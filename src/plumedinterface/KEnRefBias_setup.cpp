/*
 * KEnRefBias_setup.cpp  (KEnRef repo — kenref_plumed)
 *
 *  The VOLATILE half of KEnRefBias: its constructor (the one-time model + sub-indexing + driver setup
 *  that mirrors fillParamsStep0 and evolves with the model-abstraction). It is hosted here in the
 *  KEnRef repository so it can be updated WITHOUT re-pushing the PLUMED fork; the fork's frozen frame
 *  (src/kenref/KEnRefBias_setup.cpp) just #includes this file and compiles it within PLUMED's build.
 *
 *  Everything STABLE about the action — the class declaration, registerKeywords, the PLUMED<->Eigen
 *  glue (positions/box/force/MPI), calculate(), the boilerplate, and PLUMED_REGISTER_ACTION — lives in
 *  the fork's own KEnRefBias.{h,cpp}, so PLUMED reviewers see a real, idiomatic action there.
 */

#include "tools/Communicator.h"
#include "tools/OpenMP.h"
#include "tools/PDB.h"
#include "tools/Vector.h"

#include "core/kabsch.h"             // Kabsch_Umeyama
#include "core/IoUtils.h"            // getAtomMappingFromPdb, should_handleNames, fill_*
#include "core/buildModelIndexing.h" // the shared model + sub-indexing setup (de-dup of fillParamsStep0)

#include "plumedinterface/KEnRefBias.h"

#include <limits>

namespace PLMD::kenref {

    KEnRefBias::KEnRefBias(const ActionOptions &ao) : PLUMED_BIAS_INIT(ao), ActionAtomistic(ao) {
        /* KEnRef supports exactly ONE RANK PER REPLICA, and refuses anything else here rather than
         * producing wrong numbers. `comm` is the INTRA-replica communicator, so comm.Get_size() > 1
         * means this replica's atoms are split across ranks -- GROMACS domain decomposition, or
         * `plumed driver` without --multi, which splits atoms the same way.
         *
         * Why this is fatal rather than merely unsupported: the inter-replica communicator is only
         * set on a replica's MAIN rank (PLUMED's GROMACS patch calls "GREX setMPIIntercomm" under
         * `if (MAIN(cr))`, and GREX.cpp is what populates multi_sim_comm). An unset Communicator is
         * MPI_COMM_SELF, whose Get_size() silently returns 1 and Get_rank() 0 -- it never errors. So
         * every non-main rank would conclude "I am replica 0 of a 1-replica ensemble", compute a
         * SINGLE-REPLICA restraint, and apply it to the atoms it homes, while the main rank applied
         * the real ensemble forces to its own. The result is a spatial patchwork of two different
         * force fields, with no crash and a plausible-looking energy.
         *
         * The supported multi-replica mode is unaffected: one rank per replica leaves comm at size 1
         * on every rank, and `plumed driver --multi N` under `mpirun -np N` gives nintra == 1.
         *
         * Removing this guard is the acceptance test for real DD support, not an optimisation; the
         * fix is to gate every multi_sim_comm access on comm.Get_rank()==0 and propagate down `comm`
         * (see PLMD::function::Ensemble). Mirrors the gmx_fatal in KEnRefForceProvider.cpp. */
        if (comm.Get_size() > 1)
            error("KENREF does not support domain decomposition: this replica spans "
                  + std::to_string(comm.Get_size()) + " ranks. Run KEnRef with ONE RANK PER REPLICA "
                  "-- e.g. `mpirun -np <N> ... -multidir <N directories>`, which gives each replica a "
                  "single rank -- or reduce the rank count so no replica is decomposed. Threads are "
                  "unaffected: -ntomp / PLUMED_NUM_THREADS are free to use.");

        // ---- general parameters ----
        parse("MODEL", modelName_);
        ::kenref::bootstrapModels();
        if (!::kenref::ModelRegistry<KEnRef_Real_t>::has(modelName_))
            error("MODEL '" + modelName_ + "' is not a registered KEnRef energy model");
        parse("K", k_);
        parse("N", n_);
        parse("MAX_FORCE", maxForce_);
        parseFlag("FIT_TO_REFERENCE", fit_to_reference_);
        parseFlag("SATURATE_FORCES", saturate_forces_);

        // ---- model-tier parameters: parse every registered model keyword into modelParams_ ----
        for (const auto &[name, schema]: ::kenref::ModelRegistry<KEnRef_Real_t>::allSchemas()) {
            for (const auto &spec: schema.specs()) {
                if (modelParams_.count(spec.key))
                    continue;
                std::string v;
                parse(spec.key, v);
                if (!v.empty())
                    modelParams_[spec.key] = v;
            }
        }

        // ---- atom-name mapping (required NOW: the sub-atom list is derived from the data files) ----
        parse("ATOMNAME_MAPPING", atomname_mapping_file_);
        if (atomname_mapping_file_.empty() || atomname_mapping_file_ == "undefined")
            error("ATOMNAME_MAPPING is required");
        parse("REF", reference_pdb_);
        if (reference_pdb_.empty()) {
            reference_pdb_ = atomname_mapping_file_;
            std::cout << "  REF not specified - using ATOMNAME_MAPPING file as reference: " << reference_pdb_ << std::endl;
        }

        atomName_to_globalSerial_map_ = IoUtils::getAtomMappingFromPdb<std::string, int>(
            atomname_mapping_file_, IoUtils::fill_atomId_to_index_Map);
        if (atomName_to_globalSerial_map_.empty())
            error("No atom mapping found in ATOMNAME_MAPPING file: " + atomname_mapping_file_);
        const bool handleNames = IoUtils::should_handleNames(atomName_to_globalSerial_map_);

        // Replica count for the model's grouping (PLATEAUS). multi_sim_comm is set before construction.
        const int numSimulations = (Communicator::plumedHasMPI() && Communicator::initialized())
                                       ? multi_sim_comm.Get_size() : 1;

        // ---- shared model + sub-atom indexing via kenref::buildModelIndexing (the SAME path GROMACS
        //      and the offline tools use). The 1-based PDB serials come back as 1-based subAtomGlobalIds;
        //      buildCache + finalizeIndexing happen inside, so the model is ready after this call. ----
        auto mi = ::kenref::buildModelIndexing<KEnRef_Real_t>(
            modelName_, *this, atomName_to_globalSerial_map_, handleNames, numSimulations,
            static_cast<int>(OpenMP::getNumThreads()));
        sub0Id_to_global1Serial_ = std::move(mi.subAtomGlobalIds);
        numSubAtoms_ = static_cast<int>(sub0Id_to_global1Serial_.size());
        if (numSubAtoms_ == 0)
            error("No restrained atoms derived from experimental data - check the data files and ATOMNAME_MAPPING");
        subAtoms_.clear();
        for (int s: sub0Id_to_global1Serial_)
            subAtoms_.emplace_back(AtomNumber().setSerial(s));

        // ---- guide atoms + the merged request list (guide first, then sub) ----
        parseAtomList("GUIDE_ATOMS", guideAtoms_);
        atoms_.clear();
        atoms_.insert(atoms_.end(), guideAtoms_.begin(), guideAtoms_.end());
        atoms_.insert(atoms_.end(), subAtoms_.begin(), subAtoms_.end());
        requestAtoms(atoms_);

        // serial -> local index in getPositions() (atoms_ order)
        serial_to_localIdx_.clear();
        for (int li = 0; li < static_cast<int>(atoms_.size()); ++li)
            serial_to_localIdx_[atoms_[li].serial()] = li;

        // ---- reference structure (Angstrom) + centered copy for the driver ----
        {
            PDB pdb;
            if (!pdb.read(reference_pdb_, usingNaturalUnits(), 0.1))
                error("Cannot read reference PDB: " + reference_pdb_);
            const int nGuide = static_cast<int>(guideAtoms_.size());
            guideAtomsReferenceCoords_.resize(nGuide, 3);
            for (int i = 0; i < nGuide; ++i) {
                Vector pos = pdb.getPosition(guideAtoms_[i]); // nm
                guideAtomsReferenceCoords_(i, 0) = to<KEnRef_Real_t>(pos[0]) * 10; // -> Angstrom
                guideAtomsReferenceCoords_(i, 1) = to<KEnRef_Real_t>(pos[1]) * 10;
                guideAtomsReferenceCoords_(i, 2) = to<KEnRef_Real_t>(pos[2]) * 10;
            }
            guideAtomsReferenceCoordsCentered_ =
                    Kabsch_Umeyama<KEnRef_Real_t>::translateCenterOfMassToOrigin(guideAtomsReferenceCoords_);
        }

        // ---- construct the driver (saturate only if requested; else an infinite threshold = no clamp) ----
        const KEnRef_Real_t maxForceSquared = saturate_forces_
                ? maxForce_ * maxForce_ : std::numeric_limits<KEnRef_Real_t>::infinity();
        driver_ = std::make_unique<::kenref::KEnRefDriver<KEnRef_Real_t>>(
            std::move(mi.model), k_, n_, maxForceSquared, guideAtomsReferenceCoordsCentered_);

        // ---- gather/scatter scratch (allocated once) ----
        const size_t subSize = static_cast<size_t>(numSubAtoms_) * 3;
        const size_t nSims = static_cast<size_t>(std::max(1, numSimulations));
        allSimulationsSubAtomsX_buffer_.assign(subSize * nSims, 0);
        allDerivatives_buffer_.assign(subSize * nSims, 0);
        derivatives_buffer_.assign(subSize, 0);

        // ---- output components ----
        addComponent("energy");
        componentIsNotPeriodic("energy");
        addComponent("rmsd");
        componentIsNotPeriodic("rmsd");

        checkRead();

        std::cout << "  KEnRef bias  model=" << modelName_ << "  k=" << k_ << "  n=" << n_ << "\n";
        std::cout << "  " << guideAtoms_.size() << " guide atoms, " << subAtoms_.size() << " restrained atoms\n";

        cite("Restraining interproton angular and distance dynamics with KEnRef. "
            "Amr Alhossary & Colin Smith, J Phys Chem B. 130, 11 (2026) DOI: 10.1021/acs.jpcb.5c08554");
    }

} // namespace PLMD::kenref
