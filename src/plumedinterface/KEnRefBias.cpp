/*
 * KEnRefBias.cpp  (KEnRef repo — kenref_plumed)
 *
 *  The PLUMED EngineAdapter (see KEnRefBias.h). The old per-model SIGMA/PLATEAUS branches and the
 *  inline fit/no-jump/gather/compute/scatter/apply pipeline are gone: the model is built from the
 *  registry by name, and KEnRefDriver::step() runs the pipeline, calling back through the
 *  EngineAdapter methods implemented here for the PLUMED-specific operations.
 *
 *  PLUMED_REGISTER_ACTION lives here (this TU is compiled inside PLUMED's build by the fork's frame).
 */

#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "tools/Vector.h"
#include "tools/Tensor.h"
#include "tools/Communicator.h"
#include "tools/OpenMP.h"
#include "tools/PDB.h"

#include "core/kabsch.h"        // Kabsch_Umeyama
#include "core/IoUtils.h"       // getAtomMappingFromPdb, should_handleNames, fill_*
#include "core/ModelRegistry.h"

#include "plumedinterface/KEnRefBias.h"

#include <chrono>
#include <limits>
#include <set>

namespace PLMD::kenref {

    PLUMED_REGISTER_ACTION(KEnRefBias, "KENREF")

    // ============================================================
    //  registerKeywords
    // ============================================================
    void KEnRefBias::registerKeywords(Keywords &keys) {
        Bias::registerKeywords(keys);
        ActionAtomistic::registerKeywords(keys);
        keys.setDisplayName("KENREF");

        // ---- General (framework-owned) keywords ----
        keys.add("compulsory", "MODEL", "SIGMA", "The energy model to use (e.g. SIGMA or PLATEAUS)");
        keys.add("compulsory", "K", "1.0", "Force constant");
        keys.add("compulsory", "N", "0.25", "Power scaling factor");
        keys.add("optional", "MAX_FORCE", "Maximum force magnitude (default 9999)");
        keys.add("atoms", "GUIDE_ATOMS", "Atoms used for alignment to reference");
        keys.add("optional", "REF", "Reference structure PDB for alignment");
        keys.add("compulsory", "ATOMNAME_MAPPING", "undefined",
                 "PDB file with atom-name -> atom-index mapping (also used as reference if REF is omitted)");
        keys.addFlag("FIT_TO_REFERENCE", false, "Fit coordinates to reference before calculating energy");
        keys.addFlag("SATURATE_FORCES", false, "Clamp forces to MAX_FORCE");

        // ---- Model-specific keywords (the union across every registered model, optional; the
        //      selected model validates required-ness in buildCache). Deduplicated by key. ----
        ::kenref::bootstrapModels();
        std::set<std::string> seen;
        for (const auto &[name, schema]: ::kenref::ModelRegistry<KEnRef_Real_t>::allSchemas()) {
            for (const auto &spec: schema.specs()) {
                if (seen.insert(spec.key).second)
                    keys.add("optional", spec.key, spec.doc);
            }
        }

        // ---- Output components ----
        keys.addOutputComponent("energy", "default", "scalar", "Total KEnRef restraint energy");
        keys.addOutputComponent("rmsd", "default", "scalar", "RMSD from reference (after fitting)");

        keys.use("RESTART");
        keys.use("UPDATE_FROM");
        keys.use("UPDATE_UNTIL");
    }

    // ============================================================
    //  Constructor  (replaces the old constructor + initializeParameters/fillParamsStep0)
    // ============================================================
    KEnRefBias::KEnRefBias(const ActionOptions &ao) : PLUMED_BIAS_INIT(ao), ActionAtomistic(ao) {
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

        // ---- build the selected model from the registry and let it load its inputs ----
        auto model = ::kenref::ModelRegistry<KEnRef_Real_t>::create(modelName_);
        {
            ::kenref::ParamSchema schema; // empty: this adapter supplies every param via getRawParam
            ::kenref::InitContext<KEnRef_Real_t> initCtx{*this, schema, atomName_to_globalSerial_map_,
                                                         handleNames, numSimulations};
            model->buildCache(initCtx);
        }

        // ---- derive subAtoms_ from the model's atom-name pairs (sorted unique 1-based serials) ----
        std::set<int> subSerials;
        for (const auto &[a1, a2]: model->atomNamePairs()) {
            subSerials.insert(atomName_to_globalSerial_map_.at(a1));
            subSerials.insert(atomName_to_globalSerial_map_.at(a2));
        }
        if (subSerials.empty())
            error("No restrained atoms derived from experimental data - check the data files and ATOMNAME_MAPPING");
        sub0Id_to_global1Serial_.assign(subSerials.begin(), subSerials.end());
        numSubAtoms_ = static_cast<int>(sub0Id_to_global1Serial_.size());
        subAtoms_.clear();
        for (int s: sub0Id_to_global1Serial_)
            subAtoms_.emplace_back(AtomNumber().setSerial(s));

        // atomName -> sub0Id (for finalizeIndexing)
        std::map<int, int> serial_to_sub0Id;
        for (int i = 0; i < numSubAtoms_; ++i)
            serial_to_sub0Id[sub0Id_to_global1Serial_[i]] = i;
        std::map<std::string, int> atomName_to_sub0Id_map;
        for (const auto &[name, serial]: atomName_to_globalSerial_map_) {
            const auto it = serial_to_sub0Id.find(serial);
            if (it != serial_to_sub0Id.end())
                atomName_to_sub0Id_map[name] = it->second;
        }

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

        // ---- finalize the model's indexing (name->sub0Id pairs + primed per-data caches) ----
        {
            ::kenref::IndexingContext<KEnRef_Real_t> idxCtx{atomName_to_sub0Id_map,
                                                           static_cast<int>(OpenMP::getNumThreads())};
            model->finalizeIndexing(idxCtx);
        }

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
            std::move(model), k_, n_, maxForceSquared, guideAtomsReferenceCoordsCentered_);

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

    // ============================================================
    //  Helpers: current PLUMED positions (nm) -> Angstrom Eigen
    // ============================================================
    CoordsMatrixType<KEnRef_Real_t> KEnRefBias::getGuideAtomsX() const {
        const std::vector<Vector> &pos = getPositions();
        const int nGuide = static_cast<int>(guideAtoms_.size());
        CoordsMatrixType<KEnRef_Real_t> gx(nGuide, 3);
        for (int i = 0; i < nGuide; ++i) {
            const int li = serial_to_localIdx_.at(guideAtoms_[i].serial());
            gx(i, 0) = to<KEnRef_Real_t>(pos[li][0]) * 10;
            gx(i, 1) = to<KEnRef_Real_t>(pos[li][1]) * 10;
            gx(i, 2) = to<KEnRef_Real_t>(pos[li][2]) * 10;
        }
        return gx;
    }

    void KEnRefBias::fillSubAtomsX(CoordsMatrixType<KEnRef_Real_t> &out) const {
        const std::vector<Vector> &pos = getPositions();
        for (int i = 0; i < out.rows(); ++i) {
            const int serial = sub0Id_to_global1Serial_[i];
            const int li = serial_to_localIdx_.at(serial);
            out(i, 0) = to<KEnRef_Real_t>(pos[li][0]) * 10; // nm -> Angstrom
            out(i, 1) = to<KEnRef_Real_t>(pos[li][1]) * 10;
            out(i, 2) = to<KEnRef_Real_t>(pos[li][2]) * 10;
        }
    }

    // ============================================================
    //  kenref::EngineAdapter implementation
    // ============================================================
    std::optional<std::string> KEnRefBias::getRawParam(const std::string &key) const {
        const auto it = modelParams_.find(key);
        if (it == modelParams_.end() || it->second.empty())
            return std::nullopt;
        return it->second;
    }

    int KEnRefBias::numOmpThreads() const {
        return static_cast<int>(OpenMP::getNumThreads());
    }

    void KEnRefBias::getLocalModelX(int /*localModel*/, CoordsMatrixType<KEnRef_Real_t> &guideX,
                                    CoordsMatrixType<KEnRef_Real_t> &subX,
                                    Eigen::Matrix<KEnRef_Real_t, 3, 3> &box) const {
        guideX = getGuideAtomsX();
        subX.resize(numSubAtoms_, 3);
        fillSubAtomsX(subX);
        // PLUMED box Tensor (nm) -> Eigen 3x3 (raw; the driver scales it to Angstrom via toAngstrom).
        const Tensor &b = getPbc().getBox();
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                box(i, j) = to<KEnRef_Real_t>(b[i][j]);
    }

    void KEnRefBias::addLocalModelDerivatives(int /*localModel*/, const CoordsMatrixType<KEnRef_Real_t> &derivs) {
        // derivs are already inverse-fitted, unit-scaled and saturated by the driver. KEnRef returns
        // dE/dx; the applied force is -dE/dx, so we negate.
        for (int i = 0; i < derivs.rows(); ++i) {
            const int serial = sub0Id_to_global1Serial_[i];
            const int li = serial_to_localIdx_.at(serial);
            const Vector d(static_cast<double>(derivs(i, 0)),
                           static_cast<double>(derivs(i, 1)),
                           static_cast<double>(derivs(i, 2)));
            addForce(getValueIndices(getAbsoluteIndex(li)), -d);
        }
    }

    void KEnRefBias::gatherFittedSubAtomsX(const std::vector<CoordsMatrixType<KEnRef_Real_t>> &localFitted,
                                           std::vector<CoordsMatrixType<KEnRef_Real_t>> &all) const {
        if (!isMultiSim_) {
            all = localFitted; // single replica owns the whole ensemble
            return;
        }
        const auto &fitted = localFitted[0];
        const int subSize = static_cast<int>(fitted.size()); // rows*3
        multi_sim_comm.Gather(const_cast<KEnRef_Real_t *>(fitted.data()), subSize,
                              allSimulationsSubAtomsX_buffer_.data(), subSize, 0);
        if (simulationIndex_ == 0) {
            const Eigen::Index nSub = fitted.rows();
            all.resize(numSimulations_);
            for (int m = 0; m < numSimulations_; ++m) // zero-copy view materialised by the assignment
                all[m] = CoordsMapTypeConst<KEnRef_Real_t>(&allSimulationsSubAtomsX_buffer_[m * nSub * 3], nSub, 3);
        }
    }

    void KEnRefBias::scatterModelDerivatives(const std::vector<CoordsMatrixType<KEnRef_Real_t>> &allPerModel,
                                             std::vector<CoordsMatrixType<KEnRef_Real_t>> &localPerModel) const {
        if (!isMultiSim_) {
            localPerModel = allPerModel; // single replica
            return;
        }
        const int subSize = numSubAtoms_ * 3;
        if (simulationIndex_ == 0)
            for (int m = 0; m < numSimulations_; ++m)
                std::copy_n(allPerModel[m].data(), subSize, &allDerivatives_buffer_[m * subSize]);
        multi_sim_comm.Scatter(allDerivatives_buffer_.data(), subSize, derivatives_buffer_.data(), subSize, 0);
        localPerModel.resize(1);
        localPerModel[0] = CoordsMapType<KEnRef_Real_t>(derivatives_buffer_.data(), numSubAtoms_, 3);
    }

    // ============================================================
    //  calculate()
    // ============================================================
    void KEnRefBias::calculate() {
        const auto begin = std::chrono::high_resolution_clock::now();
        const long step = getStep();

        // Refresh per-step replica state (read by the EngineAdapter callbacks the driver invokes).
        isMultiSim_ = Communicator::plumedHasMPI() && Communicator::initialized() && multi_sim_comm.Get_size() > 1;
        numSimulations_ = isMultiSim_ ? multi_sim_comm.Get_size() : 1;
        simulationIndex_ = isMultiSim_ ? multi_sim_comm.Get_rank() : 0;

        // Run the whole engine-agnostic pipeline (no-jump -> fit -> gather -> compute -> scatter ->
        // inverse-fit -> unit-scale -> saturate -> addForce) through the adapter callbacks above.
        const KEnRef_Real_t energy = driver_->step(*this, step % 10 == 0);

        // ---- RMSD (crude diagnostic; current guide atoms vs reference, Angstrom -> nm) ----
        double rmsd = 0.0;
        if (guideAtomsReferenceCoords_.rows() > 0) {
            const CoordsMatrixType<KEnRef_Real_t> diff = getGuideAtomsX() - guideAtomsReferenceCoords_;
            rmsd = std::sqrt(diff.rowwise().squaredNorm().mean()) * 0.1;
        }

        // ---- output components ----
        const double bias_value = energy;
        getPntrToComponent("bias")->set(bias_value);
        getPntrToComponent("energy")->set(bias_value);
        getPntrToComponent("rmsd")->set(rmsd);
        setBias(bias_value);

        const auto end = std::chrono::high_resolution_clock::now();
        calculate_time_ += std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count();
        if (simulationIndex_ == 0 && step % 10 == 0)
            std::cout << "  Step " << step << "  KEnRef energy = " << energy << std::endl;
    }

    // ============================================================
    //  Ambiguity resolution (lock both Atoms and Arguments/CVs)
    // ============================================================
    void KEnRefBias::lockRequests() {
        ActionAtomistic::lockRequests();
        Bias::lockRequests();
    }

    void KEnRefBias::unlockRequests() {
        ActionAtomistic::unlockRequests();
        Bias::unlockRequests();
    }

    void KEnRefBias::calculateNumericalDerivatives(ActionWithValue *a) {
        ActionAtomistic::calculateNumericalDerivatives(a);
    }

} // namespace PLMD::kenref
