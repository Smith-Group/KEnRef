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
#include "tools/Pbc.h"
#include "tools/Tensor.h"
#include "tools/Vector.h"

#include "core/kabsch.h"             // Kabsch_Umeyama
#include "core/IoUtils.h"            // getAtomMappingFromPdb, should_handleNames, fill_*
#include "core/buildModelIndexing.h" // the shared model + sub-indexing setup (de-dup of fillParamsStep0)

#include "plumedinterface/KEnRefBias.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

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
         * (see PLMD::function::Ensemble). Mirrors the gmx_fatal in KEnRefForceProvider.cpp.
         *
         * THIS GUARD IS LOAD-BEARING FOR CODE THAT IS NOT HERE. KEnRefBias::calculate() reads
         * multi_sim_comm.Get_size() UNCONDITIONALLY (KEnRefBias.cpp, the "refresh per-step replica
         * state" line). That read is correct ONLY because this guard holds: one rank per replica
         * means every rank IS its replica's main rank, so multi_sim_comm is always set wherever that
         * line runs. Delete this guard and that line starts reading an unset communicator, which
         * reports "replica 0 of 1" in silence rather than failing -- the exact defect this guard was
         * written to prevent, reintroduced at a different site. Whoever lifts the guard owns that
         * line too. (The remaining multi_sim_comm uses -- the Gather and Scatter -- are already
         * gated behind `if (!isMultiSim_) return`, so this is the only unconditional one.) */
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

            /* ---- the spanning tree that makeRequestedAtomsWhole() walks ----
             *
             * Same reference PDB, but now for ALL of atoms_ (guide + sub), not just the guide set: the
             * repair has to place every atom the restraint reads.
             *
             * Why a Euclidean minimum spanning tree and not "nearest image to my predecessor in the
             * list" or "nearest image to a shared anchor": both of those were implemented on the
             * GROMACS side and both BROKE already-whole structures. The requested atoms are SPARSE --
             * consecutive entries of the guide list are 117 atoms apart in the molecule -- so list
             * order says nothing about proximity; and GB3 reaches further from its centroid (2.171 nm)
             * than half its shortest box vector (2.165 nm), past which "nearest image to an anchor" is
             * simply the wrong image. An MST has neither failure: every edge joins genuine spatial
             * neighbours, so every edge is short, so every minimum image along it is unambiguous.
             *
             * Prim's algorithm, O(n^2) over a few hundred atoms, once, at setup. */
            const int nAtoms = static_cast<int>(atoms_.size());
            std::vector<Vector> refPos(nAtoms);
            for (int i = 0; i < nAtoms; ++i)
                refPos[i] = pdb.getPosition(atoms_[i]); // nm

            referenceTree_.clear();
            longestReferenceTreeEdge_ = 0.0;
            if (nAtoms > 1) {
                referenceTree_.reserve(static_cast<size_t>(nAtoms) - 1);
                std::vector<char> inTree(nAtoms, 0);
                std::vector<double> best(nAtoms, std::numeric_limits<double>::max());
                std::vector<int> bestParent(nAtoms, 0);
                inTree[0] = 1;
                for (int j = 1; j < nAtoms; ++j)
                    best[j] = modulo2(delta(refPos[0], refPos[j]));
                for (int added = 1; added < nAtoms; ++added) {
                    // The unattached atom closest to the tree joins next; that keeps every edge short.
                    int next = -1;
                    for (int j = 0; j < nAtoms; ++j)
                        if (!inTree[j] && (next < 0 || best[j] < best[next]))
                            next = j;
                    inTree[next] = 1;
                    /* Emitted in attachment order, which is exactly parents-before-children: `next`
                     * attaches to bestParent[next], already in the tree and therefore already emitted
                     * (or the root). makeRequestedAtomsWhole() relies on that ordering. */
                    referenceTree_.emplace_back(bestParent[next], next);
                    longestReferenceTreeEdge_ = std::max(longestReferenceTreeEdge_, std::sqrt(best[next]));
                    for (int j = 0; j < nAtoms; ++j) {
                        if (inTree[j])
                            continue;
                        const double d2 = modulo2(delta(refPos[next], refPos[j]));
                        if (d2 < best[j]) {
                            best[j] = d2;
                            bestParent[j] = next;
                        }
                    }
                }
            }
            wholePositions_.assign(static_cast<size_t>(nAtoms), Vector(0.0, 0.0, 0.0));

            /* The health check on the reference. Every edge should join spatial neighbours, so a few
             * Angstrom is normal. A value approaching half the box means the reference is ITSELF split
             * across the boundary, in which case this tree links the wrapped fragment by a long edge,
             * minimum-image finds nothing to correct, and the repair below silently does nothing. That
             * is precisely the state res/PBC-BROKEN.md describes for the un-repaired GB3_27_10us.pdb
             * (29.17 A edge against a 30.61 A half-box), and it is why that file was made whole. */
            std::cout << "  KEnRef periodic repair: spanning tree over " << nAtoms
                      << " requested atoms, longest reference edge " << longestReferenceTreeEdge_ * 10
                      << " Angstrom" << std::endl;
        }

        // ---- construct the driver (saturate only if requested; else an infinite threshold = no clamp) ----
        const KEnRef_Real_t maxForceSquared = saturate_forces_
                ? maxForce_ * maxForce_ : std::numeric_limits<KEnRef_Real_t>::infinity();
        driver_ = std::make_unique<::kenref::KEnRefDriver<KEnRef_Real_t>>(
            std::move(mi.model), k_, n_, maxForceSquared, guideAtomsReferenceCoordsCentered_);
        /* Refuse a molecule split across a periodic boundary rather than refining a torn structure.
         *
         * This matters MORE on the PLUMED side than on the GROMACS one. GROMACS repairs the
         * coordinates from the topology before the driver sees them; PLUMED has no topology here and
         * cannot, so this check is the only thing standing between a wrapped input and a silently
         * wrong bias. PLUMED's own makeWhole() is not a substitute: it walks a spanning tree built
         * from the MOLINFO *reference* coordinates, so when the reference is broken in the same way as
         * the frames -- the usual case, both coming from one trajectory -- the fragment is linked by
         * an edge shorter than half the box and minimum-image leaves it exactly where it was.
         *
         * Enabled here, in the KEnRef-side constructor, so it needs no change to the stable frame in
         * the PLUMED repository (see the H3 split). */
        driver_->enablePeriodicSplitCheck();

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

    /* Rebuild wholePositions_ for the current step. See the header for why this exists at all.
     *
     * Hosted here, beside the tree that drives it, rather than in the fork's frozen frame: the frame
     * only needs to CALL it (getLocalModelX) and read the result (getGuideAtomsX/fillSubAtomsX). */
    void KEnRefBias::makeRequestedAtomsWhole() const {
        const long step = getStep();
        if (wholePositionsStep_ == step)
            return;   // already rebuilt this step; getGuideAtomsX() also runs after the driver, for rmsd
        const bool firstWholePositions = wholePositionsStep_ < 0;
        wholePositionsStep_ = step;

        /* The reference tree is only USABLE if every one of its edges is shorter than the reach of
         * minimum-image reasoning -- half the smallest box width. Past that the nearest image of a
         * child to its parent is simply not the right one, and the walk below re-images atoms that
         * were correctly placed instead of repairing the ones that were not.
         *
         * That is exactly what a BROKEN reference produces: its wrapped fragment has no near
         * neighbour, so the tree reaches across the box to link it, and every frame is then walked
         * through that one bad edge.
         *
         * THIS GUARD IS NOT REDUNDANT, AND THE REASON IS COUNTERINTUITIVE. A bad edge does not
         * reliably give a wrong answer -- it gives an arbitrary one, and on the committed GB3 sets it
         * happens to give the RIGHT one. Measured:
         *
         *   this tree, over the 458 REQUESTED atoms   longest edge 31.49 A  -> re-images, CORRECT
         *   PLUMED makeWhole(), over all 862 atoms    longest edge 29.17 A  -> no-op,     TORN
         *
         * Both against the same 30.61 A half-box-vector: one lands either side of it. The sparser set
         * has to reach further to link the fragment, and reaching further is what accidentally saved
         * it. So with a broken reference this code silently returned the correct energy -- and nothing
         * downstream could have told us otherwise, because the driver's split-check sees a whole
         * structure and passes. The next fixture, or the next box, flips the coin the other way.
         *
         * Refusing on the EDGE, rather than trusting the outcome, is what makes the repair honest:
         * it stops when its own precondition is violated instead of when the result happens to look
         * wrong. PLUMED's own makeWhole() has the same weakness and no such guard.
         *
         * The box is not known until a step runs -- PLUMED hands it over per step, never at setup --
         * so the setup-time log records the edge and this is where it is judged. Checked once: the
         * reference does not change, and the box would have to shrink by a factor of several for a
         * healthy structure to reach this (the repaired GB3 reference gives 2.73 A against a 21.65 A
         * half-width, a factor of 7.9 of margin). */
        if (firstWholePositions) {
            const Tensor &b = getPbc().getBox();
            double smallestWidth = std::numeric_limits<double>::max();
            for (int k = 0; k < 3; ++k)
                if (b[k][k] > 0)
                    smallestWidth = std::min(smallestWidth, b[k][k]);
            if (smallestWidth < std::numeric_limits<double>::max()
                && longestReferenceTreeEdge_ >= smallestWidth / 2) {
                error("the REFERENCE structure (" + reference_pdb_ + ") cannot be used to repair "
                      "periodic images: the spanning tree over the requested atoms has an edge of "
                      + std::to_string(longestReferenceTreeEdge_ * 10) + " Angstrom, which is not less "
                      "than half the smallest box width (" + std::to_string(smallestWidth * 10 / 2)
                      + " Angstrom), so the minimum image along it is not unambiguous.\n"
                      "This almost always means the REFERENCE is itself split across a periodic "
                      "boundary: a wrapped fragment has no near neighbour, so the tree reaches across "
                      "the box to link it. Make the reference whole -- `gmx trjconv -pbc whole` with a "
                      "topology -- and check the frames too. KEnRef restrains a WHOLE molecule (global "
                      "fit, pair distances with no cutoff), so continuing would silently give a "
                      "different energy rather than fail.");
            }
        }

        const std::vector<Vector> &pos = getPositions();
        /* Start from the raw positions and only ever overwrite an atom that genuinely needs a
         * different image. Whole input therefore comes back bit-for-bit identical, which is what makes
         * it safe to run this unconditionally instead of trying to detect breakage first. */
        wholePositions_.assign(pos.begin(), pos.end());

        for (const auto &[parent, child]: referenceTree_) {
            const Vector &anchor = wholePositions_[parent];   // already placed: parents come first
            /* PLUMED's own minimum-image, so triclinic boxes are handled by the code that owns that
             * problem rather than by hand-rolled arithmetic here. pbcDistance(a, b) is the shortest
             * b - a, and it is invariant to which image of b we hand it. */
            const Vector shortest = pbcDistance(anchor, pos[child]);
            const Vector raw = delta(anchor, pos[child]);
            /* shortest - raw is exactly a lattice vector: zero when the atom is already in the right
             * image, and otherwise at least a box vector long. The tolerance below sits many orders of
             * magnitude above the round-off in that subtraction and below any real box, so this is a
             * decision about lattice vectors, not a floating-point threshold in disguise. */
            if (modulo2(shortest - raw) > 1e-12)
                wholePositions_[child] = anchor + shortest;
        }
    }

} // namespace PLMD::kenref
