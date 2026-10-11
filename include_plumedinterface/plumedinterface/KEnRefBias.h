#ifndef PLUMED_kenref_KEnRefBias_h
#define PLUMED_kenref_KEnRefBias_h

/*
 * KEnRefBias.h  (KEnRef repo — kenref_plumed)
 *
 *  ROLE: the PLUMED consumer of kenref_core. After the model-abstraction restructure it is a thin
 *  EngineAdapter: it owns a KEnRefDriver (which runs the whole per-step pipeline) and implements the
 *  handful of PLUMED-specific callbacks (positions/box, force application, inter-replica gather/scatter).
 *  It carries NO per-model logic — the selected EnergyModel is built from the registry by name.
 *
 *  This source lives in the KEnRef repo (permissive license, single source of truth). The PLUMED fork's
 *  src/kenref/ "frozen frame" merely compiles it (within PLUMED's build) and registers the action.
 */

#include <map>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "core/KEnRef.h"          // CoordsMatrixType<>, KEnRef_Real_t
#include "core/DefaultEngineAdapter.h"
#include "core/KEnRefDriver.h"

#include "core/ActionAtomistic.h"
#include "bias/Bias.h"

// Narrowing/static-cast helper used for PLUMED(double) <-> KEnRef_Real_t conversions.
template<typename To, typename From>
To to(From v) noexcept { return static_cast<To>(v); }

namespace PLMD::kenref {

    class KEnRefBias final : public bias::Bias,
                            public ActionAtomistic,
                            public ::kenref::DefaultEngineAdapter<KEnRef_Real_t> {
        // ---- general (framework) parameters ----
        KEnRef_Real_t k_ = 1.0;
        KEnRef_Real_t n_ = 0.25;
        KEnRef_Real_t maxForce_ = 9999.0;
        bool fit_to_reference_ = true;   // kept for rmsd/back-compat; the driver always fits
        bool saturate_forces_ = false;
        std::string modelName_;          // registry model name (e.g. SIGMA / PLATEAUS)

        // ---- model-tier parameter values (backs getRawParam) ----
        std::map<std::string, std::string> modelParams_;

        // ---- atom lists / name mapping ----
        std::vector<AtomNumber> atoms_;       // full list passed to requestAtoms() = guide + sub
        std::vector<AtomNumber> guideAtoms_;  // alignment atoms
        std::vector<AtomNumber> subAtoms_;    // restrained atoms
        std::string atomname_mapping_file_, reference_pdb_;
        std::map<std::string, int> atomName_to_globalSerial_map_;  // normalised name -> 1-based serial
        std::vector<int> sub0Id_to_global1Serial_;                 // compact sub-id -> 1-based serial
        std::map<int, int> serial_to_localIdx_;                    // serial -> index in getPositions()

        // ---- reference structure (Kabsch) ----
        CoordsMatrixType<KEnRef_Real_t> guideAtomsReferenceCoords_;          // Angstrom (for rmsd)
        CoordsMatrixType<KEnRef_Real_t> guideAtomsReferenceCoordsCentered_;  // passed to the driver

        // ---- the shared engine-agnostic pipeline ----
        std::unique_ptr<::kenref::KEnRefDriver<KEnRef_Real_t>> driver_;
        int numSubAtoms_ = 0;

        /* ---- replica state, fixed for the whole run and set ONCE at setup ----
         *
         * These used to be refreshed every step from `multi_sim_comm`, which was only correct while
         * domain decomposition was refused: that communicator is populated only on a replica's MAIN
         * rank (PLUMED's GROMACS patch calls "GREX setMPIIntercomm" under `if (MAIN(cr))`), and an
         * unset Communicator is MPI_COMM_SELF -- Get_size() returns 1 and Get_rank() 0 without ever
         * erroring. Every non-main rank therefore concluded "replica 0 of 1" and restrained its own
         * domain against a single-replica ensemble. PLUMED #1443.
         *
         * They now come from Action::getNumberOfReplicas() / getReplicaIndex() (PLUMED #1481), which
         * are correct on EVERY rank and need no communication, so they can be read at setup. */
        int numSimulations_ = 1;
        int replicaIndex_ = 0;
        bool isMultiSim_ = false;
        /*! \brief Whether this rank may use `multi_sim_comm` -- true on exactly one rank per replica.
         *
         * `multi_sim_comm` is valid only on a replica's main rank, so every collective over it is
         * gated on this. Distinct from "did I compute the ensemble answer": see simulationIndex(). */
        bool isReplicaMainRank_ = true;

        /* ---- periodic repair of the requested atoms ----
         *
         * KEnRef restrains a WHOLE molecule: the Kabsch fit and the pair distances are global and have
         * no cutoff, so a molecule split across a periodic boundary does not fail -- it silently gives a
         * different energy and different forces. GROMACS repairs its coordinates from the topology's
         * bond graph; PLUMED is never handed a topology, so the tree is built instead from the REFERENCE
         * structure, which is a Euclidean minimum spanning tree over exactly the atoms we requested.
         *
         * That is only sound because the reference is itself whole -- see res/PBC-BROKEN.md. A broken
         * reference has no near neighbour for its wrapped fragment, so the tree reaches across the box
         * to link it and every frame is walked through that one bad edge. Whether the result is then
         * wrong is a COIN FLIP, not a certainty -- on the committed GB3 sets a broken reference happens
         * to give the correct answer -- so makeRequestedAtomsWhole() refuses on the EDGE rather than
         * trusting the outcome. The full measurement is with that refusal in KEnRefBias_setup.cpp.
         *
         * The driver's split-check still stands behind all of this: repair first, refuse second. */

        //! Spanning tree over atoms_, as (parent, child) LOCAL indices, parents always before children.
        std::vector<std::pair<int, int>> referenceTree_;
        //! Longest edge of referenceTree_ in the reference structure, nm. A few A is healthy.
        double longestReferenceTreeEdge_ = 0.0;
        //! getPositions(), with every atom moved to the image nearest its tree parent. Rebuilt per step.
        mutable std::vector<Vector> wholePositions_;
        //! The step wholePositions_ was last rebuilt for; -1 = never. Makes the rebuild idempotent.
        mutable long wholePositionsStep_ = -1;

        // ---- gather/scatter scratch (mutable: written by the const EngineAdapter callbacks) ----
        mutable std::vector<KEnRef_Real_t> allSimulationsSubAtomsX_buffer_;
        mutable std::vector<KEnRef_Real_t> allDerivatives_buffer_;
        mutable std::vector<KEnRef_Real_t> derivatives_buffer_;

        long long calculate_time_ = 0;

        // ---- private helpers (current PLUMED positions -> Angstrom Eigen) ----
        [[nodiscard]] CoordsMatrixType<KEnRef_Real_t> getGuideAtomsX() const;
        void fillSubAtomsX(CoordsMatrixType<KEnRef_Real_t>& out) const;

        /*! \brief Refresh wholePositions_ for the current step: every requested atom at the image
         *         nearest its parent in referenceTree_.
         *
         * Walks parents-before-children, so each atom is placed relative to one already-placed
         * neighbour and the minimum image is unambiguous. A BIT-FOR-BIT no-op on input that is already
         * whole -- no tree edge then needs a shift, and nothing is written -- which is what lets it run
         * unconditionally, exactly as the GROMACS-side repair does. Cheap to call more than once in a
         * step: it remembers which step it built. */
        void makeRequestedAtomsWhole() const;

    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        explicit KEnRefBias(const ActionOptions&);
        static void registerKeywords(Keywords&);
        void calculate() override;

        // Explicit overrides to resolve the bias::Bias / ActionAtomistic multiple-inheritance ambiguity.
        void lockRequests() override;
        void unlockRequests() override;
        void calculateNumericalDerivatives(ActionWithValue* a) override;

        // ---- kenref::EngineAdapter -------------------------------------------------------------
        [[nodiscard]] std::optional<std::string> getRawParam(const std::string& key) const override;
        [[nodiscard]] int numOmpThreads() const override;
        // numModelsInThisProcess() == 1 (one replica per process) comes from DefaultEngineAdapter.
        void getLocalModelX(int localModel, CoordsMatrixType<KEnRef_Real_t>& guideX,
                            CoordsMatrixType<KEnRef_Real_t>& subX,
                            Eigen::Matrix<KEnRef_Real_t, 3, 3>& box) const override;
        void addLocalModelDerivatives(int localModel, const CoordsMatrixType<KEnRef_Real_t>& derivs) override;
        [[nodiscard]] int numModelsTotal()  const override { return numSimulations_; }

        /*! \brief Which replica's slot this rank owns, or -1 if it owns none.
         *
         * THE -1 IS LOAD-BEARING AND IT IS NOT A STYLE CHOICE. KEnRefDriver::step() decides who runs
         * the whole-ensemble model with `adapter.simulationIndex() == 0`. getReplicaIndex() is valid
         * on every rank, so returning it unconditionally makes that test true on EVERY RANK OF
         * REPLICA 0 -- under domain decomposition the ensemble model would run once per domain
         * instead of once per job, and each of those ranks would need the gathered `all`, which a
         * Gather rooted on the replica's main rank does not give it.
         *
         * So two different questions get two different answers, and conflating them is the
         * characteristic domain-decomposition bug in this code base (KEnRefForceProvider.h says the
         * same thing for the GROMACS engine, where the predicates are named isReplicaMasterRank and
         * isEnsembleMasterRank):
         *   isReplicaMainRank_            -- "may I use multi_sim_comm", one rank per replica
         *   simulationIndex() == 0        -- "did I compute the answer", one rank per job
         * They coincide at one rank per replica, which is why mixing them up survives testing. */
        [[nodiscard]] int simulationIndex() const override {
            return isReplicaMainRank_ ? replicaIndex_ : -1;
        }
        void gatherFittedSubAtomsX(const std::vector<CoordsMatrixType<KEnRef_Real_t>>& localFitted,
                                   std::vector<CoordsMatrixType<KEnRef_Real_t>>& all) const override;
        void scatterModelDerivatives(const std::vector<CoordsMatrixType<KEnRef_Real_t>>& allPerModel,
                                     std::vector<CoordsMatrixType<KEnRef_Real_t>>& localPerModel) const override;
    };

} // namespace PLMD::kenref

#endif // PLUMED_kenref_KEnRefBias_h
