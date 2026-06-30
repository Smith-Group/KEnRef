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

#include "core/ActionAtomistic.h"
#include "bias/Bias.h"

#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "core/KEnRef.h"          // CoordsMatrixType<>, KEnRef_Real_t
#include "core/DefaultEngineAdapter.h"
#include "core/KEnRefDriver.h"

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

        // ---- per-step replica state (refreshed at the top of calculate()) ----
        int numSimulations_ = 1;
        int simulationIndex_ = 0;
        bool isMultiSim_ = false;

        // ---- gather/scatter scratch (mutable: written by the const EngineAdapter callbacks) ----
        mutable std::vector<KEnRef_Real_t> allSimulationsSubAtomsX_buffer_;
        mutable std::vector<KEnRef_Real_t> allDerivatives_buffer_;
        mutable std::vector<KEnRef_Real_t> derivatives_buffer_;

        long long calculate_time_ = 0;

        // ---- private helpers (current PLUMED positions -> Angstrom Eigen) ----
        [[nodiscard]] CoordsMatrixType<KEnRef_Real_t> getGuideAtomsX() const;
        void fillSubAtomsX(CoordsMatrixType<KEnRef_Real_t>& out) const;

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
        [[nodiscard]] int simulationIndex() const override { return simulationIndex_; }
        void gatherFittedSubAtomsX(const std::vector<CoordsMatrixType<KEnRef_Real_t>>& localFitted,
                                   std::vector<CoordsMatrixType<KEnRef_Real_t>>& all) const override;
        void scatterModelDerivatives(const std::vector<CoordsMatrixType<KEnRef_Real_t>>& allPerModel,
                                     std::vector<CoordsMatrixType<KEnRef_Real_t>>& localPerModel) const override;
    };

} // namespace PLMD::kenref

#endif // PLUMED_kenref_KEnRefBias_h
