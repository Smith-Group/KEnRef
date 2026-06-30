/*
 * DefaultEngineAdapter.h
 *
 *  ROLE: an EngineAdapter base that supplies the trivial "single process, single model, no
 *  inter-replica MPI" defaults, so each consumer overrides ONLY what it genuinely needs (the
 *  "default-implementation base" idiom). It adds no new concept — it is purely convenience over the
 *  EngineAdapter interface.
 *
 *  Who overrides what:
 *    - live single-sim MD            : override the I/O (getRawParam / getLocalModelX /
 *                                      addLocalModelDerivatives) + numOmpThreads; keep every default.
 *    - live multi-sim MD (GMX/PLUMED): additionally override the MPI cluster — numModelsTotal,
 *                                      simulationIndex, gatherFittedSubAtomsX, scatterModelDerivatives.
 *    - offline tools                 : additionally override numModelsInThisProcess / numModelsTotal
 *                                      (= N trajectories); reuse the trivial gather/scatter + index 0.
 *
 *  NOT defaulted (kept pure in EngineAdapter): getRawParam, getLocalModelX, and — deliberately —
 *  addLocalModelDerivatives, so a force-applying engine cannot silently forget to apply forces
 *  (see the note in EngineAdapter.h). Energy/metric-only consumers opt into the no-op explicitly.
 */

#ifndef KENREF_DEFAULTENGINEADAPTER_H_
#define KENREF_DEFAULTENGINEADAPTER_H_

#include <vector>

#include "EngineAdapter.h"   // EngineAdapter<>, CoordsMatrixType<>

namespace kenref {

template<typename Real>
class DefaultEngineAdapter : public EngineAdapter<Real> {
public:
    // Single model held by this process (the live-MD default; offline overrides to N).
    [[nodiscard]] int numModelsInThisProcess() const override { return 1; }

    // Single simulation / master / all threads (multi-sim engines + offline override as needed).
    [[nodiscard]] int numModelsTotal()  const override { return 1; }
    [[nodiscard]] int simulationIndex() const override { return 0; }
    [[nodiscard]] int numOmpThreads()   const override { return 0; } // 0 = use all available

    // Trivial single-process gather/scatter: this process owns the whole ensemble, so the "collect on
    // the master" / "distribute back" operations are plain moves. Multi-replica engines override.
    void gatherFittedSubAtomsX(const std::vector<CoordsMatrixType<Real>>& localFitted,
                               std::vector<CoordsMatrixType<Real>>& all) const override {
        all = localFitted;
    }
    void scatterModelDerivatives(const std::vector<CoordsMatrixType<Real>>& allPerModel,
                                 std::vector<CoordsMatrixType<Real>>& localPerModel) const override {
        localPerModel = allPerModel;
    }
};

} // namespace kenref

#endif // KENREF_DEFAULTENGINEADAPTER_H_
