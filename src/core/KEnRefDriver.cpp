/*
 * KEnRefDriver.cpp
 *
 *  The shared per-step refinement pipeline (see KEnRefDriver.h). This is the engine-agnostic half of
 *  what KEnRefForceProvider::calculateForces and KEnRefBias::calculate used to each implement; the
 *  per-model dispatch is now a single model_->compute() and the engine-specific I/O is the adapter.
 */

#include "core/KEnRefDriver.h"

#include "core/kabsch.h"
#include "core/restore_no_jump.h"

namespace kenref {

template<typename Real>
Real KEnRefDriver<Real>::step(EngineAdapter<Real>& adapter, bool printStatistics) {
    const int numLocal = adapter.numModelsInThisProcess();
    // Queried once per step: the engine may change its thread count at runtime (e.g. GROMACS).
    const int numOmpThreads = adapter.numOmpThreads();

    affines_.resize(numLocal);
    localFitted_.resize(numLocal);
    if (firstStep_) {
        lastGuide_.resize(numLocal);
        lastSub_.resize(numLocal);
    }

    // ---- per local model: fetch, no-jump, fit to the reference --------------------------------
    CoordsMatrixType<Real> guideX, subX;
    Eigen::Matrix<Real, 3, 3> box;
    for (int i = 0; i < numLocal; ++i) {
        adapter.getLocalModelX(i, guideX, subX, box);

        // No-jump correction against the previous frame. On the very first step there is no previous
        // frame, so we skip it — equivalent to the old force provider, which primed lastFrame* to the
        // step-0 frame in fillParamsStep0 (making its step-0 restoreNoJump a no-op).
        // N.B. no-jump must run BEFORE find3DAffineTransform/applyTransform/applyInverseOfTransform.
        if (!firstStep_) {
            restoreNoJump(guideX, lastGuide_[i], box, /*toAngstrom*/ true, numOmpThreads, printStatistics);
            restoreNoJump(subX,   lastSub_[i],   box, /*toAngstrom*/ true, numOmpThreads, printStatistics);
        }
        lastGuide_[i] = guideX;
        lastSub_[i]   = subX;

        affines_[i]     = Kabsch_Umeyama<Real>::find3DAffineTransform(guideX, referenceGuideAtomsCoordsCentered_, false, false, false);
        localFitted_[i] = Kabsch_Umeyama<Real>::applyTransform(affines_[i], subX);
    }
    firstStep_ = false;

    // ---- gather every model's fitted coords onto the master, compute, scatter derivatives back ---
    adapter.gatherFittedSubAtomsX(localFitted_, allModels_);

    Real energy = Real(0);
    if (adapter.simulationIndex() == 0) {
        StepContext<Real> ctx{allModels_, k_, n_, /*gradient*/ true, numOmpThreads};
        auto [e, grad] = model_->compute(ctx);
        energy = e;
        perModelDerivs_ = std::move(grad.value());
    }

    adapter.scatterModelDerivatives(perModelDerivs_, localDerivs_);

    // ---- per local model: inverse-fit, unit-scale, saturate, apply -----------------------------
    const Real scale = model_->forceUnitScale();
    for (int i = 0; i < numLocal; ++i) {
        CoordsMatrixType<Real> d = Kabsch_Umeyama<Real>::applyInverseOfTransform(affines_[i], localDerivs_[i]);
        if (scale != Real(1)) {
            d *= scale;  // SIGMA/RELAX: Å⁻¹ -> nm⁻¹ (×10); PLATEAUS: 1 (manuscript back-compat)
        }
        KEnRef<Real>::saturate(d, maxForceSquared_, numOmpThreads);
        adapter.addLocalModelDerivatives(i, d);
    }

    return energy;
}

template class KEnRefDriver<double>;
template class KEnRefDriver<float>;

} // namespace kenref
