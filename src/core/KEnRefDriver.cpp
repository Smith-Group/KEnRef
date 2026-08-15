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

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace kenref {

namespace {

/*! \brief Refuse to refine a structure that is split across a periodic boundary.
 *
 * KEnRef restrains a WHOLE molecule: the Kabsch superposition and the pair distances are global and
 * have no cutoff, so a molecule torn across a periodic boundary does not fail loudly -- it quietly
 * yields a different energy and different forces. That is the worst failure mode this project has,
 * and it is not hypothetical: every committed GB3 test set carries such a structure, and every serial
 * run against them was refining a torn protein (see res/PBC-BROKEN.md).
 *
 * This check lives in the shared driver deliberately, so BOTH engines are covered by one
 * implementation. The GROMACS side repairs the coordinates from the topology before the driver ever
 * sees them, so it passes. The PLUMED side has no topology and cannot currently repair them, so it
 * stops here rather than continuing with wrong numbers.
 *
 * The test: no atom may sit further than half the smallest box vector from the set's centroid. A
 * molecule that genuinely exceeds that cannot be treated by minimum-image reasoning at all, so the
 * limit is meaningful in its own right rather than an arbitrary tolerance.
 *
 * Note that PLUMED's own ActionAtomistic::makeWhole() does NOT solve this in general: it walks a
 * Euclidean spanning tree built from the MOLINFO *reference* coordinates, so when the reference is
 * broken in the same way as the frame -- the usual case, since both come from the same trajectory --
 * the tree links the wrapped fragment to the main body by an edge SHORTER than half the box, and
 * minimum-image then concludes there is nothing to do. Measured on the committed GB3 reference: the
 * shortest such edge is 29.17 A against a 30.61 A half-box, so the break survives by 1.44 A. */
template<typename Real>
void refuseIfSplitByPeriodicity(const CoordsMatrixType<Real>& x,
                                const Eigen::Matrix<Real, 3, 3>& box_nm, const char* what) {
    if (x.rows() < 2)
        return;
    /* UNITS: the adapters hand over coordinates already scaled to Angstrom but leave the box in nm,
     * which is why restoreNoJump() takes a toAngstrom flag and scales the box itself. Match that
     * convention here, or the comparison is wrong by a factor of ten. */
    const Eigen::Matrix<Real, 3, 3> box = box_nm * Real(10);
    Real smallestBoxVector = std::numeric_limits<Real>::max();
    for (int k = 0; k < 3; ++k)
        if (box(k, k) > Real(0))
            smallestBoxVector = std::min(smallestBoxVector, box(k, k));
    if (smallestBoxVector == std::numeric_limits<Real>::max())
        return;                                   // no periodicity declared; nothing to check
    const Real limit = smallestBoxVector / Real(2);

    const Eigen::RowVector3<Real> centroid = x.colwise().mean();
    Eigen::Index worstRow = 0;
    Real worst = Real(0);
    for (Eigen::Index i = 0; i < x.rows(); ++i) {
        const Real d = (x.row(i) - centroid).norm();
        if (d > worst) {
            worst = d;
            worstRow = i;
        }
    }
    if (worst > limit) {
        std::ostringstream msg;
        msg << "KEnRef refuses to continue: the '" << what << "' atoms are spread over "
            << worst << " (units of the input) from their centroid, which exceeds half the smallest "
            << "box vector (" << limit << "). Atom " << worstRow << " is the furthest.\n"
            << "This almost always means the restrained molecule is SPLIT across a periodic boundary, "
            << "so the superposition and the pair distances would be computed on a torn structure -- "
            << "silently giving wrong energies and forces rather than failing.\n"
            << "Fix the input rather than the symptom: make the molecule whole before refining it "
            << "(e.g. `gmx trjconv -pbc mol` on the structure AND on the reference), and check the "
            << "reference file too -- it is commonly broken in the same way as the frames. The "
            << "GROMACS engine repairs this automatically from the topology; the PLUMED engine has no "
            << "topology available and cannot.";
        throw std::runtime_error(msg.str());
    }
}

} // namespace

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

        /* Before anything is computed from these coordinates, establish that they describe a whole
         * molecule. Checked every step, not only the first: a molecule can drift across a boundary
         * mid-run, and by then nothing downstream would notice. */
        if (checkPeriodicSplit_) {
            refuseIfSplitByPeriodicity(guideX, box, "guide");
            refuseIfSplitByPeriodicity(subX, box, "sub");
        }

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
