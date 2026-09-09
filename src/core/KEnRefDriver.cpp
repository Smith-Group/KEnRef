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
#include <vector>

namespace kenref {

namespace detail {

/*! \brief Refuse to refine a structure that is split across a periodic boundary.
 *
 * KEnRef restrains a WHOLE molecule: the Kabsch superposition and the pair distances are global and
 * have no cutoff, so a molecule torn across a periodic boundary does not fail loudly -- it quietly
 * yields a different energy and different forces. That is the worst failure mode this project has,
 * and it is not hypothetical: every committed GB3 test set carries such a structure, and every serial
 * run against them was refining a torn protein (see res/PBC-BROKEN.md).
 *
 * This check lives in the shared driver deliberately, so BOTH engines are covered by one
 * implementation. Both repair first and are checked second: the GROMACS side rebuilds the coordinates
 * from the topology's bond graph, the PLUMED side from a spanning tree over the reference structure.
 * This is the backstop that catches a repair which did not work -- most importantly a PLUMED reference
 * that is itself broken, in which case the tree links the wrapped fragment by a long edge and the
 * repair silently does nothing.
 *
 * WHAT IS ACTUALLY TESTED: a TEAR, not a size. The question is whether re-imaging some atoms would
 * make this set materially more compact -- because that, and only that, is what a periodic split is.
 *
 * For each lattice direction the atoms are projected onto their fractional coordinate. Two extents are
 * then compared along that direction:
 *   * the extent AS DELIVERED, max(s) - min(s) over the raw fractional coordinates; and
 *   * the smallest extent ANY choice of periodic images could give, which is 1 minus the largest
 *     CYCLIC gap between neighbouring atoms on the unit circle.
 * For a set that is whole these are equal, exactly: with no atom wrapped, the largest cyclic gap is
 * the empty space OUTSIDE the molecule, and 1 minus it is the molecule's own extent. They differ only
 * when the largest gap is an INTERIOR void -- i.e. when the set is split, with atoms piled at both
 * ends of the box and nothing in the middle. The difference is then how much a correct re-imaging
 * would recover, which is the size of the tear.
 *
 * WHY NOT "no atom further than half the smallest box vector from the centroid", which this replaces:
 * that is a test of SIZE, and it fails on structures that are perfectly whole. Measured on the
 * committed GB3 sets after the bond-graph repair, the restrained atoms reach 21.657 A against a
 * 21.6561 A half-box -- over by 0.004% -- so a 500-step run aborted at about step 120 on a molecule
 * that was never torn. Extent is also not the property that matters: once a set is whole, nothing
 * downstream uses minimum-image reasoning over the whole molecule. The pair distances have no cutoff
 * and no PBC at all, the Kabsch fit is global, and restoreNoJump images each atom against ITS OWN
 * position in the previous frame -- a per-step displacement, not a molecular radius. Half the box
 * bounds the reach of minimum-image reasoning, and both repairs deliberately avoid needing it: they
 * walk short edges (a bond, or a nearest-neighbour link in the reference), never a molecular radius.
 *
 * Note that PLUMED's own ActionAtomistic::makeWhole() does NOT solve this in general: it walks a
 * Euclidean spanning tree built from the MOLINFO *reference* coordinates, so when the reference is
 * broken in the same way as the frame -- the usual case, since both come from the same trajectory --
 * the tree links the wrapped fragment to the main body by an edge SHORTER than half the box, and
 * minimum-image then concludes there is nothing to do. Measured on the un-repaired GB3 reference: the
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
    for (int k = 0; k < 3; ++k)
        if (!(box(k, k) > Real(0)))
            return;                                // no periodicity declared; nothing to check

    /* Row k of the box IS lattice vector k (the convention restoreNoJump uses: box(m, d) is component
     * d of vector m), so cartesian = fractional * box and fractional = cartesian * box^-1. Working in
     * fractional coordinates is what makes this correct for a triclinic cell: whatever the cell shape,
     * a periodic image differs by a whole number in exactly one fractional component. */
    const Eigen::Matrix<Real, 3, 3> toFractional = box.inverse();
    const CoordsMatrixType<Real> s = x * toFractional;

    /* The tolerance below which a difference between the two extents is not called a tear. For a whole
     * set the difference is identically zero, so this only has to clear round-off; it is set at 2% of
     * the box vector so that a marginal case -- where re-imaging would barely help because the set
     * really does span the cell -- is not turned into a hard stop over an ambiguity. A genuine tear is
     * nothing like this small: on the committed GB3 sets it recovers 28.6 A of a 61.2 A box. */
    constexpr Real kTearToleranceFraction = Real(0.02);

    std::vector<Real> wrapped(static_cast<std::size_t>(s.rows()));
    for (int k = 0; k < 3; ++k) {
        Real lo = s(0, k), hi = s(0, k);
        for (Eigen::Index i = 0; i < s.rows(); ++i) {
            lo = std::min(lo, s(i, k));
            hi = std::max(hi, s(i, k));
            // Fold onto [0,1). std::fmod keeps the sign, so a negative coordinate needs one more turn.
            Real f = std::fmod(s(i, k), Real(1));
            if (f < Real(0))
                f += Real(1);
            wrapped[static_cast<std::size_t>(i)] = f;
        }
        const Real deliveredExtent = hi - lo;

        std::sort(wrapped.begin(), wrapped.end());
        // Largest gap between neighbours ON THE CIRCLE -- the last one wraps past 1 back to the first.
        Real largestGap = wrapped.front() + Real(1) - wrapped.back();
        for (std::size_t i = 1; i < wrapped.size(); ++i)
            largestGap = std::max(largestGap, wrapped[i] - wrapped[i - 1]);
        const Real bestPossibleExtent = Real(1) - largestGap;

        const Real recoverable = deliveredExtent - bestPossibleExtent;
        if (recoverable <= kTearToleranceFraction)
            continue;

        const Real boxVectorLength = box.row(k).norm();
        std::ostringstream msg;
        msg << "KEnRef refuses to continue: the '" << what << "' atoms are SPLIT across the periodic "
            << "boundary along box vector " << k << ". They span " << deliveredExtent * boxVectorLength
            << " (units of the input) along it, but re-imaging them correctly would fit them into "
            << bestPossibleExtent * boxVectorLength << " -- there is a void of "
            << largestGap * boxVectorLength << " through the middle of the set, with atoms piled at "
            << "both ends of a " << boxVectorLength << " box vector.\n"
            << "The superposition and the pair distances would be computed on a torn structure -- "
            << "silently giving wrong energies and forces rather than failing.\n"
            << "Both engines repair this before reaching here, so this means the repair did not work. "
            << "The usual cause is a REFERENCE structure that is itself split: the PLUMED side builds "
            << "its spanning tree from the reference, and a reference broken the same way as the "
            << "frames links the wrapped fragment by an edge shorter than half the box, so "
            << "minimum-image leaves it exactly where it was. Make the reference whole -- "
            << "`gmx trjconv -pbc whole` against a topology -- and check the frames too. On the "
            << "GROMACS side the repair comes from the topology instead, so this points at the "
            << "topology being unavailable (or KENREF_NO_MAKEWHOLE being set).";
        throw std::runtime_error(msg.str());
    }
}

template void refuseIfSplitByPeriodicity<double>(const CoordsMatrixType<double>&,
                                                 const Eigen::Matrix<double, 3, 3>&, const char*);
template void refuseIfSplitByPeriodicity<float>(const CoordsMatrixType<float>&,
                                                const Eigen::Matrix<float, 3, 3>&, const char*);

} // namespace detail

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
            detail::refuseIfSplitByPeriodicity(guideX, box, "guide");
            detail::refuseIfSplitByPeriodicity(subX, box, "sub");
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
