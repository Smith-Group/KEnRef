/*
 * KEnRefDriver.h
 *
 *  ROLE: the single, engine-agnostic per-step refinement pipeline that turns engine coordinates into
 *  restraint forces. It owns the user-selected EnergyModel and is invoked once per MD step by each
 *  engine consumer; it is the shared replacement for the per-engine calculateForces/calculate bodies.
 *
 *  The shared per-step refinement pipeline, owned by kenref_core so every engine consumer runs the
 *  SAME engine-agnostic math. The driver holds the selected EnergyModel and the general params, and
 *  on each MD step it: fetches this process's model coordinates through the EngineAdapter, removes
 *  periodic jumps, fits each to the reference (Kabsch), gathers the fitted coords to the master,
 *  computes energy + per-model derivatives, scatters them back, inverse-fits + unit-scales +
 *  saturates, and applies the forces — all via the adapter's ~handful of engine-specific callbacks.
 *
 *  The model is prepared (buildCache + finalizeIndexing) by the consumer using its engine-specific
 *  atom indexing BEFORE constructing the driver; the driver only needs the centered guide-atom
 *  reference coordinates for alignment and the general k/n/maxForce params.
 */

#ifndef KENREF_KENREFDRIVER_H_
#define KENREF_KENREFDRIVER_H_

#include <memory>
#include <vector>

#include <Eigen/Geometry>

#include "KEnRef.h"
#include "EnergyModel.h"
#include "EngineAdapter.h"

namespace kenref {

template<typename Real>
class KEnRefDriver {
public:
    // numOmpThreads is NOT cached here — it is queried from the EngineAdapter once per step()
    // (the engine may change it at runtime).
    KEnRefDriver(std::unique_ptr<EnergyModel<Real>> model,
                 Real k, Real n, Real maxForceSquared,
                 CoordsMatrixType<Real> guideAtomsReferenceCoordsCentered)
        : model_(std::move(model)), k_(k), n_(n), maxForceSquared_(maxForceSquared),
          referenceGuideAtomsCoordsCentered_(std::move(guideAtomsReferenceCoordsCentered)) {}

    /*! \brief Enable the periodic-split refusal (see the check in KEnRefDriver.cpp).
     *
     * OFF by default, and deliberately opt-in rather than automatic: the check compares the spread of
     * the restrained atoms against the box, so it is only meaningful when the box really is the
     * simulation's periodic box. A live engine knows that; a unit test feeding synthetic coordinates
     * with a placeholder box does not, and would be failed spuriously by it.
     *
     * Both live engines turn it on -- GROMACS in initParamsAtSetup(), PLUMED in the KEnRefBias
     * constructor -- so any real refinement is covered. */
    void enablePeriodicSplitCheck(bool enable = true) { checkPeriodicSplit_ = enable; }

    // Run one MD step. Returns the total energy (valid on the master, simulationIndex()==0; 0 elsewhere).
    Real step(EngineAdapter<Real>& adapter, bool printStatistics = false);

    [[nodiscard]] EnergyModel<Real>& model() { return *model_; }
    [[nodiscard]] const EnergyModel<Real>& model() const { return *model_; }

private:
    std::unique_ptr<EnergyModel<Real>> model_;
    Real k_;
    Real n_;
    Real maxForceSquared_;
    CoordsMatrixType<Real> referenceGuideAtomsCoordsCentered_;
    //! Whether to refuse structures split across a periodic boundary; see enablePeriodicSplitCheck().
    bool checkPeriodicSplit_ = false;

    // Per-step scratch, reused across steps (one entry per model held by THIS process).
    std::vector<Eigen::Transform<Real, 3, Eigen::Affine>> affines_;
    std::vector<CoordsMatrixType<Real>> localFitted_;
    std::vector<CoordsMatrixType<Real>> allModels_;       // gathered on the master
    std::vector<CoordsMatrixType<Real>> perModelDerivs_;  // master's per-model derivatives
    std::vector<CoordsMatrixType<Real>> localDerivs_;     // scattered back to this process
    std::vector<CoordsMatrixType<Real>> lastGuide_;       // previous-frame refs for no-jump
    std::vector<CoordsMatrixType<Real>> lastSub_;
    bool firstStep_ = true;
};

} // namespace kenref

#endif // KENREF_KENREFDRIVER_H_
