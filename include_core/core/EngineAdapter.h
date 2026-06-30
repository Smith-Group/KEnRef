/*
 * EngineAdapter.h
 *
 *  ROLE: the ONE abstract interface that every MD engine (GROMACS, PLUMED, the offline tools)
 *  implements so that kenref_core can run its per-step refinement without any engine-specific code.
 *  It is the single core<->engine boundary: the driver reaches the engine only through this.
 *
 *  The single core<->engine boundary for the KEnRef restructure. KEnRefDriver (in kenref_core)
 *  owns the whole per-step pipeline; it calls back through this abstract interface for the handful
 *  of operations that are genuinely engine-specific. Each consumer implements EngineAdapter ONCE,
 *  generically over all models:
 *    - KEnRefForceProvider  (GROMACS, live MD : one MPI replica per process)
 *    - KEnRefBias           (PLUMED,  live MD : one MPI replica per process)
 *    - TrajectoryEngineAdapter (offline tools : one process, N trajectory files)
 *
 *  No engine types appear in these signatures (only Eigen / std), so kenref_core stays free of
 *  GROMACS and PLUMED headers and keeps building/validating standalone.
 *
 *  Multi-model assembly. KEnRef computes one energy over the coordinates of ALL models. Here a "model"
 *  (as in numModelsInThisProcess/numModelsTotal) is a SYSTEM COPY — one MD replica (live) or one input
 *  trajectory (offline) — NOT an EnergyModel (the restraint type SIGMA/PLATEAUS/RELAX); the restraint
 *  couples all the system-copy models. Each process owns some LOCAL models
 *  (live: exactly 1; offline: all N); the master (simulationIndex()==0) assembles every model's
 *  fitted coordinates, computes, then the per-model derivatives are distributed back. The driver
 *  performs the engine-agnostic math (no-jump correction, Kabsch fit, inverse-fit, unit scale,
 *  saturation); the adapter performs only engine I/O and the gather/scatter (or, single-process,
 *  the trivial move).
 */

#ifndef KENREF_ENGINEADAPTER_H_
#define KENREF_ENGINEADAPTER_H_

#include <string>
#include <vector>
#include <optional>

#include <Eigen/Dense>

#include "KEnRef.h" // CoordsMatrixType<>, KEnRef_Real_t

namespace kenref {

template<typename Real>
class EngineAdapter {
public:
    virtual ~EngineAdapter() = default;

    // ---- Parameters (schema-driven) -------------------------------------------------
    // Raw string the user supplied for `key`, or nullopt if unset. The driver/model converts to
    // the declared ParamType. PLUMED implements this with parse(); GROMACS with CLI11/Settings.
    [[nodiscard]] virtual std::optional<std::string> getRawParam(const std::string& key) const = 0;

    // Number of OpenMP threads to use for THIS step's kernels. Queried by the driver once per step
    // because the engine may change it at runtime (e.g. GROMACS gmx_omp_nthreads_get per step). The
    // KEnRef convention is that 0 means "use all available"; an engine may instead return its own
    // managed thread count.
    [[nodiscard]] virtual int numOmpThreads() const = 0;

    // ---- Per-step, per model held by THIS process (engine-native, Angstrom) ---------
    // Live MD: exactly one model in this process (this replica). Offline tools: one per input
    // trajectory (all N live in the single process).
    [[nodiscard]] virtual int numModelsInThisProcess() const = 0;

    // Fill this local model's CURRENT guide-atom coords (for alignment), restrained sub-atom coords,
    // and the simulation box (row vectors, for periodic no-jump correction). Implementations resize
    // `guideX` to [nGuide x 3] and `subX` to [nSub x 3].
    virtual void getLocalModelX(int localModel,
                                CoordsMatrixType<Real>& guideX,
                                CoordsMatrixType<Real>& subX,
                                Eigen::Matrix<Real, 3, 3>& box) const = 0;

    // Add this local model's derivatives ([nSub x 3], already inverse-fitted, unit-scaled and
    // saturated by the driver) to the engine's forces. Offline tools that only report energy may
    // leave this a no-op.
    //
    // DELIBERATELY pure (no default): unlike the single-process/no-MPI methods, this is intentionally
    // NOT given a no-op default in DefaultAdapter — a force-applying engine must not be able to
    // silently forget to apply forces. Energy/metric-only consumers opt into the no-op explicitly.
    virtual void addLocalModelDerivatives(int localModel,
                                          const CoordsMatrixType<Real>& derivs) = 0;

    // ---- Cross-process assembly (master collects every model) -----------------------
    // Total models across all processes; master is simulationIndex()==0.
    [[nodiscard]] virtual int numModelsTotal()  const = 0;
    [[nodiscard]] virtual int simulationIndex() const = 0;

    // Gather every process's local fitted sub-coords into `all` (sized numModelsTotal) on the
    // master. Single-process engines (single-sim live, or offline) simply move localFitted -> all.
    virtual void gatherFittedSubAtomsX(const std::vector<CoordsMatrixType<Real>>& localFitted,
                                       std::vector<CoordsMatrixType<Real>>& all) const = 0;

    // Scatter the master's per-model derivatives back so each process receives its local models'
    // blocks. Single-process engines simply move allPerModel -> localPerModel.
    virtual void scatterModelDerivatives(const std::vector<CoordsMatrixType<Real>>& allPerModel,
                                         std::vector<CoordsMatrixType<Real>>& localPerModel) const = 0;
};

} // namespace kenref

#endif // KENREF_ENGINEADAPTER_H_
