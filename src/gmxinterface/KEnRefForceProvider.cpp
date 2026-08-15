/*
 * KEnRefForceProvider.cpp
 *
 *      Author: amr
 */

#include <iostream>
//#include <typeinfo>
#include <cmath>
#include <algorithm>    // std::min, in the make-whole plausibility check
#include <cstdlib>      // std::getenv, for the KENREF_DD_SELFCHECK switch
#include <limits>       // std::numeric_limits, same
#include <filesystem>   // std::filesystem::absolute — pulled in transitively by libc++, NOT by libstdc++
#include <memory>
#include <vector>
#include <Eigen/Core>
#include <utility>
#include<unistd.h>

#include "gromacs/mdtypes/commrec.h"
#include "gromacs/version.h"
#include "mpi.h"

#include "gromacs/pbcutil/pbc.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/domdec/ga2la.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/arrayref.h"
#include "core/KEnRef.h"
#include "core/kabsch.h"
#include "core/IoUtils.h"
#include "core/ModelRegistry.h"
#include "core/buildModelIndexing.h"
#include "gmxinterface/KEnRefForceProvider.h"

#include <gromacs/mdrunutility/mdmodulesnotifiers.h>

#include "gmxinterface/KEnRefMDModule.h"
#include "gmxinterface/gmxwrapper.h"

#define VERBOSE false

static constexpr auto singleStr = "single";

// MPI datatype matching KEnRef_Real_t (for the per-replica gather/scatter).
static const MPI_Datatype KENREF_MPI_REAL = std::is_same_v<KEnRef_Real_t, float> ? MPI_FLOAT : MPI_DOUBLE;

/* GROMACS 2026 replaced ForceProviderInput's `const t_commrec& cr_` with an MpiComm reference plus a
 * (possibly null) gmx_domdec_t pointer, and dropped t_commrec::rankInDefaultCommunicator entirely.
 * These accessors are the only place that difference is expressed, so the body of calculateForces()
 * and friends stays version-agnostic. Each one defers to GROMACS's own predicate rather than
 * reimplementing it, so the semantics cannot drift between releases. */
#if GMX_VERSION >= 20260000
#    define KENREF_GMX_HAS_MPICOMM_FORCEPROVIDERINPUT 1
#else
#    define KENREF_GMX_HAS_MPICOMM_FORCEPROVIDERINPUT 0
#endif

namespace {

//! The domain-decomposition object for this step, or nullptr when DD is not in use.
inline const gmx_domdec_t *kenrefDd(const gmx::ForceProviderInput &in) {
#if KENREF_GMX_HAS_MPICOMM_FORCEPROVIDERINPUT
    return in.dd_;
#else
    return in.cr_.dd;
#endif
}

//! The communicator spanning this rank's group (PP or PME).
inline MPI_Comm kenrefGroupComm(const gmx::ForceProviderInput &in) {
#if KENREF_GMX_HAS_MPICOMM_FORCEPROVIDERINPUT
    return in.mpiComm_.comm();
#else
    return in.cr_.mpi_comm_mygroup;
#endif
}

/*! \brief A rank id for diagnostics only.
 *
 * Up to 2025 this reported the rank in the default (whole-run) communicator. That member is gone in
 * 2026 and ForceProviderInput exposes only the group communicator, so on 2026 this is the rank
 * within the group. It is printed for debugging and is not used in any calculation. */
inline int kenrefDiagnosticRank(const gmx::ForceProviderInput &in) {
#if KENREF_GMX_HAS_MPICOMM_FORCEPROVIDERINPUT
    return in.mpiComm_.rank();
#else
    return in.cr_.rankInDefaultCommunicator;
#endif
}

//! Whether atoms are ordered by domain decomposition (so ga2la lookups are required).
inline bool kenrefHaveDDAtomOrdering(const gmx::ForceProviderInput &in) {
#if KENREF_GMX_HAS_MPICOMM_FORCEPROVIDERINPUT
    return haveDDAtomOrdering(in.dd_);
#else
    return haveDDAtomOrdering(in.cr_);
#endif
}

//! Whether there is actual particle-particle domain decomposition.
inline bool kenrefHavePPDomainDecomposition(const gmx::ForceProviderInput &in) {
#if KENREF_GMX_HAS_MPICOMM_FORCEPROVIDERINPUT
    return havePPDomainDecomposition(in.dd_);
#else
    return havePPDomainDecomposition(&in.cr_);
#endif
}

/*! \brief Number of ranks this simulation (replica) is spread over.
 *
 * This is the INTRA-simulation communicator, so it counts the ranks sharing one replica's atoms --
 * i.e. the domain-decomposition width. It is 1 in the only configuration KEnRef supports today, which
 * is what makes every reduction below a no-op there. */
inline int kenrefSimCommSize(const gmx::ForceProviderInput &in) {
    const MPI_Comm comm = kenrefGroupComm(in);
    if (comm == MPI_COMM_NULL)
        return 1;
    int size = 1;
    MPI_Comm_size(comm, &size);
    return size;
}

/*! \brief Whether this rank is the main rank of ITS OWN simulation.
 *
 * Not to be confused with "rank 0 of the world": with `-multidir` every replica has its own main rank.
 * This is the rank that may touch `mainRanksComm_`, because that communicator is MPI_COMM_NULL
 * everywhere else (gmx_multisim_t documents it as "valid only on main ranks"). */
inline bool kenrefIsSimMainRank(const gmx::ForceProviderInput &in) {
    const MPI_Comm comm = kenrefGroupComm(in);
    if (comm == MPI_COMM_NULL)
        return true;              // no intra-simulation communicator ⇒ this rank is alone ⇒ it is main
    int rank = 0;
    MPI_Comm_rank(comm, &rank);
    return rank == 0;
}

/*! \brief Sum a coordinate block across the ranks of one simulation, in place.
 *
 * The zero-fill convention: every rank leaves the rows it does not own at zero and writes only its own,
 * so the sum reconstructs the complete, correctly-ordered block on every rank. This is EXACT, not
 * merely accurate: each row has exactly one non-zero contributor (ga2la::findHome resolves an atom on
 * its home rank only), and `0.0 + x == x` in IEEE-754 for any finite x, with addition commutative --
 * so the result does not depend on the reduction order or the rank count.
 *
 * A size-1 communicator short-circuits, which keeps the supported one-rank-per-replica configuration
 * bit-identical to the code that existed before domain-decomposition support. */
inline void kenrefSumOverSimRanks(const gmx::ForceProviderInput &in, KEnRef_Real_t *data, std::size_t count) {
    const MPI_Comm comm = kenrefGroupComm(in);
    if (comm == MPI_COMM_NULL || count == 0)
        return;
    int size = 1;
    MPI_Comm_size(comm, &size);
    if (size <= 1)
        return;                   // nothing to combine; also the supported configuration
    MPI_Allreduce(MPI_IN_PLACE, data, static_cast<int>(count), KENREF_MPI_REAL, MPI_SUM, comm);
}

/*! \brief Broadcast a block from a simulation's main rank to the rest of that simulation.
 *
 * Needed because only the main rank takes part in the cross-replica exchange (mainRanksComm_ is
 * MPI_COMM_NULL elsewhere), yet every rank has to apply forces to the atoms it owns. No-op on a
 * size-1 communicator, so the supported configuration is untouched. */
inline void kenrefBcastFromSimMain(const gmx::ForceProviderInput &in, KEnRef_Real_t *data, std::size_t count) {
    const MPI_Comm comm = kenrefGroupComm(in);
    if (comm == MPI_COMM_NULL || count == 0)
        return;
    int size = 1;
    MPI_Comm_size(comm, &size);
    if (size <= 1)
        return;
    MPI_Bcast(data, static_cast<int>(count), KENREF_MPI_REAL, 0, comm);
}


/*! \brief Whether the domain-decomposition self-check is switched on (KENREF_DD_SELFCHECK=1).
 *
 * Off by default: the check costs an extra allreduce per atom set per step. Read once. */
inline bool kenrefDdSelfCheckEnabled() {
    static const bool enabled = [] {
        const char *e = std::getenv("KENREF_DD_SELFCHECK");
        return e != nullptr && e[0] != '\0' && e[0] != '0';
    }();
    return enabled;
}

/*! \brief Assert that every row of an atom set was written by EXACTLY ONE rank.
 *
 * This is the property the whole zero-fill design rests on. `0.0 + x == x` exactly in IEEE-754 and
 * addition is commutative, so summing the per-rank buffers reproduces the full set bit-for-bit and
 * independently of the rank count -- but ONLY if each row has exactly one non-zero contributor. Two
 * ranks claiming the same atom would double it; none claiming it would silently leave it at the
 * origin. Neither shows up as a crash, and a wrong coordinate merely biases the restraint, so this
 * is precisely the kind of failure that has to be tested for rather than waited for.
 *
 * \p ownerCount holds 1 for each row this rank wrote; the same reduction used for the coordinates
 * turns it into a global writer count per row.
 *
 * Collective over the simulation's ranks: every rank calls it, from the same place, every step. */
inline void kenrefCheckExactlyOneWriter(const gmx::ForceProviderInput &in,
                                        std::vector<KEnRef_Real_t> &ownerCount,
                                        const char *what, int64_t step) {
    kenrefSumOverSimRanks(in, ownerCount.data(), ownerCount.size());
    std::size_t unowned = 0;
    std::size_t contested = 0;
    std::size_t firstBad = 0;
    bool haveBad = false;
    for (std::size_t i = 0; i < ownerCount.size(); ++i) {
        const KEnRef_Real_t c = ownerCount[i];
        if (c == KEnRef_Real_t(0))
            ++unowned;
        else if (c != KEnRef_Real_t(1))
            ++contested;
        else
            continue;
        if (!haveBad) {
            firstBad = i;
            haveBad = true;
        }
    }
    if (unowned != 0 || contested != 0) {
        gmx_fatal(FARGS,
                  "KEnRef domain-decomposition self-check FAILED for the '%s' atom set at step %lld: "
                  "of %zu rows, %zu had NO owning rank and %zu were claimed by more than one "
                  "(first offending row: %zu). The zero-fill reduction is exact only when every row "
                  "is written exactly once, so the coordinates this step are wrong.",
                  what, static_cast<long long>(step), ownerCount.size(), unowned, contested, firstBad);
    }
}

} // namespace

KEnRefForceProvider::KEnRefForceProvider() = default;

KEnRefForceProvider::~KEnRefForceProvider() = default;

[[maybe_unused]] void KEnRefForceProvider::setSimulationContext(gmx::SimulationContext *simulationContext) {
    this->simulationContext_ = simulationContext;
}

void KEnRefForceProvider::setGuideAtom0Indices(std::shared_ptr<std::vector<int> const> targetAtoms0Indices) {
    this->guideAtom0Indices_ = std::move(targetAtoms0Indices);
}

void KEnRefForceProvider::setGuideAtomsReferenceCoords(
        std::shared_ptr<const CoordsMatrixType<KEnRef_Real_t>> &guideAtomsReferenceCoords) {
    this->guideAtomsReferenceCoords_ = std::move(guideAtomsReferenceCoords);
    //keep another cashed version after moving its Center of Mass (COM) to the origin for faster processing
    this->guideAtomsReferenceCoordsCentered_ = std::make_shared<const CoordsMatrixType<KEnRef_Real_t>>(
            Kabsch_Umeyama<KEnRef_Real_t>::translateCenterOfMassToOrigin(*this->guideAtomsReferenceCoords_));
    //TODO do we need a guideAtomsReferenceCoordsCOM_ ?
}


void KEnRefForceProvider::calculateForces(const gmx::ForceProviderInput &forceProviderInput,
                                          gmx::ForceProviderOutput *forceProviderOutput) {
    auto begin = std::chrono::high_resolution_clock::now();
#if VERBOSE
    std::cout << "calculateForces() called" << std::endl;
#endif
    const auto homenr = forceProviderInput.homenr_; // total number of atoms in the system (or domain dec ?)
    GMX_ASSERT(homenr >= 0, "number of home atoms must be non-negative.");

    const auto &step = forceProviderInput.step_;

    // Stash the per-step GROMACS state + MPI info so the EngineAdapter callbacks (invoked by
    // driver_->step(*this) below) can reach the current input/output and the multi-simulation info.
    isMultiSimulation_ = this->simulationContext_->multiSimulation_ != nullptr;
    numSimulations_ = isMultiSimulation_ ? this->simulationContext_->multiSimulation_->numSimulations_ : 1;
    simulationIndex_ = isMultiSimulation_ ? this->simulationContext_->multiSimulation_->simulationIndex_ : 0;
    mainRanksComm_ = isMultiSimulation_ ? this->simulationContext_->multiSimulation_->mainRanksComm_ : MPI_COMM_NULL;
    currentInput_ = &forceProviderInput;
    currentOutput_ = forceProviderOutput;
    /* Refresh every step rather than caching at setup: it is one MPI_Comm_rank on an intra-simulation
     * communicator, and a cached value could go stale if the communicator were ever rebuilt. Always
     * true at one rank per replica. */
    isSimMainRank_ = kenrefIsSimMainRank(forceProviderInput);

    // (the stdout redirection and the whole of the one-time setup now happen in initParamsAtSetup(),
    //  called from KEnRefMDModule::initForceProviders -- see the comment on its declaration)

    /* Diagnostics, one line per replica rather than one per domain: under domain decomposition every
     * rank would otherwise print this every 10 steps. */
    if (step % 10 == 0 && isReplicaMasterRank())
        std::cout
                << "--> numSimulations " << numSimulations_ << "\n"
                << "--> rank " << kenrefDiagnosticRank(forceProviderInput) << " " << (isMultiSimulation_ ? simulationIndex_ : -1) << "\n"
                << "--> simulationIndex " << simulationIndex_ << "\tstep " << step << std::endl;

    if (!paramsInitialized) {
        volatile bool holdToDebug = Settings::debug;
        while (/*simulationIndex > 0 &&*/ holdToDebug) {
            sleep(1);
        }
    }

    if (!paramsInitialized) {
        GMX_ASSERT(check_box(PbcType::Unset, forceProviderInput.box_) == nullptr, "Invalid box.");

        /* Domain decomposition is supported, but ONLY with bond connectivity available.
         *
         * Under DD GROMACS hands out coordinates wrapped into the box, so a restrained molecule
         * crossing a periodic boundary arrives broken and the restraint would be computed on a
         * distorted structure -- silently. Repairing that requires the topology: walking bonds is the
         * only unambiguous way to choose an image for a molecule that can be wider than half the box.
         * Without the topology there is no way to repair it, so refuse rather than return wrong
         * numbers. This mirrors GROMACS's own stance for whole-molecule-requiring features. */
        if (kenrefHavePPDomainDecomposition(forceProviderInput) && !moleculeGraph_.isBuilt()) {
            gmx_fatal(FARGS,
                      "KEnRef cannot run under domain decomposition without the molecular topology: "
                      "this simulation has %d PP ranks, and no bond connectivity is available. Under "
                      "domain decomposition coordinates are wrapped into the box, so a restrained "
                      "molecule crossing a periodic boundary must be repaired using its bonds before "
                      "the restraint is computed. Run ONE RANK PER SIMULATION -- e.g. `mpirun -np <N> "
                      "... -multidir <N directories>` -- or reduce the rank count so no simulation is "
                      "decomposed. Threads are unaffected: -ntomp is free to use.",
                      kenrefDd(forceProviderInput)->nnodes);
        }

        if (isReplicaMasterRank()) {
            const gmx_domdec_t *dd = kenrefDd(forceProviderInput);
            std::cout << "Number of atoms = " << homenr << std::endl;
            std::cout << "havePPDomainDecomposition: " << kenrefHavePPDomainDecomposition(forceProviderInput) << std::endl;
            std::cout << "haveDDAtomOrdering: " << kenrefHaveDDAtomOrdering(forceProviderInput) << std::endl;
            // dd is null when domain decomposition is not in use at all, so do not dereference blindly.
            std::cout << "dd->nnodes: " << (dd != nullptr ? dd->nnodes : 0) << std::endl;
            std::cout << "ranks in this simulation: " << kenrefSimCommSize(forceProviderInput) << std::endl;
        }
        paramsInitialized = true;
    }

    // Run the shared per-step refinement pipeline: getLocalModelX -> no-jump -> Kabsch fit -> gather
    // -> model.compute -> scatter -> inverse-fit -> unit-scale -> saturate -> addLocalModelDerivatives.
    const KEnRef_Real_t energy = driver_->step(*this, /*printStatistics*/ (step % 10 == 0));

    /* Report from ONE rank only. `energy` is returned by the driver on the rank that actually ran the
     * model; every other rank gets 0, so gating on simulationIndex_ alone would print a spurious
     * "Energy: 0" once per domain. */
    if ((step % 10 == 0) && isEnsembleMasterRank())
        std::cout << "Step: " << step << " Energy: " << energy << std::endl;

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    this->calculateForces_time += elapsed.count();
    if (!(step % 10) && isEnsembleMasterRank()) {
        printf("This iteration (%ld): %.5f seconds. All walltime %.3f seconds\n", step,
               static_cast<double>(elapsed.count()) * 1e-9,
               static_cast<double>(calculateForces_time) * 1e-9);
        std::cout << "=================" << std::endl;
    }

    currentInput_ = nullptr;
    currentOutput_ = nullptr;
}

// ==================================================================================================
//  EngineAdapter implementation (GROMACS, live MD: exactly one model in this process/replica)
// ==================================================================================================

std::optional<std::string> KEnRefForceProvider::getRawParam(const std::string &key) const {
    // Model-tier params come from the GROMACS Settings (set by the CLI/TOML in gmxwrapper).
    if (key == "EXP_DATA_FOLDER")
        return Settings::experimentalDataFolder.empty() ? std::nullopt
                                                        : std::optional<std::string>(Settings::experimentalDataFolder);
    if (key == "EXP_DATA_FILE")
        return Settings::experimentalDataFileName.empty() ? std::nullopt
                                                          : std::optional<std::string>(Settings::experimentalDataFileName);
    if (key == "PROTON_MHZ")
        return std::optional<std::string>(std::to_string(Settings::proton_mhz));
    return std::nullopt;
}

int KEnRefForceProvider::numOmpThreads() const {
    // Queried fresh each step (GROMACS may change its managed thread count at runtime), matching the
    // original which passed gmx_omp_nthreads_get(...) to every kernel call.
    return gmx_omp_nthreads_get(ModuleMultiThread::Default);
}

void KEnRefForceProvider::getLocalModelX(int /*localModel*/, CoordsMatrixType<KEnRef_Real_t> &guideX,
                                         CoordsMatrixType<KEnRef_Real_t> &subX,
                                         Eigen::Matrix<KEnRef_Real_t, 3, 3> &box) const {
    /* Two ways to obtain the restrained coordinates.
     *
     * (a) Via the WHOLE MOLECULE, when connectivity is available. Gather every atom of the molecules
     *     that carry restrained atoms, repair the periodic images by walking the bond spanning tree,
     *     then slice the guide and sub sets out of the result. This is what makes domain decomposition
     *     correct: GROMACS hands out box-wrapped coordinates, so a molecule crossing a boundary arrives
     *     broken, and only connectivity can say which image each atom belongs in. Both sets come out of
     *     ONE repaired coordinate array, so they are necessarily in the same frame -- which the Kabsch
     *     fit requires, since its transform is derived from the guide atoms and applied to the sub ones.
     *
     * (b) Directly, when there is no topology. Correct without domain decomposition, where coordinates
     *     arrive exactly as the tpr holds them.
     *
     * Path (a) is a bit-for-bit no-op on already-whole input -- every bond difference is already the
     * minimum image, so nothing is rewritten -- which is why it can run unconditionally rather than only
     * under DD, and why it does not perturb the serial results. */
    /* Escape hatch, for reproducing pre-repair results and for testing the driver's refusal check.
     * Safe to offer only because that check exists: with the repair off, a torn structure is refused
     * rather than silently refined. Read once. */
    static const bool makeWholeDisabled = [] {
        const char *e = std::getenv("KENREF_NO_MAKEWHOLE");
        return e != nullptr && e[0] != '\0' && e[0] != '0';
    }();

    if (moleculeGraph_.isBuilt() && !makeWholeDisabled) {
        fillOwnedRows(moleculeX_, moleculeAtomSet_, moleculeGraph_.atoms(), *currentInput_, "molecule");
        /* Report at step 0 whether the input was actually broken. This is worth saying out loud: a
         * molecule split across a periodic boundary is invisible in the output otherwise, and it
         * silently changes the restraint. All five committed GB3 test sets are split this way. */
        const KEnRef_Real_t longestBefore = moleculeGraph_.checkBondLengths(moleculeX_);
        moleculeGraph_.makeWhole(moleculeX_, currentInput_->box_);
        if (currentInput_->step_ == 0 && isReplicaMasterRank() && longestBefore > KEnRef_Real_t(0.5)) {
            std::cout << "KEnRef: input molecule was SPLIT across a periodic boundary (longest bond "
                      << longestBefore << " nm); repaired to " << moleculeGraph_.checkBondLengths(moleculeX_)
                      << " nm using the topology's connectivity." << std::endl;
        }

        guideX.resize(static_cast<Eigen::Index>(guideRowInMolecule_.size()), 3);
        for (std::size_t i = 0; i < guideRowInMolecule_.size(); ++i)
            guideX.row(static_cast<Eigen::Index>(i)) = moleculeX_.row(guideRowInMolecule_[i]);
        subX.resize(static_cast<Eigen::Index>(subRowInMolecule_.size()), 3);
        for (std::size_t i = 0; i < subRowInMolecule_.size(); ++i)
            subX.row(static_cast<Eigen::Index>(i)) = moleculeX_.row(subRowInMolecule_[i]);

        if (kenrefDdSelfCheckEnabled()) {
            /* A bond far longer than a chemical bond is the one way this construction fails, and it
             * would otherwise be silent. 0.5 nm is generous for any bond or constraint. */
            const KEnRef_Real_t longest = moleculeGraph_.checkBondLengths(moleculeX_);
            if (longest > KEnRef_Real_t(0.5)) {
                gmx_fatal(FARGS,
                          "KEnRef make-whole left a bond of %g nm, which is far longer than any "
                          "chemical bond. The connectivity and the coordinates disagree, so the "
                          "restrained structure is not correctly re-imaged.",
                          static_cast<double>(longest));
            }
        }
    } else {
        guideX = getGuideAtomsX(*currentInput_, /*toAngstrom*/ false);
        subX.resize(static_cast<Eigen::Index>(this->sub0Id_to_global1Id_->size()), 3);
        fillSubAtomsX(subX, *currentInput_, /*toAngstrom*/ false);
    }

    guideX *= 10;   // nm -> Angstrom, deferred until after make-whole
    subX *= 10;
    // GROMACS box (nm, real[3][3]) -> Eigen 3x3 (raw; the driver scales it to Angstrom via toAngstrom).
    const matrix &b = currentInput_->box_;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            box(i, j) = static_cast<KEnRef_Real_t>(b[i][j]);
}

void KEnRefForceProvider::addLocalModelDerivatives(int /*localModel*/, const CoordsMatrixType<KEnRef_Real_t> &derivs) {
    /* The pre-apply barrier that used to be here is gone. It was a no-op in disguise: every rank
     * reaches this function through scatterModelDerivatives(), which now ends in a broadcast down the
     * simulation's own communicator (or, at one rank per simulation, in nothing at all). A broadcast
     * already orders every rank of the simulation against the root, so the barrier added a full
     * synchronisation per step and bought nothing. Forces are written to disjoint local atoms, so
     * there is nothing to order between ranks here anyway. */
    const auto &force = currentOutput_->forceWithVirial_.force_;
    //Finally, add them to corresponding atoms
    /* Each rank applies forces ONLY to the atoms it owns. There is no reduction here and there must
     * not be one: the force array is indexed by local atom, the rank sets are disjoint, and a rank
     * that does not own an atom has nowhere to put its force. Skipping is therefore the whole fix --
     * a non-owned atom is another rank's responsibility, not a lost contribution.
     *
     * Mirrors fillOwnedRows(): same two ways of deciding ownership, same re-fetch-every-step rule for
     * the ArrayRefs. It has to stay the mirror image, because a row filled from atom A but whose force
     * lands on atom B would be a silent, energy-conserving-looking corruption. */
    if (subAtomSet_.has_value()) {
        const gmx::ArrayRef<const int> localIndex = subAtomSet_->localIndex();
        const gmx::ArrayRef<const int> collectiveIndex = subAtomSet_->collectiveIndex();
        for (std::size_t k = 0; k < localIndex.size(); ++k) {
            const auto row = static_cast<Eigen::Index>(collectiveIndex[k]);
            //next line assumes that the basic type of force is **real**
            force[localIndex[k]] -= {
                    static_cast<real>(derivs(row, 0)),
                    static_cast<real>(derivs(row, 1)),
                    static_cast<real>(derivs(row, 2))
            };
        }
    } else {
        const auto &subGlobal0Ids = this->subGlobal0Ids_;
        for (int i = 0; i < derivs.rows(); i++) {
            const int *piLocal = kenrefDd(*currentInput_)->ga2la->findHome(subGlobal0Ids[i]);
            if (piLocal == nullptr)
                continue;         // not this rank's atom; its owner applies this row
            // TODO optimize this line/process
            force[*piLocal] -= {
                    static_cast<real>(derivs(i, 0)),
                    static_cast<real>(derivs(i, 1)),
                    static_cast<real>(derivs(i, 2))
            };
        }
    }
}

void KEnRefForceProvider::gatherFittedSubAtomsX(const std::vector<CoordsMatrixType<KEnRef_Real_t>> &localFitted,
                                                std::vector<CoordsMatrixType<KEnRef_Real_t>> &all) const {
    const auto &fitted = localFitted[0]; // exactly one model per replica
    if (!isMultiSimulation_) {
        all.assign(localFitted.begin(), localFitted.end()); // single replica: it is the whole ensemble
        return;
    }
    /* ONLY a simulation's main rank may touch mainRanksComm_: gmx_multisim_t documents it as "valid
     * only on main ranks" and it is MPI_COMM_NULL everywhere else, so calling a collective on it from
     * a non-main rank is undefined. With one rank per replica every rank is a main rank, which is why
     * this has always worked and why the guard below changes nothing there. */
    const int subAtomsXSize = static_cast<int>(fitted.size());
    if (!isReplicaMasterRank())
        return;                   // this rank's atoms already reached the main rank via the reduction
    MPI_Gather(const_cast<KEnRef_Real_t *>(fitted.data()), subAtomsXSize, KENREF_MPI_REAL,
               this->allSimulationsSubAtomsX_->data(), subAtomsXSize, KENREF_MPI_REAL, 0, mainRanksComm_);
    if (isEnsembleMasterRank()) {
        const Eigen::Index nSub = fitted.rows();
        all.resize(numSimulations_);
        for (int m = 0; m < numSimulations_; ++m) // zero-copy view of each model's slice, materialised by the copy
            all[m] = CoordsMapTypeConst<KEnRef_Real_t>(&allSimulationsSubAtomsX_->data()[m * nSub * 3], nSub, 3);
    }
}

void KEnRefForceProvider::scatterModelDerivatives(const std::vector<CoordsMatrixType<KEnRef_Real_t>> &allPerModel,
                                                  std::vector<CoordsMatrixType<KEnRef_Real_t>> &localPerModel) const {
    const int nSub = static_cast<int>(this->subAtomsX_->rows());
    const int subAtomsXSize = nSub * 3;
    if (!isMultiSimulation_) {
        /* Single replica. There is no cross-replica exchange, but under domain decomposition there is
         * still more than one rank here, and only the simulation's MAIN rank ran the model: the driver
         * gates that on simulationIndex()==0, which we deliberately report as -1 off the main rank so
         * the ensemble is assembled exactly once. So `allPerModel` is populated on the main rank and
         * EMPTY everywhere else, and simply copying it would leave the other ranks indexing an empty
         * vector. Size the destination and broadcast the main rank's result down the simulation.
         *
         * No-op at one rank per simulation: that rank is its own main rank, the copy below is the
         * whole operation and the broadcast short-circuits on a size-1 communicator. */
        /* The rank that fills the buffer MUST be the rank the broadcast is rooted at, which is the
         * replica master. With a single replica that is also the rank the driver let compute, so
         * allPerModel is populated exactly there -- assert it rather than rely on the coincidence. */
        GMX_ASSERT(isReplicaMasterRank() == isEnsembleMasterRank(),
                   "With one replica the replica master and the ensemble master must be the same rank");
        localPerModel.resize(1);
        if (isReplicaMasterRank()) {
            localPerModel[0] = allPerModel[0];
        } else {
            localPerModel[0].resize(nSub, 3);
            localPerModel[0].setZero();
        }
        kenrefBcastFromSimMain(*currentInput_, localPerModel[0].data(),
                               static_cast<std::size_t>(subAtomsXSize));
        return;
    }
    auto allDerivatives_buffer = this->allDerivatives_buffer_.get();
    auto derivatives_buffer = this->derivatives_buffer_.get();
    // Master packs the per-model derivatives contiguously, then scatters one block to each replica.
    // As in gatherFittedSubAtomsX, only main ranks may use mainRanksComm_.
    if (isReplicaMasterRank()) {
        if (isEnsembleMasterRank()) {
            for (int m = 0; m < numSimulations_; ++m)
                std::copy_n(allPerModel[m].data(), subAtomsXSize, &allDerivatives_buffer[m * subAtomsXSize]);
        }
        MPI_Scatter(allDerivatives_buffer, subAtomsXSize, KENREF_MPI_REAL,
                    derivatives_buffer, subAtomsXSize, KENREF_MPI_REAL, 0, mainRanksComm_);
    }
    /* Every rank needs the derivatives, because every rank applies forces to the atoms it owns -- but
     * only the main rank received them. Push them down the simulation. No-op at one rank per replica. */
    kenrefBcastFromSimMain(*currentInput_, derivatives_buffer,
                           static_cast<std::size_t>(subAtomsXSize));
    localPerModel.resize(1);
    localPerModel[0] = CoordsMapType<KEnRef_Real_t>(derivatives_buffer, nSub, 3); // copy out of the buffer
}

void KEnRefForceProvider::fillOwnedRows(CoordsMatrixType<KEnRef_Real_t> &dest,
                                        const std::optional<gmx::LocalAtomSet> &atomSet,
                                        const std::vector<int> &global0Ids,
                                        const gmx::ForceProviderInput &forceProviderInput,
                                        const char *what) const {
    const auto nRows = static_cast<std::size_t>(dest.rows());
    /* Under domain decomposition a rank owns only a subset of the atoms. Rows it does not own stay
     * ZERO here and are supplied by the reduction at the end -- a missing atom is another rank's
     * responsibility, not an error. */
    dest.setZero();

    const bool selfCheck = kenrefDdSelfCheckEnabled();
    std::vector<KEnRef_Real_t> ownerCount;
    if (selfCheck)
        ownerCount.assign(nRows, KEnRef_Real_t(0));

    if (atomSet.has_value()) {
        /* Preferred path. DD maintains these two arrays for us, so no atom is looked up at all:
         * localIndex()[k] is where the k-th atom this rank owns lives in the local coordinate array,
         * and collectiveIndex()[k] is that atom's slot in the set's own ordering -- i.e. exactly the
         * row it belongs in.
         *
         * TRAP: both underlying vectors are cleared and re-push_back()ed at every repartition, so the
         * ArrayRefs are re-fetched here on every step. Hoisting them out of the step would leave them
         * dangling the moment DD repartitions. */
        const gmx::ArrayRef<const int> localIndex = atomSet->localIndex();
        const gmx::ArrayRef<const int> collectiveIndex = atomSet->collectiveIndex();
        for (std::size_t k = 0; k < localIndex.size(); ++k) {
            const auto row = static_cast<Eigen::Index>(collectiveIndex[k]);
            const gmx::RVec atom_x = forceProviderInput.x_[localIndex[k]];
            for (int j = 0; j < 3; ++j)
                dest(row, j) = static_cast<KEnRef_Real_t>(atom_x[j]);
            if (selfCheck)
                ownerCount[static_cast<std::size_t>(row)] += KEnRef_Real_t(1);
        }
    } else {
        // Fallback when no LocalAtomSetManager was delivered: ask ga2la about every atom in the set
        // and skip the ones that are not home on this rank.
        for (std::size_t i = 0; i < nRows; ++i) {
            const int *piLocal = kenrefDd(forceProviderInput)->ga2la->findHome(global0Ids[i]);
            if (piLocal == nullptr)
                continue;         // not this rank's atom; another rank contributes this row
            const gmx::RVec atom_x = forceProviderInput.x_[*piLocal];
            for (int j = 0; j < 3; ++j)
                dest(static_cast<Eigen::Index>(i), j) = static_cast<KEnRef_Real_t>(atom_x[j]);
            if (selfCheck)
                ownerCount[i] += KEnRef_Real_t(1);
        }
    }

    if (selfCheck)
        kenrefCheckExactlyOneWriter(forceProviderInput, ownerCount, what, forceProviderInput.step_);

    // Combine the per-rank contributions BEFORE any scaling: every row has exactly one non-zero
    // contributor, so the sum is exact and rank-count independent. No-op at one rank per replica.
    kenrefSumOverSimRanks(forceProviderInput, dest.data(), nRows * 3u);
    /* Periodic re-imaging is NOT done here. Both sets must be made whole in one shared frame, so that
     * happens once in getLocalModelX() after both have been assembled -- see the comment there. */
}

void KEnRefForceProvider::fillSubAtomsX(CoordsMatrixType<KEnRef_Real_t> &subAtomsX,
                                        const gmx::ForceProviderInput &forceProviderInput,
                                        const bool toAngstrom) const {
    fillOwnedRows(subAtomsX, subAtomSet_, subGlobal0Ids_, forceProviderInput, "sub");
    if (toAngstrom)
        subAtomsX *= 10;
}

CoordsMatrixType<KEnRef_Real_t> KEnRefForceProvider::getGuideAtomsX(const gmx::ForceProviderInput &forceProviderInput,
                                                                    const bool toAngstrom) const {
    const auto guideAtom0IndicesSize = static_cast<Eigen::Index>(this->guideAtom0Indices_->size());
    auto guideAtomsX_ZEROIndexed = CoordsMatrixType<KEnRef_Real_t>(guideAtom0IndicesSize, 3);
    fillOwnedRows(guideAtomsX_ZEROIndexed, guideAtomSet_, *this->guideAtom0Indices_,
                  forceProviderInput, "guide");
#if VERBOSE
    std::cout << "guideAtomsX_ZEROIndexed shape is (" << guideAtomsX_ZEROIndexed.rows() << ", " << guideAtomsX_ZEROIndexed.cols() << ")" << std::endl;
    std::cout << "guideAtomsX_ZEROIndexed" << std::endl << guideAtomsX_ZEROIndexed << std::endl;
#endif
    if (toAngstrom)
        guideAtomsX_ZEROIndexed *= 10;
    return guideAtomsX_ZEROIndexed; //RETURN BY VALUE
}

void KEnRefForceProvider::setLocalAtomSetManager(gmx::LocalAtomSetManager *manager) {
    this->localAtomSetManager_ = manager;
}

void KEnRefForceProvider::setTopology(const gmx_mtop_t *mtop) {
    this->mtop_ = mtop;
}

void KEnRefForceProvider::initParamsAtSetup() {
    auto begin = std::chrono::high_resolution_clock::now();

    /* The multi-simulation facts are constant for the run, so read them once here. They used to be
     * read at step 0; calculateForces() still refreshes them every step, which is free and keeps a
     * single source of truth. numSimulations_ in particular is needed BELOW, by buildModelIndexing. */
    isMultiSimulation_ = this->simulationContext_->multiSimulation_ != nullptr;
    numSimulations_ = isMultiSimulation_ ? this->simulationContext_->multiSimulation_->numSimulations_ : 1;
    simulationIndex_ = isMultiSimulation_ ? this->simulationContext_->multiSimulation_->simulationIndex_ : 0;
    mainRanksComm_ = isMultiSimulation_ ? this->simulationContext_->multiSimulation_->mainRanksComm_ : MPI_COMM_NULL;
    const bool isMultiSimulation = isMultiSimulation_;
    const int numSimulations = numSimulations_;

    /* Redirect stdout before anything below prints, so the setup output still lands in the per-replica
     * file exactly as it did when this ran at step 0. */
    {
        std::string alt_out_path = Settings::alt_out_path;
        if (alt_out_path.empty()) {
            std::cout << "No alt_out_path defined. Will not redirect stdout stream." << std::endl;
        } else {
            if (isMultiSimulation_) {
                auto pos = alt_out_path.find(singleStr);
                std::string simIndexStr = std::to_string(simulationIndex_ + 1);
                if (pos != std::string::npos)
                    alt_out_path.replace(pos, strlen(singleStr), simIndexStr, 0, simIndexStr.length());
            }
            alt_out_path = std::filesystem::absolute(alt_out_path);
            std::cout << "Attempting to redirect output to " << alt_out_path << std::endl;
            if (std::freopen(alt_out_path.c_str(), "a", stdout)) {
                std::cout << "Redirected stdout stream" << std::endl;
            } else {
                std::cout << "FAILED to redirect output to " << alt_out_path << std::endl;
            }
        }
    }

    std::cout << "Energy model: " << selectedEnergyModel << std::endl;

    this->atomName_to_atomGlobalId_map_ = std::make_shared<std::map<std::string, int> >(
        IoUtils::getAtomMappingFromPdb<std::string, int>(Settings::atomName_mapping_fileName, IoUtils::fill_atomId_to_index_Map));
    GMX_ASSERT(!atomName_to_atomGlobalId_map_->empty(), "No atom mapping found");
    auto &atomName_to_atomGlobalId_map = *this->atomName_to_atomGlobalId_map_;

    const KEnRef_Real_t maxForceSquared = std::pow(Settings::max_force, 2.f);
    const KEnRef_Real_t k = Settings::k;
    const KEnRef_Real_t n = Settings::n;
    const bool handleNames = IoUtils::should_handleNames(atomName_to_atomGlobalId_map);

    std::cout << "KEnRef_Real_t type is: " << typeid(KEnRef_Real_t).name() << '\n';

    // ---- shared model + sub-atom indexing (registry create + buildCache + sub0Id indexing +
    //      finalizeIndexing), replacing the per-model SIGMA/PLATEAUS switch AND the hand-rolled
    //      sub-indexing. The 1-based global serials we pass in come back as 1-based subAtomGlobalIds. ----
    const std::string& modelName = selectedEnergyModel;
    GMX_ASSERT(kenref::ModelRegistry<KEnRef_Real_t>::has(modelName),
               "Unknown energy model. please set \"--model\" to a registered model (e.g. SIGMA, PLATEAUS, RELAX).");
    auto mi = kenref::buildModelIndexing<KEnRef_Real_t>(
        modelName, *this, atomName_to_atomGlobalId_map, handleNames, numSimulations, numOmpThreads());
    this->sub0Id_to_global1Id_ = std::make_shared<std::vector<int> >(std::move(mi.subAtomGlobalIds));

    /* ---- register the two atom sets with domain decomposition ----
     * This is the ONLY point at which registration is legal: we are inside initForceProviders(), which
     * runner.cpp calls after setupNotifier.notify(&atomSets) (so the manager exists) but before the
     * first dd_partition_system() in the MD loop (so DD will compute this set's indices from the very
     * first partition). LocalAtomSetManager::add() copies the indices, so the temporaries below may die.
     *
     * Both sets are built from ZERO-based global indices: LocalAtomSet speaks GROMACS's global atom
     * numbering, whereas sub0Id_to_global1Id_ holds ONE-based serials from the mapping PDB. */
    this->subGlobal0Ids_.clear();
    this->subGlobal0Ids_.reserve(this->sub0Id_to_global1Id_->size());
    for (const int global1Id : *this->sub0Id_to_global1Id_)
        this->subGlobal0Ids_.push_back(global1Id - 1);

    /* ---- bond connectivity for periodic make-whole ----
     * Built over every molecule containing a restrained atom, from the global topology. This is what
     * lets a molecule broken across a periodic boundary be repaired: walking bonds is unambiguous at
     * any molecule size, whereas every distance-based rule fails once the molecule approaches half the
     * box (measured here: 2.171 nm reach versus a 2.165 nm half-box). */
    if (mtop_ != nullptr) {
        std::vector<int> restrained = this->subGlobal0Ids_;
        restrained.insert(restrained.end(), this->guideAtom0Indices_->begin(),
                          this->guideAtom0Indices_->end());
        moleculeGraph_.build(*mtop_, restrained);

        if (moleculeGraph_.isBuilt()) {
            moleculeX_.resize(static_cast<Eigen::Index>(moleculeGraph_.atoms().size()), 3);
            guideRowInMolecule_.clear();
            guideRowInMolecule_.reserve(this->guideAtom0Indices_->size());
            for (const int g : *this->guideAtom0Indices_)
                guideRowInMolecule_.push_back(moleculeGraph_.rowOf(g));
            subRowInMolecule_.clear();
            subRowInMolecule_.reserve(this->subGlobal0Ids_.size());
            for (const int g : this->subGlobal0Ids_)
                subRowInMolecule_.push_back(moleculeGraph_.rowOf(g));
            std::cout << "KEnRef: molecule connectivity built -- " << moleculeGraph_.atoms().size()
                      << " atoms in " << moleculeGraph_.numFragments() << " fragment(s)" << std::endl;
        }
    } else {
        std::cout << "KEnRef: no topology available; periodic make-whole is disabled" << std::endl;
    }

    if (localAtomSetManager_ != nullptr) {
        subAtomSet_ = localAtomSetManager_->add(this->subGlobal0Ids_);
        guideAtomSet_ = localAtomSetManager_->add(*this->guideAtom0Indices_); // already zero-based
        if (moleculeGraph_.isBuilt())
            moleculeAtomSet_ = localAtomSetManager_->add(moleculeGraph_.atoms());
        std::cout << "KEnRef: registered LocalAtomSets -- " << this->subGlobal0Ids_.size()
                  << " sub atoms, " << this->guideAtom0Indices_->size() << " guide atoms" << std::endl;
    } else {
        std::cout << "KEnRef: no LocalAtomSetManager available; falling back to ga2la lookups"
                  << std::endl;
    }

    // ---- per-step buffers (allocated once, reused every step) ----
    this->subAtomsX_ = std::make_shared<CoordsMatrixType<KEnRef_Real_t> >(this->sub0Id_to_global1Id_->size(), 3); //contains needed atoms only
    this->allSimulationsSubAtomsX_ =
            isMultiSimulation
                ? std::make_shared<CoordsMatrixType<KEnRef_Real_t> >(numSimulations * this->subAtomsX_->rows(), 3)
                : this->subAtomsX_;
    this->allDerivatives_buffer_ = std::shared_ptr<KEnRef_Real_t[]>(
        new KEnRef_Real_t[this->subAtomsX_->size() * numSimulations]);
    this->derivatives_buffer_ =
        isMultiSimulation
            ? std::shared_ptr<KEnRef_Real_t[]>(new KEnRef_Real_t[this->subAtomsX_->size()])
            : this->allDerivatives_buffer_;

    // ---- construct the driver that owns the model + the shared per-step pipeline ----
    this->driver_ = std::make_unique<kenref::KEnRefDriver<KEnRef_Real_t> >(
        std::move(mi.model), k, n, maxForceSquared, *this->guideAtomsReferenceCoordsCentered_);
    // Live MD: the box really is the periodic box, so the split check is meaningful. It should never
    // fire here, because make-whole runs first -- it is the backstop for when that fails or is off.
    this->driver_->enablePeriodicSplitCheck();

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    this->calculateForces_time -= elapsed.count(); //Exclude THIS method from wall time calculation
}

#undef VERBOSE
