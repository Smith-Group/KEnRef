/*
 * KEnRefForceProvider.cpp
 *
 *      Author: amr
 */

#include <iostream>
//#include <typeinfo>
#include <cmath>
#include <filesystem>   // std::filesystem::absolute — pulled in transitively by libc++, NOT by libstdc++
#include <memory>
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

    if (!paramsInitialized) {
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

    if (step % 10 == 0)
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

        /* KEnRef supports exactly ONE RANK PER SIMULATION, and refuses anything else here rather than
         * producing wrong numbers.
         *
         * Under domain decomposition each rank owns only a subset of the atoms, and the local ordering is
         * re-sorted at every repartition. This code assumes the opposite: it resolves every sub-atom and
         * guide atom through ga2la->findHome() and dereferences the result without checking, and the
         * per-replica MPI_Gather/MPI_Scatter below run on mainRanksComm_, which is only a valid
         * communicator on a simulation's main rank. With one rank per simulation both assumptions hold --
         * every atom is home, and every rank is a main rank -- which is why this has never been hit.
         *
         * Left unguarded the failure is silent then fatal: the energy comes out `nan` at step 0, the box
         * follows, and the run dies later inside MPI or on an illegal instruction, with nothing in the
         * output naming KEnRef. Refuse up front instead. Removing this guard is the acceptance test for
         * real DD support, not an optimisation. */
        if (kenrefHavePPDomainDecomposition(forceProviderInput))
        {
            gmx_fatal(FARGS,
                      "KEnRef does not support domain decomposition: this simulation has %d PP ranks. "
                      "Run KEnRef with ONE RANK PER SIMULATION -- e.g. `mpirun -np <N> ... -multidir "
                      "<N directories>`, which gives each replica a single rank -- or add `-ntmpi 1` / "
                      "reduce the rank count so no simulation is decomposed. Threads are unaffected: "
                      "-ntomp is free to use.",
                      kenrefDd(forceProviderInput)->nnodes);
        }

        std::cout << "Number of atoms = " << homenr << std::endl;
        std::cout << "havePPDomainDecomposition: " << kenrefHavePPDomainDecomposition(forceProviderInput) << std::endl;
        std::cout << "haveDDAtomOrdering: " << kenrefHaveDDAtomOrdering(forceProviderInput) << std::endl;
        std::cout << "dd->nnodes: " << kenrefDd(forceProviderInput)->nnodes << std::endl;
        fillParamsStep0(homenr, numSimulations_, forceProviderInput); //TODO Optimize this function via OMP
        paramsInitialized = true;
    }

    // Run the shared per-step refinement pipeline: getLocalModelX -> no-jump -> Kabsch fit -> gather
    // -> model.compute -> scatter -> inverse-fit -> unit-scale -> saturate -> addLocalModelDerivatives.
    const KEnRef_Real_t energy = driver_->step(*this, /*printStatistics*/ (step % 10 == 0));

    if ((step % 10 == 0) && simulationIndex_ == 0)
        std::cout << "Step: " << step << " Energy: " << energy << std::endl;

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    this->calculateForces_time += elapsed.count();
    if (!(step % 10) && simulationIndex_ == 0) {
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
    guideX = getGuideAtomsX(*this->guideAtom0Indices_, *currentInput_, true);
    subX.resize(static_cast<Eigen::Index>(this->sub0Id_to_global1Id_->size()), 3);
    fillSubAtomsX(subX, *this->sub0Id_to_global1Id_, *currentInput_, true);
    // GROMACS box (nm, real[3][3]) -> Eigen 3x3 (raw; the driver scales it to Angstrom via toAngstrom).
    const matrix &b = currentInput_->box_;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            box(i, j) = static_cast<KEnRef_Real_t>(b[i][j]);
}

void KEnRefForceProvider::addLocalModelDerivatives(int /*localModel*/, const CoordsMatrixType<KEnRef_Real_t> &derivs) {
    // N.B. this assumes all ranks hit this line; mirrors the original pre-apply barrier.
    if (isMultiSimulation_ && kenrefHaveDDAtomOrdering(*currentInput_)) {
        gmx_barrier(kenrefGroupComm(*currentInput_));
    }
    const auto &force = currentOutput_->forceWithVirial_.force_;
    const auto &sub0Id_to_global1Id = *this->sub0Id_to_global1Id_;
    //Finally, add them to corresponding atoms
    /* Each rank applies forces ONLY to the atoms it owns. There is no reduction here and there must
     * not be one: the force array is indexed by local atom, the rank sets are disjoint, and a rank
     * that does not own an atom has nowhere to put its force. Skipping is therefore the whole fix --
     * a non-owned atom is another rank's responsibility, not a lost contribution. */
    for (int i = 0; i < derivs.rows(); i++) {
        const int global0Id = sub0Id_to_global1Id[i] - 1;
        const int *piLocal = kenrefDd(*currentInput_)->ga2la->findHome(global0Id);
        if (piLocal == nullptr)
            continue;             // not this rank's atom; its owner applies this row
        //next line assumes that the basic type of force is **real**
        // TODO optimize this line/process
        force[*piLocal] -= {
                static_cast<real>(derivs(i, 0)),
                static_cast<real>(derivs(i, 1)),
                static_cast<real>(derivs(i, 2))
        };
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
    if (!isSimMainRank_)
        return;                   // this rank's atoms already reached the main rank via the reduction
    MPI_Gather(const_cast<KEnRef_Real_t *>(fitted.data()), subAtomsXSize, KENREF_MPI_REAL,
               this->allSimulationsSubAtomsX_->data(), subAtomsXSize, KENREF_MPI_REAL, 0, mainRanksComm_);
    if (simulationIndex_ == 0) {
        const Eigen::Index nSub = fitted.rows();
        all.resize(numSimulations_);
        for (int m = 0; m < numSimulations_; ++m) // zero-copy view of each model's slice, materialised by the copy
            all[m] = CoordsMapTypeConst<KEnRef_Real_t>(&allSimulationsSubAtomsX_->data()[m * nSub * 3], nSub, 3);
    }
}

void KEnRefForceProvider::scatterModelDerivatives(const std::vector<CoordsMatrixType<KEnRef_Real_t>> &allPerModel,
                                                  std::vector<CoordsMatrixType<KEnRef_Real_t>> &localPerModel) const {
    if (!isMultiSimulation_) {
        localPerModel = allPerModel; // single replica
        return;
    }
    const int nSub = static_cast<int>(this->subAtomsX_->rows());
    const int subAtomsXSize = nSub * 3;
    auto allDerivatives_buffer = this->allDerivatives_buffer_.get();
    auto derivatives_buffer = this->derivatives_buffer_.get();
    // Master packs the per-model derivatives contiguously, then scatters one block to each replica.
    // As in gatherFittedSubAtomsX, only main ranks may use mainRanksComm_.
    if (isSimMainRank_) {
        if (simulationIndex_ == 0) {
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

void KEnRefForceProvider::fillSubAtomsX(CoordsMatrixType<KEnRef_Real_t> &subAtomsX,
                                        const std::vector<int> &sub0Id_to_global1Id,
                                        const gmx::ForceProviderInput &forceProviderInput, const bool toAngstrom) {
    /* Under domain decomposition a rank owns only a subset of the atoms, so findHome() legitimately
     * returns null for the rest -- it is NOT an error. Leave those rows at zero and let
     * kenrefSumOverSimRanks() below fill them in from the rank that does own them. With one rank per
     * replica every atom is home, the reduction short-circuits, and this is bit-identical to the
     * previous unguarded code. */
    subAtomsX.setZero();
    for (int i = 0; i < subAtomsX.rows(); i++) {
        const int global0Id = sub0Id_to_global1Id[i] - 1;
        const int *piLocal = kenrefDd(forceProviderInput)->ga2la->findHome(global0Id);
        if (piLocal == nullptr)
            continue;             // not this rank's atom; another rank contributes this row
        const gmx::RVec atom_x = forceProviderInput.x_[*piLocal];
#if VERBOSE
        std::cout << sub0Id_to_global1Id[i] << "\t" << global0Id << "\t" << *piLocal << "\t x: " << atom_x[0] << ", " << atom_x[1] << ", " << atom_x[2] << std::endl;
#endif
        // (rows this rank does not own were skipped above and stay zero until the reduction)
        if constexpr (std::is_same_v<KEnRef_Real_t, real>) {
            auto subAtomsX_buffer = subAtomsX.data();
            const auto rvec = atom_x.as_vec();
            std::copy_n(rvec, 3, &subAtomsX_buffer[i * 3]);
        } else {
            for (int j = 0; j < 3; ++j) {
                subAtomsX(i, j) = static_cast<KEnRef_Real_t>(atom_x[j]);
            }
        }
    }
    // Combine the per-rank contributions BEFORE scaling: every row now has exactly one non-zero
    // contributor, so the sum is exact and rank-count independent. No-op at one rank per replica.
    kenrefSumOverSimRanks(forceProviderInput, subAtomsX.data(),
                          static_cast<std::size_t>(subAtomsX.rows()) * 3u);
    if (toAngstrom)
        subAtomsX *= 10;
}

CoordsMatrixType<KEnRef_Real_t> KEnRefForceProvider::getGuideAtomsX(const std::vector<int> &guideAtom0Indices,
                                                                    const gmx::ForceProviderInput &forceProviderInput,
                                                                    const bool toAngstrom) {
    long guideAtom0IndicesSize = static_cast<long>(guideAtom0Indices.size());
    auto guideAtomsX_ZEROIndexed = CoordsMatrixType<KEnRef_Real_t>(guideAtom0IndicesSize, 3);
    KEnRef_Real_t *guideAtomsX_ZEROIndexed_buffer = guideAtomsX_ZEROIndexed.data();
    // Same zero-fill convention as fillSubAtomsX: rows this rank does not own stay zero and are
    // supplied by the reduction below. findHome() returning null is expected under DD, not an error.
    guideAtomsX_ZEROIndexed.setZero();
    for (auto i = 0; i < guideAtom0IndicesSize; i++) {
        const int *pi = &guideAtom0Indices[i];
        const int *piLocal = kenrefDd(forceProviderInput)->ga2la->findHome(*pi);
        if (piLocal == nullptr)
            continue;             // not this rank's atom
        const gmx::RVec atom_x = forceProviderInput.x_[*piLocal];

        auto rvec = atom_x.as_vec();
        std::copy_n(rvec, 3, &guideAtomsX_ZEROIndexed_buffer[i * 3]);
    }
    kenrefSumOverSimRanks(forceProviderInput, guideAtomsX_ZEROIndexed.data(),
                          static_cast<std::size_t>(guideAtom0IndicesSize) * 3u);
    //TODO make a unit test to validate that the value coming in rvec is equal to the value in guideAtomsX_ZEROIndexed
#if VERBOSE
    std::cout << "guideAtomsX_ZEROIndexed shape is (" << guideAtomsX_ZEROIndexed.rows() << ", " << guideAtomsX_ZEROIndexed.cols() << ")" << std::endl;
    std::cout << "guideAtomsX_ZEROIndexed" << std::endl << guideAtomsX_ZEROIndexed << std::endl;
#endif
    if (toAngstrom)
        guideAtomsX_ZEROIndexed *= 10;
    return guideAtomsX_ZEROIndexed; //RETURN BY VALUE
}

void KEnRefForceProvider::fillParamsStep0(const size_t homenr, int numSimulations, const gmx::ForceProviderInput &forceProviderInput) {
    auto begin = std::chrono::high_resolution_clock::now();
    bool isMultiSimulation = this->simulationContext_->multiSimulation_ != nullptr;
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

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    this->calculateForces_time -= elapsed.count(); //Exclude THIS method from wall time calculation
}

#undef VERBOSE
