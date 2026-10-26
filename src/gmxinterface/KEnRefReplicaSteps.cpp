/*
 * KEnRefReplicaSteps.cpp
 *
 * The collective half of the mismatched-steps refusal. See the header for why it exists and what is
 * compared.
 */

#include "gmxinterface/KEnRefReplicaSteps.h"

#include <filesystem>

#include "mpi.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/mdrunutility/multisim.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/fatalerror.h"

void kenref::refuseMismatchedReplicaSteps(const gmx_multisim_t *ms, const std::string &tprFile,
                                          std::int64_t nstepsCommandline,
                                          const std::vector<std::string> &replicaNames) {
	// Only a replica's main rank is in mainRanksComm_; on every other rank it is MPI_COMM_NULL.
	if (ms == nullptr || ms->mainRanksComm_ == MPI_COMM_NULL)
		return;

	/* Read this replica's tpr. GROMACS reads it again later, on the same rank; doing it here, before
	 * the runner exists, is what lets the refusal come before any replica takes a step. */
	t_inputrec ir;
	t_state    state;
	gmx_mtop_t mtop;
	read_tpx_state(std::filesystem::path(tprFile), &ir, &state, &mtop);
	const std::int64_t key = replicaStepKey(ir.init_step, ir.nsteps, nstepsCommandline);

	// mainRanksComm_ holds one rank per replica, in simulation order, so rank i is replica i.
	std::vector<std::int64_t> keys(ms->numSimulations_);
	MPI_Allgather(&key, 1, MPI_INT64_T, keys.data(), 1, MPI_INT64_T, ms->mainRanksComm_);

	const ReplicaStepsVerdict verdict = checkReplicaSteps(keys, nstepsCommandline > kNstepsFromTpr, replicaNames);
	if (verdict.consistent)
		return;
	// The same pattern as GROMACS's own initial-step refusal: one rank reports, the rest wait to be aborted.
	if (ms->simulationIndex_ == 0)
		gmx_fatal(FARGS, "%s", verdict.report.c_str());
	MPI_Barrier(ms->mainRanksComm_);
}
