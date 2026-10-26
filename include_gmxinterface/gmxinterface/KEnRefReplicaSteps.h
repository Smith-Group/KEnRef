/*
 * KEnRefReplicaSteps.h
 *
 * Refuse a multi-simulation whose replicas would not run the same number of steps.
 *
 * WHY. In a -multidir run KEnRefForceProvider does an MPI_Gather and an MPI_Scatter on the
 * replicas' main-rank communicator EVERY STEP. If one replica has fewer steps to run, it finishes
 * and exits, and the others block forever in the next collective. The run HANGS instead of failing.
 *
 * GROMACS does not prevent it. logInitialMultisimStatus() (mdrunutility/multisim.cpp) compares the
 * replicas' nsteps and only warns -- "The number of steps is not consistent across multi
 * simulations, but we are proceeding anyway!" -- whatever simulationsShareState says. Only a
 * mismatched INITIAL step is fatal, and only when the simulations share state, which
 * KEnRefMDModule's multisim subscription makes true. The only GROMACS check on the step count
 * itself, check_multi_int64(..., init_step + nsteps, ...), is reached through -replex alone.
 *
 * WHAT IS COMPARED. GROMACS ends a run at step_rel == nsteps, after these adjustments to the tpr's
 * values (runner.cpp, checkpoint.cpp; identical in 2025.4 and 2026.3):
 *   - restarting from a checkpoint sets  nsteps += init_step - checkpoint_step,  init_step =
 *     checkpoint_step. The LAST step, init_step + nsteps, is unchanged by a restart.
 *   - `mdrun -nsteps N` (anything but the default -2) then sets nsteps = N, AFTER the checkpoint.
 * So the last step follows from the tpr alone unless -nsteps was given, and then every replica runs N
 * steps from wherever it starts. Either way, together with the equal initial step GROMACS already
 * enforces, equal keys mean equal numbers of steps. replicaStepKey() is that value; -1 is
 * "unlimited" (nsteps = -1), which ends only on -maxh or a signal, and the multisim subscription
 * propagates those across replicas.
 *
 * The decision is pure and lives in this header so google_tests/testReplicaSteps.cpp can drive it
 * without MPI. refuseMismatchedReplicaSteps() (KEnRefReplicaSteps.cpp) reads the tpr and does the
 * collective.
 */

#ifndef KENREFREPLICASTEPS_H_
#define KENREFREPLICASTEPS_H_

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

struct gmx_multisim_t;

namespace kenref {

//! mdrun's default for -nsteps: "take nsteps from the tpr".
constexpr std::int64_t kNstepsFromTpr = -2;

/*! \brief The value every replica must agree on, or -1 for an unlimited run.
 *
 * \p nstepsCommandline is mdrun's -nsteps (kNstepsFromTpr when not given). Given, it is the number of
 * steps each replica runs from its starting step; otherwise the key is the tpr's last step. The
 * command line is shared by all replicas, so all of them take the same branch. */
[[nodiscard]] inline std::int64_t replicaStepKey(std::int64_t tprInitStep, std::int64_t tprNsteps,
                                                 std::int64_t nstepsCommandline) {
	if (nstepsCommandline > kNstepsFromTpr)
		return nstepsCommandline < 0 ? -1 : nstepsCommandline;
	if (tprNsteps < 0)
		return -1;
	return tprInitStep + tprNsteps;
}

//! What the replicas' keys say. \p report is empty when they agree.
struct ReplicaStepsVerdict {
	bool        consistent = true;
	std::string report;
};

/*! \brief Compare the keys of all replicas, one per replica in simulation order.
 *
 * \p replicaNames are the -multidir directories, used only to label the report; replicas without a
 * name are labelled by index alone. Fewer than two replicas are vacuously consistent. */
[[nodiscard]] inline ReplicaStepsVerdict checkReplicaSteps(const std::vector<std::int64_t> &keys,
                                                           bool fromCommandLine,
                                                           const std::vector<std::string> &replicaNames) {
	ReplicaStepsVerdict v;
	for (std::size_t i = 1; i < keys.size(); ++i)
		v.consistent = v.consistent && keys[i] == keys[0];
	if (v.consistent)
		return v;

	v.report = "KEnRef: the replicas of this multi-simulation would not run the same number of steps.\n"
	           "KEnRef exchanges data between the replicas at every step, so the run would hang as soon as "
	           "the first replica finished, instead of failing. GROMACS only warns about this.\n";
	v.report += fromCommandLine ? "Steps given by mdrun -nsteps, per replica:\n"
	                            : "Last step of each replica (init-step + nsteps from its tpr):\n";
	for (std::size_t i = 0; i < keys.size(); ++i) {
		v.report += "  replica " + std::to_string(i);
		if (i < replicaNames.size())
			v.report += " (" + replicaNames[i] + ")";
		v.report += ": " + (keys[i] < 0 ? std::string("unlimited (nsteps = -1)") : std::to_string(keys[i])) + "\n";
	}
	v.report += "Give every replica's tpr the same nsteps and init-step (gmx convert-tpr -nsteps), or pass "
	            "mdrun -nsteps.";
	return v;
}

/*! \brief Refuse, at setup, a multi-simulation whose replicas would run different numbers of steps.
 *
 * Collective over the replicas' main ranks: each reads its own tpr (\p tprFile, relative to the
 * replica directory mdrun has already changed into) and the keys are all-gathered. On a mismatch the
 * first replica's main rank calls gmx_fatal, which aborts every replica, and the others wait for it.
 * A no-op for a single simulation and on ranks that are not a replica's main rank. */
void refuseMismatchedReplicaSteps(const gmx_multisim_t *ms, const std::string &tprFile,
                                  std::int64_t nstepsCommandline, const std::vector<std::string> &replicaNames);

} // namespace kenref

#endif /* KENREFREPLICASTEPS_H_ */
