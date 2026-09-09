/*
 * KEnRefOwnership.h
 *
 * The pure decision behind KENREF_DD_SELFCHECK, split out so it can be tested without MPI.
 *
 * WHY THIS IS A SEPARATE HEADER. Under domain decomposition each rank fills only the rows it owns and
 * leaves the rest at zero; summing the per-rank buffers then reconstructs the set exactly, because
 * `0 + x == x` in IEEE-754 and addition is commutative. That is exact ONLY if every row has exactly
 * one non-zero contributor. Two ranks claiming a row would double it; none claiming it would leave it
 * at the origin. Neither crashes, and a wrong coordinate merely biases the restraint -- so the
 * property has to be tested for rather than waited for.
 *
 * The check that enforces it used to live entirely inside an MPI collective in
 * KEnRefForceProvider.cpp, which made it unreachable from a unit test. Everything below the reduction
 * is pure arithmetic on a vector of counts, so it lives here and `google_tests/testOwnershipCheck.cpp`
 * drives it with hand-built vectors. The reduction itself stays in the .cpp and is exercised by
 * `dd_validate.sh`'s real multi-rank runs.
 *
 * No GROMACS and no MPI, so the decision is testable without a harness. `core/KEnRef.h` comes along
 * for `KEnRef_Real_t`, which keeps the scalar type in lockstep with `KENREF_DOUBLE` rather than
 * hardcoding double; that transitively pulls Eigen, which costs nothing because every consumer of
 * this header already has it.
 *
 * HISTORY, because it is the reason this file exists at all. The self-check was originally proven to
 * FIRE -- not merely to pass -- by a temporary `KENREF_DD_FAULT` environment hook that corrupted
 * `ownerCount` just before the check. That hook was NEVER COMMITTED: it carried the comment "Remove
 * before committing", and `git log -S KENREF_DD_FAULT --all -- '*.cpp' '*.h'` is empty across the
 * whole history. It survived only inside one session transcript under `~/.claude/projects/...`, from
 * which it was recovered on 2026-09-09 -- a single JSONL file, in a directory nobody backs up, as the
 * sole evidence that the project's parallel-correctness check works at all.
 *
 * The two fault modes that hook injected are reproduced verbatim as test cases, so this file
 * demonstrably covers what that demonstration covered rather than merely resembling it:
 *     mode '2'  ownerCount[0] += 1                 -- a row claimed twice
 *     mode '0'  ownerCount[size() - 1] = 0         -- a row claimed by nobody
 */

#ifndef KENREFOWNERSHIP_H_
#define KENREFOWNERSHIP_H_

#include <cstddef>
#include <vector>

#include "core/KEnRef.h"

namespace kenref {

/*! \brief What a global owner-count vector says about who wrote each row.
 *
 * \p firstBad is the index of the first offending row and is meaningful only when !ok(); it is left
 * at 0 otherwise, which is why ok() must be consulted first rather than firstBad tested directly. */
struct OwnershipVerdict {
	std::size_t unowned   = 0;  //!< rows no rank claimed (count 0)
	std::size_t contested = 0;  //!< rows more than one rank claimed (count > 1)
	std::size_t firstBad  = 0;  //!< index of the first offending row; only meaningful when !ok()

	//! Every row was written exactly once, which is what makes the zero-fill sum exact.
	[[nodiscard]] bool ok() const { return unowned == 0 && contested == 0; }
};

/*! \brief Classify a GLOBAL owner-count vector -- one entry per row, already reduced across ranks.
 *
 * Exactly 1 is the only good value. 0 means nobody wrote the row; anything else means more than one
 * rank did. Counts arrive as reals because they ride the same floating-point reduction as the
 * coordinates, so they are compared against exact integral values -- which is safe here: small
 * integers are exact in both float and double, and the sum of a handful of exact 1.0s is exact.
 *
 * An empty vector is vacuously fine: there are no rows to own. */
[[nodiscard]] inline OwnershipVerdict verifyExactlyOneWriter(const KEnRef_Real_t *counts, std::size_t n) {
	OwnershipVerdict v;
	bool haveBad = false;
	for (std::size_t i = 0; i < n; ++i) {
		const KEnRef_Real_t c = counts[i];
		if (c == KEnRef_Real_t(0))
			++v.unowned;
		else if (c != KEnRef_Real_t(1))
			++v.contested;
		else
			continue;
		if (!haveBad) {
			v.firstBad = i;
			haveBad    = true;
		}
	}
	return v;
}

//! Convenience overload for the vector the force provider actually holds.
[[nodiscard]] inline OwnershipVerdict verifyExactlyOneWriter(const std::vector<KEnRef_Real_t> &counts) {
	return verifyExactlyOneWriter(counts.data(), counts.size());
}

} // namespace kenref

#endif /* KENREFOWNERSHIP_H_ */
