/*
 * testOwnershipCheck.cpp
 *
 * Tests for kenref::verifyExactlyOneWriter -- the pure decision behind KENREF_DD_SELFCHECK.
 *
 * WHAT THIS REPLACES, AND WHY IT EXISTS AT ALL. The domain-decomposition self-check was once proven to
 * actually FIRE -- not merely to pass, which is a much weaker statement -- by a temporary
 * `KENREF_DD_FAULT` environment hook that corrupted the owner-count vector immediately before the
 * check ran. That hook was NEVER COMMITTED. It carried the comment "Remove before committing", and
 *
 *     git log -S KENREF_DD_FAULT --all -- '*.cpp' '*.h'
 *
 * is empty across the entire history (with `-S KENREF_DD_SELFCHECK` as the positive control). It
 * survived in exactly one place: a session transcript under ~/.claude/projects/, from which it was
 * recovered on 2026-09-09. For the intervening period the evidence that this project's central
 * parallel-correctness check worked at all was a single JSONL file in a directory nobody backs up,
 * and re-running the demonstration meant re-applying a snippet by hand.
 *
 * These tests retire that dependency. THE TWO FAULT MODES BELOW ARE THE RECOVERED HOOK'S, VERBATIM:
 *
 *     if (fault[0] == '2') ownerCount[0] += KEnRef_Real_t(1);              // row claimed twice
 *     if (fault[0] == '0') ownerCount[ownerCount.size() - 1] = KEnRef_Real_t(0); // row claimed by nobody
 *
 * reproduced as RecoveredFaultMode2_RowClaimedTwice and RecoveredFaultMode0_RowClaimedByNobody. They
 * are written that way on purpose, so this file demonstrably covers what the original demonstration
 * covered rather than merely resembling it. Do not "tidy" them into a single parameterised case: the
 * point is the correspondence, not the coverage.
 *
 * TWO KINDS OF TEST LIVE HERE AND BOTH ARE LOAD-BEARING. The named cases pin the shapes someone
 * thought of, so a failure says which property broke. MatchesThePreExtractionLoopOnRandomVectors
 * compares against the pre-extraction loop over shapes nobody thought of, so it catches an
 * extraction that is wrong in a way no hand-written case anticipates -- but it can only say "you
 * changed behaviour". Neither subsumes the other; do not delete either as redundant.
 *
 * Scope: the pure decision only. The MPI reduction that produces the global counts stays in
 * KEnRefForceProvider.cpp and is exercised by dd_validate.sh's real multi-rank runs -- 1, 2 and 4
 * ranks and -npme 1, with KENREF_DD_SELFCHECK=1 on throughout.
 */

#include "gtest/gtest.h"

#include <random>
#include <vector>

#include "gmxinterface/KEnRefOwnership.h"

using kenref::OwnershipVerdict;
using kenref::verifyExactlyOneWriter;

namespace {

//! The good case: one writer per row, which is what the zero-fill reduction requires.
std::vector<KEnRef_Real_t> wellOwned(std::size_t n) {
	return std::vector<KEnRef_Real_t>(n, KEnRef_Real_t(1));
}

/*! \brief The ORIGINAL loop, lifted verbatim from KEnRefForceProvider.cpp before the extraction.
 *
 * Kept as a reference implementation so the differential test below can prove the extraction is
 * behaviour-preserving rather than merely assert it. Do not "simplify" this to call the real
 * function -- that would make the comparison vacuous, which is the one way this test can fail
 * silently. If the real decision ever legitimately changes, this copy must be updated deliberately
 * and the commit should say why. */
struct OriginalLoopResult {
	std::size_t unowned = 0, contested = 0, firstBad = 0;
};
OriginalLoopResult originalLoop(const std::vector<KEnRef_Real_t> &ownerCount) {
	OriginalLoopResult o;
	bool haveBad = false;
	for (std::size_t i = 0; i < ownerCount.size(); ++i) {
		const KEnRef_Real_t c = ownerCount[i];
		if (c == KEnRef_Real_t(0))
			++o.unowned;
		else if (c != KEnRef_Real_t(1))
			++o.contested;
		else
			continue;
		if (!haveBad) {
			o.firstBad = i;
			haveBad    = true;
		}
	}
	return o;
}

} // namespace

// ---------------------------------------------------------------------------------------------
// The property itself.
// ---------------------------------------------------------------------------------------------

TEST(OwnershipCheckTest, AcceptsExactlyOneWriterPerRow) {
	const OwnershipVerdict v = verifyExactlyOneWriter(wellOwned(64));
	EXPECT_TRUE(v.ok());
	EXPECT_EQ(v.unowned, 0u);
	EXPECT_EQ(v.contested, 0u);
}

/* An empty set is vacuously fine -- there are no rows to own. Worth pinning: the force provider calls
 * this for every registered atom set, and a set can legitimately be empty on a rank that owns none. */
TEST(OwnershipCheckTest, AcceptsAnEmptySet) {
	EXPECT_TRUE(verifyExactlyOneWriter(std::vector<KEnRef_Real_t>{}).ok());
}

// ---------------------------------------------------------------------------------------------
// The two modes the recovered KENREF_DD_FAULT hook injected, reproduced verbatim.
// ---------------------------------------------------------------------------------------------

//! Recovered hook, mode '2': ownerCount[0] += 1. Two ranks claim row 0, so its coordinate doubles.
TEST(OwnershipCheckTest, RecoveredFaultMode2_RowClaimedTwice) {
	std::vector<KEnRef_Real_t> counts = wellOwned(32);
	counts[0] += KEnRef_Real_t(1);

	const OwnershipVerdict v = verifyExactlyOneWriter(counts);
	EXPECT_FALSE(v.ok()) << "a doubly-claimed row must be caught; this is the mode the never-committed "
	                        "KENREF_DD_FAULT=2 hook injected";
	EXPECT_EQ(v.contested, 1u);
	EXPECT_EQ(v.unowned, 0u);
	EXPECT_EQ(v.firstBad, 0u);
}

//! Recovered hook, mode '0': ownerCount[size()-1] = 0. Nobody claims the last row; it stays at the origin.
TEST(OwnershipCheckTest, RecoveredFaultMode0_RowClaimedByNobody) {
	std::vector<KEnRef_Real_t> counts = wellOwned(32);
	counts[counts.size() - 1] = KEnRef_Real_t(0);

	const OwnershipVerdict v = verifyExactlyOneWriter(counts);
	EXPECT_FALSE(v.ok()) << "an unowned row must be caught; this is the mode the never-committed "
	                        "KENREF_DD_FAULT=0 hook injected";
	EXPECT_EQ(v.unowned, 1u);
	EXPECT_EQ(v.contested, 0u);
	EXPECT_EQ(v.firstBad, counts.size() - 1);
}

// ---------------------------------------------------------------------------------------------
// Beyond the recovered modes: the cases those two do not pin.
// ---------------------------------------------------------------------------------------------

/* firstBad must be the FIRST offending row, not the last one seen. The original loop latched it with a
 * flag; a rewrite that simply assigns on every mismatch would pass both fault-mode tests above (each
 * has exactly one bad row) and fail here. */
TEST(OwnershipCheckTest, ReportsTheFirstOffendingRowNotTheLast) {
	std::vector<KEnRef_Real_t> counts = wellOwned(16);
	counts[4]  = KEnRef_Real_t(0);
	counts[9] += KEnRef_Real_t(1);
	counts[12] = KEnRef_Real_t(0);

	const OwnershipVerdict v = verifyExactlyOneWriter(counts);
	EXPECT_FALSE(v.ok());
	EXPECT_EQ(v.unowned, 2u);
	EXPECT_EQ(v.contested, 1u);
	EXPECT_EQ(v.firstBad, 4u);
}

//! Both kinds at once are counted separately; a run can suffer both, and the diagnostic prints both.
TEST(OwnershipCheckTest, CountsUnownedAndContestedIndependently) {
	std::vector<KEnRef_Real_t> counts = wellOwned(10);
	counts[1] = KEnRef_Real_t(0);
	counts[2] = KEnRef_Real_t(3);   // three ranks claimed it
	counts[7] = KEnRef_Real_t(0);

	const OwnershipVerdict v = verifyExactlyOneWriter(counts);
	EXPECT_EQ(v.unowned, 2u);
	EXPECT_EQ(v.contested, 1u);
	EXPECT_EQ(v.firstBad, 1u);
}

/* Everything unowned is what a set registered but never filled looks like -- the silent failure mode
 * that leaves every restrained atom at the origin. It must not be mistaken for the empty set above. */
TEST(OwnershipCheckTest, CatchesAWhollyUnownedSet) {
	const std::vector<KEnRef_Real_t> counts(8, KEnRef_Real_t(0));
	const OwnershipVerdict v = verifyExactlyOneWriter(counts);
	EXPECT_FALSE(v.ok());
	EXPECT_EQ(v.unowned, 8u);
	EXPECT_EQ(v.firstBad, 0u);
}

/* The counts ride the same floating-point reduction as the coordinates, so they arrive as reals. Small
 * integers are exact in both float and double and a sum of exact 1.0s is exact, which is what makes
 * the equality comparisons legitimate. Pin that rather than leave it as a comment: a fractional count
 * is not a rounding artefact to be tolerated, it means something summed that should not have. */
TEST(OwnershipCheckTest, TreatsAFractionalCountAsContested) {
	std::vector<KEnRef_Real_t> counts = wellOwned(4);
	counts[2] = KEnRef_Real_t(0.5);

	const OwnershipVerdict v = verifyExactlyOneWriter(counts);
	EXPECT_FALSE(v.ok());
	EXPECT_EQ(v.contested, 1u);
	EXPECT_EQ(v.unowned, 0u);
	EXPECT_EQ(v.firstBad, 2u);
}

// ---------------------------------------------------------------------------------------------
// Differential check against the pre-extraction code.
// ---------------------------------------------------------------------------------------------

/* THIS AND THE HAND-BUILT CASES ANSWER DIFFERENT QUESTIONS. KEEP BOTH.
 *
 * The cases above pin the shapes someone thought of, by name, so a failure says which property
 * broke. This one covers the shapes nobody thought of -- degenerate combinations, unusual lengths,
 * interleavings -- and says only "you changed behaviour". Neither subsumes the other: delete the
 * named cases and a failure here tells you nothing about what broke; delete this one and an
 * extraction that is wrong in an unanticipated shape passes.
 *
 * It exists because "the loop moved verbatim" is a claim about intent. Comparing against the
 * original lifted out of git is a claim about behaviour, and that is the one worth making. First run
 * of this comparison, 2026-09-09, was 200,000 vectors with zero mismatches on all three fields and
 * on ok(); it is kept smaller here so the suite stays fast, and the seed is fixed so a failure is
 * reproducible rather than a story about a run nobody can repeat. */
TEST(OwnershipCheckTest, MatchesThePreExtractionLoopOnRandomVectors) {
	std::mt19937 rng(12345);                          // fixed: a failure must be reproducible
	std::uniform_int_distribution<int> len(0, 40);    // includes the empty set
	std::uniform_int_distribution<int> val(0, 4);     // 0 unowned, 1 good, >1 contested

	for (int trial = 0; trial < 20000; ++trial) {
		std::vector<KEnRef_Real_t> counts(static_cast<std::size_t>(len(rng)));
		for (auto &c : counts)
			c = KEnRef_Real_t(val(rng));
		if (!counts.empty() && trial % 7 == 0)
			counts[counts.size() / 2] = KEnRef_Real_t(0.5);   // a fractional count must not slip through

		const OriginalLoopResult   expected = originalLoop(counts);
		const OwnershipVerdict     actual   = verifyExactlyOneWriter(counts);

		ASSERT_EQ(actual.unowned, expected.unowned)     << "trial " << trial << ", len " << counts.size();
		ASSERT_EQ(actual.contested, expected.contested) << "trial " << trial << ", len " << counts.size();
		ASSERT_EQ(actual.firstBad, expected.firstBad)   << "trial " << trial << ", len " << counts.size();
		ASSERT_EQ(actual.ok(), expected.unowned == 0 && expected.contested == 0) << "trial " << trial;
	}
}

//! The raw-pointer overload is what the force provider's vector resolves to; keep them in agreement.
TEST(OwnershipCheckTest, PointerAndVectorOverloadsAgree) {
	std::vector<KEnRef_Real_t> counts = wellOwned(6);
	counts[3] = KEnRef_Real_t(0);

	const OwnershipVerdict a = verifyExactlyOneWriter(counts);
	const OwnershipVerdict b = verifyExactlyOneWriter(counts.data(), counts.size());
	EXPECT_EQ(a.ok(), b.ok());
	EXPECT_EQ(a.unowned, b.unowned);
	EXPECT_EQ(a.contested, b.contested);
	EXPECT_EQ(a.firstBad, b.firstBad);
}
