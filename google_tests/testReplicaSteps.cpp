/* The decision behind refusing a multi-simulation whose replicas would run different numbers of steps.
 *
 * KEnRefReplicaSteps.h explains why: KEnRef's per-step exchange between replicas hangs when one
 * finishes first, and GROMACS only warns. These tests drive the pure part with hand-built values; the
 * collective and the tpr reading are exercised by real 2-replica -multidir runs. */
#include <cstdint>
#include <string>
#include <vector>

#include "gtest/gtest.h"
#include "gmxinterface/KEnRefReplicaSteps.h"

using kenref::checkReplicaSteps;
using kenref::kNstepsFromTpr;
using kenref::replicaStepKey;

// Without -nsteps, the key is the tpr's last step, which a checkpoint restart does not change.
TEST(ReplicaStepKey, IsTheLastStepFromTheTpr) {
	EXPECT_EQ(replicaStepKey(0, 20, kNstepsFromTpr), 20);
	EXPECT_EQ(replicaStepKey(100, 20, kNstepsFromTpr), 120);
	// What GROMACS does on restart: nsteps += init_step - cpt, init_step = cpt. The key must not move.
	const std::int64_t init = 0, nsteps = 20, cpt = 10;
	EXPECT_EQ(replicaStepKey(cpt, nsteps + init - cpt, kNstepsFromTpr), replicaStepKey(init, nsteps, kNstepsFromTpr));
}

// -nsteps overrides the tpr after the checkpoint is applied, so the key is the step count itself.
TEST(ReplicaStepKey, CommandLineOverridesTheTpr) {
	EXPECT_EQ(replicaStepKey(0, 20, 5), 5);
	EXPECT_EQ(replicaStepKey(100, 20, 0), 0);
	EXPECT_EQ(replicaStepKey(0, 20, -1), -1);   // -nsteps -1 is unlimited
}

TEST(ReplicaStepKey, UnlimitedIsMinusOne) {
	EXPECT_EQ(replicaStepKey(0, -1, kNstepsFromTpr), -1);
	EXPECT_EQ(replicaStepKey(50, -1, kNstepsFromTpr), -1);   // not 49: init-step does not apply
}

TEST(CheckReplicaSteps, EqualKeysAreConsistent) {
	const auto v = checkReplicaSteps({20, 20, 20}, false, {"a", "b", "c"});
	EXPECT_TRUE(v.consistent);
	EXPECT_TRUE(v.report.empty());
	EXPECT_TRUE(checkReplicaSteps({-1, -1}, false, {}).consistent);
}

TEST(CheckReplicaSteps, FewerThanTwoReplicasAreConsistent) {
	EXPECT_TRUE(checkReplicaSteps({}, false, {}).consistent);
	EXPECT_TRUE(checkReplicaSteps({7}, false, {}).consistent);
}

// The case that hangs: nsteps 20 in one replica's tpr and 10 in the other's.
TEST(CheckReplicaSteps, MismatchIsRefusedNamingEveryReplica) {
	const auto v = checkReplicaSteps({20, 10}, false, {"rep0", "rep1"});
	EXPECT_FALSE(v.consistent);
	EXPECT_NE(v.report.find("would not run the same number of steps"), std::string::npos);
	EXPECT_NE(v.report.find("Last step of each replica"), std::string::npos);
	EXPECT_NE(v.report.find("replica 0 (rep0): 20\n"), std::string::npos);
	EXPECT_NE(v.report.find("replica 1 (rep1): 10\n"), std::string::npos);
}

// One replica differing, wherever it sits: the first replica is not special.
TEST(CheckReplicaSteps, OneOddReplicaAnywhereIsRefused) {
	EXPECT_FALSE(checkReplicaSteps({20, 20, 20, 10}, false, {}).consistent);
	EXPECT_FALSE(checkReplicaSteps({10, 20, 20, 20}, false, {}).consistent);
	EXPECT_FALSE(checkReplicaSteps({20, 10, 20}, false, {}).consistent);
}

TEST(CheckReplicaSteps, UnlimitedAgainstFiniteIsRefused) {
	const auto v = checkReplicaSteps({-1, 20}, false, {});
	EXPECT_FALSE(v.consistent);
	EXPECT_NE(v.report.find("replica 0: unlimited (nsteps = -1)\n"), std::string::npos);
	EXPECT_NE(v.report.find("replica 1: 20\n"), std::string::npos);
}

TEST(CheckReplicaSteps, ReportSaysWhereTheValuesCameFrom) {
	const auto v = checkReplicaSteps({5, 6}, true, {});
	EXPECT_NE(v.report.find("Steps given by mdrun -nsteps"), std::string::npos);
	EXPECT_EQ(v.report.find("Last step of each replica"), std::string::npos);
}

// Missing names fall back to the index rather than reading past the end.
TEST(CheckReplicaSteps, MissingNamesAreLabelledByIndex) {
	const auto v = checkReplicaSteps({1, 2, 3}, false, {"only"});
	EXPECT_NE(v.report.find("replica 0 (only): 1\n"), std::string::npos);
	EXPECT_NE(v.report.find("replica 1: 2\n"), std::string::npos);
	EXPECT_NE(v.report.find("replica 2: 3\n"), std::string::npos);
}
