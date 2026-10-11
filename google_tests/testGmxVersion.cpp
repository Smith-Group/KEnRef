/*
 * testGmxVersion.cpp
 *
 * kenref::gmxVersionNumber reads the release out of gmx_version()'s string so that main() can refuse to run
 * inside a GROMACS other than the one KEnRef was compiled against (KEnRefGmxVersion.h has why). The strings
 * below are the formats GROMACS actually produces; the last test runs the real comparison main() makes, in
 * this binary, against the libgromacs it loaded.
 */
#include "gmxpre.h"
#include "gromacs/version.h"
#include "gromacs/utility/baseversion.h"
#include <gtest/gtest.h>
#include "gmxinterface/KEnRefGmxVersion.h"

using kenref::gmxVersionNumber;

TEST(GmxVersionTestSuite, readsAReleaseInTheEncodingOfGmxVersion) {
    EXPECT_EQ(gmxVersionNumber("2025.4"), 20250004);
    EXPECT_EQ(gmxVersionNumber("2026.3"), 20260003);
    EXPECT_EQ(gmxVersionNumber("2025.10"), 20250010);
}

TEST(GmxVersionTestSuite, ignoresWhatABuildFromGitAppends) {
    // what GMX-DEST/{2025.4,2026.3} report: built from git checkouts of the release tags
    EXPECT_EQ(gmxVersionNumber("2025.4-dev-20251121-560f2601a7"), 20250004);
    EXPECT_EQ(gmxVersionNumber("2026.3-dev-20260625-1210900145"), 20260003);
}

TEST(GmxVersionTestSuite, readsADevelopmentBranchAsPatchZero) {
    // main between releases: no ".N", and GMX_VERSION is year * 10000
    EXPECT_EQ(gmxVersionNumber("2027-dev-20261005-abcdef0123"), 20270000);
    EXPECT_EQ(gmxVersionNumber("2027"), 20270000);
}

TEST(GmxVersionTestSuite, distinguishesReleasesThatShareASoname) {
    // 2025.4 and 2025.5 are both libgromacs_mpi.so.10: the case the loader lets through
    EXPECT_NE(gmxVersionNumber("2025.5-dev-20260101-0123456789"), gmxVersionNumber("2025.4-dev-20251121-560f2601a7"));
}

TEST(GmxVersionTestSuite, reportsAnUnreadableStringAsUnreadable) {
    // never as some other release: -1 matches no GMX_VERSION, so main() refuses
    EXPECT_EQ(gmxVersionNumber(""), -1);
    EXPECT_EQ(gmxVersionNumber("VERSION 2026.3"), -1);
    EXPECT_EQ(gmxVersionNumber("v2026.3"), -1);
    EXPECT_EQ(gmxVersionNumber("226.3"), -1);
    EXPECT_EQ(gmxVersionNumber("20263"), -1);
    EXPECT_EQ(gmxVersionNumber("2026."), -1);
    EXPECT_EQ(gmxVersionNumber("2026.x"), -1);
}

TEST(GmxVersionTestSuite, theLoadedLibgromacsIsTheOneThisWasCompiledAgainst) {
    // main()'s check, made here against the real library; it fails if this binary picked up another GROMACS
    EXPECT_EQ(gmxVersionNumber(gmx_version()), GMX_VERSION) << "gmx_version() = '" << gmx_version() << "'";
}
