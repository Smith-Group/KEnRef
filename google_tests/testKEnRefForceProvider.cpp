/*
 * testKEnRefForceProvider.cpp
 *
 * The previous KEnRefForceProviderTestSuite.TestRestoreNoJump migrated to the core test suite
 * (google_tests/testRestoreNoJump.cpp) when restoreNoJump moved into kenref_core during the
 * model-abstraction restructure: the GROMACS force provider no longer has a restoreNoJump member
 * (it is now a thin kenref::EngineAdapter driving a kenref::KEnRefDriver). New GROMACS-specific
 * force-provider tests belong here.
 */

#include "gtest/gtest.h"
