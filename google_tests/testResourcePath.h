/*
 * testResourcePath.h
 *
 *  Where the test fixtures live, baked in at BUILD time.
 *
 *  Why this exists: every fixture path used to be written "../../res/...", i.e. relative to the
 *  PROCESS WORKING DIRECTORY. That silently requires the test binary to be run from a directory
 *  exactly two levels below the repository root -- true for a cmake-build-X/google_tests directory
 *  inside the repo, and false for everything else. Run the suite from anywhere else and ~14 tests fail with
 *  "cannot open file", which reads exactly like a code regression. That misdiagnosis has already
 *  cost one session a near-miss regression report.
 *
 *  The fix is to stop asking where the process is and ask where the SOURCE is: CMake knows
 *  ${CMAKE_CURRENT_SOURCE_DIR} and passes it down as KENREF_TEST_RES_DIR, so a fixture path no
 *  longer depends on the build directory's position at all. See google_tests/CMakeLists.txt, where
 *  this was already a TODO.
 *
 *  Usage:  KENREF_RES("google_tests/2lum.pdb")  ->  a compile-time string literal, usable in
 *          constexpr context exactly like the raw literal it replaced.
 */

#ifndef KENREF_TESTRESOURCEPATH_H
#define KENREF_TESTRESOURCEPATH_H

#ifndef KENREF_TEST_RES_DIR
/* Deliberately a hard error rather than a "../../res" fallback. A fallback would restore the exact
 * CWD-dependence this header removes, and would do it SILENTLY -- the tests would pass from one
 * directory and fail from another, which is the trap, not a mitigation. If this fires, the
 * target_compile_definitions() in google_tests/CMakeLists.txt was lost. */
#  error "KENREF_TEST_RES_DIR is not defined. google_tests/CMakeLists.txt must pass it; see testResourcePath.h."
#endif

//! Absolute path of a fixture, as a compile-time string literal.
#define KENREF_RES(relative_path) KENREF_TEST_RES_DIR "/" relative_path

//! The fixture root itself, for the few call sites that need the bare directory.
#define KENREF_RES_ROOT KENREF_TEST_RES_DIR "/"

#endif //KENREF_TESTRESOURCEPATH_H
