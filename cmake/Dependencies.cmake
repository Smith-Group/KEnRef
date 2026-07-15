# ============================================================================
# Dependencies.cmake — locate (or fetch) KEnRef's external dependencies WITHOUT
# hardcoding machine-specific paths.
#
# Each dependency is provided by a kenref_provide_<dep>() macro that:
#   1. tries STANDARD discovery (find_package / pkg-config) — honoring the user's
#      CMAKE_PREFIX_PATH / <Pkg>_ROOT / <Pkg>_DIR hints and default system
#      locations (/usr, /usr/local, ...);
#   2. if not found AND KENREF_FETCH_MISSING=ON, downloads & builds it as a
#      sub-build — FetchContent for CMake/header-only deps; ExternalProject for
#      autotools deps such as PLUMED (see the scaffolding at the bottom);
#   3. otherwise fails with an actionable message.
#
# Macros (not functions) so the imported targets they create live in the caller's
# scope. Include this file, then call the kenref_provide_*() you need.
#
# WIRED: Eigen3 (find/fetch), MPI (find, required), OpenMP (find, optional), FFTW3 (detect),
#        GROMACS (find; ExternalProject scaffold), PLUMED (detect; build.sh owns the autotools build).
# ============================================================================
include_guard(GLOBAL)
include(FetchContent)

option(KENREF_FETCH_MISSING "Download & build missing dependencies as sub-builds instead of failing" ON)

# Fetched deps (Eigen, GROMACS) are built with the SAME build type as KEnRef (they must match); default Release.
if(CMAKE_BUILD_TYPE)
    set(KENREF_DEP_BUILD_TYPE "${CMAKE_BUILD_TYPE}")
else()
    set(KENREF_DEP_BUILD_TYPE "Release")
endif()

# Fetched-source download cache: by default CMake's own <build>/_deps (in the build folder — tied to this
# build, cleaned with it). For a shared cache across build dirs, set -DKENREF_DEPS_CACHE=DIR.
set(KENREF_DEPS_CACHE "" CACHE PATH "Shared cache dir for fetched sources (default: <build>/_deps)")
if(KENREF_DEPS_CACHE)
    set(FETCHCONTENT_BASE_DIR "${KENREF_DEPS_CACHE}" CACHE PATH "" FORCE)
endif()
message(STATUS "Fetched-deps cache (FETCHCONTENT_BASE_DIR): ${FETCHCONTENT_BASE_DIR}")

# ----------------------------------------------------------------------------
# Eigen3 (header-only)  ->  imported target Eigen3::Eigen
# Also sets EIGEN3_INCLUDE_DIR (raw path) in the caller's scope: a few KEnRef
# targets consume it directly (kenref_and_eigen3 header-copy, the .pc Cflags).
# ----------------------------------------------------------------------------
set(KENREF_EIGEN_MIN_VERSION "3.4.0"                                  CACHE STRING "Minimum acceptable Eigen3 version")
set(KENREF_EIGEN_GIT_URL     "https://gitlab.com/libeigen/eigen.git"  CACHE STRING "Eigen git URL for the fetch fallback")
set(KENREF_EIGEN_GIT_TAG     "latest"                                 CACHE STRING "Eigen tag to fetch when not found ('latest' = resolve the newest release tag)")
set(KENREF_EIGEN_INSTALL_DIR "${CMAKE_BINARY_DIR}/eigen-install"      CACHE PATH   "Prefix a FETCHED Eigen is INSTALLED to (default: build folder; set e.g. /usr/local or $ENV{HOME}/.local to persist)")

# Resolve the newest Eigen RELEASE tag from the remote (git ls-remote), highest semver. Only used by the
# fetch fallback below; sets _kn_eigen_fetch_tag. 'latest' -> newest; anything else -> use it verbatim.
function(_kenref_resolve_eigen_tag out_var)
    if(NOT KENREF_EIGEN_GIT_TAG STREQUAL "latest")
        set(${out_var} "${KENREF_EIGEN_GIT_TAG}" PARENT_SCOPE)
        return()
    endif()
    find_package(Git QUIET)
    set(_best "")
    if(GIT_EXECUTABLE)
        execute_process(COMMAND "${GIT_EXECUTABLE}" ls-remote --tags --refs "${KENREF_EIGEN_GIT_URL}"
            OUTPUT_VARIABLE _lsr RESULT_VARIABLE _rc ERROR_QUIET TIMEOUT 30)
        if(_rc EQUAL 0 AND _lsr)
            string(REGEX MATCHALL "refs/tags/[0-9]+\\.[0-9]+[.0-9]*" _reftags "${_lsr}")
            set(_vers "")
            foreach(_r IN LISTS _reftags)
                string(REGEX REPLACE "refs/tags/" "" _v "${_r}")
                list(APPEND _vers "${_v}")
            endforeach()
            if(_vers)
                list(SORT _vers COMPARE NATURAL)
                list(GET _vers -1 _best)
            endif()
        endif()
    endif()
    if(NOT _best)
        set(_best "3.4.0")   # network/git unavailable -> last-known-good
        message(WARNING "Eigen3: could not resolve the newest release tag (git ls-remote failed) -> using ${_best}")
    endif()
    set(${out_var} "${_best}" PARENT_SCOPE)
endfunction()

macro(kenref_provide_eigen)
    if(NOT TARGET Eigen3::Eigen)
        # 1) standard discovery (honors CMAKE_PREFIX_PATH, Eigen3_ROOT, Eigen3_DIR, system dirs). NO version
        #    constraint: Eigen 5.x's config uses SameMajorVersion compat and would REJECT a "3.4.0" request;
        #    KEnRef works with Eigen 3.4+ and 5.x alike, so accept whatever is found.
        find_package(Eigen3 CONFIG QUIET)
        if(Eigen3_FOUND)
            message(STATUS "Eigen3: found v${Eigen3_VERSION} (${Eigen3_DIR})  [external — not shipped with kenref]")
            set(KENREF_EIGEN_FETCHED FALSE)
        elseif(KENREF_FETCH_MISSING)
            # 2) fetch Eigen source (download only — SOURCE_SUBDIR cmake/ avoids add_subdirectory of Eigen's
            #    build), then INSTALL it with Eigen's OWN cmake so it becomes a REAL, findable Eigen3
            #    (headers + eigen3.pc + Eigen3Config), not a transient build-interface cache. Install prefix
            #    defaults to the build folder (no sudo); point KENREF_EIGEN_INSTALL_DIR elsewhere to persist.
            _kenref_resolve_eigen_tag(_kn_eigen_fetch_tag)
            message(STATUS "Eigen3: not found -> fetching ${_kn_eigen_fetch_tag} + installing to ${KENREF_EIGEN_INSTALL_DIR}")
            FetchContent_Declare(Eigen3
                GIT_REPOSITORY "${KENREF_EIGEN_GIT_URL}"
                GIT_TAG        "${_kn_eigen_fetch_tag}"
                GIT_SHALLOW    TRUE
                SOURCE_SUBDIR  cmake)
            FetchContent_MakeAvailable(Eigen3)   # populates eigen3_SOURCE_DIR (no configure of Eigen)
            if(NOT EXISTS "${KENREF_EIGEN_INSTALL_DIR}/share/eigen3/cmake/Eigen3Config.cmake")
                # Complete Eigen install = headers + Eigen3Config/Version + eigen3.pc — exactly what distro
                # packages (libeigen3-dev) ship. We turn OFF Eigen's OPTIONAL reference blas/lapack libs: they
                # take minutes to compile, kenref uses Eigen header-only, and distros don't ship them either.
                # (Set -DEIGEN_BUILD_BLAS=ON -DEIGEN_BUILD_LAPACK=ON via the build.sh `--` passthrough if you
                #  really want them.) CMAKE_EXPORT_PACKAGE_REGISTRY=OFF prevents the ~/.cmake/packages pollution.
                # built with KEnRef's build type (${KENREF_DEP_BUILD_TYPE}); same as kenref, never hardcoded.
                execute_process(COMMAND "${CMAKE_COMMAND}" -S "${eigen3_SOURCE_DIR}" -B "${CMAKE_BINARY_DIR}/eigen-build"
                                "-DCMAKE_INSTALL_PREFIX=${KENREF_EIGEN_INSTALL_DIR}" "-DCMAKE_BUILD_TYPE=${KENREF_DEP_BUILD_TYPE}"
                                -DEIGEN_BUILD_BLAS=OFF -DEIGEN_BUILD_LAPACK=OFF -DEIGEN_BUILD_TESTING=OFF
                                -DBUILD_TESTING=OFF -DEIGEN_BUILD_DOC=OFF -DEIGEN_BUILD_DEMOS=OFF
                                -DCMAKE_EXPORT_PACKAGE_REGISTRY=OFF
                                RESULT_VARIABLE _kn_ec)
                if(_kn_ec EQUAL 0)
                    execute_process(COMMAND "${CMAKE_COMMAND}" --install "${CMAKE_BINARY_DIR}/eigen-build" RESULT_VARIABLE _kn_ei)
                endif()
                if(NOT (_kn_ec EQUAL 0 AND _kn_ei EQUAL 0))
                    message(FATAL_ERROR "Eigen3 install failed (see output above).")
                endif()
            endif()
            # find the staged Eigen -> normal 'found' path (Eigen3_DIR straight at the config dir, FORCE past
            # the NOTFOUND the first find_package cached). KENREF_EIGEN_FETCHED=TRUE -> it ships with kenref
            # (KEnRefCore.cmake installs it into the kenref prefix).
            set(Eigen3_DIR "${KENREF_EIGEN_INSTALL_DIR}/share/eigen3/cmake" CACHE PATH "KEnRef-built Eigen3" FORCE)
            find_package(Eigen3 CONFIG REQUIRED NO_CMAKE_PACKAGE_REGISTRY)
            set(KENREF_EIGEN_FETCHED TRUE)
            message(STATUS "Eigen3: built (${KENREF_DEP_BUILD_TYPE}) + staged v${Eigen3_VERSION} (${Eigen3_DIR}) — will ship with kenref")
        else()
            message(FATAL_ERROR
                "Eigen3 (>= ${KENREF_EIGEN_MIN_VERSION}) not found and KENREF_FETCH_MISSING=OFF.\n"
                "  Install it (system package, or into /usr/local), OR pass a hint:\n"
                "    -DCMAKE_PREFIX_PATH=<eigen-prefix>   or   -DEigen3_DIR=<prefix>/share/eigen3/cmake\n"
                "  OR allow the sub-build:   -DKENREF_FETCH_MISSING=ON")
        endif()
    endif()

    if(NOT TARGET Eigen3::Eigen)
        message(FATAL_ERROR "Eigen3::Eigen target missing after find/fetch.")
    endif()

    # Derive EIGEN3_INCLUDE_DIR from the target: unwrap generator expressions, strip a leading -I,
    # and pick the entry that actually contains the Eigen/ header directory.
    get_target_property(_kn_eigen_inc Eigen3::Eigen INTERFACE_INCLUDE_DIRECTORIES)
    if(NOT _kn_eigen_inc)
        set(_kn_eigen_inc "${Eigen3_INCLUDE_DIRS}")
    endif()
    string(REGEX REPLACE "\\$<BUILD_INTERFACE:([^>]+)>" "\\1" _kn_eigen_inc "${_kn_eigen_inc}")
    string(REGEX REPLACE "\\$<INSTALL_INTERFACE:[^>]*>" ""    _kn_eigen_inc "${_kn_eigen_inc}")
    string(REGEX REPLACE "^-I"                          ""    _kn_eigen_inc "${_kn_eigen_inc}")
    foreach(_p IN LISTS _kn_eigen_inc)
        if(_p AND EXISTS "${_p}/Eigen")
            set(_kn_eigen_inc "${_p}")
            break()
        endif()
    endforeach()
    set(EIGEN3_INCLUDE_DIR "${_kn_eigen_inc}")
    unset(_kn_eigen_inc)
    unset(_p)
    message(STATUS "Eigen3 include dir: ${EIGEN3_INCLUDE_DIR}")
endmacro()

# ----------------------------------------------------------------------------
# MPI  ->  imported targets MPI::MPI_CXX / MPI::MPI_C  (REQUIRED)
# KEnRef's multi-sim gather/scatter needs it. No fetch (MPI is a system / HPC-module dependency).
# Compiler-agnostic: find_package(MPI) works with GCC or Clang; it detects the mpicxx/mpicc wrappers
# on PATH (or via MPI_HOME / MPI_C_COMPILER / CMAKE_PREFIX_PATH). The chosen COMPILER (Clang here) lives
# in the dev env (~/scripts/*.sh, CLion) — NOT hardcoded in this CMake.
# ----------------------------------------------------------------------------
macro(kenref_provide_mpi)
    if(NOT TARGET MPI::MPI_CXX)
        find_package(MPI QUIET COMPONENTS C CXX)
        if(MPI_CXX_FOUND)
            message(STATUS "MPI: found ${MPI_CXX_COMPILER} (v${MPI_CXX_VERSION})")
        else()
            message(FATAL_ERROR
                "MPI (C/CXX) not found, but KEnRef needs it for multi-sim gather/scatter.\n"
                "  Install OpenMPI/MPICH (or `module load` your HPC MPI), OR hint with\n"
                "  -DMPI_HOME=<prefix> / put mpicxx,mpicc on PATH / -DCMAKE_PREFIX_PATH=<mpi-prefix>.")
        endif()
    endif()
endmacro()

# ----------------------------------------------------------------------------
# OpenMP  ->  imported target OpenMP::OpenMP_CXX  (OPTIONAL — soft dependency)
# The hot kernels use OpenMP pragmas; without OpenMP they compile to serial. find_package supplies the
# compiler-correct flag automatically (-fopenmp for GCC, -fopenmp=libomp for Clang) — again, compiler-neutral.
# ----------------------------------------------------------------------------
macro(kenref_provide_openmp)
    if(NOT TARGET OpenMP::OpenMP_CXX)
        find_package(OpenMP QUIET COMPONENTS CXX)
        if(OpenMP_CXX_FOUND)
            message(STATUS "OpenMP: found (${OpenMP_CXX_FLAGS})")
        else()
            message(STATUS "OpenMP: NOT found -> KEnRef kernels compile to serial (pragmas ignored).")
        endif()
    endif()
endmacro()

# ----------------------------------------------------------------------------
# FFTW3 — KEnRef does NOT link it directly; only the GROMACS/PLUMED sub-builds need it. This just
# DETECTS a system FFTW (pkg-config fftw3) and sets KENREF_FFTW_FOUND, so a sub-build can be pointed at
# it; if absent, the GROMACS sub-build can build its own (-DGMX_BUILD_OWN_FFTW=ON) and PLUMED treats FFT
# as optional. (No standalone-FFTW fetch here — the sub-builds handle that.)
# ----------------------------------------------------------------------------
macro(kenref_provide_fftw)
    set(KENREF_FFTW_FOUND FALSE)
    find_package(PkgConfig QUIET)
    if(PkgConfig_FOUND)
        pkg_check_modules(KN_FFTW3 QUIET fftw3)
        if(KN_FFTW3_FOUND)
            set(KENREF_FFTW_FOUND TRUE)
            message(STATUS "FFTW3: found v${KN_FFTW3_VERSION} (${KN_FFTW3_LIBDIR})")
        endif()
    endif()
    if(NOT KENREF_FFTW_FOUND)
        message(STATUS "FFTW3: not found -> GROMACS/PLUMED sub-builds provide their own (GROMACS: -DGMX_BUILD_OWN_FFTW=ON).")
    endif()
endmacro()

# ----------------------------------------------------------------------------
# GROMACS — the KEnRef GMX interface is a PATCHED gmx mdrun: it needs (a) the installed CMake package
# (find_package(GROMACS) at <install>/share/cmake/gromacs[_mpi]) AND (b) the gromacs SOURCE + BUILD tree
# (internal headers). This is STOCK (upstream) gromacs — KEnRef does NOT fork or patch it; it depends on
# gromacs's source + build + library and drives it with its OWN classes (KEnRefForceProvider, the gmx
# wrapper). Discovery honours the caller-set GROMACS_DIR (from GROMACS_INSTALL_DIR) / CMAKE_PREFIX_PATH —
# concrete paths come from your dev env (~/scripts, CLion), not hardcoded here. Fetch fallback =
# ExternalProject clone+build of upstream gromacs (heavy; scaffolded for now).
# ----------------------------------------------------------------------------
# The kenref-gmx force provider's GROMACS (stock gromacs it links) — its download config. The PLUMED-batched
# GROMACS is a separate consumer, handled on the PLUMED side (build-and-batch.sh), not here. build.sh's
# --with-gmx does the actual fetch+build; these are the defaults it uses.
set(KENREF_GMX_GROMACS_GIT_URL   "https://gitlab.com/gromacs/gromacs.git" CACHE STRING
    "Upstream GROMACS git URL for the kenref-gmx force-provider fetch (stock gromacs)")
# The tag to fetch. Pinned to a known-good release that KEnRef is developed/tested against rather than a moving
# branch, so an auto-fetch is reproducible. (The force provider does not version-gate — override at will.)
set(KENREF_GMX_GROMACS_GIT_TAG   "v2025.4"                                CACHE STRING
    "GROMACS git tag/branch to fetch for the kenref-gmx force provider")
# Fetch INTO the kenref build dir (tied to this build; no cross-version collision with a system gromacs).
set(KENREF_GMX_GROMACS_FETCH_DIR "${CMAKE_BINARY_DIR}/gromacs-src"        CACHE PATH
    "Where to fetch GROMACS for the kenref-gmx force provider (inside the kenref build dir)")

# FIND a usable stock GROMACS (install cmake package + source + build tree — the force provider needs all
# three), OR FETCH the source and DELEGATE its compilation to GROMACS's OWN cmake, then use that. All the
# gromacs logic lives HERE (build manager), not in build.sh; build.sh only sets -DBUILD_KENREF_GMX=ON.
# Sets GROMACS_SRC_DIR / GROMACS_BUILD_DIR / GROMACS_INSTALL_DIR / GROMACS_DIR for KEnRefGMX.cmake.
macro(kenref_provide_gromacs)
    # "provided" = the GROMACS install CMake package is FINDABLE (find_package succeeds). That is the signal
    # to USE it — never fetch when a gromacs is already provided (a missing SRC/BUILD dir is a user setup
    # error the gmx build reports clearly, NOT a reason to silently pull a different gromacs). Only fetch when
    # no install is found at all. This keeps a CLion/dev reconfigure from ever triggering a 20-min fetch.
    find_package(GROMACS NAMES gromacs_mpi gromacs QUIET)
    if(GROMACS_FOUND)
        message(STATUS "GROMACS: using provided ${GROMACS_VERSION} (${GROMACS_CONFIG}); src=${GROMACS_SRC_DIR} build=${GROMACS_BUILD_DIR}")
    elseif(KENREF_FETCH_MISSING OR (DEFINED GROMACS_SRC_DIR AND EXISTS "${GROMACS_SRC_DIR}/CMakeLists.txt"))
        # Build GROMACS ourselves. Two ways in:
        #   * the user gave -DGROMACS_SRC_DIR=<their checkout>  -> build THAT (no download, honoured even when
        #     KENREF_FETCH_MISSING=OFF: building a source you supplied is not a download).
        #   * nothing given + downloads enabled                 -> clone a stock GROMACS first.
        # (Previously this branch unconditionally overwrote GROMACS_SRC_DIR with the fetch dir, silently
        #  discarding a user-provided checkout and cloning a different GROMACS.)
        if(NOT (DEFINED GROMACS_SRC_DIR AND EXISTS "${GROMACS_SRC_DIR}/CMakeLists.txt"))
            find_package(Git REQUIRED)
            set(GROMACS_SRC_DIR "${KENREF_GMX_GROMACS_FETCH_DIR}")
            if(NOT EXISTS "${GROMACS_SRC_DIR}/CMakeLists.txt")
                message(STATUS "GROMACS: fetching ${KENREF_GMX_GROMACS_GIT_TAG} -> ${GROMACS_SRC_DIR}")
                execute_process(COMMAND "${GIT_EXECUTABLE}" clone --depth 1 --branch "${KENREF_GMX_GROMACS_GIT_TAG}"
                                "${KENREF_GMX_GROMACS_GIT_URL}" "${GROMACS_SRC_DIR}" RESULT_VARIABLE _kg_rc)
                if(NOT _kg_rc EQUAL 0)
                    message(FATAL_ERROR "GROMACS clone failed. Provide -DGROMACS_BUILD_DIR=<an existing gromacs "
                                        "build dir> or -DGROMACS_SRC_DIR=<a gromacs source checkout>.")
                endif()
            endif()
        else()
            message(STATUS "GROMACS: building the PROVIDED source -> ${GROMACS_SRC_DIR} (no download)")
        endif()
        set(GROMACS_BUILD_DIR   "${CMAKE_BINARY_DIR}/gromacs-build")
        set(GROMACS_INSTALL_DIR "${CMAKE_BINARY_DIR}/gromacs-install")
        # FORCE the cache past the NOTFOUND the QUIET find above cached, so the final find_package succeeds.
        set(GROMACS_DIR         "${CMAKE_BINARY_DIR}/gromacs-install/share/cmake/gromacs_mpi" CACHE PATH "KEnRef-built GROMACS" FORCE)
        if(NOT EXISTS "${GROMACS_DIR}")   # not built+installed yet -> delegate the build to gromacs's own cmake
            set(_kg_simd "")
            if(DEFINED ACCEL AND NOT ACCEL STREQUAL "")
                set(_kg_simd "-DGMX_SIMD=${ACCEL}")
            endif()
            message(STATUS "GROMACS: building via its own cmake -> ${GROMACS_INSTALL_DIR}  (one-time; can take 20+ min)")
            # no -G: gromacs builds with the default generator (or $CMAKE_GENERATOR); Make and Ninja both work.
            # build type matches KEnRef (${KENREF_DEP_BUILD_TYPE}), never hardcoded.
            # GMX_INSTALL_LEGACY_API=ON is REQUIRED: gromacs defaults it OFF and then does NOT install the
            # public headers kenref-gmx compiles against (gromacs/math/vectypes.h, gromacs/topology/index.h, …).
            execute_process(COMMAND "${CMAKE_COMMAND}" -S "${GROMACS_SRC_DIR}" -B "${GROMACS_BUILD_DIR}"
                            "-DCMAKE_BUILD_TYPE=${KENREF_DEP_BUILD_TYPE}" -DGMX_MPI=ON -DGMX_INSTALL_LEGACY_API=ON
                            ${_kg_simd} "-DCMAKE_INSTALL_PREFIX=${GROMACS_INSTALL_DIR}"
                            RESULT_VARIABLE _kg_c)
            if(_kg_c EQUAL 0)
                execute_process(COMMAND "${CMAKE_COMMAND}" --build "${GROMACS_BUILD_DIR}" RESULT_VARIABLE _kg_b)
            endif()
            if(_kg_c EQUAL 0 AND _kg_b EQUAL 0)
                execute_process(COMMAND "${CMAKE_COMMAND}" --install "${GROMACS_BUILD_DIR}" RESULT_VARIABLE _kg_i)
            endif()
            if(NOT (_kg_c EQUAL 0 AND _kg_b EQUAL 0 AND _kg_i EQUAL 0))
                message(FATAL_ERROR "GROMACS build/install failed — see the output above.")
            endif()
        endif()
    else()
        message(FATAL_ERROR
            "GROMACS not found, and downloads are disabled.\n"
            "  Point KEnRef at an EXISTING gromacs build:  -DGROMACS_BUILD_DIR=<gromacs's cmake build dir>\n"
            "    (the dir you ran `cmake -B` into for GROMACS — it has CMakeCache.txt and src/include/config.h;\n"
            "     NOT the install prefix, NOT the source. Its source and install are derived from it.)\n"
            "  OR have KEnRef build YOUR gromacs source:    -DGROMACS_SRC_DIR=<gromacs source checkout>\n"
            "  OR allow a fetch+build of a stock gromacs:   -DKENREF_FETCH_MISSING=ON  (build.sh --download ON)")
    endif()
    find_package(GROMACS NAMES gromacs_mpi gromacs REQUIRED)
    message(STATUS "GROMACS: ${GROMACS_VERSION} (config ${GROMACS_CONFIG}); src=${GROMACS_SRC_DIR} build=${GROMACS_BUILD_DIR}")
endmacro()

