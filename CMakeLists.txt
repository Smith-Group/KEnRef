cmake_minimum_required(VERSION 3.24)
project(KEnRef VERSION 2.0.0)
# Single source of truth for the KEnRef version. Flows to the pkg-config files (@PROJECT_VERSION@ in
# cmake/*.pc.in) and, via those, to PLUMED's configure floor (KENREF_MIN_VERSION in plumed configure.ac,
# checked as `kenref_core >= X`). Bump here when the public API / numerical contract changes; keep the
# PLUMED floor in lockstep. v2.0.0 == the model-abstraction restructure (branch v2.0 and its subbranches).
set(CMAKE_CXX_STANDARD 17)

# Emit compile_commands.json so clangd / IDEs resolve includes for every translation unit (including
# newly added core/model files). FORCE so it wins over a stale cache value. Harmless for the build.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE BOOL "Export compile_commands.json for clangd/IDEs" FORCE)

string(TOLOWER ${CMAKE_BUILD_TYPE} build_type)
message(STATUS "the build type is: ${build_type}")

# ============================================================================
# BUILD OPTIONS
# ============================================================================
option(BUILD_KENREF_CORE "Build KENREF_CORE library" ON)
option(BUILD_KENREF_GMX "Build KENREF GMX interface" ON)
# The PLUMED interface is entirely OPT-IN — a plain kenref_core build neither COMPILES nor EXPORTS it.
# (In production the PLUMED fork compiles src/plumedinterface itself, Option A.)
#   * BUILD_KENREF_PLUMED          — COMPILE-VALIDATE plumedinterface against PLUMED's headers (a static
#                                    check lib; needs -DPLUMED_SRC_DIR). Default OFF: only when asked.
#   * KENREF_EXPORT_PLUMEDINTERFACE— INSTALL the plumedinterface headers + src/plumedinterface/*.cpp +
#                                    kenref_plumed.pc (the source a PLUMED build compiles in-tree, Option A).
#                                    Default OFF; forced ON below when a PLUMED build here needs it.
option(BUILD_KENREF_PLUMED "Compile-validate the KENREF PLUMED interface (needs -DPLUMED_SRC_DIR)" OFF)
option(KENREF_EXPORT_PLUMEDINTERFACE "Install the plumedinterface headers/source + kenref_plumed.pc (for a PLUMED build)" OFF)

# PLUMED: KEnRef does not REIMPLEMENT plumed's build; it DELEGATES to the plumed side (build-only.sh /
# build-and-batch.sh), and CMake itself invokes that script (via install(CODE) — after kenref_core is
# installed, so the plumed script reuses it). See cmake/InvokePlumed.cmake. build.sh just sets these:
option(KENREF_WITH_PLUMED          "Also build PLUMED (kenref module), delegated to the plumed side at install time" OFF)
option(KENREF_WITH_PLUMED_GROMACS  "Also build PLUMED + a batched GROMACS (delegates to build-and-batch.sh)"        OFF)
if(KENREF_WITH_PLUMED_GROMACS)
	set(KENREF_WITH_PLUMED ON CACHE BOOL "implied by KENREF_WITH_PLUMED_GROMACS" FORCE)
endif()
if(KENREF_WITH_PLUMED)
	# the delegated plumed build compiles plumedinterface from the kenref install -> it must be exported.
	set(KENREF_EXPORT_PLUMEDINTERFACE ON CACHE BOOL "required by KENREF_WITH_PLUMED" FORCE)
endif()

# Per-model compile-time enable switches (Style-2 registration). Default ON; turning one OFF drops that
# model's registration from bootstrapModels() so it is unavailable at runtime. Passed to the core
# targets as KENREF_ENABLE_<MODEL> compile definitions (see cmake/KEnRefCore.cmake).
option(KENREF_ENABLE_SIGMA    "Enable the SIGMA energy model"    ON)
option(KENREF_ENABLE_PLATEAUS "Enable the PLATEAUS energy model" ON)
option(KENREF_ENABLE_RELAX    "Enable the RELAX energy model"    ON)

### comment the next line to make it float. Uncomment to make it double.
add_compile_options(-DKENREF_DOUBLE)

if (NOT DEFINED ACCEL)
	# Auto-detect this machine's SIMD (like GROMACS): pick the highest of AVX_512 / AVX2_256 / AVX_256.
	# Override with -DACCEL=... for a fixed/portable build. Linux reads /proc/cpuinfo; elsewhere -> AVX2_256.
	set(ACCEL "AVX2_256")   # portable fallback
	if(EXISTS "/proc/cpuinfo")
		file(READ "/proc/cpuinfo" _kn_cpuinfo LIMIT 16384)
		if(_kn_cpuinfo MATCHES "avx512f")
			set(ACCEL "AVX_512")
		elseif(_kn_cpuinfo MATCHES "avx2")
			set(ACCEL "AVX2_256")
		elseif(_kn_cpuinfo MATCHES "[ \t]avx[ \t]")
			set(ACCEL "AVX_256")
		endif()
		unset(_kn_cpuinfo)
	endif()
	message(STATUS "ACCEL not set -> auto-detected ${ACCEL} for this machine (override with -DACCEL=...)")
else ()
	message(STATUS "ACCEL variable is: ${ACCEL}")
endif ()

#Important to fix incompatibility bug with PLUMED
if(ACCEL STREQUAL "AVX_512")
	add_compile_options(-march=skylake-avx512)
elseif(ACCEL STREQUAL "AVX2_256")
	add_compile_options(-march=haswell)
elseif(ACCEL STREQUAL "AVX_256")
	add_compile_options(-march=sandybridge)
endif()
# ============================================================================
# GROMACS PATHS (only needed for GMX interface)
# ============================================================================
if(BUILD_KENREF_GMX)
	set(SUFFIX_MPI _mpi)

	# Accept paths as CMake variables or environment variables with fallback defaults
	if(NOT DEFINED GROMACS_SRC_DIR)
		if(DEFINED ENV{GROMACS_SRC_DIR})
			set(GROMACS_SRC_DIR $ENV{GROMACS_SRC_DIR})
		else()
			set(GROMACS_SRC_DIR "${CMAKE_SOURCE_DIR}/../gromacs")
		endif()
	endif()

	if(NOT DEFINED GROMACS_BUILD_DIR)
		if(DEFINED ENV{GROMACS_BUILD_DIR})
			set(GROMACS_BUILD_DIR $ENV{GROMACS_BUILD_DIR})
		else()
			set(GROMACS_BUILD_DIR "${GROMACS_SRC_DIR}/cmake-build-${build_type}")
		endif()
	endif()

	if(NOT DEFINED GROMACS_INSTALL_DIR)
		if(DEFINED ENV{GROMACS_INSTALL_DIR})
			set(GROMACS_INSTALL_DIR $ENV{GROMACS_INSTALL_DIR})
		else()
			# Standard default (no machine-specific path). The KEnRef GMX interface depends on a STOCK
			# gromacs (its source + build + library) via KEnRef's own classes — not a fork. In practice
			# set GROMACS_INSTALL_DIR/SRC_DIR/BUILD_DIR in your dev env (~/scripts, CLion) to your gromacs
			# source+build (they already do).
			set(GROMACS_INSTALL_DIR "/usr/local/gromacs")
		endif()
	endif()

	message(STATUS "GROMACS_SRC_DIR: ${GROMACS_SRC_DIR}")
	message(STATUS "GROMACS_BUILD_DIR: ${GROMACS_BUILD_DIR}")
	message(STATUS "GROMACS_INSTALL_DIR: ${GROMACS_INSTALL_DIR}")
endif()

set(GROMACS_DIR ${GROMACS_INSTALL_DIR}/share/cmake/gromacs${SUFFIX_MPI}) #where libgromacs*.cmake are located

# ============================================================================
# INSTALLATION PATHS
# ============================================================================
# Default install prefix = <system folder>/kenref — flat, like GROMACS's /usr/local/gromacs (no build_type or
# ACCEL in the path; build type defaults to Release, ACCEL is auto-detected above). KENREF_INSTALL_BASE is the
# overridable system folder (default /usr/local); or pass -DCMAKE_INSTALL_PREFIX to set the full path directly.
# Installing under /usr/local needs write permission (sudo) — build.sh handles that; for a no-sudo build pass
# -DCMAKE_INSTALL_PREFIX=<somewhere-writable> or -DKENREF_INSTALL_BASE=$HOME/.local.
set(KENREF_INSTALL_BASE "/usr/local" CACHE PATH "System folder for KEnRef's default install prefix")
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
	if(BUILD_KENREF_GMX)
		set(CMAKE_INSTALL_PREFIX "${KENREF_INSTALL_BASE}/kenref-gmx" CACHE PATH "root directory to install KEnRef GMX" FORCE)
	else()
		set(CMAKE_INSTALL_PREFIX "${KENREF_INSTALL_BASE}/kenref" CACHE PATH "root directory to install KEnRef CORE" FORCE)
	endif()
endif()

# ONE set of core libraries: the kenref-gmx build must CONSUME an installed kenref_core, never install a
# duplicate into the gmx prefix. Building CORE+GMX together puts both under a single prefix (CMake has one
# CMAKE_INSTALL_PREFIX per build) -> the core libs would land in kenref-gmx too. Warn and point at the
# two-phase flow (build.sh --with-gmx does it: core -> KENREF_INSTALL_BASE/kenref, executables -> .../kenref-gmx).
if(BUILD_KENREF_CORE AND BUILD_KENREF_GMX)
	message(WARNING
		"Building CORE and GMX in ONE build installs BOTH under a single prefix (${CMAKE_INSTALL_PREFIX}), "
		"duplicating libkenref_core there. For one canonical copy of the libraries, build in two phases: "
		"(1) -DBUILD_KENREF_CORE=ON -DBUILD_KENREF_GMX=OFF  -> installs libs to ${KENREF_INSTALL_BASE}/kenref ; "
		"(2) -DBUILD_KENREF_CORE=OFF -DBUILD_KENREF_GMX=ON  -> installs executables to ${KENREF_INSTALL_BASE}/kenref-gmx "
		"(consuming the installed core). `build.sh --with-gmx` does this automatically.")
endif()

#set(CMAKE_MODULE_PATH /usr/local/gromacs-dev/share/gromacs/template/cmake:$CMAKE_MODULE_PATH)

# Dependency discovery (standard find -> fetch fallback; no hardcoded machine paths).
# Pass hints via -DCMAKE_PREFIX_PATH=... / -DEigen3_DIR=... when a dep lives outside default locations.
include(cmake/Dependencies.cmake)

# ============================================================================
# COMMON DEPENDENCIES
# ============================================================================
kenref_provide_mpi()   # MPI::MPI_CXX / MPI::MPI_C (required)

if(BUILD_KENREF_CORE)
	kenref_provide_eigen()   # find Eigen3 (or fetch it if missing); sets Eigen3::Eigen + EIGEN3_INCLUDE_DIR
endif()

if(BUILD_KENREF_GMX)
	message(STATUS "Searching for Gromacs in ${GROMACS_DIR}")
	kenref_provide_fftw()      # detect system FFTW (informational; the gromacs sub-build can build its own)
	kenref_provide_gromacs()   # find_package(GROMACS) honoring GROMACS_DIR; fetch scaffold if missing
endif()

if(BUILD_KENREF_CORE)
	kenref_provide_openmp()   # optional: OpenMP::OpenMP_CXX if present, else kernels compile to serial
endif()

# ============================================================================
# INCLUDE COMPONENT MODULES
# ============================================================================

if(BUILD_KENREF_CORE)
	include(cmake/KEnRefCore.cmake)
endif()

if(BUILD_KENREF_GMX)
	# If building GMX without CORE, need to find installed CORE
	if(NOT BUILD_KENREF_CORE AND NOT TARGET KENREF::CORE)
		message(STATUS "Searching for installed KENREF CORE...")
		# Mirror the core install default (KENREF_INSTALL_BASE/kenref); override with
		# -DKENREF_CORE_CMAKE_PATH=<prefix>/lib/cmake/KEnRef or -DCMAKE_PREFIX_PATH=<core-prefix>.
		if(NOT DEFINED KENREF_CORE_CMAKE_PATH)
			set(KENREF_CORE_CMAKE_PATH "${KENREF_INSTALL_BASE}/kenref/lib/cmake/KEnRef")
		endif()
		message(STATUS "  Looking in: ${KENREF_CORE_CMAKE_PATH}")

		# KENREF CORE depends on Eigen3, so provide it first (find or fetch).
		kenref_provide_eigen()

		if(EXISTS "${KENREF_CORE_CMAKE_PATH}/KEnRefCoreConfig.cmake")
			message(STATUS "  Found KEnRefCoreConfig.cmake, including it...")
			include("${KENREF_CORE_CMAKE_PATH}/KEnRefCoreConfig.cmake")

			# The exported target is KENREF::kenref_core, create an alias to KENREF::CORE
			if(TARGET KENREF::kenref_core AND NOT TARGET KENREF::CORE)
				add_library(KENREF::CORE ALIAS KENREF::kenref_core)
			endif()

			message(STATUS "  ✓ Loaded KENREF CORE from ${KENREF_CORE_CMAKE_PATH}")
		else()
			message(FATAL_ERROR "KENREF CORE not found at ${KENREF_CORE_CMAKE_PATH}\nPlease build and install KENREF CORE first, or build with -DBUILD_KENREF_CORE=ON")
		endif()
	endif()

	include(cmake/KEnRefGMX.cmake)
endif()

if(BUILD_KENREF_PLUMED)
	# Needs KENREF::CORE (build it here, or it must be found as above for the GMX path).
	if(NOT TARGET KENREF::CORE)
		message(FATAL_ERROR "BUILD_KENREF_PLUMED requires KENREF::CORE (build with -DBUILD_KENREF_CORE=ON).")
	endif()
	include(cmake/KEnRefPlumed.cmake)
endif()

# Post-install: env.sh + TCL modulefile (always) + a clear summary of where things went.
if(BUILD_KENREF_CORE OR BUILD_KENREF_GMX)
	include(cmake/PostInstall.cmake)
endif()

# PLUMED delegation, invoked by CMake at install time (AFTER PostInstall so env.sh is installed first).
if(KENREF_WITH_PLUMED)
	if(NOT BUILD_KENREF_CORE)
		message(FATAL_ERROR "KENREF_WITH_PLUMED needs kenref_core (build with -DBUILD_KENREF_CORE=ON).")
	endif()
	include(cmake/InvokePlumed.cmake)
endif()

# ============================================================================
# TESTING
# ============================================================================
option(BUILD_TESTS "Build Google tests" OFF)
if(BUILD_TESTS)
	enable_testing()
	add_subdirectory(google_tests)
endif()

# ============================================================================
# PKG-CONFIG GENERATION FOR PLUMED INTEGRATION
# ============================================================================

# Get MPI flags for pkg-config
if(TARGET MPI::MPI_CXX)
    get_target_property(MPI_CXX_COMPILE_OPTIONS MPI::MPI_CXX INTERFACE_COMPILE_OPTIONS)
    get_target_property(MPI_CXX_LINK_LIBRARIES MPI::MPI_CXX INTERFACE_LINK_LIBRARIES)


	message("-------->>>>>>>>>>>>>>>>MPI_CXX_COMPILE_OPTIONS is \"${MPI_CXX_COMPILE_OPTIONS}\"")


    if(MPI_CXX_COMPILE_OPTIONS)
        string(JOIN " " MPI_CXX_COMPILE_OPTIONS_STR ${MPI_CXX_COMPILE_OPTIONS})
    else()
        set(MPI_CXX_COMPILE_OPTIONS_STR "")
    endif()
    if(MPI_CXX_LINK_LIBRARIES)
        string(JOIN " " MPI_CXX_LIBRARIES_STR ${MPI_CXX_LINK_LIBRARIES})
    else()
        set(MPI_CXX_LIBRARIES_STR "")
    endif()
endif()

# ============================================================================
# COMMON INSTALLATION SETTINGS
# ============================================================================

# Add install RPATH
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)