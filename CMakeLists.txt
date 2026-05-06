cmake_minimum_required(VERSION 3.24)
project(KEnRef)
set(CMAKE_CXX_STANDARD 17)

string(TOLOWER ${CMAKE_BUILD_TYPE} build_type)
message(STATUS "the build type is: ${build_type}")

# ============================================================================
# BUILD OPTIONS
# ============================================================================
option(BUILD_KENREF_CORE "Build KENREF_CORE library" ON)
option(BUILD_KENREF_GMX "Build KENREF GMX interface" ON)

### comment the next line to make it float. Uncomment to make it double.
add_compile_options(-DKENREF_DOUBLE)

if (NOT DEFINED ACCEL)
	set(ACCEL "AVX2_256")
	message(STATUS "ACCEL variable not defined. Setting it to ${ACCEL}")
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
			set(GROMACS_INSTALL_DIR "/smithlab/opt/gromacs-4-kenref/2025/2025.4/${build_type}/${ACCEL}")
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
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
	if(BUILD_KENREF_CORE AND NOT BUILD_KENREF_GMX)
		# Core only: /smithlab/opt/kenref/${build_type}/${ACCEL}
		set(CMAKE_INSTALL_PREFIX "/smithlab/opt/kenref/${build_type}/${ACCEL}" CACHE PATH "root directory to install KEnRef CORE" FORCE)
	elseif(BUILD_KENREF_GMX)
		# GMX interface (with or without core): /smithlab/opt/kenref-gmx/${year}/${GROMACS_VERSION}/${build_type}/${ACCEL}
		# Extract GROMACS year and version from GROMACS_INSTALL_DIR
		# Expected format: /smithlab/opt/gromacs-*/2025/2025.x/${build_type}/${ACCEL}
		string(REGEX MATCH "gromacs-4-kenref/([0-9]+)/([0-9.]+)" GROMACS_VERSION_MATCH "${GROMACS_INSTALL_DIR}")
		if(CMAKE_MATCH_1 AND CMAKE_MATCH_2)
			set(GROMACS_YEAR "${CMAKE_MATCH_1}")
			set(GROMACS_VERSION "${CMAKE_MATCH_2}")
		else()
			# Fallback: try to get from GROMACS package if available
			if(GROMACS_VERSION)
				# Use existing GROMACS_VERSION variable
				string(REGEX MATCH "^([0-9]+)" YEAR_MATCH "${GROMACS_VERSION}")
				if(CMAKE_MATCH_1)
					set(GROMACS_YEAR "${CMAKE_MATCH_1}")
				else()
					set(GROMACS_YEAR "unknown")
				endif()
			else()
				message(WARNING "Could not extract GROMACS version from GROMACS_INSTALL_DIR: ${GROMACS_INSTALL_DIR}")
				set(GROMACS_YEAR "unknown")
				set(GROMACS_VERSION "unknown")
			endif()
		endif()
		set(CMAKE_INSTALL_PREFIX "/smithlab/opt/kenref-gmx/${GROMACS_YEAR}/${GROMACS_VERSION}/${build_type}/${ACCEL}" CACHE PATH "root directory to install KEnRef GMX" FORCE)
	endif()
endif()

#set(CMAKE_MODULE_PATH /usr/local/gromacs-dev/share/gromacs/template/cmake:$CMAKE_MODULE_PATH)
set(Eigen3_DIR /smithlab/opt/eigen/5.0.1/${build_type}/share/eigen3/cmake)

# ============================================================================
# COMMON DEPENDENCIES
# ============================================================================
find_package(MPI REQUIRED)

if(BUILD_KENREF_CORE)
	find_package (Eigen3 5.0 REQUIRED)
endif()
if (Eigen3_FOUND)
	message("Eigen3 FOUND")
	# Get Eigen include directory for pkg-config
	get_target_property(EIGEN3_INCLUDE_DIR Eigen3::Eigen INTERFACE_INCLUDE_DIRECTORIES)
	if(NOT EIGEN3_INCLUDE_DIR)
		set(EIGEN3_INCLUDE_DIR ${Eigen3_INCLUDE_DIRS})
	endif()
	message(STATUS "Eigen3 include directory: ${EIGEN3_INCLUDE_DIR}")
	# Extract just the directory path without -I flag for pkg-config
	string(REGEX REPLACE "^-I" "" EIGEN3_INCLUDE_DIR_CLEAN "${EIGEN3_INCLUDE_DIR}")
	set(EIGEN3_INCLUDE_DIR "${EIGEN3_INCLUDE_DIR_CLEAN}")
	message(STATUS "Eigen3 include directory (clean): ${EIGEN3_INCLUDE_DIR}")
endif ()

if(BUILD_KENREF_GMX)
	message(STATUS "Searching for Gromacs in ${GROMACS_DIR}")
	find_package(GROMACS NAMES gromacs_mpi gromacs REQUIRED)
	if (GROMACS_FOUND)
		message("GROMACS FOUND")
		message(STATUS "GROMACS_DIR folder is: ${GROMACS_DIR}")
		message(STATUS "GROMACS_CONFIG file is: ${GROMACS_CONFIG}")
	endif ()
endif()

if(BUILD_KENREF_CORE)
	find_package(OpenMP REQUIRED)
	# ============================================================================
	# TEST EIGEN MATRIX REDUCTION SUPPORT
	# ============================================================================
	message(STATUS "Testing Eigen matrix reduction support with OpenMP...")

# Try to compile and run the test
try_run(
		EIGEN_REDUCTION_TEST_EXITCODE
		EIGEN_REDUCTION_TEST_COMPILED
		${CMAKE_BINARY_DIR}
		${CMAKE_CURRENT_SOURCE_DIR}/cmake/test_eigen_reduction.cpp
		COMPILE_DEFINITIONS ${OpenMP_CXX_FLAGS}
		LINK_LIBRARIES ${OpenMP_CXX_LIBRARIES} Eigen3::Eigen
		CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${EIGEN3_INCLUDE_DIR}"
		OUTPUT_VARIABLE EIGEN_TEST_OUTPUT
)

# Check if test compiled and ran successfully
if(EIGEN_REDUCTION_TEST_COMPILED AND EIGEN_REDUCTION_TEST_EXITCODE EQUAL 0)
	set(EIGEN_MATRIX_REDUCTION_SUPPORTED TRUE)
	message(STATUS "✓ Eigen matrix reduction supported: YES")
else()
	set(EIGEN_MATRIX_REDUCTION_SUPPORTED FALSE)
	message(STATUS "✗ Eigen matrix reduction supported: NO")
	if(NOT EIGEN_REDUCTION_TEST_COMPILED)
		message(STATUS "  Test failed to compile")
	else()
		message(STATUS "  Test compiled but failed to run correctly")
	endif()
	message(STATUS "  Test output: ${EIGEN_TEST_OUTPUT}")
endif()

	# Set the preprocessor definition
	if(EIGEN_MATRIX_REDUCTION_SUPPORTED)
		add_compile_definitions(EIGEN_MATRIX_REDUCTION_SUPPORTED=1)
	else()
		add_compile_definitions(EIGEN_MATRIX_REDUCTION_SUPPORTED=0)
	endif()
	# ============================================================================
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
		set(KENREF_CORE_CMAKE_PATH "/smithlab/opt/kenref/${build_type}/${ACCEL}/lib/cmake/KEnRef")
		message(STATUS "  Looking in: ${KENREF_CORE_CMAKE_PATH}")

		# KENREF CORE depends on Eigen3, so find it first
		if(NOT TARGET Eigen3::Eigen)
			find_package(Eigen3 3.4 REQUIRED)
			if(Eigen3_FOUND)
				message(STATUS "  Found Eigen3 for KENREF CORE")
			endif()
		endif()

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