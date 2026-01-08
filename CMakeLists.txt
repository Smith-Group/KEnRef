cmake_minimum_required(VERSION 3.24)
project(KEnRef)
set(CMAKE_CXX_STANDARD 17)

string(TOLOWER ${CMAKE_BUILD_TYPE} build_type)
message(STATUS "the build type is: ${build_type}")

### comment the next line to make it float. Uncomment to make it double.
add_compile_options(-DKENREF_DOUBLE)

if (NOT DEFINED ACCEL)
	set(ACCEL "AVX2_256")
	message(STATUS "ACCEL variable not defined. Setting it to ${ACCEL}")
else ()
	message(STATUS "ACCEL variable is: ${ACCEL}")
endif ()
#TODO rearrange setting *_DIR variables after reading https://cmake.org/cmake/help/latest/command/find_package.html [then delete this line :)]
##change these if needed
set(SUFFIX_MPI _mpi)
set(GROMACS_SRC_DIR /home/amr/git/gromacs) # make sure you `git checkout` the right Gromacs branch and tag first.
set(GROMACS_BUILD_DIR /home/amr/git/gromacs/cmake-build-${build_type}) # build_type is all_small
set(GROMACS_INSTALL_DIR /smithlab/opt/gromacs-dev/2025/2025.3/${build_type}/${ACCEL})
set(GROMACS_DIR ${GROMACS_INSTALL_DIR}/share/cmake/gromacs${SUFFIX_MPI}) #where libgromacs*.cmake are located

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
	string(FIND ${GROMACS_INSTALL_DIR} "gromacs-dev" replaceable)
	if (replaceable GREATER_EQUAL -1)
		string(REPLACE "gromacs-dev" "kenref-dev" proposed_path ${GROMACS_INSTALL_DIR})
		set(CMAKE_INSTALL_PREFIX ${proposed_path} CACHE PATH "root directory to install KEnRef (default is a tree parallel to GROMACS_INSTALL_DIR)" FORCE)
	endif ()
endif()

#set(CMAKE_MODULE_PATH /usr/local/gromacs-dev/share/gromacs/template/cmake:$CMAKE_MODULE_PATH)
set(Eigen3_DIR /smithlab/opt/eigen/3.4.0/share/eigen3/cmake) #set(Eigen3_DIR /smithlab/opt/eigen/3.4.9/share/eigen3/cmake)

find_package(MPI REQUIRED)

find_package (Eigen3 3.4 REQUIRED)
if (Eigen3_FOUND)
	message("Eigen3 FOUND")
endif ()

message(STATUS "Searching for Gromacs in ${GROMACS_DIR}")
find_package(GROMACS NAMES gromacs_mpi gromacs REQUIRED)
if (GROMACS_FOUND)
	message("GROMACS FOUND")
	message(STATUS "GROMACS_DIR folder is: ${GROMACS_DIR}")
	message(STATUS "GROMACS_CONFIG file is: ${GROMACS_CONFIG}")
endif ()

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
		CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${Eigen3_INCLUDE_DIRS}"
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

# SEPARATE CORE SOURCES (Permissive license)
FILE(GLOB core_sources "src/core/*.cpp")  # Only pure math/algorithms
FILE(GLOB gmxinterface_sources "src/gmxinterface/*.cpp")
FILE(GLOB tools_sources "src/tools/*.cpp")  # Tools depend on both

# Create separate targets with clear boundaries (libraries with simple names)
add_library(kenref_core STATIC ${core_sources})           # Lowercase, no namespace
add_library(kenref_gmxinterface STATIC ${gmxinterface_sources})

# Create namespace aliases
add_library(KENREF::CORE ALIAS kenref_core)
add_library(KENREF::GMXINTERFACE ALIAS kenref_gmxinterface)

# SEPARATE DEPENDENCIES based on license compatibility
# Core: Only permissively licensed dependencies
target_link_libraries(kenref_core PUBLIC
    MPI::MPI_CXX           # Check MPI license
    Eigen3::Eigen          # MPL2 - permissive
#        Eigen3::unsupported # may need to restore it later
    # NO GROMACS here!
)

# GMXInterface: GPL-compatible dependencies
target_link_libraries(kenref_gmxinterface PUBLIC
    KENREF::CORE           # Core is permissively licensed
    Gromacs::libgromacs    # GPL - only here!
    MPI::MPI_CXX           # Check MPI license
)

# Separate include directories
target_include_directories(kenref_core PUBLIC
    #${EIGEN3_INCLUDE_DIR}
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include_core>
    $<INSTALL_INTERFACE:include/core>
)
target_include_directories(kenref_gmxinterface PUBLIC
    # GROMACS includes - THE ORDER IS CRITICAL
    "${GROMACS_BUILD_DIR}/src/include/"                    # 1. config.h
    "${GROMACS_SRC_DIR}/src/include"                       # 2. Core API
    "${GROMACS_SRC_DIR}/src/gromacs/math/include"          # 3. Math
    "${GROMACS_SRC_DIR}/src/gromacs/mdtypes"               # 4. MD types
    "${GROMACS_SRC_DIR}/src/gromacs/utility/include"       # 5. Utility
    "${GROMACS_SRC_DIR}/src/programs"                      # 6. Programs
    "${GROMACS_SRC_DIR}/src"							   # 7. I am not sure
    # Our interface (lowest priority)
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include_gmxinterface>
    $<INSTALL_INTERFACE:include/gmxinterface>
)

add_executable(KEnRef "src/gmxinterface/gmxwrapper.cpp")
target_link_libraries(KEnRef PRIVATE KENREF::CORE KENREF::GMXINTERFACE)

# Tools depend on both
add_executable(s2calc "src/tools/S2OrderParamsCalculator.cpp" "src/tools/AbstractCalculator.cpp")
add_executable(energycalc "src/tools/EnergyCalculator.cpp" "src/tools/AbstractCalculator.cpp")
target_link_libraries(s2calc PRIVATE
    KENREF::CORE           # Permissive part
    KENREF::GMXINTERFACE   # GPL part - tools become GPL
)

target_link_libraries(energycalc PRIVATE
    KENREF::CORE
    KENREF::GMXINTERFACE   # Tools inherit GPL
)

install(TARGETS KEnRef DESTINATION bin)
install(TARGETS s2calc DESTINATION bin)
install(TARGETS energycalc DESTINATION bin)

# SEPARATE INSTALLATION with license files
install(TARGETS kenref_core
    EXPORT KEnRefCoreTargets
    DESTINATION lib
    # Include permissive LICENSE file
)
install(TARGETS kenref_gmxinterface
    EXPORT KEnRefGmxTargets
    DESTINATION lib
    # Include GPL-compatible LICENSE file
)

# Separate exports with different namespaces or configs
install(EXPORT KEnRefCoreTargets
    FILE KEnRefCoreConfig.cmake
    NAMESPACE KENREF::
    DESTINATION lib/cmake/KEnRef
)

install(EXPORT KEnRefGmxTargets
    FILE KEnRefGmxConfig.cmake
    NAMESPACE KENREF::
    DESTINATION lib/cmake/KEnRef
)

install(FILES LICENSE_CORE.txt DESTINATION share/licenses/kenref)
install(FILES LICENSE_GMX.txt DESTINATION share/licenses/kenref)
install(FILES LICENSE_CLI11.txt DESTINATION share/licenses/cli11)
###########################################
## Add gtest
add_subdirectory(google_tests)
###########################################

