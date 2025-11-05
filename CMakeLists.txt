cmake_minimum_required(VERSION 3.24)
project(KEnRef)
set(CMAKE_CXX_STANDARD 17)

string(TOLOWER ${CMAKE_BUILD_TYPE} build_type)
message(STATUS "the build type is: ${build_type}")

### comment the next line to make it float. Uncomment to make it double.
#add_compile_options(-DDOUBLE)

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
## comment the next line to make it float. Uncomment to make it double.
add_compile_options(-DDOUBLE)

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


FILE(GLOB core_sources "src/core/*.cpp")
FILE(GLOB gmxinterface_sources "src/gmxinterface/*.cpp")

include_directories(include)
include_directories(
#		"${GROMACS_INSTALL_DIR}/include"  # Not neeed here anymore, but needed for libgromacs
        "${GROMACS_BUILD_DIR}/src/include/"  #  Needed for "config.h"
#        "${GROMACS_SRC_DIR}/src/gromacs/gmxlib"
        "${GROMACS_SRC_DIR}/src/gromacs/mdtypes"
        "${GROMACS_SRC_DIR}/src/gromacs/utility/include"
        "${GROMACS_SRC_DIR}/src/programs"
        "${GROMACS_SRC_DIR}/src"
        "${GROMACS_SRC_DIR}/src/include"
        "${GROMACS_SRC_DIR}/src/gromacs/math/include"
)
#include_directories("${EIGEN3_INCLUDE_DIR}")

# Create libraries with simple names
add_library(kenref_core STATIC ${core_sources})           # Lowercase, no namespace
add_library(kenref_gmxinterface STATIC ${gmxinterface_sources})
# Create namespace aliases
add_library(KENREF::CORE ALIAS kenref_core)
add_library(KENREF::GMXINTERFACE ALIAS kenref_gmxinterface)

add_executable(KEnRef "src/gmxinterface/gmxwrapper.cpp")
add_executable(s2calc "src/tools/S2OrderParamsCalculator.cpp" "src/tools/AbstractCalculator.cpp")
add_executable(energycalc "src/tools/EnergyCalculator.cpp" "src/tools/AbstractCalculator.cpp")

target_link_libraries(kenref_core PUBLIC MPI::MPI_CXX Eigen3::Eigen
		#        Eigen3::unsupported # may need to restore it later
)
target_link_libraries(kenref_gmxinterface PUBLIC KENREF::CORE Gromacs::libgromacs)
target_link_libraries(KEnRef PRIVATE KENREF::CORE KENREF::GMXINTERFACE)
target_link_libraries(s2calc PRIVATE KENREF::CORE KENREF::GMXINTERFACE)
target_link_libraries(energycalc PRIVATE KENREF::CORE KENREF::GMXINTERFACE)

install(TARGETS KEnRef DESTINATION bin)
install(TARGETS s2calc DESTINATION bin)
install(TARGETS energycalc DESTINATION bin)
install(TARGETS kenref_core kenref_gmxinterface DESTINATION lib)
# Optional: Export with namespaces
install(TARGETS kenref_core kenref_gmxinterface
		EXPORT KEnRefTargets
		LIBRARY DESTINATION lib
		ARCHIVE DESTINATION lib
)
install(EXPORT KEnRefTargets
		FILE KEnRefTargets.cmake
		NAMESPACE KENREF::
		DESTINATION lib/cmake/KEnRef
)
###########################################
## Add gtest
add_subdirectory(google_tests)
###########################################

