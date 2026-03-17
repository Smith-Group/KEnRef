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

# SEPARATE CORE SOURCES (Permissive license)
FILE(GLOB core_sources "src/core/*.cpp")  # Only pure math/algorithms
FILE(GLOB gmxinterface_sources "src/gmxinterface/*.cpp")
#FILE(GLOB tools_sources "src/tools/*.cpp")  # Tools depend on both
#FILE(GLOB core_headers "src/core/*.h")
#FILE(GLOB gmxinterface_headers "src/gmxinterface/*.h")

# Create separate targets with clear boundaries (libraries with simple names)
add_library(kenref_core STATIC ${core_sources})           # Lowercase, no namespace
add_library(kenref_gmxinterface STATIC ${gmxinterface_sources})

# Create namespace aliases
add_library(KENREF::CORE ALIAS kenref_core)
add_library(KENREF::GMXINTERFACE ALIAS kenref_gmxinterface)

# ============================================================================
# CREATE KENREF_AND_EIGEN3 LIBRARY (FOR PLUMED WHEN EIGEN3 IS NOT AVAILABLE)
# ============================================================================

message(STATUS "Preparing kenref_and_eigen3 library...")

# Create a target that includes Eigen headers
add_library(kenref_and_eigen3 STATIC ${core_sources})

# Copy Eigen headers to our source tree for distribution
set(EIGEN_HEADERS_DEST "${CMAKE_CURRENT_BINARY_DIR}/eigen_headers")
file(MAKE_DIRECTORY "${EIGEN_HEADERS_DEST}")

# Function to copy Eigen headers
function(copy_eigen_headers eigen_source_dir dest_dir)
    file(GLOB_RECURSE eigen_headers
        "${eigen_source_dir}/Eigen/*"
        "${eigen_source_dir}/unsupported/Eigen/*"
    )
    foreach(header IN LISTS eigen_headers)
        file(RELATIVE_PATH rel_path "${eigen_source_dir}" "${header}")
        set(dest_file "${dest_dir}/${rel_path}")
        get_filename_component(dest_dir_path "${dest_file}" DIRECTORY)
        file(MAKE_DIRECTORY "${dest_dir_path}")
        configure_file("${header}" "${dest_file}" COPYONLY)
    endforeach()
endfunction()

# Copy Eigen headers if they exist
if(EXISTS "${EIGEN3_INCLUDE_DIR}/Eigen")
    copy_eigen_headers("${EIGEN3_INCLUDE_DIR}" "${EIGEN_HEADERS_DEST}")
    message(STATUS "Eigen headers copied to: ${EIGEN_HEADERS_DEST}")
else()
    message(WARNING "Eigen headers not found at ${EIGEN3_INCLUDE_DIR}/Eigen")
endif()

# Create the kenref_and_eigen3 library
target_include_directories(kenref_and_eigen3 PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include_core>
    $<BUILD_INTERFACE:${EIGEN_HEADERS_DEST}>
    $<INSTALL_INTERFACE:include/core>
    $<INSTALL_INTERFACE:include/eigen>
)

target_link_libraries(kenref_and_eigen3 PUBLIC
    MPI::MPI_CXX
)

# Create namespace alias
add_library(KENREF::CORE_WITH_EIGEN ALIAS kenref_and_eigen3)

# ============================================================================
# SEPARATE DEPENDENCIES based on license compatibility
# ============================================================================

# Core: Only permissively licensed dependencies
target_link_libraries(kenref_core PUBLIC
    MPI::MPI_CXX           # Check MPI license
    Eigen3::Eigen          # MPL2 - permissive
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
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include_core>
    $<INSTALL_INTERFACE:include/core>
    ${EIGEN3_INCLUDE_DIR} # TODO check its license
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
    ${EIGEN3_INCLUDE_DIR}
)

# ============================================================================
# EXECUTABLES
# ============================================================================

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

# Create pkg-config file for kenref_core (with external Eigen)
configure_file(
		"${CMAKE_CURRENT_SOURCE_DIR}/cmake/kenref_core.pc.in"
		"${CMAKE_CURRENT_BINARY_DIR}/kenref_core.pc"
		@ONLY
)

# Create pkg-config file for kenref_and_eigen3 (self-contained)
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/kenref_and_eigen3.pc.in"
    "${CMAKE_CURRENT_BINARY_DIR}/kenref_and_eigen3.pc"
    @ONLY
)

# Create pkg-config file for kenref_gmxinterface (with GROMACS dependency)
configure_file(
		"${CMAKE_CURRENT_SOURCE_DIR}/cmake/kenref_gmxinterface.pc.in"
		"${CMAKE_CURRENT_BINARY_DIR}/kenref_gmxinterface.pc"
		@ONLY
)

# ============================================================================
# INSTALLATION - Install in BOTH locations
# ============================================================================

# Install executables
install(TARGETS KEnRef s2calc energycalc
		RUNTIME DESTINATION bin
		LIBRARY DESTINATION lib
		ARCHIVE DESTINATION lib
)

# Install libraries
install(TARGETS kenref_core kenref_gmxinterface kenref_and_eigen3
		EXPORT KEnRefTargets
		DESTINATION lib
		INCLUDES DESTINATION include
)

# Install headers
install(DIRECTORY include_core/
		DESTINATION include/core
		FILES_MATCHING PATTERN "*.h"
)

install(DIRECTORY include_gmxinterface/
		DESTINATION include/gmxinterface
		FILES_MATCHING PATTERN "*.h"
)

#------------------------
# Install pkg-config files
#------------------------

# For self-contained version, install Eigen headers
install(DIRECTORY "${EIGEN_HEADERS_DEST}/Eigen"
    DESTINATION include/eigen
    PATTERN "*.h"
)

if(EXISTS "${EIGEN_HEADERS_DEST}/unsupported")
    install(DIRECTORY "${EIGEN_HEADERS_DEST}/unsupported/Eigen"
        DESTINATION include/unsupported/eigen
        PATTERN "*.h"
    )
endif()

# Install pkg-config files in BOTH standard locations
install(FILES
		"${CMAKE_CURRENT_BINARY_DIR}/kenref_core.pc"
    "${CMAKE_CURRENT_BINARY_DIR}/kenref_and_eigen3.pc"
		"${CMAKE_CURRENT_BINARY_DIR}/kenref_gmxinterface.pc"
		DESTINATION lib/pkgconfig
)

# Also install to share/pkgconfig for compatibility
install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/kenref_core.pc"
    "${CMAKE_CURRENT_BINARY_DIR}/kenref_and_eigen3.pc"
    "${CMAKE_CURRENT_BINARY_DIR}/kenref_gmxinterface.pc"
    DESTINATION share/pkgconfig
)

# Install export configuration
install(EXPORT KEnRefTargets
		FILE KEnRefConfig.cmake
		NAMESPACE KENREF::
		DESTINATION lib/cmake/KEnRef
)

# Install license files
install(FILES LICENSE_CORE.txt DESTINATION share/licenses/kenref)
install(FILES LICENSE_GMX.txt DESTINATION share/licenses/kenref)
install(FILES LICENSE_CLI11.txt DESTINATION share/licenses/cli11)


# ============================================================================
# TESTING
# ============================================================================
add_subdirectory(google_tests)

# ============================================================================
# CREATE A SIMPLE TEST FOR PLUMED INTEGRATION
# ============================================================================
add_executable(test_kenref_plumed
		${CMAKE_CURRENT_SOURCE_DIR}/cmake/test_plumed_integration.cpp
)
target_link_libraries(test_kenref_plumed PRIVATE kenref_and_eigen3)
target_include_directories(test_kenref_plumed PRIVATE
		${CMAKE_CURRENT_SOURCE_DIR}/include_core
)

# Test that both libraries work
add_executable(test_kenref_core
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/test_plumed_integration.cpp
)
target_link_libraries(test_kenref_core PRIVATE kenref_core)
target_include_directories(test_kenref_core PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include_core
		${EIGEN3_INCLUDE_DIR}
)

# ============================================================================
# ADDITIONAL TARGET PROPERTIES
# ============================================================================

# Set VERSION and SOVERSION properties for shared libraries (if any)
set_target_properties(kenref_core kenref_gmxinterface kenref_and_eigen3
    PROPERTIES
    VERSION 1.0.0
    SOVERSION 1
)

# Add install RPATH
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)