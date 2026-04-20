# ============================================================================
# KENREF CORE MODULE
# Permissive license - no GROMACS dependencies
# ============================================================================

message(STATUS "Configuring KENREF CORE...")

# SEPARATE CORE SOURCES (Permissive license)
FILE(GLOB core_sources "${CMAKE_CURRENT_SOURCE_DIR}/src/core/*.cpp")  # Only pure math/algorithms

# Create separate targets with clear boundaries (libraries with simple names)
add_library(kenref_core STATIC ${core_sources})           # Lowercase, no namespace

# Create namespace aliases
add_library(KENREF::CORE ALIAS kenref_core)

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

# Separate include directories
target_include_directories(kenref_core PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include_core>
    $<INSTALL_INTERFACE:include/core>
    ${EIGEN3_INCLUDE_DIR} # TODO check its license
)

# ============================================================================
# PKG-CONFIG GENERATION FOR PLUMED INTEGRATION
# ============================================================================

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

# ============================================================================
# INSTALLATION
# ============================================================================

# Install libraries
install(TARGETS kenref_core kenref_and_eigen3
		EXPORT KEnRefCoreTargets
		DESTINATION lib
		INCLUDES DESTINATION include
)

# Install headers
install(DIRECTORY include_core/
		DESTINATION include/core
		FILES_MATCHING PATTERN "*.h"
)

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
		DESTINATION lib/pkgconfig
)

# Also install to share/pkgconfig for compatibility
install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/kenref_core.pc"
    "${CMAKE_CURRENT_BINARY_DIR}/kenref_and_eigen3.pc"
    DESTINATION share/pkgconfig
)

# Install export configuration
install(EXPORT KEnRefCoreTargets
		FILE KEnRefCoreConfig.cmake
		NAMESPACE KENREF::
		DESTINATION lib/cmake/KEnRef
)

# Install license files
install(FILES LICENSE_CORE.txt DESTINATION share/licenses/kenref)
install(FILES LICENSE_CLI11.txt DESTINATION share/licenses/cli11)

# ============================================================================
# ADDITIONAL TARGET PROPERTIES
# ============================================================================

# Set VERSION and SOVERSION properties for shared libraries (if any)
set_target_properties(kenref_core kenref_and_eigen3
    PROPERTIES
    VERSION 1.0.0
    SOVERSION 1
)

# ============================================================================
# TESTING
# ============================================================================

# CREATE A SIMPLE TEST FOR PLUMED INTEGRATION
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
