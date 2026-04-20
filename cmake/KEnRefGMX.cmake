# ============================================================================
# KENREF GMX INTERFACE MODULE
# GPL-compatible - depends on GROMACS
# ============================================================================

message(STATUS "Configuring KENREF GMX Interface...")

# Check if KENREF CORE is available
if(NOT TARGET KENREF::CORE)
    message(FATAL_ERROR "KENREF GMX Interface requires KENREF CORE. Either build with BUILD_KENREF_CORE=ON or ensure KENREF CORE is installed.")
endif()

FILE(GLOB gmxinterface_sources "${CMAKE_CURRENT_SOURCE_DIR}/src/gmxinterface/*.cpp")

add_library(kenref_gmxinterface STATIC ${gmxinterface_sources})

# Create namespace aliases
add_library(KENREF::GMXINTERFACE ALIAS kenref_gmxinterface)

# ============================================================================
# SEPARATE DEPENDENCIES based on license compatibility
# ============================================================================

# GMXInterface: GPL-compatible dependencies
target_link_libraries(kenref_gmxinterface PUBLIC
    KENREF::CORE           # Core is permissively licensed
    Gromacs::libgromacs    # GPL - only here!
    MPI::MPI_CXX           # Check MPI license
)

# Separate include directories
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
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include_core>  # For CLI11 and KENREF headers
    $<INSTALL_INTERFACE:include/gmxinterface>
    $<INSTALL_INTERFACE:include/core>  # For CLI11 in installed version
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
# PKG-CONFIG GENERATION
# ============================================================================

# Create pkg-config file for kenref_gmxinterface (with GROMACS dependency)
configure_file(
		"${CMAKE_CURRENT_SOURCE_DIR}/cmake/kenref_gmxinterface.pc.in"
		"${CMAKE_CURRENT_BINARY_DIR}/kenref_gmxinterface.pc"
		@ONLY
)

# ============================================================================
# INSTALLATION
# ============================================================================

# Install executables
install(TARGETS KEnRef s2calc energycalc
		RUNTIME DESTINATION bin
		LIBRARY DESTINATION lib
		ARCHIVE DESTINATION lib
)

# Install libraries
install(TARGETS kenref_gmxinterface
		EXPORT KEnRefGMXTargets
		DESTINATION lib
		INCLUDES DESTINATION include
)

# Install headers
install(DIRECTORY include_gmxinterface/
		DESTINATION include/gmxinterface
		FILES_MATCHING PATTERN "*.h"
)

# Install pkg-config files
install(FILES
		"${CMAKE_CURRENT_BINARY_DIR}/kenref_gmxinterface.pc"
		DESTINATION lib/pkgconfig
)

# Also install to share/pkgconfig for compatibility
install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/kenref_gmxinterface.pc"
    DESTINATION share/pkgconfig
)

# Install export configuration
install(EXPORT KEnRefGMXTargets
		FILE KEnRefGMXConfig.cmake
		NAMESPACE KENREF::
		DESTINATION lib/cmake/KEnRef
)

# Install license files
install(FILES LICENSE_GMX.txt DESTINATION share/licenses/kenref)

# ============================================================================
# ADDITIONAL TARGET PROPERTIES
# ============================================================================

# Set VERSION and SOVERSION properties for shared libraries (if any)
set_target_properties(kenref_gmxinterface
    PROPERTIES
    VERSION 1.0.0
    SOVERSION 1
)

# ============================================================================
# TESTING
# ============================================================================
add_subdirectory(google_tests)
