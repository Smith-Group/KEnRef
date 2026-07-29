# ============================================================================
# KENREF GMX INTERFACE MODULE
# GPL-compatible - depends on GROMACS
# ============================================================================

message(STATUS "Configuring KENREF GMX Interface...")

# Check if KENREF CORE is available
if(NOT TARGET KENREF::CORE)
    message(FATAL_ERROR "KENREF GMX Interface requires KENREF CORE. Either build with BUILD_KENREF_CORE=ON or ensure KENREF CORE is installed.")
endif()

# CONFIGURE_DEPENDS: re-glob at build time so new files are auto-detected without a manual reconfigure.
FILE(GLOB gmxinterface_sources CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/gmxinterface/*.cpp")

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
    # GROMACS includes - THE ORDER IS CRITICAL.
    # BUILD_INTERFACE-only: these are absolute paths into THIS machine's gromacs source/build trees. They are
    # needed to COMPILE kenref_gmxinterface here, but must not leak into the INSTALLED target's interface —
    # they are meaningless on another machine, and when gromacs is auto-fetched INTO this build dir CMake
    # rightly refuses to export a build-dir path ("...is prefixed in the build directory", generate fails).
    # A downstream consumer supplies its own gromacs include paths.
    $<BUILD_INTERFACE:${GROMACS_BUILD_DIR}/src/include/>                    # 1. config.h
    $<BUILD_INTERFACE:${GROMACS_SRC_DIR}/src/include>                       # 2. Core API
    $<BUILD_INTERFACE:${GROMACS_SRC_DIR}/src/gromacs/math/include>          # 3. Math
    $<BUILD_INTERFACE:${GROMACS_SRC_DIR}/src/gromacs/mdtypes>               # 4. MD types
    $<BUILD_INTERFACE:${GROMACS_SRC_DIR}/src/gromacs/utility/include>       # 5. Utility
    $<BUILD_INTERFACE:${GROMACS_SRC_DIR}/src/programs>                      # 6. Programs
    $<BUILD_INTERFACE:${GROMACS_SRC_DIR}/src>                               # 7. I am not sure
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
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
)

# ============================================================================
# SIMD / EIGEN-ABI CONSUMER CHECK, FORCED AT INSTALL
# ============================================================================
# Compiled with the gmx build's flags and linked against KENREF::CORE (in-tree, or a pre-installed core for a
# gmx-only build). This is the check that actually bites: if the linked core was built with a different ACCEL
# than this gmx build, its baked alignment differs from this TU's -> non-zero -> install FATALs.
add_executable(kenref_gmx_simd_check "${CMAKE_CURRENT_SOURCE_DIR}/cmake/simd_check_main.cpp")
target_link_libraries(kenref_gmx_simd_check PRIVATE KENREF::CORE)
target_include_directories(kenref_gmx_simd_check PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include_core")
if(NOT CMAKE_CROSSCOMPILING)
    # rc 1 = genuine mismatch; any other non-zero = the check could not RUN (see KEnRefCore.cmake) -> warn only.
    install(CODE "
        message(STATUS \"KEnRef: verifying kenref-gmx vs kenref_core SIMD/Eigen ABI ...\")
        execute_process(COMMAND \"$<TARGET_FILE:kenref_gmx_simd_check>\" RESULT_VARIABLE _kn_gmx_simd_rc)
        if(_kn_gmx_simd_rc EQUAL 1)
            message(FATAL_ERROR \"kenref-gmx / kenref_core SIMD/Eigen ABI MISMATCH — rebuild both with the same ACCEL and Eigen.\")
        elseif(NOT _kn_gmx_simd_rc EQUAL 0)
            message(WARNING \"KEnRef: could not RUN the kenref-gmx SIMD/Eigen ABI check (rc=\${_kn_gmx_simd_rc}); skipping it \"
                            \"(not a mismatch — commonly the compiler's runtime libs are not loadable here).\")
        endif()
    ")
endif()

# ============================================================================
# TESTING
# ============================================================================
add_subdirectory(google_tests)
