# ============================================================================
# KENREF CORE MODULE
# Permissive license - no GROMACS dependencies
# ============================================================================

message(STATUS "Configuring KENREF CORE...")

# SEPARATE CORE SOURCES (Permissive license)
# CONFIGURE_DEPENDS: re-glob at build time so newly added files are auto-detected (and registered with
# the targets / IDE) without a manual cmake reconfigure. Plain GLOB is evaluated only at configure time.
FILE(GLOB core_sources CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/core/*.cpp")  # Only pure math/algorithms
# Public headers — listed on the targets purely so IDEs (CLion/clangd) treat them as project files
# with the target's include flags (CMake does not compile .h). Keeps new headers recognised.
FILE(GLOB core_headers CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/include_core/core/*.h")

# Create separate targets with clear boundaries (libraries with simple names)
add_library(kenref_core STATIC ${core_sources} ${core_headers})  # Lowercase, no namespace

# Create namespace aliases
add_library(KENREF::CORE ALIAS kenref_core)

# ============================================================================
# CREATE KENREF_AND_EIGEN3 LIBRARY (FOR PLUMED WHEN EIGEN3 IS NOT AVAILABLE)
# ============================================================================

message(STATUS "Preparing kenref_and_eigen3 library...")

# Create a target that includes Eigen headers
add_library(kenref_and_eigen3 STATIC ${core_sources} ${core_headers})

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

# NOTE: no MPI here. kenref_core contains no MPI code; MPI belongs to the gmx interface + tools (which link it
# themselves), and the PLUMED bias uses PLUMED's own Communicator. See CMakeLists.txt (kenref_provide_mpi).

# Create namespace alias
add_library(KENREF::CORE_WITH_EIGEN ALIAS kenref_and_eigen3)

# ============================================================================
# SEPARATE DEPENDENCIES based on license compatibility
# ============================================================================

# Core: Only permissively licensed dependencies
target_link_libraries(kenref_core PUBLIC
    Eigen3::Eigen          # MPL2 - permissive
    # NO GROMACS here!  and NO MPI: the core has no MPI code (consumers that need it link it themselves).
)

# OpenMP (optional): when found, enables -fopenmp so the core's `#pragma omp` directives are active.
# Linked conditionally so a machine without OpenMP can still configure and build (serial) the core.
if(TARGET OpenMP::OpenMP_CXX)
    target_link_libraries(kenref_core PUBLIC OpenMP::OpenMP_CXX)
    target_link_libraries(kenref_and_eigen3 PUBLIC OpenMP::OpenMP_CXX)
endif()

# Separate include directories
target_include_directories(kenref_core PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include_core>
    $<INSTALL_INTERFACE:include/core>
    # BUILD-only: keeps the (possibly in-build-tree, e.g. FetchContent) Eigen path out of the exported
    # interface. Installed consumers resolve Eigen from the bundled kenref_and_eigen3 headers instead.
    $<BUILD_INTERFACE:${EIGEN3_INCLUDE_DIR}>
)

# ============================================================================
# SHARED VARIANT OF THE CORE (Phase E)
# ============================================================================
# Build kenref_core also as a SHARED library so distributors / PLUMED can link it dynamically. This is
# ADDITIVE: the static kenref_core above is untouched, so every existing consumer keeps linking the .a.
# Same sources + usage requirements; OUTPUT_NAME kenref_core => libkenref_core.so beside libkenref_core.a.
add_library(kenref_core_shared SHARED ${core_sources} ${core_headers})
add_library(KENREF::CORE_SHARED ALIAS kenref_core_shared)
target_link_libraries(kenref_core_shared PUBLIC Eigen3::Eigen)   # no MPI (see kenref_core above)
if(TARGET OpenMP::OpenMP_CXX)
    target_link_libraries(kenref_core_shared PUBLIC OpenMP::OpenMP_CXX)
endif()
target_include_directories(kenref_core_shared PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include_core>
    $<INSTALL_INTERFACE:include/core>
    $<BUILD_INTERFACE:${EIGEN3_INCLUDE_DIR}>   # build-only; keeps a fetched (in-tree) Eigen out of the export
)
set_target_properties(kenref_core_shared PROPERTIES
    OUTPUT_NAME kenref_core
    POSITION_INDEPENDENT_CODE ON
)

# ============================================================================
# PER-MODEL ENABLE DEFINITIONS (Style-2 registration; options in top CMakeLists.txt)
# ============================================================================
# bootstrapModels() guards each register*Model() with #if KENREF_ENABLE_<MODEL>; drive those from the
# CMake options so a model can be dropped at compile time. Applied to all core targets (PUBLIC so
# consumers that include the headers see the same switches).
foreach(kenref_lib kenref_core kenref_core_shared kenref_and_eigen3)
    target_compile_definitions(${kenref_lib} PUBLIC
        KENREF_ENABLE_SIGMA=$<BOOL:${KENREF_ENABLE_SIGMA}>
        KENREF_ENABLE_PLATEAUS=$<BOOL:${KENREF_ENABLE_PLATEAUS}>
        KENREF_ENABLE_RELAX=$<BOOL:${KENREF_ENABLE_RELAX}>
    )
endforeach()

# ============================================================================
# COMPILE-TIME KENREF <-> EIGEN ABI GUARD
# ============================================================================
# kenref_core stores Eigen objects in its own containers, so Eigen's alignment is baked into libkenref_core's
# ABI. Generate a header recording the EIGEN_MAX_ALIGN_BYTES / version kenref_core is built with (via a tiny
# generator compiled with THIS build's flags + Eigen), which core/EigenAbiCheck.h (included from core/KEnRef.h)
# static_asserts against — so a consumer TU whose Eigen disaligns from kenref_core FAILS TO COMPILE. Skipped
# when cross-compiling (the generator can't run on the build host); the runtime SIMD gate still covers it.
set(_kn_eigen_abi_hdr "${CMAKE_BINARY_DIR}/generated/core/KEnRefEigenAbi.h")
if(NOT CMAKE_CROSSCOMPILING)
    add_executable(kenref_eigen_abi_gen "${CMAKE_CURRENT_SOURCE_DIR}/cmake/eigen_abi_gen.cpp")
    target_link_libraries(kenref_eigen_abi_gen PRIVATE Eigen3::Eigen)
    # Run the generator through a wrapper that degrades gracefully: executing a just-built binary can fail in
    # perfectly valid environments (e.g. -stdlib=libc++ with the toolchain's runtime libs not on
    # LD_LIBRARY_PATH -> exec fails), and the compile-time guard is only a belt over the runtime check — it
    # must never break the build. See cmake/RunEigenAbiGen.cmake.
    add_custom_command(OUTPUT "${_kn_eigen_abi_hdr}"
        COMMAND "${CMAKE_COMMAND}" -E make_directory "${CMAKE_BINARY_DIR}/generated/core"
        COMMAND "${CMAKE_COMMAND}" "-DKN_GEN=$<TARGET_FILE:kenref_eigen_abi_gen>" "-DKN_OUT=${_kn_eigen_abi_hdr}"
                -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/RunEigenAbiGen.cmake"
        DEPENDS kenref_eigen_abi_gen
        COMMENT "Recording kenref_core's Eigen ABI (alignment/version) -> KEnRefEigenAbi.h"
        VERBATIM)
    add_custom_target(kenref_eigen_abi DEPENDS "${_kn_eigen_abi_hdr}")
    foreach(kenref_lib kenref_core kenref_core_shared kenref_and_eigen3)
        add_dependencies(${kenref_lib} kenref_eigen_abi)
        target_include_directories(${kenref_lib} PUBLIC $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/generated>)
    endforeach()
    # Install alongside the other core headers (include/core/core/*.h) so external consumers get the guard too.
    install(FILES "${_kn_eigen_abi_hdr}" DESTINATION include/core/core)
endif()

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

# Create pkg-config file for the PLUMED interface (Option A: exposes the repo's plumedinterface header
# + source include paths so the PLUMED fork's frozen frame compiles src/plumedinterface within PLUMED's
# build without hardcoding a checkout path). Installed/headers+source below.
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/kenref_plumed.pc.in"
    "${CMAKE_CURRENT_BINARY_DIR}/kenref_plumed.pc"
    @ONLY
)

# Self-contained counterpart, mirroring kenref_and_eigen3.pc on the LIBRARY side. Without it the two
# halves fall out of step: on a machine with no eigen3.pc the library check succeeds via
# kenref_and_eigen3 while the include check fails, because kenref_plumed.pc `Requires: eigen3`.
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/kenref_plumed_and_eigen3.pc.in"
    "${CMAKE_CURRENT_BINARY_DIR}/kenref_plumed_and_eigen3.pc"
    @ONLY
)

# ============================================================================
# INSTALLATION
# ============================================================================

# Install libraries
install(TARGETS kenref_core kenref_core_shared kenref_and_eigen3
		EXPORT KEnRefCoreTargets
		DESTINATION lib
		INCLUDES DESTINATION include
)

# Install headers
install(DIRECTORY include_core/
		DESTINATION include/core
		FILES_MATCHING PATTERN "*.h"
)

# Install the PLUMED interface headers + SOURCE — ONLY when asked (KENREF_EXPORT_PLUMEDINTERFACE). Option A:
# the fork compiles src/plumedinterface in its own build, so a PLUMED build needs these exported; a plain
# core build does not, and should not carry them. Headers -> include/plumedinterface/plumedinterface/*.h ;
# source -> src/plumedinterface/*. kenref_plumed.pc's Cflags (-I include/plumedinterface, -I src) make the
# frame's #include forwarders resolve.
if(KENREF_EXPORT_PLUMEDINTERFACE)
	install(DIRECTORY include_plumedinterface/
			DESTINATION include/plumedinterface
			FILES_MATCHING PATTERN "*.h"
	)
	install(DIRECTORY src/plumedinterface
			DESTINATION src
			FILES_MATCHING PATTERN "*.cpp"
	)
endif()

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

# Ship the KEnRef-BUILT Eigen (external flavor) INTO the kenref prefix, so the install is self-sufficient
# (external kenref_core Requires: eigen3 -> resolves from the same prefix). Only when WE built it
# (KENREF_EIGEN_FETCHED); an external/system Eigen is left where the user has it. Copies the staged Eigen
# tree (include/eigen3 + share/{pkgconfig,eigen3}) into the kenref prefix at `cmake --install`.
option(KENREF_INSTALL_EIGEN "Install the KEnRef-built Eigen into the kenref prefix (ships with the project)" ON)
if(KENREF_EIGEN_FETCHED AND KENREF_INSTALL_EIGEN)
    install(DIRECTORY "${KENREF_EIGEN_INSTALL_DIR}/" DESTINATION ".")
endif()

# Install pkg-config files in BOTH standard locations (core ones always; kenref_plumed.pc only when the
# plumedinterface is exported — it points at plumedinterface paths that aren't installed otherwise).
install(FILES
		"${CMAKE_CURRENT_BINARY_DIR}/kenref_core.pc"
    "${CMAKE_CURRENT_BINARY_DIR}/kenref_and_eigen3.pc"
		DESTINATION lib/pkgconfig
)
install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/kenref_core.pc"
    "${CMAKE_CURRENT_BINARY_DIR}/kenref_and_eigen3.pc"
    DESTINATION share/pkgconfig
)
if(KENREF_EXPORT_PLUMEDINTERFACE)
	install(FILES "${CMAKE_CURRENT_BINARY_DIR}/kenref_plumed.pc"
	              "${CMAKE_CURRENT_BINARY_DIR}/kenref_plumed_and_eigen3.pc" DESTINATION lib/pkgconfig)
	install(FILES "${CMAKE_CURRENT_BINARY_DIR}/kenref_plumed.pc"
	              "${CMAKE_CURRENT_BINARY_DIR}/kenref_plumed_and_eigen3.pc" DESTINATION share/pkgconfig)
endif()

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

# Set VERSION / SOVERSION from the single source of truth (project(KEnRef VERSION ...) in CMakeLists.txt),
# so the shared-library soname tracks the same version as the pkg-config files and the PLUMED floor.
#
# POSITION_INDEPENDENT_CODE applies to the STATIC archives too, not just kenref_core_shared. Our static
# libraries are not end products: consumers link them *into a shared object* -- PLUMED embeds one in
# libplumedKernel.so, and the GROMACS force provider does the same. A non-PIC .a cannot go into a .so at
# all; ld rejects it outright with
#     relocation R_X86_64_PC32 against symbol ... can not be used when making a shared object
# Setting it only on the shared target hid the bug locally, because where libkenref_core.so sits beside
# libkenref_core.a the linker silently prefers the .so. The self-contained kenref_and_eigen3 flavour has
# no shared counterpart, so on any machine without eigen3.pc -- which is what selects that flavour -- the
# non-PIC archive is the only candidate and the consumer's link fails.
set_target_properties(kenref_core kenref_core_shared kenref_and_eigen3
    PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
    POSITION_INDEPENDENT_CODE ON
)

# ----------------------------------------------------------------------------
# Version-consistency manifest — written by CMake at install time (single source of truth; replaces the
# hand-rolled one build.sh used to emit). Records what this kenref_core actually is, so a downstream
# PLUMED/GROMACS build (or a later run) can confirm it is pinned to the expected version/ISA.
# ----------------------------------------------------------------------------
execute_process(COMMAND git -C "${CMAKE_SOURCE_DIR}" rev-parse --short HEAD
    OUTPUT_VARIABLE _kn_git_sha OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
if(NOT _kn_git_sha)
    set(_kn_git_sha "n/a")
endif()
install(CODE "file(WRITE \"\${CMAKE_INSTALL_PREFIX}/kenref-build-manifest.txt\"
\"# KEnRef build manifest (generated by CMake at install)
kenref_version = ${PROJECT_VERSION}
kenref_git     = ${_kn_git_sha}
build_type     = ${CMAKE_BUILD_TYPE}
accel          = ${ACCEL}
cxx            = ${CMAKE_CXX_COMPILER}
\")")

# ============================================================================
# SIMD / EIGEN-ABI SELF-CHECK, FORCED AT INSTALL
# ============================================================================
# A tiny executable (core flags, links kenref_core) that returns non-zero if the caller's Eigen alignment
# differs from the linked core's. Run at `cmake --install` via install(CODE) so an ISA/ABI-mismatched build
# cannot be installed. Trivially passes for a single-configure core build; the real bite is the gmx/plumed
# consumer check (KEnRefGMX.cmake) against a pre-installed core of a different ACCEL. See core/SimdSignature.h.
add_executable(kenref_simd_check "${CMAKE_CURRENT_SOURCE_DIR}/cmake/simd_check_main.cpp")
target_link_libraries(kenref_simd_check PRIVATE kenref_core)
target_include_directories(kenref_simd_check PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include_core")

# Skip the runtime gate when cross-compiling (the check exe can't run on the build host).
if(NOT CMAKE_CROSSCOMPILING)
    # rc 0 = match, rc 1 = genuine mismatch (simd_check_main's only failure return). ANY OTHER rc means the
    # check could not RUN (e.g. 127: built with -stdlib=libc++ but libc++.so is not loadable here) — that is
    # not evidence of a mismatch, so warn instead of falsely accusing the build and refusing to install.
    install(CODE "
        message(STATUS \"KEnRef: verifying kenref_core SIMD/Eigen ABI ...\")
        execute_process(COMMAND \"$<TARGET_FILE:kenref_simd_check>\" RESULT_VARIABLE _kn_simd_rc)
        if(_kn_simd_rc EQUAL 1)
            message(FATAL_ERROR \"kenref_core SIMD/Eigen ABI self-check FAILED — refusing to install an ISA-mismatched build.\")
        elseif(NOT _kn_simd_rc EQUAL 0)
            message(WARNING \"KEnRef: could not RUN the SIMD/Eigen ABI self-check (rc=\${_kn_simd_rc}); skipping it. \"
                            \"This is not a mismatch — commonly the compiler's runtime libs are not loadable \"
                            \"(e.g. -stdlib=libc++ without LD_LIBRARY_PATH).\")
        endif()
    ")
endif()

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
