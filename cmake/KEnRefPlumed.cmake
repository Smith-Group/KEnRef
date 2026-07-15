# ============================================================================
# KENREF PLUMED INTERFACE (kenref_plumed)
# ============================================================================
# KEnRefBias — the PLUMED consumer of kenref_core (a kenref::EngineAdapter driving KEnRefDriver).
#
# Integration model (Option A, model-abstraction restructure): in PRODUCTION the PLUMED fork's frozen
# frame (src/kenref/) compiles src/plumedinterface/KEnRefBias.cpp *within PLUMED's own build* (where
# PLUMED's config + headers live). This repo target therefore exists for OFFLINE COMPILE-VALIDATION:
# it builds libkenref_plumed.a against PLUMED's headers so CI / local builds can catch API breakage
# without a full PLUMED rebuild. Undefined PLUMED symbols are expected (resolved when linked into the
# plumed binary), so this is a STATIC/archive target only — never linked into an executable here.
#
# Requires PLUMED headers (KEnRefBias derives from PLMD::bias::Bias / ActionAtomistic). Set
# -DPLUMED_SRC_DIR=/path/to/plumed2/src ; only built when BUILD_KENREF_PLUMED=ON.

message(STATUS "Configuring KENREF PLUMED interface (compile-validation lib)...")

# No hardcoded PLUMED location: compile-validation needs a PLUMED source. This file is only included when
# BUILD_KENREF_PLUMED RESOLVED to ON — and AUTO resolves to OFF when the headers are absent (see
# cmake/Components.cmake). So reaching here without headers means the user FORCED =ON: fail firmly rather than
# silently skipping what they explicitly asked for. (Use =AUTO for the graceful "build it only if possible".)
if(NOT DEFINED PLUMED_SRC_DIR OR PLUMED_SRC_DIR STREQUAL "")
    message(FATAL_ERROR
        "BUILD_KENREF_PLUMED=ON requires PLUMED headers, but no -DPLUMED_SRC_DIR=<plumed>/src was given.\n"
        "  Provide them:  -DPLUMED_SRC_DIR=/path/to/plumed2/src\n"
        "  OR use AUTO:   -DBUILD_KENREF_PLUMED=AUTO  (skips when PLUMED headers are absent)")
endif()

if(NOT EXISTS "${PLUMED_SRC_DIR}/bias/Bias.h")
    message(FATAL_ERROR
        "BUILD_KENREF_PLUMED=ON but PLUMED headers were not found at PLUMED_SRC_DIR='${PLUMED_SRC_DIR}' "
        "(no bias/Bias.h).\n  Point it at <plumed>/src, or use -DBUILD_KENREF_PLUMED=AUTO to skip gracefully.")
endif()

file(GLOB plumedinterface_sources CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/plumedinterface/*.cpp")
file(GLOB plumedinterface_headers CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/include_plumedinterface/plumedinterface/*.h")

add_library(kenref_plumed STATIC ${plumedinterface_sources} ${plumedinterface_headers})

target_include_directories(kenref_plumed PRIVATE
    "${PLUMED_SRC_DIR}"                                           # PLMD headers (bias/, core/, tools/, config/)
    "${CMAKE_CURRENT_SOURCE_DIR}/include_core"
    "${CMAKE_CURRENT_SOURCE_DIR}/include_plumedinterface"
)

# Pull in kenref_core's usage requirements (Eigen, KENREF_ENABLE_* defs, include dirs). KENREF_DOUBLE is
# already applied globally (top CMakeLists add_compile_options). PLUMED runtime symbols stay unresolved
# by design (archive only), so don't add an executable/link step.
target_link_libraries(kenref_plumed PUBLIC KENREF::CORE)

add_library(KENREF::PLUMED ALIAS kenref_plumed)
message(STATUS "  kenref_plumed compiles against PLUMED headers at: ${PLUMED_SRC_DIR}")
