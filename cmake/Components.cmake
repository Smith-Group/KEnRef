# ============================================================================
# Components.cmake — tri-state component resolution (ON / OFF / AUTO) and the
# download policy shared by every KEnRef component.
#
# A component (CORE, GMX, PLUMED interface, the delegated PLUMED build, …) is
# BUILT when it resolves to ON. The switches are tri-state STRINGs, default AUTO:
#
#   ON    — always build; if its heavy external engine is missing it is fetched
#           (when KENREF_FETCH_MISSING) or the configure FATALs (when not).
#   OFF   — never build.
#   AUTO  — build ONLY if the component's heavy external engine (GROMACS / PLUMED)
#           is ALREADY available (found or user-provided). AUTO NEVER triggers a
#           download on its own — no surprise multi-minute GROMACS / PLUMED clone.
#
# The download switch (KENREF_FETCH_MISSING, default ON) governs fetching a
# *sub-dependency* (Eigen; a stock GROMACS for a committed gmx; the PLUMED fork
# for a committed --with-plumed) once a component is committed — it does NOT make
# AUTO build a component whose engine is absent.
# ============================================================================

# kenref_read_foreign_cache(<other-project-build-dir> <CACHE_KEY> <out>)
#   Reads one entry out of ANOTHER CMake project's CMakeCache.txt (e.g. GROMACS's build dir), so we can learn
#   what that project was configured with instead of guessing. Sets ${out} to "" when unavailable.
function(kenref_read_foreign_cache build_dir key out)
    set(${out} "" PARENT_SCOPE)
    if(NOT EXISTS "${build_dir}/CMakeCache.txt")
        return()
    endif()
    file(STRINGS "${build_dir}/CMakeCache.txt" _kn_line REGEX "^${key}:[A-Za-z]+=")
    if(_kn_line)
        list(GET _kn_line 0 _kn_line)
        string(REGEX REPLACE "^${key}:[A-Za-z]+=" "" _kn_line "${_kn_line}")
        set(${out} "${_kn_line}" PARENT_SCOPE)
    endif()
endfunction()

# kenref_read_cpp_define(<header> <MACRO> <out>)
#   Reads `#define <MACRO> <value>` out of a generated C/C++ header (e.g. GROMACS's config.h), so flavor facts
#   come from the authoritative generated file rather than from a directory-name heuristic.
function(kenref_read_cpp_define header macro out)
    set(${out} "" PARENT_SCOPE)
    if(NOT EXISTS "${header}")
        return()
    endif()
    file(STRINGS "${header}" _kn_line REGEX "^#define[ \t]+${macro}[ \t]+")
    if(_kn_line)
        list(GET _kn_line 0 _kn_line)
        string(REGEX REPLACE "^#define[ \t]+${macro}[ \t]+" "" _kn_line "${_kn_line}")
        string(STRIP "${_kn_line}" _kn_line)
        set(${out} "${_kn_line}" PARENT_SCOPE)
    endif()
endfunction()

# Validate + upper-case a tri-state value. Call on the raw cache string.
function(kenref_normalize_tristate value out)
    string(TOUPPER "${value}" _v)
    if(NOT _v MATCHES "^(ON|OFF|AUTO)$")
        message(FATAL_ERROR "tri-state option must be ON, OFF or AUTO (got '${value}').")
    endif()
    set(${out} "${_v}" PARENT_SCOPE)
endfunction()

# kenref_resolve_component(<PRETTY_NAME> <tri-state> <available-bool> <out_bool>)
#   Sets ${out_bool} to ON/OFF (a real boolean usable in if()).
#   For a forced ON it also sets ${out_bool}_FORCED=ON so the provide step knows
#   the user demanded it (→ fetch-or-FATAL rather than a silent skip).
function(kenref_resolve_component pretty tri available out)
    kenref_normalize_tristate("${tri}" tri)
    unset(${out}_FORCED PARENT_SCOPE)
    if(tri STREQUAL "OFF")
        set(${out} OFF PARENT_SCOPE)
        message(STATUS "  ${pretty}: OFF")
    elseif(tri STREQUAL "ON")
        set(${out} ON PARENT_SCOPE)
        set(${out}_FORCED ON PARENT_SCOPE)
        message(STATUS "  ${pretty}: ON (forced)")
    else() # AUTO
        if(available)
            set(${out} ON PARENT_SCOPE)
            message(STATUS "  ${pretty}: AUTO -> ON (engine available)")
        else()
            set(${out} OFF PARENT_SCOPE)
            message(STATUS "  ${pretty}: AUTO -> OFF (engine not present; pass =ON to force build + provide/download it)")
        endif()
    endif()
endfunction()
