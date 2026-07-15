# ============================================================================
# PostInstall.cmake — after any KEnRef install, ALWAYS:
#   * write ${prefix}/env.sh (source-able: PATH / LD_LIBRARY_PATH / PKG_CONFIG_PATH / CMAKE_PREFIX_PATH),
#   * generate a TCL environment-modules modulefile UNDER the prefix (${prefix}/modulefiles/kenref/<version>),
#     unless KENREF_MODULEFILE_DIR is emptied, and OPTIONALLY symlink it into a system modulefiles dir,
#   * print a clear summary of WHERE libs/executables went and WHAT to add to the environment.
#
# ONE common prefix: core libs + gmx executables + the plumedinterface export all live under
# CMAKE_INSTALL_PREFIX, so there is a single env.sh / modulefile for the whole install.
# Included from CMakeLists.txt (for core and/or gmx builds).
# ============================================================================

# --- component identity: one prefix -> one module name ------------------------
set(_kn_comp "kenref")

# --- assemble the environment path lists (everything under the one prefix) ----
set(_kn_bins "")
set(_kn_libs "${CMAKE_INSTALL_PREFIX}/lib")
set(_kn_pcs  "${CMAKE_INSTALL_PREFIX}/lib/pkgconfig" "${CMAKE_INSTALL_PREFIX}/share/pkgconfig")
set(_kn_prefixes "${CMAKE_INSTALL_PREFIX}")
if(BUILD_KENREF_GMX)
    list(APPEND _kn_bins "${CMAKE_INSTALL_PREFIX}/bin")       # KEnRef, energycalc, s2calc
endif()
# external eigen3.pc (if this build used a system/fetched Eigen that ships one), so a downstream PLUMED build
# that sources env.sh prefers "kenref_core + external eigen3" over the self-contained kenref_and_eigen3 flavor.
foreach(_ed "${Eigen3_DIR}" "${Eigen3_DIR}/share/pkgconfig" "${Eigen3_DIR}/lib/pkgconfig"
            "${Eigen3_DIR}/../../pkgconfig"              # <prefix>/share/eigen3/cmake -> <prefix>/share/pkgconfig
            "${EIGEN3_INCLUDE_DIR}/../share/pkgconfig" "${EIGEN3_INCLUDE_DIR}/../../share/pkgconfig")
    if(_ed AND EXISTS "${_ed}/eigen3.pc")
        list(APPEND _kn_pcs "${_ed}")
        break()
    endif()
endforeach()
list(REMOVE_DUPLICATES _kn_libs)
list(REMOVE_DUPLICATES _kn_pcs)

# --- env.sh (colon-joined, source-able) --------------------------------------
string(REPLACE ";" ":" _kn_bins_s "${_kn_bins}")
string(REPLACE ";" ":" _kn_libs_s "${_kn_libs}")
string(REPLACE ";" ":" _kn_pcs_s  "${_kn_pcs}")
string(REPLACE ";" ":" _kn_pref_s "${_kn_prefixes}")

# NB: every expansion uses ${VAR:-} (empty default). env.sh is routinely sourced by scripts running under
# `set -u` (the plumed-side build scripts use `set -euo pipefail`), where a bare ${PKG_CONFIG_PATH} on a
# machine that has never set it is an "unbound variable" error that ABORTS the caller's build.
set(_kn_env "#!/bin/sh\n# KEnRef ${_kn_comp} ${PROJECT_VERSION} — `source` this file to use the installed KEnRef.\n")
if(_kn_bins)
    string(APPEND _kn_env "export PATH=\"${_kn_bins_s}:\${PATH:-}\"\n")
endif()
string(APPEND _kn_env "export LD_LIBRARY_PATH=\"${_kn_libs_s}:\${LD_LIBRARY_PATH:-}\"\n")
string(APPEND _kn_env "export PKG_CONFIG_PATH=\"${_kn_pcs_s}:\${PKG_CONFIG_PATH:-}\"\n")
string(APPEND _kn_env "export CMAKE_PREFIX_PATH=\"${_kn_pref_s}:\${CMAKE_PREFIX_PATH:-}\"\n")
file(WRITE "${CMAKE_BINARY_DIR}/kenref-env.sh" "${_kn_env}")
install(FILES "${CMAKE_BINARY_DIR}/kenref-env.sh" DESTINATION "." RENAME "env.sh")

# --- TCL modulefile: UNDER the prefix by default -----------------------------
# Everything from the KEnRef project lives in one prefix, so the modulefile does too:
# ${prefix}/modulefiles/<comp>/<version>. Override with KENREF_MODULEFILE_DIR (empty = skip generation).
# NOT a CACHE default: a cached derived value would go STALE when the same build dir is re-configured with a
# different --prefix (the modulefile would then be written under the OLD prefix). Deriving it as a normal
# variable makes it track CMAKE_INSTALL_PREFIX; an explicit -DKENREF_MODULEFILE_DIR=... is still honoured
# (CMake defines it, so the `if(NOT DEFINED ...)` below leaves it alone). Set it empty to skip generation.
if(NOT DEFINED KENREF_MODULEFILE_DIR)
    set(KENREF_MODULEFILE_DIR "${CMAKE_INSTALL_PREFIX}/modulefiles/${_kn_comp}")
endif()
# Optionally ALSO expose it from a shared system modulefiles dir via a symlink, so `module use <link-dir>`
# finds it without adding the in-prefix path. Empty = no link (default; build.sh offers to set it to
# /usr/local/modulefiles). The link is skipped when it would coincide with the real dir (e.g. prefix=/usr/local).
set(KENREF_MODULEFILE_LINK_DIR "" CACHE PATH "System modulefiles dir to symlink the modulefile into (empty = none)")

set(_kn_mod "#%Module1.0\n")
string(APPEND _kn_mod "## KEnRef ${_kn_comp} ${PROJECT_VERSION}\n")
string(APPEND _kn_mod "module-whatis \"KEnRef ${_kn_comp} ${PROJECT_VERSION} (kinetic ensemble refinement)\"\n")
foreach(_d IN LISTS _kn_bins)
    string(APPEND _kn_mod "prepend-path PATH ${_d}\n")
endforeach()
foreach(_d IN LISTS _kn_libs)
    string(APPEND _kn_mod "prepend-path LD_LIBRARY_PATH ${_d}\n")
endforeach()
foreach(_d IN LISTS _kn_pcs)
    string(APPEND _kn_mod "prepend-path PKG_CONFIG_PATH ${_d}\n")
endforeach()
foreach(_d IN LISTS _kn_prefixes)
    string(APPEND _kn_mod "prepend-path CMAKE_PREFIX_PATH ${_d}\n")
endforeach()
file(WRITE "${CMAKE_BINARY_DIR}/kenref-modulefile" "${_kn_mod}")

if(KENREF_MODULEFILE_DIR)
    install(FILES "${CMAKE_BINARY_DIR}/kenref-modulefile"
            DESTINATION "${KENREF_MODULEFILE_DIR}" RENAME "${PROJECT_VERSION}")

    # Optional system-dir symlink: <link-dir>/<comp> -> <modulefile-dir>. Best-effort (a permission failure
    # warns, doesn't abort the install). Skipped when the link target equals the real dir.
    if(KENREF_MODULEFILE_LINK_DIR)
        set(_kn_link_comp_dir "${KENREF_MODULEFILE_LINK_DIR}/${_kn_comp}")
        if(NOT _kn_link_comp_dir STREQUAL KENREF_MODULEFILE_DIR)
            install(CODE "
                set(_lnk  \"${_kn_link_comp_dir}\")
                set(_real \"${KENREF_MODULEFILE_DIR}\")
                get_filename_component(_lnk_parent \"\${_lnk}\" DIRECTORY)
                file(MAKE_DIRECTORY \"\${_lnk_parent}\")
                if(EXISTS \"\${_lnk}\" OR IS_SYMLINK \"\${_lnk}\")
                    file(REMOVE \"\${_lnk}\")
                endif()
                execute_process(COMMAND \"${CMAKE_COMMAND}\" -E create_symlink \"\${_real}\" \"\${_lnk}\"
                                RESULT_VARIABLE _kn_lnrc)
                if(_kn_lnrc EQUAL 0)
                    message(STATUS \"modulefiles linked: \${_lnk} -> \${_real}\")
                else()
                    message(WARNING \"could not link modulefiles into \${_lnk} (permission?). Manually: \"
                        \"ln -sfn \${_real} \${_lnk}\")
                endif()
            ")
        endif()
    endif()
endif()

# --- clear post-install summary (printed at `cmake --install` time) ----------
if(BUILD_KENREF_GMX)
    set(_kn_exe_line "  executables -> ${CMAKE_INSTALL_PREFIX}/bin   (KEnRef, energycalc, s2calc)\n")
else()
    set(_kn_exe_line "")
endif()
if(KENREF_MODULEFILE_DIR)
    get_filename_component(_kn_mod_parent "${KENREF_MODULEFILE_DIR}" DIRECTORY)
    set(_kn_mod_line "  or:     module use ${_kn_mod_parent} && module load ${_kn_comp}/${PROJECT_VERSION}\n")
    if(KENREF_MODULEFILE_LINK_DIR AND NOT "${KENREF_MODULEFILE_LINK_DIR}/${_kn_comp}" STREQUAL KENREF_MODULEFILE_DIR)
        string(APPEND _kn_mod_line "  or:     module use ${KENREF_MODULEFILE_LINK_DIR} && module load ${_kn_comp}/${PROJECT_VERSION}\n")
    endif()
else()
    set(_kn_mod_line "")
endif()
install(CODE "message(STATUS \"\n\
================ KEnRef ${_kn_comp} ${PROJECT_VERSION} installed ================\n\
  libraries   -> ${CMAKE_INSTALL_PREFIX}/lib   (libkenref_core.{a,so}, libkenref_and_eigen3.a)\n\
${_kn_exe_line}\
To use it:\n\
  source ${CMAKE_INSTALL_PREFIX}/env.sh\n\
${_kn_mod_line}\
  or add manually:  PATH+=[${_kn_bins_s}]  LD_LIBRARY_PATH+=${_kn_libs_s}  PKG_CONFIG_PATH+=${_kn_pcs_s}  CMAKE_PREFIX_PATH+=${_kn_pref_s}\n\
==========================================================================\")")
