# ============================================================================
# PostInstall.cmake — after any KEnRef install, ALWAYS:
#   * write ${prefix}/env.sh (source-able: PATH / LD_LIBRARY_PATH / PKG_CONFIG_PATH / CMAKE_PREFIX_PATH),
#   * generate a TCL environment-modules modulefile (unless KENREF_MODULEFILE_DIR is empty),
#   * print a clear summary of WHERE libs vs executables went and WHAT to add to the environment.
#
# ONE set of core libraries: the kenref-gmx component installs no core libs (they live under the kenref
# prefix); its env therefore also references that core prefix so the executables resolve libkenref_core.
# Included from CMakeLists.txt (for core and/or gmx builds).
# ============================================================================

# --- component identity ------------------------------------------------------
if(BUILD_KENREF_GMX)
    set(_kn_comp "kenref-gmx")
else()
    set(_kn_comp "kenref")
endif()

# --- where the CORE library actually lives (never duplicated in the gmx prefix) ---------------------
if(BUILD_KENREF_CORE)
    set(_kn_core_prefix "${CMAKE_INSTALL_PREFIX}")            # core built+installed here
elseif(DEFINED KENREF_CORE_CMAKE_PATH)
    get_filename_component(_kn_core_prefix "${KENREF_CORE_CMAKE_PATH}/../../.." ABSOLUTE)   # <prefix>/lib/cmake/KEnRef
else()
    set(_kn_core_prefix "${KENREF_INSTALL_BASE}/kenref")
endif()

# --- assemble the environment path lists (dedup, keep order) ------------------
set(_kn_bins "")
set(_kn_libs "")
set(_kn_pcs  "")
set(_kn_prefixes "")

# this component's own prefix (share/pkgconfig covers a KEnRef-built Eigen shipped into the prefix)
list(APPEND _kn_prefixes "${CMAKE_INSTALL_PREFIX}")
list(APPEND _kn_libs     "${CMAKE_INSTALL_PREFIX}/lib")
list(APPEND _kn_pcs      "${CMAKE_INSTALL_PREFIX}/lib/pkgconfig" "${CMAKE_INSTALL_PREFIX}/share/pkgconfig")
if(BUILD_KENREF_GMX)
    list(APPEND _kn_bins "${CMAKE_INSTALL_PREFIX}/bin")       # KEnRef, energycalc, s2calc
    # + the CORE prefix (libkenref_core lives there, not here)
    if(NOT _kn_core_prefix STREQUAL "${CMAKE_INSTALL_PREFIX}")
        list(APPEND _kn_prefixes "${_kn_core_prefix}")
        list(APPEND _kn_libs     "${_kn_core_prefix}/lib")
        list(APPEND _kn_pcs      "${_kn_core_prefix}/lib/pkgconfig")
    endif()
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
list(REMOVE_DUPLICATES _kn_prefixes)
list(REMOVE_DUPLICATES _kn_libs)
list(REMOVE_DUPLICATES _kn_pcs)

# --- env.sh (colon-joined, source-able) --------------------------------------
string(REPLACE ";" ":" _kn_bins_s "${_kn_bins}")
string(REPLACE ";" ":" _kn_libs_s "${_kn_libs}")
string(REPLACE ";" ":" _kn_pcs_s  "${_kn_pcs}")
string(REPLACE ";" ":" _kn_pref_s "${_kn_prefixes}")

set(_kn_env "#!/bin/sh\n# KEnRef ${_kn_comp} ${PROJECT_VERSION} — `source` this file to use the installed KEnRef.\n")
if(_kn_bins)
    string(APPEND _kn_env "export PATH=\"${_kn_bins_s}:\${PATH}\"\n")
endif()
string(APPEND _kn_env "export LD_LIBRARY_PATH=\"${_kn_libs_s}:\${LD_LIBRARY_PATH}\"\n")
string(APPEND _kn_env "export PKG_CONFIG_PATH=\"${_kn_pcs_s}:\${PKG_CONFIG_PATH}\"\n")
string(APPEND _kn_env "export CMAKE_PREFIX_PATH=\"${_kn_pref_s}:\${CMAKE_PREFIX_PATH}\"\n")
file(WRITE "${CMAKE_BINARY_DIR}/kenref-env.sh" "${_kn_env}")
install(FILES "${CMAKE_BINARY_DIR}/kenref-env.sh" DESTINATION "." RENAME "env.sh")

# --- TCL modulefile (always generated unless the dir is emptied) --------------
set(KENREF_MODULEFILE_DIR "${KENREF_INSTALL_BASE}/modulefiles/${_kn_comp}"
    CACHE PATH "Directory for the generated TCL modulefile (set empty to skip module-file generation)")
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
else()
    set(_kn_mod_line "")
endif()
install(CODE "message(STATUS \"\n\
================ KEnRef ${_kn_comp} ${PROJECT_VERSION} installed ================\n\
  libraries   -> ${_kn_core_prefix}/lib   (libkenref_core.{a,so}, libkenref_and_eigen3.a)\n\
${_kn_exe_line}\
To use it:\n\
  source ${CMAKE_INSTALL_PREFIX}/env.sh\n\
${_kn_mod_line}\
  or add manually:  PATH+=[${_kn_bins_s}]  LD_LIBRARY_PATH+=${_kn_libs_s}  PKG_CONFIG_PATH+=${_kn_pcs_s}  CMAKE_PREFIX_PATH+=${_kn_pref_s}\n\
==========================================================================\")")
