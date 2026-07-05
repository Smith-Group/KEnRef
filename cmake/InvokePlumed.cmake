# ============================================================================
# InvokePlumed.cmake — CMake DELEGATES the PLUMED build to the plumed side, and
# invokes its script itself (execute_process), instead of build.sh doing it.
#
# The plumed build needs kenref_core INSTALLED (its `ensure_kenref` finds it via
# pkg-config and reuses it), so it must run AFTER kenref's files land at the
# install prefix. The natural CMake hook for that is install(CODE): it runs at
# `cmake --install`, after the kenref install rules. There it clones PLUMED into
# the build dir (or uses KENREF_PLUMED_SRC_DIR) and runs build-only.sh /
# build-and-batch.sh, sourcing the installed env.sh so the script sees kenref
# (and the preferred eigen flavor). No loop: kenref is already installed.
#
# Included from CMakeLists.txt when KENREF_WITH_PLUMED. (build.sh only sets -D.)
# ============================================================================
find_package(Git REQUIRED)

set(KENREF_PLUMED_GIT_URL     "https://github.com/Smith-Group/kenref-plumed2" CACHE STRING
    "PLUMED git URL to clone for the delegated build (TEMPORARY PR fork; upstream once merged)")
set(KENREF_PLUMED_GIT_TAG     "kenref-plumed-master" CACHE STRING "PLUMED branch/tag carrying the kenref module")
set(KENREF_PLUMED_SRC_DIR     ""  CACHE PATH   "Local PLUMED tree (else cloned into the build dir)")
set(KENREF_PLUMED_INSTALL_DIR ""  CACHE PATH   "PLUMED install prefix (else the plumed script's default)")
# for --with-plumed-gromacs: a provided GROMACS 2025.x source (else the plumed side fetches one)
set(KENREF_GROMACS_SRC_DIR    ""  CACHE PATH   "Provided GROMACS 2025.x source for the batched build")

if(KENREF_WITH_PLUMED_GROMACS)
    set(_kn_pl_script "build-and-batch.sh")
else()
    set(_kn_pl_script "build-only.sh")
endif()

# where the PLUMED tree will be (clone happens at install time so configure stays light)
if(KENREF_PLUMED_SRC_DIR)
    set(_kn_plumed_src "${KENREF_PLUMED_SRC_DIR}")
else()
    set(_kn_plumed_src "${CMAKE_BINARY_DIR}/plumed-src")
endif()

# args forwarded to the plumed script (build type / accel / prefixes). Space-joined for the bash command.
set(_kn_pl_args "-y --build-type ${CMAKE_BUILD_TYPE}")
if(DEFINED ACCEL AND NOT ACCEL STREQUAL "")
    string(APPEND _kn_pl_args " --accel ${ACCEL}")
endif()
if(KENREF_PLUMED_INSTALL_DIR)
    string(APPEND _kn_pl_args " --plumed-prefix ${KENREF_PLUMED_INSTALL_DIR}")
endif()
if(KENREF_WITH_PLUMED_GROMACS AND KENREF_GROMACS_SRC_DIR)
    string(APPEND _kn_pl_args " --gromacs-src ${KENREF_GROMACS_SRC_DIR}")
endif()

# SKIP-IF-ALREADY-BUILT: the plumed script installs to KENREF_PLUMED_INSTALL_DIR (or <plumed>/install/plumed,
# its default). If a plumed binary is already there, don't re-run on every `cmake --install`.
option(KENREF_PLUMED_FORCE_REBUILD "Re-run the delegated PLUMED build even if it is already installed" OFF)
if(KENREF_PLUMED_INSTALL_DIR)
    set(_kn_plumed_prefix "${KENREF_PLUMED_INSTALL_DIR}")
else()
    set(_kn_plumed_prefix "${_kn_plumed_src}/install/plumed")
endif()
if(KENREF_PLUMED_FORCE_REBUILD)
    set(_kn_force TRUE)
else()
    set(_kn_force FALSE)
endif()

message(STATUS "KENREF_WITH_PLUMED: CMake will delegate to ${_kn_plumed_src}/src/kenref/${_kn_pl_script} at install time")

# Runs at `cmake --install`, after kenref is installed. ${..} = configure-time; \${..} = install-time.
install(CODE "
if(EXISTS \"${_kn_plumed_prefix}/bin/plumed\" AND NOT ${_kn_force})
    message(STATUS \"PLUMED: already built at ${_kn_plumed_prefix} — skipping (set -DKENREF_PLUMED_FORCE_REBUILD=ON to rebuild).\")
else()
    set(_pl_src \"${_kn_plumed_src}\")
    set(_pl_script \"\${_pl_src}/src/kenref/${_kn_pl_script}\")
    if(NOT EXISTS \"\${_pl_script}\")
        message(STATUS \"PLUMED: cloning ${KENREF_PLUMED_GIT_URL}@${KENREF_PLUMED_GIT_TAG} -> \${_pl_src}\")
        execute_process(COMMAND \"${GIT_EXECUTABLE}\" clone --branch \"${KENREF_PLUMED_GIT_TAG}\" \"${KENREF_PLUMED_GIT_URL}\" \"\${_pl_src}\" RESULT_VARIABLE _kn_c)
        if(NOT _kn_c EQUAL 0)
            message(FATAL_ERROR \"PLUMED clone failed; set -DKENREF_PLUMED_SRC_DIR=<plumed tree>.\")
        endif()
    endif()
    message(STATUS \"PLUMED: delegating build -> \${_pl_script}  (reuses the just-installed kenref)\")
    execute_process(
        COMMAND bash -c \"source '\${CMAKE_INSTALL_PREFIX}/env.sh' 2>/dev/null; exec '\${_pl_script}' ${_kn_pl_args}\"
        RESULT_VARIABLE _kn_r)
    if(NOT _kn_r EQUAL 0)
        message(FATAL_ERROR \"PLUMED delegated build failed (see output above).\")
    endif()
endif()
")
