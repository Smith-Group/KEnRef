#!/usr/bin/env bash
# =============================================================================
# build.sh — KEnRef build orchestrator (Linux/macOS).
#
# KEnRef builds KEnRef; each MD engine builds itself. So this script builds
# kenref_core (and, optionally, the kenref-gmx force provider) with KEnRef's own
# CMake, and for anything PLUMED it DELEGATES to the plumed side — it never
# reimplements plumed's or gromacs's build.
#
#   (default)              kenref_core (the library)          -> --prefix (/usr/local/kenref)
#   --with-gmx             kenref-gmx force provider (exes)    -> --gmx-prefix (/usr/local/kenref-gmx)
#                          (fetches+builds a stock GROMACS into the build dir if none is provided)
#   --with-plumed          + PLUMED (kenref module)  — DELEGATES to <plumed>/src/kenref/build-only.sh
#   --with-plumed-gromacs  + PLUMED + a batched GROMACS — DELEGATES to build-and-batch.sh
#
# For --with-plumed[-gromacs] the script only sets CMake options; KEnRef's CMake itself clones PLUMED and
# invokes its build script (at install time, so it reuses the just-installed kenref). This script stays a
# shallow shell: it maps flags to `cmake -D…`. Toolchain from the environment: CXX / CC / CXXFLAGS.
# =============================================================================
set -euo pipefail
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ---- defaults ---------------------------------------------------------------
WITH_GMX=0
WITH_PLUMED=0
WITH_PLUMED_GROMACS=0
ASSUME_YES=0
BUILD_TYPE="Release"
ACCEL=""                                  # empty => CMake/GROMACS auto-detect this machine's SIMD
JOBS="$(command -v nproc >/dev/null 2>&1 && nproc || echo 4)"
PREFIX="/usr/local/kenref"                # core/library prefix (the ONE home of libkenref_*)
GMX_PREFIX="/usr/local/kenref-gmx"        # kenref-gmx executables prefix
MODULEFILE_DIR=""                         # empty => CMake default (<install-base>/modulefiles/<comp>)
# PLUMED (delegation): clone coords + optional local tree / install prefix
PLUMED_SRC=""       PLUMED_PREFIX=""
GROMACS_SRC=""                             # for --with-plumed-gromacs: a provided 2025.x source (else fetched by plumed)
# kenref-gmx force-provider GROMACS (stock): provide the 3 dirs, else KEnRef's CMake fetches+builds it.
GMX_GROMACS_SRC=""  GMX_GROMACS_BUILD=""  GMX_GROMACS_INSTALL=""
EXTRA_CMAKE=()      # anything after `--` -> forwarded to cmake (e.g. -G Ninja, extra -D...)

say()  { printf '\n\033[1;34m==> %s\033[0m\n' "$*"; }
warn() { printf '\033[1;33mWARNING: %s\033[0m\n' "$*" >&2; }
die()  { printf '\033[1;31mERROR: %s\033[0m\n' "$*" >&2; exit 1; }

usage() {
    cat <<EOF
Usage: ./build.sh [options]        (run with no options for interactive mode)

Components:
  (default)              build kenref_core (the library) only
  --with-gmx             also build the kenref-gmx force provider (executables), consuming the core
  --with-plumed          also build PLUMED with the kenref module   (delegated to the plumed side)
  --with-plumed-gromacs  also build PLUMED + a batched GROMACS       (delegated to build-and-batch.sh)

Config:
  --build-type T         Release | Debug | RelWithDebInfo      (default: ${BUILD_TYPE})
  --accel A              AVX_512 | AVX_256 | AVX2_256           (default: auto-detect this machine)
  --jobs N               parallel build jobs                    (default: ${JOBS})
  --prefix DIR           core/library install prefix            (default: ${PREFIX}; sudo if needed)
  --gmx-prefix DIR       kenref-gmx executables prefix          (default: ${GMX_PREFIX})
  --modulefile-dir DIR   where TCL modulefiles are written      (default: <install-base>/modulefiles/<comp>)
  --gmx-gromacs-src DIR      stock GROMACS source   (for --with-gmx; omit all 3 to auto fetch+build)
  --gmx-gromacs-build DIR    stock GROMACS build dir (kenref-gmx needs src+build+install)
  --gmx-gromacs-install DIR  stock GROMACS install
  --plumed-src DIR       PLUMED source tree (omit to clone the PR fork into the build dir)
  --plumed-prefix DIR    PLUMED install prefix
  --gromacs-src DIR      provided GROMACS 2025.x source for --with-plumed-gromacs (else plumed fetches it)
  -y | --yes             non-interactive; take defaults
  -h | --help
  --                     forward everything after it to cmake, e.g.  -- -G Ninja -DKENREF_EIGEN_INSTALL_DIR=/usr/local

Toolchain via environment: CXX, CC, CXXFLAGS (e.g. CXX=mpicxx CXXFLAGS="-stdlib=libc++").
EOF
}

# ---- arg parsing ------------------------------------------------------------
while [ $# -gt 0 ]; do
    case "$1" in
        --with-gmx)             WITH_GMX=1 ;;
        --with-plumed)          WITH_PLUMED=1 ;;
        --with-plumed-gromacs)  WITH_PLUMED_GROMACS=1 ;;
        --build-type)           BUILD_TYPE="$2"; shift ;;
        --accel)                ACCEL="$2"; shift ;;
        --jobs)                 JOBS="$2"; shift ;;
        --prefix)               PREFIX="$2"; shift ;;
        --gmx-prefix)           GMX_PREFIX="$2"; shift ;;
        --modulefile-dir)       MODULEFILE_DIR="$2"; shift ;;
        --gmx-gromacs-src)      GMX_GROMACS_SRC="$2"; shift ;;
        --gmx-gromacs-build)    GMX_GROMACS_BUILD="$2"; shift ;;
        --gmx-gromacs-install)  GMX_GROMACS_INSTALL="$2"; shift ;;
        --plumed-src)           PLUMED_SRC="$2"; shift ;;
        --plumed-prefix)        PLUMED_PREFIX="$2"; shift ;;
        --gromacs-src)          GROMACS_SRC="$2"; shift ;;
        -y|--yes)               ASSUME_YES=1 ;;
        -h|--help)              usage; exit 0 ;;
        --)                     shift; EXTRA_CMAKE=( "$@" ); break ;;   # everything after -- goes to cmake
        *) die "unknown option '$1' (try --help)" ;;
    esac
    shift
done

interactive() { [ "$ASSUME_YES" = 0 ] && [ -t 0 ]; }
ask_yesno()   { local a; read -rp "$1 [y/N] " a; [[ "$a" =~ ^[Yy] ]]; }

if interactive && [ "$WITH_GMX" = 0 ] && [ "$WITH_PLUMED" = 0 ] && [ "$WITH_PLUMED_GROMACS" = 0 ]; then
    say "Interactive mode (pass flags or -y to skip). Building kenref_core; add components?"
    ask_yesno "Also build the kenref-gmx force provider (executables)?" && WITH_GMX=1
    if ask_yesno "Also build PLUMED with the kenref module?"; then
        if ask_yesno "  ... and also batch a GROMACS with it?"; then WITH_PLUMED_GROMACS=1; else WITH_PLUMED=1; fi
    fi
fi
[ "$WITH_PLUMED_GROMACS" = 1 ] && WITH_PLUMED=0   # the -gromacs variant supersedes plain --with-plumed
DELEGATE_PLUMED=0; { [ "$WITH_PLUMED" = 1 ] || [ "$WITH_PLUMED_GROMACS" = 1 ]; } && DELEGATE_PLUMED=1

bt="$(echo "$BUILD_TYPE" | tr '[:upper:]' '[:lower:]')"

# install with sudo only when the target prefix isn't writable
do_install() { # do_install <build-dir> <prefix>
    local bd=$1 pfx=$2 parent="$2"
    while [ -n "$parent" ] && [ ! -e "$parent" ]; do parent="$(dirname "$parent")"; done
    if [ -w "$parent" ]; then cmake --install "$bd"; else say "prefix ${pfx} not writable -> installing with sudo"; sudo cmake --install "$bd"; fi
}

# No -G here: CMake uses your default generator (or $CMAKE_GENERATOR — e.g. "Unix Makefiles" or "Ninja").
# Anything after `--` on the command line (e.g. `-- -G Ninja -DFOO=bar`) is forwarded verbatim to cmake.
COMMON=( -DCMAKE_BUILD_TYPE="$BUILD_TYPE" )
[ -n "$ACCEL" ]          && COMMON+=( -DACCEL="$ACCEL" )
[ -n "${CXX:-}" ]        && COMMON+=( -DCMAKE_CXX_COMPILER="$CXX" )
[ -n "${CC:-}" ]         && COMMON+=( -DCMAKE_C_COMPILER="$CC" )
[ -n "${CXXFLAGS:-}" ]   && COMMON+=( -DCMAKE_CXX_FLAGS="$CXXFLAGS" )
[ -n "$MODULEFILE_DIR" ] && COMMON+=( -DKENREF_MODULEFILE_DIR="$MODULEFILE_DIR" )
[ ${#EXTRA_CMAKE[@]} -gt 0 ] && COMMON+=( "${EXTRA_CMAKE[@]}" )   # user passthrough (after --)

# =============================================================================
# PHASE 1 — kenref_core -> PREFIX (the one home of libkenref_*)
# =============================================================================
say "kenref_core -> ${PREFIX}"
CORE_DIR="${REPO_ROOT}/cmake-build-${bt}-core-orch"
CORE_ARGS=( "${COMMON[@]}" -DBUILD_KENREF_CORE=ON -DBUILD_KENREF_GMX=OFF -DBUILD_KENREF_PLUMED=OFF -DCMAKE_INSTALL_PREFIX="$PREFIX" )
# PLUMED delegation is driven by CMake (InvokePlumed.cmake invokes the plumed script at install time); we
# just set the options. CMake forces KENREF_EXPORT_PLUMEDINTERFACE and does the clone + build itself.
if [ "$WITH_PLUMED_GROMACS" = 1 ]; then
    CORE_ARGS+=( -DKENREF_WITH_PLUMED_GROMACS=ON )
    [ -n "$GROMACS_SRC" ] && CORE_ARGS+=( -DKENREF_GROMACS_SRC_DIR="$GROMACS_SRC" )
elif [ "$WITH_PLUMED" = 1 ]; then
    CORE_ARGS+=( -DKENREF_WITH_PLUMED=ON )
fi
if [ "$DELEGATE_PLUMED" = 1 ]; then
    [ -n "$PLUMED_SRC" ]    && CORE_ARGS+=( -DKENREF_PLUMED_SRC_DIR="$PLUMED_SRC" )
    [ -n "$PLUMED_PREFIX" ] && CORE_ARGS+=( -DKENREF_PLUMED_INSTALL_DIR="$PLUMED_PREFIX" )
fi
cmake -S "$REPO_ROOT" -B "$CORE_DIR" "${CORE_ARGS[@]}"
cmake --build "$CORE_DIR" -j "$JOBS"
# NOTE: if --with-plumed[-gromacs], this install ALSO runs the delegated PLUMED build (CMake install(CODE)).
do_install "$CORE_DIR" "$PREFIX"

# =============================================================================
# PHASE 2a (optional) — kenref-gmx force provider -> GMX_PREFIX  (scenario 2)
# GROMACS (find or fetch+build) is handled by KEnRef's CMake (kenref_provide_gromacs); this is just a
# cmake call. Provide a stock gromacs via --gmx-gromacs-src/-build/-install, or CMake fetches+builds one.
# =============================================================================
if [ "$WITH_GMX" = 1 ]; then
    say "kenref-gmx executables -> ${GMX_PREFIX}  (consuming core at ${PREFIX})"
    GMX_DIR="${REPO_ROOT}/cmake-build-${bt}-gmx-orch"
    GMX_ARGS=( "${COMMON[@]}" -DBUILD_KENREF_CORE=OFF -DBUILD_KENREF_GMX=ON
               -DCMAKE_INSTALL_PREFIX="$GMX_PREFIX" -DKENREF_CORE_CMAKE_PATH="${PREFIX}/lib/cmake/KEnRef" -DCMAKE_PREFIX_PATH="$PREFIX" )
    [ -n "$GMX_GROMACS_SRC" ]     && GMX_ARGS+=( -DGROMACS_SRC_DIR="$GMX_GROMACS_SRC" )
    [ -n "$GMX_GROMACS_BUILD" ]   && GMX_ARGS+=( -DGROMACS_BUILD_DIR="$GMX_GROMACS_BUILD" )
    [ -n "$GMX_GROMACS_INSTALL" ] && GMX_ARGS+=( -DGROMACS_INSTALL_DIR="$GMX_GROMACS_INSTALL" )
    cmake -S "$REPO_ROOT" -B "$GMX_DIR" "${GMX_ARGS[@]}"
    cmake --build "$GMX_DIR" -j "$JOBS"
    do_install "$GMX_DIR" "$GMX_PREFIX"
fi

say "DONE."
[ -f "${PREFIX}/kenref-build-manifest.txt" ] && cat "${PREFIX}/kenref-build-manifest.txt"
