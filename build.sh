#!/usr/bin/env bash
# =============================================================================
# build.sh — KEnRef build orchestrator (Linux/macOS).
#
# KEnRef builds KEnRef; each MD engine builds itself. This script builds kenref_core (and, optionally, the
# kenref-gmx force provider) with KEnRef's own CMake, and for anything PLUMED it DELEGATES to the plumed side
# — it never reimplements plumed's or gromacs's build. EVERYTHING installs under ONE prefix (--prefix,
# default /usr/local/kenref).
#
# Components are TRI-STATE: ON / OFF / AUTO (default AUTO). AUTO builds a component only if its heavy engine
# (GROMACS / PLUMED) is ALREADY present; it never downloads one. Use =ON to force a build (fetching sub-deps
# when --download ON). See INSTALL.md.
#
#   (default)                 kenref_core (the library)
#   --with-gmx[=ON|OFF|AUTO]   kenref-gmx force provider (executables), consuming the core
#   --with-plumed[=ON|OFF|AUTO] + PLUMED (kenref module) — DELEGATES to <plumed>/src/kenref/build-only.sh
#   --with-plumed-gromacs      + PLUMED + a batched GROMACS — DELEGATES to build-and-batch.sh
#
# A single CMake configure builds whatever resolves ON; a single `cmake --install` writes the one prefix
# (and runs the delegated PLUMED build via install(CODE) when requested). This script stays a shallow shell:
# it maps flags to `cmake -D…`. Toolchain from the environment: CXX / CC / CXXFLAGS.
# =============================================================================
set -euo pipefail
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ---- defaults ---------------------------------------------------------------
WITH_GMX="AUTO"                           # ON | OFF | AUTO
WITH_PLUMED="AUTO"                        # ON | OFF | AUTO
WITH_PLUMED_GROMACS=0
EXPORT_PLUMEDINTERFACE=0                   # -DKENREF_EXPORT_PLUMEDINTERFACE (plumed-side delegation sets this)
DOWNLOAD="ON"                             # ON | OFF  (KENREF_FETCH_MISSING)
ASSUME_YES=0
BUILD_TYPE="Release"
ACCEL=""                                  # empty => CMake/GROMACS auto-detect this machine's SIMD
JOBS="$(command -v nproc >/dev/null 2>&1 && nproc || echo 4)"
PREFIX="/usr/local/kenref"                # the ONE install prefix for everything
MODULEFILE_DIR=""                         # empty => CMake default (<prefix>/modulefiles/kenref)
MODULEFILE_LINK=""                        # non-empty => also symlink modulefiles into this system dir
# PLUMED (delegation): clone coords + optional local tree / install prefix
PLUMED_SRC=""       PLUMED_PREFIX=""
GROMACS_SRC=""                             # for --with-plumed-gromacs: a provided 2025.x source (else plumed fetches)
# kenref-gmx force-provider GROMACS: give its BUILD dir (source+install+MPI flavor are derived from it), or a
# SOURCE for KEnRef to build; give neither and KEnRef fetches+builds a stock gromacs (--download ON).
GMX_GROMACS_BUILD=""  GMX_GROMACS_SRC=""
EXTRA_CMAKE=()      # anything after `--` -> forwarded to cmake (e.g. -G Ninja, extra -D...)

say()  { printf '\n\033[1;34m==> %s\033[0m\n' "$*"; }
warn() { printf '\033[1;33mWARNING: %s\033[0m\n' "$*" >&2; }
die()  { printf '\033[1;31mERROR: %s\033[0m\n' "$*" >&2; exit 1; }

# normalize a tri-state token to upper-case; validate ON/OFF/AUTO
tri() { local v; v="$(printf '%s' "$1" | tr '[:lower:]' '[:upper:]')"; case "$v" in ON|OFF|AUTO) printf '%s' "$v";; *) die "expected ON|OFF|AUTO, got '$1'";; esac; }
onoff() { local v; v="$(printf '%s' "$1" | tr '[:lower:]' '[:upper:]')"; case "$v" in ON|OFF) printf '%s' "$v";; *) die "expected ON|OFF, got '$1'";; esac; }

usage() {
    cat <<EOF
Usage: ./build.sh [options]        (run with no options for interactive mode)

Components (tri-state ON/OFF/AUTO, default AUTO; bare flag = ON):
  (default)                  build kenref_core (the library)
  --with-gmx[=ON|OFF|AUTO]   also build the kenref-gmx force provider (executables), consuming the core
  --with-plumed[=ON|OFF|AUTO] also build PLUMED with the kenref module   (delegated to the plumed side)
  --with-plumed-gromacs      also build PLUMED + a batched GROMACS        (delegated to build-and-batch.sh)
  --export-plumedinterface   install the plumedinterface source + kenref_plumed.pc (the plumed side sets this
                             when it delegates the kenref build; enables compiling the module in-tree)

  AUTO builds a component only if its engine (GROMACS/PLUMED) is already present; it never downloads one.
  Use =ON to force it (missing sub-deps are then fetched when --download ON, else the configure FATALs).

Config:
  --download ON|OFF      fetch missing sub-deps (Eigen/GROMACS/PLUMED)   (default: ${DOWNLOAD})
  --build-type T         Release | Debug | RelWithDebInfo                (default: ${BUILD_TYPE})
  --accel A              GROMACS GMX_SIMD tier (AVX_512/AVX2_256/AVX_256/AVX2_128/AVX_128_FMA/SSE4.1/…)
                                                                          (default: auto-detect this machine)
  --jobs N               parallel build jobs                             (default: ${JOBS})
  --prefix DIR           the ONE install prefix for everything           (default: ${PREFIX}; sudo if needed)
  --modulefile-dir DIR   where the TCL modulefile is written             (default: <prefix>/modulefiles/kenref)
  --link-modulefiles[=DIR]  also symlink the modulefile into a system dir (default DIR: /usr/local/modulefiles)
  --gmx-gromacs-build DIR    GROMACS's own CMake BUILD dir -- the dir you ran `cmake -B` into for GROMACS
                             (it holds CMakeCache.txt + src/include/config.h). NOT the install prefix and NOT
                             the source: KEnRef derives both of those, and the MPI flavor, from it.
                             This is normally the ONLY gromacs option you need.
  --gmx-gromacs-src DIR      GROMACS source checkout, to have KEnRef BUILD gromacs for you (no build dir yet).
                             Omit both to fetch+build a stock gromacs (needs --download ON).
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
        --with-gmx)             WITH_GMX="ON" ;;
        --with-gmx=*)           WITH_GMX="$(tri "${1#*=}")" ;;
        --with-plumed)          WITH_PLUMED="ON" ;;
        --with-plumed=*)        WITH_PLUMED="$(tri "${1#*=}")" ;;
        --with-plumed-gromacs)  WITH_PLUMED_GROMACS=1 ;;
        --export-plumedinterface) EXPORT_PLUMEDINTERFACE=1 ;;
        --download)             DOWNLOAD="$(onoff "$2")"; shift ;;
        --download=*)           DOWNLOAD="$(onoff "${1#*=}")" ;;
        --build-type)           BUILD_TYPE="$2"; shift ;;
        --accel)                ACCEL="$2"; shift ;;
        --jobs)                 JOBS="$2"; shift ;;
        --prefix)               PREFIX="$2"; shift ;;
        --modulefile-dir)       MODULEFILE_DIR="$2"; shift ;;
        --link-modulefiles)     MODULEFILE_LINK="/usr/local/modulefiles" ;;
        --link-modulefiles=*)   MODULEFILE_LINK="${1#*=}" ;;
        --gmx-gromacs-build)    GMX_GROMACS_BUILD="$2"; shift ;;
        --gmx-gromacs-src)      GMX_GROMACS_SRC="$2"; shift ;;
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
ask_tri()     { local __v=$1 __p=$2 a; read -rp "${__p} [on/off/AUTO] " a; [ -n "$a" ] && printf -v "$__v" '%s' "$(tri "$a")"; }

if interactive; then
    say "Interactive mode (pass flags or -y to skip). Building kenref_core; component choices (ON/OFF/AUTO):"
    ask_tri WITH_GMX    "  kenref-gmx force provider (executables)?"
    ask_tri WITH_PLUMED "  PLUMED with the kenref module?"
    if [ "$WITH_PLUMED" != "OFF" ] && ask_yesno "  ... and also batch a GROMACS with it?"; then WITH_PLUMED_GROMACS=1; fi
    ask_yesno "  Allow downloading missing sub-deps (Eigen/GROMACS/PLUMED)?" && DOWNLOAD="ON" || DOWNLOAD="OFF"
    # Modulefile system-dir link: only offer when the prefix is NOT /usr/local (there the modulefile is
    # already under /usr/local/modulefiles, so a link would be a self-referential no-op).
    if [ "$PREFIX" != "/usr/local" ] && [ -z "$MODULEFILE_LINK" ]; then
        ask_yesno "  Symlink modulefiles into /usr/local/modulefiles?" && MODULEFILE_LINK="/usr/local/modulefiles"
    fi
fi
[ "$WITH_PLUMED_GROMACS" = 1 ] && WITH_PLUMED="ON"   # batching implies building plumed

bt="$(echo "$BUILD_TYPE" | tr '[:upper:]' '[:lower:]')"

# install with sudo only when the target prefix isn't writable
do_install() { # do_install <build-dir> <prefix>
    local bd=$1 pfx=$2 parent="$2"
    while [ -n "$parent" ] && [ ! -e "$parent" ]; do parent="$(dirname "$parent")"; done
    if [ -w "$parent" ]; then cmake --install "$bd"; else say "prefix ${pfx} not writable -> installing with sudo"; sudo cmake --install "$bd"; fi
}

# =============================================================================
# ONE configure / build / install — CMake builds whatever resolves ON and writes the single prefix.
# =============================================================================
say "kenref -> ${PREFIX}   (gmx=${WITH_GMX} plumed=${WITH_PLUMED}$([ "$WITH_PLUMED_GROMACS" = 1 ] && echo '+gromacs') download=${DOWNLOAD})"
BUILD_DIR="${REPO_ROOT}/cmake-build-${bt}-orch"
ARGS=( -DCMAKE_BUILD_TYPE="$BUILD_TYPE"
       -DBUILD_KENREF_CORE=ON
       -DBUILD_KENREF_GMX="$WITH_GMX"
       -DKENREF_WITH_PLUMED="$WITH_PLUMED"
       -DKENREF_FETCH_MISSING="$DOWNLOAD"
       -DCMAKE_INSTALL_PREFIX="$PREFIX" )
[ "$WITH_PLUMED_GROMACS" = 1 ]      && ARGS+=( -DKENREF_WITH_PLUMED_GROMACS=ON )
[ "$EXPORT_PLUMEDINTERFACE" = 1 ]   && ARGS+=( -DKENREF_EXPORT_PLUMEDINTERFACE=ON )
[ -n "$ACCEL" ]              && ARGS+=( -DACCEL="$ACCEL" )
[ -n "${CXX:-}" ]           && ARGS+=( -DCMAKE_CXX_COMPILER="$CXX" )
[ -n "${CC:-}" ]            && ARGS+=( -DCMAKE_C_COMPILER="$CC" )
[ -n "${CXXFLAGS:-}" ]      && ARGS+=( -DCMAKE_CXX_FLAGS="$CXXFLAGS" )
[ -n "$MODULEFILE_DIR" ]    && ARGS+=( -DKENREF_MODULEFILE_DIR="$MODULEFILE_DIR" )
[ -n "$MODULEFILE_LINK" ]   && ARGS+=( -DKENREF_MODULEFILE_LINK_DIR="$MODULEFILE_LINK" )
# kenref-gmx force-provider GROMACS (provide all three, or CMake fetches+builds when --with-gmx is committed)
[ -n "$GMX_GROMACS_BUILD" ]   && ARGS+=( -DGROMACS_BUILD_DIR="$GMX_GROMACS_BUILD" )   # src+install derived from it
[ -n "$GMX_GROMACS_SRC" ]     && ARGS+=( -DGROMACS_SRC_DIR="$GMX_GROMACS_SRC" )
# PLUMED delegation coords (CMake clones + invokes the plumed script at install time)
[ -n "$PLUMED_SRC" ]    && ARGS+=( -DKENREF_PLUMED_SRC_DIR="$PLUMED_SRC" )
[ -n "$PLUMED_PREFIX" ] && ARGS+=( -DKENREF_PLUMED_INSTALL_DIR="$PLUMED_PREFIX" )
[ -n "$GROMACS_SRC" ]   && ARGS+=( -DKENREF_GROMACS_SRC_DIR="$GROMACS_SRC" )
[ ${#EXTRA_CMAKE[@]} -gt 0 ] && ARGS+=( "${EXTRA_CMAKE[@]}" )   # user passthrough (after --)

cmake -S "$REPO_ROOT" -B "$BUILD_DIR" "${ARGS[@]}"
cmake --build "$BUILD_DIR" -j "$JOBS"
# NOTE: if --with-plumed[-gromacs], this install ALSO runs the delegated PLUMED build (CMake install(CODE)).
do_install "$BUILD_DIR" "$PREFIX"

say "DONE."
[ -f "${PREFIX}/kenref-build-manifest.txt" ] && cat "${PREFIX}/kenref-build-manifest.txt"
