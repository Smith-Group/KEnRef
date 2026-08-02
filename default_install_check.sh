#!/usr/bin/env bash
# =============================================================================
# default_install_check.sh — does the OUT-OF-THE-BOX install actually work?
#
# Not the /smithlab/opt matrix. This answers a different question: if someone
# follows the documented quick start with NO options, do they end up with
# something that runs an MD simulation — from EITHER side, and do the two
# coexist?
#
# Only documented commands are used (`./build.sh -y`, `build-and-batch.sh -y`,
# `source env.sh`, `source GMXRC`). Every smoke run happens under `env -i`, so
# nothing is propped up by this shell's PATH/LD_LIBRARY_PATH — that is exactly
# the mistake that makes a broken install look fine on the developer's machine.
# A documented environment that is NOT sufficient to run is recorded as a
# FINDING; the script never quietly patches around it.
#
# PREFIXES — each product where its own build system puts it (both sides agree
# on these since the fork's prefix unification):
#     kenref  -> /usr/local/kenref     plumed -> /usr/local
#     gromacs -> /usr/local/gromacs    (kenref-gmx's own gromacs stays in its build tree)
# Because both sides install kenref to the SAME path, the machine ends up with
# exactly one kenref — which is what makes pass 4 below a re-run and not an
# uninstall.
#
# PASSES
#   0. `kenref-plumed-master` builds at all           (one-off; that branch is
#      PR-minimal and needs --kenref-src, so it is checked once, not in the loop)
#   1. KEnRef side alone, clean system   -> smoke -> REMOVED again
#   2. PLUMED side alone, clean system   -> must BUILD its own kenref -> smoke -> kept
#   3. KEnRef side ON TOP of #2          -> smoke native AND re-smoke plumed
#   4. PLUMED side ON TOP of #3's kenref -> must REUSE it (asserted) -> smoke both
#
# Passes 2 and 4 are the same command with different starting states, and the
# script asserts they take DIFFERENT paths through ensure_kenref (build vs
# reuse). That assertion is the whole point of running it twice.
#
# Usage:
#   ./default_install_check.sh                # passes 0-4
#   ./default_install_check.sh 2 3            # selected passes
#   NSTEPS=100 ./default_install_check.sh     # shorter MD
#   KENREF_SRC=... PLUMED_GIT=... ./default_install_check.sh
#
# Needs sudo (the prefixes are system paths). Primes it once, up front.
# =============================================================================
set -uo pipefail

KENREF_SRC="${KENREF_SRC:-/home/amr/CLionProjects/KEnRef}"
PLUMED_GIT="${PLUMED_GIT:-https://github.com/Smith-Group/kenref-plumed2.git}"
PLUMED_BRANCH_MIN="${PLUMED_BRANCH_MIN:-kenref-plumed-master}"      # PR-minimal, needs --kenref-src
PLUMED_BRANCH_DEF="${PLUMED_BRANCH_DEF:-kenref-plumed-downloads}"   # the zero-option branch
KENREF_PREFIX=/usr/local/kenref
KENREF_GMX_PREFIX=/usr/local/kenref-gmx
PLUMED_PREFIX=/usr/local
GROMACS_PREFIX=/usr/local/gromacs
NSTEPS="${NSTEPS:-200}"

WORK="$(mktemp -d /tmp/kenref_definst.XXXX)"; LOGS="$WORK/logs"; mkdir -p "$LOGS"

RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'; BLUE='\033[0;34m'; NC='\033[0m'
say()  { printf "\n${BLUE}======== %s ========${NC}\n" "$*"; }
ok()   { printf "${GREEN}✓ %s${NC}\n" "$*"; }
warn() { printf "${YELLOW}WARNING: %s${NC}\n" "$*" >&2; }
bad()  { printf "${RED}✗ %s${NC}\n" "$*" >&2; }

PASSES=(0 1 2 3 4)
[ $# -gt 0 ] && PASSES=("$@")

FINDINGS=()
finding() { FINDINGS+=("$1"); bad "FINDING: $1"; }

step() {                                   # step <label> <cmd...>
    local label="$1"; shift
    printf "\n${BLUE}--- %s ---${NC}\n" "$label"
    if "$@" >"$LOGS/$label.log" 2>&1; then ok "$label"; return 0; fi
    bad "$label FAILED — $LOGS/$label.log"
    tail -n 25 "$LOGS/$label.log" | sed 's/^/    | /'
    return 1
}

# ---------------------------------------------------------------------------
# Smoke input: the smallest committed set — PLATEAUS, single member, 1 rank —
# the same fixture smoke_test_installs.sh drives, so a failure here is about the
# INSTALL and not about the data.
# ---------------------------------------------------------------------------
SET="$KENREF_SRC/res/plateau/md_single"
REF="$KENREF_SRC/res/sigma/md_double/GB3_27_10us.pdb"
EXP="$SET/singleton_data_GB3_3-5_all.csv"
GUIDE_ATOMS=39,60,82,101,117,234,241,256,270,284,358,368,383,397,407,422,444,454,474,496,513,534,544,558,626,642,666,680,701,771,785

prepare_case() {                           # prepare_case <dir>
    local d="$1"; mkdir -p "$d"; cp "$SET/md-00.tpr" "$d/md.tpr" || return 1
    { echo "model = PLATEAUS"; echo 'guide = "guideC-alpha"'
      echo "index = \"$SET/KEnRefAtomIndex.ndx\""; echo "exp-data-file = \"$EXP\""
      echo "atomname-mapping = \"$REF\""; echo "ref = \"$REF\""; } > "$d/kenref.toml"
    # NO `ARG=` — KENREF biases on coordinates, not a CV. See smoke_test_installs.sh for the full story.
    { echo "KENREF ..."; echo "  LABEL=kenref"; echo "  MODEL=PLATEAUS"
      echo "  K=1e8"; echo "  N=0.25"; echo "  EXP_DATA_FILE=$EXP"
      echo "  GUIDE_ATOMS=$GUIDE_ATOMS"; echo "  REF=$REF"; echo "  ATOMNAME_MAPPING=$REF"
      echo "  MAX_FORCE=999"; echo "  FIT_TO_REFERENCE"; echo "  SATURATE_FORCES"; echo "... KENREF"
      echo "PRINT ARG=kenref.bias,kenref.energy FILE=kenref.out STRIDE=50"; } > "$d/kenref.dat"
}

# `env -i` guarantees the install stands on its own: only the one documented
# setup line under test is allowed to contribute to the environment.
in_clean_env() { env -i HOME="$HOME" TERM="${TERM:-dumb}" bash -lc "set -e; $1; $2"; }

MODINIT='source /usr/share/Modules/init/bash 2>/dev/null || true'
env_kenref() { echo "$MODINIT; source ${KENREF_PREFIX}/env.sh"; }
# PLUMED at /usr/local puts `plumed` on the default PATH already; GROMACS ships GMXRC.
env_plumed() { echo "$MODINIT; source ${GROMACS_PREFIX}/bin/GMXRC"; }

smoke_kenref() {                           # smoke_kenref <tag>
    # SEPARATE statements: `local a="$1" b="...$a"` does NOT work — bash expands every argument to the
    # `local` builtin BEFORE running it, so $a is still unbound (fatal under `set -u`).
    local tag="$1"; local d="$WORK/kn-$tag"
    prepare_case "$d" || { finding "could not stage inputs ($tag)"; return 1; }
    step "smoke-kenref-$tag" in_clean_env "$(env_kenref)" \
        "cd '$d' && KEnRef --k 1e8 --n 0.25 --max-force 999 --params '$d/kenref.toml' -- -nsteps $NSTEPS -ntomp 1 -deffnm md"
}
smoke_plumed() {                           # smoke_plumed <tag>
    local tag="$1"; local d="$WORK/pl-$tag"        # separate statements — see smoke_kenref
    prepare_case "$d" || { finding "could not stage inputs ($tag)"; return 1; }
    step "smoke-plumed-$tag" in_clean_env "$(env_plumed)" \
        "cd '$d' && gmx_mpi mdrun -plumed '$d/kenref.dat' -nsteps $NSTEPS -ntomp 1 -deffnm md"
}

# On failure, answer the two questions that matter for an install: is the binary
# even on PATH, and are its shared libraries resolvable?
diagnose() {                               # diagnose <setup> <exe>
    printf "  ${YELLOW}diagnosis for '%s':${NC}\n" "$2"
    in_clean_env "$1" "command -v $2 || echo '(not on PATH)'"                              2>&1 | sed 's/^/    /'
    in_clean_env "$1" "ldd \$(command -v $2) 2>/dev/null | grep -i 'not found' || echo '(no missing libraries)'" 2>&1 | sed 's/^/    /'
}

report_layout() {                          # report_layout <label> <dir...>
    printf "\n${BLUE}  what landed (%s):${NC}\n" "$1"; shift
    for d in "$@"; do
        if [ -e "$d" ]; then printf "    %-46s %s\n" "$d" "$(du -sh "$d" 2>/dev/null | cut -f1)"
        else printf "    %-46s ${YELLOW}(absent)${NC}\n" "$d"; fi
    done
}

# ---------------------------------------------------------------------------
clone_plumed() {                           # clone_plumed <branch> -> echoes path
    local br="$1"; local dst="$WORK/plumed-$br"    # separate statements — see smoke_kenref
    [ -d "$dst" ] && { echo "$dst"; return 0; }
    git clone --quiet --branch "$br" "$PLUMED_GIT" "$dst" >&2 || return 1
    echo "$dst"
}

build_kenref_side() {                      # the documented KEnRef quick start
    step "build-kenref-side-$1" bash -c "cd '$KENREF_SRC' && ./build.sh --with-gmx -y"
}
build_plumed_side() {                      # <tag> <plumed tree> [extra args...]
    local tag="$1" tree="$2"; shift 2
    step "build-plumed-side-$tag" bash -c "cd '$tree' && ./src/kenref/build-and-batch.sh -y $*"
}

remove_kenref_install() {
    say "Removing the KEnRef install (so the next pass genuinely starts from nothing)"
    for d in "$KENREF_PREFIX" "$KENREF_GMX_PREFIX"; do
        [ -e "$d" ] && { sudo rm -rf "$d" && ok "removed $d"; }
    done
    # The BUILD tree is left alone: it is not an install, and wiping it would
    # force a 20-minute GROMACS rebuild in the next KEnRef pass.
    [ -e "$KENREF_PREFIX" ] && finding "$KENREF_PREFIX survived removal"
    return 0
}

# ensure_kenref() either reuses an installed kenref or builds one. Which path it
# took is the thing passes 2 and 4 disagree about, so assert it from the log.
assert_kenref_path() {                     # assert_kenref_path <log label> reuse|build
    local log="$LOGS/$1.log" want="$2" got
    grep -q "reusing it" "$log" 2>/dev/null && got=reuse || got=build
    if [ "$got" = "$want" ]; then ok "kenref was ${want}d, as expected"
    else finding "expected the PLUMED side to ${want} kenref, but it ${got}ed it (see $log)"; fi
}

# ---------------------------------------------------------------------------
say "Default end-user install check"
printf "  KEnRef source : %s\n  PLUMED repo   : %s\n" "$KENREF_SRC" "$PLUMED_GIT"
printf "  prefixes      : kenref=%s  plumed=%s  gromacs=%s\n" "$KENREF_PREFIX" "$PLUMED_PREFIX" "$GROMACS_PREFIX"
printf "  build type    : Release (the default on both sides)\n  MD length     : %s steps\n" "$NSTEPS"
printf "  passes        : %s\n  logs          : %s\n" "${PASSES[*]}" "$LOGS"

say "Preflight"
# Both routes can CLONE from the remotes, so a local-only commit would be invisible here.
if [ "$(git -C "$KENREF_SRC" rev-parse master 2>/dev/null)" = "$(git -C "$KENREF_SRC" rev-parse origin/master 2>/dev/null)" ]; then
    ok "KEnRef master == origin/master"
else
    warn "KEnRef master differs from origin/master — a clone-based route would fetch the REMOTE state."
fi
if [ -e "$GROMACS_PREFIX" ] && [ ! -d "$GROMACS_PREFIX/bin" ]; then
    finding "$GROMACS_PREFIX exists but is not a GROMACS install (old buildtype-keyed tree?) — archive and remove it first"
fi
sudo -v || { bad "sudo is required (the prefixes are system paths)"; exit 1; }

for p in "${PASSES[@]}"; do
case "$p" in
0)  say "PASS 0 — does the PR-minimal branch ($PLUMED_BRANCH_MIN) build?"
    # That branch deliberately has no auto-download, so it needs --kenref-src. Checked ONCE, with
    # build-only.sh: proving PLUMED+kenref compiles there does not require a GROMACS source too.
    tree="$(clone_plumed "$PLUMED_BRANCH_MIN")" || { finding "pass 0: could not clone $PLUMED_BRANCH_MIN"; continue; }
    step "build-plumed-min" bash -c "cd '$tree' && ./src/kenref/build-only.sh --kenref-src '$KENREF_SRC' -y" \
        || finding "pass 0: $PLUMED_BRANCH_MIN does not build"
    step "plumed-has-kenref" in_clean_env "$MODINIT" "plumed --version && plumed manual --action=KENREF >/dev/null" \
        || finding "pass 0: the installed plumed does not expose the KENREF action"
    ;;
1)  say "PASS 1 — KEnRef side alone, clean system"
    build_kenref_side p1 || finding "pass 1: the default KEnRef build failed"
    report_layout "KEnRef side" "$KENREF_PREFIX" "$KENREF_GMX_PREFIX" "$KENREF_SRC/cmake-build-release-orch/gromacs-install"
    smoke_kenref p1 || { finding "pass 1: a default KEnRef install cannot run an MD from 'source env.sh' alone"; diagnose "$(env_kenref)" KEnRef; }
    remove_kenref_install
    ;;
2)  say "PASS 2 — PLUMED side alone, clean system (must build its OWN kenref)"
    [ -e "$KENREF_PREFIX" ] && warn "$KENREF_PREFIX exists — this pass is meant to start without it."
    tree="$(clone_plumed "$PLUMED_BRANCH_DEF")" || { finding "pass 2: could not clone $PLUMED_BRANCH_DEF"; continue; }
    build_plumed_side p2 "$tree" || finding "pass 2: the zero-option PLUMED build failed"
    assert_kenref_path build-plumed-side-p2 build
    report_layout "PLUMED side" "$KENREF_PREFIX" "$PLUMED_PREFIX/bin/plumed" "$GROMACS_PREFIX"
    smoke_plumed p2 || { finding "pass 2: a default PLUMED install cannot run an MD from GMXRC alone"; diagnose "$(env_plumed)" gmx_mpi; }
    ;;
3)  say "PASS 3 — KEnRef side ON TOP of the PLUMED install"
    build_kenref_side p3 || finding "pass 3: the KEnRef build failed with a PLUMED install already present"
    smoke_kenref p3 || { finding "pass 3: the KEnRef install does not work alongside PLUMED"; diagnose "$(env_kenref)" KEnRef; }
    # The coexistence check: installing KEnRef over the shared /usr/local/kenref must not break PLUMED,
    # whose kernel links libkenref_core from exactly that prefix.
    smoke_plumed p3 || { finding "pass 3: PLUMED STOPPED WORKING after the KEnRef side was installed over its kenref"; diagnose "$(env_plumed)" gmx_mpi; }
    ;;
4)  say "PASS 4 — PLUMED side ON TOP of the existing kenref (must REUSE it)"
    [ -e "$KENREF_PREFIX" ] || finding "pass 4: $KENREF_PREFIX is missing — nothing to reuse (run pass 3 first)"
    tree="$(clone_plumed "$PLUMED_BRANCH_DEF")" || { finding "pass 4: could not clone $PLUMED_BRANCH_DEF"; continue; }
    build_plumed_side p4 "$tree" || finding "pass 4: the PLUMED build failed against an existing kenref"
    assert_kenref_path build-plumed-side-p4 reuse
    smoke_plumed p4 || { finding "pass 4: PLUMED does not run when built against the pre-existing kenref"; diagnose "$(env_plumed)" gmx_mpi; }
    smoke_kenref p4 || { finding "pass 4: the KEnRef side broke after PLUMED was rebuilt over the shared kenref"; diagnose "$(env_kenref)" KEnRef; }
    ;;
esac
done

say "Result"
report_layout "final state" "$KENREF_PREFIX" "$KENREF_GMX_PREFIX" "$PLUMED_PREFIX/bin/plumed" "$GROMACS_PREFIX"
echo
if [ ${#FINDINGS[@]} -eq 0 ]; then
    ok "no findings — the documented default install works from both sides, and they share one kenref."
else
    bad "${#FINDINGS[@]} finding(s):"
    for f in "${FINDINGS[@]}"; do echo "    - $f"; done
fi
echo "  logs + run dirs: $WORK"
[ ${#FINDINGS[@]} -eq 0 ]
