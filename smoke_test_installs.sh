#!/usr/bin/env bash
# =============================================================================
# KEnRef / GROMACS / PLUMED install smoke test  (/smithlab/opt, AVX_512, release)
# Drives the COMMITTED in-repo test sets res/{sigma,plateau}/{md_single,md_double}
# and res/relax/md_double. 500-step live MD per (model x member-count), run both
# from the KEnRef native force provider and from gmx mdrun -plumed.
#   single set -> 1 replica (num_models=1) ;  double set -> 2 replicas (num_models=2)
#
# Environment is set the SANCTIONED way (no manual LD_LIBRARY_PATH):
#   * openmpi + plumed  -> `module load`
#   * gromacs           -> `source .../GMXRC`
#   * kenref            -> `source .../env.sh`
# The two sides use DIFFERENT gromacs (the kenref one vs the plumed-patched one),
# so each run gets its env in a fresh subshell.
#
# TWO TIERS, one script — they differ ONLY in where things are installed:
#   staging (default)  /smithlab/opt/<pkg>-dev/...   hierarchical: a directory
#                      level per dimension, all build types, all accels
#   deploy             /smithlab/opt/<pkg>/...       flat, release-only: the whole
#                      config encoded in ONE directory NAME
# Running both is the point: the tiers are built from the same sources, so a
# result that differs between them is about the INSTALL LAYOUT, not the science.
#
# THE PLUMED SIDE'S PERIODIC REPAIR IS NOW KEnRef'S OWN. Since plumed/plumed2#1515
# (PLUMED >= 1e03000ea) KEnRefBias calls makeRequestedAtomsWhole() and reads the
# repaired coordinates, so the generated .dat deliberately carries NO MOLINFO and
# NO WHOLEMOLECULES: what is under test is OUR repair, and PLUMED's own make-whole
# above the KENREF line would mask it. Until #1515 those two lines WERE required --
# the repair was compiled but never called, and without them the SIGMA and RELAX
# sets refused on the PLUMED side because md-00_PBC-BROKEN.tpr is torn. (PLATEAUS
# passed either way, because its restrained set contains none of the 12 displaced
# atoms, which makes its green uninformative about periodicity -- so it is not the
# set to watch here.) Against a PLUMED older than #1515 the SIGMA/RELAX sets will
# refuse: that is the correct answer, and the fix is a newer PLUMED, NOT adding
# MOLINFO/WHOLEMOLECULES back. rt_validate_plumed.sh asserts the kernel version
# outright; see its header, and res/PBC-BROKEN.md.
#
# WHAT THIS ASSERTS. Each (model x member-count) is checked THREE ways on the
# step-0 energy: native-vs-reference, plumed-vs-reference, and the two engines
# against each other. Until 2026-08 this script compared NOTHING -- it scored
# runs on exit status alone, so "8/8 pass" meant eight runs finished, not that
# anything agreed. The summary now separates `runs` from `checks` precisely so
# that number can never be read as engine agreement again.
#
# Only the STEP-0 energy is a valid assertion. md.trr is not (kenref-dd: 0/4-4/4
# over 5 repeats per configuration, worse under DD), and neither is the step-N
# energy -- it moves in lockstep with md.trr. The lockstep is the informative
# part: KEnRef is a deterministic function of the coordinates handed to it, so
# that drift arrives through GROMACS's positions, not ours. Step 0 was 5/5 in
# every configuration tested.
#
# Usage:  ./smoke_test_installs.sh [staging|deploy]
#         TIER=deploy KENREF_VER=2.0.0 ./smoke_test_installs.sh
# =============================================================================
# The repository THIS script came from, so a copy run from a clone drives that clone's
# res/ sets. It used to be hardcoded to /home/amr/CLionProjects/KEnRef, which meant an
# edited copy in a temp clone silently tested the MAIN CHECKOUT's fixtures -- the same
# class of mistake as a stale install faking correctness. REPO= still overrides, which is
# what a frozen copy outside the tree needs (never run a long job from a file you are editing).
REPO="${REPO:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
# 500, not 10: nstlist is 40, so a 10-step run never leaves the initial pair-list
# partition -- which is exactly how a step-120 abort stayed invisible through every
# earlier DD/PBC check.
NSTEPS=${NSTEPS:-500}
REF=$REPO/res/sigma/md_double/GB3_27_10us.pdb   # same GB3 system for every set


# GUIDE_ATOMS used to be a hardcoded list of 31 serials here. The `[ guideC-alpha ]`
# group that the KEnRef side reads BY NAME out of the same .ndx has 33 (805 and 819
# were missing), so the two engines were not even fitting on the same atoms -- any
# energy comparison would have failed for that reason alone. Derive it, so the two
# sides cannot drift apart again. Same extraction as rt_validate_plumed.sh.
guide_group() {   # ndx-file -> comma-separated serials
    awk '/^\[ guideC-alpha \]/{f=1;next} /^\[/{f=0} f' "$1" \
        | tr -s ' \n' '\n' | grep -v '^$' | paste -sd,
}

TIER="${1:-${TIER:-staging}}"
OPT=/smithlab/opt
ACCEL="${ACCEL:-AVX_512}"
KENREF_VER="${KENREF_VER:-2.0.0}"        # the kenref version dimension / name token
# GROMACS versions to exercise, IN ORDER. A list, not a single default, because the
# tiers are append-only and hold more than one series at a time — testing one and
# calling it "the" version hid the fact that nothing else was ever checked.
#
# Versions with no matching prefix on disk are SKIPPED LOUDLY (named, with the path
# that was missing) rather than failing pre-flight. That turns "this tree has no 2026"
# from a crash into data, and makes it harmless to leave future versions listed here
# before they exist. A SILENT skip would be the very failure mode this script has been
# cleaned of — never make it quiet.
#
# !! DO NOT BUILD A 2026 gromacs-4-plumed YET. PLUMED's 2026 patch silently drops the
# !! GREX/multisim block (hunk #3 fails, patch.sh never checks patch's exit status), so
# !! a 2026 tree built today is quietly wrong. Listing 2026 here is fine — it only ever
# !! skips until such a tree exists — but do not create one to satisfy the list.
GMX_VERS="${GMX_VERS:-2025.4 2026.0 2026.3}"
# A single GMX_VER= still works, and wins: the one-version invocation stays available.
[ -n "${GMX_VER:-}" ] && GMX_VERS="$GMX_VER"
PLUMED_VER="${PLUMED_VER:-}"             # empty = DISCOVER it; see pick_plumed_ver below.
# There is no sane hardcoded default: both tiers name PLUMED by a describe token
# (`main.1e03000ea`), which changes at every re-pin -- the old default "master" named
# nothing in either tier once the pipeline landed, and the script could only be run by
# passing PLUMED_VER= by hand. So discover it, and REFUSE when the answer is ambiguous:
# picking the newest by mtime is an accident of write order, not a version decision.
pick_plumed_ver() {   # <glob-dir> <pattern> -> echoes the one match, or fails
    local dir="$1" pat="$2" n m
    m=$(cd "$dir" 2>/dev/null && ls -d $pat 2>/dev/null)
    n=$(printf '%s\n' "$m" | grep -c .)
    if [ "$n" -eq 0 ]; then
        echo "no plumed install matching $dir/$pat" >&2; return 1
    elif [ "$n" -gt 1 ]; then
        echo "several plumed installs match $dir/$pat:" >&2
        printf '  %s\n' $m >&2
        echo "  set PLUMED_VER=<token> to choose; it is not this script's call." >&2
        return 1
    fi
    printf '%s\n' "$m"
}

MODINIT='source /usr/share/Modules/init/bash'
OMPI='module load openmpi/5.0.7_clang20.1.1_cuda12.4.1; module load llvm/20.1.1'

# Resolve the four paths for ONE gromacs version. Returns non-zero (with a reason on
# stderr) if this version is not present in this tier, so the caller can skip it.
tier_paths() {   # tier_paths <gmxver>
  local GMX_VER="$1"
  case "$TIER" in
  staging)
    # Hierarchical: <base>-dev/<version dims>/<buildtype>/<accel>
    KENREF_GMXRC="$OPT/gromacs-4-kenref-dev/${GMX_VER%%.*}/$GMX_VER/release/$ACCEL/bin/GMXRC"
    KENREF_ENV="$OPT/kenref-gmx-dev/$KENREF_VER/${GMX_VER%%.*}/$GMX_VER/release/$ACCEL/env.sh"
    PLUMED_GMXRC="$OPT/gromacs-4-plumed-dev/${GMX_VER%%.*}/$GMX_VER/release/$ACCEL/bin/GMXRC"
    local pv="$PLUMED_VER"
    [ -n "$pv" ] || pv=$(pick_plumed_ver "$OPT/modulefiles/plumed-dev" \
                             "*/$KENREF_VER/release/$ACCEL") || return 1
    pv="${pv%%/*}"
    PLUMED_MODULE="plumed-dev/$pv/$KENREF_VER/release/$ACCEL"
    ;;
  deploy)
    # Flat: the config is the directory NAME. Identity is carried by PACKAGE + name, not by a
    # composite token inside one name (the pre-2026-10-19 tree put <kenref-ver>+<gmx-ver> in a
    # single `kenref/` prefix; that tree is now under $OPT/bkp/pre-pipeline-20261019/):
    #   kenref/<k>_<accel>                      the core alone
    #   kenref-gmx/<k>_<accel>-gmx-<g>          the interface, which names its GROMACS
    #   plumed/<p>_<accel>-kenref<k>            one PLUMED build, which names its KEnRef
    # So do NOT reintroduce a `<k>+<g>` token: the interface already says which GROMACS it was
    # compiled against, and one PLUMED serves any GROMACS through its gromacs-4-plumed sibling.
    # The interface prefix has NO GMXRC of its own -- GROMACS lives in gromacs-4-kenref/<g>_<accel>,
    # which is what its modulefile loads, so this script sources that GMXRC and the interface env.sh.
    # The plumed name carries a describe token, so it is DISCOVERED rather than constructed:
    # `2.10.1-787-g3eec8d7ae` is not something to hardcode, and the commit-count is not even monotonic
    # across a history rewrite. Exactly one match is required -- see pick_plumed_ver.
    KENREF_GMXRC="$OPT/gromacs-4-kenref/${GMX_VER}_${ACCEL}/bin/GMXRC"
    KENREF_ENV="$OPT/kenref-gmx/${KENREF_VER}_${ACCEL}-gmx-${GMX_VER}/env.sh"
    PLUMED_GMXRC="$OPT/gromacs-4-plumed/${GMX_VER}_${ACCEL}/bin/GMXRC"
    local pm
    if [ -n "$PLUMED_VER" ]; then
        pm="${PLUMED_VER}_${ACCEL}-kenref${KENREF_VER}"
        [ -e "$OPT/modulefiles/plumed/$pm" ] || { echo "no such plumed modulefile: $pm" >&2; return 1; }
    else
        pm=$(pick_plumed_ver "$OPT/modulefiles/plumed" "*_${ACCEL}-kenref${KENREF_VER}") || return 1
    fi
    PLUMED_MODULE="plumed/$pm"
    ;;
  *)  echo "unknown tier '$TIER' (expected: staging | deploy)" >&2; exit 2 ;;
  esac

  # Report the FIRST missing path by name. Ten confusing GROMACS errors later is the
  # alternative, and that is how the absent tree looked like a science problem.
  local p
  for p in "$KENREF_ENV" "$KENREF_GMXRC" "$PLUMED_GMXRC"; do
      [ -e "$p" ] || { echo "no such path: $p" >&2; return 1; }
  done
  [ "${PLUMED_MODULE#plumed/}" = "" ] && { echo "no plumed modulefile for $ACCEL/$KENREF_VER" >&2; return 1; }
  ENV_KENREF="$MODINIT; $OMPI; source $KENREF_GMXRC; source $KENREF_ENV; export OMP_NUM_THREADS=1"
  ENV_PLUMED="$MODINIT; $OMPI; source $PLUMED_GMXRC; module load $PLUMED_MODULE; export OMP_NUM_THREADS=1"
  return 0
}

WORK=$(mktemp -d /tmp/kenref_smoke.XXXX); echo "workdir = $WORK"; echo

# Two independent tallies. `runs` only ever meant "the process exited 0"; keeping it
# under its own name is what stops it being quoted as agreement.
runs_ok=0; runs_bad=0; checks_ok=0; checks_bad=0

runlog() { local label=$1; shift
  echo "===== $label ====="
  if env -i HOME="$HOME" bash -c "$*" >"$WORK/$label.log" 2>&1; then echo "  ran OK"; runs_ok=$((runs_ok+1))
  else echo "  RUN FAILED -- tail:"; tail -n 6 "$WORK/$label.log" | sed 's/^/    | /'; runs_bad=$((runs_bad+1)); fi
}

# The two engines report step 0 in DIFFERENT formats, from different code:
#   native  KEnRefForceProvider.cpp  ->  "Step: 0 Energy: <v>"      (ensemble master rank, every 10 steps)
#   plumed  KEnRefBias::calculate    ->  "  Step 0  KEnRef energy = <v>"
# Both are std::cout at default precision (6 significant digits). That shared formatting
# is what makes an EXACT string comparison meaningful between the two engines, and the
# "engines agree" check below stays exact: when these two disagree it is a real integration
# defect, not float noise. Against the TABLE, every set is exact too, except the two SIGMA
# sets: that energy is the repaired-GB3 residual (3.65e-05), whose 6th figure moves with the
# machine, so it must be within a relative 1e-05 instead. See the note in rt_validate.sh for
# the measurements: identical CI configurations gave 3.6515e-05 and 3.65151e-05 (2.7e-06 apart,
# the runner's hardware); the torn protein's 0.00125246 is 33 relative away; PLUMED's
# put-in-box rounding under DD is ~3e-08. Same quantity, same tolerance, in all four scripts.
e0_native() { grep -a -m1 -E '^Step: 0 Energy:'    "$1" | awk '{print $4}'; }
e0_plumed() { grep -a -m1 -E 'Step 0  KEnRef energy' "$1" | awk '{print $NF}'; }

KN_GB3_RTOL="1e-05"

# kn_close GOT WANT RTOL|exact -> prints |GOT-WANT|/|WANT| ("n/a" when either is not a number);
# status 0 iff that is <= RTOL, or, for "exact", iff the two strings are identical.
kn_close() {
    local out rc
    out="$(awk -v g="$1" -v w="$2" -v t="$3" 'BEGIN {
        num = "^[-+]?([0-9]+[.]?[0-9]*|[.][0-9]+)([eE][-+]?[0-9]+)?$"
        if (g !~ num || w !~ num || w + 0 == 0) { print "n/a"; exit 1 }
        d = (g - w) / w; if (d < 0) d = -d
        printf "%.1e\n", d
        exit (t == "exact" ? 0 : !(d <= t + 0)) }')"
    rc=$?
    echo "$out"
    if [ "$3" = exact ]; then [ -n "$1" ] && [ "$1" = "$2" ]; return; fi
    return $rc
}

check() {   # label got want rtol|exact
  local label=$1 got=$2 want=$3 rtol=${4:-exact} rel ok
  rel="$(kn_close "$got" "$want" "$rtol")"; ok=$?
  if [ "$ok" -eq 0 ]; then
      printf '  %-24s %-14s OK  (rel %s, %s)\n' "$label" "$got" "$rel" "$rtol"; checks_ok=$((checks_ok+1))
  else
      printf '  %-24s %-14s *** MISMATCH (expected %s; rel %s, %s)\n' "$label" "${got:-none}" "$want" "$rel" "$rtol"; checks_bad=$((checks_bad+1))
  fi
}

# args: MODEL(SIGMA|PLATEAUS|RELAX)  SET(res path)  NREP(1|2)  EXPKIND(file|folder)  EXPVAL(abs)  K  E0  RTOL(exact|<rel>)
do_set() {
  local MODEL=$1 SET=$2 NREP=$3 EXPKIND=$4 EXPVAL=$5 K=$6 E0=$7 RTOL=${8:-exact}
  local name=$(basename "$SET"); local tag="${MODEL,,}_${name}_gmx${GMX_VER_CUR}"
  # md-00_PBC-BROKEN.tpr: the protein is SPLIT across a periodic boundary (see res/PBC-BROKEN.md).
  # Deliberate: it makes this smoke test cover the make-whole path on real data.
  local d="$WORK/$tag"; local NDX="$SET/KEnRefAtomIndex.ndx" TPR="$SET/md-00_PBC-BROKEN.tpr"
  local GUIDE_ATOMS; GUIDE_ATOMS=$(guide_group "$NDX")
  [ -n "$GUIDE_ATOMS" ] || { echo "  SKIP $tag: no [ guideC-alpha ] group in $NDX" >&2; checks_bad=$((checks_bad+3)); return; }
  local reps=""; for i in $(seq 1 "$NREP"); do r=$(printf "repl_%02d" "$i"); mkdir -p "$d/$r"; cp "$TPR" "$d/$r/md.tpr"; reps="$reps $r"; done
  local MULTI=""; [ "$NREP" -gt 1 ] && MULTI="-multidir$reps"
  local wd="$d"; [ "$NREP" -eq 1 ] && wd="$d/repl_01"

  # ---- KEnRef-side params (v2 TOML, absolute paths) ----
  { echo "model = $MODEL"; echo 'guide = "guideC-alpha"'; echo "index = \"$NDX\"";
    if [ "$EXPKIND" = folder ]; then echo "exp-data-folder = \"$EXPVAL\""; [ "$MODEL" = SIGMA ] && echo "proton-mhz = 700.0";
    else echo "exp-data-file = \"$EXPVAL\""; fi
    echo "atomname-mapping = \"$REF\""; echo "ref = \"$REF\""; } > "$d/kenref.toml"

  # ---- PLUMED-side .dat ----
  # NO `ARG=`. KENREF biases on coordinates, never on a collective variable, so it takes no ARG at all.
  # The empty `ARG=` this used to emit was a workaround for PLUMED 2.11-dev briefly making ARG compulsory
  # for every Bias -- which meant the same input could not serve 2.10 (ARG must be ABSENT) and that
  # 2.11-dev (ARG must be PRESENT). The module's Bias no longer declares it, so passing it now fails with
  #   ERROR in input to action KENREF ... cannot understand the following words from the input line : ARG
  # NO MOLINFO and NO WHOLEMOLECULES, deliberately: KEnRef repairs the requested
  # atoms itself since plumed/plumed2#1515, and PLUMED's own make-whole here would
  # mask it. See the header -- do not add them back to make an older PLUMED pass.
  { echo "KENREF ..."; echo "  LABEL=kenref"; echo "  MODEL=$MODEL"; echo "  K=$K"; echo "  N=0.25";
    if [ "$EXPKIND" = folder ]; then echo "  EXP_DATA_FOLDER=$EXPVAL"; [ "$MODEL" = SIGMA ] && echo "  PROTON_MHZ=700.0";
    else echo "  EXP_DATA_FILE=$EXPVAL"; fi
    echo "  GUIDE_ATOMS=$GUIDE_ATOMS"; echo "  REF=$REF"; echo "  ATOMNAME_MAPPING=$REF";
    echo "  MAX_FORCE=999"; echo "  FIT_TO_REFERENCE"; echo "  SATURATE_FORCES"; echo "... KENREF";
    echo "PRINT ARG=kenref.bias,kenref.energy,kenref.rmsd FILE=kenref.out STRIDE=100"; } > "$d/kenref.dat"

  runlog "${tag}_kenref" "$ENV_KENREF; cd '$wd' && mpirun -np $NREP KEnRef --k $K --n 0.25 --max-force 999 --params '$d/kenref.toml' -- $MULTI -nsteps $NSTEPS -ntomp 1 -deffnm md"
  runlog "${tag}_plumed" "$ENV_PLUMED; cd '$wd' && mpirun -np $NREP gmx_mpi mdrun $MULTI -plumed '$d/kenref.dat' -nsteps $NSTEPS -ntomp 1 -deffnm md"

  # ---- the actual assertion ----
  local en ep
  en=$(e0_native "$WORK/${tag}_kenref.log")
  ep=$(e0_plumed "$WORK/${tag}_plumed.log")
  echo "----- $tag: step-0 energy -----"
  check "native vs reference" "$en" "$E0" "$RTOL"
  check "plumed vs reference" "$ep" "$E0" "$RTOL"
  # Checked separately from the two above: if both engines move together the reference is
  # stale (a build or physics change); if only one moves the two paths have diverged.
  if [ -n "$en" ] && [ "$en" = "$ep" ]; then
      printf '  %-24s %-14s OK\n' "engines agree" "$en"; checks_ok=$((checks_ok+1))
  else
      printf '  %-24s native=%s plumed=%s  *** ENGINES DISAGREE\n' "engines agree" "${en:-none}" "${ep:-none}"
      checks_bad=$((checks_bad+1))
  fi
  echo
}

# ---- the matrix ----
# The trailing column is the expected step-0 energy, WITH THE MAKE-WHOLE REPAIR ACTIVE.
# The older recorded agreements (SIGMA single 0.00125246 vs 0.001252, etc.) were both
# engines agreeing on a TORN protein -- see res/PBC-BROKEN.md. These supersede them.
# Expected step-0 energies, PER SIMD WIDTH. A narrower vector changes the order the energy
# terms are summed and the 6th significant figure moves, so one row of literals cannot serve all
# three widths -- the same reason rt_validate.sh and dd_validate.sh carry a table instead of one
# set of numbers. Measured per width, every value from a run AT that width across all four
# GROMACS releases, with one distinct value per (width,set) across those four -- which is what
# justifies keying on width alone. Nothing is carried over between widths.
#
#   set                  AVX_512      AVX2_256     AVX_256
#   plateaus_md_single   0.0018212    0.00182122   0.00182126
#   plateaus_md_double   2.87581      2.87575      2.87584
#   sigma_md_single      3.65155e-05  3.6515e-05   3.65161e-05
#   sigma_md_double      3.65155e-05  3.6515e-05   3.65161e-05
#   relax_md_double      34.5337      34.5337      34.5337
#
# The two sigma rows being equal at every width is a property of the fixtures, not a copied cell,
# and RELAX does not move with width at six printed figures. The sigma values are deliberately the
# SAME literals rt_validate.sh and rt_validate_plumed.sh assert for this fixture at each width:
# same fixture, same quantity, so they must stay in lockstep, and a divergence between the two
# tables means one of them is WRONG rather than merely different.
#
# Keyed on $ACCEL -- and here that is sound, unlike rt_validate_plumed.sh where $BIN and $PLUMED
# are independent inputs and the width has to be read off the binary that actually computes the
# number. In this script $ACCEL is the single source for BOTH engines' prefixes (KENREF_ENV,
# KENREF_GMXRC and PLUMED_MODULE are all built from it above), so the two halves cannot be
# constructed at different widths. If this script is ever changed to take the two prefixes
# separately, this table needs the same treatment as that gate.
#
# An unknown (width, set) is a hard failure, never a silent skip and never another width's
# number: borrowing AVX_512's row is how 32 of 60 smoke checks read as a product failure on
# 2026-10-14 when the references were the only thing wrong.
#
# WHEN YOU EDIT THIS TABLE, diff it MECHANICALLY against whatever you are copying from -- do not
# read it. Adjacent widths differ in the 6th significant figure only, so a cell copied from the
# neighbouring column is invisible on review: writing this table the first time put AVX2_256's
# 3.6515e-05 into the AVX_256 sigma_md_double cell, and it was caught by generating the 15 cells
# from the source table and diffing, not by anyone looking at it. The run would also have caught
# it -- that is the point of the table -- but an hour later.
expect_e0() {   # set-key -> expected step-0 energy on $ACCEL
    case "$ACCEL:$1" in
        AVX_512:plateaus_md_single)  echo 0.0018212   ;;
        AVX_512:plateaus_md_double)  echo 2.87581     ;;
        AVX_512:sigma_md_single)     echo 3.65155e-05 ;;
        AVX_512:sigma_md_double)     echo 3.65155e-05 ;;
        AVX_512:relax_md_double)     echo 34.5337     ;;
        AVX2_256:plateaus_md_single) echo 0.00182122  ;;
        AVX2_256:plateaus_md_double) echo 2.87575     ;;
        AVX2_256:sigma_md_single)    echo 3.6515e-05  ;;
        AVX2_256:sigma_md_double)    echo 3.6515e-05  ;;
        AVX2_256:relax_md_double)    echo 34.5337     ;;
        AVX_256:plateaus_md_single)  echo 0.00182126  ;;
        AVX_256:plateaus_md_double)  echo 2.87584     ;;
        AVX_256:sigma_md_single)     echo 3.65161e-05 ;;
        AVX_256:sigma_md_double)     echo 3.65161e-05 ;;
        AVX_256:relax_md_double)     echo 34.5337     ;;
        *) echo "expect_e0: no reference energy for ($ACCEL, $1)." >&2
           echo "  MEASURE it on that width -- do not borrow another width's number." >&2
           return 1 ;;
    esac
}

# Resolved BEFORE any MD runs, all five, so a missing reference stops the run instead of
# surfacing as a mismatch an hour in.
E0_PLATEAUS_SINGLE=$(expect_e0 plateaus_md_single) || exit 1
E0_PLATEAUS_DOUBLE=$(expect_e0 plateaus_md_double) || exit 1
E0_SIGMA_SINGLE=$(expect_e0 sigma_md_single)       || exit 1
E0_SIGMA_DOUBLE=$(expect_e0 sigma_md_double)       || exit 1
E0_RELAX_DOUBLE=$(expect_e0 relax_md_double)       || exit 1

run_matrix() {
do_set PLATEAUS "$REPO/res/plateau/md_single" 1 file   "$REPO/res/plateau/md_single/singleton_data_GB3_3-5_all.csv" 1e8 "$E0_PLATEAUS_SINGLE"
do_set PLATEAUS "$REPO/res/plateau/md_double" 2 file   "$REPO/res/plateau/md_double/singleton_data_GB3_all.csv"     1e8 "$E0_PLATEAUS_DOUBLE"
do_set SIGMA    "$REPO/res/sigma/md_single"   1 folder "$REPO/res/sigma/md_single/"                                 1.0 "$E0_SIGMA_SINGLE" "$KN_GB3_RTOL"
do_set SIGMA    "$REPO/res/sigma/md_double"   2 folder "$REPO/res/sigma/md_double/"                                 1.0 "$E0_SIGMA_DOUBLE" "$KN_GB3_RTOL"
do_set RELAX    "$REPO/res/relax/md_double"   2 folder "$REPO/res/relax/md_double/"                                 1   "$E0_RELAX_DOUBLE"
}

# ---- drive it once per gromacs version --------------------------------------
ran=""; skipped=""
for GMX_VER_CUR in $GMX_VERS; do
    echo "==================================================="
    if ! tier_paths "$GMX_VER_CUR" 2>"$WORK/resolve_$GMX_VER_CUR.err"; then
        echo "SKIP gromacs=$GMX_VER_CUR  (not present in the '$TIER' tier)"
        sed 's/^/    | /' "$WORK/resolve_$GMX_VER_CUR.err"
        skipped="$skipped $GMX_VER_CUR"; echo; continue
    fi
    echo "tier=$TIER  accel=$ACCEL  kenref=$KENREF_VER  gromacs=$GMX_VER_CUR"
    echo "  kenref env : $KENREF_ENV"
    echo "  kenref gmx : $KENREF_GMXRC"
    echo "  plumed gmx : $PLUMED_GMXRC"
    echo "  plumed mod : $PLUMED_MODULE"
    echo
    ran="$ran $GMX_VER_CUR"
    run_matrix
done

echo; echo "==================================================="
echo "gromacs: ran[${ran:- none}]  skipped[${skipped:- none}]"
echo "runs   : $runs_ok completed, $runs_bad failed   <- only that the process exited 0"
echo "checks : $checks_ok passed,    $checks_bad failed   <- the step-0 energy assertions"
echo "Logs+outputs under $WORK/"
# chain-smoke-then-deploy.sh already branches on this script's exit status; until now the
# script always ended on an echo and so always returned 0, leaving that branch dead.
#
# Skipping EVERY version must fail. Otherwise an empty run reports "0 failed" and exits 0
# -- a green light for having tested nothing, which is the same bug class as scoring
# PASS on exit status alone.
if [ -z "$ran" ]; then
    echo "NOTHING RAN: no listed gromacs version is present in the '$TIER' tier." >&2
    echo "  listed: $GMX_VERS" >&2
    exit 1
fi
if [ "$runs_bad" -eq 0 ] && [ "$checks_bad" -eq 0 ]; then exit 0; else exit 1; fi
