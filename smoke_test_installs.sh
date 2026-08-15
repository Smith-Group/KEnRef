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
# The two sides use DIFFERENT gromacs (the kenref one vs the plumed-batched one),
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
# Usage:  ./smoke_test_installs.sh [staging|deploy]
#         TIER=deploy KENREF_VER=2.0.0 ./smoke_test_installs.sh
# =============================================================================
REPO=/home/amr/CLionProjects/KEnRef
NSTEPS=${NSTEPS:-500}
GUIDE_ATOMS=39,60,82,101,117,234,241,256,270,284,358,368,383,397,407,422,444,454,474,496,513,534,544,558,626,642,666,680,701,771,785
REF=$REPO/res/sigma/md_double/GB3_27_10us.pdb   # same GB3 system for every set

TIER="${1:-${TIER:-staging}}"
OPT=/smithlab/opt
ACCEL="${ACCEL:-AVX_512}"
KENREF_VER="${KENREF_VER:-2.0.0}"        # the kenref version dimension / name token
GMX_VER="${GMX_VER:-2025.4}"
PLUMED_VER="${PLUMED_VER:-master}"       # staging only; deploy uses a describe token (see below)

MODINIT='source /usr/share/Modules/init/bash'
OMPI='module load openmpi/5.0.7_clang20.1.1_cuda12.4.1; module load llvm/20.1.1'

case "$TIER" in
staging)
    # Hierarchical: <base>-dev/<version dims>/<buildtype>/<accel>
    KENREF_GMXRC="$OPT/gromacs-4-kenref-dev/${GMX_VER%%.*}/$GMX_VER/release/$ACCEL/bin/GMXRC"
    KENREF_ENV="$OPT/kenref-gmx-dev/$KENREF_VER/${GMX_VER%%.*}/$GMX_VER/release/$ACCEL/env.sh"
    PLUMED_GMXRC="$OPT/gromacs-4-plumed-dev/${GMX_VER%%.*}/$GMX_VER/release/$ACCEL/bin/GMXRC"
    PLUMED_MODULE="plumed-dev/$PLUMED_VER/$KENREF_VER/release/$ACCEL"
    ;;
deploy)
    # Flat: the config is the directory NAME. The kenref prefix holds BOTH the core and the gmx
    # interface, and its paired GROMACS is the `-gmx-<ver>` sibling — that pairing is why they are
    # named to sort together. The plumed name carries a describe token, so it is discovered rather
    # than constructed: `2.10.1-778-gb1dc713b2` is not something to hardcode.
    KENREF_GMXRC="$OPT/kenref/${KENREF_VER}_${ACCEL}-gmx-${GMX_VER}/bin/GMXRC"
    KENREF_ENV="$OPT/kenref/${KENREF_VER}_${ACCEL}/env.sh"
    PLUMED_GMXRC="$OPT/gromacs-4-plumed/${GMX_VER}_${ACCEL}/bin/GMXRC"
    PLUMED_MODULE="plumed/$(ls -t "$OPT/modulefiles/plumed/"*"_${ACCEL}-kenref${KENREF_VER}" 2>/dev/null | head -1 | xargs -r basename)"
    ;;
*)  echo "unknown tier '$TIER' (expected: staging | deploy)" >&2; exit 1 ;;
esac

echo "tier=$TIER  accel=$ACCEL  kenref=$KENREF_VER  gromacs=$GMX_VER"
echo "  kenref env : $KENREF_ENV"
echo "  kenref gmx : $KENREF_GMXRC"
echo "  plumed gmx : $PLUMED_GMXRC"
echo "  plumed mod : $PLUMED_MODULE"
# Fail NOW with the missing path, rather than 10 confusing GROMACS errors later.
for p in "$KENREF_ENV" "$KENREF_GMXRC" "$PLUMED_GMXRC"; do
    [ -e "$p" ] || { echo "MISSING: $p  — is the '$TIER' tier built?" >&2; exit 1; }
done
[ "${PLUMED_MODULE#plumed/}" = "" ] && { echo "MISSING: no plumed modulefile for $ACCEL/$KENREF_VER in the '$TIER' tier" >&2; exit 1; }

ENV_KENREF="$MODINIT; $OMPI; source $KENREF_GMXRC; source $KENREF_ENV; export OMP_NUM_THREADS=1"
ENV_PLUMED="$MODINIT; $OMPI; source $PLUMED_GMXRC; module load $PLUMED_MODULE; export OMP_NUM_THREADS=1"

WORK=$(mktemp -d /tmp/kenref_smoke.XXXX); echo "workdir = $WORK"; echo
pass=0; fail=0
runlog() { local label=$1; shift
  echo "===== $label ====="
  if env -i HOME="$HOME" bash -c "$*" >"$WORK/$label.log" 2>&1; then echo "  PASS"; pass=$((pass+1))
  else echo "  FAIL -- tail:"; tail -n 6 "$WORK/$label.log" | sed 's/^/    | /'; fail=$((fail+1)); fi
}

# args: MODEL(SIGMA|PLATEAUS)  SET(res path)  NREP(1|2)  EXPKIND(file|folder)  EXPVAL(abs)  K
do_set() {
  local MODEL=$1 SET=$2 NREP=$3 EXPKIND=$4 EXPVAL=$5 K=$6
  local name=$(basename "$SET"); local tag="${MODEL,,}_${name}"
  # md-00_PBC-BROKEN.tpr: the protein is SPLIT across a periodic boundary (see res/PBC-BROKEN.md).
  # Deliberate: it makes this smoke test cover the make-whole path on real data.
  local d="$WORK/$tag"; local NDX="$SET/KEnRefAtomIndex.ndx" TPR="$SET/md-00_PBC-BROKEN.tpr"
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
  { echo "KENREF ..."; echo "  LABEL=kenref"; echo "  MODEL=$MODEL"; echo "  K=$K"; echo "  N=0.25";
    if [ "$EXPKIND" = folder ]; then echo "  EXP_DATA_FOLDER=$EXPVAL"; [ "$MODEL" = SIGMA ] && echo "  PROTON_MHZ=700.0";
    else echo "  EXP_DATA_FILE=$EXPVAL"; fi
    echo "  GUIDE_ATOMS=$GUIDE_ATOMS"; echo "  REF=$REF"; echo "  ATOMNAME_MAPPING=$REF";
    echo "  MAX_FORCE=999"; echo "  FIT_TO_REFERENCE"; echo "  SATURATE_FORCES"; echo "... KENREF";
    echo "PRINT ARG=kenref.bias,kenref.energy,kenref.rmsd FILE=kenref.out STRIDE=100"; } > "$d/kenref.dat"

  runlog "${tag}_kenref" "$ENV_KENREF; cd '$wd' && mpirun -np $NREP KEnRef --k $K --n 0.25 --max-force 999 --params '$d/kenref.toml' -- $MULTI -nsteps $NSTEPS -ntomp 1 -deffnm md"
  runlog "${tag}_plumed" "$ENV_PLUMED; cd '$wd' && mpirun -np $NREP gmx_mpi mdrun $MULTI -plumed '$d/kenref.dat' -nsteps $NSTEPS -ntomp 1 -deffnm md"
}

# ---- the matrix ----
do_set PLATEAUS "$REPO/res/plateau/md_single" 1 file   "$REPO/res/plateau/md_single/singleton_data_GB3_3-5_all.csv" 1e8
do_set PLATEAUS "$REPO/res/plateau/md_double" 2 file   "$REPO/res/plateau/md_double/singleton_data_GB3_all.csv"     1e8
do_set SIGMA    "$REPO/res/sigma/md_single"   1 folder "$REPO/res/sigma/md_single/"                                 1.0
do_set SIGMA    "$REPO/res/sigma/md_double"   2 folder "$REPO/res/sigma/md_double/"                                 1.0
do_set RELAX    "$REPO/res/relax/md_double"   2 folder "$REPO/res/relax/md_double/"                               1

echo; echo "==================================================="
echo "SUMMARY: $pass passed, $fail failed. Logs+outputs under $WORK/"
