#!/usr/bin/env bash
# =============================================================================
# KEnRef / GROMACS / PLUMED install smoke test  (/smithlab/opt, AVX_512, release)
# Drives the COMMITTED in-repo test sets res/sigma/{md_single,md_double} and
# res/plateau/{md_single,md_double}. 500-step live MD per (model x member-count),
# run both from the KEnRef native force provider and from gmx mdrun -plumed.
#   single set -> 1 replica (num_models=1) ;  double set -> 2 replicas (num_models=2)
#
# Environment is set the SANCTIONED way (no manual LD_LIBRARY_PATH):
#   * openmpi + plumed  -> `module load`
#   * gromacs           -> `source .../GMXRC`
#   * kenref            -> `source .../env.sh`
# The two sides use DIFFERENT gromacs (gromacs-4-kenref vs gromacs-4-plumed), so
# each run gets its env in a fresh subshell.
# =============================================================================
REPO=/home/amr/CLionProjects/KEnRef
NSTEPS=${NSTEPS:-500}
GUIDE_ATOMS=39,60,82,101,117,234,241,256,270,284,358,368,383,397,407,422,444,454,474,496,513,534,544,558,626,642,666,680,701,771,785
REF=$REPO/res/sigma/md_double/GB3_27_10us.pdb   # same GB3 system for every set

MODINIT='source /usr/share/Modules/init/bash'
OMPI='module load openmpi/5.0.7_clang20.1.1_cuda12.4.1; module load llvm/20.1.1'
# KEnRef side env (gromacs-4-kenref GMXRC + kenref-gmx env.sh):
ENV_KENREF="$MODINIT; $OMPI; source /smithlab/opt/gromacs-4-kenref/2025/2025.4/release/AVX_512/bin/GMXRC; source /smithlab/opt/kenref-gmx/2.0/2025/2025.4/release/AVX_512/env.sh; export OMP_NUM_THREADS=1"
# PLUMED side env (gromacs-4-plumed GMXRC + module load plumed):
ENV_PLUMED="$MODINIT; $OMPI; source /smithlab/opt/gromacs-4-plumed/2025/2025.4/release/AVX_512/bin/GMXRC; module load plumed/master/2.0/release/AVX_512; export OMP_NUM_THREADS=1"

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
  local d="$WORK/$tag"; local NDX="$SET/KEnRefAtomIndex.ndx" TPR="$SET/md-00.tpr"
  local reps=""; for i in $(seq 1 "$NREP"); do r=$(printf "repl_%02d" "$i"); mkdir -p "$d/$r"; cp "$TPR" "$d/$r/md.tpr"; reps="$reps $r"; done
  local MULTI=""; [ "$NREP" -gt 1 ] && MULTI="-multidir$reps"
  local wd="$d"; [ "$NREP" -eq 1 ] && wd="$d/repl_01"

  # ---- KEnRef-side params (v2 TOML, absolute paths) ----
  { echo "model = $MODEL"; echo 'guide = "guideC-alpha"'; echo "index = \"$NDX\"";
    if [ "$EXPKIND" = folder ]; then echo "exp-data-folder = \"$EXPVAL\""; [ "$MODEL" = SIGMA ] && echo "proton-mhz = 700.0";
    else echo "exp-data-file = \"$EXPVAL\""; fi
    echo "atomname-mapping = \"$REF\""; echo "ref = \"$REF\""; } > "$d/kenref.toml"

  # ---- PLUMED-side .dat ----
  { echo "KENREF ..."; echo "  LABEL=kenref"; echo "  ARG="; echo "  MODEL=$MODEL"; echo "  K=$K"; echo "  N=0.25";
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
