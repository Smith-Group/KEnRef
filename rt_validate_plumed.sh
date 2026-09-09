#!/usr/bin/env bash
# =============================================================================
# KEnRef runtime gate: periodic wholeness, PLUMED side (and engine agreement)
#
# WHY THIS EXISTS. rt_validate.sh pins the GROMACS engine. This is its PLUMED
# counterpart, and it exists because the two engines repair periodic breaks by
# DIFFERENT mechanisms and must still agree to the last digit:
#
#   GROMACS  walks the topology's BOND graph          (KEnRefMoleculeGraph)
#   PLUMED   walks a spanning tree over the REFERENCE (makeRequestedAtomsWhole)
#
# PLUMED is never handed a topology, so the reference structure is all it has.
# That makes the reference load-bearing in a way it is not for GROMACS, and it
# is why res/*/GB3_27_10us.pdb had to be made whole before this could work at
# all -- see res/PBC-BROKEN.md.
#
# Three cases, and the contrast between them is the whole point:
#
#   A  broken GB3 + repaired reference  -> PLUMED runs and MATCHES GROMACS
#   B  whole ubiquitin                  -> the repair must be a no-op here
#   C  broken GB3 + BROKEN reference    -> PLUMED must REFUSE, not guess
#
# A alone proves nothing: a repair that did nothing would still "run". B is
# what proves the repair does not damage good input. C is what proves the
# reference is really carrying the repair -- without it, A could be luck.
#
# NOTE ON THE FRAME. plumed driver is fed a .trr, not a .gro, deliberately.
# GRO carries 0.001 nm; on these fixtures the step-0 SIGMA energy is a near-zero
# residual (3.65e-05) and that rounding alone moves it to 3.4e-03 -- a factor of
# 90, with nothing wrong. The frame is produced by running the tpr for zero
# steps, so it is the tpr's own single-precision coordinates.
#
# Usage:  PLUMED=/path/to/plumed ./rt_validate_plumed.sh
#         BIN=/path/to/KEnRef MPIRUN=/path/to/mpirun PLUMED=... ./rt_validate_plumed.sh
# =============================================================================
set -u

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN="${BIN:-$REPO/cmake-build-cli-gmx/KEnRef}"
PLUMED="${PLUMED:-plumed}"
MPIRUN="${MPIRUN:-mpirun}"
WORK="${WORK:-${TMPDIR:-/tmp}/kenref-rt-validate-plumed.$$}"

export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export PLUMED_NUM_THREADS="${PLUMED_NUM_THREADS:-1}"

# The same step-0 energies rt_validate.sh asserts on the GROMACS side. That they
# are the SAME numbers here is the point: two engines, two repair mechanisms.
EXPECT_GB3="3.65155e-05"
EXPECT_UBQ="7.18869e-08"

GB3_SET="$REPO/res/sigma/md_single"
GB3_TPR="$GB3_SET/md-00_PBC-BROKEN.tpr"
UBQ_SET="$REPO/res/10nsstart+fitting"
UBQ_TPR="$UBQ_SET/t000.00.tpr"

[ -x "$BIN" ] || { echo "rt_validate_plumed: no KEnRef binary at $BIN (set BIN=)" >&2; exit 2; }
command -v "$PLUMED" >/dev/null || { echo "rt_validate_plumed: no plumed (set PLUMED=)" >&2; exit 2; }
command -v "$MPIRUN" >/dev/null || { echo "rt_validate_plumed: no mpirun (set MPIRUN=)" >&2; exit 2; }
"$PLUMED" driver --help 2>&1 | grep -q -- --itrr || {
    echo "rt_validate_plumed: this plumed has no --itrr driver support" >&2; exit 2; }

mkdir -p "$WORK"

# Zero-based-to-1-based is not our problem here: read the group straight out of the
# .ndx so PLUMED's GUIDE_ATOMS can never drift from what the KEnRef side uses.
guide_group() {   # ndx-file
    awk '/^\[ guideC-alpha \]/{f=1;next} /^\[/{f=0} f' "$1" \
        | tr -s ' \n' '\n' | grep -v '^$' | paste -sd,
}

# Frame 0 in FULL precision: run the tpr for zero steps and keep the trr it writes.
make_frame() {   # label tpr -> $WORK/<label>.trr
    local label="$1" tpr="$2"
    local d="$WORK/frame_$label"
    rm -rf "$d"; mkdir -p "$d"; cp "$tpr" "$d/md.tpr"
    ( cd "$d" && "$MPIRUN" --oversubscribe -np 1 "$BIN" --k 1.0 --n 0.25 --max-force 999 \
          --params "$WORK/$label.toml" -- -nsteps 0 -ntomp 1 -deffnm md ) > "$WORK/frame_$label.log" 2>&1
    [ -s "$d/md.trr" ] || { echo "rt_validate_plumed: no frame produced for $label" >&2; return 1; }
    cp "$d/md.trr" "$WORK/$label.trr"
}

fails=0
run_case() {   # label dat trr expect|REFUSE
    local label="$1" dat="$2" trr="$3" expect="$4"
    local d="$WORK/$label"; rm -rf "$d"; mkdir -p "$d"
    ( cd "$d" && "$PLUMED" driver --plumed "$dat" --itrr "$trr" ) > "$WORK/$label.log" 2>&1
    local rc=$? energy refused=no
    grep -a -q "ERROR in input to action KENREF\|KEnRef refuses to continue" "$WORK/$label.log" && refused=yes
    energy=$(grep -a -m1 -E "Step 0  KEnRef energy" "$WORK/$label.log" | awk '{print $NF}')

    printf '%-28s rc=%-4s ' "$label" "$rc"
    if [ "$expect" = REFUSE ]; then
        if [ "$refused" = yes ] && [ "$rc" -ne 0 ]; then
            echo "REFUSED as required"
            grep -a -m1 -o 'cannot be used to repair periodic images:.*' "$WORK/$label.log" \
                | cut -c1-160 | sed 's/^/      | /'
        else
            echo "*** FAIL: expected a refusal, got rc=$rc energy=${energy:-none}"; fails=$((fails+1))
        fi
    else
        if [ "$rc" -eq 0 ] && [ "$energy" = "$expect" ]; then
            echo "E0=$energy  (GROMACS gives $expect)  OK"
        else
            echo "*** FAIL: E0=${energy:-none} expected $expect"; fails=$((fails+1))
        fi
    fi
}

# ---- inputs -----------------------------------------------------------------
cat > "$WORK/gb3.toml" <<EOF
model = SIGMA
guide = "guideC-alpha"
index = "$GB3_SET/KEnRefAtomIndex.ndx"
exp-data-folder = "$GB3_SET"
proton-mhz = 700.0
atomname-mapping = "$GB3_SET/GB3_27_10us.pdb"
ref = "$GB3_SET/GB3_27_10us.pdb"
EOF
cat > "$WORK/ubq.toml" <<EOF
model = PLATEAUS
guide = "guideC-alpha"
index = "$UBQ_SET/KEnRefAtomIndex.ndx"
exp-data-file = "$UBQ_SET/singleton_data_10nsstart+fit_3-5_1977pairs_80_A.csv"
atomname-mapping = "$UBQ_SET/00ns.pdb"
ref = "$UBQ_SET/00ns.pdb"
EOF

# The A and C inputs differ in ONE thing: which reference PLUMED is given.
for kind in whole broken; do
    case $kind in
        whole)  ref="$GB3_SET/GB3_27_10us.pdb" ;;
        broken) ref="$GB3_SET/GB3_27_10us_PBC-BROKEN.pdb" ;;
    esac
    cat > "$WORK/gb3_$kind.dat" <<EOF
KENREF ...
  LABEL=kenref
  MODEL=SIGMA
  K=1.0
  N=0.25
  EXP_DATA_FOLDER=$GB3_SET/
  PROTON_MHZ=700.0
  GUIDE_ATOMS=$(guide_group "$GB3_SET/KEnRefAtomIndex.ndx")
  REF=$ref
  ATOMNAME_MAPPING=$ref
  MAX_FORCE=999
  FIT_TO_REFERENCE
... KENREF
EOF
done
cat > "$WORK/ubq.dat" <<EOF
KENREF ...
  LABEL=kenref
  MODEL=PLATEAUS
  K=1.0
  N=0.25
  EXP_DATA_FILE=$UBQ_SET/singleton_data_10nsstart+fit_3-5_1977pairs_80_A.csv
  GUIDE_ATOMS=$(guide_group "$UBQ_SET/KEnRefAtomIndex.ndx")
  REF=$UBQ_SET/00ns.pdb
  ATOMNAME_MAPPING=$UBQ_SET/00ns.pdb
  MAX_FORCE=999
  FIT_TO_REFERENCE
... KENREF
EOF

echo "KEnRef PBC wholeness gate -- PLUMED side"
echo "  plumed   $PLUMED"
echo "  KEnRef   $BIN   (used only to emit the step-0 frames)"
echo "  workdir  $WORK"
echo
make_frame gb3 "$GB3_TPR" || exit 2
make_frame ubq "$UBQ_TPR" || exit 2

echo "A/B) PLUMED must repair the torn fixture, agree with GROMACS, and leave the whole one alone:"
run_case gb3_broken_repaired  "$WORK/gb3_whole.dat"  "$WORK/gb3.trr" "$EXPECT_GB3"
run_case ubiquitin_whole      "$WORK/ubq.dat"        "$WORK/ubq.trr" "$EXPECT_UBQ"
echo
echo "C) with the reference itself split, PLUMED cannot repair -- and must say so:"
run_case gb3_broken_reference "$WORK/gb3_broken.dat" "$WORK/gb3.trr" REFUSE

echo
if [ "$fails" -eq 0 ]; then
    echo "rt_validate_plumed: PASS (3/3)"
else
    echo "rt_validate_plumed: FAIL ($fails of 3)"
fi
exit $((fails == 0 ? 0 : 1))
