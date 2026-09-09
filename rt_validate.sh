#!/usr/bin/env bash
# =============================================================================
# KEnRef runtime gate: periodic wholeness  (the whole/broken contrast)
#
# WHY THIS EXISTS. KEnRef restrains a *whole* molecule: the Kabsch superposition
# and the pair distances are global and have no cutoff. A molecule split across a
# periodic boundary therefore does not fail -- it quietly returns a different
# energy and different forces. Every pre-2026-08 serial run against the committed
# GB3 sets was refining a torn protein, and engine-vs-engine, byte-identity and
# R-ground-truth checks all passed while it did. See res/PBC-BROKEN.md.
#
# So this gate never validates on the broken fixture alone. It runs THREE cases,
# and the contrast between them is the whole point:
#
#   A  broken GB3, repair on   -> the bond-connectivity make-whole fixes it, PASS
#   B  whole ubiquitin         -> the repair must be a no-op here, PASS
#   C  broken GB3, repair off  -> KENREF_NO_MAKEWHOLE=1, the split-check REFUSES
#
# A alone proves nothing (it would also "pass" if the repair were a no-op and the
# energy simply wrong). B is what proves the repair does not damage good input.
# C is what proves the check can still see a tear when the repair is not there.
#
# Both fixtures are committed in-repo, so this is self-contained.
#
# Usage:  ./rt_validate.sh
#         BIN=/path/to/KEnRef MPIRUN=/path/to/mpirun NSTEPS=2 ./rt_validate.sh
# =============================================================================
set -u

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN="${BIN:-$REPO/cmake-build-cli-gmx/KEnRef}"
MPIRUN="${MPIRUN:-mpirun}"
NSTEPS="${NSTEPS:-2}"
WORK="${WORK:-${TMPDIR:-/tmp}/kenref-rt-validate.$$}"

export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"

# Expected step-0 energies. These are properties of the committed fixtures, not of
# the machine: they were bit-identical across 1, 2 and 4 ranks and across repeats.
EXPECT_GB3="3.65155e-05"
EXPECT_UBQ="7.18869e-08"

GB3_SET="$REPO/res/sigma/md_single"
GB3_TPR="$GB3_SET/md-00_PBC-BROKEN.tpr"          # the tpr is named for what it IS
UBQ_SET="$REPO/res/10nsstart+fitting"
UBQ_TPR="$UBQ_SET/t000.00.tpr"

[ -x "$BIN" ] || { echo "rt_validate: no KEnRef binary at $BIN (set BIN=)" >&2; exit 2; }
command -v "$MPIRUN" >/dev/null || { echo "rt_validate: no mpirun (set MPIRUN=)" >&2; exit 2; }

mkdir -p "$WORK"
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

fails=0

run_case() {   # label tpr toml expect|REFUSE [env=val ...]
    local label="$1" tpr="$2" toml="$3" expect="$4"; shift 4
    local d="$WORK/$label"
    rm -rf "$d"; mkdir -p "$d"; cp "$tpr" "$d/md.tpr"
    ( cd "$d" && env "$@" "$MPIRUN" --oversubscribe -np 1 "$BIN" \
          --k 1.0 --n 0.25 --max-force 999 --params "$toml" \
          -- -nsteps "$NSTEPS" -ntomp 1 -deffnm md ) > "$WORK/$label.log" 2>&1
    local rc=$?
    local refused="no" energy
    grep -a -q "KEnRef refuses to continue" "$WORK/$label.log" && refused="yes"
    energy=$(grep -a -m1 -E "^Step: 0 Energy:" "$WORK/$label.log" | awk '{print $4}')

    printf '%-26s rc=%-3s ' "$label" "$rc"
    if [ "$expect" = REFUSE ]; then
        # A refusal is the pass here: a non-zero exit AND the diagnostic.
        if [ "$refused" = yes ] && [ "$rc" -ne 0 ]; then
            echo "REFUSED as required"
            sed -n '/KEnRef refuses to continue/,+4p' "$WORK/$label.log" | sed 's/^/      | /'
        else
            echo "*** FAIL: expected a refusal, got rc=$rc energy=${energy:-none}"
            fails=$((fails + 1))
        fi
    else
        if [ "$rc" -eq 0 ] && [ "$energy" = "$expect" ]; then
            echo "E0=$energy  (expected $expect)  OK"
        else
            echo "*** FAIL: E0=${energy:-none} expected $expect"
            fails=$((fails + 1))
        fi
    fi
}

echo "KEnRef PBC wholeness gate"
echo "  binary   $BIN"
echo "  workdir  $WORK"
echo
echo "A/B) the repair must fix the torn fixture AND leave the whole one alone:"
run_case gb3_broken_repaired "$GB3_TPR" "$WORK/gb3.toml" "$EXPECT_GB3"
run_case ubiquitin_whole     "$UBQ_TPR" "$WORK/ubq.toml" "$EXPECT_UBQ"
echo
echo "C) repair disabled -- the split-check must still catch the tear:"
run_case gb3_repair_disabled "$GB3_TPR" "$WORK/gb3.toml" REFUSE KENREF_NO_MAKEWHOLE=1

echo
if [ "$fails" -eq 0 ]; then
    echo "rt_validate: PASS (3/3)"
else
    echo "rt_validate: FAIL ($fails of 3)"
fi
exit $((fails == 0 ? 0 : 1))
