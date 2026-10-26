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
# So this gate never validates on the broken fixture alone. It runs FOUR cases,
# and the contrast between them is the whole point:
#
#   A  broken GB3, repair on   -> the bond-connectivity make-whole fixes it, PASS
#   B  whole ubiquitin         -> the repair must be a no-op here, PASS
#   C  broken GB3, repair off  -> KENREF_NO_MAKEWHOLE=1, the split-check REFUSES
#   D  broken GB3, no repair   -> its ENERGY must fail A's comparison
#
# A alone proves nothing (it would also "pass" if the repair were a no-op and the
# energy simply wrong). B is what proves the repair does not damage good input.
# C is what proves the check can still see a tear when the repair is not there.
# D is what proves A's tolerance (below) cannot accept the torn protein's energy.
#
# Both fixtures are committed in-repo, so this is self-contained.
#
# Usage:  ./rt_validate.sh
#         BIN=/path/to/KEnRef MPIRUN=/path/to/mpirun NSTEPS=2 ./rt_validate.sh
# D needs energycalc, taken from beside $BIN (set ENERGYCALC= otherwise); without
# it D is reported as SKIPPED, never as passed.
# =============================================================================
set -u

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN="${BIN:-$REPO/cmake-build-cli-gmx/KEnRef}"
MPIRUN="${MPIRUN:-mpirun}"
ENERGYCALC="${ENERGYCALC:-$(dirname "$BIN")/energycalc}"
NSTEPS="${NSTEPS:-2}"
WORK="${WORK:-${TMPDIR:-/tmp}/kenref-rt-validate.$$}"

export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"

# Expected step-0 energies, PER SIMD WIDTH.
#
# The previous comment here claimed these were "properties of the committed fixtures, not
# of the machine: bit-identical across 1, 2 and 4 ranks and across repeats". That was
# measured on AVX_512 and is NOT true in general. A narrower vector changes the order the
# energy terms are summed, and the 6th significant figure moves. Measured 2026-10-14 on the
# v2.0.3 tier (kenref 087d753, gromacs 2026.4), and constant across all four GROMACS
# releases within a width -- which is what establishes that the SIMD width, not GROMACS, is
# the variable:
#
#   fixture      AVX_512       AVX2_256      AVX_256
#   gb3_broken   3.65155e-05   3.6515e-05    3.65161e-05
#   ubiquitin    7.18869e-08   7.18869e-08   7.18869e-08     (width-insensitive so far)
#
# AVX_256 CAVEAT: there, and only there, gb3_broken's step-0 energy is 3.65161e-05 at np=1
# but 3.6516e-05 at np=2 and np=4 -- a reproducible 2.7e-06 relative serial-vs-DD gap (four
# independent confirmations, deterministic within each configuration). np=2 and np=4 agree
# with EACH OTHER exactly, so the decomposition axis is self-consistent and it is the serial
# path that is the outlier. The value below is the np=1 one, because this script runs serial.
# dd_validate.sh, which sweeps rank counts, therefore reports three genuinely-failing cells
# for gb3 on AVX_256. That red is the check working, not a stale number -- see dd_validate.sh.
#
# DERIVED FROM $BIN'S OWN PATH, not taken as a knob. The path says which build is actually
# under test; a knob can disagree with the binary it claims to describe. A dev-tree binary has
# no width in its path, so ACCEL= must be set explicitly there.
# NO FALLBACK to AVX_512's numbers for an unrecognised width: silently reusing them is exactly
# how 32 of 60 smoke checks read as a product failure on 2026-10-14 when the references were
# the only thing wrong.
# Both tier LAYOUTS, because the width sits in a different place in each: staging spells it
# .../release/AVX_512/..., the deployment tree spells it 2.0.3+2026.4_AVX_512 (and the plumed
# prefix 2.10.1-830-g684c3c255_AVX_512-kenref2.0.3+2026.4). Matching only the staging spelling
# made every deployment-tier run exit 2. AVX2_256 cannot be mistaken for AVX_256: the literal
# "AVX_256" does not occur inside "AVX2_256".
kn_width_of() {   # path -> AVX_512|AVX2_256|AVX_256, or "" when the path names no width
    case "$1" in
        *AVX_512*)  echo AVX_512  ;;
        *AVX2_256*) echo AVX2_256 ;;
        *AVX_256*)  echo AVX_256  ;;
        *)          echo ""       ;;
    esac
}
KN_WIDTH="$(kn_width_of "$BIN")"
KN_WIDTH="${KN_WIDTH:-${ACCEL:-}}"
case "$KN_WIDTH" in
    AVX_512)  EXPECT_GB3="3.65155e-05"; EXPECT_UBQ="7.18869e-08" ;;
    AVX2_256) EXPECT_GB3="3.6515e-05";  EXPECT_UBQ="7.18869e-08" ;;
    AVX_256)  EXPECT_GB3="3.65161e-05"; EXPECT_UBQ="7.18869e-08" ;;
    *) echo "rt_validate: no reference energies for this build: $BIN${ACCEL:+ (ACCEL=$ACCEL)}" >&2
       echo "  The expected energies differ per width, so there is nothing safe to compare against." >&2
       echo "  Its SIMD width is not one of AVX_512|AVX2_256|AVX_256, or no path names a width." >&2
       echo "  Set ACCEL to one of those -- and if it IS a width with no recorded" >&2
       echo "  reference, MEASURE it rather than borrowing another width's number." >&2
       exit 2 ;;
esac

# HOW CLOSE IS "EQUAL". Every energy here is compared as an exact string of std::cout's six
# significant figures, EXCEPT ONE: the repaired-GB3 SIGMA energy (case A), which must be within
# a relative 1e-05 of the table above.
#
# Why that one. It is a near-zero residual -- 3.65e-05, left over from terms of order 1 -- so
# its 6th figure is the first place summation order shows, and the order depends on more than
# the SIMD width the table is keyed on. One step in that printed figure is 2.7e-06 relative.
# Measured: two CI runs of the IDENTICAL configuration (GitHub ubuntu-24.04 runner, clang,
# AVX2_256, GROMACS 2026.4, same commit) gave 3.6515e-05 and 3.65151e-05, 2.7e-06 apart,
# in this script and in rt_validate_plumed.sh alike. Same code, same input: the runner's
# hardware. An exact string made that a coin toss.
#
# Why 1e-05 cannot let the failure this gate exists for through. With no repair -- the torn
# protein -- the same frame gives 0.00125246, 34x the repaired value: a relative difference of
# 33, six orders of magnitude outside the tolerance. Case D measures that on every run rather
# than trusting this comment. The smallest effect anyone has measured on this energy, PLUMED's
# put-in-box rounding under domain decomposition (~3e-08 relative), is well inside it.
#
# What it does NOT replace. The widths' values differ from one another by 1.4e-05 (AVX_512 vs
# AVX2_256) to 3.0e-05 (AVX2_256 vs AVX_256) relative -- outside the tolerance, but not by
# much, so the per-width table above still decides what is expected. Ubiquitin's energy is not
# a residual and stays exact. Comparisons BETWEEN two engines or two rank counts on the same
# machine stay exact too (rt_validate_plumed.sh, dd_validate.sh): those are not hardware noise.
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

run_case() {   # label tpr toml expect|REFUSE rtol|exact [env=val ...]
    local label="$1" tpr="$2" toml="$3" expect="$4" rtol="$5"; shift 5
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
        local rel ok
        rel="$(kn_close "$energy" "$expect" "$rtol")"; ok=$?
        if [ "$rc" -eq 0 ] && [ "$ok" -eq 0 ]; then
            echo "E0=$energy  expected $expect  rel $rel  ($rtol)  OK"
        else
            echo "*** FAIL: E0=${energy:-none}  expected $expect  rel $rel  ($rtol)"
            fails=$((fails + 1))
        fi
    fi
}

echo "KEnRef PBC wholeness gate"
echo "  binary   $BIN"
echo "  workdir  $WORK"
echo
# D: the torn protein's energy, put through case A's comparison, must FAIL it. The engine never
# reports that energy -- it repairs, and with the repair off the split-check refuses (case C) --
# so it is computed by energycalc, which does neither, from the step-0 frame exactly as GROMACS
# holds it: KEnRef repairs its own copy of the coordinates, never the state it was handed.
skipped=0
torn_case() {
    local label=gb3_torn_energy d="$WORK/gb3_torn"
    printf '%-26s ' "$label"
    if [ ! -x "$ENERGYCALC" ]; then
        echo "SKIPPED: no energycalc at $ENERGYCALC (set ENERGYCALC=)"
        skipped=1; return
    fi
    rm -rf "$d"; mkdir -p "$d"; cp "$GB3_TPR" "$d/md.tpr"
    ( cd "$d" && "$MPIRUN" --oversubscribe -np 1 "$BIN" --k 1.0 --n 0.25 --max-force 999 \
          --params "$WORK/gb3.toml" -- -nsteps 0 -ntomp 1 -deffnm md ) > "$d/frame.log" 2>&1
    ( cd "$d" && "$ENERGYCALC" -m SIGMA -g guideC-alpha -d "$GB3_SET/KEnRefAtomIndex.ndx" \
          -r "$GB3_SET/GB3_27_10us.pdb" -x "$GB3_SET" -z 700 -k 1.0 -n 0.25 -i md.trr -o torn.out ) \
          > "$d/energycalc.log" 2>&1
    local torn rel ok
    torn="$(awk 'NR == 1 {print $2}' "$d/torn.out" 2>/dev/null)"
    rel="$(kn_close "$torn" "$EXPECT_GB3" "$KN_GB3_RTOL")"; ok=$?
    if [ "$rel" = n/a ]; then
        echo "*** FAIL: no torn energy computed, see $d/energycalc.log"; fails=$((fails + 1))
    elif [ "$ok" -ne 0 ]; then
        echo "E=$torn  vs $EXPECT_GB3  rel $rel > $KN_GB3_RTOL  REJECTED as required"
    else
        echo "*** FAIL: torn E=$torn is within $KN_GB3_RTOL of $EXPECT_GB3: case A cannot tell"
        echo "                           repaired from torn (see $d)"; fails=$((fails + 1))
    fi
}

echo "A/B) the repair must fix the torn fixture AND leave the whole one alone:"
run_case gb3_broken_repaired "$GB3_TPR" "$WORK/gb3.toml" "$EXPECT_GB3" "$KN_GB3_RTOL"
run_case ubiquitin_whole     "$UBQ_TPR" "$WORK/ubq.toml" "$EXPECT_UBQ" exact
echo
echo "C) repair disabled -- the split-check must still catch the tear:"
run_case gb3_repair_disabled "$GB3_TPR" "$WORK/gb3.toml" REFUSE exact KENREF_NO_MAKEWHOLE=1
echo
echo "D) the torn protein's energy must fail A's tolerance:"
torn_case

echo
total=$((4 - skipped))
if [ "$fails" -eq 0 ]; then
    if [ "$skipped" -eq 1 ]; then
        echo "rt_validate: PASS ($total/$total, D SKIPPED)"
    else
        echo "rt_validate: PASS ($total/$total)"
    fi
else
    echo "rt_validate: FAIL ($fails of $total)"
fi
exit $((fails == 0 ? 0 : 1))
