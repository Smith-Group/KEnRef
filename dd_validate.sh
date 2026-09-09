#!/usr/bin/env bash
# =============================================================================
# KEnRef runtime gate: domain decomposition
#
# WHAT THIS IS ASSERTING. Under DD each rank holds only its own domain, so KEnRef
# locates its atoms through gmx::LocalAtomSet, has every rank fill ONLY the rows it
# owns (zeroes elsewhere), and sums the buffers over the simulation's ranks. That
# reconstructs the set exactly -- 0 + x == x in IEEE-754, and addition is
# commutative -- but ONLY if every row has exactly one writer. Two ranks claiming a
# row would double it; none would leave it at the origin. Neither crashes; both
# merely bias the restraint. KENREF_DD_SELFCHECK=1 turns that property into a
# runtime assertion (it is silent on success and gmx_fatal on failure), so a run
# that EXITS 0 with the self-check on is the self-check passing.
#
# The observable consequence is rank-count independence: the step-0 KEnRef energy
# must be bit-identical at 1, 2 and 4 ranks, and with a separate PME rank
# (mpi_comm_mygroup excludes PME-only ranks, so no reduction involves them).
#
# NOT asserted: anything downstream of a step of MD. This harness is not bitwise
# reproducible run-to-run -- repeating an IDENTICAL command gives a different
# md.trr/md.gro/md.edr a good fraction of the time. The step-N KEnRef energy
# inherits that drift, and measurably so: over 5 repeats the runs whose md.trr
# matched had one step-10 energy and the runs whose md.trr differed had another,
# in exact lockstep. That is GROMACS's nondeterminism arriving through the
# coordinates, not KEnRef's -- KEnRef is a deterministic function of the positions
# it is handed, which is why the STEP-0 energy is 5/5 everywhere: step 0 is
# evaluated on the tpr's own coordinates, before any MD has run.
#
# So the pass condition is the step-0 energy alone. The step-N energy and the
# trajectory files are REPORTED as match counts over N repeats -- they measure the
# harness, and a low count there is expected, not a regression.
#
# Both fixtures are run, torn and whole -- never the broken one alone. See
# res/PBC-BROKEN.md and rt_validate.sh.
#
# WHY THE LENGTH IS NOT NEGOTIABLE -- this harness has already earned it once. Run at 10
# steps, everything below passes. Run at 500, and on 2026-08-23 the GB3 rows aborted at
# ~step 120, identically at 1 rank and at 4: the split-check of the day tested SIZE (extent
# from the centroid vs half the box) rather than TEARING, and whole GB3 sits 0.004% on the
# wrong side of that line -- 21.657 against a 21.6561 half-box -- so MD breathing pushed a
# perfectly whole molecule over it. The check was rewritten to compare the delivered extent
# against the smallest extent any choice of periodic images could give, which can only fire
# on an interior void, and google_tests/testPeriodicSplitCheck.cpp pins that as a regression
# test. Everything now passes at 500. Do not lower the default back to save wall-clock.
#
# Usage:  ./dd_validate.sh
#         BIN=/path/to/KEnRef MPIRUN=/path/to/mpirun N=5 NSTEPS=10 ./dd_validate.sh
# =============================================================================
set -u

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN="${BIN:-$REPO/cmake-build-cli-gmx/KEnRef}"
MPIRUN="${MPIRUN:-mpirun}"
# 500, NOT a handful. nstlist is 40 (GROMACS raises it to 50), so a short run only ever
# uses the INITIAL partition -- which is precisely the case that already worked before any
# of the DD work. 500 steps gives ~10 repartitions; grep md.log for "DD  step" to confirm.
NSTEPS="${NSTEPS:-500}"
N="${N:-5}"                                       # repeats, for the match counts
WORK="${WORK:-${TMPDIR:-/tmp}/kenref-dd-validate.$$}"

export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"

EXPECT_GB3="3.65155e-05"
EXPECT_UBQ="7.18869e-08"

GB3_SET="$REPO/res/sigma/md_single"
GB3_TPR="$GB3_SET/md-00_PBC-BROKEN.tpr"
UBQ_SET="$REPO/res/10nsstart+fitting"
UBQ_TPR="$UBQ_SET/t000.00.tpr"

[ -x "$BIN" ] || { echo "dd_validate: no KEnRef binary at $BIN (set BIN=)" >&2; exit 2; }
command -v "$MPIRUN" >/dev/null || { echo "dd_validate: no mpirun (set MPIRUN=)" >&2; exit 2; }

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

# Run once. Echoes "<rc> <E0> <E_last>"; the self-check is on, so rc!=0 includes it firing.
one_run() {   # dir tpr toml np extra-mdrun-args...
    local d="$1" tpr="$2" toml="$3" np="$4"; shift 4
    rm -rf "$d"; mkdir -p "$d"; cp "$tpr" "$d/md.tpr"
    ( cd "$d" && KENREF_DD_SELFCHECK=1 "$MPIRUN" --oversubscribe -np "$np" "$BIN" \
          --k 1.0 --n 0.25 --max-force 999 --params "$toml" \
          -- -nsteps "$NSTEPS" -ntomp 1 "$@" -deffnm md ) > "$d/run.log" 2>&1
    local rc=$?
    echo "$rc" \
         "$(grep -a -m1 -E "^Step: 0 Energy:"        "$d/run.log" | awk '{print $4}')" \
         "$(grep -a -m1 -E "^Step: $NSTEPS Energy:"  "$d/run.log" | awk '{print $4}')"
}

# ---- part 1: rank-count independence -------------------------------------------------
ranks_independent() {   # label tpr toml expect
    local label="$1" tpr="$2" toml="$3" expect="$4"
    echo "  $label  (expected step-0 energy $expect)"
    local cfg np extra
    for cfg in "1:" "2:-npme 0" "4:-npme 0" "4:-npme 1"; do
        np="${cfg%%:*}"; extra="${cfg#*:}"
        local tag="${label}_np${np}${extra:+_${extra// /}}"
        # shellcheck disable=SC2086
        read -r rc e0 elast <<<"$(one_run "$WORK/ranks/$tag" "$tpr" "$toml" "$np" $extra)"
        local grid
        grid=$(grep -a -m1 "Domain decomposition grid" "$WORK/ranks/$tag/md.log" 2>/dev/null \
               | sed 's/^ *//')
        printf '    %-24s rc=%-3s %-34s E0=%-13s E%-4s=%-13s ' \
               "np=$np ${extra:-(serial)}" "$rc" "${grid:-no DD (single rank)}" "$e0" "$NSTEPS" "$elast"
        if [ "$rc" -eq 0 ] && [ "$e0" = "$expect" ]; then
            echo OK
        else
            echo "*** FAIL"
            fails=$((fails + 1))
        fi
    done
}

# ---- part 2: repeat the SAME configuration, count what actually matches ---------------
repeatability() {   # label tpr toml np extra...
    local label="$1" tpr="$2" toml="$3" np="$4"; shift 4
    local base="$WORK/repro/$label" i rc e0 elast
    local e0_ref="" el_ref="" e0_same=0 el_same=0 trr=0 gro=0 edr=0
    rm -rf "$base"; mkdir -p "$base"
    for i in $(seq 1 "$N"); do
        read -r rc e0 elast <<<"$(one_run "$base/r$i" "$tpr" "$toml" "$np" "$@")"
        [ "$rc" -eq 0 ] || { fails=$((fails + 1)); }
        if [ "$i" -eq 1 ]; then
            e0_ref="$e0"; el_ref="$elast"
        else
            cmp -s "$base/r1/md.trr" "$base/r$i/md.trr" && trr=$((trr + 1))
            cmp -s "$base/r1/md.gro" "$base/r$i/md.gro" && gro=$((gro + 1))
            cmp -s "$base/r1/md.edr" "$base/r$i/md.edr" && edr=$((edr + 1))
        fi
        [ "$e0"    = "$e0_ref" ] && e0_same=$((e0_same + 1))
        [ "$elast" = "$el_ref" ] && el_same=$((el_same + 1))
    done
    printf '    %-22s np=%-2s  E0 %d/%d  E%s %d/%d  |  md.trr %d/%d  md.gro %d/%d  md.edr %d/%d\n' \
           "$label" "$np" "$e0_same" "$N" "$NSTEPS" "$el_same" "$N" \
           "$trr" "$((N - 1))" "$gro" "$((N - 1))" "$edr" "$((N - 1))"
    # Only the STEP-0 energy is a pass condition; see the header. E$NSTEPS and the
    # trajectory files ride on GROMACS's run-to-run drift and are information only.
    #
    # DO NOT "improve" this by also asserting E$NSTEPS. Measured on GB3 at 500 steps, the
    # step-500 energy reads 4.22 / 4.26 / 4.28 / 4.23 at 1, 2, 4 and 4-with-PME ranks --
    # while the step-0 energy is BIT-IDENTICAL at every one of those rank counts and 5/5
    # stable across repeats. Both facts come from the same runs. The divergence is
    # GROMACS's summation order, amplified by 500 steps of MD; it is not KEnRef moving,
    # and an E$NSTEPS assertion would fail for a reason that has nothing to do with the
    # code under test.
    if [ "$e0_same" -ne "$N" ]; then
        echo "      *** FAIL: the step-0 KEnRef energy is not reproducible at a fixed rank count"
        fails=$((fails + 1))
    fi
}

echo "KEnRef domain-decomposition gate"
echo "  binary   $BIN"
echo "  workdir  $WORK"
echo "  repeats  $N   steps $NSTEPS   self-check ON"
echo
echo "1) rank-count independence -- the step-0 energy must not move:"
ranks_independent gb3_broken "$GB3_TPR" "$WORK/gb3.toml" "$EXPECT_GB3"
ranks_independent ubiquitin  "$UBQ_TPR" "$WORK/ubq.toml" "$EXPECT_UBQ"
echo
echo "2) repeatability at a FIXED rank count (E0 asserts; the rest is reported):"
repeatability ubiquitin_serial "$UBQ_TPR" "$WORK/ubq.toml" 1
repeatability ubiquitin_dd4    "$UBQ_TPR" "$WORK/ubq.toml" 4 -npme 0
repeatability gb3_serial       "$GB3_TPR" "$WORK/gb3.toml" 1
repeatability gb3_dd4          "$GB3_TPR" "$WORK/gb3.toml" 4 -npme 0

echo
if [ "$fails" -eq 0 ]; then
    echo "dd_validate: PASS"
else
    echo "dd_validate: FAIL ($fails checks)"
fi
exit $((fails == 0 ? 0 : 1))
