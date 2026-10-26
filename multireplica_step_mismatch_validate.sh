#!/usr/bin/env bash
# =============================================================================
# KEnRef runtime gate: multi-replica step agreement  (-multidir)
#
# WHY THIS EXISTS. In a -multidir run KEnRef exchanges data between the replicas
# at EVERY step (a gather and a scatter on the replicas' main ranks). A replica
# with fewer steps to run finishes and exits, and the others then block forever
# in the next exchange: the run HANGS instead of failing. GROMACS does not stop
# it -- a mismatched nsteps only draws "The number of steps is not consistent
# across multi simulations, but we are proceeding anyway!". KEnRef therefore
# refuses such a run at setup (gmxinterface/KEnRefReplicaSteps.h).
#
# Every other gate copies ONE tpr into every replica, so none of them can reach
# this: the replicas here are deliberately given DIFFERENT tprs. Five cases,
# two replicas each, all from the committed res/sigma/md_double fixture:
#
#   A  nsteps 20 / 20                          -> runs to step 20
#   B  nsteps 20 / 10                          -> KEnRef REFUSES      (hung before the fix)
#   C  both restarted (-cpi) at step 10, one
#      tpr changed to end at 15 instead of 20  -> KEnRef REFUSES      (hung before the fix)
#   D  restarted from checkpoints at steps
#      10 and 5                                -> GROMACS REFUSES (initial step)
#   E  tprs 20 / 10, but mdrun -nsteps 20      -> runs to step 20
#
# A and E are the contrast: they prove the refusal is not simply "multi-replica
# is refused". C is the restart path, where the comparison has to use the LAST
# step (init-step + nsteps) because a checkpoint restart rewrites both. D is
# GROMACS's own check, which applies only when the simulations share state, and
# KEnRef's multisim subscription is what makes that true -- so D failing means
# that subscription regressed, not that GROMACS changed its mind.
#
# A refusal must come quickly, and must carry its message: a hang is caught by
# a timeout (TO, default 120 s; the refusals measured under 1 s, and the hangs
# they replace sat until killed), and a crash with no message is a FAIL, not a
# refusal. Measured 2026-10-25 against release/1.x with and without the fix, on
# GROMACS 2025.4 and 2026.3: without it, B and C were killed at the timeout.
#
# This gate asserts no energy: a replica that reaches its last step is the
# pass condition. (The fixture's tpr and reference are the torn GB3 ones on
# this branch; whether the energy is right is not what is tested here.)
#
# Usage:  GMX=/path/to/gmx_mpi ./multireplica_step_mismatch_validate.sh
#         BIN=/path/to/KEnRef MPIRUN=/path/to/mpirun GMX=... ./multireplica_step_mismatch_validate.sh
# GMX is any GROMACS gmx binary of the release KEnRef was built against; it only
# writes the per-replica tprs (convert-tpr) and the checkpoints C and D restart from.
# =============================================================================
set -u

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN="${BIN:-$REPO/cmake-build-cli-gmx/KEnRef}"
MPIRUN="${MPIRUN:-mpirun}"
GMX="${GMX:-$(command -v gmx_mpi || command -v gmx || echo gmx_mpi)}"
TO="${TO:-120}"
WORK="${WORK:-${TMPDIR:-/tmp}/kenref-multireplica-step-mismatch-validate.$$}"

export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"

SET="$REPO/res/sigma/md_double"
TPR="$SET/md-00.tpr"

[ -x "$BIN" ] || { echo "multireplica_step_mismatch_validate: no KEnRef binary at $BIN (set BIN=)" >&2; exit 2; }
command -v "$MPIRUN" >/dev/null || { echo "multireplica_step_mismatch_validate: no mpirun (set MPIRUN=)" >&2; exit 2; }
command -v "$GMX" >/dev/null || { echo "multireplica_step_mismatch_validate: no gmx binary (set GMX=)" >&2; exit 2; }
command -v timeout >/dev/null || { echo "multireplica_step_mismatch_validate: needs timeout(1) to tell a hang from a run" >&2; exit 2; }

mkdir -p "$WORK"
cat > "$WORK/gb3.toml" <<EOF
model = SIGMA
guide = "guideC-alpha"
index = "$SET/KEnRefAtomIndex.ndx"
exp-data-folder = "$SET"
proton-mhz = 700.0
atomname-mapping = "$SET/GB3_27_10us.pdb"
ref = "$SET/GB3_27_10us.pdb"
EOF

fails=0

mk_tpr() {   # dir nsteps
    mkdir -p "$1"
    "$GMX" convert-tpr -s "$TPR" -nsteps "$2" -o "$1/md.tpr" > "$1/convert.log" 2>&1 \
        || { echo "multireplica_step_mismatch_validate: convert-tpr failed, see $1/convert.log" >&2; exit 2; }
}

mk_cpt() {   # dir step -- plain GROMACS, one simulation, leaves md.cpt at <step>
    ( cd "$1" && "$GMX" mdrun -nsteps "$2" -ntomp 1 -deffnm md -cpo md.cpt ) > "$1/plain.log" 2>&1 \
        || { echo "multireplica_step_mismatch_validate: plain mdrun failed, see $1/plain.log" >&2; exit 2; }
}

run_case() {   # case expect(RUN|KENREF|GMX) extra-mdrun-args...
    local c="$1" expect="$2"; shift 2
    local log="$WORK/$c.log" t0=$SECONDS rc
    ( cd "$WORK/$c" && timeout -k 10 "$TO" "$MPIRUN" --oversubscribe -np 2 "$BIN" \
          --k 1.0 --n 0.25 --max-force 999 --params "$WORK/gb3.toml" \
          -- -multidir rep0 rep1 -ntomp 1 -deffnm md "$@" ) > "$log" 2>&1
    rc=$?
    local secs=$((SECONDS - t0))
    # GROMACS wraps fatal-error text across lines, so match on the flattened log.
    local flat; flat="$(tr -d '\000' < "$log" | tr '\n' ' ' | tr -s ' ')"
    local by=none
    case "$flat" in *"would not run the same number of steps"*) by=KENREF ;; esac
    case "$flat" in *"not consistent across multi simulations which share the state"*) by=GMX ;; esac
    # The ensemble's energy is printed every 10 steps; step 20 printed means the run reached its end.
    local last
    last="$(grep -a -h -c -E '^Step: 20 Energy:' "$log" "$WORK/$c"/rep*/md.log 2>/dev/null \
            | awk '{s += $1} END {print s + 0}')"

    printf '%-2s rc=%-3s %3ss  ' "$c" "$rc" "$secs"
    if [ "$rc" -eq 124 ] || [ "$rc" -eq 137 ]; then
        echo "*** FAIL: HUNG -- killed after ${TO}s. This is the deadlock the refusal exists to prevent."
        fails=$((fails + 1)); return
    fi
    case "$expect" in
        RUN)
            if [ "$rc" -eq 0 ] && [ "$by" = none ] && [ "$last" -gt 0 ]; then
                echo "ran to step 20 as required"
            else
                echo "*** FAIL: expected a complete run (refusal: $by, step 20 reached: $last), see $log"
                fails=$((fails + 1))
            fi ;;
        KENREF|GMX)
            if [ "$rc" -ne 0 ] && [ "$by" = "$expect" ]; then
                echo "REFUSED by $expect as required"
                if [ "$expect" = KENREF ]; then
                    grep -a -A4 'Last step of each replica\|Steps given by mdrun' "$log" | head -4 | sed 's/^/      | /'
                fi
            else
                echo "*** FAIL: expected a refusal by $expect, got rc=$rc and refusal by $by, see $log"
                fails=$((fails + 1))
            fi ;;
    esac
}

echo "KEnRef multi-replica step-agreement gate"
echo "  binary   $BIN"
echo "  gmx      $GMX"
echo "  workdir  $WORK"
echo
for c in A B C D E; do rm -rf "${WORK:?}/$c"; mkdir -p "$WORK/$c"; done
mk_tpr "$WORK/A/rep0" 20; mk_tpr "$WORK/A/rep1" 20
mk_tpr "$WORK/B/rep0" 20; mk_tpr "$WORK/B/rep1" 10
mk_tpr "$WORK/C/rep0" 20; mk_tpr "$WORK/C/rep1" 20
mk_cpt "$WORK/C/rep0" 10; mk_cpt "$WORK/C/rep1" 10
mk_tpr "$WORK/C/rep1.end15" 15; mv "$WORK/C/rep1.end15/md.tpr" "$WORK/C/rep1/md.tpr"
mk_tpr "$WORK/D/rep0" 20; mk_tpr "$WORK/D/rep1" 20
mk_cpt "$WORK/D/rep0" 10; mk_cpt "$WORK/D/rep1" 5
mk_tpr "$WORK/E/rep0" 20; mk_tpr "$WORK/E/rep1" 10

echo "Equal steps must run, and an mdrun -nsteps override makes them equal:"
run_case A RUN
run_case E RUN -nsteps 20
echo
echo "Different steps must be refused at setup, not hang:"
run_case B KENREF
run_case C KENREF -cpi md.cpt
echo
echo "A mismatched initial step is GROMACS's to refuse (needs KEnRef's multisim subscription):"
run_case D GMX -cpi md.cpt

echo
if [ "$fails" -eq 0 ]; then
    echo "multireplica_step_mismatch_validate: PASS (5/5)"
else
    echo "multireplica_step_mismatch_validate: FAIL ($fails of 5)"
fi
exit $((fails == 0 ? 0 : 1))
