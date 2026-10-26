#!/usr/bin/env bash
# =============================================================================
# KEnRef runtime gate: periodic wholeness, PLUMED side (and engine agreement)
#
# WHY THIS EXISTS. rt_validate.sh pins the GROMACS engine. This is its PLUMED
# counterpart, and it exists because the two engines repair periodic breaks by
# DIFFERENT mechanisms and must still agree to the last digit:
#
#   GROMACS  walks the topology's BOND graph            (KEnRefMoleculeGraph)
#   PLUMED   walks a spanning tree over the REFERENCE
#
# WHICH PLUMED-SIDE REPAIR THIS TESTS (2026-10-19): KEnRef's OWN
# makeRequestedAtomsWhole(), and nothing else. REQUIRES PLUMED >= 1e03000ea
# (plumed/plumed2#1515), which added the call in getLocalModelX() and switched
# both readers -- getGuideAtomsX and fillSubAtomsX -- from getPositions() to the
# repaired wholePositions_. The generated inputs below therefore carry NO
# MOLINFO and NO WHOLEMOLECULES: PLUMED's own make-whole would repair the frame
# first and this gate would measure PLUMED's repair while claiming to measure
# ours. That is exactly what it did between 2026-10-13 and this commit, which is
# why those two lines were in it and why they are gone now.
#
# AGAINST AN OLDER PLUMED THIS GATE REFUSES TO RUN, by checking the kernel
# rather than a version string. No version string can tell the two apart:
# `plumed info --long-version` is 2.11.0-dev on both sides of #1515, and
# --git-version reports a describe token that cannot be ordered without a
# plumed2 clone and says "-dirty" for any locally modified build. Measured:
#   684c3c255 (pre)   makeRequestedAtomsWhole defined, 0 call sites
#   1e03000ea (post)  defined, 1 call site
# so the predicate is "the kernel defines the repair and calls it", which is the
# condition that makes these assertions mean what they say. It is not a proxy.
#
# WHY A REFUSAL AND NOT A FALLBACK: post-#1515 a regression in OUR repair also
# ends at the split-check backstop with the same refusal text an old kernel
# produces. So "this PLUMED predates #1515" and "our repair broke" cannot be
# told apart from behaviour, and a gate that adapted between them would file the
# second as the first -- agreement mistaken for correctness, rebuilt inside the
# check meant to catch it. Refuse, name the kernel, and let a human re-pin.
#
# AND NOTE WHY THIS MATTERED: before 2026-10-13 this script asserted the same
# three energies with no MOLINFO/WHOLEMOLECULES lines, and passed -- but only
# against a local PLUMED build carrying the unshipped patch. From committed
# sources it could not pass. Same trap as the never-committed KENREF_DD_FAULT
# hook: a gate whose green lived outside the repository. It is in PLUMED now.
#
# PLUMED is never handed a topology, so the reference structure is all it has --
# for its own repair exactly as for ours. That makes the reference load-bearing
# in a way it is not for GROMACS, and it is why res/*/GB3_27_10us.pdb had to be
# made whole before this could work at all -- see res/PBC-BROKEN.md.
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
# A and B each assert TWO things: PLUMED's energy against the per-width table
# below, and PLUMED's energy against GROMACS's own step-0 energy for the same
# frame, EXACTLY. GROMACS's comes from the zero-step run that writes the frame
# (see NOTE ON THE FRAME), so the two engines are compared on the same machine
# in the same job, and the agreement check does not rest on the table at all.
#
# C ASSERTS THE MECHANISM, not merely that the run died. A crash and a refusal
# both abort the driver -- both give rc=134 -- so matching on exit status alone
# lets a segfault pass as a refusal. C therefore requires one of the two
# refusal TEXTS, and prints which one fired:
#
#   * our reference guard, "cannot be used to repair periodic images" -- only
#     reachable once makeRequestedAtomsWhole() is actually called;
#   * the KEnRefDriver split-check backstop, "the 'sub' atoms are SPLIT" --
#     what fires today, when PLUMED's own makeWhole is handed a torn reference.
#
# Measured 2026-10-13: a torn reference leaves the set WORSE than no repair at
# all (span 67.7628 with MOLINFO against 64.6615 without), which is the
# reference being load-bearing, visibly, for PLUMED's repair too.
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
# are the SAME numbers here is the point: two engines, two repair mechanisms. Since
# the table check for gb3 became a tolerance (below), "the same" is asserted
# directly instead, against GROMACS's energy for the same frame -- see the header.
# Expected step-0 energies, PER SIMD WIDTH. See the long note in rt_validate.sh: a narrower
# vector changes the order the energy terms are summed and the 6th significant figure moves,
# so one set of literals cannot serve all three widths. Measured 2026-10-14 on the v2.0.3
# tier (kenref 087d753, gromacs 2026.4):
#
#   fixture      AVX_512       AVX2_256      AVX_256
#   gb3_broken   3.65155e-05   3.6515e-05    3.65161e-05
#   ubiquitin    7.18869e-08   7.18869e-08   7.18869e-08
#
# This script runs serial, so AVX_256's np=1 value is the right one here. (On AVX_256 the DD
# runs differ from serial by 2.7e-06; that is dd_validate.sh's business, documented there.)
#
# Derived from the computing binary's OWN PATH, not taken as a knob: the path says which build
# is actually under test. On this script that binary is $PLUMED, not $BIN -- see the note below
# the layout one. NO FALLBACK for an unrecognised width -- borrowing AVX_512's numbers is how 32
# of 60 smoke checks read as a product failure on 2026-10-14 when only the references were wrong.
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

# WHICH binary's width, though? On this script the asserted energy is computed by $PLUMED --
# by the libkenref_core that the plumed kernel links -- while $BIN only runs the tpr for zero
# steps to emit a frame, and a zero-step frame is just the tpr's own coordinates, so $BIN's
# width cannot reach PLUMED's number. ($BIN's width DOES reach the GROMACS energy that the same
# zero-step run prints, which A and B compare against PLUMED's exactly; a BIN of another width
# therefore fails there, as engines disagreeing, which is what that run would be.)
# Keying on $BIN was therefore keying on the wrong build: BIN from an AVX_512 tier with PLUMED
# from an AVX_256 one would compare AVX_512's references against AVX_256's energies. $BIN stays
# as the fallback for a plumed whose path names no width, and when both name one and they
# DISAGREE there is no set of references that describes the run, so stop.
# $PLUMED is resolved through the PATH first, so the usual `module load plumed-dev/...` form --
# where $PLUMED is the bare word "plumed" -- still yields a path with the width in it.
KN_PLUMED_PATH="$(command -v "$PLUMED" 2>/dev/null || echo "$PLUMED")"
KN_W_PLUMED="$(kn_width_of "$KN_PLUMED_PATH")"
KN_W_BIN="$(kn_width_of "$BIN")"
if [ -n "$KN_W_PLUMED" ] && [ -n "$KN_W_BIN" ] && [ "$KN_W_PLUMED" != "$KN_W_BIN" ]; then
    echo "rt_validate_plumed: this run mixes SIMD widths, so no reference fits it." >&2
    echo "  plumed ($KN_W_PLUMED): $KN_PLUMED_PATH" >&2
    echo "  KEnRef ($KN_W_BIN): $BIN" >&2
    echo "  The energy asserted here is PLUMED's; run both from the same width." >&2
    exit 2
fi
KN_WIDTH="${KN_W_PLUMED:-${KN_W_BIN:-${ACCEL:-}}}"
case "$KN_WIDTH" in
    AVX_512)  EXPECT_GB3="3.65155e-05"; EXPECT_UBQ="7.18869e-08" ;;
    AVX2_256) EXPECT_GB3="3.6515e-05";  EXPECT_UBQ="7.18869e-08" ;;
    AVX_256)  EXPECT_GB3="3.65161e-05"; EXPECT_UBQ="7.18869e-08" ;;
    *) echo "$(basename "$0"): no reference energies for this build: $KN_PLUMED_PATH${ACCEL:+ (ACCEL=$ACCEL)}" >&2
       echo "  The expected energies differ per width, so there is nothing safe to compare against." >&2
       echo "  Its SIMD width is not one of AVX_512|AVX2_256|AVX_256, or no path names a width." >&2
       echo "  Set ACCEL to one of those -- and if it IS a width with no recorded" >&2
       echo "  reference, MEASURE it rather than borrowing another width's number." >&2
       exit 2 ;;
esac

# The repaired-GB3 energy against the table: within a relative 1e-05, not an exact string. See
# the note in rt_validate.sh for the measurements: identical CI configurations gave 3.6515e-05
# and 3.65151e-05 in BOTH gates (2.7e-06 apart, the runner's hardware); the torn protein's
# 0.00125246 is 33 relative away; PLUMED's put-in-box rounding under DD is ~3e-08. Ubiquitin
# stays exact, and so does PLUMED against GROMACS: two engines on one machine is not noise.
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
GB3_TPR="$GB3_SET/md-00_PBC-BROKEN.tpr"
UBQ_SET="$REPO/res/10nsstart+fitting"
UBQ_TPR="$UBQ_SET/t000.00.tpr"

[ -x "$BIN" ] || { echo "rt_validate_plumed: no KEnRef binary at $BIN (set BIN=)" >&2; exit 2; }
command -v "$PLUMED" >/dev/null || { echo "rt_validate_plumed: no plumed (set PLUMED=)" >&2; exit 2; }
command -v "$MPIRUN" >/dev/null || { echo "rt_validate_plumed: no mpirun (set MPIRUN=)" >&2; exit 2; }
"$PLUMED" driver --help 2>&1 | grep -q -- --itrr || {
    echo "rt_validate_plumed: this plumed has no --itrr driver support" >&2; exit 2; }

# REQUIRE PLUMED >= 1e03000ea (#1515), by asking the KERNEL whether the repair is actually called.
# Rationale and the measured 0-vs-1 call-site figures are in the header. Three ways this can fail and
# all three EXIT, never continue: no kernel found, no binutils to inspect it, or the repair present
# but never called. Refusing is the point -- the assertions below are only meaningful if our repair
# runs, and no version string distinguishes the two kernels.
# CONDITION, not a fact: this reads a direct call instruction, so a kernel built with -flto could
# inline the call away and read as 0 on a CORRECT build. The staging and deployment builds are -O3
# with no -flto (checked via `plumed info --configuration`), so it is sound for them; if that ever
# changes, this predicate needs replacing rather than relaxing.
require_our_repair() {
    local kernel sym n
    kernel="${PLUMED_KERNEL:-}"
    [ -n "$kernel" ] || kernel="$(dirname "$(command -v "$PLUMED")")/../lib/libplumedKernel.so"
    [ -f "$kernel" ] || {
        echo "rt_validate_plumed: cannot find libplumedKernel.so (tried '$kernel')." >&2
        echo "  This gate must inspect the kernel to know whether KEnRef's own repair is called." >&2
        echo "  Set PLUMED_KERNEL=/path/to/libplumedKernel.so." >&2; exit 2; }
    command -v nm >/dev/null 2>&1 && command -v objdump >/dev/null 2>&1 || {
        echo "rt_validate_plumed: needs nm and objdump (binutils) to inspect $kernel." >&2
        echo "  Refusing rather than assuming: see the header." >&2; exit 2; }
    sym="$(nm -D "$kernel" 2>/dev/null | awk '/makeRequestedAtomsWhole/{print $3; exit}')"
    [ -n "$sym" ] || {
        echo "rt_validate_plumed: $kernel does not define makeRequestedAtomsWhole at all." >&2
        echo "  That is not a pre-#1515 PLUMED, it is a kernel without KEnRef's repair -- check" >&2
        echo "  that this plumed was built against KEnRef >= 2.0.0." >&2; exit 2; }
    n="$(objdump -d --no-show-raw-insn "$kernel" 2>/dev/null | grep -cE "(call|jmp).*<${sym}>")"
    [ "${n:-0}" -gt 0 ] || {
        echo "rt_validate_plumed: this gate needs PLUMED >= 1e03000ea (plumed/plumed2#1515)." >&2
        echo "  kernel      : $kernel" >&2
        # The KERNEL's own embedded git version, not `$PLUMED info --git-version`: the exe and the
        # kernel can disagree (PLUMED_KERNEL overrides which kernel is loaded), and reporting the
        # exe's version beside the kernel's path is how a diagnostic ends up describing the wrong
        # artefact. Cosmetic only -- the refusal above is decided by the call-site count.
        echo "  kernel built from : $(strings -a "$kernel" 2>/dev/null \
            | grep -oE '^[0-9a-f]{9}(-dirty)?$' | head -1 || echo unknown)" >&2
        echo "  plumed exe        : $("$PLUMED" info --git-version 2>/dev/null || echo unknown)" >&2
        echo "  It DEFINES makeRequestedAtomsWhole but never calls it, so the per-step readers are" >&2
        echo "  still reading getPositions() and the torn fixture would reach the numerics torn." >&2
        echo "  Re-pin PLUMED past #1515 and rebuild. Do NOT add MOLINFO/WHOLEMOLECULES to make" >&2
        echo "  this pass: that measures PLUMED's repair while claiming to measure ours." >&2; exit 2; }
    echo "  repair   KEnRef's own makeRequestedAtomsWhole, $n call site(s) in the kernel"
}
require_our_repair

mkdir -p "$WORK"

# Zero-based-to-1-based is not our problem here: read the group straight out of the
# .ndx so PLUMED's GUIDE_ATOMS can never drift from what the KEnRef side uses.
guide_group() {   # ndx-file
    awk '/^\[ guideC-alpha \]/{f=1;next} /^\[/{f=0} f' "$1" \
        | tr -s ' \n' '\n' | grep -v '^$' | paste -sd,
}

# WHOLEMOLECULES needs the molecule's atom range. Derive it from the reference
# PDB -- the last atom that is not solvent or a counter-ion -- rather than
# hardcoding a number that would rot the moment a fixture changed. The protein
# is a contiguous prefix in every fixture here (GB3's reference is protein-only,
# 862 atoms; ubiquitin's 00ns.pdb is 1-1231 protein followed by 24720 SOL
# atoms), and this REFUSES if that stops being true instead of silently chaining
# water into the entity, which would make WHOLEMOLECULES both wrong and slow.
protein_entity() {   # reference-pdb -> "1-N"
    awk '
        /^ATOM|^HETATM/ {
            r = substr($0, 18, 4); gsub(/ /, "", r); serial = substr($0, 7, 5) + 0
            if (r == "SOL" || r == "HOH" || r == "WAT" || r == "TIP3" ||
                r == "NA" || r == "CL" || r == "NA+" || r == "CL-" ||
                r == "K" || r == "MG" || r == "ZN" || r == "CA") {
                if (first_solvent == 0) first_solvent = serial
            } else last_protein = serial
        }
        END {
            if (last_protein == 0) {
                print "protein_entity: no non-solvent atoms in the reference" > "/dev/stderr"; exit 1
            }
            if (first_solvent != 0 && first_solvent < last_protein) {
                printf "protein_entity: solvent interleaves the molecule (first solvent %d, last protein %d) -- 1-N is not the right entity\n", \
                       first_solvent, last_protein > "/dev/stderr"
                exit 1
            }
            print "1-" last_protein
        }' "$1"
}

# Frame 0 in FULL precision: run the tpr for zero steps and keep the trr it writes. The same run
# prints GROMACS's step-0 KEnRef energy for that frame; keep it for the engine-agreement check.
make_frame() {   # label tpr -> $WORK/<label>.trr, $WORK/<label>.gmx_e0
    local label="$1" tpr="$2"
    local d="$WORK/frame_$label"
    rm -rf "$d"; mkdir -p "$d"; cp "$tpr" "$d/md.tpr"
    ( cd "$d" && "$MPIRUN" --oversubscribe -np 1 "$BIN" --k 1.0 --n 0.25 --max-force 999 \
          --params "$WORK/$label.toml" -- -nsteps 0 -ntomp 1 -deffnm md ) > "$WORK/frame_$label.log" 2>&1
    [ -s "$d/md.trr" ] || { echo "rt_validate_plumed: no frame produced for $label" >&2; return 1; }
    cp "$d/md.trr" "$WORK/$label.trr"
    grep -a -m1 -E "^Step: 0 Energy:" "$WORK/frame_$label.log" | awk '{print $4}' > "$WORK/$label.gmx_e0"
}

fails=0
run_case() {   # label dat trr expect|REFUSE [rtol|exact gromacs-energy]
    local label="$1" dat="$2" trr="$3" expect="$4" rtol="${5:-exact}" gmx="${6:-}"
    local d="$WORK/$label"; rm -rf "$d"; mkdir -p "$d"
    ( cd "$d" && "$PLUMED" driver --plumed "$dat" --itrr "$trr" ) > "$WORK/$label.log" 2>&1
    local rc=$? energy guard=no backstop=no
    # WHICH message, not just "it died": see the header. rc alone cannot tell a
    # refusal from a crash, because the driver aborts either way.
    grep -a -q 'cannot be used to repair periodic images' "$WORK/$label.log" && guard=yes
    grep -a -q "the 'sub' atoms are SPLIT"                "$WORK/$label.log" && backstop=yes
    energy=$(grep -a -m1 -E "Step 0  KEnRef energy" "$WORK/$label.log" | awk '{print $NF}')

    printf '%-28s rc=%-4s ' "$label" "$rc"
    if [ "$expect" = REFUSE ]; then
        if [ "$rc" -eq 0 ]; then
            echo "*** FAIL: expected a refusal, but it ran: E0=${energy:-none}"; fails=$((fails+1))
        elif [ "$guard" = yes ]; then
            echo "REFUSED by OUR reference guard"
            grep -a -m1 -o 'cannot be used to repair periodic images:.*' "$WORK/$label.log" \
                | cut -c1-160 | sed 's/^/      | /'
            echo "      | This is the expected refusal since #1515: the guard lives inside"
            echo "      | makeRequestedAtomsWhole(), so reaching it proves the repair RAN."
        elif [ "$backstop" = yes ]; then
            echo "*** the SPLIT-CHECK BACKSTOP fired, not the reference guard."
            echo "    Since #1515 the repair runs first, so the guard should have caught this."
            echo "    Reaching the backstop instead means the repair ran and left the set torn."
            grep -a -m1 -o "the 'sub' atoms are SPLIT.*" "$WORK/$label.log" \
                | cut -c1-160 | sed 's/^/      | /'
        else
            echo "*** FAIL: aborted with rc=$rc and NEITHER refusal message -- that is a crash,"
            echo "                                  not a refusal (see $WORK/$label.log)"; fails=$((fails+1))
        fi
    else
        local rel_t ok_t rel_g ok_g
        rel_t="$(kn_close "$energy" "$expect" "$rtol")"; ok_t=$?
        rel_g="$(kn_close "$energy" "$gmx" exact)"; ok_g=$?
        if [ "$rc" -eq 0 ] && [ "$ok_t" -eq 0 ] && [ "$ok_g" -eq 0 ]; then
            echo "E0=$energy  table $expect (rel $rel_t, $rtol)  GROMACS $gmx (exact)  OK"
        else
            echo "*** FAIL: E0=${energy:-none}  table $expect (rel $rel_t, $rtol)  GROMACS ${gmx:-none} (rel $rel_g, exact)"
            fails=$((fails+1))
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
    ent="$(protein_entity "$ref")" || exit 2
    cat > "$WORK/gb3_$kind.dat" <<EOF
# KENREF now repairs the requested atoms itself (plumed/plumed2#1515), so this input
# deliberately carries NO MOLINFO and NO WHOLEMOLECULES: what is under test is OUR
# makeRequestedAtomsWhole(), and PLUMED's own make-whole would mask it.
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
ubq_ent="$(protein_entity "$UBQ_SET/00ns.pdb")" || exit 2
cat > "$WORK/ubq.dat" <<EOF
# The repair lines belong here too, and that is the point of case B: it exists
# to prove the repair does not damage input that is already whole, so the repair
# must be ENABLED for the case to prove anything. Ubiquitin is whole, so these
# two lines must change the answer by nothing at all.
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
echo "  KEnRef   $BIN   (emits the step-0 frames, and GROMACS's energy for them)"
echo "  workdir  $WORK"
echo
make_frame gb3 "$GB3_TPR" || exit 2
make_frame ubq "$UBQ_TPR" || exit 2

echo "A/B) PLUMED must repair the torn fixture, agree with GROMACS, and leave the whole one alone:"
run_case gb3_broken_repaired  "$WORK/gb3_whole.dat"  "$WORK/gb3.trr" "$EXPECT_GB3" "$KN_GB3_RTOL" "$(cat "$WORK/gb3.gmx_e0")"
run_case ubiquitin_whole      "$WORK/ubq.dat"        "$WORK/ubq.trr" "$EXPECT_UBQ" exact          "$(cat "$WORK/ubq.gmx_e0")"
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
