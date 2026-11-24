#!/usr/bin/env bash
# KEnRef GROMACS-plugin gate: a stock `gmx mdrun` that loads kenref.so must give the energies of bin/KEnRef.
#
# A GROMACS built with GMX_MDMODULE_PLUGINS=ON loads the plugin when the tpr names it, as grompp records
# "mdmodule-plugins = kenref" from the mdp. The plugin reads its settings from the KENREF_* environment
# variables and KENREF_PARAMS instead of the KEnRef command line.
#
#   A/B/C) rt_validate.sh's single-simulation cases, through the plugin, against its reference energies.
#   M)     a two-replica -multidir run of the md_double fixture: the plugin must print the step-0 energy
#          that bin/KEnRef prints for the same replicas.
#
# rt_validate's GB3 fixture has only a tpr, so kenref_tpr_add_plugin writes the copies that name the plugin.
#
# Environment:
#   GMX      the stock gmx_mpi of the GROMACS KEnRef was built against
#   PLUGIN   kenref.so                         (default: next to BIN)
#   BIN      bin/KEnRef, for case M            (default: cmake-build-cli-gmx/KEnRef)
#   ADDPLUGIN kenref_tpr_add_plugin             (default: next to BIN)
#   ACCEL    the SIMD width, when no path names it
set -u
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN="${BIN:-$REPO/cmake-build-cli-gmx/KEnRef}"
PLUGIN="${PLUGIN:-$(dirname "$BIN")/kenref.so}"
ADDPLUGIN="${ADDPLUGIN:-$(dirname "$BIN")/kenref_tpr_add_plugin}"
GMX="${GMX:-$(command -v gmx_mpi || echo gmx_mpi)}"
MPIRUN="${MPIRUN:-mpirun}"
NSTEPS="${NSTEPS:-2}"
WORK="${WORK:-${TMPDIR:-/tmp}/kenref-plugin-validate.$$}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"

# The reference energies of rt_validate.sh, which differ per SIMD width
case "${ACCEL:-}${BIN}" in
    *AVX_512*)  EXPECT_GB3="3.65155e-05"; EXPECT_UBQ="7.18869e-08" ;;
    *AVX2_256*) EXPECT_GB3="3.6515e-05";  EXPECT_UBQ="7.18869e-08" ;;
    *AVX_256*)  EXPECT_GB3="3.65161e-05"; EXPECT_UBQ="7.18869e-08" ;;
    *) echo "plugin_validate: set ACCEL to AVX_512, AVX2_256 or AVX_256" >&2; exit 2 ;;
esac
GB3_RTOL="1e-05"

for f in "$PLUGIN" "$ADDPLUGIN" "$BIN"; do
    [ -e "$f" ] || { echo "plugin_validate: $f does not exist (set PLUGIN=, ADDPLUGIN= or BIN=)" >&2; exit 2; }
done
command -v "$GMX" >/dev/null || { echo "plugin_validate: no gmx_mpi (set GMX=)" >&2; exit 2; }
command -v "$MPIRUN" >/dev/null || { echo "plugin_validate: no mpirun (set MPIRUN=)" >&2; exit 2; }

GB3_SET="$REPO/res/sigma/md_single"
UBQ_SET="$REPO/res/10nsstart+fitting"
DOUBLE_SET="$REPO/res/sigma/md_double"
mkdir -p "$WORK/plugins"
# GROMACS loads <name>.so from GMX_MDMODULE_PLUGIN_PATH
ln -sf "$(cd "$(dirname "$PLUGIN")" && pwd)/$(basename "$PLUGIN")" "$WORK/plugins/kenref.so"
export GMX_MDMODULE_PLUGIN_PATH="$WORK/plugins"

write_toml() {   # file set model data-key data-value reference-pdb
    cat > "$1" <<EOF
model = $3
guide = "guideC-alpha"
index = "$2/KEnRefAtomIndex.ndx"
$4 = "$5"
proton-mhz = 700.0
atomname-mapping = "$6"
ref = "$6"
EOF
}
write_toml "$WORK/gb3.toml" "$GB3_SET" SIGMA exp-data-folder "$GB3_SET" "$GB3_SET/GB3_27_10us.pdb"
write_toml "$WORK/ubq.toml" "$UBQ_SET" PLATEAUS exp-data-file \
    "$UBQ_SET/singleton_data_10nsstart+fit_3-5_1977pairs_80_A.csv" "$UBQ_SET/00ns.pdb"
write_toml "$WORK/double.toml" "$DOUBLE_SET" SIGMA exp-data-folder "$DOUBLE_SET" "$DOUBLE_SET/GB3_27_10us.pdb"

# The settings that rt_validate gives bin/KEnRef on its command line
KENREF_ENV=(KENREF_K=1.0 KENREF_N=0.25 KENREF_MAX_FORCE=999)

add_plugin() {   # in.tpr out.tpr
    "$ADDPLUGIN" "$1" "$2" kenref > "$2.log" 2>&1 \
        || { echo "plugin_validate: kenref_tpr_add_plugin failed, see $2.log" >&2; exit 2; }
}

close() {   # got want rtol|exact -> prints the relative difference, fails when out of tolerance
    awk -v g="$1" -v w="$2" -v t="$3" 'BEGIN {
        num = "^[-+]?([0-9]+[.]?[0-9]*|[.][0-9]+)([eE][-+]?[0-9]+)?$"
        if (g !~ num || w !~ num || w + 0 == 0) { print "n/a"; exit 1 }
        d = (g - w) / w; if (d < 0) d = -d
        printf "%.1e\n", d
        exit (t == "exact" ? (g != w) : !(d <= t + 0)) }'
}

step0_energy() {   # log
    grep -a -m1 -E "^Step: 0 Energy:" "$1" | awk '{print $4}'
}

fails=0
run_case() {   # label tpr toml expect|REFUSE rtol|exact [env=val ...]
    local label="$1" tpr="$2" toml="$3" expect="$4" rtol="$5"; shift 5
    local d="$WORK/$label"
    rm -rf "$d"; mkdir -p "$d"; add_plugin "$tpr" "$d/md.tpr"
    ( cd "$d" && env "${KENREF_ENV[@]}" KENREF_PARAMS="$toml" "$@" "$MPIRUN" --oversubscribe -np 1 \
          -x GMX_MDMODULE_PLUGIN_PATH -x KENREF_K -x KENREF_N -x KENREF_MAX_FORCE -x KENREF_PARAMS \
          "$GMX" mdrun -nsteps "$NSTEPS" -ntomp 1 -deffnm md ) > "$WORK/$label.log" 2>&1
    local rc=$? energy rel
    energy="$(step0_energy "$WORK/$label.log")"
    printf '%-24s rc=%-3s ' "$label" "$rc"
    if ! grep -a -q "KEnRefMDModule::initForceProviders" "$WORK/$label.log"; then
        echo "*** FAIL: GROMACS did not run the plugin, see $WORK/$label.log"
        fails=$((fails + 1))
    elif [ "$expect" = REFUSE ]; then
        if grep -a -q "KEnRef refuses to continue" "$WORK/$label.log" && [ "$rc" -ne 0 ]; then
            echo "REFUSED as required"
        else
            echo "*** FAIL: expected a refusal, got rc=$rc energy=${energy:-none}"
            fails=$((fails + 1))
        fi
    elif rel="$(close "$energy" "$expect" "$rtol")" && [ "$rc" -eq 0 ]; then
        echo "E0=$energy  expected $expect  rel $rel  ($rtol)  OK"
    else
        echo "*** FAIL: E0=${energy:-none}  expected $expect  rel ${rel:-n/a}  ($rtol), see $WORK/$label.log"
        fails=$((fails + 1))
    fi
}

echo "KEnRef plugin gate"
echo "  gmx      $GMX"
echo "  plugin   $PLUGIN"
echo "  workdir  $WORK"
echo
echo "A/B/C) rt_validate's cases through the plugin:"
run_case gb3_broken_repaired "$GB3_SET/md-00_PBC-BROKEN.tpr" "$WORK/gb3.toml" "$EXPECT_GB3" "$GB3_RTOL"
run_case ubiquitin_whole "$UBQ_SET/t000.00.tpr" "$WORK/ubq.toml" "$EXPECT_UBQ" exact
run_case gb3_repair_disabled "$GB3_SET/md-00_PBC-BROKEN.tpr" "$WORK/gb3.toml" REFUSE exact KENREF_NO_MAKEWHOLE=1

echo
echo "M) two replicas with -multidir, the plugin against bin/KEnRef:"
for kind in bin plugin; do
    for r in rep0 rep1; do
        mkdir -p "$WORK/multi_$kind/$r"
        if [ "$kind" = plugin ]; then
            add_plugin "$DOUBLE_SET/md-00_PBC-BROKEN.tpr" "$WORK/multi_$kind/$r/md.tpr"
        else
            cp "$DOUBLE_SET/md-00_PBC-BROKEN.tpr" "$WORK/multi_$kind/$r/md.tpr"
        fi
    done
done
( cd "$WORK/multi_bin" && "$MPIRUN" --oversubscribe -np 2 "$BIN" --k 1.0 --n 0.25 --max-force 999 \
      --params "$WORK/double.toml" -- -multidir rep0 rep1 -nsteps "$NSTEPS" -ntomp 1 -deffnm md ) \
    > "$WORK/multi_bin.log" 2>&1
bin_rc=$?
( cd "$WORK/multi_plugin" && env "${KENREF_ENV[@]}" KENREF_PARAMS="$WORK/double.toml" "$MPIRUN" --oversubscribe \
      -np 2 -x GMX_MDMODULE_PLUGIN_PATH -x KENREF_K -x KENREF_N -x KENREF_MAX_FORCE -x KENREF_PARAMS \
      "$GMX" mdrun -multidir rep0 rep1 -nsteps "$NSTEPS" -ntomp 1 -deffnm md ) > "$WORK/multi_plugin.log" 2>&1
plugin_rc=$?
bin_e="$(step0_energy "$WORK/multi_bin.log")"
plugin_e="$(step0_energy "$WORK/multi_plugin.log")"
printf '%-24s rc=%s/%s  ' multidir_two_replicas "$bin_rc" "$plugin_rc"
if [ "$bin_rc" -eq 0 ] && [ "$plugin_rc" -eq 0 ] && [ -n "$plugin_e" ] && [ "$plugin_e" = "$bin_e" ]; then
    echo "E0=$plugin_e  bin/KEnRef $bin_e  OK"
else
    echo "*** FAIL: plugin E0=${plugin_e:-none}  bin/KEnRef E0=${bin_e:-none}, see $WORK/multi_*.log"
    fails=$((fails + 1))
fi

echo
if [ "$fails" -eq 0 ]; then
    echo "plugin_validate: PASS (4/4)"
else
    echo "plugin_validate: FAIL ($fails of 4)"
fi
exit $((fails == 0 ? 0 : 1))
