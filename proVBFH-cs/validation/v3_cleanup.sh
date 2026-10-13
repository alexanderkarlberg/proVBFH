#!/bin/bash
# V3: clean-up regression. The default example card (default binary, user
# analysis of proVBFH) and the same card with the p1506 binary (analysis
# analysis/p1506_analysis.f) are run at cs_order 1, 2 and 3 with a small
# ncall and a fixed seed, with the current build and with a reference
# build made before the clean-up (REFDIR). Every bin of every .top file and
# the xsct file must agree to <= 1e-12 relative (v3_compare.py); files that
# the clean-up no longer writes (cs_spikes.dat, fort.78, ...) are ignored.
# Login node, about 5 minutes for N = 10000.
#
#   v3_cleanup.sh [work dir]
# Environment: REFDIR (default /ptmp/mpp/akarlber/csprod/v3ref: proVBFH-cs
# and proVBFH-cs-p1506 of the pre-clean-up build, commit 5905f1e), BINDIR
# (default ..), N (ncall1 = ncall2 of the exclusive part; default 10000),
# SEED (default 4711), CARD (default ../example/powheg.input).
here=$(cd "$(dirname "$0")" && pwd)
W=${1:-/ptmp/mpp/akarlber/csprod/v3}
REFDIR=${REFDIR:-/ptmp/mpp/akarlber/csprod/v3ref}
BINDIR=${BINDIR:-$here/..}
N=${N:-10000}; SEED=${SEED:-4711}
card=${CARD:-$here/../example/powheg.input}
mk() {
    local d=$1 o=$2; mkdir -p "$d"; rm -f "$d"/*
    sed -E "s/^iseed .*/iseed $SEED/; s/^ncall1 [0-9]+/ncall1 $N/; s/^ncall2 [0-9]+/ncall2 $N/; s/^itmx1 .*/itmx1 2/; s/^itmx2 .*/itmx2 2/; s/^cs_order .*/cs_order $o/" "$card" > "$d/powheg.input"
    cp "$here/../example/vbfnlo.input" "$d/"
}
mkdir -p "$W"
pids=()
for v in default p1506; do
    exe=proVBFH-cs; [ $v = p1506 ] && exe=proVBFH-cs-p1506
    for o in 1 2 3; do
        mk $W/ref-$v-$o $o; ( cd $W/ref-$v-$o && $REFDIR/$exe > run.log 2>&1 ) & pids+=($!)
        mk $W/new-$v-$o $o; ( cd $W/new-$v-$o && $BINDIR/$exe > run.log 2>&1 ) & pids+=($!)
    done
done
wait "${pids[@]}"
rc=0
for v in default p1506; do
    for o in 1 2 3; do
        out=$(python3 $here/v3_compare.py $W/new-$v-$o $W/ref-$v-$o 2>&1); r=$?
        [ $r -eq 0 ] && s=PASS || { s=FAIL; rc=1; }
        echo "V3 $v cs_order $o: $s  ($out)"
    done
done
exit $rc
