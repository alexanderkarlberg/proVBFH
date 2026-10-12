#!/bin/bash
# V1: the combined run (cs_order 1, 2, 3) reproduces the old binary's
# cs_part 1 (inclusive, qcd_order = order, ncall mapped) + cs_part 2
# (exclusive, cs_order 1 for NLO, 3 for NNLO; qcd_order 2 as in production),
# same seed and card, bin by bin (<= 1e-12 relative for the parts, see
# v1_compare.py for the totals). Login node, a few minutes.
#
#   v1_regression.sh [work dir]
# Environment: REF (reference binary, built from tag cs-study-2026-10, default
# /ptmp/mpp/akarlber/csprod/ref/proVBFH-cs/proVBFH-cs-p1506), BIN (default
# ../proVBFH-cs-p1506), N1, N2 (ncall1, ncall2 of the exclusive part; default
# 20000), SEED (default 4711).
here=$(cd "$(dirname "$0")" && pwd)
W=${1:-/ptmp/mpp/akarlber/csprod/v1}
REF=${REF:-/ptmp/mpp/akarlber/csprod/ref/proVBFH-cs/proVBFH-cs-p1506}
BIN=${BIN:-$here/../proVBFH-cs-p1506}
N1=${N1:-20000}; N2=${N2:-20000}; SEED=${SEED:-4711}
card=$here/../example/powheg.input
mapn() { local n=$1; if [ $n -le 10000 ]; then echo $n; elif [ $n -le 10000000 ]; then echo $((n/10)); else echo $((n/100)); fi; }
# card <dir> key=value ...: example card with the given keys replaced (added if absent)
mk() {
    local d=$1; shift; mkdir -p "$d"; rm -f "$d"/*.top "$d"/*.dat
    sed -E "s/^iseed .*/iseed $SEED/; s/^ncall1 [0-9]+/ncall1 $N1/; s/^ncall2 [0-9]+/ncall2 $N2/; s/^itmx1 .*/itmx1 2/; s/^itmx2 .*/itmx2 2/" "$card" > "$d/powheg.input"
    for kv in "$@"; do
        k=${kv%%=*}; v=${kv#*=}
        if grep -q "^$k " "$d/powheg.input"; then sed -i -E "s/^$k .*/$k $v/" "$d/powheg.input"; else echo "$k $v" >> "$d/powheg.input"; fi
    done
    cp "$here/../example/vbfnlo.input" "$d/"
}
mkdir -p "$W"
pids=()
for o in 1 2 3; do
    mk $W/new-$o cs_order=$o
    ( cd $W/new-$o && $BIN > run.log 2>&1 ) & pids+=($!)
    mk $W/ref-incl-$o qcd_order=$o cs_part=1 ncall1=$(mapn $N1) ncall2=$(mapn $N2)
    ( cd $W/ref-incl-$o && $REF > run.log 2>&1 ) & pids+=($!)
    if [ $o -ge 2 ]; then
        mk $W/ref-excl-$o qcd_order=2 cs_part=2 cs_order=$((o == 2 ? 1 : 3))
        ( cd $W/ref-excl-$o && $REF > run.log 2>&1 ) & pids+=($!)
    fi
done
wait "${pids[@]}"
names=(LO NLO NNLO); rc=0
for o in 1 2 3; do
    ex=$W/ref-excl-$o; [ $o -eq 1 ] && ex=-
    out=$(python3 $here/v1_compare.py $W/new-$o $W/ref-incl-$o $ex 2>&1); r=$?
    [ $r -eq 0 ] && s=PASS || { s=FAIL; rc=1; }
    echo "V1 ${names[$((o-1))]}: $s  ($out)"
done
exit $rc
