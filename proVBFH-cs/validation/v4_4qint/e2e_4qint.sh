#!/bin/bash
# V4 end-to-end check of cs_4qint 2 (the interference alone): with the same
# VEGAS grids (readingrid 1) and one main-pass iteration, the three modes see
# the same phase-space points, so bin by bin
#   EXCL(cs_4qint 1) - EXCL(cs_4qint 0) = EXCL(cs_4qint 2)
# up to rounding (e2e_compare.py). p1506 binary, NNLO, default card with a
# small ncall. Login node, a few minutes.   e2e_4qint.sh [work dir] [ncall]
set -e
here=$(cd "$(dirname "$0")" && pwd)
CS=$here/../..
W=${1:-/ptmp/mpp/akarlber/csint/e2e}
N=${2:-20000}
mkdir -p $W/grid
card() { sed -E "s/^iseed .*/iseed 4711/; s/^ncall1 [0-9]+/ncall1 $N/; s/^ncall2 [0-9]+/ncall2 $N/; s/^itmx1 .*/itmx1 2/; s/^itmx2 .*/itmx2 $1/; s/^cs_4qint .*/cs_4qint $2/" $CS/example/powheg.input; }
# grids from a cs_4qint 0 run
card 1 0 > $W/grid/powheg.input; cp $CS/example/vbfnlo.input $W/grid/
( cd $W/grid && $CS/proVBFH-cs-p1506 > run.log 2>&1 )
for m in 0 1 2; do
    d=$W/m$m; mkdir -p $d
    card 1 $m > $d/powheg.input; echo "readingrid 1" >> $d/powheg.input
    cp $CS/example/vbfnlo.input $W/grid/grids*.dat $d/
    ( cd $d && $CS/proVBFH-cs-p1506 > run.log 2>&1 ) &
done
wait
python3 $here/e2e_compare.py $W
