#!/bin/bash
# usage: run_stages.sh <rundir> <seed index> <stages, e.g. "1 2">
# Runs pwhg_main for the given parallel stages and seed, in <rundir>.
set -e
EXE=${EXE:-/ptmp/mpp/akarlber/h3j/powheg/POWHEG-BOX-V2/VBF_HJJJ/pwhg_main}
cd "$1"
seed=$2
for st in $3; do
  grep -q "^parallelstage  *$st\b" powheg.input || sed -i "s/^parallelstage .*/parallelstage  $st/" powheg.input
  echo "=== stage $st seed $seed start $(date +%s) $(date)"
  echo $seed | $EXE > run-st${st}-$(printf %04d $seed).log 2>&1
  echo "=== stage $st seed $seed end   $(date +%s) $(date)"
done
