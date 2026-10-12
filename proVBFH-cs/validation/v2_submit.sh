#!/bin/bash
# V2 (physics check of the default card): jobs of the default p1506 card at
# NLO or NNLO, one job = inclusive + exclusive part, for the comparison with
# production/results/final-2026-10-07 (v2_check.py). Submits one Slurm array
# on alma (not on et07-et40, see ~/.claude/cluster-notes.md). The same script
# makes the 2-job pilot (NJOBS=2) and the production of the check.
#
#   NJOBS=2 v2_submit.sh <nlo|nnlo> <first seed> [key=value ...]
#
# Environment: NJOBS (default 100), TIME (default 8:00:00), MEM (3000MB), ROOT
# (default /ptmp/mpp/akarlber/csprod/v2), BIN (default ../proVBFH-cs-p1506).
# key=value replaces keys of the card (e.g. ncall2=3200000). Job directories
# job-<seed> under $ROOT/<order>-<tag>, where tag is the first seed; a task
# whose directory is finished exits at once (resubmitting is safe).
set -e
order=$1; first=$2; shift 2
here=$(cd "$(dirname "$0")" && pwd)
prod=$here/../production
NJOBS=${NJOBS:-100}; TIME=${TIME:-8:00:00}; MEM=${MEM:-3000MB}
ROOT=${ROOT:-/ptmp/mpp/akarlber/csprod/v2}
BIN=${BIN:-$here/../proVBFH-cs-p1506}
case $order in nlo) co=2 ;; nnlo) co=3 ;; *) echo "order: nlo or nnlo"; exit 1 ;; esac
run=$ROOT/$order-$first
$prod/slurm/make_jobs.sh p1506 $first $NJOBS $run cs_order=$co "$@"
mkdir -p $run/logs
sbatch --partition=alma --exclude='et[07-40]' --mem=$MEM --time=$TIME -J v2-$order \
    --array=1-$NJOBS --export=ALL,BIN=$BIN,LIST=$run/jobs.list \
    -o $run/logs/%x.%A_%a.out $prod/slurm/run_array.sh
