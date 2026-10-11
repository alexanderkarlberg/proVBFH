#!/bin/bash
# submit.sh <name> <array a-b> <time> [BIN]   one single-core job per line of /ptmp/mpp/akarlber/smallx/<name>/jobs.list
R=/ptmp/mpp/akarlber/smallx
EXCL="$(paste -sd, /ptmp/mpp/akarlber/cs-production/bad_nodes),et[07-40]"
exec sbatch --parsable --nice=50 --exclude="$EXCL" --array=$2 --time=$3 --mem=1000MB -J smx-$1 \
  --export=ALL,BIN=$R/bin/smx-run,SMX_BIN=${4:-new},LIST=$R/$1/jobs.list,OFFSET=${5:-0} \
  -o $R/logs/%x.%A_%a.out /u/akarlber/work/proVBFH-smallx/proVBFH-cs/production/slurm/run_array.sh
