#!/bin/bash
# Submit an array over lines a..b of a jobs.list with the proven runner.
#   submit.sh cs     <jobs.list> <a-b> [time]
#   submit.sh vbfnlo <jobs.list> <a-b> [time]
set -e
H=$(dirname "$(readlink -f "$0")")
R=/u/akarlber/work/proVBFH/proVBFH-cs/production/slurm/run_array.sh
L=/ptmp/mpp/akarlber/chsplit/logs; mkdir -p $L
EXCL="$(paste -sd, /ptmp/mpp/akarlber/cs-production/bad_nodes),et[07-40]"
if [ "$1" = cs ]; then
  sbatch --parsable --nice=50 --exclude="$EXCL" --array=$3 --time=${4:-18:00:00} --mem=1000MB -J chsplit-cs \
    --export=ALL,CHAN_MULTI=1,BIN=/ptmp/mpp/akarlber/chsplit/bin/proVBFH-cs-p1506-chan,LIST=$2 \
    -o $L/%x.%A_%a.out $R
else
  sbatch --parsable --nice=50 --exclude="$EXCL" --array=$3 --time=${4:-5:00:00} --mem=2000MB -J chsplit-vbfnlo \
    --export=ALL,BIN=/ptmp/mpp/akarlber/chsplit/bin/vbfnlo-chan-run,LIST=$2 \
    -o $L/%x.%A_%a.out $R
fi
