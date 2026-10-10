#!/bin/bash
# submit.sh <cs|vb> <slice> <array a-b> [time]   (job dirs from mkset.sh)
R=/ptmp/mpp/akarlber/xslice; B=/ptmp/mpp/akarlber/chsplit/bin
A=/u/akarlber/work/proVBFH/proVBFH-cs/production/slurm/run_array.sh
EXCL="$(paste -sd, /ptmp/mpp/akarlber/cs-production/bad_nodes),et[07-40]"
[ $1 = cs ] && { bin=$B/cs-xslice-run; mem=1000MB; t=${4:-8:00:00}; } || { bin=$B/vbfnlo-chan-run; mem=2000MB; t=${4:-6:00:00}; }
exec sbatch --parsable --nice=50 --exclude="$EXCL" --array=$3 --time=$t --mem=$mem -J xs-$1-$2 \
  --export=ALL,BIN=$bin,LIST=$R/$1-$2/jobs.list -o $R/logs/%x.%A_%a.out $A
