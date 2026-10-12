#!/bin/bash
# submit.sh <jobs.list> <array a-b> <time h:mm:ss> <cs|vbfnlo> [mem]   (single-core alma array, bad nodes + et07-40 excluded)
D=/u/akarlber/work/proVBFH/proVBFH-cs/production
EXCL="$(paste -sd, /ptmp/mpp/akarlber/cs-production/bad_nodes),et[07-40]"
BIN=$D/h3j-crosscheck/pieces/$([ "$4" = cs ] && echo cs-pieces-run || echo vbfnlo-pieces-run)
exec sbatch --parsable --nice=50 --exclude="$EXCL" --array=$2 --time=$3 --mem=${5:-1000MB} -J pc-$4 \
  --export=ALL,BIN=$BIN,LIST=$1 -o $(dirname $1)/logs/%x.%A_%a.out $D/slurm/run_array.sh
