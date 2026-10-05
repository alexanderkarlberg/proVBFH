#!/bin/bash -l
# usage: sbatch [--array=..] job.sbatch <rundir> <stages> [seed]   (seed defaults to array index)
#SBATCH --partition=alma
#SBATCH --exclude=et36
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2000MB
export OMP_NUM_THREADS=1
seed=${3:-$SLURM_ARRAY_TASK_ID}
echo "host $(hostname) rundir $1 stages $2 seed $seed EXE=${EXE:-default}"
/ptmp/mpp/akarlber/h3j/powheg/runs/run_stages.sh "$1" "$seed" "$2"
