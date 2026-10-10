#!/bin/bash -l
#SBATCH -J dg3j
#SBATCH -o /ptmp/mpp/akarlber/diag/r3/slurm.%A_%a.out
#SBATCH --partition=alma
#SBATCH --ntasks=1
#SBATCH --mem=2000MB
#SBATCH --time=16:00:00
#SBATCH --array=0-11
JOBS=(1002198 1001652 1001012 1000297 1001387 1001262 1000975 1001364 1000215 1001513 1002128 1000860)
j=${JOBS[$SLURM_ARRAY_TASK_ID]}
cd /ptmp/mpp/akarlber/diag/r3/job-$j && CHAN_MULTI=1 DG_THR=0.01 /ptmp/mpp/akarlber/diag/build3/proVBFH-cs/proVBFH-cs-p1506 > run.log 2>&1; touch done
