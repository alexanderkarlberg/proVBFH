#!/bin/bash -l
# Slurm array job: task i runs the binary $BIN in line (i + OFFSET) of $LIST.
#
#   sbatch --array=1-200 --time=8:00:00 -J p1506-excl \
#          --export=ALL,BIN=/path/proVBFH-cs-p1506,LIST=/path/jobs.list[,OFFSET=0] \
#          -o /path/logs/%x.%A_%a.out run_array.sh
#
# OFFSET lets one jobs.list be submitted as several arrays (MaxArraySize and
# the queue limit). A task whose directory already has a finished run.log
# exits without running (resubmitting an array does not repeat work).
#SBATCH --partition=alma
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1000MB
export OMP_NUM_THREADS=1
i=$((SLURM_ARRAY_TASK_ID + ${OFFSET:-0}))
d=$(sed -n "${i}p" "$LIST")
[ -n "$d" ] && [ -d "$d" ] || { echo "task $i: no directory in $LIST"; exit 1; }
cd "$d"
if [ -f run.log ] && grep -q "CPU" run.log && [ -f done ]; then
    echo "task $i: $d already done"; exit 0
fi
echo "task $i: $d on $(hostname), $(date)"
/usr/bin/time -v "$BIN" > run.log 2> time.log
status=$?
[ $status -eq 0 ] && touch done
echo "task $i: exit $status, $(date)"
exit $status
