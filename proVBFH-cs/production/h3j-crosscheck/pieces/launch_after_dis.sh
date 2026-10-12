#!/bin/bash
# Self-renewing tiny job: once no NNLOJET/disorder production (dis-nj*, dis-*-d2) is queued, submit the
# proVBFH-cs per-piece production (AK, 12 Oct; VBFNLO M2MIN runs dropped). Otherwise check again in 30 min.
#SBATCH --partition=alma
#SBATCH --ntasks=1
#SBATCH --mem=500MB
#SBATCH --time=0:20:00
P=/u/akarlber/work/proVBFH/proVBFH-cs/production/h3j-crosscheck/pieces
B=/ptmp/mpp/akarlber/pieces/prod
log=$B/launch.log
mkdir -p $B
if squeue -u $USER -h -o %j | grep -qE '^dis-(nj|[EZW]-d2)'; then
    echo "$(date) DIS still queued: $(squeue -u $USER -h -o %j | grep -cE '^dis-')" >> $log
    [ -z "$(squeue -u $USER -h -n pc-launch -t PD -o %i)" ] && \
        sbatch --parsable -J pc-launch --begin=now+30minutes --exclude='et[07-40]' -o $B/launch.%j.out $P/launch_after_dis.sh >> $log
    exit 0
fi
[ -e $B/launched ] && exit 0
touch $B/launched
mkdir -p $B/{lo,v,rd,rd-c4,rd-c6}/logs
{
$P/mkjobs.sh cs $B/lo 7200 20 1 1 4800000 200000 && $P/submit.sh $B/lo/jobs.list 1-20 2:00:00 cs
$P/mkjobs.sh cs $B/v 7300 20 3 6 1200000 100000 && $P/submit.sh $B/v/jobs.list 1-20 4:00:00 cs
$P/mkjobs.sh cs $B/rd 7400 330 3 120 4800000 200000 && $P/submit.sh $B/rd/jobs.list 1-330 10:00:00 cs
$P/mkjobs.sh cs $B/rd-c4 7400 40 3 120 4800000 200000 1d-4 && $P/submit.sh $B/rd-c4/jobs.list 1-40 10:00:00 cs
$P/mkjobs.sh cs $B/rd-c6 7400 40 3 120 4800000 200000 1d-6 && $P/submit.sh $B/rd-c6/jobs.list 1-40 10:00:00 cs
} >> $log 2>&1
echo "$(date) submitted" >> $log
