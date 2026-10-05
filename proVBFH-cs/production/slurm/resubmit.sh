#!/bin/bash
# Resubmit the unfinished lines of a jobs.list (no "done" marker), once no task
# of that list is pending or running any more (tasks stuck in COMPLETING on a
# broken node do not count: they never write "done").
#
#   resubmit.sh <name> <binary> <jobs.list> <time limit>
#
# <name> must be the job name used for the list's arrays (the check uses it).
set -e
name=$1; bin=$2; list=$3; tlim=$4
here=$(cd "$(dirname "$0")" && pwd)
# any array of this list counts, also ones on other partitions named <name>-<suffix>
active=$(squeue -u "$USER" -h -r -t PD,R,CF -o "%j" | awk -v n="$name" '$1 == n || index($1, n "-") == 1' | wc -l)
if [ "$active" -gt 0 ]; then echo "$name: $active tasks still pending/running, not resubmitting"; exit 0; fi
missing=$(awk '{ if (system("test -f " $1 "/done")) printf "%s%d", (n++ ? "," : ""), NR }' "$list")
if [ -z "$missing" ]; then echo "$name: all $(wc -l < "$list") done"; exit 0; fi
logs=$(dirname "$(dirname "$(dirname "$list")")")/../logs/$(basename "$(dirname "$(dirname "$(dirname "$list")")")")
# stay under the per-user queue limit (25,000, with a margin); the rest goes
# in a later round
room=$((24800 - $(squeue -u "$USER" -h -r | wc -l)))
[ "$room" -le 0 ] && { echo "$name: no room in the queue"; exit 0; }
missing=$(echo "$missing" | tr ',' '\n' | head -n "$room" | paste -sd,)
n=$(echo "$missing" | tr ',' '\n' | wc -l)
bad=$(paste -sd, /ptmp/mpp/akarlber/cs-production/bad_nodes 2>/dev/null)
id=$(sbatch --parsable ${bad:+--exclude=$bad} --array="$missing" --time="$tlim" -J "$name" \
    --export=ALL,BIN="$bin",LIST="$list",OFFSET=0 \
    -o "$logs/%x.%A_%a.out" "$here/run_array.sh")
echo "$name: resubmitted $n unfinished lines as array $id" | tee -a "$list.submitted"
