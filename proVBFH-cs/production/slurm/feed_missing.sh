#!/bin/bash
# Keep a bounded number of the never-started lines of a jobs.list in the
# queue: submit lines with neither a "done" marker nor a time.log (i.e. never
# started) and not submitted before by this script, as long as the user's
# total queue stays below <cap>. Lines already fed are kept in <list>.fed.
#
#   feed_missing.sh <name> <binary> <jobs.list> <time limit> <cap> [nice]
set -e
name=$1; bin=$2; list=$3; tlim=$4; cap=$5; nice=${6:-0}
here=$(cd "$(dirname "$0")" && pwd)
fed=$list.fed
touch "$fed"
room=$((cap - $(squeue -u "$USER" -h -r | wc -l)))
[ "$room" -le 0 ] && { echo "$name: queue at cap ($cap), nothing submitted"; exit 0; }
lines=$(awk -v fedf="$fed" 'BEGIN { while ((getline l < fedf) > 0) f[l] = 1 }
    { if (!(NR in f) && system("test -f " $1 "/done -o -f " $1 "/time.log")) print NR }' "$list" \
    | head -n "$room" | paste -sd,)
[ -z "$lines" ] && { echo "$name: no never-started lines left"; exit 0; }
logs=$(dirname "$(dirname "$(dirname "$list")")")/../logs/$(basename "$(dirname "$(dirname "$(dirname "$list")")")")
bad=$(paste -sd, /ptmp/mpp/akarlber/cs-production/bad_nodes 2>/dev/null)
id=$(sbatch --parsable --nice="$nice" ${bad:+--exclude=$bad} --array="$lines" --time="$tlim" -J "$name" \
    --export=ALL,BIN="$bin",LIST="$list",OFFSET=0 -o "$logs/%x.%A_%a.out" "$here/run_array.sh")
echo "$lines" | tr ',' '\n' >> "$fed"
n=$(echo "$lines" | tr ',' '\n' | wc -l)
echo "$name: fed $n never-started lines as array $id" | tee -a "$list.submitted"
