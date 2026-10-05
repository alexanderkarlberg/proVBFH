#!/bin/bash
# Submit the not-yet-submitted lines of a jobs.list as an array, as far as the
# per-user queue limit allows (keeps a margin). Run it repeatedly (e.g. hourly).
#
#   feed.sh <name> <binary> <jobs.list> <time limit> [limit 25000] [margin 200]
#
# State: <jobs.list>.next holds the first line not yet submitted.
set -e
name=$1; bin=$2; list=$3; tlim=$4; limit=${5:-25000}; margin=${6:-200}
here=$(cd "$(dirname "$0")" && pwd)
logs=$(dirname "$(dirname "$(dirname "$list")")")/../logs/$(basename "$(dirname "$(dirname "$(dirname "$list")")")")
mkdir -p "$logs"
state=$list.next
next=$(cat "$state" 2>/dev/null || echo 1)
total=$(wc -l < "$list")
room=$((limit - margin - $(squeue -u "$USER" -h -r | wc -l)))
n=$((total - next + 1)); [ $n -gt $room ] && n=$room
if [ $n -le 0 ]; then echo "$name: next $next of $total, room $room: nothing submitted"; exit 0; fi
bad=$(paste -sd, /ptmp/mpp/akarlber/cs-production/bad_nodes 2>/dev/null)
id=$(sbatch --parsable --nice=100 ${bad:+--exclude=$bad} --array=1-$n --time="$tlim" -J "$name" \
    --export=ALL,BIN="$bin",LIST="$list",OFFSET=$((next - 1)) \
    -o "$logs/%x.%A_%a.out" "$here/run_array.sh")
echo $((next + n)) > "$state"
echo "$name: submitted lines $next-$((next + n - 1)) of $total as array $id" | tee -a "$list.submitted"
