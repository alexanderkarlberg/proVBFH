#!/bin/bash
# Create one directory per job for a Slurm array (see CLUSTER-INSTRUCTIONS.md).
#
#   make_jobs.sh <setup> <excl|incl> <first seed> <njobs> <run dir> [key=value ...]
#
# Writes <run dir>/job-<seed>/{powheg.input,vbfnlo.input} from the card
# production/<setup>/powheg-<part>.input with "iseed <seed>", and appends the
# job directories to <run dir>/jobs.list (line i = array task i). Extra
# key=value arguments replace that key in the card (e.g. ncall2=3200000).
# Refuses to overwrite an existing job directory, so seeds are never reused.
set -e
setup=$1; part=$2; first=$3; n=$4; run=$5; shift 5
here=$(cd "$(dirname "$0")/.." && pwd)
card=$here/$setup/powheg-$part.input
[ -f "$card" ] || { echo "no card $card"; exit 1; }
mkdir -p "$run"
run=$(cd "$run" && pwd)
for ((s = first; s < first + n; s++)); do
    d=$run/job-$s
    [ -e "$d" ] && { echo "$d exists, not overwriting"; exit 1; }
    mkdir "$d"
    sed -E "s/^iseed .*/iseed $s/" "$card" > "$d/powheg.input"
    for kv in "$@"; do
        k=${kv%%=*}; v=${kv#*=}
        grep -q "^$k " "$d/powheg.input" || { echo "key $k not in card"; exit 1; }
        sed -i -E "s/^$k .*/$k $v/" "$d/powheg.input"
    done
    cp "$here/$setup/vbfnlo.input" "$d/"
    echo "$d" >> "$run/jobs.list"
done
echo "$run/jobs.list: $(wc -l < "$run/jobs.list") jobs"
