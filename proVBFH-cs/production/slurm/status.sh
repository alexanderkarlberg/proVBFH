#!/bin/bash
# Summary of the production: jobs finished per part, Slurm states, cluster nodes.
#   status.sh [run root, default /ptmp/mpp/akarlber/cs-production]
root=${1:-/ptmp/mpp/akarlber/cs-production}
date
for list in "$root"/prod-*/*/*/jobs.list; do
    d=$(dirname "$list")
    n=$(wc -l < "$list")
    done=$(find "$d" -maxdepth 2 -name done | wc -l)
    printf "%-40s %6d / %6d done\n" "${d#$root/}" "$done" "$n"
done
echo "my jobs by state:"
squeue -u "$USER" -h -r -o "%T" | sort | uniq -c
echo "alma nodes by state:"
sinfo -p alma -h -N -o "%T" | sort | uniq -c
