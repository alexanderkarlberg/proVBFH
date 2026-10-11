#!/bin/bash
# mkset.sh <name> <first seed> <n> <ncall2> [extra powheg.input line] [CHAN_XSLICE value]
#   job directories /ptmp/mpp/akarlber/smallx/<name>/job-<seed> from the fixmh card (ncall2 as given), jobs.list
set -e
R=/ptmp/mpp/akarlber/smallx; base=$R/$1; s0=$2; n=$3
src=/ptmp/mpp/akarlber/cs-production/fixmh/p1506/excl/job-1000201
mkdir -p "$base"
for ((s = s0; s < s0 + n; s++)); do
  d=$base/job-$s
  [ -e "$d" ] && { echo "exists: $d"; exit 1; }
  mkdir "$d"; cp $src/vbfnlo.input "$d/"
  sed -e "s/^iseed .*/iseed $s/" -e "s/^ncall2 .*/ncall2 $4/" $src/powheg.input > "$d/powheg.input"
  [ -n "$5" ] && echo "$5" >> "$d/powheg.input"
  : > "$d/chan.env"
  [ -n "$6" ] && echo "CHAN_XSLICE=\"$6\"" >> "$d/chan.env"
  echo "$d" >> "$base/jobs.list"
done
