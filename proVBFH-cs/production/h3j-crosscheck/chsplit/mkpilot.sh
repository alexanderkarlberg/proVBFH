#!/bin/bash
# Pilot job directories: proVBFH-cs seeds 1000201-1000300 (the fixmh seeds, so
# W1 reproduces those jobs), VBFNLO 8 channels x seeds 2001-2020 (one jobs.list).
H=$(dirname "$(readlink -f "$0")")
bash $H/mkjobs.sh cs /ptmp/mpp/akarlber/chsplit/cs 1000201 100
V=/ptmp/mpp/akarlber/chsplit/vbfnlo/runs
for b in NC CC; do
  for i in qq qg gq gg; do
    bash $H/mkjobs.sh vbfnlo $V/$b-$i 2001 20 $b $i
    cat $V/$b-$i/jobs.list >> $V/jobs.list
  done
done
wc -l /ptmp/mpp/akarlber/chsplit/cs/jobs.list $V/jobs.list
