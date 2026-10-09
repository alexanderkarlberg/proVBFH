#!/bin/bash
# Login-node validation of the channel-split VBFNLO: one iteration of 2^17 points
# (LO and NLO), seed 2001, so all runs see the same phase-space points.
#   vbfnlo_tests.sh setup | run <names...>
# V0: the original install (vbfnlo-run); ALL, NC, CC, qq, qg, gq, gg: the copy.
H=$(dirname "$(readlink -f "$0")")
B=/ptmp/mpp/akarlber/chsplit/test-vbfnlo
if [ "$1" = setup ]; then
  bash $H/mkjobs.sh vbfnlo $B/V0 2001 1 ALL ALL 17 1
  rm $B/V0/job-2001/chan.env
  bash $H/mkjobs.sh vbfnlo $B/ALL 2001 1 ALL ALL 17 1
  bash $H/mkjobs.sh vbfnlo $B/NC 2001 1 NC ALL 17 1
  bash $H/mkjobs.sh vbfnlo $B/CC 2001 1 CC ALL 17 1
  for c in qq qg gq gg; do bash $H/mkjobs.sh vbfnlo $B/$c 2001 1 ALL $c 17 1; done
elif [ "$1" = run ]; then
  shift
  for t in "$@"; do
    cd $B/$t/job-2001
    if [ $t = V0 ]; then w=/ptmp/mpp/akarlber/h3j/vbfnlo/vbfnlo-run; else w=/ptmp/mpp/akarlber/chsplit/bin/vbfnlo-chan-run; fi
    nohup /usr/bin/time -v $w > run.log 2> time.log < /dev/null &
    echo "started $t"
  done
fi
