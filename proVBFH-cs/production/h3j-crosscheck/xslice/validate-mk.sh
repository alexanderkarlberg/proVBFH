#!/bin/bash
# validation dirs: cs {old,off,hi,lo1,lo2}, vbfnlo {off,hi,lo1,lo2}
B=/ptmp/mpp/akarlber/chsplit/xtest
for t in old off hi lo1 lo2; do
  mkdir -p $B/cs-$t; cp /ptmp/mpp/akarlber/chsplit/test-cs/A/{powheg.input,vbfnlo.input} $B/cs-$t/
done
for t in off hi lo1 lo2; do
  d=$B/vb-$t; mkdir -p $d; cp /ptmp/mpp/akarlber/chsplit/test-vbfnlo/ALL/job-2001/* $d/ 2>/dev/null
  rm -f $d/chan.env $d/done; rm -rf $d/run.log $d/time.log
  [ $t = off ] && x= || x="CHAN_XSLICE=0.02 $t"
  printf 'CHAN_BOSON=ALL\nCHAN_INIT=ALL\n%s\n' "${x:+CHAN_XSLICE=\"$x\"}" > $d/chan.env
done
