#!/bin/bash
# mkset.sh <cs|vb> <slice hi|lo1|lo2> <first seed> <n>  -> /ptmp/mpp/akarlber/xslice/<code>-<slice>/job-*, jobs.list
H=/u/akarlber/work/proVBFH/proVBFH-cs/production/h3j-crosscheck/chsplit
R=/ptmp/mpp/akarlber/xslice; code=$1; sl=$2; base=$R/$code-$sl
if [ $code = cs ]; then bash $H/mkjobs.sh cs $base $3 $4; else bash $H/mkjobs.sh vbfnlo $base $3 $4 ALL ALL; fi
for d in $(tail -n $4 $base/jobs.list); do
  [ $code = cs ] && : > $d/chan.env
  echo "CHAN_XSLICE=\"0.02 $sl\"" >> $d/chan.env
done
