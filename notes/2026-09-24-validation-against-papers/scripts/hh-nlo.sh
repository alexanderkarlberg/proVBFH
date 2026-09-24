#!/bin/bash
# Wait for the H NLO run to finish, combine it, then run HH NLO with reduced stats
set -e
SP=/tmp/claude-1000/-home-karlberg-cernbox-proVBFH-github/60e1988e-5395-46ce-8aa4-480ce2ec5103/scratchpad
WS=$SP/ws
export PATH=$SP/deps/bin:$PATH LHAPDF_DATA_PATH=$SP/deps/share/LHAPDF
: skip wait
: skip H combine
echo "$(date +%T) (H already done)"
d=$SP/valruns/proVBFHH-NLO; rm -rf $d; mkdir -p $d; cp $WS/proVBFHH/example/* $d
sed -e "s/^lhans\([12]\) .*/lhans\1 91500/; s/^qcd_order .*/qcd_order 2/; s/^inclusive_only .*/inclusive_only 0/; s/^runningscales .*/runningscales 1/" \
    -e "s/^ebeam\([12]\) .*/ebeam\1 7000d0/; s/^ncall2 .*/ncall2 1000000/" $WS/proVBFHH/example/powheg.input > $d/powheg.input-save
sed -i "s/^HMASS .*/HMASS 125.0d0/; s/^WMASS .*/WMASS 80.379/; s/^WWIDTH .*/WWIDTH 2.141d0/; s/^ZWIDTH .*/ZWIDTH 2.4952d0/; s/^HWIDTH .*/HWIDTH 4.030d-3/" $d/vbfnlo.input
rm -f $d/powheg.input
sed "s|^EXEC=.*|EXEC=$WS/proVBFHH/proVBFHH|" $WS/proVBFH/aux/runpar.sh > $d/runpar.sh
echo "$(date +%T) proVBFHH NLO start"
(cd $d; bash runpar.sh)
(cd $d; $WS/proVBFH/aux/combine_runs pwg-*-NLO.top > combine.log 2>&1)
echo "$(date +%T) proVBFHH NLO done"
