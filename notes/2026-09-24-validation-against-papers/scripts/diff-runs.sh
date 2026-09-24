#!/bin/bash
# Differential LO/NLO runs with the setups of 1506.02660 (H) and 1811.07918 (HH)
set -e
SP=/tmp/claude-1000/-home-karlberg-cernbox-proVBFH-github/60e1988e-5395-46ce-8aa4-480ce2ec5103/scratchpad
WS=$SP/ws
export PATH=$SP/deps/bin:$PATH LHAPDF_DATA_PATH=$SP/deps/share/LHAPDF
setup() { # dir pkg qcd_order
  local d=$1 pkg=$2 ord=$3
  rm -rf $d; mkdir -p $d; cp $WS/$pkg/example/* $d
  if [ $pkg = proVBFH ]; then lhaid=261000; mw=80.398; else lhaid=91500; mw=80.379; fi
  sed -i "s/^lhans\([12]\) .*/lhans\1 $lhaid/; s/^qcd_order .*/qcd_order $ord/; s/^inclusive_only .*/inclusive_only 0/; s/^runningscales .*/runningscales 1/" $d/powheg.input
  [ $pkg = proVBFH ] && sed -i 's/^ebeam\([12]\) .*/ebeam\1 6500d0/' $d/powheg.input
  [ $pkg = proVBFHH ] && sed -i 's/^ebeam\([12]\) .*/ebeam\1 7000d0/' $d/powheg.input
  sed -i "s/^HMASS .*/HMASS 125.0d0/; s/^WMASS .*/WMASS $mw/; s/^WWIDTH .*/WWIDTH 2.141d0/; s/^ZWIDTH .*/ZWIDTH 2.4952d0/; s/^HWIDTH .*/HWIDTH 4.030d-3/" $d/vbfnlo.input
}
for pkg in proVBFH proVBFHH; do
  exe=$WS/$pkg/$pkg
  d=$SP/valruns/$pkg-LO; setup $d $pkg 1
  sed -i 's/^ncall1 .*/ncall1 1000000/; s/^ncall2 .*/ncall2 10000000/' $d/powheg.input
  echo "$(date +%T) $pkg LO start"
  (cd $d; for i in $(seq 8); do echo $i | $exe > run-$i.log 2>&1 & done; wait)
  (cd $d; $WS/proVBFH/aux/combine_runs pwg-LO-*.top > combine.log 2>&1)
  echo "$(date +%T) $pkg LO done"
done
for pkg in proVBFH proVBFHH; do
  exe=$WS/$pkg/$pkg
  d=$SP/valruns/$pkg-NLO; setup $d $pkg 2
  mv $d/powheg.input $d/powheg.input-save
  sed "s|^EXEC=.*|EXEC=$exe|" $WS/proVBFH/aux/runpar.sh > $d/runpar.sh
  echo "$(date +%T) $pkg NLO start"
  (cd $d; bash runpar.sh)
  (cd $d; $WS/proVBFH/aux/combine_runs pwg-*-NLO.top > combine.log 2>&1)
  echo "$(date +%T) $pkg NLO done"
done
echo ALLDONE
