# setup_run <dir> <pkg> [extra powheg.input lines...]
# 13 TeV, NNPDF30_nnlo_as_0118, EW parameters and scale of 2005.11334,
# VBF cuts from example/, LO (qcd_order 1) differential mode
SP=/tmp/claude-1000/-home-karlberg-cernbox-proVBFH-github/60e1988e-5395-46ce-8aa4-480ce2ec5103/scratchpad
setup_run() {
  local d=$1 pkg=$2; shift 2
  rm -rf $d; mkdir -p $d; cp $SP/ws/$pkg/example/* $d
  sed -i "s/^lhans\([12]\) .*/lhans\1 261000/; s/^qcd_order .*/qcd_order 1/; s/^inclusive_only .*/inclusive_only 0/; s/^runningscales .*/runningscales 1/; s/^ebeam\([12]\) .*/ebeam\1 6500d0/; s/^nonfact .*/nonfact 0/" $d/powheg.input
  sed -i "s/^HMASS .*/HMASS 125.0d0/; s/^WMASS .*/WMASS 80.398d0/; s/^WWIDTH .*/WWIDTH 2.141d0/; s/^ZWIDTH .*/ZWIDTH 2.4952d0/; s/^HWIDTH .*/HWIDTH 4.030d-3/" $d/vbfnlo.input
  sed -i -e '$a\' $d/powheg.input   # ensure trailing newline
  for l in "$@"; do k=${l%% *}; sed -i "/^$k /d" $d/powheg.input; echo "$l" >> $d/powheg.input; done
}
