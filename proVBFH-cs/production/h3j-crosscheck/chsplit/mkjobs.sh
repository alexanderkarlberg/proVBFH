#!/bin/bash
# Job directories for the channel split (appends to <base>/jobs.list; refuses
# to reuse a directory).
#   mkjobs.sh cs <base> <first seed> <n> [ncall2]
#        proVBFH-cs, the fixmh card (runningscales 0, cs_scales 1) with iseed;
#        run with CHAN_MULTI=1 (W1 all, W2-W9 channels)
#   mkjobs.sh vbfnlo <base> <first seed> <n> <CHAN_BOSON> <CHAN_INIT> [points] [iterations]
#        VBFNLO process 110, the fixmh card, channel in chan.env
set -e
mode=$1; base=$2; s0=$3; n=$4
mkdir -p "$base"
if [ "$mode" = cs ]; then
  src=/ptmp/mpp/akarlber/cs-production/fixmh/p1506/excl/job-1000201
  for ((s = s0; s < s0 + n; s++)); do
    d=$base/job-$s
    [ -e "$d" ] && { echo "exists: $d"; exit 1; }
    mkdir "$d"
    cp $src/vbfnlo.input "$d/"
    sed -e "s/^iseed .*/iseed $s/" ${5:+-e "s/^ncall2 .*/ncall2 $5/"} $src/powheg.input > "$d/powheg.input"
    echo "$d" >> "$base/jobs.list"
  done
elif [ "$mode" = vbfnlo ]; then
  src=/ptmp/mpp/akarlber/h3j/vbfnlo/prod/nlo-fixmh/job-2001
  for ((s = s0; s < s0 + n; s++)); do
    d=$base/job-$s
    [ -e "$d" ] && { echo "exists: $d"; exit 1; }
    mkdir "$d"
    for f in anom_HVV.dat anomV.dat cuts.dat ggflo.dat histograms.dat kk_coupl_inp.dat kk_input.dat spin2coupl.dat; do
      cp $src/$f "$d/"
    done
    sed -e "s/^LO_POINTS = .*/LO_POINTS = ${7:-23}/; s/^NLO_POINTS = .*/NLO_POINTS = ${7:-23}/" \
        -e "s/^LO_ITERATIONS = .*/LO_ITERATIONS = ${8:-5}/; s/^NLO_ITERATIONS = .*/NLO_ITERATIONS = ${8:-5}/" \
        $src/vbfnlo.dat > "$d/vbfnlo.dat"
    echo "SEED = $s" > "$d/random.dat"
    printf 'CHAN_BOSON=%s\nCHAN_INIT=%s\n' "$5" "$6" > "$d/chan.env"
    echo "$d" >> "$base/jobs.list"
  done
else
  echo "usage: see header"; exit 1
fi
