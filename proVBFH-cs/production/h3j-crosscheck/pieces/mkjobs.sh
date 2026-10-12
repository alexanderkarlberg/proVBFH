#!/bin/bash
# Job directories for the per-piece study; appends to <base>/jobs.list.
#   mkjobs.sh cs <base> <first seed> <n> <cs_order> <cs_piece mask> <ncall2> [ncall1] [cutoff]
#        proVBFH-cs: the fixmh p1506 card (cs_part 2, cs_scales 1) with cs_order, cs_piece
#        (bits: 1 w, 2 wv, 4 e12 virtual, 8 e12 dipoles, 16 wr, 32 wd, 64 e3; 127 = all)
#   mkjobs.sh vbfnlo <base> <first seed> <n> <log2 points> [m2min GeV]
#        VBFNLO process 110 fixmh card with 2^points per iteration (production: 23)
set -e
mode=$1; base=$2; s0=$3; n=$4
mkdir -p "$base"
if [ "$mode" = cs ]; then
  src=/ptmp/mpp/akarlber/cs-production/fixmh/p1506/excl/job-1000201
  for ((s = s0; s < s0 + n; s++)); do
    d=$base/job-$s; [ -e "$d" ] && { echo "exists: $d"; exit 1; }
    mkdir "$d"; cp $src/vbfnlo.input "$d/"
    sed -e "s/^iseed .*/iseed $s/" -e "s/^ncall2 .*/ncall2 $7/" -e "s/^cs_order .*/cs_order $5/" \
        ${8:+-e "s/^ncall1 .*/ncall1 $8/"} ${9:+-e "s/^cs_cutoff .*/cs_cutoff $9/"} $src/powheg.input > "$d/powheg.input"
    echo "cs_piece $6" >> "$d/powheg.input"
    echo "$d" >> "$base/jobs.list"
  done
elif [ "$mode" = vbfnlo ]; then
  src=/ptmp/mpp/akarlber/h3j/vbfnlo/prod/nlo-fixmh/job-2001
  for ((s = s0; s < s0 + n; s++)); do
    d=$base/job-$s; [ -e "$d" ] && { echo "exists: $d"; exit 1; }
    mkdir "$d"
    for f in anom_HVV.dat anomV.dat cuts.dat ggflo.dat histograms.dat kk_coupl_inp.dat kk_input.dat spin2coupl.dat; do cp $src/$f "$d/"; done
    sed -e "s/^LO_POINTS = .*/LO_POINTS = $5/; s/^NLO_POINTS = .*/NLO_POINTS = $5/" $src/vbfnlo.dat > "$d/vbfnlo.dat"
    echo "SEED = $s" > "$d/random.dat"
    [ -n "$6" ] && echo "VBFNLO_M2MIN=$6" > "$d/piece.env"
    echo "$d" >> "$base/jobs.list"
  done
else echo "usage: see header"; exit 1; fi
