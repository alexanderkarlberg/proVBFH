#!/bin/bash
# Combine the H+3j cross-check runs (VBFNLO, public POWHEG with the scale
# fix) and make the 3-/4-jet plots, ratio to proVBFH-cs (1506.02660 set-up).
#   h3j_update.sh   -> /ptmp/.../combined/p1506/vbfnlo/plots-h3j-{nlo,lo}
set -e
T=$(cd "$(dirname "$0")/../../tools" && pwd)
O=/ptmp/mpp/akarlber/cs-production/combined/p1506/vbfnlo
V=/ptmp/mpp/akarlber/h3j/vbfnlo/prod
W=/ptmp/mpp/akarlber/h3j/powheg/runs
cd $O
python3 - <<'PY'
import re
keep=['sig(all VBF cuts 3 jets)','sig(all VBF cuts 4 jets)','ptj3','yj3','y*j3','min{rap(j1,j3),rap(j3,j2)}','ptj4','yj4']
def sub(src,dst):
    out=[];on=False
    for l in open(src):
        m=re.match(r"\s*#\s*(.*?)\s+index\s+\d+\s*$",l)
        if m: on = m.group(1).strip() in keep
        if on: out.append(l)
    open(dst,'w').write(''.join(out))
import os
R='/u/akarlber/work/proVBFH/proVBFH-cs/production/reference/p1506/'
C='/ptmp/mpp/akarlber/cs-production/combined/p1506/'
for f in ['nnlo-W1','nnlo-W2','nnlo-W3','nlo-W1','nlo-W2','nlo-W3']: sub(C+f+'.top','cs-'+f+'.top')
for f in ['11','HH','22']: sub(R+f+'.top','old-'+f+'.top')
PY
comb() { # out glob
  local files=$(ls $2 2>/dev/null | wc -l)
  [ $files -ge 5 ] || { echo "$1: only $files files, skipped"; return 1; }
  python3 $T/combine_parts.py --error scatter --part "$2" --out $1.tmp | grep -v "^part" || true
  python3 - $1.tmp $1 <<'PY'
import re,sys
keep=['sig(all VBF cuts 3 jets)','sig(all VBF cuts 4 jets)','ptj3','yj3','y*j3','min{rap(j1,j3),rap(j3,j2)}','ptj4','yj4']
out=[];on=False
for l in open(sys.argv[1]):
    m=re.match(r"\s*#\s*(.*?)\s+index\s+\d+\s*$",l)
    if m: on = m.group(1).strip() in keep
    if on: out.append(l)
open(sys.argv[2],'w').write(''.join(out))
PY
  echo "$1: $files files"; }
done_glob() { for d in $1; do [ -f $d/done ] && echo $d/$2; done; }
comb vbfnlo-nlo-3j.top "$V/nlo/job-*/p1506_nlo.top"
comb vbfnlo-lo-3j.top "$V/lo/job-*/p1506_lo.top"
alts=()
comb powheg-nlo-pt1.top "$W/nlo-pt1/pwg-[0-9][0-9][0-9][0-9]-NLO.top" && alts+=(--alt "POWHEG VBF_HJJJ NLO, ptcut 1 GeV" powheg-nlo-pt1.top powheg-nlo-pt1.top powheg-nlo-pt1.top)
comb powheg-nlo-nf5-pt1.top "$W/nlo-nf5nc-nw-pt1/pwg-[0-9][0-9][0-9][0-9]-NLO.top" && alts+=(--alt "POWHEG VBF_HJJJ NLO, like-for-like (5 fl., narrow H, neg. PDFs), ptcut 1 GeV" powheg-nlo-nf5-pt1.top powheg-nlo-nf5-pt1.top powheg-nlo-nf5-pt1.top)
# the plot has four alternative styles: the 0.1 GeV ptcut run (no effect) is
# left out to show the like-for-like run with fixes 1-3 and st_nlight 5
comb powheg-nlo-fix123.top "$W/nlo-lfl-fix123-nl5-pt1/pwg-[0-9][0-9][0-9][0-9]-NLO.top" && alts+=(--alt "POWHEG VBF_HJJJ NLO, like-for-like + fixes 1–3, st_nlight 5" powheg-nlo-fix123.top powheg-nlo-fix123.top powheg-nlo-fix123.top)
loalt=()
comb powheg-lo.top "$W/lo-prod/pwg-[0-9][0-9][0-9][0-9]-NLO.top" && loalt=(--alt "POWHEG VBF_HJJJ LO, 4 flavours" powheg-lo.top powheg-lo.top powheg-lo.top)
comb powheg-lo-nf5.top "$W/lo-nf5fix/pwg-[0-9][0-9][0-9][0-9]-NLO.top" && loalt+=(--alt "POWHEG VBF_HJJJ LO, 5 flavours" powheg-lo-nf5.top powheg-lo-nf5.top powheg-lo-nf5.top)
comb powheg-lo-nf5nc.top "$W/lo-nf5nc-nw/pwg-[0-9][0-9][0-9][0-9]-NLO.top" && loalt+=(--alt "POWHEG VBF_HJJJ LO, like-for-like (5 fl., narrow H, neg. PDFs)" powheg-lo-nf5nc.top powheg-lo-nf5nc.top powheg-lo-nf5nc.top)
rm -rf plots-h3j-nlo plots-h3j-lo
nice python3 $T/plot_compare.py --new vbfnlo-nlo-3j.top vbfnlo-nlo-3j.top vbfnlo-nlo-3j.top \
    --old cs-nnlo-W1.top cs-nnlo-W2.top cs-nnlo-W3.top \
    --alt "proVBFH (old), NNLO" old-11.top old-HH.top old-22.top "${alts[@]}" \
    --oldshort proVBFH-cs --reflabel proVBFH-cs --ratio-range 0.85 1.15 --outdir plots-h3j-nlo \
    --title "1506.02660, O(αs²), ≥3 jets" --json summary-h3j-nlo.json \
    --newlabel "VBFNLO 3.0, NLO H+3j" --oldlabel "proVBFH-cs, NNLO" 2>&1 | grep -v Warn
nice python3 $T/plot_compare.py --new vbfnlo-lo-3j.top vbfnlo-lo-3j.top vbfnlo-lo-3j.top \
    --old cs-nlo-W1.top cs-nlo-W2.top cs-nlo-W3.top "${loalt[@]}" \
    --oldshort proVBFH-cs --reflabel proVBFH-cs --ratio-range 0.9 1.1 --outdir plots-h3j-lo \
    --title "1506.02660, tree-level H+3j" --json summary-h3j-lo.json \
    --newlabel "VBFNLO 3.0, LO H+3j" --oldlabel "proVBFH-cs, NLO (≥3 jets: tree H+3j)" 2>&1 | grep -v Warn
for f in vbfnlo-nlo-3j powheg-nlo-pt1 powheg-nlo-nf5-pt1 powheg-nlo-fix123 vbfnlo-lo-3j powheg-lo powheg-lo-nf5 powheg-lo-nf5nc cs-nnlo-W1 cs-nlo-W1 old-11; do
  [ -f $f.top ] && echo "$f: $(grep -A1 'VBF cuts 3 jets' $f.top | tail -1 | awk '{printf "%.5f +- %.5f", $3, $4}')"; done
