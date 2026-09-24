#!/bin/bash
# Non-factorisable validation runs against 2005.11334 (13 TeV, VBF cuts)
SP=/tmp/claude-1000/-home-karlberg-cernbox-proVBFH-github/60e1988e-5395-46ce-8aa4-480ce2ec5103/scratchpad
export LHAPDF_DATA_PATH=$SP/lhapdf:/usr/local/share/LHAPDF
source $SP/nf-setup.sh
R=$SP/nfruns
HSTAT=("ncall1 200000" "ncall2 2000000" "itmx1 3" "itmx2 3")
HHSTAT=("ncall1 200000" "ncall2 1000000" "itmx1 2" "itmx2 2")
BSTAT=("ncall1 500000" "ncall2 5000000" "itmx1 3" "itmx2 3")
TT=("box_t_off 1" "box_u_off 1"); BB=("tri_off 1"); ALL=()
ONE=("tri2_off 1" "box2_off 1"); TWO=("tri1_off 1" "box1_off 1")
run() { local d=$1 pkg=$2; shift 2; setup_run $R/$d $pkg "$@"
        (cd $R/$d; s=$SECONDS; echo 1 | nice -n19 $SP/ws/$pkg/$pkg > run.log 2>&1; echo "$(date +%T) $d rc=$? $((SECONDS-s)) s") & }
run H-NF proVBFH "nonfact 1" "${HSTAT[@]}"
for set in TT BB ALL; do
  eval "sw=(\"\${$set[@]}\")"
  run HH-born-$set proVBFHH "tensorME 1" "${sw[@]}" "${BSTAT[@]}"
  run HH-NF1-$set proVBFHH "nonfact 1" "${sw[@]}" "${ONE[@]}" "${HHSTAT[@]}"
  run HH-NF2-$set proVBFHH "nonfact 1" "${sw[@]}" "${TWO[@]}" "${HHSTAT[@]}"
done
run HH-NF-ALL proVBFHH "nonfact 1" "${HHSTAT[@]}"
wait
echo ALLDONE
