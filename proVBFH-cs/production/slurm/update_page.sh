#!/bin/bash
# Combine everything finished so far and rebuild the status/results page.
#   update_page.sh <page dir>
# Runs on the login node (niced): the batch queue may be slow. Orders with
# fewer than 20 finished exclusive (NLO, NNLO) or inclusive (LO) jobs are skipped.
set -e
page=$1
R=/ptmp/mpp/akarlber/cs-production; C=$R/combined; S=$(cd "$(dirname "$0")" && pwd)
T=$S/../../tools; src=$R/page-src
ndone() { ls "$1"/job-*/done 2>/dev/null | wc -l; }
for s in p1506 hxswg136; do
    [ "$(ndone $R/prod-nnlo/$s/excl)" -ge 20 ] && [ "$(ndone $R/prod-nnlo/$s/incl)" -ge 10 ] && nice bash "$S/combine_and_plot.sh" $s nnlo $C
    [ "$(ndone $R/prod-lonlo/$s/nlo-excl)" -ge 20 ] && [ "$(ndone $R/prod-lonlo/$s/nlo-incl)" -ge 10 ] && nice bash "$S/combine_and_plot.sh" $s nlo $C
    [ "$(ndone $R/prod-lonlo/$s/lo-incl)" -ge 20 ] && nice bash "$S/combine_and_plot.sh" $s lo $C
    opts=(); [ $s = hxswg136 ] && opts=(--ratio-range 0.85 1.15 $(cat "$S/hxswg136-xcuts.txt"))
    nice python3 "$T/plot_orders.py" --dir $C/$s --outdir $C/$s/plots-orders --title "$s" "${opts[@]}" 2>&1 | grep -v Warning || true
done
rm -rf "$page"; mkdir -p "$page"
python3 "$T/build_page.py" --data $C --runs $R --style $src/style.html --static $src/static.html \
    --notes $src/findings.html --out "$page"
