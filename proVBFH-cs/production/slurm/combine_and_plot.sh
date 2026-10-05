#!/bin/bash -l
# Combine the runs of one set-up and order per weight and make the comparison
# plots with the old proVBFH results (where they exist).
#
#   combine_and_plot.sh <setup> <order> <out dir> [run root]
#
# <order>: nnlo (prod-nnlo/<setup>/{excl,incl}), nlo (prod-lonlo/<setup>/
# {nlo-excl,nlo-incl}) or lo (prod-lonlo/<setup>/lo-incl). Writes
# <out dir>/<setup>/<order>-W{1,2,3}.top (exclusive + inclusive, error of the
# mean from the seed scatter, no trimming) and, where an old reference
# exists (1506: NNLO; HXSWG: LO, NLO, NNLO), plots-<order>/ and
# summary-<order>.json. Only jobs with a "done" marker are used.
#SBATCH --partition=alma
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4000MB
#SBATCH --time=4:00:00
set -e
setup=$1; order=$2; out=$3/$setup; root=${4:-/ptmp/mpp/akarlber/cs-production}
here=${PROD_DIR:-/u/akarlber/work/proVBFH/proVBFH-cs/production}
tools=$here/../tools
case $order in
    nnlo) excl=$root/prod-nnlo/$setup/excl; incl=$root/prod-nnlo/$setup/incl ;;
    nlo)  excl=$root/prod-lonlo/$setup/nlo-excl; incl=$root/prod-lonlo/$setup/nlo-incl ;;
    lo)   excl=; incl=$root/prod-lonlo/$setup/lo-incl ;;
    *) echo "order must be nnlo, nlo or lo"; exit 1 ;;
esac
mkdir -p "$out"
# finished jobs, as links in one directory per part (combine_parts takes globs)
parts=()
for p in excl incl; do
    src=${!p}; [ -n "$src" ] || continue
    l=$out/.done-$order-$p; rm -rf "$l"; mkdir "$l"
    for d in "$src"/job-*; do [ -f "$d/done" ] && ln -s "$d" "$l/$(basename "$d")"; done
    echo "$order $p: $(ls "$l" | wc -l) finished jobs"
    prefix=EXCL; [ $p = incl ] && prefix=LO
    parts+=("$l/job-*/pwg-$prefix-WEIGHT.top")
done
for w in W1 W2 W3; do
    args=(); for g in "${parts[@]}"; do args+=(--part "${g/WEIGHT/$w}"); done
    python3 "$tools/combine_parts.py" --error scatter "${args[@]}" --out "$out/$order-$w.top" | { grep -v "^part" || true; }
done
ref=$here/reference/$setup
old=; extra=(--oldlabel "proVBFH (old) ${order^^}")
case $setup-$order in
    p1506-nnlo) old="$ref/11.top $ref/HH.top $ref/22.top"; title="1506.02660, NNLO" ;;
    hxswg136-*) old="$ref/$order-central.top $ref/$order-min.top $ref/$order-max.top"
                title="HXSWG 13.6 TeV, ${order^^}"
                if [ $order = nnlo ]; then
                    # the study's trimmed files, plus two other merges of the same seeds
                    extra=(--oldlabel "proVBFH (old) NNLO, study trim" --oldshort study
                           --alt "proVBFH (old) NNLO, plain" "$ref/plain/nnlo-central-plain.top"
                                 "$ref/plain/nnlo-min-plain.top" "$ref/plain/nnlo-max-plain.top"
                           --alt "proVBFH (old) NNLO, sym. trim" "$ref/symtrim/nnlo-central-symtrim.top"
                                 "$ref/symtrim/nnlo-min-symtrim.top" "$ref/symtrim/nnlo-max-symtrim.top")
                fi ;;
esac
# HXSWG plots: fixed ratio range (AK) and the noisy far tails of the *-STXS
# histograms cut (shown and in the chi2)
[ $setup = hxswg136 ] && extra+=(--ratio-range 0.85 1.15 $(cat "$here/slurm/hxswg136-xcuts.txt"))
if [ -n "$old" ]; then
    python3 "$tools/plot_compare.py" --new "$out/$order"-W{1,2,3}.top --old $old \
        --outdir "$out/plots-$order" --title "$title" --json "$out/summary-$order.json" \
        --newlabel "proVBFH-cs ${order^^}" "${extra[@]}" 2>&1 | { grep -v Warning || true; }
fi
