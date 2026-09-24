#!/usr/bin/env bash
#
# Short, low-statistics runs checking that an already-built package runs
# end to end and produces finite output. These are NOT physics
# validations: the numbers are far too imprecise for that.
#
# Usage: ci/smoke-tests.sh <proVBFH|proVBFHH|proVBFH-inclusive>
#
# Must be run from the repository root. Uses the PDF set PDF4LHC21_40
# (LHAPDF ID 93100), which must be installed. Run directories are
# created under $RUNROOT [default: a fresh temporary directory].
#
set -euo pipefail

PKG=${1:?usage: $0 <package>}
ROOT=$(pwd)
RUNDIR="${RUNROOT:-$(mktemp -d)}/$PKG"
PDFNAME=PDF4LHC21_40
PDFID=93100
NCALL=20000

rm -rf "$RUNDIR"
mkdir -p "$RUNDIR"

# fail if a file is missing, empty, or contains NaN/Inf
check_output() {
    local f=$1
    if [ ! -s "$f" ]; then
        echo "ERROR: expected output $f is missing or empty" >&2
        exit 1
    fi
    if grep -Eiq '(^|[^a-z])(nan|inf|infinity)([^a-z]|$)' "$f"; then
        echo "ERROR: $f contains NaN or Inf:" >&2
        cat "$f" >&2
        exit 1
    fi
}

# run a command, sending its output to a log that is printed on failure
run_logged() {
    local log=$1; shift
    local start=$SECONDS
    if ! "$@" > "$log" 2>&1; then
        echo "ERROR: '$*' failed, last lines of $log:" >&2
        tail -50 "$log" >&2
        exit 1
    fi
    echo "  ok ($((SECONDS - start)) s)"
}

# ----------------------------------------------------------------------
# proVBFH-inclusive: N3LO runs (includes LO, NLO, NNLO) for H and HH,
# plus a 7-point scale variation
# ----------------------------------------------------------------------
smoke_inclusive() {
    local exe args
    args=(-pdf $PDFNAME -ncall1 $NCALL -ncall2 $NCALL -itmx1 1 -itmx2 1)
    for exe in provbfh_incl provbfhh_incl; do
        mkdir -p "$RUNDIR/$exe" "$RUNDIR/$exe-scales"

        echo "$exe -n3lo"
        (cd "$RUNDIR/$exe" && run_logged run.log \
            "$ROOT/proVBFH-inclusive/$exe" -n3lo "${args[@]}" -iseed 1)
        check_output "$RUNDIR/$exe/xsct_n3lo_seed0001.dat"
        cat "$RUNDIR/$exe/xsct_n3lo_seed0001.dat"

        echo "$exe -nlo -7scaleuncert"
        (cd "$RUNDIR/$exe-scales" && run_logged run.log \
            "$ROOT/proVBFH-inclusive/$exe" -nlo -7scaleuncert "${args[@]}" -iseed 2)
        check_output "$RUNDIR/$exe-scales/xsct_nlo_seed0002.dat"
    done
}

# ----------------------------------------------------------------------
# proVBFH / proVBFHH: POWHEG-style runs based on example/
# ----------------------------------------------------------------------

# create a run directory from example/ with reduced statistics;
# extra sed expressions can be passed as arguments
setup_powheg_run() {
    local dir=$1; shift
    mkdir -p "$dir"
    cp -r "$ROOT/$PKG/example/." "$dir"
    sed -i.bak \
        -e "s/^lhans\([12]\) .*/lhans\1 $PDFID/" \
        -e "s/^ncall1 .*/ncall1 $NCALL/" \
        -e "s/^ncall2 .*/ncall2 $NCALL/" \
        -e "s/^itmx1 .*/itmx1 1/" \
        -e "s/^itmx2 .*/itmx2 1/" \
        "$@" "$dir/powheg.input"
}

smoke_powheg() {
    local exe="$ROOT/$PKG/$PKG" dir

    # inclusive-only mode at N3LO
    dir="$RUNDIR/inclusive-n3lo"
    setup_powheg_run "$dir" -e "s/^qcd_order .*/qcd_order 4/" \
                            -e "s/^inclusive_only .*/inclusive_only 1/"
    echo "$PKG inclusive_only N3LO"
    (cd "$dir" && echo 1 | run_logged run.log "$exe")
    check_output "$dir/xsct-n3lo-0001.dat"
    cat "$dir/xsct-n3lo-0001.dat"

    # fully differential NNLO: one grid iteration (stage 1) then stage 2
    dir="$RUNDIR/differential-nnlo"
    setup_powheg_run "$dir" -e "s/^qcd_order .*/qcd_order 3/" \
                            -e "s/^inclusive_only .*/inclusive_only 0/" \
                            -e "s/^xgriditeration .*/xgriditeration 1/" \
                            -e "s/^parallelstage .*/parallelstage 1/"
    echo "$PKG differential NNLO, stage 1"
    (cd "$dir" && echo 1 | run_logged run-st1.log "$exe")
    check_output "$dir/pwggridinfo-btl-xg1-0001.dat"

    sed -i.bak "s/^parallelstage .*/parallelstage 2/" "$dir/powheg.input"
    echo "$PKG differential NNLO, stage 2"
    (cd "$dir" && echo 1 | run_logged run-st2.log "$exe")
    check_output "$dir/pwg-st2-0001-stat.dat"
    check_output "$dir/pwg-0001-NNLO.top"
    cat "$dir/pwg-st2-0001-stat.dat"
}

case "$PKG" in
    proVBFH-inclusive) smoke_inclusive ;;
    proVBFH|proVBFHH)  smoke_powheg ;;
    *) echo "unknown package $PKG" >&2; exit 1 ;;
esac
echo "Run directories are in $RUNDIR"
echo "All smoke tests for $PKG passed."
