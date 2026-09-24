#!/usr/bin/env bash
#
# Checks of the non-factorisable (NF) corrections, after proVBFH and
# proVBFHH have been built: short NF runs (nonfact 1, from example/)
# with the gluon-mass regulator lambda = MV^2 (default) and
# lambda = 0.01 MV^2 (nf_regfact 0.01). The NF correction does not
# depend on lambda: the log(lambda/MV^2) terms of the 1-loop square and
# of the 2-loop interference cancel point by point. Since the 1-loop
# triangle is computed analytically and the log coefficient of the
# 2-loop triangle from the numerical angular integral, the check also
# tests the angular integration. The two runs see the same phase-space
# points, so the results must agree to rounding (tolerance 1e-6).
#
# Usage: ci/nf-checks.sh [proVBFH|proVBFHH ...]   [default: both]
#
# Must be run from the repository root. Uses PDF4LHC21_40 (LHAPDF ID
# 93100). Run directories are created under $RUNROOT [default: a fresh
# temporary directory].
#
set -euo pipefail

ROOT=$(pwd)
RUNDIR="${RUNROOT:-$(mktemp -d)}/nf-checks"
PKGS=("$@")
[ ${#PKGS[@]} -gt 0 ] || PKGS=(proVBFH proVBFHH)
NCALL=5000
TOL=1e-6
mkdir -p "$RUNDIR"

# sig(all VBF cuts 2 jets) of an LO .top file
vbf_sigma() {
    awk '/all VBF cuts 2 jets/ { getline; print $3; exit }' "$1"
}

for pkg in "${PKGS[@]}"; do
    for reg in 1 0.01; do
        dir="$RUNDIR/$pkg-regfact$reg"
        rm -rf "$dir"
        mkdir -p "$dir"
        cp -r "$ROOT/$pkg/example/." "$dir"
        sed -i.bak \
            -e "s/^lhans\([12]\) .*/lhans\1 93100/" \
            -e "s/^qcd_order .*/qcd_order 1/" \
            -e "s/^inclusive_only .*/inclusive_only 0/" \
            -e "s/^nonfact .*/nonfact 1/" \
            -e "s/^ncall1 .*/ncall1 $NCALL/" \
            -e "s/^ncall2 .*/ncall2 $NCALL/" \
            -e "s/^itmx1 .*/itmx1 1/" \
            -e "s/^itmx2 .*/itmx2 1/" \
            "$dir/powheg.input"
        printf '\nnf_regfact %s\n' "$reg" >> "$dir/powheg.input"
        echo "$pkg NF, nf_regfact $reg"
        if ! (cd "$dir" && echo 1 | "$ROOT/$pkg/$pkg" > run.log 2>&1); then
            echo "ERROR: $pkg NF run (nf_regfact $reg) failed:" >&2
            tail -50 "$dir/run.log" >&2
            exit 1
        fi
    done
    a=$(vbf_sigma "$RUNDIR/$pkg-regfact1/pwg-LO-0001.top")
    b=$(vbf_sigma "$RUNDIR/$pkg-regfact0.01/pwg-LO-0001.top")
    echo "$pkg NF, sigma with VBF cuts: lambda = MV^2: $a, lambda = 0.01 MV^2: $b"
    if ! awk -v a="$a" -v b="$b" -v tol="$TOL" \
            'BEGIN { if (a + 0 == 0) exit 1; d = (a - b) / a; if (d < 0) d = -d; exit !(d <= tol) }'; then
        echo "ERROR: $pkg NF correction depends on the regulator (or is zero)" >&2
        exit 1
    fi
done
echo "All NF checks passed."
