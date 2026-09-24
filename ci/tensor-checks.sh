#!/usr/bin/env bash
#
# Checks of the tensor code, after proVBFH-inclusive has been built:
#  1. "make check": unit test of tensor.f90 against its previous
#     implementation (tests/tensor_legacy.f90) and explicit formulae;
#  2. the three copies of tensor.f90 are identical;
#  3. with zero widths, provbfhh_incl gives the same cross section with
#     the tensor (-tensorME) and the analytic matrix element, at LO and
#     NNLO (they agree to machine precision pointwise there, so the
#     results agree to all printed digits; the tolerance only allows
#     for rounding in the VEGAS grid adaptation).
#
# Must be run from the repository root. Uses PDF4LHC21_40. Run
# directories are created under $RUNROOT [default: a fresh temporary
# directory].
#
set -euo pipefail

ROOT=$(pwd)
RUNDIR="${RUNROOT:-$(mktemp -d)}/tensor-checks"
TOL=1e-6
rm -rf "$RUNDIR"
mkdir -p "$RUNDIR"

echo "make check (proVBFH-inclusive)"
(cd proVBFH-inclusive && make check)

echo "tensor.f90 copies"
for f in proVBFH/src/inclusive/tensor.f90 proVBFHH/src/inclusive/tensor.f90; do
    if ! cmp -s proVBFH-inclusive/src/tensor.f90 "$f"; then
        echo "ERROR: $f differs from proVBFH-inclusive/src/tensor.f90" >&2
        exit 1
    fi
done
echo "  ok"

args=(-pdf PDF4LHC21_40 -sqrts 14000 -ncall1 20000 -ncall2 50000 -itmx1 1 -itmx2 1
      -iseed 1 -wwidth 0 -zwidth 0 -hwidth 0)
for order in lo nnlo; do
    for me in analytic tensor; do
        dir="$RUNDIR/$order-$me"
        mkdir -p "$dir"
        extra=()
        [ "$me" = tensor ] && extra=(-tensorME)
        if ! (cd "$dir" && "$ROOT/proVBFH-inclusive/provbfhh_incl" -$order "${args[@]}" \
                  ${extra[@]+"${extra[@]}"} > run.log 2>&1); then
            echo "ERROR: provbfhh_incl -$order ($me) failed:" >&2
            tail -50 "$dir/run.log" >&2
            exit 1
        fi
    done
    a=$(grep -v '^#' "$RUNDIR/$order-analytic/xsct_${order}_seed0001.dat" | head -1 | awk '{print $1}')
    t=$(grep -v '^#' "$RUNDIR/$order-tensor/xsct_${order}_seed0001.dat" | head -1 | awk '{print $1}')
    echo "provbfhh_incl -$order, zero widths: analytic $a, tensor $t"
    if ! awk -v a="$a" -v t="$t" -v tol="$TOL" \
            'BEGIN { d = (a - t) / a; if (d < 0) d = -d; exit !(a == a && d <= tol) }'; then
        echo "ERROR: tensor and analytic matrix elements differ by more than $TOL" >&2
        exit 1
    fi
done
echo "All tensor checks passed."
