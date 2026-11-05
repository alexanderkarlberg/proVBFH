#!/usr/bin/env bash
#
# Single-unresolved limit checks of the identical-quark interference of the
# four-quark reals (cs_4qint 1, proVBFH-cs), after proVBFH-cs has been built.
# validation/v4_4qint/limits_4qint.f90 measures the scalings of the
# interference I and of the square N in the collinear and soft limits, and
# checks the Cauchy-Schwarz bound on I, and that the code's interference
# term agrees with N + I/2 (criteria in that program; it prints
# 'all checks pass' and exits non-zero otherwise). The reference is
# validation/v4_4qint/ref.f90, an explicit Dirac-matrix calculation.
#
# About 4 minutes: run weekly and on demand, not on every push.
#
# Usage: ci/limit-checks.sh          [after make in proVBFH-cs]
#
# Must be run from the repository root. Uses the default example card
# (cs_4qint 1). Run directories are created under $RUNROOT [default: a
# fresh temporary directory].
#
set -euo pipefail

ROOT=$(pwd)
CS=$ROOT/proVBFH-cs
V4=$CS/validation/v4_4qint
RUNDIR="${RUNROOT:-$(mktemp -d)}/limit-checks"
rm -rf "$RUNDIR"
mkdir -p "$RUNDIR"
cd "$RUNDIR"

# the program links the package's objects, except the main program and the
# analysis of the run (as in validation/v4_4qint/run_v4.sh)
OBJS=$(ls "$CS"/obj/*.o | grep -v 'cs_main.o\|_analysis.o\|user_analysis.o')
OBJS="$OBJS $CS/obj/user_analysis.o"
LIBS="-Wl,-rpath,$(lhapdf-config --prefix)/lib -L$(lhapdf-config --prefix)/lib -lLHAPDF"
LIBS="$LIBS -Wl,-rpath,$(hoppet-config --prefix)/lib $(hoppet-config --libs)"
LIBS="$LIBS $(fastjet-config --libs --plugins) -lstdc++"

gfortran -O2 -ffree-line-length-none -fno-automatic -c "$V4/ref.f90" -o ref.o
gfortran -O2 -ffree-line-length-none -fno-automatic -I. "$V4/limits_4qint.f90" ref.o \
    $OBJS $LIBS -o limits_4qint

mkdir -p run
sed -E 's/^cs_4qint .*/cs_4qint 1/' "$CS/example/powheg.input" > run/powheg.input
grep -q '^cs_4qint' run/powheg.input || echo 'cs_4qint 1' >> run/powheg.input
cp "$CS/example/vbfnlo.input" run/

if ! (cd run && ../limits_4qint > limits.log 2>&1); then
    echo "ERROR: limits_4qint failed, last lines of $RUNDIR/run/limits.log:" >&2
    tail -30 run/limits.log >&2
    exit 1
fi
sed -n '/---- summary/,$p' run/limits.log
if ! grep -q 'all checks pass' run/limits.log; then
    echo "ERROR: single-unresolved limits of the interference failed" >&2
    exit 1
fi
echo "All limit checks passed."
