#!/bin/bash
# V4 (cs_4qint): builds check_4qint against the proVBFH-cs objects (obj/,
# after "make") and runs it with cs_4qint 1, 0 and 2 (default example card), and the single-
# unresolved limit test limits_4qint (cs_4qint 1).
# Login node, about 25 minutes (4 processes).   run_v4.sh [work dir] [npt per mode] [permutations]
set -e
here=$(cd "$(dirname "$0")" && pwd)
CS=$here/../..
W=${1:-/ptmp/mpp/akarlber/csint/v4}
NPT=${2:-20}; NPERM=${3:-24}
mkdir -p "$W"; cd "$W"
OBJS=$(ls $CS/obj/*.o | grep -v "cs_main.o\|_analysis.o\|user_analysis.o")
OBJS="$OBJS $CS/obj/user_analysis.o"
HP=hoppet-config; LH=lhapdf-config; FJ=fastjet-config
LIBS="-Wl,-rpath,$($LH --prefix)/lib -L$($LH --prefix)/lib -lLHAPDF -Wl,-rpath,$($HP --prefix)/lib $($HP --libs) $($FJ --libs --plugins) -lstdc++"
gfortran -O2 -ffree-line-length-none -fno-automatic -c $here/ref.f90 -o ref.o
gfortran -O2 -ffree-line-length-none -fno-automatic -I. $here/check_4qint.f90 ref.o $OBJS $LIBS -o check_4qint
gfortran -O2 -ffree-line-length-none -fno-automatic -I. $here/limits_4qint.f90 ref.o $OBJS $LIBS -o limits_4qint
for f in 1 0 2; do
    mkdir -p flag$f
    sed -E "s/^cs_4qint .*/cs_4qint $f/" $CS/example/powheg.input > flag$f/powheg.input
    grep -q "^cs_4qint" flag$f/powheg.input || echo "cs_4qint $f" >> flag$f/powheg.input
    cp $CS/example/vbfnlo.input flag$f/
    ( cd flag$f && ../check_4qint $NPT $NPERM > check.log 2>&1 ) &
done
( cd flag1 && ../limits_4qint > limits.log 2>&1 ) &
wait
for f in 1 0 2; do echo "== cs_4qint $f"; grep -v "^ rflav\|^ alr_tag\|^ flst_alrtags\|^ flav+tag" flag$f/check.log | grep -A40 "V4: cs_4qint" | grep -v "LHAPDF\|Thanks\|Eur.Phys\|^ *$"; done
echo "== single-unresolved limits (cs_4qint 1)"; sed -n '/---- summary/,$p' flag1/limits.log | grep -v "LHAPDF\|Thanks\|Eur.Phys\|^ *$" | grep "summary\|cc  limit\|collinear\|soft\|Cauchy\|caller"
