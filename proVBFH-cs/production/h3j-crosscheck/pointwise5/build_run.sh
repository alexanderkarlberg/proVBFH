#!/bin/bash
# Build the LD_PRELOAD library (chsplit VBFNLO 3.0 copy + patches) and run the limit sequences. Login node only.
# usage: build_run.sh WORKDIR NPOINTS     (WORKDIR must hold a vbfnlo.dat with PROCESS=110, NLO only, NLO_POINTS ~16)
H=/ptmp/mpp/akarlber/chsplit/vbfnlo; SRC=$H/src; L=$H/install/lib/VBFNLO; P=$(dirname $(readlink -f $0))
W=${1:-.}; N=${2:-60}; mkdir -p $W/src; cd $W/src
cp $SRC/amplitudes/hjjj/m2s_qqh4j.F v2.F; patch -s v2.F $P/m2s_qqh4j_lim5.patch
cp $SRC/amplitudes/hjjj/qqh4q_sub.F q4.F; patch -s q4.F $P/qqh4q_sub_dump.patch
FL="-fPIC -shared -fno-automatic -ffixed-line-length-none -ffree-line-length-none -fimplicit-none -fbacktrace -DQREAL=REALPART -DQIMAG=IMAGPART -DQCMPLX=COMPLEX -g -O2 -DHAVE_CONFIG_H -DWITH_NLO -I$SRC/include -I$SRC/utilities -I$SRC/amplitudes/hjjj -I$SRC/.. -J."
gfortran $FL -o libt.so v2.F q4.F -L$L -lVBFNLOUtilities -lVBFNLOLoops -lVBFNLOAmplitudes || exit 1
mkdir -p ../run; cd ../run
export LD_LIBRARY_PATH=$L:$LD_LIBRARY_PATH
rm -f fort.9*
# PW_LIM=N: limit sequences for the first N points with lok(1); PW_NOCUTS=1 reproduces the old (no jet cuts on mapped configs) behaviour
PW_LIM=$N LD_PRELOAD=../src/libt.so timeout 1800 $H/install/bin/vbfnlo > run.log 2>&1
python3 $P/agg_limits.py fort.99     # summary table; details: fort.99 (LIM lines), fort.98 (per-subprocess LK, per-dipole LE), fort.96 (jets per config)
