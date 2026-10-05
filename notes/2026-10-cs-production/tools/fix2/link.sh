#!/bin/bash
# usage: link.sh BUILD_DIR prog.f out [extra .f files compiled before linking]
P=$1; SRC=$2; OUT=$3; shift 3
FF="gfortran -fno-automatic -ffixed-line-length-none -std=legacy -fallow-argument-mismatch -O1"
OBJS=$(ls $P/obj-gfortran/*.o | grep -v '/pwhg_main.o$')
LHL=$(lhapdf-config --libdir)
$FF -I$P -I$P/../include -I$P/vbfnlo-files "$SRC" "$@" $OBJS -L$P/vbfnlo-files -lvbfnlo \
  -Wl,-rpath,$LHL -L$LHL -lLHAPDF $(fastjet-config --libs --plugins) -lstdc++ -lz -o "$OUT"
