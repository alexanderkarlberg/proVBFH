#!/bin/bash
# usage: buildpw.sh driver.f90 exe   (in a scratch copy of proVBFH-cs with objects built; links all objects but cs_main.o)
set -e
D=$1; X=$2
INC="-I/u/akarlber/.local/include/hoppet -Isrc -I../proVBFH/src -I../proVBFH/src/inclusive -I../proVBFH/src/exclusive -I../proVBFH/src/powheg-files -I../proVBFH/src/powheg-files/include -I../proVBFH/analysis -I../proVBFH/src/powheg-files/hjjj-files -I../proVBFH/src/powheg-files/hjjj-files/vbfnlo-files -Jobj"
OBJS=$(ls obj/*.o | grep -v "obj/cs_main.o" | tr '\n' ' ')
# keep the original link order: take it from make.log
ORD=$(sed -n 56p make.log | tr ' ' '\n' | grep '^obj/.*\.o$' | grep -v cs_main.o | tr '\n' ' ')
gfortran -O3 -fno-automatic -ffree-line-length-none $INC -c -o obj/$(basename $D .f90).o $D
gfortran -O3 -fno-automatic -ffree-line-length-none -I/u/akarlber/.local/include/hoppet $ORD obj/$(basename $D .f90).o $(cat libs.txt | sed 's/^LIBS=//') -o $X
