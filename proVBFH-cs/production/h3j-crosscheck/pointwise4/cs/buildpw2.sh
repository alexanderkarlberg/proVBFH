#!/bin/bash
# like buildpw.sh, but links the p1506 analysis (obj/p1506_analysis.o) instead of user_analysis.o
set -e
D=$1; X=$2
INC="-I/u/akarlber/.local/include/hoppet -Isrc -I../proVBFH/src -I../proVBFH/src/inclusive -I../proVBFH/src/exclusive -I../proVBFH/src/powheg-files -I../proVBFH/src/powheg-files/include -I../proVBFH/analysis -I../proVBFH/src/powheg-files/hjjj-files -I../proVBFH/src/powheg-files/hjjj-files/vbfnlo-files -Jobj"
ORD=$(sed -n 56p make.log | tr ' ' '\n' | grep '^obj/.*\.o$' | grep -v cs_main.o | sed 's#obj/user_analysis.o#obj/p1506_analysis.o#' | tr '\n' ' ')
gfortran -O3 -fno-automatic -ffree-line-length-none $INC -c -o obj/$(basename $D .f90).o $D
gfortran -O3 -fno-automatic -ffree-line-length-none -I/u/akarlber/.local/include/hoppet $ORD obj/$(basename $D .f90).o $(cat libs.txt | sed 's/^LIBS=//') -o $X
