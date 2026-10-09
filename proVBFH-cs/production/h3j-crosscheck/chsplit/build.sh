#!/bin/bash -l
# build the channel-split VBFNLO copy (16 cores on alma)
#SBATCH --partition=alma
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16000MB
#SBATCH --time=2:00:00
#SBATCH -J chsplit-vbfnlo-build
#SBATCH -o /ptmp/mpp/akarlber/chsplit/vbfnlo/build.%j.out
set -e
cd /ptmp/mpp/akarlber/chsplit/vbfnlo/src
make distclean > /dev/null 2>&1 || true
./configure --prefix=/ptmp/mpp/akarlber/chsplit/vbfnlo/install --enable-quad \
  --with-LHAPDF=/u/akarlber/.local --with-gsl=/u/akarlber/.local \
  FC=gfortran F77=gfortran CC=gcc CXX=g++ > ../configure.log 2>&1
make -j16 > ../make.log 2>&1
make install > ../install.log 2>&1
echo BUILD_OK
