#!/bin/bash
# usage: mkrun.sh <name> <bornonly> <ptcut> <ncall1> <itmx1> <ncall2> <itmx2>
set -e
R=/ptmp/mpp/akarlber/h3j/powheg/runs
mkdir -p $R/$1
cp $R/template/vbfnlo.input $R/template/pwgseeds.dat $R/$1/
sed -e "s/BORNONLY/$2/; s/PTCUT/$3/; s/NCALL1/$4/; s/ITMX1/$5/; s/NCALL2/$6/; s/ITMX2/$7/; s/STAGE/1/" \
    $R/template/powheg.input > $R/$1/powheg.input
