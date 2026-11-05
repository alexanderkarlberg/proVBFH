#!/bin/bash

mergedata 4 11.top HH.top 22.top
mv fort.12 max.top
mergedata 5 11.top HH.top 22.top
mv fort.12 min.top

mergedata 4 ratio-11.top ratio-HH.top ratio-22.top
mv fort.12 max-ratio.top
mergedata 5 ratio-11.top ratio-HH.top ratio-22.top
mv fort.12 min-ratio.top

cd ..
gnuplot NLO-NNLO-thin.gp
cd -

pdfjam ../NLO-NNLO-thin.pdf 1-3,11 -o bugfix.pdf

