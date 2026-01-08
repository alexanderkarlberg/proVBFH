#!/bin/bash

for en in {14,100}; do
    for ord in {LO,NLO,NNLO,N3LO}; do
	lowcase_ord=`echo $ord | tr '[:upper:]' '[:lower:]'`
	for lam in {-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,0.5,0.8,0.9,1,1.1,1.2,1.25,1.5,1.75,2,3,4,5,6,7,8,9,10}; do
	    floatlam100=`echo "$lam * 100" | bc`
	    lam100=`printf '%03.0f\n' $floatlam100`
	    echo "python2 create_lambda_band.py -l ${lam} -e ${en} -${lowcase_ord} >! ${en}tev_${lowcase_ord}_hist_lambda${lam100}.dat" >> create_hists_exec.sh
	done
    done
done
