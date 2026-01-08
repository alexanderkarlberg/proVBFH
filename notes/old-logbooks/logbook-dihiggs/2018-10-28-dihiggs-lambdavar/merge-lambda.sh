#!/bin/bash

echo '# ' $0 $*
echo '# {sqrts, center, min, max}'
for ord in {LO,NLO,NNLO,N3LO}; do
    lowcase_ord=`echo $ord | tr '[:upper:]' '[:lower:]'`
    echo '# '${ord}' cross-section'
    for lam in {-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,0.5,1,1.25,1.5,1.75,2,3,4,5,6,7,8,9,10}; do
	floatlam100=`echo "$lam * 100" | bc`
	lam100=`printf '%03.0f\n' $floatlam100`
	echo ${lam}.0 $(mergeidx.pl -f ${lowcase_ord}_hist_lambda${lam100}.dat ${ord}\.cross | sed -n '2p')
    done
    echo
    echo
done
