#!/bin/bash

energ=('7' '10' '13' '20' '27' '30' '50' '70' '100')

echo '# ' $0 $*
echo '# {sqrts, center, min, max}'
for ord in {LO,NLO,NNLO,N3LO}; do
    echo '# '${ord}' cross-section'
    for sqrts in {7,10,14,20,30,50,70,100}; do
	echo $sqrts $(mergeidx.pl -f n3lo_hist_${sqrts}tev.dat ${ord}\.cross | sed -n '2p')
    done
    echo
    echo
done
