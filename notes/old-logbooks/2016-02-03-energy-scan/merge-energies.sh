#!/bin/bash

energ=('1' '7' '10' '13' '14' '20' '30' '40' '60' '80' '100')

echo '# ' $0 $*
echo '# {sqrts, center, min, max}'
for ord in {LO,NLO,NNLO,N3LO}; do
    echo '# '${ord}' cross-section'
    for sqrts in {1,7,10,13,14,20,30,40,60,80,100}; do
	echo $sqrts $(mergeidx.pl -f n3lo_hist_${sqrts}tev.dat ${ord}\.cross | sed -n '2p')
    done
    echo
    echo
done
