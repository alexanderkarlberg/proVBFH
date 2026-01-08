#!/bin/bash

energ=('1' '7' '10' '13' '20' '30' '40' '60' '80' '100')

echo '# ' $0 $*
echo '# {sqrts, pdf-uncert}'
echo '# N3LO cross-section'
for sqrts in {1,7,10,13,20,30,40,60,80,100}; do
	echo $sqrts $(mergeidx.pl -f n3lo_pdfuncert_${sqrts}tev_schemeB_Q10.dat N3LO\.cross | sed -n '2p')
done
echo
echo
echo '# NNLO cross-section'
for sqrts in {1,7,10,13,20,30,40,60,80,100}; do
	echo $sqrts $(mergeidx.pl -f nnlo_pdfuncert_${sqrts}tev_schemeB_Q10.dat NNLO\.cross | sed -n '2p')
done
echo
echo
