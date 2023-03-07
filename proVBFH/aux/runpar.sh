#!/bin/bash

> Timings.txt
NCORES=8
ITMX=3
EXEC=../proVBFH

# First compile the proVBFH executable in the ../ directory
#

# two stages of importance sampling grid calculation
for igrid in $(seq $ITMX)
do
    
    (echo -n st1 xg$igrid ' ' ; date ) >> Timings.txt
    
    cat powheg.input-save | sed "s/xgriditeration.*/xgriditeration $igrid/ ; s/parallelstage.*/parallelstage 1/ ; s/fakevirt.*/fakevirt 1/ " > powheg.input
    
    for i in $(seq 1 $NCORES)
    do
	echo $i | $EXEC > run-st1-xg$igrid-$i.log 2>&1 &
#echo $i 
    done
    wait
    
done

# compute NLO and upper bounding envelope for underlying born comfigurations
cat powheg.input-save | sed 's/parallelstage.*/parallelstage 2/ ; s/fakevirt.*/fakevirt 0/' > powheg.input
(echo -n st2 ' ' ; date ) >> Timings.txt
for i in $(seq 1 $NCORES)
do
    echo $i | ../EXEC > run-st2-$i.log 2>&1 &
done
wait
(echo -n end ' ' ; date ) >> Timings.txt
exit 0;

