#!/bin/bash

function char {
	case $1 in
	    [1-9])  echo 000$1 ;;
	    [1-9][0-9])  echo 00$1 ;;
	    [1-9][0-9][0-9])  echo 0$1 ;;
	    [1-9][0-9][0-9][0-9])  echo $1 ;;
	esac
}

for i in $(seq 1 1000)
do
    ch=`char $i`
    if [ ! -f provbfh-$ch-NLO-W1.top ]
    then
	echo "echo $i | ../proVBFH-cs"
    fi
done > jobs.txt
thsub -a jobs.txt

#mergedata 1 provbfh-0*-NLO-W1.top; mv fort.12 W1.top
#mergedata 1 provbfh-0*-NLO-W2.top; mv fort.12 W2.top
#mergedata 1 provbfh-0*-NLO-W3.top; mv fort.12 W3.top
