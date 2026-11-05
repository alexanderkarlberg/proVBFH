#!/bin/bash

function char {
	case $1 in
	    [1-9])  echo 000$1 ;;
	    [1-9][0-9])  echo 00$1 ;;
	    [1-9][0-9][0-9])  echo 0$1 ;;
	    [1-9][0-9][0-9][0-9])  echo $1 ;;
	esac
}

for i in $(seq 1 9999)
do
    ch=`char $i`
    if [ ! -f pwg-$ch-W1.top ]
    then
	echo "echo $i | ../proVBFH-cs"
    fi
done > jobs.txt
thsub -a jobs.txt

#mergedata 1 pwg-0*W1*; mv fort.12 W1.top
#mergedata 1 pwg-0*W2*; mv fort.12 W2.top
#mergedata 1 pwg-0*W3*; mv fort.12 W3.top
