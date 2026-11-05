#!/bin/bash

for i in {1..12}; do echo $i | ../proVBFH-cs& done
wait

mergedata 1 pwg-LO-*W1*; mv fort.12 W1.top
mergedata 1 pwg-LO-*W2*; mv fort.12 W2.top
mergedata 1 pwg-LO-*W3*; mv fort.12 W3.top
