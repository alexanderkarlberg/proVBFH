#!/bin/bash

for i in {1..12}; do echo $i | ../proVBFH-cs& done
wait

mergedata 1 provbfh-*-LO-W1.top; mv fort.12 W1.top
mergedata 1 provbfh-*-LO-W2.top; mv fort.12 W2.top
mergedata 1 provbfh-*-LO-W3.top; mv fort.12 W3.top
