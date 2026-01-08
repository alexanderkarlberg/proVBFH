#!/bin/zsh
# get various hadron/parton ratios

dir=../../gavins-code/results/



obs=vbfcut-stdjets.cross

for obs in vbfcut-stdjets.cross nocut-stdjets.cross vbfcut-nomass-alljets.cross
do
    echo
    echo Observable is $obs
    for mc in pythia6-P11 pythia8-4C herwig-AUET2
    do
        echo $mc: hadron/parton and hadr+UE/parton
        paste <(grep $obs $dir/lhc13-Hgaga-$mc-parton.res) <(grep $obs $dir/lhc13-Hgaga-$mc-hadron.res) <(grep $obs $dir/lhc13-Hgaga-$mc-hadron+UE.res) | awk '{print $21/$7,$35/$7}'
    done
done
