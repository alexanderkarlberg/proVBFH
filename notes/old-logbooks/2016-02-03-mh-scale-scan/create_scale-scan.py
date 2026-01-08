#!/usr/bin/env python
#
# Set up multiple folders with different seeds for
# parallel runs of vbfnnlo
#
# Assumes existing vbfnlo.input and powheg.input files
# in current directory, which are used as base for all runs
# (changing only the seed in powheg.input)
#
# Usage example on lxplus in NNLOruns:
# cd [...]/VBFH/vbf_nnlo/
# [edit powheg.input as needed]
# python create_runs.py -NNLO -order-min 3
# cd NNLOruns
# for se in {1..10}; do cd seed_${se};mybsub -q 2nd ../../vbf_w;cd ..; done
#
# .. or equivalently
# python create_runs.py -NNLO -order-min 3 -mybsub
#

# Define which seeds should be used (eg 1,..,10)

order=['LO','NLO','NNLO','N3LO']

import os, shutil, numpy, cmdline, subprocess, hfile, sys

maxseed = cmdline.value('-maxseed',1)
seeds=range(1,maxseed+1)
ordermin = cmdline.value('-order-min',1)

#scales=[('1d0','1d0'),('0.5d0','0.5d0'),('2d0','2d0'),('1d0','0.5d0'),('0.5d0','1d0'),('2d0','1d0'),('1d0','2d0'),('0.25d0','0.25d0'),('0.25d0','1d0'),('1d0','0.25d0'),('4d0','4d0'),('4d0','1d0'),('1d0','4d0')]

scales=[('0.1d0','0.1d0'),('0.25d0','0.25d0'),('0.5d0','0.5d0'),('1d0','1d0'),('2d0','2d0'),('4d0','4d0')]

#scales=[('0.05d0','0.05d0'),('0.1d0','0.1d0'),('0.25d0','0.25d0'),('0.5d0','0.5d0'),('1d0','1d0'),('2d0','2d0'),('4d0','4d0')]


index = 1
allxsc = numpy.zeros((len(scales),5))
for o in order:
    indscale = 0
    for scale in scales:
        directory=o+'runs-scalevar/muf'+scale[0]+'mur'+scale[1]+'/seed_1'
        mu = float(scale[0][:-2])
        fn = directory+'/analysis_'+o.lower()+'_seed0001.top'
        xsc = hfile.get_array(fn,'sig incl cuts index')
        allxsc[indscale,0] = mu
        allxsc[indscale,index] = xsc[:,2]
        indscale += 1

    index += 1
print ('# '+cmdline.cmdline())
print ('# '+o+' scale-scan_cross-section')
print ('# mu LO NLO NNLO N3LO')
print (hfile.reformat(allxsc))
