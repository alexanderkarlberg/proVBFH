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

ordermin = cmdline.value('-order-min',1)

#energy in TeV
energy=cmdline.value('-sqrts',27)
sqrts='%itev' % (energy)

scales=[('E','E'),('Q','Q'),('H','H'),('1','1'),('2','2'),('4','4'),('8','8')]
scdic={'E':1.0/8.0,
       'Q':1.0/4.0,
       'H':1.0/2.0,
       '1':1.0,
       '2':2.0,
       '4':4.0,
       '8':8.0}

index = 1
allxsc = numpy.zeros((len(scales),5))
for o in order:
    indscale = 0
    for scale in scales:
        directory=sqrts+'/run_'+o.lower()+'_'+scale[0]+scale[1]
        fn = directory+'/pwg-LO.top'
        mu = scdic[scale[0]]
        xsc = hfile.get_array(fn,'sig incl cuts index')
        allxsc[indscale,0] = mu
        allxsc[indscale,index] = xsc[:,2]
        indscale += 1

    index += 1
print ('# '+cmdline.cmdline())
print ('# '+o+' scale-scan_cross-section')
print ('# mu LO NLO NNLO N3LO')
print (hfile.reformat(allxsc))
