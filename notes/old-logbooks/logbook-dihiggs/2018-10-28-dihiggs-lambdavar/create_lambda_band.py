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

order=['N3LO']

import os, shutil, numpy, cmdline, subprocess, hfile

if cmdline.present('-lo'):
    order=['LO']

if cmdline.present('-nlo'):
    order=['NLO']

if cmdline.present('-nnlo'):
    order=['NNLO']
    
if cmdline.present('-n3lo'):
    order=['N3LO']

scales=[('1','1'),('H','H'),('2','2'),('1','H'),('H','1'),('2','1'),('1','2')]

if cmdline.present('-3point'):
    scales=scales[:3]
if cmdline.present('-novar'):
    scales=scales[0]
#print scales

#energy in TeV
sqrts=''
if cmdline.present('-e'):
    en=cmdline.value('-e',100)
    sqrts='%itev_'%en
lval=cmdline.value('-l',1.0)
lstr='%slambda%03i' % (sqrts, lval*100)

print ('# '+cmdline.cmdline())
for o in order:
    allxsc   = {}
    allptHH  = {}
    allyHH   = {}
    allmHH   = {}
    allphiHH = {}
    allptH_h = {}
    allyH_h  = {}
    allptH_s = {}
    allyH_s  = {}

    for scale in scales:
        sc=scale[0]+scale[1]
        directory=lstr+'/run_'+o.lower()+'_'+scale[0]+scale[1]
        #directory=o+'runs-'+sqrts+'scale'+sc_choice+'/muf'+scale[0]+'mur'+scale[1]+'/seed_1'
        fn = directory+'/pwg-LO.top'
        xsc = hfile.get_array(fn,'sig incl cuts index')
        ptHH = hfile.get_array(fn,'ptHH')
        yHH  = hfile.get_array(fn,'yHH')
        mHH  = hfile.get_array(fn,'mHH')
        phiHH  = hfile.get_array(fn,'index  18')
        ptH_h = hfile.get_array(fn,'ptH_h')
        yH_h  = hfile.get_array(fn,'yH_h')
        ptH_s = hfile.get_array(fn,'ptH_s')
        yH_s  = hfile.get_array(fn,'yH_s')
        allxsc[sc] = xsc[:,2]
        allptHH[sc] = ptHH[:,2]
        allyHH[sc]  = yHH[:,2]
        allmHH[sc]  = mHH[:,2]
        allphiHH[sc]  = phiHH[:,2]
        allptH_h[sc] = ptH_h[:,2]
        allyH_h[sc]  = yH_h[:,2]
        allptH_s[sc] = ptH_s[:,2]
        allyH_s[sc]  = yH_s[:,2]
        ptHHbin = ptHH[:,:2]
        yHHbin  = yHH[:,:2]
        mHHbin  = mHH[:,:2]
        phiHHbin  = phiHH[:,:2]
        ptH_hbin = ptH_h[:,:2]
        yH_hbin  = yH_h[:,:2]
        ptH_sbin = ptH_s[:,:2]
        yH_sbin  = yH_s[:,:2]
    
    finalxsc=numpy.empty((1,3))
    finalptHH=numpy.empty((len(ptHHbin[:,0]),5))
    finalyHH=numpy.empty((len(yHHbin[:,0]),5))
    finalmHH=numpy.empty((len(mHHbin[:,0]),5))
    finalphiHH=numpy.empty((len(phiHHbin[:,0]),5))
    finalptH_h=numpy.empty((len(ptH_hbin[:,0]),5))
    finalyH_h=numpy.empty((len(yH_hbin[:,0]),5))
    finalptH_s=numpy.empty((len(ptH_sbin[:,0]),5))
    finalyH_s=numpy.empty((len(yH_sbin[:,0]),5))
    finalxsc[:,0]=allxsc['11']
    finalxsc[:,1]=hfile.minimum(allxsc)
    finalxsc[:,2]=hfile.maximum(allxsc)
    finalptHH[:,:2]=ptHHbin
    finalptHH[:,2]=allptHH['11']
    finalptHH[:,3]=hfile.minimum(allptHH)
    finalptHH[:,4]=hfile.maximum(allptHH)
    finalyHH[:,:2]=yHHbin
    finalyHH[:,2] =allyHH['11']
    finalyHH[:,3] =hfile.minimum(allyHH)
    finalyHH[:,4] =hfile.maximum(allyHH)
    finalmHH[:,:2]=mHHbin
    finalmHH[:,2] =allmHH['11']
    finalmHH[:,3] =hfile.minimum(allmHH)
    finalmHH[:,4] =hfile.maximum(allmHH)
    finalphiHH[:,:2]=phiHHbin
    finalphiHH[:,2] =allphiHH['11']
    finalphiHH[:,3] =hfile.minimum(allphiHH)
    finalphiHH[:,4] =hfile.maximum(allphiHH)
    finalptH_h[:,:2]=ptH_hbin
    finalptH_h[:,2]=allptH_h['11']
    finalptH_h[:,3]=hfile.minimum(allptH_h)
    finalptH_h[:,4]=hfile.maximum(allptH_h)
    finalyH_h[:,:2]=yH_hbin
    finalyH_h[:,2] =allyH_h['11']
    finalyH_h[:,3] =hfile.minimum(allyH_h)
    finalyH_h[:,4] =hfile.maximum(allyH_h)
    finalptH_s[:,:2]=ptH_sbin
    finalptH_s[:,2]=allptH_s['11']
    finalptH_s[:,3]=hfile.minimum(allptH_s)
    finalptH_s[:,4]=hfile.maximum(allptH_s)
    finalyH_s[:,:2]=yH_sbin
    finalyH_s[:,2] =allyH_s['11']
    finalyH_s[:,3] =hfile.minimum(allyH_s)
    finalyH_s[:,4] =hfile.maximum(allyH_s)
    print ('# '+o+' cross-section {center, min, max}')
    print (hfile.reformat(finalxsc))
    print ('')
    print ('# '+o+' ptHH {ptmin, ptmax, center, min, max}')
    print (hfile.reformat(finalptHH))
    print ('')
    print ('# '+o+' yHH {ptmin, ptmax, center, min, max}')
    print (hfile.reformat(finalyHH))
    print ('')
    print ('# '+o+' mHH {ptmin, ptmax, center, min, max}')
    print (hfile.reformat(finalmHH))
    print ('')
    print ('# '+o+' phiHH {ptmin, ptmax, center, min, max}')
    print (hfile.reformat(finalphiHH))
    print ('')
    print ('# '+o+' ptH_h {ptmin, ptmax, center, min, max}')
    print (hfile.reformat(finalptH_h))
    print ('')
    print ('# '+o+' yH_h {ptmin, ptmax, center, min, max}')
    print (hfile.reformat(finalyH_h))
    print ('')
    print ('# '+o+' ptH_s {ptmin, ptmax, center, min, max}')
    print (hfile.reformat(finalptH_s))
    print ('')
    print ('# '+o+' yH_s {ptmin, ptmax, center, min, max}')
    print (hfile.reformat(finalyH_s))
    print ('')

