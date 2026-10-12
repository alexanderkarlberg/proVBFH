#!/usr/bin/env python3
"""summ.py <jobs.list or dirs...>: per job, sigma(>=3j) +- err (fb) from pwg-EXCL.top (proVBFH-cs) or the
VBFNLO log (LO, Born+V+I+KP and real-dip, p1506hist estimator = last iteration), with CPU time."""
import sys,re,os
for a in sys.argv[1:]:
    ds=[l.strip() for l in open(a)] if a.endswith('.list') else [a]
    for d in ds:
        n=os.path.basename(d)
        try:
            if os.path.exists(d+'/pwg-EXCL.top'):
                t=open(d+'/pwg-EXCL.top').read()
                m=re.search(r"sig\(all VBF cuts 3 jets\) index\s+\d+\n\s*\S+\s+\S+\s+(\S+)\s+(\S+)",t)
                c=[None,re.findall(r"CPU\s+(\S+) s",open(d+"/run.log").read())[-1]]
                o=re.search(r"cs_order (\d+)",open(d+'/powheg.input').read())[1]; p=re.search(r"cs_piece (\d+)",open(d+'/powheg.input').read())[1]
                print("%s cs order %s piece %s  3j = %.2f +- %.2f fb  cpu %.2f h"%(n,o,p,1e3*float(m[1].replace('D','E')),1e3*float(m[2].replace('D','E')),float(c[1])/3600))
            else:
                t=open(d+'/run.log').read(); r={}
                for m in re.finditer(r"integration nlo =\s+(-?\d), Ncall =\s+(\d+), sig\(all VBF cuts 3 jets\) =\s+(\S+) \+-\s+(\S+) fb",t): r[m[1]]=(float(m[3]),float(m[4]))
                s=int(re.search(r"Runtime: .*\((\d+) seconds",t)[1])
                print("%s VBFNLO  LO %.2f+-%.2f  B+V+I+KP %.2f+-%.2f  real-dip %.2f+-%.2f  cpu %.2f h"%((n,)+r['-1']+r['1']+r['2']+(s/3600,)))
        except Exception as e: print(n,"incomplete",type(e).__name__)
