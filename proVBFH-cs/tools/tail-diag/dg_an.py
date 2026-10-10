#!/usr/bin/env python3
"""Analysis of the 0.01 pb stage-2 reruns (build3: cs_dg.dat, cs_dgacc.dat).
cs_dg.dat columns: 0 point(excl_stats(1)) 1 type(1 line piece w+wv+e12, 2 real, 3 dipole, 4 e3) 2 line 3 m 4 n
 5:14 c(kk=1..9) pb [W1 all; 2-5 NC qq,qg,gq,gg; 6-9 CC] 14:17 (w,wv,e12)*vw (type 1) 17:23 x1 x2 Q1 Q2 1-xp z
 23:25 dipole (y or 1-x, z) 25 vw 26:54 p(0:3,1:7) 54:74 xrand.  Weights: contribution to sigma(>=3j) before /itmx2(3)."""
import sys, glob, os, numpy as np
R='/ptmp/mpp/akarlber/diag/r3/'
CH=['all','NCqq','NCqg','NCgq','NCgg','CCqq','CCqg','CCgq','CCgg']
TY={1:'line(w+wv+e12)',2:'real',3:'dipole',4:'e3'}
REF=[None,29.848,3.894,3.988,-0.135,71.929,8.322,8.661,-0.233]  # pilot channel means, fb (8 Oct)
def top(f,name='sig(all VBF cuts 3 jets)'):
    t=open(f).read().split('\n')
    for i,l in enumerate(t):
        if l.startswith('# '+name+' index'):
            for m in t[i+1:i+4]:
                w=m.split()
                if len(w)>=4: return 1000*float(w[2].replace('D','E'))
def sig(d):
    out=[]
    for k in range(1,10):
        f=d+'/pwg-EXCL-W%d.top'%k
        out.append(top(f) if os.path.exists(f) else np.nan)
    return np.array(out)
def acc(d):
    L=open(d+'/cs_dgacc.dat').read().split('\n'); i=max(j for j,l in enumerate(L) if l.startswith('points'))
    npt=int(L[i].split()[1]); A=np.zeros((4,9,2)); C=np.zeros((4,2))
    for g in range(4):
        a=L[i+1+2*g].split(':')[1].split(); C[g]=[float(a[1]),float(a[2])]; A[g,:,0]=[float(x) for x in a[3:12]]
        A[g,:,1]=[float(x) for x in L[i+2+2*g].split()[2:11]]
    return npt,A,C
def load(d):
    r=[l.split() for l in open(d+'/cs_dg.dat') if l.strip()]
    return np.array(r,float)
def kin(row):
    p=row[26:54].reshape(7,4); n=int(row[4]); out=[]
    for i in range(3,n):   # partons 4..n (index 3 = H? no: 0,1 in; 2 H; 3.. partons)
        E,px,py,pz=p[i]; pt=np.hypot(px,py); y=0.5*np.log((E+pz)/(E-pz)) if E>abs(pz) else np.sign(pz)*99
        out.append('E%.0f pt%.2g y%+.1f'%(E,pt,y))
    return '; '.join(out)
if __name__=='__main__':
    jobs=sorted(glob.glob(R+'job-*')); S={}
    print('job   sigma3j W1  | rerun - ref per channel (fb)')
    for d in jobs:
        j=d.split('-')[-1]
        if not os.path.exists(d+'/done'): print(j,'not done'); continue
        s=sig(d); prod=top('/ptmp/mpp/akarlber/cs-production/fixmh/p1506/excl/job-%s/pwg-EXCL.top'%j)
        npt,A,C=acc(d); D=load(d); S[j]=(s,prod,A,C,D,npt)
        print('%s rerun W1 %.4f production %.4f diff %.2e | pts %d logged %d'%(j,s[0],prod,s[0]-prod,npt,len(D)))
    print('\nPer job: excess (fb) over pilot channel mean, and logged (>=thr in any channel) sum/3 (fb), by channel')
    for j,(s,prod,A,C,D,npt) in S.items():
        print('job',j,'W1 %.1f'%s[0]); 
        print('  %-5s'%'', ' '.join('%9s'%c for c in CH))
        print('  %-5s'%'excess',' '.join('%9.2f'%(s[k]-(REF[k] if k else 125.52)) for k in range(9)))
        for g in range(4):
            print('  %-5s'%TY[g+1][:5],' '.join('%9.2f'%(1000*A[g,k,1]/3) for k in range(9)),' above; n=%d'%C[g,1])
        print('  %-5s'%'below',' '.join('%9.2f'%(1000*A[:,k,0].sum()/3) for k in range(9)))
        tot=1000*A[:,:,:].sum((0,2))/3
        print('  %-5s'%'total',' '.join('%9.2f'%t for t in tot),' (total of all contributions = sigma(>=3j) incl. stage-2 points only)')
