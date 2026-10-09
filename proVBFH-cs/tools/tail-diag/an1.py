import glob,collections,numpy as np
nm={2:'NCqq',3:'NCqg',4:'NCgq',5:'NCgg',6:'CCqq',7:'CCqg',8:'CCgq',9:'CCgg'}
pn=['l1','l2','e3','r1','r2']
R=[]
fs=glob.glob('study/job-*/cs_chspikes.dat')
for f in fs:
    for l in open(f):
        w=l.split()
        if len(w)<30 or w[7]!='F': continue
        R.append((int(w[0]),float(w[1]),np.array([float(x) for x in w[2:7]]),f.split('/')[1],w))
print(len(fs),'files',len(R),'stage1 rows')
for thr in (0.1,0.3,1,3,10):
    c=collections.Counter(r[0] for r in R if r[1]>thr)
    print('tot>%g'%thr,' '.join('%s:%d'%(nm[k],c[k]) for k in nm))
print('sum tot^2 (>0.1):',' '.join('%s:%.0f'%(nm[k],sum(r[1]**2 for r in R if r[0]==k)) for k in nm))
for k in (3,4,7,8,2,6):
    c=collections.Counter(pn[int(np.argmax(abs(r[2])))] for r in R if r[0]==k and r[1]>0.3)
    print(nm[k],'dominant piece (>0.3):',dict(c))
# per-job counts of gq vs qg > 1
for k in (3,4,7,8):
    cj=collections.Counter(r[3] for r in R if r[0]==k and r[1]>1); print(nm[k],'jobs with events>1:',len(cj),'max per job',max(cj.values()) if cj else 0)
# top events gq/qg
for k in (4,8,3,7):
    t=sorted([r for r in R if r[0]==k],key=lambda r:-r[1])[:3]
    for r in t: print(nm[k],r[3],'tot %.2f'%r[1],'pieces',np.round(r[2],2),'flags',r[4][8:12])
