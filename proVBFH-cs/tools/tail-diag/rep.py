import sys,re,numpy as np
j=sys.argv[1]
pts=[];cur=None
for l in open(f'rep-{j}/replay.log'):
    if '===== replay point' in l: cur={'chk':{},'kin':[],'I':None}; pts.append(cur)
    elif cur is None: continue
    elif ' CHK ks' in l:
        f=l.split(); cur['chk'][int(f[2])]=[float(x) for x in f[-5:]]
    elif ' KIN' in l: cur['kin'].append(l.rstrip())
    elif 'integrand (vegas weight 1)' in l: cur['I']=float(l.split()[-1])
c1=[float(l.split()[0]) for l in open(f'job-{j}/cs_spikes.dat')]
for p,c in zip(pts,c1): p['vw']=c/p['I'] if p['I'] else 0; p['c1']=c
vws=np.array([p['vw'] for p in pts]); print('vw quantiles',np.percentile(vws,[0,25,50,75,100]))
st2=[p for p in pts if p['vw']<np.median(vws)]  # stage 2 has smaller weights
print(len(pts),len(st2))
st2.sort(key=lambda p:-p['c1'])
for p in st2[:int(sys.argv[2])]:
    print('c1 %.3f vw %.2e I %.3e'%(p['c1'],p['vw'],p['I']))
    for k in (1,2,3,4,5,6,7,8,9):
        if k in p['chk']: print('   ks',k,'l1 l2 e3 r1 r2',['%.3f'%(x*1) for x in p['chk'][k]])
    if len(sys.argv)>3: print('\n'.join(p['kin']))
