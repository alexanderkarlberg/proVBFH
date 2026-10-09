import sys,numpy as np
exec(open('rep.py').read().split("vws=")[0])
s2=[p for p in pts if p['vw']<3e-6]
print(j,'points',len(pts),'with vw<3e-6:',len(s2))
for p in sorted(s2,key=lambda p:-p['c1'])[:3]:
    print(' c1 %.3f vw %.2e'%(p['c1'],p['vw']))
    for k in (1,2,3,4,6,7,8):
        print('   ks',k,['%.4f'%(x*p['vw']) for x in p['chk'][k]])
