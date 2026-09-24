import numpy as np, sys
exec(open('check.py').read().split("q1,q2,P1,P2,k1,k2=kin()")[0])
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv)>1 else 1)
q1,q2,P1,P2,k1,k2=kin()
np.savetxt('point.txt', np.vstack([q1,q2,P1,P2,k1,k2]), fmt='%.17e')
def G(P,q,k,lowerit):
    L=(lambda v: g@v) if lowerit else (lambda v: v)
    if k==1: return np.outer(L(q),L(q))/dot(q,q)-g
    if k==2:
        Ph=P-dot(P,q)/dot(q,q)*q; return np.outer(L(Ph),L(Ph))/dot(P,q)
    return E_lower(P,q)/(2*dot(P,q))
def coeffs(k,l,lowerit):
    W1=G(P1,q1,k,lowerit); W2=G(P2,q2,l,lowerit)
    s=lambda A,B,C: S(W1,W2,Mtensor(k1,k2,q1,A,B,C,lowerit)).real
    d={'AA':s(1,0,0),'BB':s(0,1,0),'CC':s(0,0,1)}
    d['AB']=s(1,1,0)-d['AA']-d['BB']; d['AC']=s(1,0,1)-d['AA']-d['CC']; d['BC']=s(0,1,1)-d['BB']-d['CC']
    return d
for (k,l,name) in [(1,1,'F1F1'),(3,3,'F3F3')]:
    for lowerit,lab in [(True,'correct'),(False,'code   ')]:
        d=coeffs(k,l,lowerit)
        print(name, lab, ' '.join(f"{t}={d[t]: .8e}" for t in ['AA','AB','AC','BB','BC','CC']))
