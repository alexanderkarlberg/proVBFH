import numpy as np
d=np.loadtxt('me_terms.dat')
names=['A','B','C','A+B','A+C','B+C']
v={n:(d[:,3*i],d[:,3*i+1],d[:,3*i+2]) for i,n in enumerate(names)}
def terms(k):
    A,B,C=v['A'][k],v['B'][k],v['C'][k]
    return {'AA':A,'BB':B,'CC':C,'2Re AB':v['A+B'][k]-A-B,'2Re AC':v['A+C'][k]-A-C,'2Re BC':v['B+C'][k]-B-C}
N,S,F=terms(0),terms(1),terms(2)
tot=sum(F[t] for t in F).sum()
print("%-8s %14s %14s %14s %11s %11s" % ("term","nontensor","tensor s+a","tensor full","(N-S)/tot","(F-S)/tot"))
for t in N:
    print("%-8s %14.6e %14.6e %14.6e %11.3e %11.3e" % (t,N[t].sum()/tot,S[t].sum()/tot,F[t].sum()/tot,(N[t]-S[t]).sum()/tot,(F[t]-S[t]).sum()/tot))
