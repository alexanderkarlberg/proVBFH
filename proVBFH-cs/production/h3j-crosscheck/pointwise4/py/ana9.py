import collections, sys
fn='pw4_ab.txt'
rows=[]
for line in open(fn):
    t=line.split()
    ev=int(t[0]); a=int(t[1]); b=int(t[2]); vb=float(t[3]); cs=float(t[4])*97.2708; cu=float(t[5])*97.2708; pp=float(t[6])
    rows.append((ev,a,b,vb,cs,cu,pp))
evs=sorted({r[0] for r in rows}); print('events',len(evs))
def cls(a,b): return (a==0)+(b==0)
tot=collections.defaultdict(lambda:[0.0]*6)
for ev,a,b,vb,cs,cu,pp in rows:
    t=tot[cls(a,b)]; t[0]+=vb; t[1]+=cs; t[2]+=cu; t[3]+=vb*pp; t[4]+=cs*pp; t[5]+=cu*pp
print('class(n gluons in)  ME-only: VB/CS_all  VB/CS_used   pdf-weighted: VB/CS_all  VB/CS_used')
for k,t in sorted(tot.items()):
    print(k,'%.5f %.5f      %.5f %.5f'%(t[0]/t[1],t[0]/t[2],t[3]/t[4],t[3]/t[5]))
S=[sum(t[i] for t in tot.values()) for i in range(6)]
print('ALL','%.5f %.5f      %.5f %.5f'%(S[0]/S[1],S[0]/S[2],S[3]/S[4],S[3]/S[5]))
# per pair worst
pa=collections.defaultdict(lambda:[0.0]*6)
for ev,a,b,vb,cs,cu,pp in rows:
    t=pa[(a,b)]; t[0]+=vb; t[1]+=cs; t[2]+=cu; t[3]+=vb*pp; t[4]+=cs*pp; t[5]+=cu*pp
bad=sorted(pa.items(),key=lambda z:-abs(z[1][3]-z[1][4]))[:12]
print('largest pdf-weighted (VB-CS_all) pairs:')
for ab,t in bad: print(ab,'VB/CS_all %.5f  VB/CS_used %.5f  pdfw %.5f %.5f'%(t[0]/t[1],t[0]/t[2],t[3]/t[4],t[3]/t[5]))
