import math, sys, itertools, collections
# parse dump
def parse(fn, maxev=None):
    ev=None; out=[]
    with open(fn) as f:
        for line in f:
            if line.startswith('EV'):
                t=line.split(); iev=int(t[1]); nk=int(t[2])
                lok=[x=='T' for x in t[7:7+nk]]
                ev={'iev':iev,'lok':lok,'p':{}}; out.append(ev)
                if maxev and len(out)>maxev: break
            elif line.startswith('P ') and ev is not None:
                t=line.split(); nd=int(t[1]); v=[float(x) for x in t[3:27]]
                ev['p'][nd]=[v[4*i:4*i+4] for i in range(6)]
    return out
def rap(p):
    return 0.5*math.log((p[0]+p[3])/(p[0]-p[3])) if p[0]>abs(p[3]) else 1e9*(1 if p[3]>0 else -1)
def eta(p):
    pp=math.sqrt(p[1]**2+p[2]**2+p[3]**2); 
    return 0.5*math.log((pp+p[3])/(pp-p[3])) if pp>abs(p[3]) else 1e9*(1 if p[3]>0 else -1)
def pt(p): return math.hypot(p[1],p[2])
def phi(p): return math.atan2(p[2],p[1])
def dphi(a,b):
    d=abs(a-b)
    return 2*math.pi-d if d>math.pi else d
def cluster(parts, R, use_eta, pn=-1.0):
    js=[list(p) for p in parts]
    def coord(p): return (eta(p) if use_eta else rap(p)), phi(p), pt(p)
    while True:
        n=len(js)
        if n==0: return js
        best=None; bd=1e300
        c=[coord(p) for p in js]
        for i in range(n):
            dib=c[i][2]**(2*pn)   # kt^2p, R^2 factor: compare with dij*(1/R^2) convention: d_iB = kt^2p ; d_ij = min*dR^2/R^2
            if dib<bd: bd=dib; best=('B',i)
        for i in range(n):
            for j in range(i+1,n):
                dr2=(c[i][0]-c[j][0])**2+dphi(c[i][1],c[j][1])**2
                dij=min(c[i][2]**(2*pn),c[j][2]**(2*pn))*dr2/R**2
                if dij<bd: bd=dij; best=('P',i,j)
        if best[0]=='B':
            # i is a final jet; remove it from active list but keep
            return finish(js,best[1],R,use_eta,pn)
        i,j=best[1],best[2]
        js[i]=[a+b for a,b in zip(js[i],js[j])]
        del js[j]
def finish(js,ib,R,use_eta,pn):
    # recursive: move jet ib to final list and continue with the rest
    final=[js[ib]]; rest=js[:ib]+js[ib+1:]
    return final+cluster(rest,R,use_eta,pn)
def passes(parts, use_eta_cluster, use_eta_cut, ptmin=25., ymax=4.5, mjjmin=600., dymin=4.5, opp=True, tag='hardest_after_cuts', ymaxpart=None, njmin=2):
    if ymaxpart is not None:
        parts=[p for p in parts if abs(eta(p))<ymaxpart]
    jets=cluster(parts,0.4,use_eta_cluster)
    yv=(lambda p: eta(p)) if use_eta_cut else (lambda p: rap(p))
    cand=[j for j in jets if pt(j)>ptmin and abs(yv(j))<ymax]
    cand.sort(key=lambda j:-pt(j))
    if len(cand)<njmin: return False
    j1,j2=cand[0],cand[1]
    m2=(j1[0]+j2[0])**2-(j1[1]+j2[1])**2-(j1[2]+j2[2])**2-(j1[3]+j2[3])**2
    ok = pt(j2)>ptmin and m2>mjjmin**2 and abs(yv(j1)-yv(j2))>dymin
    if opp: ok = ok and yv(j1)*yv(j2)<0
    return ok
if __name__=='__main__':
    evs=parse(sys.argv[1], int(sys.argv[2]) if len(sys.argv)>2 else None)
    variants=[('cluster y, cut y',False,False),('cluster eta, cut eta',True,True),('cluster y, cut eta',False,True),('cluster eta, cut y',True,False)]
    res={v[0]:collections.Counter() for v in variants}
    for e in evs:
        for nd,p in e['p'].items():
            if nd==1: parts=[p[2],p[3],p[4],p[5]]
            else:
                parts=[p[2],p[3],p[4]]
                h=[p[0][i]+p[1][i]-sum(x[i] for x in parts) for i in range(4)]
                m2=h[0]**2-h[1]**2-h[2]**2-h[3]**2
                if abs(m2-125.**2)>1.0: continue
                if any(x[0]<=0 or abs(x[0]**2-x[1]**2-x[2]**2-x[3]**2)>1e-6*max(x[0]**2,1) for x in parts): continue
            for name,ce,cc in variants:
                r=passes(parts,ce,cc,njmin=(3 if nd==1 else 3))
                res[name][(e['lok'][nd-1],r,nd==1)]+=1
    for name,c in res.items():
        print(name)
        for isreal in (True,False):
            print('  ','real' if isreal else 'mapped','[VB0,mine0] %d [VB0,mine1] %d [VB1,mine0] %d [VB1,mine1] %d'%(c[(False,False,isreal)],c[(False,True,isreal)],c[(True,False,isreal)],c[(True,True,isreal)]))
