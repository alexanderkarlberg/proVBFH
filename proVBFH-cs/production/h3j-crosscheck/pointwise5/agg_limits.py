import re,sys,collections,statistics as st
f=sys.argv[1] if len(sys.argv)>1 else 'run/fort.99'
ev=None;d=collections.OrderedDict()
for l in open(f):
    if l.startswith('LIMEVENT'): ev=int(l.split()[1])
    elif l.startswith('LIM ') and not l.startswith('LIM skip'):
        name=l[4:18].strip(); w=l[18:].split()
        d.setdefault((name,ev),{})[float(w[0])]=(float(w[1]),float(w[2]))
out=collections.OrderedDict()
for (name,ev),v in d.items():
    try: sing = v[1e-6][0]>0 and v[1e-6][0]/v[1e-4][0]>20 and v[1e-4][0]/v[1e-2][0]>5
    except Exception: sing=False
    out.setdefault(name,[0,0,{}])
    out[name][0]+=1
    if sing:
        out[name][1]+=1
        for lam in (1e-2,1e-4,1e-6,1e-8):
            R,S=v[lam]
            if R>0: out[name][2].setdefault(lam,[]).append(abs((R+S)/R))
print('%-14s %s'%('limit','n_sing/n   lam: median / max |R+D|/R'))
for name,(n,ns,dd) in out.items():
    s=' | '.join('%.0e: %.1e / %.1e'%(lam,st.median(x),max(x)) for lam,x in sorted(dd.items(),reverse=True))
    print('%-13s %2d/%2d  %s'%(name,ns,n,s))
