import collections, sys
good=set(int(x) for x in open('csrun/goodev.txt'))
N=97.2708
cs=collections.defaultdict(lambda: collections.defaultdict(float))
for line in open('csrun/pw4_real_details.txt'):
    t=line.split()
    if len(t)<9: continue
    ev=int(t[0])
    if ev in good: cs[ev][(int(t[1]),int(t[2]))]+=float(t[8])*N
def load(fn):
    vb=collections.defaultdict(lambda: collections.defaultdict(float)); ev=None
    for line in open(fn):
        if line.startswith('EV'): ev=int(line.split()[1])
        elif line.startswith('A ') and ev in good:
            t=line.split(); vb[ev][(int(t[1]),int(t[2]))]+=float(t[7])
    return vb
V={n:load('vbrun_%s/pw4_vbfnlo.dat'%n) for n in 'abc'}
pairs=sorted({ab for e in good for ab in V['a'][e]})
def ratio(vb,ab):
    v=sum(vb[e][ab] for e in good); c=sum(cs[e][ab] for e in good)
    return v/c if c else float('nan')
print('pair   full  noint  nocc34   (VB/CS, ME only, 250 events)')
for ab in pairs:
    if ab[0]<=0 or ab[1]==0: continue   # a>0 quarks, b nonzero
    if ab[1]<0 and abs(ab[1])>5: continue
    if abs(ab[0])>4 or abs(ab[1])>4: continue
    print(ab, ' '.join('%.4f'%ratio(V[n],ab) for n in 'abc'))
