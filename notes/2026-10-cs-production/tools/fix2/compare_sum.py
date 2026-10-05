import sys
def load(f):
    d={}
    for l in open(f):
        w=l.replace('****',' 99').split()
        ip,j=int(w[0]),int(w[1]); fl=tuple(map(int,w[2:9])); e2,e3=int(w[9]),int(w[10])
        r1,r2,r3,s=map(float,w[11:15])
        d[(ip,j)]=(fl,e2,e3,r1,r2,r3,s)
    return d
old=load('fix2_sum-fix13.dat'); new=load('fix2_sum-fix123.dat')
maxo=0; maxn=0; nnc=0; noth=0; nzero=0; fr_out=[]; fr_in1=[]; fr_in2=[]; missing=0
for k in old:
    fl,e2,e3,r1o,_,_,so=old[k]; fln,e2n,e3n,r1,r2,r3,sn=new[k]
    assert fl==fln
    isnc = (e2n>0 or e3n>0)
    if isnc:
        nnc+=1
        if e2n==0 or e3n==0: missing+=1
        if so!=0: maxn=max(maxn,abs(sn/so-1)); fr_out.append(r1/so); fr_in1.append(r2/so); fr_in2.append(r3/so)
        else: nzero+=1
    else:
        noth+=1
        if so!=0 or sn!=0: maxo=max(maxo,abs(sn-so)/max(abs(so),abs(sn)))
print('points x entries: NC outgoing-pair %d (missing partner entries: %d, zero ME: %d), other entries %d'%(nnc,missing,nzero,noth))
print('other entries: max relative difference fix123 vs fix13 = %.3e'%maxo)
print('NC flavour structures: max |sum_entries(fix123)/R(fix13) - 1| = %.3e'%maxn)
import statistics as st
for n,a in (('out',fr_out),('in1',fr_in1),('in2',fr_in2)):
    print('  share of entry %s: min %.4f median %.4f max %.4f'%(n,min(a),st.median(a),max(a)))
