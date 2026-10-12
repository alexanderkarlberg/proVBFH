# VBFNLO real-dip (nlo=2 piece) per-job tails analysis; details -> rd_tails_detail.txt
import re,glob,sys,numpy as np
D='/ptmp/mpp/akarlber/h3j/vbfnlo/prod/nlo-fixmh/job-*/run.log'
R=[];S=[];seeds=[]
for f in sorted(glob.glob(D)):
    t=open(f).read(); i=t.find('ndim = 15')
    if i<0: continue
    it=re.findall(r"iteration\s+(\d+):\s*\n\s*integral =\s*(\S+),\s+sigma =\s*(\S+)",t[i:])
    if len(it)!=5: continue
    R.append([float(a[1].replace('E','e')) for a in it]); S.append([float(a[2].replace('E','e')) for a in it]); seeds.append(int(f.split('job-')[1].split('/')[0]))
R=np.array(R);S=np.array(S);seeds=np.array(seeds);n=len(R)
o=open('rd_tails_detail.txt','w')
def P(*a):
    s=' '.join(str(x) for x in a); print(s) if P.v else None; o.write(s+'\n')
P.v=True
def sk(x): d=x-x.mean(); return (d**3).mean()/(d**2).mean()**1.5
def boot(x,B=4000,rng=np.random.default_rng(1)): return np.std([x[rng.integers(0,len(x),len(x))].mean() for _ in range(B)])
def jk(x): m=len(x); l=(x.sum()-x)/(m-1); return np.sqrt((m-1)/m*((l-l.mean())**2).sum())
print("njobs",n)
for k in (2,3,4):
    x=R[:,k];e=S[:,k];m=x.mean();se=x.std(ddof=1)/np.sqrt(n);d=np.abs(x-np.median(x)); srt=np.sort(x-m)
    o_=np.argsort(-np.abs(x-m)); c1=np.sort(x)[::-1]
    # contribution of top 1%/5% by |deviation from median|
    top=lambda f: x[np.argsort(-np.abs(x-np.median(x)))[:int(f*n)]].sum()/n
    P("iter",k+1,"mean %.2f+-%.2f jk %.2f boot %.2f median %.2f skew %.2f  mean-median %.2f"%(m,se,jk(x),boot(x),np.median(x),sk(x),m-np.median(x)))
    P("   top1%% |dev| jobs contribute %.2f fb to mean, top5%% %.2f; largest devs:"%(top(.01),top(.05)),np.round(x[o_[:6]]-m,1),"quoted err",np.round(e[o_[:6]],1))
    h=seeds<np.median(seeds); P("   halves: %.2f+-%.2f vs %.2f+-%.2f"%(x[h].mean(),x[h].std(ddof=1)/np.sqrt(h.sum()),x[~h].mean(),x[~h].std(ddof=1)/np.sqrt((~h).sum())))
    pull=(x-m)/e; P("   pull wrt plain mean: std %.2f, frac|pull|>3 %.3f; median quoted err %.2f vs scatter(MAD*1.48) %.2f vs std %.2f; corr(x,err) %.2f spearman-ish %.2f"%(pull.std(),(abs(pull)>3).mean(),np.median(e),1.48*np.median(abs(x-np.median(x))),x.std(ddof=1),np.corrcoef(x,e)[0,1],np.corrcoef(np.argsort(np.argsort(x)),np.argsort(np.argsort(e)))[0,1]))
# inverse-variance per job, all 5 iterations, and 3-5
for lab,sl in (("ivw1-5",slice(0,5)),("ivw3-5",slice(2,5))):
    w=1/S[:,sl]**2; v=(R[:,sl]*w).sum(1)/w.sum(1); P(lab,"mean %.2f+-%.2f"%(v.mean(),v.std(ddof=1)/np.sqrt(n)),"jk %.2f"%jk(v),"corr(v,err) %.2f"%np.corrcoef(v,1/np.sqrt(w.sum(1)))[0,1])
    # global ivw over jobs
    wj=1/ (1/w.sum(1)); P("   job-weighted by combined err: %.2f"%((v*wj).sum()/wj.sum()))
eq=R[:,2:5].mean(1); P("equal-weight it3-5 per job: mean %.2f+-%.2f jk %.2f boot %.2f median %.2f skew %.2f"%(eq.mean(),eq.std(ddof=1)/np.sqrt(n),jk(eq),boot(eq),np.median(eq),sk(eq)))
nw=np.array([2,4,8])/14.;en=(R[:,2:5]*nw).sum(1); P("Ncall-weighted it3-5: mean %.2f+-%.2f"%(en.mean(),en.std(ddof=1)/np.sqrt(n)))
eq45=R[:,3:5].mean(1); P("equal it4-5: %.2f+-%.2f"%(eq45.mean(),eq45.std(ddof=1)/np.sqrt(n)))
# quantile of per-job err vs value in last iteration
x=R[:,4];e=S[:,4];q=np.quantile(x,[0,.25,.5,.75,1]);
for a,b in zip(q[:-1],q[1:]):
    s_=(x>=a)&(x<=b); P("last-it value quartile [%.1f,%.1f]: median quoted err %.2f"%(a,b,np.median(e[s_])))
# iteration-wise scatter vs quoted error: is err understated for big |dev|?
bias=(R[:,4]-R[:,4].mean()); P("sum ivw-weight shift: weight share of lowest-err 10%% jobs in ivw3-5: ")
