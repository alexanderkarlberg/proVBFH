import glob,os,re,numpy as np
def val(f,name):
    t=open(f).read().split('\n'); 
    for i,l in enumerate(t):
        if l.startswith('# '+name+' index'):
            for m in t[i+1:i+4]:
                w=m.split()
                if len(w)>=4: return float(w[2].replace('D','E'))
names=['sig(all VBF cuts 3 jets)','sig(all VBF cuts 4 jets)','sig(all VBF cuts 2 jets)']
R={}
for f in sorted(glob.glob('/ptmp/mpp/akarlber/cs-production/fixmh/p1506/excl/job-*/pwg-EXCL.top')):
    if not os.path.exists(os.path.dirname(f)+'/done'): continue
    R[f.split('/')[-2]]=[val(f,n) for n in names]
J=list(R); A=1000*np.array([R[j] for j in J],float)
np.save('fix.npy',A); open('fix.jobs','w').write('\n'.join(J))
n=len(J);print('jobs',n)
for k,nm in enumerate(names):
    x=A[:,k];m=x.mean();sd=x.std(ddof=1);o=np.argsort(-x)
    print(nm,'mean %.3f +- %.3f sd/job %.2f'%(m,sd/np.sqrt(n),sd),'median %.3f'%np.median(x))
    for N in (1,3,5,10,20,50):
        y=np.delete(x,o[:N]);print('  without top %d: %.3f +- %.3f'%(N,y.mean(),y.std(ddof=1)/np.sqrt(len(y))))
    print('  top:',[(J[i],round(x[i],1)) for i in o[:8]])
    print('  frac jobs>mean+5sd_med:',(x>np.median(x)*2).sum(), 'excess share from top10: %.2f'%((x[o[:10]]-np.median(x)).sum()/(x-np.median(x)).sum()))
