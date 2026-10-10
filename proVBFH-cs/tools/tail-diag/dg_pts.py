import glob,numpy as np,dg_an as a
for d in sorted(glob.glob(a.R+'job-*')):
    D=a.load(d); j=d[-7:]
    if len(D)==0: continue
    print('job',j)
    for r in D[np.argsort(-abs(D[:,5]))]:
        ch=a.CH[1+int(np.argmax(abs(r[6:10])))] if abs(r[6:10]).max()>abs(r[10:14]).max() else a.CH[5+int(np.argmax(abs(r[10:14])))]
        print('  pt%d %s l%d m%d c=%+.1f fb(all) %s x1=%.3f x2=%.3f Q=%.0f,%.0f 1-xp=%.2g z=%.2g dm=%.2g,%.2g comp=%s'%(r[0],a.TY[int(r[1])][:4],r[2],r[3],1000*r[5]/3,ch,*r[17:21],r[21],r[22],r[23],r[24],np.round(1000*r[14:17]/3,1)))
        print('       ',a.kin(r))
