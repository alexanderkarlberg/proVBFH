import re,glob,math,statistics as s
P={}
for f in sorted(glob.glob('job-*/run.log')):
    t=open(f).read()
    for m in re.finditer(r"integration nlo =\s+(\d), Ncall =\s+(\d+), sig\(all VBF cuts 3 jets\) =\s+(\S+) \+-\s+(\S+) fb",t):
        P.setdefault(m[1],[]).append((float(m[3]),float(m[4])))
    m=re.search(r"final result at NLO.*?sigma =\s+(\S+)",t,re.S); P.setdefault('tot',[]).append((float(m[1]),0)) if m else None
for k,x in P.items():
    v=[a for a,_ in x];print(k,len(v),"mean %.3f +- %.3f"%(s.mean(v),s.stdev(v)/math.sqrt(len(v))))
