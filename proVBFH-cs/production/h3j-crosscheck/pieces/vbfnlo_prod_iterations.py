import re,glob,statistics as s,math
B={}
for f in sorted(glob.glob('job-*/run.log')):
    t=open(f).read()
    it=re.findall(r"iteration\s+(\d+):\s*\n\s*integral =\s*(\S+),\s+sigma =\s*(\S+)",t)
    if len(it)!=11: continue
    for i,(k,v,e) in enumerate(it): B.setdefault(i,[]).append(float(v.replace('E','e')))
    m=re.search(r"QCD real emission:\s+(\S+)",t); B.setdefault('int',[]).append(float(m[1]))
for i in list(range(5,11))+['int']:
    v=B[i];print(i,len(v),"mean %.2f +- %.2f"%(s.mean(v),s.stdev(v)/math.sqrt(len(v))))
