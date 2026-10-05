import sys,re,collections
# parse pwhg_checklimits: blocks "<label>  emitter e, process f..." / tag line / "alr N" / ratios r0c/r0 + flag
f=sys.argv[1]; lo=int(sys.argv[2]); hi=int(sys.argv[3])
lines=open(f).read().split('\n')
res=[]; i=0
while i<len(lines):
    l=lines[i]
    m=re.match(r'\s*(\S.*?)\s+emitter\s+(\d+), process(.*)',l)
    if m and i+2<len(lines) and lines[i+2].strip().startswith('alr'):
        label=m.group(1); em=int(m.group(2)); proc=m.group(3).split()
        alr=int(lines[i+2].split()[1]); rat=[]; j=i+3
        while j<len(lines):
            w=lines[j].split()
            try: x=float(w[0])
            except: break
            rat.append((x,' '.join(w[1:]))); j+=1
        res.append((label,em,alr,proc,rat)); i=j
    else: i+=1
sel=[r for r in res if lo<=r[2]<=hi]
print('blocks in file: %d; with alr in [%d,%d]: %d'%(len(res),lo,hi,len(sel)))
bylab=collections.defaultdict(list)
for r in sel: bylab[(r[0],r[1])].append(r)
for k,v in sorted(bylab.items()):
    last=[abs(r[4][-1][0]-1) for r in v if r[4]]
    warn=sum(1 for r in v if r[4] and 'WARN' in r[4][-1][1])
    print('%-16s emitter %d: %4d regions, |r0c/r0-1| at the smallest distance: median %.2e max %.2e, WARN flags at last point: %d'%(k[0],k[1],len(v),sorted(last)[len(last)//2],max(last),warn))
for r in sel[:3]:
    print(r[0],'emitter',r[1],'alr',r[2],'process',' '.join(r[3]))
    for x,fl in r[4]: print('    %.10f %s'%(x,fl))
