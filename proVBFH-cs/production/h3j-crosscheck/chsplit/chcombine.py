#!/usr/bin/env python3
"""Plain equal-weight mean over jobs of .top files (column 3), error = seed scatter / sqrt(n) (column 4).
  chcombine.py OUT.top in1.top in2.top ..."""
import math, sys
out, ins = sys.argv[1], sys.argv[2:]
lines = [open(f).read().split('\n') for f in ins]
n = len(ins)
res = []
for i, l in enumerate(lines[0]):
    p = l.split()
    if len(p) == 4 and not l.startswith('#'):
        v = [float(x[i].split()[2].replace('D', 'E')) for x in lines]
        m = sum(v) / n
        e = math.sqrt(sum((a - m) ** 2 for a in v) / (n - 1) / n)
        res.append(f' {p[0]} {p[1]} {m:.8E} {e:.8E}')
    else:
        res.append(l)
open(out, 'w').write('\n'.join(res))
