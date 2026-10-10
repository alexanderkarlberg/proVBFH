#!/usr/bin/env python3
"""Validation: hi + lo1 + lo2 = all in every bin of a .top pair set.
  xsumcheck.py all.top hi.top lo1.top lo2.top   (column 3 values; max relative deviation)"""
import sys
def rd(f):
    return [float(l.split()[2].replace('D','E')) for l in open(f) if len(l.split()) == 4 and not l.startswith('#')]
a, *s = [rd(f) for f in sys.argv[1:]]
assert all(len(x) == len(a) for x in s)
m = max(abs(sum(x[i] for x in s) - a[i]) for i in range(len(a)))
scale = max(abs(v) for v in a)
w = max((abs(sum(x[i] for x in s) - a[i]) / abs(a[i]), i) for i in range(len(a)) if abs(a[i]) > 1e-3 * scale)
print(f'bins {len(a)}  max |sum-all| {m:.3e} (largest bin {scale:.3e});  max rel dev (bins > 1e-3 max) {w[0]:.2e}')
