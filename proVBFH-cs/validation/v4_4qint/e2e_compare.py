#!/usr/bin/env python3
"""e2e_4qint.sh: bin by bin, EXCL(m1) - EXCL(m0) against EXCL(m2), for W1..W3.

  e2e_compare.py <work dir>
The .top files carry 8 significant digits, so every bin must satisfy
|(m1 - m0) - m2| <= 1e-7 (|m0| + |m1| + |m2|) + 1e-12 |sig(2 jets, m0)|
(rounding of the printed numbers, and of bins that cancel to zero by
construction, such as 'sig incl cuts').
Prints, per scale point, the largest ratio of the two sides over all bins and
the cross sections with cuts; exit code 0 if it is <= 1 everywhere.
"""
import os, sys
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, '..', '..', 'tools'))
from combine_parts import read_top

W = sys.argv[1]
rc = 0
for w in (1, 2, 3):
    h = [read_top(os.path.join(W, 'm%d' % m, 'pwg-EXCL-W%d.top' % w))[0] for m in (0, 1, 2)]
    worst, wname, nb = 0.0, '', 0
    floor = 1e-12 * abs(h[0]['sig(all VBF cuts 2 jets)'][0][2])
    for name in h[2]:
        for b0, b1, b2 in zip(h[0][name], h[1][name], h[2][name]):
            tol = 1e-7 * (abs(b0[2]) + abs(b1[2]) + abs(b2[2])) + floor
            dd = abs((b1[2] - b0[2]) - b2[2])
            d = dd / tol if tol > 0 else (0.0 if dd == 0 else float('inf'))
            nb += 1
            if d > worst:
                worst, wname = d, name
    ok = worst <= 1.0
    rc |= not ok
    print('W%d: %d bins, max |(m1-m0)-m2| / tol = %.2f (%s)  %s' % (w, nb, worst, wname, 'PASS' if ok else 'FAIL'))
    for name in ('sig(all VBF cuts 2 jets)', 'sig(all VBF cuts 3 jets)', 'sig(all VBF cuts 4 jets)'):
        a0, a1, a2 = (x[name][0][2] for x in h)
        print('    %-26s m1-m0 %13.6e  m2 %13.6e' % (name, a1 - a0, a2))
sys.exit(rc)
