#!/usr/bin/env python3
"""V2: the default card (v2_submit.sh) against the final p1506 production.

  v2_check.py <nlo|nnlo> '<run dir>/job-*/pwg-W%d.top'

Per scale point W1..W3: the cross sections 'sig incl cuts' and with all VBF
cuts for 2, 3 and 4 jets, mean of the finished jobs with the error of the mean
from the seed scatter, against production/results/final-2026-10-07/p1506
(pull = difference / quadrature sum of the errors). Plain means, no trimming.
Also the per-job relative scatter, to size a larger run.
"""
import glob, math, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'tools'))
from combine_parts import read_top, average

order, pat = sys.argv[1], sys.argv[2]
ref_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'production', 'results', 'final-2026-10-07', 'p1506')
names = ['sig incl cuts', 'sig(all VBF cuts 2 jets)', 'sig(all VBF cuts 3 jets)', 'sig(all VBF cuts 4 jets)']
print('%-28s %2s %5s %12s %10s %12s %10s %6s %8s' % ('histogram', 'W', 'jobs', 'new [pb]', 'err', 'ref [pb]', 'err', 'pull', 'job rel.'))
for w in (1, 2, 3):
    files = sorted(glob.glob(pat % w))
    if not files:
        continue
    new, _ = average(files, 'scatter')
    ref, _ = read_top(os.path.join(ref_dir, '%s-W%d.top' % (order, w)))
    for n in names:
        a, b = new[n][0], ref[n][0]
        e = math.hypot(a[3], b[3])
        job = a[3] * math.sqrt(len(files)) / a[2] if a[2] else 0
        print('%-28s %2d %5d %12.6g %10.2g %12.6g %10.2g %6.2f %8.3f' % (n, w, len(files), a[2], a[3], b[2], b[3],
              (a[2] - b[2]) / e if e > 0 else 0, job))
