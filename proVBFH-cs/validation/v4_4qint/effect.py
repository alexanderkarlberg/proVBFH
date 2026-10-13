#!/usr/bin/env python3
"""cs_4qint 2 runs (the four-quark interference alone): its effect on the p1506
set-up, against the final production (production/results/final-2026-10-07).

  effect.py '<run dir>/job-*/pwg-EXCL-W%d.top' [target]

Per scale point W1..W3: the interference contribution to sig(>= 2, 3, 4 jets)
(mean of the jobs, error of the mean from the seed scatter), relative to the
NNLO result and to sig(>= 2 jets); the per-job scatter and the number of jobs
for an error of target x sig(>= 2 jets) (default 1e-4). For the distributions
(W1): per histogram, the largest |effect| / error over the bins, and the median
over the bins of error / (NNLO bin value), now and with the sized number of jobs.
Plain means, no trimming.
"""
import glob, math, os, sys
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, '..', '..', 'tools'))
from combine_parts import read_top, average

pat = sys.argv[1]
target = float(sys.argv[2]) if len(sys.argv) > 2 else 1e-4
ref_dir = os.path.join(here, '..', '..', 'production', 'results', 'final-2026-10-07', 'p1506')
names = ['sig(all VBF cuts 2 jets)', 'sig(all VBF cuts 3 jets)', 'sig(all VBF cuts 4 jets)']
nneed = 0
for w in (1, 2, 3):
    files = sorted(glob.glob(pat % w))
    if not files:
        continue
    n = len(files)
    new, _ = average(files, 'scatter')
    ref, _ = read_top(os.path.join(ref_dir, 'nnlo-W%d.top' % w))
    s2 = ref[names[0]][0][2]
    print('W%d, %d jobs; NNLO sig(>=2j) = %.6g pb' % (w, n, s2))
    print('  %-26s %12s %10s %10s %10s %10s %8s' % ('', 'effect [pb]', 'err', 'rel. NNLO', 'rel. 2j', 'err rel.2j', 'jobs'))
    for nm in names:
        a = new[nm][0]
        b = ref[nm][0][2]
        job = a[3] * math.sqrt(n)              # per-job scatter
        need = math.ceil((job / (target * s2))**2)
        if w == 1:
            nneed = max(nneed, need)
        print('  %-26s %12.4e %10.2e %10.2e %10.2e %10.2e %8d' % (nm, a[2], a[3], a[2] / b, a[2] / s2, a[3] / s2, need))
# distributions at W1
files = sorted(glob.glob(pat % 1))
n = len(files)
new, _ = average(files, 'scatter')
ref, _ = read_top(os.path.join(ref_dir, 'nnlo-W1.top'))
print('distributions (W1): max |effect|/err over bins; max |effect|/NNLO over the bins with |effect| > 3 err;')
print('  median err/NNLO bin now (%d jobs) and with %d jobs' % (n, nneed))
for nm in new:
    if nm.startswith('sig') or nm not in ref:
        continue
    sig, rel, eff = [], [], [0.0]
    for a, b in zip(new[nm], ref[nm]):
        if a[3] > 0:
            sig.append(abs(a[2]) / a[3])
            if abs(a[2]) > 3 * a[3] and b[2] != 0:
                eff.append(abs(a[2] / b[2]))
        if b[2] != 0 and a[3] > 0:
            rel.append(a[3] / abs(b[2]))
    if not sig:
        continue
    rel.sort()
    med = rel[len(rel) // 2] if rel else 0
    print('  %-34s %6.2f %10.2e %10.2e %10.2e' % (nm, max(sig), max(eff), med, med * math.sqrt(n / max(nneed, 1))))
