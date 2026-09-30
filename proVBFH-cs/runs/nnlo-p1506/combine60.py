#!/usr/bin/env python3
"""The 2-jet total from the 60 exclusive seeds (nnlo-full with the current
analysis, nnlo-p1506 with the paper's: same bin name, same cuts) plus the
inclusive part, and the 3-/4-jet rates from nnlo-p1506 alone, against the
1506.02660 files (runs/ref-1506.02660/11.top)."""
import glob, math, os, sys
sys.path.insert(0, os.path.expanduser('~/cernbox/proVBFH-github/proVBFH-cs/tools'))
from combine_parts import read_top, average
R = os.path.expanduser('~/cernbox/proVBFH-github/proVBFH-cs/runs')
ref, _ = read_top(R + '/ref-1506.02660/11.top')
inc, _ = average(sorted(glob.glob(R + '/nnlo-incl-p1506/incl-s*/pwg*-LO*.top')), 'max')
def seeds(g, h):
    return [read_top(f)[0][h][0][2] for f in sorted(glob.glob(g))]
for h in ('sig(all VBF cuts 2 jets)', 'sig(all VBF cuts 3 jets)', 'sig(all VBF cuts 4 jets)'):
    a = seeds(R + '/nnlo-p1506/s*/pwg*-EXCL*.top', h)
    b = seeds(R + '/nnlo-full/s*/pwg*-EXCL*.top', h) if h.endswith('2 jets)') else []
    for lab, v in (('p1506 (30)', a), ('nnlo-full (30)', b), ('all 60', a + b)):
        if not v: continue
        n = len(v); m = sum(v) / n; e = math.sqrt(sum((x - m) ** 2 for x in v) / (n - 1) / n)
        ib = inc[h][0]; t = m + ib[2]; te = math.hypot(e, ib[3]); rb = ref[h][0]
        print(f'{h:26s} {lab:15s}: exclusive {m: .5f} +- {e:.5f}, total {t:.5f} +- {te:.5f}; '
              f'paper {rb[2]:.5f} +- {rb[3]:.5f}; diff {t-rb[2]:+.5f} ({(t-rb[2])/math.hypot(te, rb[3]):+.1f} sigma, {100*(t/rb[2]-1):+.2f}%)')
