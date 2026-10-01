#!/usr/bin/env python3
"""Fixed-scale (mu_R = mu_F = m_H) check at the 1506.02660 set-up: the 2-,
3- and 4-jet rates of proVBFH-cs (runs/fixmh-{nlo,nnlo,incl}-p1506), the
old proVBFH (runs/fixmh-old1506{-nlo,}: plain mean, median and the
combine_runs.f trimming) and VBFNLO 3.0 process 110 (>= 3 jets only; LO =
H+3j tree, NLO = the O(alpha_s^2) 3-jet rate). With --dyn, the same for
the dynamic scale mu_0(ptH): nlo-p1506, nnlo-p1506 + nnlo-p1506-h03, the
paper files and VBFNLO with scale ID 20."""
import glob, math, os, re, sys
import numpy as np
sys.path.insert(0, os.path.expanduser('~/cernbox/proVBFH-github/proVBFH-cs/tools'))
sys.path.insert(0, os.path.expanduser('~/cernbox/proVBFH-github/notes/2026-09-30-hxswg-comparison/tools'))
from combine_parts import read_top, average
from rawload import combine_runs
R = os.path.expanduser('~/cernbox/proVBFH-github/proVBFH-cs/runs')
V = os.path.expanduser('~/work/disorder-comparisons/vbfnlo/runs')
dyn = '--dyn' in sys.argv
H = ['sig(all VBF cuts 2 jets)', 'sig(all VBF cuts 3 jets)', 'sig(all VBF cuts 4 jets)']

def seeds(globs, h):
    out = []
    for g in globs:
        for f in sorted(glob.glob(g)):
            t = read_top(f)[0]
            if h in t: out.append((t[h][0][2], t[h][0][3]))
    return np.array(out).reshape(-1, 2)

def mean(v):
    n = len(v)
    if n < 2: return float('nan'), float('nan'), n
    return v.mean(), v.std(ddof=1) / math.sqrt(n), n

def fmt(m, e): return f'{m:.5f} +- {e:.5f}'

def vbfnlo(d):
    lo, nlo = [], []
    for f in sorted(glob.glob(d + '/j*/run.log')):
        s = open(f).read()
        a = re.search(r'TOTAL result \(LO\):\s+(\S+)\s+\+-\s+(\S+)', s)
        b = re.search(r'TOTAL result \(NLO\):\s+(\S+)\s+\+-\s+(\S+)', s)
        if a and b:
            lo.append((float(a.group(1)), float(a.group(2)))); nlo.append((float(b.group(1)), float(b.group(2))))
    return np.array(lo).reshape(-1, 2) / 1000, np.array(nlo).reshape(-1, 2) / 1000   # fb -> pb

def wmean(x):
    if len(x) == 0: return float('nan'), float('nan'), 0
    w = 1 / x[:, 1] ** 2; m = (w * x[:, 0]).sum() / w.sum()
    return m, max(1 / math.sqrt(w.sum()), x[:, 0].std(ddof=1) / math.sqrt(len(x)) if len(x) > 1 else 0), len(x)

if dyn:
    our_nlo = [R + '/nlo-p1506/s*/pwg*-EXCL*.top']; our_nnlo = [R + '/nnlo-p1506/s*/pwg*-EXCL*.top', R + '/nnlo-p1506-h03/s*/pwg*-EXCL*.top']
    inc_nlo = R + '/nlo-incl-p1506/incl-s*/pwg*-LO*.top'; inc_nnlo = R + '/nnlo-incl-p1506/incl-s*/pwg*-LO*.top'
    old_nlo = R + '/old1506-nlo/pwg-*-NLO.top'; old_nnlo = None; vdir = V + '/dyn-1506'
else:
    our_nlo = [R + '/fixmh-nlo-p1506/s*/pwg*-EXCL*.top']; our_nnlo = [R + '/fixmh-nnlo-p1506/s*/pwg*-EXCL*.top']
    inc_nlo = None; inc_nnlo = R + '/fixmh-incl-p1506/incl-s*/pwg*-LO*.top'
    old_nlo = R + '/fixmh-old1506-nlo/pwg-*-NLO.top'; old_nnlo = R + '/fixmh-old1506/pwg-*-NNLO.top'; vdir = V + '/fixmh-1506'

print('scale:', 'mu_0(ptH)' if dyn else 'fixed m_H', ' [pb]')
for order, ours, inc, old in (('NLO', our_nlo, inc_nlo, old_nlo), ('NNLO', our_nnlo, inc_nnlo, old_nnlo)):
    for h in H:
        line = f'{order:4s} {h:26s}'
        m, e, n = mean(seeds(ours, h)[:, 0]) if len(seeds(ours, h)) else (float('nan'),) * 2 + (0,)
        if h.endswith('2 jets)'):
            if inc and glob.glob(inc):
                ib = average(sorted(glob.glob(inc)), 'max')[0][h][0]; m, e = m + ib[2], math.hypot(e, ib[3])
            else:
                m = e = float('nan')
        line += f' ours({n}) {fmt(m, e)}'
        if old and glob.glob(old):
            s = seeds([old], h)
            if len(s) > 2:
                pm, pe, k = mean(s[:, 0]); med = np.median(s[:, 0])
                hw = (np.percentile(s[:, 0], 84) - np.percentile(s[:, 0], 16)) / 2
                cr, ce, frac = combine_runs(s[:, 0], s[:, 1])
                line += f' | old({k}) plain {fmt(pm, pe)}, median {med:.5f} +- {1.2533*hw/math.sqrt(k):.5f}, combine_runs {fmt(cr, ce)} (kept {frac:.3f})'
        if dyn and order == 'NNLO':
            rb = read_top(R + '/ref-1506.02660/11.top')[0][h][0]
            line += f' | paper {fmt(rb[2], rb[3])}'
        print(line)
lo, nlo = vbfnlo(vdir)
for lab, x in (('LO (H+3j tree)', lo), ('NLO (O(as^2) >= 3 jets)', nlo)):
    m, e, n = wmean(x)
    print(f'VBFNLO {lab:26s} ({n} jobs): {fmt(m, e)}')
