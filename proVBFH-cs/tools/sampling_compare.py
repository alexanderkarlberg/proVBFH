#!/usr/bin/env python3
"""Compare sampling variants of the same exclusive-part run.

Usage: sampling_compare.py REF VAR [VAR ...] [--hist H ...] [--range H:LO:HI ...]
Each argument is a glob of seed directories (with run.log and pwg*-EXCL*.top).
Per bin: mean over seeds, error from the seed scatter, and the CPU from the
"exclusive part: CPU" lines. The efficiency gain of VAR over REF in a bin is
  (err_ref^2 cpu_ref) / (err_var^2 cpu_var),
i.e. how much less CPU VAR needs for the same error. Printed: the gain for
the selected histograms' bins (default: the fiducial totals) and its median
over all bins with a nonzero error in both, per histogram and overall, and
chi2/n of VAR against REF (consistency). The gain is also given with the
VEGAS errors (sqrt(sum err^2)/N; more stable with few seeds, but blind to
rare large weights). --range H:LO:HI integrates histogram H over [LO, HI]
per seed and compares the integrals (seed-scatter errors), which is more
robust than single tail bins.
"""
import argparse
import glob
import math
import os
import re
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from combine_parts import read_top  # noqa: E402


def load(pattern):
    dirs = sorted(d for d in glob.glob(pattern) if glob.glob(d + '/pwg*-EXCL*.top'))
    cpu, data, order, seeds = 0.0, [], None, []
    for d in dirs:
        m = re.findall(r'exclusive part: CPU\s+([\d.]+)', open(d + '/run.log').read())
        if not m:
            continue
        cpu += float(m[-1])
        h, o = read_top(glob.glob(d + '/pwg*-EXCL*.top')[0])
        data.append(h); order = order or o
    n = len(data)
    out = {}
    for h in order:
        rows = []
        for i in range(len(data[0][h])):
            v = [x[h][i][2] for x in data]
            mu = sum(v) / n
            e = math.sqrt(sum((y - mu) ** 2 for y in v) / (n - 1) / n) if n > 1 else 0
            ev = math.sqrt(sum(x[h][i][3] ** 2 for x in data)) / n
            rows.append((data[0][h][i][0], data[0][h][i][1], mu, e, ev))
        out[h] = rows
    return out, order, cpu, n, data


def integral(data, h, lo, hi):
    """per-seed integral of histogram h over [lo, hi]: mean and scatter error"""
    v = [sum(b[2] * (b[1] - b[0]) for b in x[h] if b[0] >= lo and b[1] <= hi) for x in data]
    n = len(v); mu = sum(v) / n
    return mu, math.sqrt(sum((y - mu) ** 2 for y in v) / (n - 1) / n)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('runs', nargs='+')
    ap.add_argument('--hist', action='append')
    ap.add_argument('--range', action='append', default=[])
    a = ap.parse_args()
    ref, order, cref, nref, dref = load(a.runs[0])
    print(f'REF {a.runs[0]}: {nref} seeds, {cref/3600:.1f} CPU-h')
    sel = a.hist or [h for h in order if h.startswith('sig')]
    for pat in a.runs[1:]:
        var, _, cvar, nvar, dvar = load(pat)
        print(f'\nVAR {pat}: {nvar} seeds, {cvar/3600:.1f} CPU-h')
        allg, allv, c2, nc = [], [], 0, 0
        per = []
        for h in order:
            gs = []
            for r, v in zip(ref[h], var[h]):
                if r[3] > 0 and v[3] > 0:
                    gs.append((r[3] ** 2 * cref) / (v[3] ** 2 * cvar))
                    if r[4] > 0 and v[4] > 0:
                        allv.append((r[4] ** 2 * cref) / (v[4] ** 2 * cvar))
                    c2 += (r[2] - v[2]) ** 2 / (r[3] ** 2 + v[3] ** 2); nc += 1
            if gs:
                allg += gs; per.append((h, statistics.median(gs), len(gs)))
        for h in sel:
            for r, v in zip(ref[h], var[h]):
                if r[3] > 0 and v[3] > 0:
                    print(f'  {h:34s} [{r[0]:g},{r[1]:g}]: REF {r[2]: .5e} +- {r[3]:.1e}  VAR {v[2]: .5e} +- {v[3]:.1e}'
                          f'  gain {(r[3]**2*cref)/(v[3]**2*cvar):6.2f}')
        print(f'  gain, median over {len(allg)} bins: {statistics.median(allg):.2f}; '
              f'16-84%: {sorted(allg)[len(allg)*16//100]:.2f}-{sorted(allg)[len(allg)*84//100]:.2f}; '
              f'chi2/n VAR vs REF {c2:.1f}/{nc}')
        print(f'  gain from the VEGAS errors, median over {len(allv)} bins: {statistics.median(allv):.2f}')
        for rg in a.range:
            h, lo, hi = rg.rsplit(':', 2); lo, hi = float(lo), float(hi)
            mr, er = integral(dref, h, lo, hi); mv, ev = integral(dvar, h, lo, hi)
            g = (er ** 2 * cref) / (ev ** 2 * cvar) if ev > 0 else float('inf')
            print(f'  {h} in [{lo:g},{hi:g}]: REF {mr: .4e} +- {er:.1e}  VAR {mv: .4e} +- {ev:.1e}'
                  f'  pull {(mv-mr)/math.hypot(er, ev):+.1f}  gain {g:.2f}')
        per.sort(key=lambda t: t[1])
        print('  lowest median gains:  ' + '; '.join(f'{h} {g:.2f}' for h, g, _ in per[:4]))
        print('  highest median gains: ' + '; '.join(f'{h} {g:.2f}' for h, g, _ in per[-4:]))


if __name__ == '__main__':
    main()
