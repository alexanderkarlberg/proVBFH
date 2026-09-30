#!/usr/bin/env python3
"""Compare proVBFH-cs runs of the exclusive part at different cutoffs.

Usage: cutoff_compare.py DIR [DIR ...]
Each DIR holds seed directories s*/ with pwg*-EXCL*.top. Per DIR and
histogram bin: the mean over seeds, the error from the scatter of the
seeds (std/sqrt(N)) and the mean of the VEGAS errors/sqrt(N). Prints the
VBF-cut cross section and, for each pair of DIRs, chi2/n and the largest
pull per histogram (scatter errors). Bins that are zero up to rounding in
both DIRs (|value| < ZERO times the largest |value| of the run, e.g. the
Higgs-only histograms, for which the exclusive part vanishes pointwise)
are left out: their pulls measure only rounding noise.
"""
import glob
import math
import re
import sys

sys.path.insert(0, __file__.rsplit('/', 1)[0])
from combine_parts import read_top  # noqa: E402

ZERO = 1e-12


def load(d):
    files = sorted(glob.glob(f'{d}/s*/pwg*-EXCL*.top'))
    data = [read_top(f)[0] for f in files]
    order = read_top(files[0])[1]
    out = {}
    for h in order:
        rows = [x[h] for x in data]
        n = len(rows)
        res = []
        for i in range(len(rows[0])):
            v = [r[i][2] for r in rows]
            m = sum(v) / n
            sd = math.sqrt(sum((x - m) ** 2 for x in v) / (n - 1)) if n > 1 else 0
            ev = math.sqrt(sum(r[i][3] ** 2 for r in rows)) / n
            res.append((rows[0][i][0], rows[0][i][1], m, sd / math.sqrt(n), ev))
        out[h] = res
    return out, order, len(files)


def main():
    dirs = sys.argv[1:]
    runs = [load(d) for d in dirs]
    order = runs[0][1]
    sig = 'sig(all VBF cuts 2 jets)'
    print('sigma(VBF cuts) of the exclusive part [pb]:')
    for d, (h, _, n) in zip(dirs, runs):
        b = h[sig][0]
        print(f'  {d:12s} {n:3d} seeds: {b[2]: .5f} +- {b[3]:.5f} (seed scatter), VEGAS {b[4]:.5f}')
    scale = max(abs(b[2]) for h, _, _ in runs for r in h.values() for b in r)
    for a in range(len(dirs)):
        for c in range(a + 1, len(dirs)):
            print(f'\n{dirs[a]} vs {dirs[c]} (scatter errors):')
            tot, ntot = 0, 0
            for hname in order:
                ha, hc = runs[a][0][hname], runs[c][0][hname]
                chi, n, pmax = 0, 0, 0
                for x, y in zip(ha, hc):
                    e2 = x[3] ** 2 + y[3] ** 2
                    if e2 == 0 or max(abs(x[2]), abs(y[2])) < ZERO * scale:
                        continue
                    p = (x[2] - y[2]) / math.sqrt(e2)
                    chi += p * p; n += 1; pmax = max(pmax, abs(p))
                if n:
                    tot += chi; ntot += n
                    print(f'  {hname:28s} chi2/n = {chi:7.1f}/{n:3d}   max|pull| {pmax:5.2f}')
            print(f'  {"all":28s} chi2/n = {tot:7.1f}/{ntot:3d}')


if __name__ == '__main__':
    main()
