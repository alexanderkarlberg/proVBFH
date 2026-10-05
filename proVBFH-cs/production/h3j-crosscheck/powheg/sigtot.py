#!/usr/bin/env python3
"""Combine the 'sig(...)' one-bin histograms of pwg-*-NLO.top over seeds.
usage: sigtot.py <rundir>   -> per-seed values, mean, error (seed scatter and
quadrature of the per-seed errors)."""
import sys, glob, math, re
files = sorted(glob.glob(sys.argv[1] + '/pwg-[0-9][0-9][0-9][0-9]-NLO.top'))
names = ['sig incl cuts', 'sig(all VBF cuts 2 jets)', 'sig(all VBF cuts 3 jets)', 'sig(all VBF cuts 4 jets)']
data = {n: [] for n in names}
for f in files:
    lines = open(f).read().splitlines()
    for i, l in enumerate(lines):
        for n in names:
            if l.startswith('# ' + n + ' index'):
                v, e = map(float, lines[i + 1].split()[2:4])
                data[n].append((v, e))
print(f'{len(files)} seeds in {sys.argv[1]}')
for n in names:
    d = data[n]
    if not d: continue
    k = len(d); m = sum(v for v, _ in d) / k
    eq = math.sqrt(sum(e * e for _, e in d)) / k
    es = math.sqrt(sum((v - m) ** 2 for v, _ in d) / (k - 1) / k) if k > 1 else float('nan')
    print(f'{n:28s} mean {m*1e3:10.4f} fb  +- {eq*1e3:.4f} (stat, quad)  +- {es*1e3:.4f} (seed scatter)')
    if '3 jets' in n:
        print('   per seed [fb]:', ' '.join(f'{v*1e3:.2f}' for v, _ in d))
