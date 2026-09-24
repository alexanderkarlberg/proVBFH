#!/usr/bin/env python3
# Table 2 of 2005.11334 (VBF HH non-factorisable corrections) from the
# HH-born-*, HH-NF1-*, HH-NF2-* run directories made by nf-runs.sh:
# sig(all VBF cuts 2 jets) of pwg-LO-0001.top, TB = ALL - TT - BB.
import sys, math, os
d = sys.argv[1] if len(sys.argv) > 1 else 'runs/nonfact'
def sig(run):
    lines = open(os.path.join(d, run, 'pwg-LO-0001.top')).read().splitlines()
    for i, l in enumerate(lines):
        if 'all VBF cuts 2 jets' in l:
            v = lines[i+1].split()
            return 1e3*float(v[2]), 1e3*float(v[3])   # pb -> fb
    raise RuntimeError(run)
def comb(a, b, c):   # a - b - c with errors added in quadrature
    return a[0]-b[0]-c[0], math.sqrt(a[1]**2+b[1]**2+c[1]**2)
rows = {}
for kind in ['born', 'NF1', 'NF2']:
    t, b, a = (sig(f'HH-{kind}-{s}') for s in ('TT', 'BB', 'ALL'))
    rows[kind] = {'TT': t, 'BB': b, 'TB': comb(a, t, b), 'Sum': a}
paper = {'born': (10.393, 14.172, -23.904, 0.662), 'NF1': (0.339, 0.300, 0.318, 0.286),
         'NF2': (-0.667, -0.621, -0.644, -0.516), 'full': (-0.327, -0.320, -0.326, -0.230)}
cols = ['TT', 'BB', 'TB', 'Sum']
print('| | ' + ' | '.join(cols) + ' |')
print('|---|' + '---|'*4)
print('| Born [fb] | ' + ' | '.join(f'{rows["born"][c][0]:.3f} ({p})' for c, p in zip(cols, paper['born'])) + ' |')
for kind, lab in [('NF1', '1-loop NF'), ('NF2', '2-loop NF')]:
    out = []
    for c, p in zip(cols, paper[kind]):
        v, e = rows[kind][c]; B = rows['born'][c][0]
        out.append(f'{100*v/B:.3f}% ± {100*e/abs(B):.3f} ({p}%)')
    print(f'| {lab} | ' + ' | '.join(out) + ' |')
out = []
for c, p in zip(cols, paper['full']):
    v = rows['NF1'][c][0] + rows['NF2'][c][0]; e = math.hypot(rows['NF1'][c][1], rows['NF2'][c][1]); B = rows['born'][c][0]
    out.append(f'{100*v/B:.3f}% ± {100*e/abs(B):.3f} ({p}%)')
print('| full NF | ' + ' | '.join(out) + ' |')
if os.path.exists(os.path.join(d, 'HH-NF-ALL')):
    v, e = sig('HH-NF-ALL'); B = rows['born']['Sum'][0]
    print(f'\ndirect full NF, all diagrams: {v:.5f} ± {e:.5f} fb = {100*v/B:.3f}% ± {100*e/B:.3f} (paper -0.0015 fb, -0.230%)')
