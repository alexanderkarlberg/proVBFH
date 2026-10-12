#!/usr/bin/env python3
"""Bin-by-bin sigma(>=3j) NLO: proVBFH-cs per-piece production (lo + v + rd, plain means over jobs of
pwg-EXCL.top) against VBFNLO's 800 fixmh jobs (p1506_nlo.top, last-iteration estimator). Errors = scatter
over jobs / sqrt(N). Prints, per 3-jet observable, the bins and the summary chi2; details to bins_vs_vbfnlo.txt."""
import glob, math, statistics as st, sys
OBS = ['sig(all VBF cuts 3 jets)', 'ptj3', 'yj3', 'y*j3', 'min{rap(j1,j3),rap(j3,j2)}']
def read(f):
    h = {}; cur = None
    for l in open(f):
        if l.startswith('#'):
            cur = l[1:].rsplit('index', 1)[0].strip(); h[cur] = []
        elif l.strip() and cur is not None:
            x = l.split(); h[cur].append((float(x[0]), float(x[1]), float(x[2])))
    return h
def mean(files):
    hs = [read(f) for f in files]; out = {}
    for o in OBS:
        nb = len(hs[0][o]); out[o] = []
        for b in range(nb):
            v = [h[o][b][2] for h in hs]
            out[o].append((hs[0][o][b][0], hs[0][o][b][1], st.mean(v), st.stdev(v) / math.sqrt(len(v))))
    return out, len(hs)
B = '/ptmp/mpp/akarlber/pieces/prod'
cs = {}
for p in ['lo', 'v', 'rd']:
    cs[p], n = mean(sorted(glob.glob(B + '/%s/job-*/pwg-EXCL.top' % p)))
vb, nv = mean(sorted(glob.glob('/ptmp/mpp/akarlber/h3j/vbfnlo/prod/nlo-fixmh/job-*/p1506_nlo.top')))
det = open('bins_vs_vbfnlo.txt', 'w')
for o in OBS:
    chi2 = 0; nb = 0; rows = []
    for b, (lo, hi, v, ev) in enumerate(vb[o]):
        c = sum(cs[p][o][b][2] for p in cs); ec = math.sqrt(sum(cs[p][o][b][3] ** 2 for p in cs))
        w = hi - lo if o.startswith(('pt', 'y', 'min')) else 1
        d = c - v; e = math.hypot(ec, ev)
        if v == 0 and c == 0: continue
        pull = d / e if e > 0 else 0; chi2 += pull ** 2; nb += 1
        rows.append('%9.3g %9.3g  cs %10.4g  vb %10.4g  diff*width %8.3f +- %6.3f fb  rel %+6.2f%%  pull %+5.1f'
                    % (lo, hi, c, v, d * w, e * w, 100 * d / v if v else 0, pull))
    det.write('# %s\n' % o + '\n'.join(rows) + '\n\n')
    print('%-30s bins %2d  chi2 %6.1f' % (o, nb, chi2))
    if o in ('ptj3', 'yj3'): print('\n'.join(rows))
