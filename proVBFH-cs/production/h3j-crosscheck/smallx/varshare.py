#!/usr/bin/env python3
"""For several variables: share of sum c^2 (all points, estimated from big points + 97x subsample) per bin / sampled fraction per bin.
A ratio >> 1 flags a region where the variance sits but few points are sampled.  Usage: varshare.py 'glob' [thr]"""
import sys, io, contextlib, runpy, numpy as np
sys.argv = ['x', sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else '0.002']
with contextlib.redirect_stdout(io.StringIO()): g = runpy.run_path('wlog_an.py')
D, c, sub, big = g['D'], g['c'], g['sub'], g['big']
x1, x2, Q1, Q2 = D[:,3], D[:,4], D[:,5], D[:,6]
w2 = np.where(sub & ~big, 97*c**2, np.where(big, c**2, 0.)); wn = np.where(sub, 97., 0.)
def tab(name, v, edges):
    print(name)
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (v >= lo) & (v < hi)
        f = wn[m].sum()/wn.sum(); s = w2[m].sum()/w2.sum()
        print('   [%7.3g,%7.3g) sampled %.4f  c2 share %.4f  ratio %5.2f  (big n=%d)' % (lo, hi, f, s, s/f if f > 0 else np.nan, (big & m).sum()))
xs = np.minimum(x1, x2)
tab('x_small', xs, [0, .005, .012, .02, .035, .06, .1, .2, .4, 1])
tab('tau=x1x2', x1*x2, [0, 3e-4, 1e-3, 2e-3, 3e-3, 5e-3, 1e-2, 3e-2, .1, 1])
tab('cthj-sampling r6 (xrand(6))', D[:,12], [0, .1, .2, .4, .6, .8, .9, 1.0001])
Qs = np.where(x1 < x2, Q1, Q2); Qo = np.where(x1 < x2, Q2, Q1)
tab('Q of small-x line', Qs, [0, 20, 40, 80, 120, 200, 400, 1e4])
tab('Q of other line', Qo, [0, 20, 40, 80, 120, 200, 400, 1e4])
xp_s = np.where(x1 < x2, D[:,21], D[:,25]); wr_s = np.where(x1 < x2, D[:,23], D[:,27])
tab('xp of small-x line radiation', xp_s, [0, .02, .05, .1, .2, .4, .7, .9, .99, 1.0001])
