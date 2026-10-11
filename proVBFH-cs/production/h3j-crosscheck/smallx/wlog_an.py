#!/usr/bin/env python3
"""Variance attribution of the stage-2 per-point log (cs_wlog.dat, scratch diagnostic of 11 Oct).
Columns: idx c vw x1 x2 Q1 Q2 xrand(1:13) [ok xp z wrad]x2 dw4(2) okr1 passr1 okr2 passr2.
Every point with c > thr is logged, plus every 97th point (unbiased subsample)."""
import sys, glob, numpy as np
def load(pat):
    rows = []
    for f in sorted(glob.glob(pat)):
        for l in open(f):
            t = l.replace(' T', ' 1').replace(' F', ' 0').split()
            if len(t) == 34: rows.append([float(x) for x in t])
    return np.array(rows)
D = load(sys.argv[1]); thr = float(sys.argv[2]) if len(sys.argv) > 2 else 0.002
idx, c, vw, x1, x2, Q1, Q2 = D[:,0], D[:,1], D[:,2], D[:,3], D[:,4], D[:,5], D[:,6]
sub = (idx.astype(np.int64) % 97 == 0)
big = c > thr
print('rows', len(D), 'subsample', sub.sum(), 'big', big.sum())
# unbiased estimator of sum c and sum c^2 : subsample*97 (all points) -- big points are in the subsample w.p. 1/97 only,
# so use: S = 97*sum_sub(c, c<=thr) + sum_big(c)
def est(p):
    return 97*np.sum((c**p)[sub & ~big]) + np.sum((c**p)[big])
S1, S2 = est(1), est(2)
print('sum c %.4g  sum c^2 %.4g  (big share of c: %.3f, of c^2: %.3f)' % (S1, S2, np.sum(c[big])/S1, np.sum(c[big]**2)/S2))
# regions
xs = np.minimum(x1, x2); line = np.where(x1 < x2, 1, 2)
def share(mask, name):
    m = lambda p: (97*np.sum((c**p)[sub & ~big & mask]) + np.sum((c**p)[big & mask]))
    n = 97*np.sum(sub & mask)  # points sampled (estimate)
    print('%-26s sampled frac %.4f  share c %.3f  share c^2 %.3f  (big n=%d)' % (name, n/ (97*sub.sum()), m(1)/S1, m(2)/S2, np.sum(big & mask)))
for lo, hi in ((0, .002), (.002, .005), (.005, .012), (.012, .02), (.02, .05), (.05, .1), (.1, .3), (.3, 1)):
    share((xs >= lo) & (xs < hi), 'x_small in [%g,%g)' % (lo, hi))
share((xs < .012) & (np.maximum(x1, x2) > .1), 'x_small<.012, x_big>.1')
for ln in (1, 2): share((xs < .012) & (line == ln), 'x_small<.012, beam %d small' % ln)
print('big points: x_small median %.4f  x_big median %.3f  Q_small-x-line / other median %.0f / %.0f' % (np.median(xs[big]), np.median(np.maximum(x1,x2)[big]),
      np.median(np.where(line==1, Q1, Q2)[big]), np.median(np.where(line==1, Q2, Q1)[big])))
print('subsample: x_small median %.4f' % np.median(xs[sub]))
