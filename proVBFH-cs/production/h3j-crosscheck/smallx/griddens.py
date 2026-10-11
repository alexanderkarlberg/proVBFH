#!/usr/bin/env python3
"""Per-dimension VEGAS densities g_j(xrand_j) of the logged points, from grids-excl.top of the same job.
Usage: griddens.py 'diag/job-*' [thr] [ntop]"""
import sys, glob, numpy as np
def read_top(f):
    g = {}; cur = None
    for l in open(f):
        t = l.split()
        if len(t) == 4 and t[1] == 'XDIM': cur = int(t[3]); g[cur] = []
        elif len(t) == 2 and cur is not None:
            try: g[cur].append(float(t[1].replace('D', 'E')))
            except ValueError: pass
    return {k: np.array(v) for k, v in g.items()}
def dens(edges, x):
    e = np.concatenate([[0.], edges]); i = np.clip(np.searchsorted(e, x, side='left') - 1, 0, len(e) - 2)
    return 1/(len(edges)*(e[i+1] - e[i]))
if __name__ == '__main__':
    thr = float(sys.argv[2]) if len(sys.argv) > 2 else 0.002; ntop = int(sys.argv[3]) if len(sys.argv) > 3 else 15
    out = []
    for d in sorted(glob.glob(sys.argv[1])):
        try: G = read_top(d + '/grids-excl.top')
        except OSError: continue
        for l in open(d + '/cs_wlog.dat'):
            t = l.replace(' T', ' 1').replace(' F', ' 0').split()
            if len(t) != 34 or float(t[1]) <= thr: continue
            r = np.array([float(x) for x in t]); xr = r[7:20]
            g = np.array([dens(G[j], xr[j]) for j in range(13)])
            out.append((r, g))
    out.sort(key=lambda a: -a[0][1])
    print('n big', len(out), ' dims 1..13: Born(1 mH, 2 tau, 3 y, 4 mjj, 5 cthV, 6 cthj, 7 phi) line1(8-10) line2(11-13)')
    for r, g in out[:ntop]:
        print('c=%.3f vw=%.1e prod(1/g)/N=%.1e x=%.3f,%.3f | g:' % (r[1], r[2], np.prod(1/g)/(3*1.2e6), r[3], r[4]), ' '.join('%.2g' % x for x in g))
    G_ = np.array([g for _, g in out]); c = np.array([r[1] for r, _ in out])
    print('c^2-weighted mean log10 g per dim:', ' '.join('%.2f' % (np.sum(c**2*np.log10(G_[:, j]))/np.sum(c**2)) for j in range(13)))
    print('unweighted mean log10 g (big points) :', ' '.join('%.2f' % np.mean(np.log10(G_[:, j])) for j in range(13)))
