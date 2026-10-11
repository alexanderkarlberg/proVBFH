#!/usr/bin/env python3
"""Final statistics of the small-x channel 'on' arm (b-on pilot seeds 1302001-150 + c-on seeds 1303001-600) vs VBFNLO fixmh (800 seeds).
Plain means, seed-scatter errors. Prints sigma(>=3j), (>=4j), exactly 3j, chi2 for ptj3, yj3, ptj4."""
import glob, math, os, sys
sys.path.insert(0, '/u/akarlber/work/proVBFH/proVBFH-cs/production/h3j-crosscheck/chsplit')
import chsum
R = '/ptmp/mpp/akarlber/smallx'; VB = '/ptmp/mpp/akarlber/h3j/vbfnlo/prod/nlo-fixmh'
HN = ('ptj3', 'yj3', 'ptj4')
def hists(path):
    o, name = {}, None
    for l in open(path):
        if l.startswith('#'):
            t = l[1:].split(' index')[0].strip(); name = t if t in HN else None
            if name: o[name] = []
        elif name and l.strip():
            p = l.split()
            if len(p) == 4: o[name].append((float(p[0]), float(p[1]), 1000*float(p[2].replace('D','E'))))
    return o
def cs_jobs():
    out = {}
    for s in ('b-on', 'c-on'):
        for d in glob.glob(f'{R}/{s}/job-*'):
            if os.path.exists(d+'/done') and os.path.exists(d+'/pwg-EXCL.top'): out[int(d.split('-')[-1])] = d+'/pwg-EXCL.top'
    return out
def vb_jobs():
    return {int(d.split('-')[-1]): d+'/p1506_nlo.top' for d in glob.glob(VB+'/job-*') if os.path.exists(d+'/done') and os.path.exists(d+'/p1506_nlo.top')}
def mom(v):
    n = len(v); m = sum(v)/n; sd = math.sqrt(sum((x-m)**2 for x in v)/(n-1)); s = math.sqrt(sum((x-m)**2 for x in v)/n)
    return m, sd, sum((x-m)**3 for x in v)/n/s**3 if s else 0
def load(jobs):
    D = {}
    for s, f in jobs.items():
        try: D[s] = chsum.read_sig(f)
        except Exception: pass
    return D
C, V = cs_jobs(), vb_jobs(); cd, vd = load(C), load(V)
res = {}
for lab, D in (('cs', cd), ('vbfnlo', vd)):
    seeds = sorted(D)
    for nm, f in (('>=3j', lambda x: x[0]), ('>=4j', lambda x: x[1]), ('=3j', lambda x: x[0]-x[1])):
        v = [f(D[s]) for s in seeds]; n = len(v); m, sd, sk = mom(v); e = sd/math.sqrt(n)
        res[(lab, nm)] = (m, e)
        top = sorted(zip(v, seeds), reverse=True)[:4]
        extra = ''
        if lab == 'cs' and nm != '=3j' or nm == '=3j' and lab == 'cs':
            extra = ' | max ' + ' '.join(f'{x:.2f}({s})' for x, s in top) + ' | running ' + ' '.join(f'{k}:{sum(v[:k])/k:.2f}' for k in (50, 100, 200, 400, 600, n) if k <= n)
        print(f'{lab:6s} {nm:5s} n={n} mean {m:.3f}+-{e:.3f} sd {sd:.2f} skew {sk:.2f}{extra}')
for nm in ('>=3j', '>=4j', '=3j'):
    (a, ea), (b, eb) = res[('cs', nm)], res[('vbfnlo', nm)]
    d = a-b; e = math.hypot(ea, eb); print(f'diff cs-VBFNLO {nm}: {d:+.3f} +- {e:.3f} fb, pull {d/e:+.2f}, rel {d/b*100:+.2f}%')
# cs-only error on the full-stat 'on' set when subsetting is not needed; distributions
Hc = {}; Hv = {}
for lab, J, H in (('cs', C, Hc), ('vb', V, Hv)):
    for s, f in J.items():
        try: h = hists(f)
        except Exception: continue
        for k in HN:
            if k in h: H.setdefault(k, []).append([b[2] for b in h[k]]); H.setdefault(k+'_e', h[k])
for k in HN:
    A, B = Hc[k], Hv[k]; edges = Hc[k+'_e']; nb = len(edges); chi = 0; nd = 0; lines = []
    for i in range(nb):
        a = [x[i] for x in A]; b = [x[i] for x in B]
        ma, sa, _ = mom(a) if len(a) > 2 else (0, 0, 0); mb, sb, _ = mom(b)
        ea, eb = sa/math.sqrt(len(a)), sb/math.sqrt(len(b)); e = math.hypot(ea, eb)
        if e <= 0: continue
        p = (ma-mb)/e; chi += p*p; nd += 1; lines.append(f'  [{edges[i][0]:.3g},{edges[i][1]:.3g}] cs {ma:.4g} vb {mb:.4g} pull {p:+.1f}')
    print(f'chi2 {k}: {chi:.1f} / {nd} bins; sum of pulls {sum(float(l.split("pull")[1]) for l in lines):+.1f}')
    for l in lines:
        if abs(float(l.split('pull')[1])) >= 2: print(l)
