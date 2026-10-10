#!/usr/bin/env python3
"""Low-slice statistics: per-job sigma(>=3j,>=4j) for lo1, lo2, lo=lo1+lo2 (same seed) for cs and vb (pilot + stat jobs)."""
import glob, math, os, sys
sys.path.insert(0, '/u/akarlber/work/proVBFH/proVBFH-cs/production/h3j-crosscheck/chsplit')
import chsum
R = '/ptmp/mpp/akarlber/xslice'
def load(code, sl):
    f = 'pwg-EXCL.top' if code == 'cs' else 'p1506_nlo.top'
    o = {}
    for d in glob.glob(f'{R}/{code}-{sl}/job-*'):
        if os.path.exists(d+'/done') and os.path.exists(f'{d}/{f}'):
            try: r = chsum.read_sig(f'{d}/{f}')
            except Exception: continue
            if len(r) == 2: o[int(d.split('-')[-1])] = list(r)
    return o
def mom(v):
    n = len(v); m = sum(v)/n; sd = math.sqrt(sum((x-m)**2 for x in v)/(n-1)) if n > 1 else 0
    s = math.sqrt(sum((x-m)**2 for x in v)/n)
    return m, sd, sum((x-m)**3 for x in v)/n/s**3 if s else 0
def med(v): s = sorted(v); n = len(s); return s[n//2] if n % 2 else (s[n//2-1]+s[n//2])/2
for code in ('cs', 'vb'):
    A, B = load(code, 'lo1'), load(code, 'lo2')
    seeds = sorted(set(A) & set(B))
    sets = {'lo1': A, 'lo2': B, 'lo': {s: [A[s][i]+B[s][i] for i in (0, 1)] for s in seeds}}
    for nj, i in ((3, 0), (4, 1)):
        for k, D in sets.items():
            sd_ = sorted(D); v = [D[s][i] for s in sd_]; n = len(v)
            m, sd, sk = mom(v); e = sd/math.sqrt(n)
            top = sorted(zip(v, sd_), reverse=True)
            tr = lambda k_: sum(x for x, _ in top[k_:])/(n-k_)
            print(f'{code} {k:3s} >={nj}j n={n:3d} mean {m:.4f}+-{e:.4f} med {med(v):.4f} sd {sd:.4f} skew {sk:5.2f} | max '
                  + ' '.join(f'{x:.2f}({s})' for x, s in top[:3]) + f' | mean w/o top1/3/5: {tr(1):.4f} {tr(3):.4f} {tr(5):.4f}')
    if code == 'cs' or '--all' in sys.argv:
        for nj, i in ((3, 0), (4, 1)):
            v = [sets['lo'][s][i] for s in seeds]; n = len(v)
            if code == 'cs':
                lo, hi = min(v), sorted(v)[int(.98*n)]
                nb = 10; w = (max(v)-lo)/nb
                h = [0]*nb
                for x in v: h[min(nb-1, int((x-lo)/w))] += 1
                print(f'  hist cs lo >={nj}j (bin width {w:.3f} fb from {lo:.3f}): ' + ' '.join(map(str, h)))
            print(f'  running mean {code} lo >={nj}j: ' + ' '.join(f'{k}:{sum(v[:k])/k:.3f}' for k in (12, 25, 50, 75, 100, 112) if k <= n))
