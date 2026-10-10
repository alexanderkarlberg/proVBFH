#!/usr/bin/env python3
"""x-slice analysis: sigma(>=3j), sigma(>=4j) [fb] per slice and code, plain means, seed-scatter errors.
  analyze.py [--prod]   (production totals included with --prod)"""
import glob, math, os, sys
sys.path.insert(0, '/u/akarlber/work/proVBFH/proVBFH-cs/production/h3j-crosscheck/chsplit')
import chsum
R = '/ptmp/mpp/akarlber/xslice'
def load(code, sl):
    f = 'pwg-EXCL.top' if code == 'cs' else 'p1506_nlo.top'
    if sl == 'prod':
        dirs = glob.glob('/ptmp/mpp/akarlber/cs-production/fixmh/p1506/excl/job-*' if code == 'cs' else '/ptmp/mpp/akarlber/h3j/vbfnlo/prod/nlo-fixmh/job-*')
    else:
        dirs = glob.glob(f'{R}/{code}-{sl}/job-*')
    out = []
    for d in sorted(dirs):
        if os.path.exists(d + '/done') and os.path.exists(f'{d}/{f}'):
            try: r = chsum.read_sig(f'{d}/{f}')
            except Exception: continue
            if len(r) == 2: out.append((d, list(r)))
    return out
def skew(v):
    n = len(v); m = sum(v)/n; s = math.sqrt(sum((x-m)**2 for x in v)/n)
    return sum((x-m)**3 for x in v)/n/s**3 if s > 0 else 0.0
def main():
    sets = {}
    for code in ('cs', 'vb'):
        for sl in ('hi', 'lo1', 'lo2') + (('prod',) if '--prod' in sys.argv else ()):
            L = load(code, sl)
            if L: sets[(code, sl)] = L
    res = {}
    for (code, sl), L in sets.items():
        for i, nj in enumerate((3, 4)):
            v = [x[1][i] for x in L]
            m, e, sd, n = chsum.stats(v)
            j = max(range(n), key=lambda k: v[k])
            res[(code, sl, nj)] = (m, e, n)
            print(f'{code} {sl:5s} >={nj}j  n={n:4d}  mean {m:9.3f} +- {e:6.3f}  job sd {sd:7.2f}  skew {skew(v):5.1f}  max {v[j]:8.2f} ({os.path.basename(L[j][0])}, +{(v[j]-m)/n:.2f} fb on mean)')
    def comb(code, nj, parts):
        t = [res.get((code, p, nj)) for p in parts]
        if any(x is None for x in t): return None
        return sum(x[0] for x in t), math.sqrt(sum(x[1]**2 for x in t))
    print()
    for nj in (3, 4):
        for code in ('cs', 'vb'):
            lo = comb(code, nj, ['lo1', 'lo2']); hi = comb(code, nj, ['hi']); tot = comb(code, nj, ['lo1', 'lo2', 'hi'])
            pr = res.get((code, 'prod', nj))
            print(f'{code} >={nj}j  lo {lo and "%.3f+-%.3f" % lo}  hi {hi and "%.3f+-%.3f" % hi}  lo+hi {tot and "%.3f+-%.3f" % tot}  prod {pr and "%.3f+-%.3f" % pr[:2]}')
        for nm, p in (('lo', ['lo1', 'lo2']), ('hi', ['hi']), ('lo+hi', ['lo1', 'lo2', 'hi'])):
            a, b = comb('cs', nj, p), comb('vb', nj, p)
            if a and b:
                d = a[0]-b[0]; e = math.hypot(a[1], b[1]); print(f'   >={nj}j {nm:5s} cs - vb = {d:+.3f} +- {e:.3f}  pull {d/e:+.1f}')
        for code in ('cs', 'vb'):
            lo = comb(code, nj, ['lo1', 'lo2']); hi = comb(code, nj, ['hi']); pr = res.get((code, 'prod', nj))
            if lo and hi and pr:
                d = pr[0]-hi[0]; e = math.hypot(pr[1], hi[1]); print(f'   {code} >={nj}j: prod - hi = {d:.3f} +- {e:.3f};  lo(dedicated) = {lo[0]:.3f} +- {lo[1]:.3f};  lo - (prod-hi) = {lo[0]-d:+.3f} +- {math.hypot(lo[1], e):.3f}')
main()
