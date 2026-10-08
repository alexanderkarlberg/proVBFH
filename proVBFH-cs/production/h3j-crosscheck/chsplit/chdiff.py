#!/usr/bin/env python3
"""Per-channel proVBFH-cs - VBFNLO for sigma(all VBF cuts N jets), N=2,3,4; plain means, seed-scatter errors."""
import glob, math, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chsum
chsum.KEYS = tuple(f'sig(all VBF cuts {n} jets)' for n in (2, 3, 4))
R = '/ptmp/mpp/akarlber/chsplit'
st = chsum.stats
def rd(p): return chsum.read_sig(p)
csd = [d for d in sorted(glob.glob(R + '/cs/job-*')) if os.path.exists(d + '/done')]
cs = {n: [rd(f'{d}/pwg-EXCL-W{k}.top') for d in csd] for k, n in enumerate(chsum.CS_NAMES, 1)}
cs['sum'] = [[sum(cs[n][j][i] for n in chsum.CS_NAMES[1:]) for i in range(3)] for j in range(len(csd))]
vb = {}
for c in ['NC-qq','NC-qg','NC-gq','NC-gg','CC-qq','CC-qg','CC-gq','CC-gg']:
    vb[c.replace('-', ' ')] = [rd(j + '/p1506_nlo.top') for j in sorted(glob.glob(f'{R}/vbfnlo/runs/{c}/job-*')) if os.path.exists(j + '/done')]
vb['sum'] = None
names = chsum.CS_NAMES[1:]
nv = min(len(vb[n]) for n in names)
vb['sum'] = [[sum(vb[n][j][i] for n in names) for i in range(3)] for j in range(nv)]
print('cs jobs', len(csd), 'vbfnlo jobs/chan', {n: len(vb[n]) for n in names})
for i, nj in enumerate((2, 3, 4)):
    print(f'\n sigma(all VBF cuts {nj} jets) [fb]')
    for n in chsum.CS_NAMES[1:] + ['sum']:
        a = st([x[i] for x in cs[n]]); b = st([x[i] for x in vb[n]])
        d = a[0] - b[0]; e = math.hypot(a[1], b[1])
        print(f'{n:7s} cs {a[0]:10.3f} +- {a[1]:6.3f}  vb {b[0]:10.3f} +- {b[1]:6.3f}  diff {d:+8.3f} +- {e:6.3f} pull {d/e:+5.1f}')
    if i == 1:
        a = st([x[i] for x in cs['all']]); print(f'cs W1 (all) {a[0]:.3f} +- {a[1]:.3f}')
        # per-job W1 vs sum check
        print('max |W1-sum| per job', max(abs(x[i]-y[i]) for x, y in zip(cs['all'], cs['sum'])))
