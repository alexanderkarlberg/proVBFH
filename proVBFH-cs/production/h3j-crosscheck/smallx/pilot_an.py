#!/usr/bin/env python3
"""Pilot statistics: per-job sigma(>=3j), sigma(>=4j) [fb] (plain means, seed-scatter errors), skewness, largest jobs, CPU,
error x sqrt(CPU) = sd/sqrt(n) x sqrt(mean CPU * n) = sd x sqrt(mean CPU) (fb sqrt(core-h)).
  pilot_an.py a      low slices lo1, lo2 with the channel (and lo = lo1 + lo2, same seed)
  pilot_an.py b      full setup, off and on"""
import glob, math, os, re, sys
sys.path.insert(0, '/u/akarlber/work/proVBFH-smallx/proVBFH-cs/production/h3j-crosscheck/chsplit')
import chsum
R = '/ptmp/mpp/akarlber/smallx'
def load(name):
    o = {}
    for d in glob.glob(f'{R}/{name}/job-*'):
        if not os.path.exists(d + '/done') or not os.path.exists(d + '/pwg-EXCL.top'): continue
        try: r = chsum.read_sig(d + '/pwg-EXCL.top')
        except Exception: continue
        cpu = 0
        for l in open(d + '/run.log'):
            if 'exclusive part: CPU' in l: cpu = float(l.split()[-2])/3600
        if len(r) == 2: o[int(d.split('-')[-1])] = (r[0], r[1], cpu)
    return o
def mom(v):
    n = len(v); m = sum(v)/n; sd = math.sqrt(sum((x-m)**2 for x in v)/(n-1)); s = math.sqrt(sum((x-m)**2 for x in v)/n)
    return m, sd, (sum((x-m)**3 for x in v)/n/s**3 if s else 0)
def med(v): s = sorted(v); n = len(s); return s[n//2] if n % 2 else (s[n//2-1]+s[n//2])/2
def report(label, D):
    for nj, i in ((3, 0), (4, 1)):
        seeds = sorted(D); v = [D[s][i] for s in seeds]; n = len(v)
        m, sd, sk = mom(v); cpu = sum(D[s][2] for s in seeds)/n
        top = sorted(zip(v, seeds), reverse=True)[:3]
        print(f'{label:7s} >={nj}j n={n:3d} mean {m:.4f}+-{sd/math.sqrt(n):.4f} med {med(v):.4f} sd {sd:.4f} skew {sk:5.2f} | max '
              + ' '.join(f'{x:.2f}({s})' for x, s in top) + f' | CPU/job {cpu:.2f} h  sd*sqrt(CPU) {sd*math.sqrt(cpu):.3f}')
if sys.argv[1] == 'a':
    A, B = load('a-lo1'), load('a-lo2'); seeds = sorted(set(A) & set(B))
    report('lo1', A); report('lo2', B)
    report('lo', {s: tuple(A[s][i] + B[s][i] for i in (0, 1)) + (A[s][2] + B[s][2],) for s in seeds})
else:
    for k in ('off', 'on'):
        D = load('b-' + k)
        if len(D) > 2: report(k, D)
