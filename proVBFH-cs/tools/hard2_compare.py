#!/usr/bin/env python3
"""Compare sets of proVBFH-cs exclusive jobs (e.g. cs_hardfrac2 pilots):
plain means of the sig(...) totals, per-job scatter, the largest single-job
share, the number of bins dominated by one job, and CPU per job.
  hard2_compare.py LABEL=DIR [LABEL=DIR ...]   (DIR holds job-*/pwg-EXCL.top)"""
import sys, glob, os, re, math, collections
def read(f):
    out = {}; name = None
    for l in open(f):
        m = re.match(r'\s*#\s*(.*?)\s+index\s+\d+', l)
        if m: name = m.group(1); ib = 0; continue
        w = l.split()
        if name is None or len(w) < 4: continue
        try: out[(name, ib)] = float(w[2].replace('D', 'E'))
        except ValueError: continue
        ib += 1
    return out
for arg in sys.argv[1:]:
    lab, d = arg.split('=', 1)
    jobs = [j for j in sorted(glob.glob(d + '/job-*')) if os.path.exists(j + '/done')]
    data = [read(j + '/pwg-EXCL.top') for j in jobs]
    n = len(data)
    cpu = []
    for j in jobs:
        for l in open(j + '/time.log'):
            if 'User time' in l: cpu.append(float(l.split()[-1]) / 3600)
    keys = data[0].keys()
    dom = 0; nb = 0
    for k in keys:
        x = [dd.get(k, 0) for dd in data]; m = sum(x) / n
        v = sum((t - m) ** 2 for t in x) / (n - 1)
        if v <= 0: continue
        nb += 1; e = math.sqrt(v / n)
        if max(abs(t - m) for t in x) / (n - 1) / e > 0.5: dom += 1
    print(f'{lab}: {n} jobs, CPU {sum(cpu)/len(cpu):.2f} h/job; bins dominated by one job (shift > 0.5 sigma): {dom} of {nb}')
    for name in ['sig(all VBF cuts 2 jets)', 'sig(all VBF cuts 3 jets)', 'sig(all VBF cuts 4 jets)']:
        x = [dd[(name, 0)] * 1e3 for dd in data]; m = sum(x) / n
        sd = math.sqrt(sum((t - m) ** 2 for t in x) / (n - 1))
        mx = max(x, key=lambda t: abs(t - m))
        print(f'   {name:26s} {m:9.3f} +- {sd/math.sqrt(n):.3f} fb   per-job sd {sd:8.3f}  largest job {mx:9.2f} ({(mx-m)/sd:+.0f} sd)'
              f'   error x sqrt(CPU) {sd/math.sqrt(n)*math.sqrt(n*sum(cpu)/len(cpu)):.3f}')
