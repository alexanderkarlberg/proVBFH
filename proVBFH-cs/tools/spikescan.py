# For every bin of the p1506 NNLO exclusive W1 histograms: the job that moves
# the plain mean most when removed, in units of the seed-scatter error.
import glob, os, re, math, sys, collections
files = sorted(glob.glob('/ptmp/mpp/akarlber/cs-production/prod-nnlo/p1506/excl/job-*/pwg-EXCL-W1.top'))
files = [f for f in files if os.path.exists(os.path.dirname(f) + '/done')]
S = collections.defaultdict(float); S2 = collections.defaultdict(float)
top = collections.defaultdict(list)   # key -> list of (|x|, x, job), keep 3 largest
n = 0
for f in files:
    job = f.split('/')[-2]; name = None
    for l in open(f):
        m = re.match(r'\s*#\s*(.*?)\s+index\s+\d+', l)
        if m: name = m.group(1); ib = 0; continue
        w = l.split()
        if name is None or len(w) < 4: continue
        try: x = float(w[2].replace('D', 'E'))
        except ValueError: continue
        k = (name, ib, w[0], w[1]); ib += 1
        S[k] += x; S2[k] += x * x
        t = top[k]; t.append((abs(x), x, job)); t.sort(reverse=True); del t[3:]
    n += 1
out = []
for k in S:
    m = S[k] / n; var = max(S2[k] / n - m * m, 0); err = math.sqrt(var / (n - 1)) if n > 1 else 0
    if err == 0: continue
    a, x, job = top[k][0]
    shift = (x - m) / (n - 1)          # mean(all) - mean(without job)
    out.append((abs(shift) / err, k, m, err, x, job, abs(x - m) / math.sqrt(var)))
out.sort(reverse=True)
print('jobs', n)
for r in out[:40]:
    s, k, m, err, x, job, z = r
    print(f'{s:6.2f} sigma  {k[0]:28s} bin {k[1]:3d} [{k[2]},{k[3]}]  mean {m:.4e} +- {err:.2e}  job {job} value {x:.3e} ({z:.0f} sd)')
jobs = collections.Counter(r[5] for r in out if r[0] > 0.5)
print('jobs most often responsible (shift > 0.5 sigma):', jobs.most_common(15))
