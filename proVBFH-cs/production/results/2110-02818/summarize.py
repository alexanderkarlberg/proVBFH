#!/usr/bin/env python3
"""Cross sections (fb) of the >= 2-jet VBF-cut histogram of proVBFH-cs jobs.

  summarize.py <run dir> [<run dir of the other selection> for the per-seed difference]

Plain mean over the finished jobs (marker file "done") of pwg-W<k>.top, error of
the mean from the seed scatter. W1, W2, W3 = (1,1), (1/2,1/2), (2,2) times mu0.
"""
import glob, math, os, re, sys

NAME = "sig(all VBF cuts 2 jets)"

def val(path):
    on = False
    for line in open(path):
        m = re.match(r"\s*#\s*(.*?)\s+index", line)
        if m:
            on = m.group(1).strip() == NAME
            continue
        w = line.split()
        if on and len(w) >= 4:
            return 1000 * float(w[2].replace("D", "E"))
    raise ValueError(path)

def load(run):
    out = {}
    for d in sorted(glob.glob(run + "/job-*")):
        if not os.path.exists(d + "/done"):
            continue
        out[d.rsplit("-", 1)[1]] = [val("%s/pwg-W%d.top" % (d, k)) for k in (1, 2, 3)]
    return out

def stats(x):
    n = len(x); m = sum(x) / n
    s = math.sqrt(sum((v - m) ** 2 for v in x) / (n - 1))
    return m, s / math.sqrt(n)

a = load(sys.argv[1])
print("%s: %d jobs" % (sys.argv[1], len(a)))
for k, lab in enumerate(["mu0", "mu0/2", "2mu0"]):
    m, e = stats([v[k] for v in a.values()])
    print("  %-6s %9.3f +- %6.3f fb  (%.3f%%)" % (lab, m, e, 100 * e / m))
if len(sys.argv) > 2:
    b = load(sys.argv[2])
    common = sorted(set(a) & set(b))
    print("difference %s - %s (%d common seeds)" % (sys.argv[1], sys.argv[2], len(common)))
    for k, lab in enumerate(["mu0", "mu0/2", "2mu0"]):
        m, e = stats([a[s][k] - b[s][k] for s in common])
        print("  %-6s %9.3f +- %6.3f fb" % (lab, m, e))
