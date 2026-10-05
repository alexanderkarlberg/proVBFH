#!/usr/bin/env python3
"""CPU needed per bin to reach the errors of a reference, from a pilot.

Usage:
  pilot_estimate.py --excl 'pilot/p1506/excl/job-*' --incl 'pilot/p1506/incl/job-*'
                    --ref reference/p1506/11.top [--weight W1] [--groups p1506]
                    [--out needed.dat]

For each part, the per-bin error of the mean is taken from the scatter of
the pilot runs (std/sqrt(N)), and the part's pilot CPU from the jobs'
time.log (user + system time). With a = CPU_e err_e^2, b = CPU_i err_i^2,
the total error at CPU (C_e, C_i) is sqrt(a/C_e + b/C_i); at the optimal
split the CPU needed for a target error t is (sqrt(a) + sqrt(b))^2 / t^2.
The target t is the reference's quoted error in that bin, or without --ref
a relative error (--relerr) of the pilot's own value; --excl may then be
omitted (LO: inclusive part only). Bins with zero
value or error in the reference, or with no pilot entries, are skipped.
CPU is in CPU-hours of the cluster on which the pilot ran.
"""
import argparse
import glob
import math
import os
import re
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from combine_parts import read_top, average  # noqa: E402


def cpu_hours(dirs):
    tot = 0.0
    for d in dirs:
        with open(os.path.join(d, "time.log")) as f:
            t = f.read()
        u = float(re.search(r"User time \(seconds\): ([\d.]+)", t).group(1))
        s = float(re.search(r"System time \(seconds\): ([\d.]+)", t).group(1))
        tot += u + s
    return tot / 3600.0


def part(pattern, weight, prefix):
    dirs = sorted(d for d in glob.glob(pattern) if os.path.exists(os.path.join(d, "done")))
    files = [os.path.join(d, f"pwg-{prefix}-{weight}.top") for d in dirs]
    files = [f for f in files if os.path.exists(f)]
    return files, cpu_hours([os.path.dirname(f) for f in files])


# histogram groups for the summary: (label, test on (index, name))
GROUPS = {
    "p1506": [
        ("sigma(2 jets)", lambda i, n: n == "sig(all VBF cuts 2 jets)"),
        ("sigma(>=3 jets)", lambda i, n: n == "sig(all VBF cuts 3 jets)"),
        ("sigma(>=4 jets)", lambda i, n: n == "sig(all VBF cuts 4 jets)"),
        ("2-jet distributions", lambda i, n: 4 <= i <= 16),
        ("3-jet distributions", lambda i, n: 17 <= i <= 20),
        ("4-jet distributions", lambda i, n: 21 <= i <= 22),
    ],
    "hxswg136": [
        ("sigma incl cuts (ptj > 20 GeV)", lambda i, n: i == 0),
        ("all bins", lambda i, n: True),
    ],
}


def pct(v, p):
    v = sorted(v)
    return v[min(len(v) - 1, max(0, int(round(p / 100 * (len(v) - 1)))))]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--excl", help="omit for an inclusive-only run (LO)")
    ap.add_argument("--incl", required=True)
    ap.add_argument("--ref", help="target errors: the quoted errors of this file")
    ap.add_argument("--relerr", type=float, default=0.002,
                    help="without --ref: target relative error per bin")
    ap.add_argument("--weight", default="W1")
    ap.add_argument("--groups", default="p1506")
    ap.add_argument("--out")
    a = ap.parse_args()

    fe, ce = part(a.excl, a.weight, "EXCL") if a.excl else ([], 0.0)
    fi, ci = part(a.incl, a.weight, "LO")
    inc, order = average(fi, "scatter")
    ex = average(fe, "scatter")[0] if fe else {}
    if a.ref:
        ref, order = read_top(a.ref)
    else:
        # target: a fixed relative error of the pilot's own central value
        ref = {}
        for h in order:
            ref[h] = []
            for k, r in enumerate(inc[h]):
                v = r[2] + (ex[h][k][2] if h in ex else 0.0)
                ref[h].append([r[0], r[1], v, a.relerr * abs(v)])
    print(f"pilot: {len(fe)} exclusive runs, {ce:.1f} CPU-h; "
          f"{len(fi)} inclusive runs, {ci:.1f} CPU-h")

    need = {}  # name -> list of (xlo, xhi, CPU-h needed, rel err target)
    out = open(a.out, "w") if a.out else None
    for idx, h in enumerate(order):
        if h not in inc:
            continue
        need[h] = []
        for k, r in enumerate(ref[h]):
            if k >= len(inc[h]):
                break
            t = r[3]
            ee = ex[h][k][3] if h in ex and k < len(ex[h]) else 0.0
            ei = inc[h][k][3] if h in inc and k < len(inc[h]) else 0.0
            if not (math.isfinite(r[2]) and math.isfinite(t)) or r[2] == 0 or t <= 0 or (ee == 0 and ei == 0):
                continue
            c = (math.sqrt(ce) * ee + math.sqrt(ci) * ei) ** 2 / t ** 2
            need[h].append((r[0], r[1], c, t / abs(r[2])))
            if out:
                out.write(f"{idx:3d} {r[0]:12.5g} {r[1]:12.5g} {c:12.4g}  {h}\n")
    if out:
        out.close()

    print(f"{'group':34s} {'bins':>5s} {'median':>10s} {'16%':>10s} {'84%':>10s} {'worst':>10s}  (CPU-h)")
    for label, test in GROUPS[a.groups]:
        v = [c for i, h in enumerate(order) if h in need and test(i, h) for (_, _, c, _) in need[h]]
        if not v:
            continue
        print(f"{label:34s} {len(v):5d} {statistics.median(v):10.3g} {pct(v, 16):10.3g} "
              f"{pct(v, 84):10.3g} {max(v):10.3g}")
    if a.groups == "p1506":
        print("per histogram (median, worst):")
        for i, h in enumerate(order):
            if need.get(h) and len(need[h]) > 1:
                v = [c for (_, _, c, _) in need[h]]
                print(f"  {i:3d} {h:34s} {statistics.median(v):10.3g} {max(v):10.3g}")


if __name__ == "__main__":
    main()
