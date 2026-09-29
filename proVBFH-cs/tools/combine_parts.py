#!/usr/bin/env python3
"""Combine proVBFH-cs runs and compare with a reference .top file.

Usage:
  combine_parts.py --part 'incl-s*/pwg*-LO*.top' --part 'excl-s*/pwg*-EXCL*.top'
                   [--ref total_distrib.top] [--out total.top] [--hist NAME ...]

Each --part is a glob of .top files from runs of one part with equal
statistics: they are averaged (mean, error sqrt(sum err^2)/N). The parts
are then added (errors in quadrature). With --ref, every histogram is
compared bin by bin: chi2/n and the largest pull, and for selected
histograms (--hist, default the VBF-cut ones) the values are printed.
The .top format is that of pwhg_bookhist-multi: blocks started by
'# <title> index <n>' followed by 'xlo xhi value error' rows.
"""
import argparse
import glob
import math
import re
import sys


def read_top(path):
    hists, order, name = {}, [], None
    with open(path) as f:
        for line in f:
            m = re.match(r"\s*#\s*(.*?)\s+index\s+\d+\s*$", line)
            if m:
                name = m.group(1).strip()
                hists[name] = []
                order.append(name)
                continue
            w = line.split()
            if name is None or len(w) < 4 or line.lstrip().startswith("#"):
                continue
            try:
                hists[name].append([float(v.replace('D', 'E').replace('d', 'e')) for v in w[:4]])
            except ValueError:
                pass
    return hists, order


def average(files):
    data = [read_top(f)[0] for f in files]
    order = read_top(files[0])[1]
    out = {}
    for h in order:
        rows = [d[h] for d in data if h in d]
        n = len(rows)
        out[h] = [[r[0][0], r[0][1], sum(x[i][2] for x in rows) / n,
                   math.sqrt(sum(x[i][3] ** 2 for x in rows)) / n] for i, r in enumerate(zip(*rows))]
    return out, order


def add(parts):
    base, order = parts[0]
    out = {h: [row[:] for row in base[h]] for h in order}
    for p, _ in parts[1:]:
        for h in order:
            for i, row in enumerate(p[h]):
                out[h][i][2] += row[2]
                out[h][i][3] = math.hypot(out[h][i][3], row[3])
    return out, order


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--part", action="append", required=True)
    ap.add_argument("--ref")
    ap.add_argument("--out")
    ap.add_argument("--hist", action="append")
    a = ap.parse_args()
    parts = []
    for g in a.part:
        files = sorted(glob.glob(g))
        if not files:
            sys.exit(f"no files for {g}")
        print(f"part {g}: {len(files)} files")
        parts.append(average(files))
    tot, order = add(parts)
    if a.out:
        with open(a.out, "w") as f:
            for k, h in enumerate(order):
                f.write(f"# {h} index {k}\n")
                for r in tot[h]:
                    f.write(f"{r[0]:14.6e} {r[1]:14.6e} {r[2]:14.6e} {r[3]:14.6e}\n")
                f.write("\n\n")
    if not a.ref:
        return
    ref, _ = read_top(a.ref)
    sel = a.hist or [h for h in order if "vbf" in h.lower() or "VBF" in h]
    print(f"{'histogram':32s} {'chi2/n':>10s} {'max|pull|':>9s} {'<err new/err ref>':>18s}")
    for h in order:
        if h not in ref:
            continue
        c2, n, mp, er = 0.0, 0, 0.0, []
        for r, s in zip(tot[h], ref[h]):
            e = math.hypot(r[3], s[3])
            if e <= 0 or (r[2] == 0 and s[2] == 0):
                continue
            p = (r[2] - s[2]) / e
            c2 += p * p; n += 1; mp = max(mp, abs(p))
            if s[3] > 0:
                er.append(r[3] / s[3])
        if n:
            print(f"{h:32s} {c2:6.1f}/{n:<3d} {mp:9.2f} {sum(er)/len(er) if er else 0:18.2f}")
    for h in sel:
        if h in ref and h in tot:
            print(f"\n{h}:  bin, new, ref, (new-ref)/err")
            for r, s in zip(tot[h], ref[h]):
                e = math.hypot(r[3], s[3])
                print(f"  [{r[0]:8.3g},{r[1]:8.3g}]  {r[2]:12.5e} +- {r[3]:9.2e}   {s[2]:12.5e} +- {s[3]:9.2e}"
                      f"   {((r[2]-s[2])/e) if e > 0 else 0:6.2f}")


if __name__ == "__main__":
    main()
