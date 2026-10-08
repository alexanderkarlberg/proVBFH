#!/usr/bin/env python3
"""Per-channel sigma(>= 3 jets), sigma(>= 4 jets): plain means over jobs,
errors from the seed scatter (never inverse-variance weights).

  chsum.py cs <jobs.list>              proVBFH-cs CHAN_MULTI jobs (pwg-EXCL-W1..W9.top)
  chsum.py vbfnlo <dir> [<dir> ...]    VBFNLO, one directory of job-* per channel
                                       (p1506_nlo.top = last iteration only)
Values in fb. With --json FILE the per-job values are written as well.
"""
import glob
import json
import math
import os
import sys

KEYS = ('sig(all VBF cuts 3 jets)', 'sig(all VBF cuts 4 jets)')
CS_NAMES = ['all', 'NC qq', 'NC qg', 'NC gq', 'NC gg', 'CC qq', 'CC qg', 'CC gq', 'CC gg']


def read_sig(path):
    out, name = {}, None
    for line in open(path):
        if line.startswith('#'):
            name = line[1:].split(' index')[0].strip()
        elif name in KEYS and line.strip():
            out[name] = 1000 * float(line.split()[2].replace('D', 'E'))   # pb -> fb
            name = None
    return [out[k] for k in KEYS]


def stats(v):
    n = len(v)
    m = sum(v) / n
    var = sum((x - m) ** 2 for x in v) / (n - 1) if n > 1 else 0.0
    return m, math.sqrt(var / n), math.sqrt(var), n


def main():
    args = sys.argv[1:]
    jout = None
    if '--json' in args:
        i = args.index('--json')
        jout = args[i + 1]
        del args[i:i + 2]
    mode = args[0]
    per = {}
    if mode == 'cs':
        dirs = [l.strip() for l in open(args[1]) if l.strip()]
        done = [d for d in dirs if os.path.exists(os.path.join(d, 'done'))]
        for k, name in enumerate(CS_NAMES, 1):
            per[name] = [read_sig(os.path.join(d, f'pwg-EXCL-W{k}.top')) for d in done]
        per['sum of channels'] = [[sum(per[n][j][i] for n in CS_NAMES[1:]) for i in range(2)]
                                  for j in range(len(done))]
    else:
        for d in args[1:]:
            tops = [os.path.join(j, 'p1506_nlo.top') for j in sorted(glob.glob(os.path.join(d, 'job-*')))
                    if os.path.exists(os.path.join(j, 'done'))]
            per[os.path.basename(d.rstrip('/'))] = [read_sig(t) for t in tops]
    print(f'{"channel":18s} {"jobs":>5s} {"sigma(>=3j) [fb]":>22s} {"sigma(>=4j) [fb]":>22s} {"job scatter 3j":>15s}')
    for name, v in per.items():
        if not v:
            continue
        s3 = stats([x[0] for x in v])
        s4 = stats([x[1] for x in v])
        print(f'{name:18s} {s3[3]:5d} {s3[0]:12.3f} +- {s3[1]:6.3f} {s4[0]:12.3f} +- {s4[1]:6.3f} {s3[2]:15.3f}')
    if jout:
        json.dump(per, open(jout, 'w'))


if __name__ == '__main__':
    main()
