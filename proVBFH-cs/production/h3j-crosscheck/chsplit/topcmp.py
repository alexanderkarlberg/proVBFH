#!/usr/bin/env python3
"""Compare .top files bin by bin (channel-split validation).

  topcmp.py A.top B.top              largest |A - B| and |A - B|/|A| over all bins
  topcmp.py A.top B1.top B2.top ...  A against the sum B1 + B2 + ...

Only the central values are compared (column 3).
"""
import sys


def read_top(path):
    hists, name = {}, None
    for line in open(path):
        if line.startswith('#'):
            name = line[1:].split(' index')[0].strip()
            hists[name] = []
        elif line.strip() and name is not None:
            f = line.split()
            hists[name].append(float(f[2].replace('D', 'E')))
    return hists


def main():
    a = read_top(sys.argv[1])
    bs = [read_top(p) for p in sys.argv[2:]]
    worst_abs, worst_rel, where, n = 0.0, 0.0, '', 0
    for name, va in a.items():
        for i, x in enumerate(va):
            y = sum(b[name][i] for b in bs)
            d = abs(x - y)
            n += 1
            worst_abs = max(worst_abs, d)
            scale = max(abs(x), max(abs(b[name][i]) for b in bs))
            if scale > 0 and d / scale > worst_rel:
                worst_rel, where = d / scale, f'{name} bin {i}: {x:.10e} vs {y:.10e}'
    print(f'{n} bins; max |diff| {worst_abs:.3e}; max |diff|/max|term| {worst_rel:.3e} ({where})')
    for key in ('sig(all VBF cuts 3 jets)', 'sig(all VBF cuts 4 jets)'):
        if key in a:
            print(f'  {key}: {a[key][0]:.12e}  vs  {sum(b[key][0] for b in bs):.12e}')


if __name__ == '__main__':
    main()
