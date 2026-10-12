#!/usr/bin/env python3
"""V1 comparer: new combined run against the old cs_part 1 + cs_part 2 runs.

  v1_compare.py <new dir> <ref incl dir> <ref excl dir or ->  [--tol 1e-12]

Parts: pwg-LO-W<k>.top and pwg-EXCL-W<k>.top of the new run against the
reference parts, every bin (value and error); same print format, so these
agree to the tolerance (default 1e-12 relative). Totals: pwg-W<k>.top
against the sum of the reference parts (errors in quadrature); the
files carry 8 digits, so these are compared at the print precision (2e-7).
Prints one line, exit code 0 if everything agrees.
"""
import glob, math, os, re, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'tools'))
from combine_parts import read_top


def rel(a, b, scale=None):
    s = scale if scale is not None else max(abs(a), abs(b))
    return 0.0 if s == 0 else abs(a - b) / s


def compare(fa, fb, tol, ncmp):
    ha, oa = read_top(fa)
    hb, ob = read_top(fb)
    if oa != ob:
        return None, 'histogram lists differ'
    worst = 0.0
    for h in oa:
        if len(ha[h]) != len(hb[h]):
            return None, 'bin numbers differ in ' + h
        for r, s in zip(ha[h], hb[h]):
            ncmp[0] += 1
            for i in (2, 3):
                worst = max(worst, rel(r[i], s[i]))
    return worst, ''


def main():
    a = [x for x in sys.argv[1:] if not x.startswith('--')]
    tol = 1e-12
    if '--tol' in sys.argv:
        tol = float(sys.argv[sys.argv.index('--tol') + 1])
        a = [x for x in a if x != str(tol) and x != sys.argv[sys.argv.index('--tol') + 1]]
    new, rin, rex = a[0], a[1], a[2]
    ncmp, worst_part, worst_tot, ok, msg = [0], 0.0, 0.0, True, []
    ws = sorted(glob.glob(os.path.join(new, 'pwg-LO-W*.top')))
    if not ws:
        ws = sorted(glob.glob(os.path.join(new, 'pwg-LO.top')))   # single weight
    for fn in ws:
        w = re.search(r'LO(.*)\.top', os.path.basename(fn)).group(1)
        parts = [(fn, os.path.join(rin, 'pwg-LO%s.top' % w))]
        if rex != '-':
            parts.append((os.path.join(new, 'pwg-EXCL%s.top' % w), os.path.join(rex, 'pwg-EXCL%s.top' % w)))
        for fa, fb in parts:
            x, m = compare(fa, fb, tol, ncmp)
            if x is None:
                ok = False; msg.append(m); continue
            worst_part = max(worst_part, x)
        # totals: new total against the sum of the reference parts
        tn, _ = read_top(os.path.join(new, 'pwg%s.top' % w))
        ti, oi = read_top(parts[0][1])
        te = read_top(parts[1][1])[0] if rex != '-' else None
        for h in oi:
            for k, r in enumerate(tn[h]):
                v = ti[h][k][2] + (te[h][k][2] if te else 0.0)
                e = math.hypot(ti[h][k][3], te[h][k][3] if te else 0.0)
                sc = abs(ti[h][k][2]) + (abs(te[h][k][2]) if te else 0.0)
                ec = ti[h][k][3] + (te[h][k][3] if te else 0.0)
                worst_tot = max(worst_tot, rel(r[2], v, sc), rel(r[3], e, ec))
    ok = ok and worst_part <= tol and worst_tot <= 2e-7
    print('%d bins compared, parts: max rel diff %.1e (tol %.0e), totals vs sum of parts: %.1e (print precision 2e-7) %s'
          % (ncmp[0], worst_part, tol, worst_tot, ' '.join(msg)))
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
