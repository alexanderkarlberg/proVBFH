#!/usr/bin/env python3
"""V3 comparer: every .top file and the xsct file of two run directories,
token by token; numbers must agree to <= tol relative (default 1e-12), all
other text exactly. Every .top and xsct file of the reference run must
exist in the new one (other files, e.g. dumps no longer written, are
ignored).

  v3_compare.py <new dir> <reference dir> [--tol 1e-12]
Prints one line, exit code 0 if everything agrees.
"""
import glob, os, sys


def num(t):
    try:
        return float(t.replace('D', 'E').replace('d', 'e'))
    except ValueError:
        return None


def main():
    a = [x for x in sys.argv[1:] if not x.startswith('--')]
    tol = 1e-12
    if '--tol' in sys.argv:
        tol = float(sys.argv[sys.argv.index('--tol') + 1])
        a = [x for x in a if x != sys.argv[sys.argv.index('--tol') + 1]]
    new, ref = a[0], a[1]
    names = sorted(os.path.basename(f) for f in glob.glob(os.path.join(ref, '*.top')) + glob.glob(os.path.join(ref, 'xsct-*.dat')))
    if not names:
        print('no reference files'); sys.exit(1)
    nfile, ntok, worst, bad = 0, 0, 0.0, []
    for n in names:
        fn = os.path.join(new, n)
        if not os.path.exists(fn):
            bad.append('missing ' + n); continue
        la = open(fn).read().split('\n')
        lb = open(os.path.join(ref, n)).read().split('\n')
        if len(la) != len(lb):
            bad.append('%s: %d lines against %d' % (n, len(la), len(lb))); continue
        nfile += 1
        for x, y in zip(la, lb):
            tx, ty = x.split(), y.split()
            if len(tx) != len(ty):
                bad.append('%s: token count differs in "%s"' % (n, x[:40])); break
            for s, t in zip(tx, ty):
                u, v = num(s), num(t)
                if u is None or v is None:
                    if s != t:
                        bad.append('%s: "%s" != "%s"' % (n, s, t))
                    continue
                ntok += 1
                sc = max(abs(u), abs(v))
                if sc > 0:
                    worst = max(worst, abs(u - v) / sc)
    ok = not bad and worst <= tol
    print('%d files, %d numbers compared, max rel diff %.1e (tol %.0e) %s' % (nfile, ntok, worst, tol, '; '.join(bad[:3])))
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
