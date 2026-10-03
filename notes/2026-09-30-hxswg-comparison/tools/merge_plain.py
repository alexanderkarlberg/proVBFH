#!/usr/bin/env python3
"""Direct (untrimmed) merge of the HXSWG 13.6 TeV NNLO per-seed files.

For each scale choice (HH = mu0/2, 11 = mu0, 22 = 2 mu0; directories with the
files pwg-NNNN-NNLO.top of the study's 13.6TeV_NNLO/{HH,11,22}.tgz):
  - plain mean over all seeds per bin, error = seed scatter / sqrt(N)
    (written as nnlo-{HH,11,22}-plain.top);
  - the study's trimmed combination (combine_runs.f, exact re-implementation
    in rawload.py) as a check against the published files.
The scale band as in the study (to be verified, see --check): nnlo-central =
11, nnlo-min/max = the per-bin minimum/maximum over HH, 11, 22, written as
nnlo-{central,min,max}-plain.top; the errors of min/max are those of the
scale that gives the extremum.
Format as the study's files: xlow xhigh value error fraction (per bin width;
fraction = 1, all seeds kept).

Also a robust third estimate (nnlo-*-symtrim.top): per bin the symmetric
trimmed mean dropping the lowest and highest 0.5% of the seeds (fixed
fraction, both tails equally, no off-by-one), error from 400 bootstrap
resamplings of the seeds (the error of this estimator, unlike the study's
combined VEGAS errors); the band built the same way.

Usage: merge_plain.py RAWDIR_HH RAWDIR_11 RAWDIR_22 OUTDIR [STUDY_RESULTS_DIR]
"""
import glob, os, re, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rawload import combine_runs


def load(d):
    files = sorted(glob.glob(d + '/pwg-*-NNLO.top'))
    names, edges = [], []
    with open(files[0]) as f:
        cur = None
        for line in f:
            m = re.match(r'\s*#\s*(.*?)\s+index\s+\d+', line)
            if m:
                cur = m.group(1); continue
            w = line.split()
            if cur and len(w) >= 4:
                names.append(cur); edges.append((float(w[0]), float(w[1])))
    nb = len(names)
    rows_v, rows_e, used, bad = [], [], [], []
    for fn in files:
        v, e = [], []
        with open(fn) as f:
            for line in f:
                if line.lstrip().startswith('#'):
                    continue
                w = line.split()
                if len(w) >= 4:
                    v.append(float(w[2].replace('D', 'E'))); e.append(float(w[3].replace('D', 'E')))
        if len(v) != nb:
            bad.append(os.path.basename(fn)); continue      # incomplete file
        rows_v.append(v); rows_e.append(e); used.append(fn)
    if bad:
        print('%s: %d incomplete files skipped: %s' % (d, len(bad), ' '.join(bad[:10])))
    return np.array(rows_v), np.array(rows_e), names, np.array(edges), used


def symtrim(val, f=0.005, nboot=400, seed=1):
    """symmetric trimmed mean per bin and its bootstrap error"""
    n = val.shape[0]; k = int(f*n)
    t = np.sort(val, axis=0)[k:n-k].mean(axis=0)
    rng = np.random.default_rng(seed)
    bt = np.empty((nboot, val.shape[1]))
    for b in range(nboot):
        bt[b] = np.sort(val[rng.integers(0, n, n)], axis=0)[k:n-k].mean(axis=0)
    return t, bt.std(axis=0), 1 - 2*k/n


def write(fn, names, edges, val, err, frac, header):
    with open(fn, 'w') as f:
        f.write('# ' + header + '\n# column 5 is the fraction of seeds used\n')
        prev, idx = None, -1
        for n, (lo, hi), v, e, fr in zip(names, edges, val, err, frac):
            if n != prev:
                idx += 1
                f.write('\n\n# %s index %3d\n' % (n, idx))
                prev = n
            f.write('  %16.8E %16.8E %16.8E %16.8E %16.8E\n' % (lo, hi, v, e, fr))


def read_top(fn):
    vals = []
    for line in open(fn):
        if line.lstrip().startswith('#'):
            continue
        w = line.split()
        if len(w) >= 4:
            vals.append((float(w[2].replace('D', 'E')), float(w[3].replace('D', 'E'))))
    return np.array(vals)


def main():
    dirs = dict(zip(['HH', '11', '22'], sys.argv[1:4]))
    out = sys.argv[4]
    study = sys.argv[5] if len(sys.argv) > 5 else None
    os.makedirs(out, exist_ok=True)
    plain, perr, trim, terr, robust, rerr = {}, {}, {}, {}, {}, {}
    names = edges = None
    for sc, d in dirs.items():
        val, err, nm, ed, files = load(d)
        if names is None:
            names, edges = nm, ed
        assert nm == names
        n = len(files)
        plain[sc] = val.mean(axis=0)
        perr[sc] = val.std(axis=0, ddof=1)/np.sqrt(n)
        t = np.array([combine_runs(val[:, k], err[:, k]) for k in range(val.shape[1])])
        trim[sc], terr[sc] = t[:, 0], t[:, 1]
        st, se, fk = symtrim(val)
        robust[sc], rerr[sc] = st, se
        write(os.path.join(out, 'nnlo-%s-symtrim.top' % sc), names, edges, st, se, np.full(len(names), fk),
              'symmetric 0.5%% trimmed mean of %d seeds (%s), bootstrap error; merge_plain.py' % (n, sc))
        print('%s: symmetric 0.5%% trim %.5f +- %.5f' % (sc, st[0], se[0]))
        write(os.path.join(out, 'nnlo-%s-plain.top' % sc), names, edges, plain[sc], perr[sc],
              np.ones(len(names)), 'plain mean of %d seeds (%s), error = seed scatter/sqrt(N); merge_plain.py' % (n, sc))
        print('%s: %d seeds; sig incl cuts (ptj > 20): plain %.5f +- %.5f, trimmed %.5f +- %.5f'
              % (sc, n, plain[sc][0], perr[sc][0], trim[sc][0], terr[sc][0]))
    P = np.array([plain[s] for s in ('HH', '11', '22')])
    PE = np.array([perr[s] for s in ('HH', '11', '22')])
    imin, imax = P.argmin(axis=0), P.argmax(axis=0)
    k = np.arange(P.shape[1])
    one = np.ones(len(names))
    write(os.path.join(out, 'nnlo-central-plain.top'), names, edges, P[1], PE[1], one,
          'plain merge, central scale (11); merge_plain.py')
    write(os.path.join(out, 'nnlo-min-plain.top'), names, edges, P[imin, k], PE[imin, k], one,
          'plain merge, per-bin minimum over HH, 11, 22; merge_plain.py')
    write(os.path.join(out, 'nnlo-max-plain.top'), names, edges, P[imax, k], PE[imax, k], one,
          'plain merge, per-bin maximum over HH, 11, 22; merge_plain.py')
    R = np.array([robust[s] for s in ('HH', '11', '22')])
    RE = np.array([rerr[s] for s in ('HH', '11', '22')])
    jmin, jmax = R.argmin(axis=0), R.argmax(axis=0)
    fr = np.full(len(names), 0.99)
    write(os.path.join(out, 'nnlo-central-symtrim.top'), names, edges, R[1], RE[1], fr,
          'symmetric 0.5% trim, central scale (11); merge_plain.py')
    write(os.path.join(out, 'nnlo-min-symtrim.top'), names, edges, R[jmin, k], RE[jmin, k], fr,
          'symmetric 0.5% trim, per-bin minimum over HH, 11, 22; merge_plain.py')
    write(os.path.join(out, 'nnlo-max-symtrim.top'), names, edges, R[jmax, k], RE[jmax, k], fr,
          'symmetric 0.5% trim, per-bin maximum over HH, 11, 22; merge_plain.py')
    if study:
        # check: the trimmed combinations reproduce the study's files, and
        # its min/max are the envelope of HH, 11, 22
        T = np.array([trim[s] for s in ('HH', '11', '22')])
        for name, ref in (('central', T[1]), ('min', np.nanmin(T, axis=0)), ('max', np.nanmax(T, axis=0))):
            st = read_top(os.path.join(study, 'nnlo-%s.top' % name))[:, 0]
            ok = np.isfinite(ref) & (st != 0)
            print('study nnlo-%s.top vs trimmed %s: largest relative difference %.2e (%d bins)'
                  % (name, 'envelope' if name != 'central' else '11', np.max(np.abs(ref[ok]/st[ok] - 1)), ok.sum()))
        rel = (P[1] - T[1])/np.where(T[1] != 0, T[1], np.nan)
        pull = (P[1] - T[1])/np.where(PE[1] > 0, PE[1], np.nan)
        print('central, plain - trimmed: median relative %+.4f, median in plain errors %+.2f, positive in %.0f%% of %d bins'
              % (np.nanmedian(rel), np.nanmedian(pull), 100*np.mean(rel[np.isfinite(rel)] > 0), np.isfinite(rel).sum()))


if __name__ == '__main__':
    main()
