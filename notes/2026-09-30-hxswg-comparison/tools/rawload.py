"""Load the study's per-seed .top files into numpy arrays (cached in raw11.npz)."""
import glob, os, re, sys
import numpy as np
SP = os.environ.get('RAWDIR', os.path.dirname(os.path.abspath(__file__)))  # directory with raw11/ (the unpacked 11.tgz)
def load(d=SP + '/raw11', cache=SP + '/raw11.npz'):
    if os.path.exists(cache):
        z = np.load(cache, allow_pickle=True)
        return z['val'], z['err'], list(z['names']), z['edges'], list(z['files'])
    files = sorted(glob.glob(d + '/pwg-*-NNLO.top'))
    names, edges, rows = [], [], []
    with open(files[0]) as f:
        cur = None
        for line in f:
            m = re.match(r'\s*#\s*(.*?)\s+index\s+\d+', line)
            if m: cur = m.group(1); continue
            w = line.split()
            if cur and len(w) >= 4:
                names.append(cur); edges.append((float(w[0]), float(w[1])))
    nb = len(names)
    val = np.zeros((len(files), nb)); err = np.zeros((len(files), nb))
    for i, fn in enumerate(files):
        k = 0
        with open(fn) as f:
            for line in f:
                if line.lstrip().startswith('#'): continue
                w = line.split()
                if len(w) >= 4:
                    val[i, k] = float(w[2].replace('D', 'E')); err[i, k] = float(w[3].replace('D', 'E')); k += 1
        assert k == nb, (fn, k, nb)
    np.savez(cache, val=val, err=err, names=np.array(names, dtype=object), edges=np.array(edges), files=np.array(files, dtype=object))
    return val, err, names, np.array(edges), files

def combine_runs(v, e, limit=10.0, offbyone=True):
    """combine_runs.f imethod 0 for one bin: v, e arrays over files."""
    n = len(v)
    idx = np.argsort(v, kind='stable'); s = v[idx]; se = e[idx]
    med = (s[n//2 - 1] + s[n//2]) / 2 if n % 2 == 0 else s[n//2]
    iqr = s[(84*n)//100 - 1] - s[(16*n)//100 - 1]       # Fortran 1-based indices
    lo, hi = med - iqr/2*limit, med + iqr/2*limit
    first = int(np.sum(s < lo)); last = int(np.sum(s < hi))
    nk = last - first
    if offbyone:
        j0 = first            # Fortran loop from first_index (1-based) = 0-based first-1
        if first == 0:        # element 0: out of bounds, reads 0 in practice
            ys = s[:last].sum(); es = (se[:last]**2).sum()
        else:
            ys = s[first-1:last].sum(); es = (se[first-1:last]**2).sum()
    else:
        ys = s[first:last].sum(); es = (se[first:last]**2).sum()
    if nk == 0: return np.nan, np.nan, 0.0
    return ys/nk, np.sqrt(es)/nk, nk/n
