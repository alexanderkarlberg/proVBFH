#!/usr/bin/env python3
"""Physics checks of proVBFH-inclusive, proVBFH and proVBFH-cs.

Usage (from anywhere; the packages must be built):
  ci/physics-checks.py [options] [check ...]
  ci/physics-checks.py generate [options] [part ...]     (frozen references)

Checks (default: all):
  inclusive   The inclusive cross section (no cuts) of the three
              implementations - provbfh_incl, proVBFH (inclusive_only 1) and
              the inclusive part of proVBFH-cs - at LO, NLO, NNLO (N3LO for
              the first two) and at the scale points (1,1), (1/2,1/2), (2,2),
              each against the frozen reference (ci/reference/inclusive.ref) and
              the implementations against each other. In proVBFH-cs the
              scale points are the on-the-fly ones of one run, in the others
              separate runs.
  projection  proVBFH-cs at NLO and NNLO: the exclusive part (the
              difference between the real radiation and its projection to
              Born) must not change the Higgs-only observables (the cross
              section without cuts, ptH-incl, yH-incl): it integrates to zero
              point by point, so its histograms must vanish to rounding
              precision, at all scale points.
  cs          proVBFH-cs at LO, all histograms at the three scale points
              against the frozen references ci/reference/cs-lo-W<k>.top
              (chi2 over the bins, pull for the single-bin cross sections).
              At NLO and NNLO the exact check replaces it.
  exact       proVBFH-cs at NLO and NNLO: short, fixed runs (fixed seed,
              fixed ncall) whose outputs (histograms, cross-section lines)
              must be reproduced byte for byte against the frozen ones in
              <reference-dir>/exact. The runs are deterministic, so any
              change of the code, the compiler or a library version (hoppet,
              LHAPDF) shifts at least one digit and fails; a dependency or
              compiler update requires the deliberate regeneration
              'ci/physics-checks.py generate exact'. Ignores --seed and
              --stat, so that it stays in sync with the frozen outputs.

All runs use the hoppet benchmark PDF (ci/benchmark-pdf.sh), so no PDF set is
downloaded. The statistics are chosen so that the checks take a few minutes;
--stat F multiplies the number of points of all runs by F (the errors, and so
the tolerance in absolute terms, shrink as 1/sqrt(F)).

Options:
  --runroot DIR   directory for the runs [default: $RUNROOT or a temporary one]
  --stat F        statistics factor [default 1]
  --nsigma N      tolerance of the pulls [default 4]
  --jobs N        runs in parallel [default: number of cores]
  --seed S        seed of the runs [default 1]; the references use other seeds
  --reference-dir DIR   [default: ci/reference]
  --cs-extra K=V  extra key of the proVBFH-cs card (repeatable); to test the
                  sensitivity of the checks, e.g. --cs-extra renscfact=1.2
generate options:
  --seeds N       number of seeds of the exclusive references [default 12]
Parts of generate: inclusive, cs-lo, exact (default: all).
"""
import argparse
import concurrent.futures as cf
import filecmp
import glob
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
PDFNAME = "hoppet_benchmark"
PDFID = 900000
ORDERS = {"lo": 1, "nlo": 2, "nnlo": 3, "n3lo": 4}
MUS = (1.0, 0.5, 2.0)                 # (mu_R, mu_F) = (mu, mu) relative to the central scale
HIGGS_ONLY = ("sig incl cuts", "ptH-incl", "yH-incl")
SIGMA_VBF = "sig(all VBF cuts 2 jets)"

# statistics of the checks (before --stat); the inclusive code takes the
# number of points per iteration, cs maps ncall of the card (n <= 1e4: n,
# n <= 1e7: n/10) for its inclusive part
N_INCL = 200000
N_PROJ = 10000
MAXRELERR = 0.15      # bins with a larger relative error of the run are not part of the chi2

# the exact runs: short, fixed, so that their outputs can be frozen and
# compared byte for byte; they must stay in sync with the frozen outputs
# in <reference-dir>/exact, so they do not follow --seed and --stat
EXACT_SEED = 1
N_EXACT = 2000        # ncall1 = ncall2 of the exact runs, itmx1 = itmx2 = 1


# ----------------------------------------------------------------------------
# .top files
# ----------------------------------------------------------------------------
def read_top(path):
    """{title: [[xlo, xhi, value, error], ...]} in file order."""
    hists, name = {}, None
    with open(path) as f:
        for line in f:
            m = re.match(r"\s*#\s*(.*?)\s+index\s+\d+\s*$", line)
            if m:
                name = m.group(1).strip()
                hists[name] = []
                continue
            w = line.split()
            if name is None or len(w) < 4 or line.lstrip().startswith("#"):
                continue
            hists[name].append([float(v.replace("D", "E").replace("d", "e")) for v in w[:4]])
    return hists


def write_top(path, hists, header=""):
    with open(path, "w") as f:
        for h in header.splitlines():
            f.write("# %s\n" % h)
        for i, (name, rows) in enumerate(hists.items()):
            f.write("\n\n# %s index %3d\n" % (name, i))
            for r in rows:
                f.write(" %.8E %.8E %.8E %.8E\n" % tuple(r))


def combine(files):
    """Mean over runs of equal statistics, bin by bin; the error is the larger
    of the VEGAS error of the mean and the scatter error of the mean."""
    data = [read_top(f) for f in files]
    out = {}
    for name in data[0]:
        rows = []
        for r in zip(*[d[name] for d in data]):
            n = len(r)
            m = sum(x[2] for x in r) / n
            ev = math.sqrt(sum(x[3] ** 2 for x in r)) / n
            es = math.sqrt(sum((x[2] - m) ** 2 for x in r) / (n - 1) / n) if n > 1 else ev
            rows.append([r[0][0], r[0][1], m, max(ev, es)])
        out[name] = rows
    return out


# ----------------------------------------------------------------------------
# running the programs
# ----------------------------------------------------------------------------
class Failure(Exception):
    pass


def copy_inputs(ex, d):
    """Copy only the input files of an example directory (never run
    products: a manual run in the example directory must not break the
    checks, and the runs must not see its leftovers)."""
    for f in os.listdir(ex):
        if f.endswith(".input") or f == "pwgseeds.dat":
            shutil.copy(os.path.join(ex, f), d)


def set_card(src, dst, **kv):
    """Copy a powheg.input, replacing (or appending) the keys in kv."""
    with open(src) as f:
        lines = f.read().split("\n")
    for k, v in kv.items():
        for i, l in enumerate(lines):
            if re.match(r"\s*%s(\s|$)" % re.escape(k), l):
                lines[i] = "%s %s" % (k, v)
                break
        else:
            lines.append("%s %s" % (k, v))
    with open(dst, "w") as f:
        f.write("\n".join(lines) + "\n")


def run(cmd, cwd, stdin=None, label=None):
    """Run a program, log to run.log in cwd; raise Failure with the tail of the log if it fails."""
    os.makedirs(cwd, exist_ok=True)
    log = os.path.join(cwd, "run.log")
    t0 = time.time()
    with open(log, "w") as lf:
        p = subprocess.run(cmd, cwd=cwd, input=stdin, stdout=lf, stderr=subprocess.STDOUT, text=True)
    if p.returncode != 0:
        with open(log) as f:
            tail = "".join(f.readlines()[-30:])
        raise Failure("%s failed (exit %d), last lines of %s:\n%s" % (label or cmd[0], p.returncode, log, tail))
    return time.time() - t0


VAL = r"([-+]?[\d.]+(?:[eEdD][-+]?\d+)?)"


def xsct_values(path):
    """[(value, error), ...] of the cross-section lines of an xsct file."""
    out = []
    with open(path) as f:
        for l in f:
            if l.lstrip().startswith("#") and "+/-" not in l.split("#", 1)[0]:
                continue
            m = re.match(r"\s*%s\s*(?:\+/-)?\s*%s" % (VAL, VAL), l.split("#")[0] if "+/-" in l else l)
            if m and not l.lstrip().startswith("#"):
                out.append((float(m.group(1).replace("D", "E")), float(m.group(2).replace("D", "E"))))
    return out


class Setup:
    def __init__(self, a):
        self.a = a
        self.pdfdir = None

    def prepare(self):
        """benchmark PDF, and the environment for the programs"""
        a = self.a
        self.pdfdir = subprocess.run([os.path.join(HERE, "benchmark-pdf.sh"), os.path.join(a.runroot, "pdf")],
                                     check=True, stdout=subprocess.PIPE, text=True).stdout.strip().splitlines()[-1]
        path = [self.pdfdir]
        if os.environ.get("LHAPDF_DATA_PATH"):
            path.append(os.environ["LHAPDF_DATA_PATH"])
        else:
            path.append(subprocess.run(["lhapdf-config", "--datadir"], check=True, stdout=subprocess.PIPE,
                                       text=True).stdout.strip())
        os.environ["LHAPDF_DATA_PATH"] = ":".join(path)

    # --- provbfh_incl ---------------------------------------------------------
    def incl_exe(self, tag, order, mu, n, seed):
        exe = os.path.join(ROOT, "proVBFH-inclusive", "provbfh_incl")
        d = os.path.join(self.a.runroot, "inclusive", "provbfh_incl", "%s-mu%g-s%d" % (order, mu, seed))
        shutil.rmtree(d, ignore_errors=True)
        cmd = [exe, "-pdf", PDFNAME, "-" + order, "-sqrts", "13000", "-scale-choice", "3",
               "-xmur", str(mu), "-xmuf", str(mu), "-ncall1", str(n), "-ncall2", str(n),
               "-itmx1", "2", "-itmx2", "2", "-iseed", str(seed)]
        run(cmd, d, label="provbfh_incl -" + order)
        v = xsct_values(os.path.join(d, "xsct_%s_seed%04d.dat" % (order, seed)))
        return v[0]

    # --- proVBFH, inclusive_only ----------------------------------------------
    def provbfh_incl_only(self, order, mu, n, seed):
        d = os.path.join(self.a.runroot, "inclusive", "proVBFH", "%s-mu%g-s%d" % (order, mu, seed))
        shutil.rmtree(d, ignore_errors=True)
        os.makedirs(d)
        ex = os.path.join(ROOT, "proVBFH", "example")
        copy_inputs(ex, d)
        set_card(os.path.join(ex, "powheg.input"), os.path.join(d, "powheg.input"),
                 lhans1=PDFID, lhans2=PDFID, qcd_order=ORDERS[order], inclusive_only=1,
                 ncall1=n, ncall2=n, itmx1=2, itmx2=2, iseed=seed, renscfact=mu, facscfact=mu)
        run([os.path.join(ROOT, "proVBFH", "proVBFH")], d, stdin="1\n", label="proVBFH inclusive_only -" + order)
        files = glob.glob(os.path.join(d, "xsct-*.dat"))
        if len(files) != 1:
            raise Failure("proVBFH: expected one xsct file in %s, found %s" % (d, files))
        return xsct_values(files[0])[0]

    # --- proVBFH-cs -----------------------------------------------------------
    def cs_run(self, name, order, seed, ncall, itmx=1, inclusive_only=False, **kv):
        """Run proVBFH-cs in a fresh directory (returned); ncall1 = ncall2 = ncall."""
        d = os.path.join(self.a.runroot, "cs", name)
        shutil.rmtree(d, ignore_errors=True)
        os.makedirs(d)
        ex = os.path.join(ROOT, "proVBFH-cs", "example")
        copy_inputs(ex, d)
        set_card(os.path.join(ex, "powheg.input"), os.path.join(d, "powheg.input"),
                 lhans1=PDFID, lhans2=PDFID, qcd_order=ORDERS[order], inclusive_only=int(inclusive_only),
                 ncall1=ncall, ncall2=ncall, itmx1=itmx, itmx2=itmx, iseed=seed, cs_scales=3,
                 **dict(kv, **self.a.cs_extra))
        t = run([os.path.join(ROOT, "proVBFH-cs", "proVBFH-cs")], d, stdin="1\n", label="proVBFH-cs qcd_order %d" % ORDERS[order])
        return d, t


def reference_inclusive(refdir):
    ref = {}
    with open(os.path.join(refdir, "inclusive.ref")) as f:
        for l in f:
            if l.strip() and not l.startswith("#"):
                k, v, e = l.split()
                ref[k] = (float(v), float(e))
    return ref


# ----------------------------------------------------------------------------
# comparisons
# ----------------------------------------------------------------------------
class Report:
    def __init__(self, nsigma):
        self.nsigma = nsigma
        self.failed = []
        self.uncompared = []
        self.n = 0

    def pull(self, label, a, b, what=("run", "ref")):
        (va, ea), (vb, eb) = a, b
        err = math.hypot(ea, eb)
        p = (va - vb) / err if err > 0 else (0.0 if va == vb else math.inf)
        ok = abs(p) <= self.nsigma
        self.n += 1
        print("  %-4s %-42s %s %.6g +- %.2g, %s %.6g +- %.2g, pull %+.2f" % (
            "ok" if ok else "FAIL", label, what[0], va, ea, what[1], vb, eb, p), flush=True)
        if not ok:
            self.failed.append(label)

    def hist(self, label, run_h, ref_h):
        """all histograms of ref_h against run_h: single-bin ones by pull, the others by chi2 over the bins
        in which the run has a relative error below MAXRELERR. Sparse bins (a handful of points with very
        different weights) have non-Gaussian errors, which the run's own error estimate underestimates,
        so they are not compared; the cut was chosen from runs with several seeds against the references."""
        for name in ref_h:
            if name not in run_h:
                self.n += 1
                self.failed.append("%s %s" % (label, name))
                print("  FAIL %s: histogram '%s' missing in the run" % (label, name))
                continue
            r, g = run_h[name], ref_h[name]
            if len(r) != len(g):
                self.n += 1
                self.failed.append("%s %s" % (label, name))
                print("  FAIL %s %s: %d bins in the run, %d in the reference" % (label, name, len(r), len(g)))
                continue
            if len(r) == 1:
                self.pull("%s %s" % (label, name), (r[0][2], r[0][3]), (g[0][2], g[0][3]))
                continue
            chi2, ndf = 0.0, 0
            for a, b in zip(r, g):
                if a[2] != 0 and 0 < a[3] <= MAXRELERR * abs(a[2]):
                    chi2 += (a[2] - b[2]) ** 2 / (a[3] ** 2 + b[3] ** 2)
                    ndf += 1
            # the integral (sum of the bins) of every histogram: a single number with
            # an error, also where no bin passes the chi2 cut
            # weighted by the bin widths, so that it is the cross section of the histogram
            self.pull("%s %s integral" % (label, name),
                      (sum(a[2] * (a[1] - a[0]) for a in r), math.sqrt(sum((a[3] * (a[1] - a[0])) ** 2 for a in r))),
                      (sum(b[2] * (b[1] - b[0]) for b in g), math.sqrt(sum((b[3] * (b[1] - b[0])) ** 2 for b in g))))
            if ndf == 0:
                self.uncompared.append("%s %s" % (label, name))
                print("  --   %-42s no bin with a relative error below %.0f%%: shape not compared" % (
                    "%s %s" % (label, name), 100 * MAXRELERR), flush=True)
                continue
            limit = ndf + self.nsigma * math.sqrt(2.0 * ndf) + self.nsigma ** 2 / 2   # ~ nsigma of the chi2 distribution
            ok = chi2 <= limit
            self.n += 1
            print("  %-4s %-42s chi2/ndf = %.1f/%d (limit %.1f)" % ("ok" if ok else "FAIL", "%s %s" % (label, name), chi2, ndf, limit), flush=True)
            if not ok:
                self.failed.append("%s %s" % (label, name))


# ----------------------------------------------------------------------------
# the checks
# ----------------------------------------------------------------------------
def check_inclusive(S, rep):
    a = S.a
    ref = reference_inclusive(a.refdir)
    n = int(N_INCL * a.stat)
    ncs = 10 * n if n > 10000 else n          # cs maps n -> n/10 above 1e4
    jobs, res = {}, {}
    with cf.ThreadPoolExecutor(a.jobs) as ex:
        for o in ("lo", "nlo", "nnlo", "n3lo"):
            for mu in MUS:
                jobs[("A", o, mu)] = ex.submit(S.incl_exe, "A", o, mu, n, a.seed)
            jobs[("B", o, 1.0)] = ex.submit(S.provbfh_incl_only, o, 1.0, n, a.seed)
        for o in ("lo", "nlo", "nnlo"):
            jobs[("C", o, 1.0)] = ex.submit(S.cs_run, "inclusive-" + o, o, a.seed, ncs, 2, True)
        for k, j in jobs.items():
            res[k] = j.result()
    names = {"A": "provbfh_incl", "B": "proVBFH", "C": "proVBFH-cs"}
    for o in ("lo", "nlo", "nnlo", "n3lo"):
        print("%s:" % o.upper())
        for mu in MUS:
            r = ref["%s_%g" % (o, mu)]
            rep.pull("%s (mu=%g) provbfh_incl vs ref" % (o, mu), res[("A", o, mu)], r)
        rep.pull("%s proVBFH inclusive_only vs ref" % o, res[("B", o, 1.0)], ref["%s_1" % o])
        rep.pull("%s provbfh_incl vs proVBFH" % o, res[("A", o, 1.0)], res[("B", o, 1.0)], ("A", "B"))
        if o != "n3lo":
            d, _ = res[("C", o, 1.0)]
            vals = xsct_values(os.path.join(d, "xsct-%s.dat" % o))
            for mu, v in zip(MUS, vals):      # W1 (1,1), W2 (1/2,1/2), W3 (2,2)
                rep.pull("%s (mu=%g) proVBFH-cs vs ref" % (o, mu), v, ref["%s_%g" % (o, mu)])
            rep.pull("%s provbfh_incl vs proVBFH-cs" % o, res[("A", o, 1.0)], vals[0], ("A", "C"))


def check_projection(S, rep):
    a = S.a
    n = int(N_PROJ * a.stat)
    jobs = {}
    with cf.ThreadPoolExecutor(a.jobs) as ex:
        for o in ("nlo", "nnlo"):
            jobs[o] = ex.submit(S.cs_run, "projection-" + o, o, a.seed, n)
        res = {o: j.result() for o, j in jobs.items()}
    for o in ("nlo", "nnlo"):
        d, t = res[o]
        print("%s (%.0f s):" % (o.upper(), t))
        for k in (1, 2, 3):
            ex_h = read_top(os.path.join(d, "pwg-EXCL-W%d.top" % k))
            in_h = read_top(os.path.join(d, "pwg-LO-W%d.top" % k))
            for name in HIGGS_ONLY:
                scale = max(abs(r[2]) for r in in_h[name])
                worst = max(abs(r[2]) for r in ex_h[name])
                ok = scale > 0 and worst <= 1e-9 * scale
                rep.n += 1
                print("  %-4s W%d %-14s exclusive part max |bin| = %.2e, inclusive part max |bin| = %.2e" % (
                    "ok" if ok else "FAIL", k, name, worst, scale), flush=True)
                if not ok:
                    rep.failed.append("projection %s W%d %s" % (o, k, name))
            # the exclusive part must not be trivially zero for jet observables
            jet = sum(abs(r[2]) for r in ex_h["ptj1-vbf"])
            rep.n += 1
            ok = jet > 0 and math.isfinite(jet)
            print("  %-4s W%d exclusive part of ptj1-vbf is non-zero and finite (sum |bin| = %.2e)" % ("ok" if ok else "FAIL", k, jet))
            if not ok:
                rep.failed.append("projection %s W%d jets" % (o, k))


def check_cs(S, rep):
    a = S.a
    jobs = {}
    with cf.ThreadPoolExecutor(a.jobs) as ex:
        n = int(N_INCL * a.stat)
        jobs["lo"] = ex.submit(S.cs_run, "ref-lo", "lo", a.seed, 10 * n if n > 10000 else n, 2)
        res = {o: j.result() for o, j in jobs.items()}
    for o in ("lo",):
        d, t = res[o]
        print("%s (%.0f s):" % (o.upper(), t))
        for k in (1, 2, 3):
            rep.hist("%s W%d" % (o, k), read_top(os.path.join(d, "pwg-W%d.top" % k)),
                     read_top(os.path.join(a.refdir, "cs-%s-W%d.top" % (o, k))))


def exact_files(d):
    """The outputs of a cs run that the exact check freezes and compares:
    the histograms and the cross-section lines (not the grids)."""
    files = []
    for pattern in ("xsct-*.dat", "pwg-*.top"):
        files += sorted(glob.glob(os.path.join(d, pattern)))
    return files


def first_diff(path_a, path_b):
    """Description of the first differing line of two text files."""
    with open(path_a) as fa, open(path_b) as fb:
        for i, (la, lb) in enumerate(zip(fa, fb)):
            if la != lb:
                return "line %d: '%s' vs '%s'" % (i + 1, la.strip()[:40], lb.strip()[:40])
    na = sum(1 for _ in fa)
    nb = sum(1 for _ in fb)
    return "%s has %d lines, %s has %d" % (os.path.basename(path_a), na, os.path.basename(path_b), nb)


def check_exact(S, rep):
    """proVBFH-cs at NLO and NNLO: short, fixed runs whose outputs must be
    reproduced byte for byte against the frozen ones in
    <reference-dir>/exact (see the 'exact' check of the module docstring)."""
    a = S.a
    for o in ("nlo", "nnlo"):
        d, t = S.cs_run("exact-" + o, o, EXACT_SEED, N_EXACT, 1)
        print("%s (%.0f s):" % (o.upper(), t), flush=True)
        produced = exact_files(d)
        frozen = sorted(glob.glob(os.path.join(a.refdir, "exact", o, "*")))
        names_p = {os.path.basename(f) for f in produced}
        names_f = {os.path.basename(f) for f in frozen}
        for name in sorted(names_p ^ names_f):
            rep.n += 1
            rep.failed.append("exact %s %s" % (o, name))
            print("  FAIL %s: %s (frozen outputs and the run disagree; regenerate with"
                  " 'ci/physics-checks.py generate exact')" % (name, "not produced by the run" if name in names_f else "not frozen"), flush=True)
        for f in produced:
            name = os.path.basename(f)
            ref = os.path.join(a.refdir, "exact", o, name)
            rep.n += 1
            if not os.path.exists(ref):
                rep.failed.append("exact %s %s" % (o, name))
                print("  FAIL %-22s no frozen output %s" % (name, ref), flush=True)
            elif filecmp.cmp(f, ref, shallow=False):
                print("  ok   %-22s byte-identical" % name, flush=True)
            else:
                rep.failed.append("exact %s %s" % (o, name))
                print("  FAIL %-22s differs from the frozen output: %s" % (name, first_diff(f, ref)), flush=True)


# ----------------------------------------------------------------------------
# frozen references
# ----------------------------------------------------------------------------
def generate(S, parts):
    a = S.a
    os.makedirs(a.refdir, exist_ok=True)
    info = subprocess.run(["git", "-C", ROOT, "rev-parse", "--short", "HEAD"], stdout=subprocess.PIPE, text=True).stdout.strip()
    dirty = subprocess.run(["git", "-C", ROOT, "status", "--porcelain", "--untracked-files=no"], stdout=subprocess.PIPE, text=True).stdout.strip()
    header = ("frozen reference of ci/physics-checks.py, hoppet benchmark PDF (ci/benchmark-pdf.sh),\n"
              "13 TeV, scale of 1506.02660 (runningscales 1), generated by 'ci/physics-checks.py generate'\n"
              "with the code of commit %s%s on %s" % (info, " (+ uncommitted changes)" if dirty else "", time.strftime("%Y-%m-%d")))
    seed0 = 1000                                   # the checks use seed 1
    if "inclusive" in parts:
        print("inclusive reference (provbfh_incl, 4 x 2e6 points per iteration, 8 seeds)", flush=True)
        n = int(2000000 * a.stat)
        out = {}
        with cf.ThreadPoolExecutor(a.jobs) as ex:
            fut = {}
            for o in ("lo", "nlo", "nnlo", "n3lo"):
                for mu in MUS:
                    for s in range(8):
                        fut[(o, mu, s)] = ex.submit(S.incl_exe, "A", o, mu, n, seed0 + s)
            for k, f in fut.items():
                out.setdefault(k[:2], []).append(f.result())
        with open(os.path.join(a.refdir, "inclusive.ref"), "w") as f:
            for h in header.splitlines():
                f.write("# %s\n" % h)
            f.write("# key (order_mu, mu_R = mu_F relative to the central scale)  cross section [pb]  error\n"
                    "# provbfh_incl, mean of 8 seeds, the error is the larger of the MC error and the scatter\n")
            for (o, mu), v in out.items():
                m = sum(x[0] for x in v) / len(v)
                ev = math.sqrt(sum(x[1] ** 2 for x in v)) / len(v)
                es = math.sqrt(sum((x[0] - m) ** 2 for x in v) / (len(v) - 1) / len(v))
                f.write("%s_%g %.8e %.3e\n" % (o, mu, m, max(ev, es)))
    for o, n, itmx in (("lo", 10000000, 2),):
        if "cs-" + o not in parts:
            continue
        ns = 8
        n = int(n * a.stat)
        print("cs-%s reference: %d seeds, ncall1 = ncall2 = %d" % (o, ns, n), flush=True)
        with cf.ThreadPoolExecutor(a.jobs) as ex:
            fut = [ex.submit(S.cs_run, "gen-%s-s%d" % (o, seed0 + s), o, seed0 + s, n, itmx) for s in range(ns)]
            dirs = [f.result()[0] for f in fut]
        for k in (1, 2, 3):
            h = combine([os.path.join(d, "pwg-W%d.top" % k) for d in dirs])
            write_top(os.path.join(a.refdir, "cs-%s-W%d.top" % (o, k)), h,
                      header + "\n%s: mean of %d seeds, ncall1 = ncall2 = %d" % (o, ns, n))
    if "exact" in parts:
        # the frozen outputs must stay byte-identical, so their metadata
        # (code, dependencies, settings) goes into a README next to them
        versions = {}
        for t in ("gfortran", "hoppet-config", "lhapdf-config"):
            versions[t] = subprocess.run([t, "--version"], stdout=subprocess.PIPE,
                                         stderr=subprocess.DEVNULL, text=True).stdout.strip().splitlines()[0]
        for o in ("nlo", "nnlo"):
            d, t = S.cs_run("gen-exact-" + o, o, EXACT_SEED, N_EXACT, 1)
            out = os.path.join(a.refdir, "exact", o)
            os.makedirs(out, exist_ok=True)
            for f in exact_files(d):
                shutil.copy(f, os.path.join(out, os.path.basename(f)))
            print("exact %s reference written (%.0f s)" % (o, t), flush=True)
        with open(os.path.join(a.refdir, "exact", "README.md"), "w") as f:
            f.write("# Exact-check references\n\n"
                    "Frozen outputs of short, fixed proVBFH-cs runs, compared byte for byte\n"
                    "by `ci/physics-checks.py exact` (see its docstring). Regenerate with\n"
                    "`ci/physics-checks.py generate exact` after a deliberate change of the\n"
                    "code or of the dependencies or compiler.\n\n"
                    "- Generated: %s\n- Code: commit %s%s\n- gfortran: %s\n- hoppet: %s\n- LHAPDF: %s\n"
                    "- PDF: the hoppet benchmark (ci/benchmark-pdf.sh)\n"
                    "- Runs: qcd_order 2 and 3, iseed %d, ncall1 = ncall2 = %d, itmx1 = itmx2 = 1,\n"
                    "  13 TeV, scale of 1506.02660 (cs_scales 3), default card otherwise\n"
                    "- Compared: xsct-*.dat and pwg-*.top of each run (not the grids)\n"
                    % (time.strftime("%Y-%m-%d"), info, " (+ uncommitted changes)" if dirty else "",
                       versions["gfortran"], versions["hoppet-config"], versions["lhapdf-config"],
                       EXACT_SEED, N_EXACT))
    print("references written to", a.refdir)


# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(usage=__doc__.split("Options:")[0], add_help=True)
    ap.add_argument("checks", nargs="*")
    ap.add_argument("--runroot", default=None)
    ap.add_argument("--stat", type=float, default=1.0)
    ap.add_argument("--nsigma", type=float, default=4.0)
    ap.add_argument("--jobs", type=int, default=os.cpu_count() or 2)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--seeds", type=int, default=12)
    ap.add_argument("--cs-extra", dest="cs_extra", action="append", default=[])
    ap.add_argument("--reference-dir", dest="refdir", default=os.path.join(HERE, "reference"))
    a = ap.parse_args()
    a.cs_extra = dict(kv.split("=", 1) for kv in a.cs_extra)
    if a.runroot is None:
        a.runroot = os.path.join(os.environ.get("RUNROOT") or tempfile.mkdtemp(), "physics-checks")
    os.makedirs(a.runroot, exist_ok=True)
    a.runroot = os.path.abspath(a.runroot)
    S = Setup(a)
    S.prepare()
    if a.checks and a.checks[0] == "generate":
        parts = a.checks[1:] or ["inclusive", "cs-lo", "exact"]
        generate(S, parts)
        return 0
    checks = a.checks or ["inclusive", "projection", "cs", "exact"]
    funcs = {"inclusive": check_inclusive, "projection": check_projection, "cs": check_cs,
             "exact": check_exact}
    rep = Report(a.nsigma)
    t0 = time.time()
    for c in checks:
        if c not in funcs:
            sys.exit("unknown check %s" % c)
        print("== %s" % c, flush=True)
        try:
            funcs[c](S, rep)
        except Failure as e:
            print("ERROR: %s" % e, file=sys.stderr)
            rep.failed.append(c)
    print("%d comparisons in %.0f s, %d failed" % (rep.n, time.time() - t0, len(rep.failed)))
    if rep.uncompared:
        print("%d histograms with no bin for the shape comparison (their integrals were compared): %s" % (
            len(rep.uncompared), ", ".join(rep.uncompared)))
    for f in rep.failed:
        print("ERROR: physics check failed: %s" % f, file=sys.stderr)
    return 1 if rep.failed else 0


if __name__ == "__main__":
    sys.exit(main())
