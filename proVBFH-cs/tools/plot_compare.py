#!/usr/bin/env python3
"""Comparison plots new (proVBFH-cs) vs old (proVBFH), one per histogram.

Usage:
  plot_compare.py --new W1.top W2.top W3.top --old central.top lo.top hi.top
                  --outdir plots/p1506 [--title '1506.02660 set-up'] [--json out.json]

--new: the combined files per weight, central first (W1 = (1,1)); the scale
band is the envelope of all of them, bin by bin. --old: the reference
central file first, then the files whose envelope (with the central) is the
old band (HH and 22 for 1506.02660; min and max for the HXSWG study).
Upper panel: new central with its band and statistical error, old central
with its band; each --alt set (e.g. the plain and the symmetric-trim merges
of the HXSWG study's seeds) as a dashed or dash-dotted line, band dotted. Lower panel: ratio to the old central. Writes one PDF and
one PNG per histogram, all.pdf with every page, and with --json the numbers
(central, band, errors, chi2 of new vs old central) for a summary page.
"""
import argparse
import json
import math
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from combine_parts import read_top  # noqa: E402


XCUT = {}  # histogram name -> upper x limit shown (bins starting at or above it dropped)


def arrays(hists, name):
    r = np.array(hists[name], dtype=float)
    if name in XCUT:
        r = r[r[:, 0] < XCUT[name]]
    lo, hi, v, e = r[:, 0].copy(), r[:, 1].copy(), r[:, 2].copy(), r[:, 3].copy()
    # overflow/underflow bins with edges like +-1e50 (HXSWG analysis): draw them
    # with the width of the neighbouring bin, values rescaled to that width
    # (so the cross section in the bin is kept)
    if len(lo) > 1:
        for k, nb in ((-1, -2), (0, 1)):
            w = hi[k] - lo[k]
            wn = hi[nb] - lo[nb]
            if wn > 0 and w > 1e6 * wn:
                f = w / wn
                v[k] *= f
                e[k] *= f
                if k == -1:
                    hi[k] = lo[k] + wn
                else:
                    lo[k] = hi[k] - wn
    return lo, hi, v, e


def safe(name):
    return re.sub(r"[^A-Za-z0-9_.+-]+", "_", name).strip("_")


def step(ax, lo, hi, y, **kw):
    x = np.append(lo, hi[-1])
    ax.stairs(y, x, baseline=None, **kw)


def band(ax, lo, hi, ylo, yhi, **kw):
    x = np.append(lo, hi[-1])
    ax.stairs(yhi, x, baseline=ylo, fill=True, **kw)


ALT_STYLE = [("tab:red", "--"), ("tab:green", "-."), ("tab:purple", ":"), ("tab:brown", (0, (5, 1, 1, 1)))]


def old_set(old, h, n):
    """Central, error, band (envelope) of one set of old files; non-finite
    bins (the old combine_runs.f writes +-inf or garbage where every run is
    zero) set to 0 and flagged in fin."""
    _, _, oc, oe = arrays(old[0], h)
    if len(oc) != n:
        return None
    ov = np.array([oc] + [arrays(o, h)[2] for o in old[1:]])
    fin = np.isfinite(oc) & np.isfinite(oe) & np.all(np.isfinite(ov), 0)
    oc, oe = np.where(fin, oc, 0), np.where(fin, oe, 0)
    ov = np.where(fin, ov, 0)
    return oc, oe, ov.min(0), ov.max(0), fin


def chi2_of(nc, ne, oc, oe, fin):
    ok = fin & ((oe > 0) | (ne > 0)) & ((oc != 0) | (nc != 0))
    pulls = np.where(ok, (nc - oc) / np.sqrt(ne ** 2 + oe ** 2 + 1e-300), 0)
    nb = int(ok.sum())
    return float((pulls[ok] ** 2).sum()), nb, (float(np.abs(pulls).max()) if nb else 0.0), ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--new", nargs="+", required=True)
    ap.add_argument("--old", nargs="+", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--title", default="")
    ap.add_argument("--json")
    ap.add_argument("--newlabel", default="proVBFH-cs NNLO")
    ap.add_argument("--oldlabel", default="proVBFH (old) NNLO")
    ap.add_argument("--oldshort", default="study", help="short name of --old in the chi2 text")
    ap.add_argument("--reflabel", help="label of the ratio axis (default: ratio to old)")
    ap.add_argument("--ratio-range", nargs=2, type=float, metavar=("LO", "HI"),
                    help="fixed y range of the ratio panels (curves outside are cut)")
    ap.add_argument("--xcut", nargs=2, action="append", default=[], metavar=("NAME", "XMAX"),
                    help="show histogram NAME only below XMAX (also in the chi2); may be repeated")
    ap.add_argument("--alt", nargs=4, action="append", default=[],
                    metavar=("LABEL", "CENTRAL", "MIN", "MAX"),
                    help="a further old set (e.g. another combination of the same seeds), "
                         "drawn as a line with its band outlined; may be repeated")
    a = ap.parse_args()

    XCUT.update({n: float(x) for n, x in a.xcut})
    new = [read_top(f)[0] for f in a.new]
    old = [read_top(f)[0] for f in a.old]
    alts = [(lab, [read_top(f)[0] for f in files]) for lab, *files in a.alt]
    order = read_top(a.new[0])[1]
    os.makedirs(a.outdir, exist_ok=True)
    summary = []
    with PdfPages(os.path.join(a.outdir, "all.pdf")) as allpdf:
        for i, h in enumerate(order):
            if h not in old[0] or len(new[0][h]) == 0:
                continue
            lo, hi, nc, ne = arrays(new[0], h)
            nv = np.array([arrays(n, h)[2] for n in new])
            nmin, nmax = nv.min(0), nv.max(0)
            o1 = old_set(old, h, len(nc))
            if o1 is None:
                print(f"skip {h}: bin numbers differ")
                continue
            oc, oe, omin, omax, fin = o1
            chi2, nb, maxpull, ok = chi2_of(nc, ne, oc, oe, fin)
            o2s = []
            for (lab, files), (col, ls) in zip(alts, ALT_STYLE):
                o = old_set(files, h, len(nc)) if h in files[0] else None
                if o:
                    c2, n2, m2, _ = chi2_of(nc, ne, o[0], o[1], o[4])
                    o2s.append(dict(label=lab, col=col, ls=ls, o=o, chi2=c2, nbins=n2, maxpull=m2))
            entry = dict(index=i, name=h, lo=lo.tolist(), hi=hi.tolist(),
                         new=nc.tolist(), new_err=ne.tolist(),
                         new_min=nmin.tolist(), new_max=nmax.tolist(),
                         old=oc.tolist(), old_err=oe.tolist(),
                         old_min=omin.tolist(), old_max=omax.tolist(),
                         chi2=chi2, nbins=nb, maxpull=maxpull)
            if o2s:
                entry["alts"] = [dict(label=x["label"], val=x["o"][0].tolist(), err=x["o"][1].tolist(),
                                      min=x["o"][2].tolist(), max=x["o"][3].tolist(),
                                      chi2=x["chi2"], nbins=x["nbins"], maxpull=x["maxpull"])
                                 for x in o2s]
            summary.append(entry)

            fig, (ax, rx) = plt.subplots(2, 1, sharex=True, figsize=(6, 5.5),
                                         gridspec_kw=dict(height_ratios=[2.2, 1], hspace=0.05))
            band(ax, lo, hi, omin, omax, color="tab:gray", alpha=0.35, label=f"{a.oldlabel}, scale band")
            step(ax, lo, hi, oc, color="k", lw=1.2, label=a.oldlabel)
            band(ax, lo, hi, nmin, nmax, color="tab:blue", alpha=0.3,
                 label=f"{a.newlabel}, scale band" if np.any(nmax > nmin) else None)
            step(ax, lo, hi, nc, color="tab:blue", lw=1.4, label=a.newlabel)
            xc = 0.5 * (lo + hi)
            ax.errorbar(xc, nc, ne, fmt="none", ecolor="tab:blue", lw=1)
            ax.errorbar(xc, oc, oe, fmt="none", ecolor="k", lw=1, alpha=0.7)
            for k, x in enumerate(o2s, 1):
                o = x["o"]
                for y in (o[2], o[3]):
                    step(ax, lo, hi, y, color=x["col"], lw=0.6, ls=":")
                step(ax, lo, hi, o[0], color=x["col"], lw=1.1, ls=x["ls"], label=f"{x['label']} (band dotted)")
                ax.errorbar(xc + 0.08 * k * (hi - lo), o[0], o[1], fmt="none", ecolor=x["col"], lw=0.9)
            pos = np.concatenate([nc, oc])
            if np.all(pos > 0) and pos.max() / pos.min() > 30:
                ax.set_yscale("log")
            ax.set_ylabel("dσ/dx  [pb / unit]" if len(lo) > 1 else "σ  [pb]")
            ax.set_title(f"{a.title}: {h}" if a.title else h, fontsize=10)
            ax.legend(fontsize=7, frameon=False)
            txt = f"χ²/n = {chi2:.1f}/{nb}"
            if o2s:
                txt = "χ²/n: " + ", ".join([f"{a.oldshort} {chi2:.1f}/{nb}"] +
                                           [f"{x['label'].split(',')[-1].strip()} {x['chi2']:.1f}/{x['nbins']}" for x in o2s])
            ax.text(0.98, 0.03, txt, transform=ax.transAxes, ha="right", va="bottom", fontsize=8)

            r = np.where(oc != 0, 1 / np.where(oc != 0, oc, 1), 0)
            band(rx, lo, hi, omin * r, omax * r, color="tab:gray", alpha=0.35)
            band(rx, lo, hi, nmin * r, nmax * r, color="tab:blue", alpha=0.3)
            step(rx, lo, hi, nc * r, color="tab:blue", lw=1.4)
            rx.errorbar(xc, nc * r, ne * np.abs(r), fmt="none", ecolor="tab:blue", lw=1)
            rx.errorbar(xc, np.ones_like(xc), oe * np.abs(r), fmt="none", ecolor="k", lw=1, alpha=0.7)
            for k, x in enumerate(o2s, 1):
                o = x["o"]
                for y in (o[2], o[3]):
                    step(rx, lo, hi, y * r, color=x["col"], lw=0.6, ls=":")
                step(rx, lo, hi, o[0] * r, color=x["col"], lw=1.1, ls=x["ls"])
                rx.errorbar(xc + 0.08 * k * (hi - lo), o[0] * r, o[1] * np.abs(r), fmt="none", ecolor=x["col"], lw=0.9)
            rx.axhline(1, color="k", lw=0.8)
            vals = np.concatenate([(nc * r)[ok], (nmin * r)[ok], (nmax * r)[ok]])
            vals = vals[np.isfinite(vals) & (vals != 0)]
            if a.ratio_range:
                rx.set_ylim(*a.ratio_range)
            elif len(vals):
                d = max(0.05, min(0.5, np.percentile(np.abs(vals - 1), 95) * 1.3))
                rx.set_ylim(1 - d, 1 + d)
            rx.set_ylabel(f"ratio to {a.reflabel}" if a.reflabel else
                          (f"ratio to old ({a.oldshort})" if o2s else "ratio to old"))
            rx.set_xlabel(h)
            for ext in ("pdf", "png"):
                fig.savefig(os.path.join(a.outdir, f"{i:03d}-{safe(h)}.{ext}"),
                            dpi=110, bbox_inches="tight")
            allpdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)
    if a.json:
        with open(a.json, "w") as f:
            json.dump(summary, f)
    print(f"{len(summary)} histograms -> {a.outdir}")


if __name__ == "__main__":
    main()
