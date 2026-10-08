#!/usr/bin/env python3
"""LO, NLO and NNLO of proVBFH-cs in one plot per histogram.

Usage:
  plot_orders.py --dir combined/p1506 --outdir combined/p1506/plots-orders [--title ...]

Reads <dir>/{lo,nlo,nnlo}-W{1,2,3}.top (those orders that exist; W1 is the
central scale, the band is the envelope of W1-W3). Upper panel: each order
with its scale band and statistical error; lower panel: ratio to the NLO
central value (or the highest order below NNLO that exists).
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from combine_parts import read_top  # noqa: E402
import plot_compare  # noqa: E402
from plot_compare import arrays, safe, step, band  # noqa: E402

COLORS = {"lo": "tab:green", "nlo": "tab:orange", "nnlo": "tab:blue"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--title", default="")
    ap.add_argument("--ratio-range", nargs=2, type=float, metavar=("LO", "HI"))
    ap.add_argument("--xcut", nargs=2, action="append", default=[], metavar=("NAME", "XMAX"))
    ap.add_argument("--ratio-frac", type=float, default=1 / 3.2,
                    help="fraction of the figure height taken by the ratio panel")
    a = ap.parse_args()
    plot_compare.XCUT.update({n: float(x) for n, x in a.xcut})
    orders = [o for o in ("lo", "nlo", "nnlo")
              if all(os.path.exists(os.path.join(a.dir, f"{o}-W{w}.top")) for w in (1, 2, 3))]
    if len(orders) < 2:
        sys.exit(f"need at least two orders in {a.dir}, found {orders}")
    data = {o: [read_top(os.path.join(a.dir, f"{o}-W{w}.top"))[0] for w in (1, 2, 3)] for o in orders}
    order_names = read_top(os.path.join(a.dir, f"{orders[-1]}-W1.top"))[1]
    refo = "nlo" if "nlo" in orders else orders[0]
    os.makedirs(a.outdir, exist_ok=True)
    with PdfPages(os.path.join(a.outdir, "all.pdf")) as allpdf:
        for i, h in enumerate(order_names):
            if any(h not in data[o][0] or len(data[o][0][h]) == 0 for o in orders):
                continue
            lo_, hi_, rc, _ = arrays(data[refo][0], h)
            if any(len(arrays(data[o][0], h)[2]) != len(rc) for o in orders):
                continue
            r = np.where(rc != 0, 1 / np.where(rc != 0, rc, 1), 0)
            fig, (ax, rx) = plt.subplots(2, 1, sharex=True, figsize=(6, 5.5),
                                         gridspec_kw=dict(height_ratios=[1 - a.ratio_frac, a.ratio_frac],
                                                          hspace=0.05))
            xc = 0.5 * (lo_ + hi_)
            vals = []
            for o in orders:
                _, _, c, e = arrays(data[o][0], h)
                v = np.array([arrays(d, h)[2] for d in data[o]])
                mn, mx = v.min(0), v.max(0)
                col = COLORS[o]
                band(ax, lo_, hi_, mn, mx, color=col, alpha=0.25)
                step(ax, lo_, hi_, c, color=col, lw=1.3, label=o.upper())
                ax.errorbar(xc, c, e, fmt="none", ecolor=col, lw=0.9)
                band(rx, lo_, hi_, mn * r, mx * r, color=col, alpha=0.25)
                step(rx, lo_, hi_, c * r, color=col, lw=1.3)
                rx.errorbar(xc, c * r, e * np.abs(r), fmt="none", ecolor=col, lw=0.9)
                vals.append(np.concatenate([mn * r, mx * r]))
            pos = np.concatenate([arrays(data[o][0], h)[2] for o in orders])
            if np.all(pos > 0) and pos.max() / pos.min() > 30:
                ax.set_yscale("log")
            ax.set_ylabel("dσ/dx  [pb / unit]" if len(lo_) > 1 else "σ  [pb]")
            ax.set_title(f"{a.title}: {h}" if a.title else h, fontsize=10)
            ax.legend(fontsize=8, frameon=False)
            rx.axhline(1, color="k", lw=0.8)
            v = np.concatenate(vals)
            v = v[np.isfinite(v) & (v != 0)]
            if a.ratio_range:
                rx.set_ylim(*a.ratio_range)
            elif len(v):
                d = max(0.1, min(0.8, np.percentile(np.abs(v - 1), 95) * 1.2))
                rx.set_ylim(1 - d, 1 + d)
            rx.set_ylabel(f"ratio to {refo.upper()}")
            rx.set_xlabel(h)
            for ext in ("pdf", "png"):
                fig.savefig(os.path.join(a.outdir, f"{i:03d}-{safe(h)}.{ext}"), dpi=110, bbox_inches="tight")
            allpdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)
    print(f"orders {orders} -> {a.outdir}")


if __name__ == "__main__":
    main()
