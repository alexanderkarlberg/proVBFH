#!/usr/bin/env python3
"""Plots in the style of Fig. 2 of 1506.02660 from the .top files of this directory.

  plot_hists.py new       -> plots-new/      all histograms of the new proVBFH-cs output
                                              (LO, NLO, NNLO with scale bands, ratio to NLO)
  plot_hists.py overlap   -> plots-overlap/  histograms present in both the new output
                                              and old/ (old proVBFH NNLO), same layout plus
                                              the old NNLO; ratio to the NEW NLO central

New files: {lo,nlo,nnlo}-{central,min,max}.top (4 columns: xlo xhi val err).
Old files: old/11.top (central), old/max.top, old/min.top (envelope of 11, HH, 22).
Columns 5 of the old 11/HH/22 files (fraction of untrimmed seeds) is a diagnostic
and is ignored.

Needs numpy and matplotlib.
"""
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from matplotlib.ticker import MaxNLocator  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))

# new histogram name -> (old name or None, x label, y label, log y)
OBS = {
    "sig incl cuts": (None, "", r"$\sigma$ (inclusive cuts) [pb]", False),
    "sig(all VBF cuts 2 jets)": ("sig(all VBF cuts 2 jets)", "", r"$\sigma$ (VBF cuts) [pb]", False),
    "ptH-incl": (None, r"$p_{t,H}$ [GeV]", r"d$\sigma$/d$p_{t,H}$ [pb/GeV]", True),
    "yH-incl": (None, r"$y_H$", r"d$\sigma$/d$y_H$ [pb]", False),
    "ptj1-vbf": ("ptj1", r"$p_{t,j_1}$ [GeV]", r"d$\sigma$/d$p_{t,j_1}$ [pb/GeV]", True),
    "ptj2-vbf": ("ptj2", r"$p_{t,j_2}$ [GeV]", r"d$\sigma$/d$p_{t,j_2}$ [pb/GeV]", True),
    "ptH-vbf": ("ptH", r"$p_{t,H}$ [GeV]", r"d$\sigma$/d$p_{t,H}$ [pb/GeV]", True),
    "yj1-vbf": ("yj1", r"$y_{j_1}$", r"d$\sigma$/d$y_{j_1}$ [pb]", False),
    "yj2-vbf": ("yj2", r"$y_{j_2}$", r"d$\sigma$/d$y_{j_2}$ [pb]", False),
    "yH-vbf": ("yH", r"$y_H$", r"d$\sigma$/d$y_H$ [pb]", False),
    "M_jj-vbf": ("M_jj", r"$m_{jj}$ [GeV]", r"d$\sigma$/d$m_{jj}$ [pb/GeV]", True),
    "rap(j1,j2)-vbf": ("rap(j1,j2)", r"$\Delta y_{j_1,j_2}$", r"d$\sigma$/d$\Delta y_{j_1,j_2}$ [pb]", False),
    "phi(j1,j2)-vbf": ("phi(j1,j2)", r"$\Delta\phi_{j_1,j_2}$", r"d$\sigma$/d$\Delta\phi_{j_1,j_2}$ [pb]", False),
    "R(j1,j2)-vbf": ("R(j1,j2)", r"$R_{j_1,j_2}$", r"d$\sigma$/d$R_{j_1,j_2}$ [pb]", False),
}

# colours / styles as in Fig. 2 of the paper
LO_C, NLO_C, NNLO_C, OLD_C = "#f5b82e", "#00a8e8", "#e8410f", "#2a7a2a"


def read_top(path):
    """{name: ndarray(nrows, ncols)} in file order."""
    out, name, rows = {}, None, []
    with open(path) as f:
        for line in f:
            s = line.strip()
            m = re.match(r"#\s*(.*?)\s+index\s+\d+", s)
            if m:
                if name is not None:
                    out[name] = np.array(rows)
                name, rows = m.group(1), []
            elif s and not s.startswith("#") and name is not None:
                rows.append([float(x.replace("D", "E")) for x in s.split()])
    if name is not None:
        out[name] = np.array(rows)
    # old files have +-Infinity in empty bins (fifth column, untrimmed fraction, is 0 there)
    for k, a in out.items():
        a[~np.isfinite(a)] = 0.0
    return out


def load_new():
    d = {}
    for o in ("lo", "nlo", "nnlo"):
        for k in ("central", "min", "max"):
            d[(o, k)] = read_top(os.path.join(HERE, f"{o}-{k}.top"))
    return d


def load_old():
    return {k: read_top(os.path.join(HERE, "old", f"{f}.top"))
            for k, f in (("central", "11"), ("min", "min"), ("max", "max"))}


def band_of(h_c, h_min, h_max):
    """(lo, hi) envelope of the three scale choices, so the band always brackets the line."""
    # old min/max files contain +-Infinity in empty bins (mergedata artefact): use central there
    v = np.vstack([h_c[:, 2], h_min[:, 2], h_max[:, 2]])
    v = np.where(np.isfinite(v), v, h_c[:, 2])
    return v.min(0), v.max(0)


def stepify(edges, v):
    return np.repeat(edges, 2)[1:-1], np.repeat(v, 2)


def draw_band(ax, edges, lo, hi, color, hatch=None, alpha=0.45):
    x, l = stepify(edges, lo)
    _, h = stepify(edges, hi)
    if hatch:
        ax.fill_between(x, l, h, facecolor="none", edgecolor=color, hatch=hatch, lw=0)
    else:
        ax.fill_between(x, l, h, color=color, alpha=alpha, lw=0)


def draw_line(ax, edges, v, color, lw=1.1, ls="-"):
    x, y = stepify(edges, v)
    ax.plot(x, y, color=color, lw=lw, ls=ls)


def safe_ratio(a, b):
    with np.errstate(divide="ignore", invalid="ignore"):
        r = a / b
    return np.where(b > 0, r, np.nan)


def make_series(name, new, old=None):
    """Return list of dicts (label, central, lo, hi, err, color, hatch) in drawing order."""
    oldname = OBS[name][0]
    ser = []
    for o, lab, col, hat in (("lo", "LO", LO_C, None), ("nlo", "NLO", NLO_C, None),
                             ("nnlo", "NNLO", NNLO_C, "////")):
        c, mn, mx = (new[(o, k)][name] for k in ("central", "min", "max"))
        lo, hi = band_of(c, mn, mx)
        ser.append(dict(label=lab, h=c, lo=lo, hi=hi, color=col, hatch=hat))
    if old is not None:
        c, mn, mx = (old[k][oldname] for k in ("central", "min", "max"))
        edges_new = new[("nnlo", "central")][name][:, :2]
        if c.shape[0] != edges_new.shape[0] or not np.allclose(c[:, :2], edges_new):
            raise RuntimeError(f"bins of {name} differ between new and old")
        lo, hi = band_of(c, mn, mx)
        ser.append(dict(label="NNLO (old)", h=c, lo=lo, hi=hi, color=OLD_C, hatch="\\\\\\\\"))
    return ser


def title_text(name, overlap):
    vbf = ("vbf" in name) or ("VBF" in name)
    return ("VBF CUTS" if vbf else "INCLUSIVE CUTS") + "\nLHC 13 TeV"


def plot_hist(name, new, old=None):
    _, xlab, ylab, logy = OBS[name]
    ser = make_series(name, new, old)
    h0 = ser[1]["h"]
    edges = np.append(h0[:, 0], h0[-1, 1])
    nlo_c = new[("nlo", "central")][name][:, 2]

    fig, (ax, rx) = plt.subplots(2, 1, figsize=(5.2, 6.4), sharex=True,
                                 gridspec_kw=dict(height_ratios=[1.6, 1], hspace=0.0))
    for s in ser:
        draw_band(ax, edges, s["lo"], s["hi"], s["color"], s["hatch"])
        draw_line(ax, edges, s["h"][:, 2], s["color"])
        rr = lambda v: safe_ratio(v, nlo_c)
        draw_band(rx, edges, rr(s["lo"]), rr(s["hi"]), s["color"], s["hatch"])
        draw_line(rx, edges, rr(s["h"][:, 2]), s["color"])
    # statistical errors on the NNLO(s)
    ctr = 0.5 * (edges[:-1] + edges[1:])
    for s in ser[2:]:
        rx.errorbar(ctr, safe_ratio(s["h"][:, 2], nlo_c), yerr=safe_ratio(s["h"][:, 3], nlo_c),
                    fmt="none", ecolor=s["color"], elinewidth=0.8, capsize=0)

    ax.set_ylabel(ylab)
    if logy:
        ax.set_yscale("log")
        pos = np.concatenate([s["h"][:, 2][s["h"][:, 2] > 0] for s in ser])
        ax.set_ylim(pos.min() * 0.5, pos.max() * 5)
    else:
        top = max(np.nanmax(s["hi"]) for s in ser)
        bot = min(0, min(np.nanmin(s["lo"]) for s in ser))
        ax.set_ylim(bot, top * 1.45)
    nz = np.where(np.any([s["h"][:, 2] > 0 for s in ser], axis=0))[0]
    ax.set_xlim(edges[nz[0]], edges[nz[-1] + 1])   # populated bins only
    if not logy:
        ax.yaxis.set_major_locator(MaxNLocator(prune="lower"))
    handles = [Patch(color=s["color"], alpha=0.5) if s["hatch"] is None
               else Patch(facecolor="none", edgecolor=s["color"], hatch=s["hatch"]) for s in ser]
    ax.legend(handles, [s["label"] for s in ser], loc="upper right", frameon=False, fontsize=9)
    ax.text(0.97, 0.62 if not logy else 0.55, title_text(name, old is not None),
            transform=ax.transAxes, ha="right", va="top", fontsize=8)
    ax.text(0.03, 0.04, "NNPDF30_nnlo_as_118\n" r"$\mu_0/2<\mu_R=\mu_F<2\mu_0$ (dynamic scale)",
            transform=ax.transAxes, fontsize=7, va="bottom")

    # ratio panel range: from the NLO and NNLO(s) bands; LO may be clipped
    sel = [s for s in ser if s["label"] != "LO"]
    r = np.concatenate([safe_ratio(s[k], nlo_c) for s in sel for k in ("lo", "hi")])
    r = r[np.isfinite(r)]
    lo_r, hi_r = (np.percentile(r, 2), np.percentile(r, 98)) if r.size else (0.9, 1.1)
    pad = 0.05 * max(hi_r - lo_r, 0.1)
    rx.set_ylim(max(0.5, min(lo_r - pad, 0.9)), min(1.5, max(hi_r + pad, 1.1)))
    rx.yaxis.set_major_locator(MaxNLocator(prune="upper", nbins=6))
    rx.axhline(1, color="gray", lw=0.4)
    rx.set_ylabel("ratio to NLO" + (" (new)" if old is not None else ""))
    rx.set_xlabel(xlab)
    return fig


def plot_scalar(name, new, old=None):
    """Total cross sections: points with scale band per order."""
    ser = make_series(name, new, old)
    fig, ax = plt.subplots(figsize=(4.6, 4.2))
    for i, s in enumerate(ser):
        v, lo, hi, e = s["h"][0, 2], s["lo"][0], s["hi"][0], s["h"][0, 3]
        ax.errorbar([i], [v], yerr=[[v - lo], [hi - v]], fmt="o", color=s["color"], capsize=5, ms=4)
        ax.errorbar([i + 0.08], [v], yerr=[e], fmt="none", ecolor="k", elinewidth=0.8, capsize=2)
        ax.text(i, hi, f"{v:.4f}", ha="center", va="bottom", fontsize=7)
    ax.set_xticks(range(len(ser)))
    ax.set_xticklabels([s["label"] for s in ser])
    ax.set_ylabel(OBS[name][2])
    vals = np.array([[s["lo"][0], s["hi"][0]] for s in ser])
    pad = 0.15 * (vals.max() - vals.min())
    ax.set_ylim(vals.min() - pad, vals.max() + 2 * pad)
    ax.set_xlim(-0.6, len(ser) - 0.4)
    ax.set_title("LHC 13 TeV, " + title_text(name, old is not None).split("\n")[0].lower()
                 + " (bars: scale band; black: MC error)", fontsize=8)
    return fig


def safe(name):
    return re.sub(r"[^A-Za-z0-9]+", "_", name).strip("_")


def checks(new):
    """Cross-check the -nlo-ratio files against ratios of the main files."""
    for o in ("lo", "nlo", "nnlo"):
        for k in ("central", "min", "max"):
            rf = read_top(os.path.join(HERE, f"{o}-{k}-nlo-ratio.top"))
            worst = 0.0
            for n, h in new[(o, k)].items():
                den = new[("nlo", "central")][n][:, 2]
                m = den > 0
                worst = max(worst, np.max(np.abs(rf[n][m, 2] - h[m, 2] / den[m])))
            print(f"  ratio file {o}-{k}-nlo-ratio.top vs recomputed: max |diff| = {worst:.2e}")


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else ""
    if mode not in ("new", "overlap"):
        sys.exit(__doc__)
    new = load_new()
    old = load_old() if mode == "overlap" else None
    outdir = os.path.join(HERE, "plots-" + mode)
    os.makedirs(outdir, exist_ok=True)
    if mode == "new":
        checks(new)
    names = list(new[("nnlo", "central")].keys())
    if mode == "overlap":
        names = [n for n in names if OBS[n][0] is not None and OBS[n][0] in old["central"]]
        skipped = [n for n in new[("nnlo", "central")] if n not in names]
        print("  not in old, skipped:", ", ".join(skipped))
    with PdfPages(os.path.join(outdir, f"all-{mode}.pdf")) as pdf:
        for n in names:
            scalar = new[("nnlo", "central")][n].shape[0] == 1
            fig = plot_scalar(n, new, old) if scalar else plot_hist(n, new, old)
            fig.savefig(os.path.join(outdir, safe(n) + ".pdf"), bbox_inches="tight")
            fig.savefig(os.path.join(outdir, safe(n) + ".png"), dpi=110, bbox_inches="tight")
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)
            print("  plotted", n)


if __name__ == "__main__":
    main()
