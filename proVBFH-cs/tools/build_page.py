#!/usr/bin/env python3
"""Build the production status/results page (an HTML artifact) from the
combined results.

Usage:
  build_page.py --data /ptmp/.../combined --runs /ptmp/.../cs-production
                --style old_index.html --notes findings.html --out page/

--data: the output root of combine_and_plot.sh (<setup>/{summary-*.json,
plots-*/}, <setup>/plots-orders/). --style: an HTML file whose <style>
block is reused. --notes: an HTML fragment with the current findings,
inserted as is. Writes <out>/index.html, copies the plots it shows to
<out>/plots/..., and writes <out>/files.json (the map for the publish).
PNGs for the 1506 plots and the HXSWG NNLO comparison; multipage PDFs for
the rest (the publish limit is 255 files).
"""
import argparse
import datetime
import glob
import html
import json
import os
import re
import shutil

SETUPS = {"p1506": "1506.02660 (13 TeV, VBF cuts)", "hxswg136": "HXSWG 13.6 TeV"}
PARTS = [  # label, run dir relative to --runs
    ("p1506 NNLO exclusive", "prod-nnlo/p1506/excl"),
    ("p1506 NNLO inclusive", "prod-nnlo/p1506/incl"),
    ("p1506 NLO exclusive", "prod-lonlo/p1506/nlo-excl"),
    ("p1506 NLO inclusive", "prod-lonlo/p1506/nlo-incl"),
    ("p1506 LO inclusive", "prod-lonlo/p1506/lo-incl"),
    ("hxswg136 NNLO exclusive", "prod-nnlo/hxswg136/excl"),
    ("hxswg136 NNLO inclusive", "prod-nnlo/hxswg136/incl"),
    ("hxswg136 NLO exclusive", "prod-lonlo/hxswg136/nlo-excl"),
    ("hxswg136 NLO inclusive", "prod-lonlo/hxswg136/nlo-incl"),
    ("hxswg136 LO inclusive", "prod-lonlo/hxswg136/lo-incl"),
]
SIGMA = [  # p1506 cross sections for the key-number table
    ("σ(VBF cuts, ≥ 2 jets)", "sig(all VBF cuts 2 jets)"),
    ("σ(VBF cuts, ≥ 3 jets)", "sig(all VBF cuts 3 jets)"),
    ("σ(VBF cuts, ≥ 4 jets)", "sig(all VBF cuts 4 jets)"),
]


def count(runs, rel):
    d = os.path.join(runs, rel)
    lst = os.path.join(d, "jobs.list")
    n = sum(1 for _ in open(lst)) if os.path.exists(lst) else 0
    done = len(glob.glob(os.path.join(d, "job-*", "done")))
    return done, n


def fig(src, name, chi=None, alts=None):
    cap = f"<span>χ² {chi[0]:.1f}/{chi[1]}</span>" if chi else ""
    if chi and alts:
        parts = [f"study {chi[0]:.0f}"] + [f"{x['label'].split(',')[-1].strip()} {x['chi2']:.0f}" for x in alts]
        cap = f"<span>χ²/{chi[1]}: {', '.join(parts)}</span>"
    return (f'<figure><img src="{src}" alt="{html.escape(name)}" loading="lazy">'
            f'<figcaption><b>{html.escape(name)}</b>{cap}</figcaption></figure>')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--runs", required=True)
    ap.add_argument("--style", required=True)
    ap.add_argument("--notes", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--static", help="HTML fragment (checks, pilot) inserted after the status")
    a = ap.parse_args()

    style = re.search(r"<style>.*?</style>", open(a.style).read(), re.S).group(0)
    files = {}
    os.makedirs(a.out, exist_ok=True)

    def copy(src, rel):
        dst = os.path.join(a.out, rel)
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        shutil.copy(src, dst)
        files[rel] = rel
        return rel

    now = datetime.datetime.now().strftime("%-d %b %Y, %H:%M")
    rows = []
    for lab, rel in PARTS:
        done, n = count(a.runs, rel)
        if n == 0:
            continue
        pct = 100.0 * done / n
        rows.append(f'<tr><td>{lab}</td><td class="num">{done:,}</td><td class="num">{n:,}</td>'
                    f'<td><div class="bar"><span style="width:{pct:.1f}%"></span></div></td>'
                    f'<td class="num">{pct:.0f}%</td></tr>')
    status = ('<div class="tablebox"><table><thead><tr><th>Part</th><th class="num">Done</th>'
              '<th class="num">Jobs</th><th>Progress</th><th class="num"></th></tr></thead><tbody>'
              + "".join(rows) + "</tbody></table></div>")

    sections = []
    # key numbers, p1506 NNLO
    sj = os.path.join(a.data, "p1506", "summary-nnlo.json")
    if os.path.exists(sj):
        s = {h["name"]: h for h in json.load(open(sj))}
        trs = []
        for lab, key in SIGMA:
            if key not in s:
                continue
            h = s[key]
            f = 1000.0  # pb -> fb
            trs.append(f'<tr><td>{lab}</td><td class="num">{h["new"][0]*f:.2f} ± {h["new_err"][0]*f:.2f}</td>'
                       f'<td class="num">{h["new_min"][0]*f:.1f} – {h["new_max"][0]*f:.1f}</td>'
                       f'<td class="num">{h["old"][0]*f:.2f} ± {h["old_err"][0]*f:.2f}</td>'
                       f'<td class="num">{h["new"][0]/h["old"][0]:.4f}</td></tr>')
        a2, a3 = s.get(SIGMA[0][1]), s.get(SIGMA[1][1])
        if a2 and a3:
            n2 = (a2["new"][0] - a3["new"][0]) * 1000
            o2 = (a2["old"][0] - a3["old"][0]) * 1000
            trs.insert(2, f'<tr><td>σ(VBF cuts, exactly 2 jets)</td><td class="num">{n2:.2f}</td><td></td>'
                          f'<td class="num">{o2:.2f}</td><td class="num">{n2/o2:.4f}</td></tr>')
        sections.append('<h2>1506.02660 at NNLO: cross sections</h2><div class="tablebox"><table><thead><tr>'
                        '<th>[fb]</th><th class="num">proVBFH-cs</th><th class="num">scale band</th>'
                        '<th class="num">old proVBFH</th><th class="num">new / old</th></tr></thead><tbody>'
                        + "".join(trs) + "</tbody></table></div>")
    sections.append(open(a.notes).read())

    galleries = []
    legend = ('<div class="legend col"><span><i style="background:var(--new)"></i>proVBFH-cs, with scale band '
              '(ξ = ½, 1, 2)</span><span><i style="background:var(--old)"></i>old proVBFH with its band</span>'
              '<span>Lower panels: ratio to the old central value.</span></div>')
    for setup, title in SETUPS.items():
        # comparison with the old code, NNLO
        sj = os.path.join(a.data, setup, "summary-nnlo.json")
        if os.path.exists(sj):
            summ = json.load(open(sj))
            figs = []
            for h in summ:
                if h["name"] == "sig incl cuts":  # 1506 only: not comparable (see the notes)
                    continue
                pngs = glob.glob(os.path.join(a.data, setup, "plots-nnlo", f"{h['index']:03d}-*.png"))
                if not pngs:
                    continue
                rel = copy(pngs[0], f"plots/{setup}/nnlo/{os.path.basename(pngs[0])}")
                figs.append(fig(rel, h["name"], (h["chi2"], h["nbins"]), h.get("alts")))
            leg = legend
            if any("alts" in h for h in summ):
                leg = legend.replace('old proVBFH with its band', 'old proVBFH, the study\'s trimmed merge, with its band')
                leg = leg.replace('</div>', '<span><i style="background:#d62728"></i>old, plain (untrimmed) merge, '
                                  'dashed</span><span><i style="background:#2ca02c"></i>old, symmetric 0.5% trim, '
                                  'dash-dotted (bands dotted)</span></div>')
            galleries.append(f'<section><h2 class="col">{title}: NNLO, new vs old</h2>{leg}'
                             f'<div class="grid">{"".join(figs)}</div></section>')
        # LO and NLO comparisons with the old code (HXSWG): PDFs
        links = []
        for o in ("lo", "nlo"):
            pdf = os.path.join(a.data, setup, f"plots-{o}", "all.pdf")
            if os.path.exists(pdf):
                rel = copy(pdf, f"plots/{setup}/{o}-vs-old.pdf")
                links.append(f'<li><a href="{rel}">{o.upper()}, new vs old (PDF, one page per histogram)</a></li>')
        # orders
        od = os.path.join(a.data, setup, "plots-orders")
        if os.path.isdir(od):
            if setup == "p1506":
                figs = []
                for png in sorted(glob.glob(os.path.join(od, "*.png"))):
                    if "sig_incl_cuts" in png:
                        continue
                    name = os.path.basename(png)[4:-4]
                    rel = copy(png, f"plots/{setup}/orders/{os.path.basename(png)}")
                    figs.append(fig(rel, name))
                galleries.append(f'<section><h2 class="col">{title}: LO, NLO and NNLO</h2>'
                                 f'<p class="col">proVBFH-cs only, each order with its scale band; '
                                 f'lower panels: ratio to NLO.</p><div class="grid">{"".join(figs)}</div></section>')
            else:
                rel = copy(os.path.join(od, "all.pdf"), f"plots/{setup}/orders.pdf")
                links.append(f'<li><a href="{rel}">LO, NLO and NNLO of proVBFH-cs (PDF)</a></li>')
        if links:
            galleries.append(f'<section class="col"><h2>{title}: more plots</h2><ul>{"".join(links)}</ul></section>')

    # NLO H+3j cross-check (VBFNLO), ratio to proVBFH-cs
    for sub, title, txt in (
            ("plots-h3j-nlo", "H+3j cross-check: O(αs²) ≥ 3 jets, VBFNLO 3.0 vs proVBFH-cs",
             "VBFNLO 3.0 NLO H+3j (200 seeds, patched to fill the p1506 observables; single scale μ0) and the old "
             "proVBFH, as ratios to proVBFH-cs NNLO (its ≥ 3-jet bins are O(αs²), i.e. NLO H+3j). The 4-jet "
             "histograms are tree-level H+4j in all codes."),
            ("plots-h3j-lo", "H+3j cross-check: tree level, VBFNLO 3.0 vs proVBFH-cs",
             "VBFNLO 3.0 LO H+3j (50 seeds) against the ≥ 3-jet bins of the proVBFH-cs NLO run, which are "
             "tree-level H+3j. Ratio panels 0.95-1.05.")):
        d = os.path.join(a.data, "p1506", "vbfnlo", sub)
        if not os.path.isdir(d):
            continue
        figs = []
        for png in sorted(glob.glob(os.path.join(d, "*.png"))):
            if sub.endswith("-lo") and ("4_jets" in png or "j4" in png):  # no 4 jets at tree-level H+3j
                continue
            rel = copy(png, f"plots/p1506/{sub}/{os.path.basename(png)}")
            figs.append(fig(rel, os.path.basename(png)[4:-4]))
        galleries.append(f'<section><h2 class="col">{title}</h2><p class="col">{txt}</p>'
                         f'<div class="grid">{"".join(figs)}</div></section>')

    static = open(a.static).read() if a.static else ""
    page = f"""<title>proVBFH-cs NNLO Production</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Literata:opsz,wght@7..72,500;7..72,650&family=Public+Sans:wght@400;600&family=JetBrains+Mono:wght@400;500&display=swap">
{style}
<div class="wrap">
<header class="col">
  <h1>proVBFH-cs NNLO production</h1>
  <p class="sub">VBF Higgs with line-by-line projection-to-Born, 1506.02660 and HXSWG 13.6 TeV set-ups, scale points μ<sub>R</sub> = μ<sub>F</sub> = ξ μ<sub>0</sub>(p<sub>T,H</sub>), ξ = 1, ½, 2. MPCDF cluster, partition <code>alma</code>. Updated {now}.</p>
</header>
<section class="col">
  <h2>Where the runs stand</h2>
  {status}
  {static}
  {"".join(sections)}
</section>
{"".join(galleries)}
<footer class="col"><p>Errors are from the seed scatter of the finished jobs, no trimming. Raw per-job outputs stay in <code>/ptmp/mpp/akarlber/cs-production</code>. Log: <code>notes/2026-10-cs-production/README.md</code>.</p></footer>
</div>
"""
    open(os.path.join(a.out, "index.html"), "w").write(page)
    json.dump(files, open(os.path.join(a.out, "files.json"), "w"))
    print(f"{len(files)} files -> {a.out}")


if __name__ == "__main__":
    main()
