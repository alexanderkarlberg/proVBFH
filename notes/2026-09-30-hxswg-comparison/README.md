# proVBFH-cs vs the HXSWG VBF study at 13.6 TeV (2026-09-30)

Goal: an apples-to-apples comparison of proVBFH-cs with the NNLO results
of the LHC Higgs WG VBF study (`~/cernbox/LHCHXSWG/vbf-higgs-wg`), for
which all run cards exist, and a timing comparison on thserv18. Runs are
in `proVBFH-cs/runs/hxswg-*` and `runs/oldcode-timing`.

## The study

- Draft: `Draft/vbf-higgs-wg.tex`. Set-up (Sec. "Computational set-up"):
  13.6 TeV, PDF4LHC21_40 (lhans 93100), G_mu scheme with MW 80.379,
  MZ 91.1876, GW 2.085, GZ 2.4952, Gmu 1.16638e-5, narrow Higgs,
  mu_0^2 = mH/2 sqrt(mH^2/4 + ptH^2) (`runningscales 1`), anti-kt R = 0.4.
  Cuts in `setup/binning.txt`: STXS (ptj > 30, |yH| < 2.5, mjj > 350 for
  the VBF bins), fiducial with |yj| < 4.7, mjj > 300, |dy_jj| > 2 and
  ptj > 20 or 30.
- proVBFH results: `proVBFH/results/{lo,nlo,nnlo}-{central,min,max}.top`
  (95 histograms; the NNLO file from `aux/combine_runs.f`). NNLO card and
  `vbfnlo.input`: `proVBFH/13.6TeV_NNLO/` (added by AK today): ncall1 1M x
  itmx1 3 for the grid, stage 2 with ncall2 5M x itmx2 3 in 9999 jobs.
  AK will add the raw per-seed files in the same directory.
- Code: the study's copy `proVBFH/src` is identical to the proVBFH tag
  v2.1.0 (apart from an extra copy of `setlocalscales.f`). v2.1.0 uses
  Hoppet 1.3 (Fortran module `hoppet_v1`); the port to Hoppet 2 came after
  v2.1.0. Its exclusive part is the same code as the current proVBFH, so
  it has the three issues of the POWHEG bug report (kl index, nf 4/5, NC
  gg pair type). Earlier today I called the study copy "2023 code" and
  listed its differences with the current proVBFH as if they were
  unexplained; they are the Hoppet-2 port, the NF-correction rewrite and
  the fixed-scale default (commits after v2.1.0).
- Analysis: `proVBFH/analysis/13.6TeV_analysis.f` of the study (the
  NNLO file has its extra STXS histograms). The copy in the proVBFH repo
  is an older version (no STXS fills, early return). Copied without
  changes to `proVBFH-cs/analysis/hxswg136_analysis.f`; `make
  ANALYSIS=hxswg136` builds `proVBFH-cs-hxswg136`.
- Generation cuts: the analysis and the generation cuts use the same
  `setup_vbf_cuts`/`buildjets`/`vbfcuts` and common block; the analysis
  sets `deltay_jjmin` 2 and `yjetmax` 4.7 for the fiducial part and then
  resets them to 0 and 4.7d10. The NNLO card has those values anyway
  (ptalljetmin 20, ptjetmin 20, mjjmin 300, no dy or y cut), which are
  looser than every histogram. proVBFH-cs skips a point only if neither
  an event nor its Born counterevent passes (`cs_passes`), so these cuts
  are safe.
- The draft's definitions of fiducial (a) and (b) are swapped relative to
  its own table: its equations give (a) ptj > 30, (b) ptj > 20, but the
  (a) table (LO 2.479 pb) is ptj > 20 (`setup/binning.txt`, spreadsheet
  set-up 2a; set-up 2b, ptj > 30, is 2.022 pb). Reported to AK.

## LO check (`runs/hxswg-lo`, inclusive part at LO, 4 seeds x 6.6M points)

| | proVBFH-cs [pb] | study [pb] | pull |
|---|---|---|---|
| sig incl cuts (ptj > 20) | 2.47789 +- 0.00124 | 2.47899 +- 0.00010 | -0.88 |
| sig incl cuts (ptj > 30) | 2.02175 +- 0.00112 | 2.02213 +- 0.00009 | -0.34 |

STXS and 2D bins: chi2 about 1 per bin. ptH-STXS and ptj1-STXS have
chi2/n 220/100 and 172/99, from bins above 1.3 TeV (1e-8 pb, a handful
of events in 4 seeds, 30-55% low with underestimated errors). Integrated
over pt ranges (0-200, 200-500, 500-1000, >1000 GeV) all agree to 1%,
|pull| <= 2.3. The inputs match.

## combine_runs.f (used for the study and for 1506.02660)

- Per bin, runs outside median +- 10 x (half the central 68% range) are
  dropped (column 5 = fraction kept, about 0.99 for the study's NNLO).
  AK: the trimming is there because of huge outliers that destroy the
  variance ("probably due to the singularity you found").
- Off by one: the sum runs from `first_index` instead of
  `first_index+1`. When there are low outliers, the largest of them is
  added back (divided by the number kept); when there are none, element 0
  (out of bounds) is read. Synthetic test (38 values near 1, one at -100,
  one at +100): result -1.630 instead of 1.0014 (`scratchpad/combtest`).
  Expected effect on the study: small, the value added back is the least
  extreme low outlier, divided by about 9900.
- In bins where all runs are zero, both indices stay 0 and the result is
  element 0 divided by 0: 1506.02660's `11.top` has -inf +- inf (yH in
  [-5,-4.5]) and -4.3e-58 +- 4.3e-58. `combine_parts.py` now skips
  non-finite reference bins.
- Errors: sqrt(sum of the VEGAS errors^2 of the kept runs)/N.

## Timing (approved by AK)

- The study's binary (Feb 2024) aborts with `std::length_error` in
  FastJet: built against other FastJet headers than today's
  `/lib64/libfastjet.so`. For the timing the current proVBFH (exclusive
  part identical to v2.1.0) was built in a scratch copy with the study's
  analysis and FastJet headers from /usr/include:
  `runs/oldcode-timing/proVBFH-old.bin` (see `.info`). A 20k-point
  stage 1 and stage 2 run cleanly.
- Stage-1 grid with the NNLO card (1M x 3), one job, thA371a:
  `runs/oldcode-timing/grid`.
- `runs/hxswg-timing/run_thserv18.sh` (thserv18, nice 10, 16 cores)
  starts when nnlo-full and the grid are done: old code stage 2 with
  200k and 1M points (2 jobs each; the difference gives the time per
  point), proVBFH-cs NNLO exclusive part (8 seeds x 2.2M points) and NLO
  exclusive part (4 seeds x 6.6M points), all at the same time.

## 1506.02660 results

`~/cernbox/proVBFH-github/1506.02660.tgz` (from AK): `HH.top`, `11.top`,
`22.top` = mu_R = mu_F = {0.5, 1, 2} x mu_0, unpacked in
`runs/ref-1506.02660`. Central sig(all VBF cuts 2 jets) = 0.84383 +-
0.00046 pb (combined with the same trimming, fraction kept 0.988). The
ten VBF-cut distributions have the same binning as ours (names without
`-vbf`: `combine_parts.py --strip=-vbf`). The comparison with nnlo-full
+ nnlo-incl runs when nnlo-full is done.
