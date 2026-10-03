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

## The raw per-seed data (AK, 16:30): 9940 of the 9999 stage-2 seeds

`proVBFH/13.6TeV_NNLO/{HH,11,22}.tgz`, each 9940 files `pwg-NNNN-NNLO.top`
(95 histograms, 988 bins). Unpacked outside cernbox; loader and an exact
re-implementation of `combine_runs.f` in `tools/rawload.py` (set RAWDIR).

- **Reproduction:** `combine_runs` (limit 10, with the off-by-one) on
  the 9940 central files gives `results/nnlo-central.top` in all 946
  non-empty bins to 5e-8 (the precision of the files).
- **sig incl cuts (ptj > 20):** plain mean 2.09614 +- 0.00649 (seed
  scatter); study (trimmed) 2.08980 +- 0.00075; trimmed without the
  off-by-one 2.08964. Scatter of the seeds 0.65 pb against a VEGAS error
  of 0.072 per seed. 94 seeds dropped (41 low, 53 high). ptj > 30:
  plain 1.65531 +- 0.00671, trimmed 1.65130 +- 0.00068.
- **All bins:** plain - trimmed is positive in 59% of the bins, median
  +0.5% (+0.43 plain errors); the quoted errors are a median 4.2 times
  smaller than the plain seed-scatter errors.
- **Tails:** power laws on both sides, Hill index alpha = 1.0-1.6 (top
  20-100 seeds): infinite variance, and a mean that barely exists. The
  largest seed (61.3 pb, median 2.09) moves the plain mean by 0.3%, as
  much as the trimming; the running plain mean wanders by 0.5%.
- **The spikes are single events:** all of a big seed's excess sits in
  one bin of every histogram: 2 jets, mjj 300-700 GeV, dy_jj 2-5, ptH
  80-120 GeV, ptHjj < 20 GeV (the third parton soft or collinear), both
  signs (+59, +15, -11 pb). An unsubtracted singular region, as expected
  from issue 2 (missing ISR region) and the NC gg pair type.
- **Consequences:** the study's (and 1506.02660's, combined the same way)
  quoted errors do not include the tails, and the trimmed mean is biased
  by an unknown amount (here -0.3% relative to the plain mean). The
  off-by-one is negligible (+0.00016 pb, 0.2 quoted errors).

## Timing on thserv18 (15:30-18:20, 16 jobs at the same time, nice 10)

| run | points | CPU [s] | per point |
|---|---|---|---|
| old code stage 2, study card, study grid | 200k | 2034, 2038 | |
| old code stage 2 | 1M | 10064, 10023 | 10.0 ms (difference), 34 s set-up |
| proVBFH-cs NNLO exclusive (8 seeds) | 2.2M | 8215-8416 | 3.8 ms |
| proVBFH-cs NLO exclusive (4 seeds) | 6.6M | 2306-2466 | 0.36 ms |

- Old code, one study job (ncall2 5M x itmx2 3 = 15M points): 34 + 15M x
  10.0 ms = 1.50e5 s = 42 CPU-h; 9999 jobs: 1.50e9 s = 417,000 CPU-h
  (thserv18 equivalent; stage 1 not counted). The old code runs with
  `testplots 1` (the analysis is called for every real flavour region).
- proVBFH-cs, 8 NNLO seeds (66,580 CPU-s): per bin, the CPU needed to
  reach the study's errors, C = 66,580 s x (err_ours/err_study)^2:
  - quoted (trimmed) errors: median over 946 bins 3,600 CPU-h, speed-up
    median 116 (16-84%: 7-486), total bin (ptj > 20) 139;
  - plain seed-scatter errors: median 170 CPU-h, speed-up median 2,500
    (16-84%: 130-35,000), total bin 10,500.
  - With 8 seeds the errors are uncertain by about 25% per bin (50% in
    CPU); the inclusive part (a few CPU-min) is not included.

## 13.6 TeV NLO and NNLO from the timing runs (18:30)

Inclusive part `runs/hxswg-incl{2,3}` (4 seeds each, thA371a) + exclusive
part from `runs/hxswg-timing` (NLO 4 seeds, NNLO 8 seeds); `--error max`.

| | proVBFH-cs [pb] | study [pb] | pull |
|---|---|---|---|
| NLO, ptj > 20 | 2.18262 +- 0.00313 | 2.17937 +- 0.00026 (`nlo-central`) | +1.0 |
| NLO, ptj > 30 | 1.74104 +- 0.00288 | 1.73978 +- 0.00017 | +0.4 |
| NNLO, ptj > 20 | 2.0977 +- 0.0096 | 2.0898 (trimmed) / 2.0961 +- 0.0065 (plain) | +0.8 / +0.1 |
| NNLO, ptj > 30 | 1.6507 +- 0.0115 | 1.6513 / 1.6553 +- 0.0067 | -0.05 / -0.35 |

(The spreadsheet's NLO 2176.39 fb is an older run; `nlo-central.top`,
July 2024, is 2179.37 fb.) With 8 NNLO seeds the NNLO check is at the
0.5% level only.

Distributions: about 1 per bin except in the far tails: NLO ptj2-STXS
chi2 4234/99 and ptHjj 759/50, from bins with ptj2 > 540 GeV or ptHjj >
880 GeV. There the exclusive part dominates (hard real emission makes the
second jet; the study's K factor is 1.2 at ptj2 400-600, 1.8 at 600-1000,
3.0 at 1-2 TeV), and our exclusive part is badly under-sampled with 4
seeds (errors 40-80%, and underestimated): ours 1.5, 1.6, 0.83 for the
same K factors, while our inclusive part agrees with the study's LO to
0.1-0.5% there. VEGAS adapts to sum |w|, dominated by the bulk, and the
radiation variables are sampled logarithmically towards the soft and
collinear limits, so hard emissions get few points. Not a physics
difference as far as can be told, but a sampling weakness to fix (hard-
radiation sampling in the warm-up / channel study) before production, and
the old code's radiation sampling is better per point in these bins.

## Hard-emission sampling (AK: "improve the hard-emission sampling while we wait", 18:45)

Cause of the under-sampled tails: a line's two partons have pT^2 =
Q^2 z (1-z) (1-xp)/xp, so a 500 GeV jet at Q ~ MW needs xp of a few
1e-2, near the lower end of [xB, 1]; the logarithmic map in 1-xp puts
about 0.3% of the points at xp < 0.04, and about 6% have z in [0.25,
0.75].

Change (commits 04525fd, 93ab2e9; default off, bitwise unchanged):
`cs_hardfrac h` adds a second channel to `line_radiation` with
probability h: ln xp uniform in [ln xB, 0], z uniform; weight 1/g with
g = (1-h) g_log(xp) g_log(z) + h g_hard(xp) g_hard(z). The same in the
first step of `gen_four` and in `four_weight` (module variable
`four_hard`). Checks:
- `tests/test_kinematics` now covers npow 2, logarithmic, and logarithmic
  with h = 0.3 (pointwise and phase-space volume). The pointwise checks
  failed at 3e-9 with h = 0.3 against a tolerance of 1e-10 relative to
  Q^2: round-off, as the hard channel reaches xp ~ xB where pin = pB/xp
  is large compared with sqrt(Q^2) (relative to E_max^2 the errors are
  1e-15 for both samplings); the mass, xp and z checks are now relative
  to kappa = E_max^2/Q^2.
- `tests/test_four` with four_hard = 0 and 0.3: weight of gen_four =
  four_weight exactly, integrals of the test functions against the flat
  RAMBO reference within |pull| < 1.5.
- With the option off, the NNLO smoke run (incl. (2,0)) is bitwise
  identical to the earlier build.

Tests running (13.6 TeV study set-up): NLO exclusive h = 0.2, 0.5 (4
seeds, as `hxswg-timing/nlo-s*`, thserv22), NLO h = 0 and 0.3 with 12
seeds x 1.65M points each (thserv19; more seeds for stable tail errors),
NNLO exclusive h = 0.3 (8 seeds, as `hxswg-timing/nnlo-s*`, thserv21).
Compared with `tools/sampling_compare.py` (per-bin and integrated-range
gains in error^2 x CPU, and consistency).

### NLO result: h = 0.3 vs h = 0 (12 seeds x 1.65M each, thserv19, 19:35)

`hxswg-hard/nlo12/compare.txt`, gains in error^2 x CPU (same CPU,
2.6 CPU-h each):
- fiducial totals (NLO exclusive part): 2.3 (ptj > 20), 2.7 (ptj > 30);
- ptj2 0-200: 2.9; 600-1000: 4.7; ptHjj 100-1000: 3.7;
- ptj2 1000-2000: h = 0 finds almost nothing (-2.5e-8 +- 1.2e-8), h = 0.3
  gives +2.5e-7 +- 0.6e-7 (the study implies about +4e-7 for the
  exclusive part there);
- worse: mjj 3-5 TeV 0.8, ptH > 500 GeV 0.15 (Born-level high pT with
  soft/collinear cancellations; 30% fewer points in the logarithmic
  channel gives 0.7);
- median over all bins 1.27 (seed scatter), 1.22 (VEGAS errors).
- Consistency: chi2 1201/934 between the variants, all of the excess in
  the hard-tail histograms (ptj2-STXS 285/99, ptHjj 88/50 each, ptj1-STXS
  120/99), where h = 0 misses rare configurations and underestimates its
  errors; the other histograms 410/499.

### Second emission of gen_four (19:45, commit 7dc122a)

The 4-jet histograms at the 1506.02660 set-up (chi2 107/34; most bins
10-30% low with small errors, a few spikes) show that the second
emission also needs a hard channel: a 4th jet above 25 GeV needs it hard,
and z of the second step is logarithmic towards both ends (about 6% of
the points at z in [0.25, 0.75]). With probability four_hard the second
step now takes y uniform (FF) or ln x uniform in [ln xi3, 0] (FI) and z
uniform; four_weight uses the combined density; the cutoff checks are
applied in the hard branch only (the logarithmic maps cannot go below
the cutoff), so the default stays bitwise identical (checked).
test_four with four_hard = 0.3: weight = four_weight exactly, integrals
|pull| < 2.6 (reruns 1.0, 1.9). With h = 0.3 fewer points fail all cuts
(12,346 vs 15,755 of 25k in the smoke run), so the CPU per point rises
by about 36%; the gains are per CPU.

NNLO tests with the full channel (h = 0.3), started 19:45: 13.6 TeV set-up
(8 seeds as hxswg-timing/nnlo-s*, thserv19) and 1506.02660 set-up with the
paper's analysis (16 seeds, iseed 7301-7316, thserv18; against the 30
h = 0 seeds of nnlo-p1506, including the 3- and 4-jet rates). The
h = 0.2/0.5 NLO test on thserv22 runs at about half speed (the machine
is loaded to 75 by others), so its CPU times are not comparable with the
h = 0 run on thserv18; the 12 + 12 seed test on thserv19 is the clean one.

### NLO h = 0.2, 0.5 (4 seeds, thserv22 at half speed) and a bias check (20:20)

- h = 0.2 and 0.5 against h = 0 (4 seeds each; the CPU of the h > 0 runs
  is inflated by about 40% by the load on thserv22): hard tails improve
  as with h = 0.3 (ptHjj 100-1000 gain 6-9, mjj 3-5 TeV 3-11, ptj2
  400-600 37-88, ptj2 1-2 TeV found); median gains 0.99 and 0.92 before
  the CPU correction.
- NLO exclusive total (ptj > 20) against the value the study's NLO
  implies with our inclusive part (-0.20441 +- 0.00128): h = 0 +1.6
  sigma, h = 0.2 +1.9, h = 0.3 (12 seeds) -1.4, h = 0.5 -1.1 (ptj > 30:
  +0.8, +1.4, -2.2, -1.3). No sign of a bias of the hard channel, but the
  4-seed errors are not reliable.
- Exact check started: the line-NLO validation (cs_order 11, no cuts,
  30 seeds x 3.3M, as `stage3-validnlo`: -0.102618 +- 0.00033 against
  the structure functions' -0.102623 +- 0.000068) with h = 0.3:
  `runs/stage3-validnlo-h03`, thserv09.
- **Bias check passed (20:55):** line NLO (cs_order 11, no cuts, 30 seeds
  x 3.3M) with h = 0.3: -0.102827 +- 0.000297 pb against the structure
  functions' -0.102623 +- 0.000068 (-0.67 sigma); h = 0 gave -0.102618
  +- 0.000330 (+0.02 sigma). The hard channel of the line radiation is
  unbiased at the 0.3% level of the full inclusive integral; error^2 x
  CPU 3.88e-3 against 4.77e-3 (gain 1.23 for this total without cuts).
  (The watcher's parser failed on VEGAS's "integral =-0.99E-01+/-"
  without a space; parsed from the "sum |w1|+|w2|" line instead.)

### NNLO tests at 13.6 TeV (23:35) and a separate fraction for the second emission

`runs/hxswg-hard/compare-nnlo.txt`, 8 seeds each against h = 0
(`hxswg-timing/nnlo-s*`, 18.5 CPU-h); gains in error^2 x CPU:

| | first emission only (h = 0.3; 20.3 CPU-h) | both emissions (h = 0.3; 24.0 CPU-h) |
|---|---|---|
| total, ptj > 20 | 1.8 | 0.28 |
| total, ptj > 30 | 0.7 | 0.42 |
| ptj2 400-600 / 600-1000 | 4.1 / 7.5 | 9.3 / 3.3 |
| ptHjj 100-1000 | 41 | 13 |
| median over 946 bins | 0.73 | 0.62 |

At NNLO the variance comes mostly from the double-unresolved corners of
(2,0)+(0,2); a hard channel for the second emission takes points from
there and costs a factor 2.5-3.5 on the totals. The first emission's
channel helps the tails and is neutral on the totals (the 8-seed factors
are uncertain by about 50%). Both variants agree with h = 0 within the
errors. The second step now has its own fraction, `cs_hardfrac2`
(default 0; `four_hard2`); `cs_hardfrac` alone is bitwise identical to
the first-emission-only build (checked), and the default to the original.
test_four covers (0, 0), (0.3, 0), (0.3, 0.3). The 1506.02660 run with
both emissions at 0.3 (`nnlo-p1506-h03`) will show whether the 4-jet
observables need a small second-step fraction.

### 1506.02660 set-up with both emissions at h = 0.3 (23:47)

`nnlo-p1506-h03` (16 seeds, 62.6 CPU-h) against nnlo-p1506 (30 seeds,
h = 0, 84.0 CPU-h): gains 2.06 (2-jet total), 1.53 (3-jet), 5.5 (4-jet),
yj4 5.7, min{rap(j1,j3),rap(j3,j2)} 2.4; median over 351 bins 1.12;
chi2 378/351 between the two. So with the tighter 1506.02660 cuts the
second emission's channel helps everything, while with the study's
looser cuts it cost 2.5-3.5 on the totals (8 seeds): the best
cs_hardfrac2 depends on the set-up. Overnight tuning (`runs/hard-tune`,
16 jobs each on thserv09, 21 (13.6 TeV) and 19, 22 (1506.02660), mixed
variants per machine so the CPU is comparable): 13.6 TeV (0,0) 8 seeds,
(0.3,0) 8, (0.3,0.1) 16 (plus the earlier 8+8); 1506.02660 (0.3,0) 16,
(0.3,0.1) 16.

### Tuning results (2026-10-01, 04:20) and an open question

`runs/hard-tune` (overnight), gains in error^2 x CPU against h = 0:
- 13.6 TeV, 16 seeds per variant (new 8 + earlier 8 for A and B):
  (0.3, 0): totals 1.8 / 0.8, ptj2 400-600 / 600-1000 3.4 / 5.1, ptHjj
  100-1000 4.8, median 0.80; (0.3, 0.1): totals 0.8 / 0.5, median 0.75.
- 1506.02660 set-up, 16 seeds per variant against nnlo-p1506 (h = 0, 30):
  (0.3, 0): 2-jet 1.1, 3-jet 3.7, 4-jet 39, median 1.07; (0.3, 0.1): 1.2,
  1.7, 11, 1.18; (0.3, 0.3): 2.1, 1.5, 5.5, 1.12.
- So the second emission's channel does not pay; the first emission's
  does (the 4-jet rate even more than with the second-step channel).
  Recommended: cs_hardfrac 0.3, cs_hardfrac2 0.

**Open: ptHjj 100-1000 at 13.6 TeV differs between the samplings.**
NNLO, ptHjj-ptj20 integrated over [100, 1000]: h = 0 (A16) 9.23e-3 +-
0.43e-3 pb, (0.3, 0) (B16) 7.53e-3 +- 0.19e-3: -3.7 sigma (and -3.3 and
-1.8 in the two independent 8-seed comparisons). Per seed, h = 0 spreads
over 6.5-12.0 (median 9.7) and h = 0.3 clusters at 7.6 +- 0.5: not a
single outlier. NLO parts agree (3.81e-3 vs 3.77e-3), so the difference
is in the O(alpha_s^2) pieces (about 5.4e-3 vs 3.7e-3). The study gives
8.41e-3 (plain), between the two (its NLO 3.78e-3 agrees with ours; its
>= 3-jet region has the missing-ISR-region excess). Correction: the
statement above that the NNLO variants "agree with h = 0 within the
errors" holds for the totals but not for this tail.
K+P convolutions use their own Gauss quadrature (not the sampling). To
locate it: `runs/hard-diag`, (2,0)+(0,2) only (cs_order 2, cs_only2) and
(1,1) only (cs_order 3, cs_only2, cs_no20), h = 0 and 0.3, 16 seeds
each, 13.6 TeV set-up (thserv09, 21, 19, 22).

### The open question, per piece (2026-10-01, 07:00)

`runs/hard-diag`, 16 seeds per piece and sampling, 13.6 TeV set-up;
ptHjj-ptj20 over [100, 1000] in 1e-3 pb:

| | h = 0 | h = 0.3 | pull |
|---|---|---|---|
| (2,0)+(0,2) only | 3.75 +- 0.32 | 3.46 +- 0.15 | -0.8 |
| (1,1) only | 0.22 +- 0.11 | 0.24 +- 0.05 | +0.1 |
| NLO part (nlo12) | 3.81 +- 0.05 | 3.77 +- 0.03 | -0.8 |
| sum | 7.78 +- 0.34 | 7.47 +- 0.16 | |
| full NNLO runs | 9.23 +- 0.43 (A16) | 7.53 +- 0.19 (B16) | |

Each piece agrees between the samplings, and the sum of the h = 0 pieces
agrees with h = 0.3. The outlier is the full NNLO run at h = 0, 2.7
sigma above the sum of its own pieces: its seeds spread over 6.5-12; the
full run with (0.3, 0.1) also has one large seed (16.0). In this tail the
full-run estimates have heavy-tailed seed distributions and the 16-seed
errors are not reliable; no sign of a bias of the hard channel. Gains of
the channel per piece: (2,0)+(0,2) totals 3.7 / 4.0, ptHjj ranges 4-7;
(1,1) totals 2.6 / 3.6, ranges 2.3-6.

(2026-10-01, before the first push: commit 3e57a64 had swept in five
untracked files with `git add notes` - the build products
`tools/find_regions.o` and `tools/list_regions` of the stage-2 notes and
three run outputs of AK's `notes/scale-setting-in-vbf-hh`. The branch was
rebuilt without them (the files stay on disk, untracked); the commits from
there on have new hashes, e.g. f215582 -> 7dc122a, daf218b -> c785573.)

## Direct merge of all three scales (3 Oct)

AK: "could you also do a direct merge of the HXSWG files and commit and
push those? ... to try and estimate the bias from the trimming". Done with
`tools/merge_plain.py`, output in
`proVBFH-cs/production/reference/hxswg136/plain/` (README there):
- plain means per scale with seed-scatter errors;
- the central/min/max band.
The trimmed re-combination reproduces the study's
`nnlo-{central,min,max}.top` to 5e-8, which confirms that min/max are the
per-bin envelope of HH, 11, 22.

σ(ptj > 20), plain − trimmed:
- HH: +0.0082 pb (+0.40%);
- 11: +0.0063 pb (+0.30%);
- 22: +0.0048 pb (+0.23%).
One HH seed file (pwg-7524) is empty and is skipped.

## A more robust third reference (3 Oct)

AK: "Can you come up with a more robust way of trimming them as well as a
third reference?" Estimators compared on the seed distributions (central
scale): plain mean, the study's trimming, median-of-means, symmetric
quantile trimming, winsorising, Peng's tail-corrected mean.
- Peng: unstable for Hill index ≈ 1.
- Median-of-means: drifts with the number of groups in skewed bins.
- Symmetric 0.5% trim with bootstrap errors: stable against the trim
  fraction (0.25%/1%/2%: median shifts 0.27/0.28/0.49 errors), errors
  1.2× the study's. Chosen.
Added to `tools/merge_plain.py`; output and README in
`proVBFH-cs/production/reference/hxswg136/symtrim/`.

σ(ptj > 20): HH 2.07465 ± 0.00111, 11 2.09022 ± 0.00088, 22 2.11129 ±
0.00072 pb, i.e. 0.02–0.04% above the study's trimmed values. Over all
bins the symmetric trim is a median +0.05% above the study's values, but
more than 3 errors away in 81 of 946 bins. The plain merge stays a median
+0.34 plain errors above. The trimmed estimators are biased low by the
upper tail and the plain mean has unreliable errors, so for a comparison
the symmetric trim and the plain merge bracket the reference.
