# proVBFH-cs NNLO production on the MPCDF cluster (Oct 2026)

Instructions: `proVBFH-cs/production/CLUSTER-INSTRUCTIONS.md`. Set-ups
1506.02660 (`p1506`) and HXSWG 13.6 TeV (`hxswg136`), NNLO, `cs_scales 3`
(W1 = (1,1), W2 = (1/2,1/2), W3 = (2,2)).

Cluster: MPCDF t2 (login `mppui1.t2.mpcdf.mpg.de`), Slurm partition `alma`.
Raw per-job outputs stay on the cluster in `/ptmp/mpp/akarlber/cs-production/`
(not committed). Layout there:

- `bin/`: the binaries used, frozen copies (`GIT_COMMIT` gives the commit).
- `test/`, `regr/`: timing test and the regression check (section 6).
- `pilot/<setup>/{excl,incl}/job-<seed>/`: pilot jobs.
- `prod/<setup>/{excl,incl}/job-<seed>/`: production jobs.
- `logs/`: Slurm stdout per array task.

## Build (2 Oct 2026)

Branch `2026-09-cs-p2b` at 52686fd. gfortran 11.5 (system), hoppet 2.3.0,
LHAPDF 6.5.6, FastJet from `~/.local`; `proVBFH/Makefile.inc` from
`./configure`; `make ANALYSIS=p1506`, `make ANALYSIS=hxswg136` in
`proVBFH-cs`. PDF sets `NNPDF30_nnlo_as_0118` (261000) and `PDF4LHC21_40`
(93100; the HXSWG card's `lhans1 93100`) downloaded from
lhapdfsets.web.cern.ch into `~/.local/share/LHAPDF`.

## Scripts

- `proVBFH-cs/production/slurm/make_jobs.sh`: one directory per job from
  the card, unique `iseed`, appended to `jobs.list`; refuses to reuse a
  directory (so a seed is never run twice).
- `proVBFH-cs/production/slurm/run_array.sh`: array task i runs line
  i + OFFSET of `jobs.list`; skips finished jobs (`done` marker), so an
  array can be resubmitted to refill failed tasks.
- `proVBFH-cs/tools/pilot_estimate.py`: CPU needed per bin for a target
  error from the pilot's seed scatter.
- `proVBFH-cs/tools/plot_compare.py`: the comparison plots.

Seeds: exclusive 1000001..., inclusive 2000001... per set-up (pilot first,
production continues after the last pilot seed).

## Checks (section 6)

- Timing test (50k points, `cs_scales 3`): p1506 91 s, hxswg136 about 160 s
  CPU, 140 MB RSS.
- Regression: one exclusive job per set-up on a fixed grid (`readingrid 1`
  on the grid of a run without `cs_scales`), without `cs_scales` vs with
  `cs_scales 3`: W1 equals the run without variations in every bin of every
  histogram (largest absolute difference 4e-16 for p1506 and 1e-17 for
  hxswg136, i.e. rounding). A first attempt compared a warm-up run with a
  `readingrid 1` run, which see different points (statistical differences
  only): both runs of the comparison must read the same grid.

## Pilot

Cards as given (p1506: 200k x 2 warm-up, 1.6M x 3; hxswg136: 200k x 2,
600k x 3; inclusive 200k x 3, 2M x 3). 200 exclusive + 20 inclusive jobs
per set-up, submitted 2 Oct 17:35, all running within a minute.

- Inclusive: 145 s CPU per job. σ_tot (no cuts, inclusive part, 20 seeds):
  p1506 3.88774 ± 0.00035 pb, hxswg136 4.30914 ± 0.00039 pb.
- Exclusive: 184 of 200 (p1506) and 198 of 200 (hxswg136) completed, the
  rest NODE_FAIL (nodes kt12, ct18; not rerun). CPU per job (`cs_scales 3`):
  p1506 2.64 h (5.2M points, 1.83 ms per point), hxswg136 1.97 h (2.2M
  points, 3.2 ms per point); wall-time spread small (max/median 1.12). No
  NaN, `cs_spikes.dat` empty, RSS 65-135 MB.

### CPU needed per bin for the quoted (trimmed) errors of the old results

`pilot_estimate.py`, central scale (W1), seed-scatter errors, cluster
CPU-h (these include `cs_scales 3`). Targets: the quoted errors of
`reference/p1506/11.top` and `reference/hxswg136/nnlo-central.top`.

| | bins | median | 16% | 84% | worst |
|---|---|---|---|---|---|
| 1506: σ(2 jets) | 1 | 930 | | | |
| 1506: σ(≥ 3 jets) | 1 | 19,500 | | | |
| 1506: σ(≥ 4 jets) | 1 | 302,000 | | | |
| 1506: 2-jet distributions | 211 | 1,400 | 730 | 3,800 | (outliers) |
| 1506: 3-jet distributions | 102 | 64,000 | 36,000 | 273,000 | |
| 1506: 4-jet distributions | 33 | 492,000 | 193,000 | 2.0·10⁶ | |
| HXSWG: σ(incl cuts, p_T,j > 20) | 1 | 2,960 | | | |
| HXSWG: all bins | 946 | 14,900 | 2,300 | 271,000 | |

Compared with the thA371a estimates in the instructions (thserv18 CPU-h,
16 and 8 seeds): the 1506 numbers are in line once the faster cluster
cores are taken into account, σ(≥ 3 jets) is a factor 3 cheaper. The
HXSWG median is 4x higher than the 8-seed estimate (3,600): with 8 seeds
the scatter of heavy-tailed bins is underestimated. The "worst" bins
(up to 1e95) are bins where the old file's quoted error is zero to
rounding or garbage, not physics.

## NNLO production (submitted 2 Oct 22:2x)

AK (by message, away): "you are authorised to submit jobs for the runs as
long as the timings come in as expected". Sizing, about 7.5 h per job
(well under 24 h, few enough jobs for the 25k queue limit):

| | jobs | ncall2 x itmx2 | per job | total |
|---|---|---|---|---|
| p1506 exclusive | 11,000 | 4.8M x 3 | 7.5 h | 83k CPU-h |
| hxswg136 exclusive | 11,000 | 2.4M x 3 | 6.8 h | 75k CPU-h |
| inclusive, each set-up | 1,000 | 2M x 3 | 2.5 min | 40 CPU-h |

Expected: p1506 all 2-jet bins, σ(≥ 3 jets) and the median 3-jet bins
at or below the quoted errors; hxswg136 about two thirds of the bins.
Seeds: exclusive 1000201-1011200, inclusive 2000021-2001020. Array ids in
`prod-nnlo/array_ids.txt`. The two exclusive arrays have `nice=500` so
that the LO/NLO jobs submitted later start first.

## LO and NLO (AK, 2 Oct: "do the full LO and NLO runs as well with the scale variation ... the accuracy target should be much better ... should not visibly fluctuate in plots")

- LO = inclusive part with `qcd_order 1` (no exclusive part at LO).
- NLO = inclusive part with `qcd_order 2` + exclusive part with
  `cs_order 1` (the exclusive card's `qcd_order 2` kept; alpha_s and PDFs
  do not depend on it).
- Pilot (`pilot-lonlo`, 2 Oct 22:3x): 50 NLO exclusive jobs (5.2M points)
  and 20 + 20 inclusive jobs per set-up. NLO exclusive CPU per job:
  p1506 0.175 h (0.12 ms per point), hxswg136 0.385 h (0.27 ms per
  point). 19 tasks hung on et07/et08 (see below) and were cancelled.
- CPU needed for a relative error per bin (`pilot_estimate.py --relerr`),
  cluster CPU-h:

| | target | σ | median bin | 84% of bins |
|---|---|---|---|---|
| p1506 LO | 1e-3 | 0.01 | 0.8 | 13 |
| p1506 NLO 2-jet | 1e-3 | 2.0 | 71 | 2,150 |
| p1506 NLO 3-jet | 1e-3 | 31 (σ ≥ 3j) | 481 | 26,800 |
| hxswg136 LO | 1e-3 | 0.014 | 7.4 | 1,170 |
| hxswg136 NLO | 1e-3 | 19 | 2,640 | 299,000 |

- Production (`prod-lonlo`, submitted 3 Oct 00:45; sized to fit the 25k
  queue limit next to the NNLO production, long jobs because of the start
  rate): NLO exclusive p1506 450 jobs x 10.4 h (`ncall2` 28.6M), 4.7k
  CPU-h; hxswg136 1,800 x 10.5 h (`ncall2` 13M), 19k CPU-h; inclusive NLO
  150 jobs (`ncall2` 40M) and LO 150 jobs (`ncall2` 80M) per set-up.
  Seeds: NLO exclusive 3000051..., NLO inclusive 4000021..., LO
  5000021....
- **Correction (3 Oct 06:30).** The NLO exclusive job times above are
  wrong by a factor 3.6: I converted seconds to hours wrongly. 28.6M x 3
  points at 0.12 ms is 2.9 h, not 10.4 h (the first jobs finished in 3.0 h
  with the full 81.7M points). So the submitted NLO exclusive budget was
  1.3k (p1506) and 5.2k (hxswg136) CPU-h, not 4.7k and 19k. Fix: more jobs
  of the same size, to the intended budget: p1506 +1,200 (seeds
  3000501-3001700, 1,650 jobs in all), hxswg136 +4,800 (3001851-3006650,
  6,600 in all), fed into the queue as the 25k limit allows
  (`slurm/feed.sh`, state in `jobs.list.next`). With that the expected
  relative errors are the ones planned: p1506 NLO 2-jet bins about 3e-4
  (median) and 7e-4 (84%), 3-jet median 3e-4; hxswg136 NLO median 4e-4,
  84% 4e-3; LO below 1e-3 in 84% of the bins.
- NNLO sizing check: the production p1506 jobs take 4.4-5.4 CPU-h (median
  5.2), not 7.5 h, because the production nodes (et, gt) are faster per
  point than the pilot's (ct, kt). The statistics are set by the number of
  points: 11,000 x 14.8M = 170x the pilot's points, as planned.

## Interim result: 1506.02660 NNLO (3 Oct 08:00)

1,362 exclusive + 988 inclusive jobs, `combine_and_plot.sh` into
`/ptmp/.../interim/p1506` (scatter errors, no trimming). Page:
https://claude.ai/artifact/Wpkbr9jhcBCcoTWMRZnGrU

| [fb] | proVBFH-cs | old (11.top) |
|---|---|---|
| σ(≥ 2 jets) | 840.85 ± 0.34 | 843.83 ± 0.46 |
| σ(≥ 3 jets) | 127.20 ± 0.31 | 133.24 ± 0.06 |
| exactly 2 jets (difference) | 713.65 | 710.60 |
| σ(≥ 4 jets) | 17.07 ± 0.23 | 16.88 ± 0.01 |

**Correction to the earlier conclusion** (instructions section 1: "2-jet
σ: −0.56%, −3.4σ ..., of which the ≥ 3-jet region carries all"): the
≥ 3-jet deficit (−6.0 fb, −4.5%, as expected) is now larger than the
total deficit (−3.0 fb); the exactly-2-jet cross section is 0.43% above
the old code (about 4-5σ against the trimmed old errors). The 2-jet
distributions show χ²/n of 1-10 (p_T,H, y_H, H_T worst). The old errors
are trimmed (too small), so the pulls are overstated, but the shift is
not explained yet. Open, to be revisited with the full statistics.

"sig incl cuts" is not comparable: the new code fills it with the total
inclusive cross section (3.887 pb), the old file with that of its
`phspcuts`-restricted generation (0.923 pb).

## Cluster problems (night 2-3 Oct)

From about 19:30 on 2 Oct about 18 et* nodes stopped responding and about
57 nodes hung in "completing" (mostly another user's jobs). My pilot ran on
ct/kt nodes and was not involved. Effect on the production: job starts
dropped from about 90/min to a few per 10 min; about 180 tasks on
et07-et11 hung before their job script started (no Slurm output file). Those
will be resubmitted after the arrays end (`run_array.sh` skips finished
tasks).

## The `standard` partition is the old cluster (3 Oct 19:20-19:50)

AK: "given how badly the cluster is behaving you are free to try
experimenting with other queues". I moved the last 1,000 hxswg136 NNLO
exclusive tasks to `standard` (ct41-56) and reported them as running: that
was wrong. They showed RUNNING in squeue but all failed within ~20 s
(ExitCode 0:53, no output): `standard`/`long`/`extralong`/`special` belong
to the old cluster (CentOS 7), whose /u and /ptmp are not the new
cluster's (/u/akarlber and /ptmp/mpp/akarlber do not exist there; probed
with exit codes). Using it would need the old login node mppui4 (separate
home, CentOS 7 rebuild), which this session cannot reach. The 1,000 lines
were resubmitted on alma (array 48603568, nice 500); none had run.
`resubmit.sh` now counts arrays named `<name>-<suffix>` as well, and
keeps under the queue limit.

## NLO VBF H+3j cross-checks (section 9 of the instructions, from 3 Oct)

- VBFNLO 3.0 and the public POWHEG-BOX-V2 VBF_HJJJ (r4135): build and LO
  checks in progress, in `/ptmp/mpp/akarlber/h3j/{vbfnlo,powheg}`.
- **Correction (AK, 3 Oct ~19:45):** no fixed-scale runs. VBFNLO and the
  public POWHEG run in exactly the 1506.02660 set-up and analysis
  (μ0(p_T,H), the 1506 cuts, the p1506 observables), LO first and then
  NLO, so that their curves go into the existing 1506 plots; VBFNLO's
  histogramming is to be patched to fill our observables. No proVBFH runs
  are needed (the production is the reference). The 1,000 proVBFH-cs
  μ = m_H jobs (array 48603304, none started) were cancelled and the
  deferred p1506 NNLO tail (lines 10001-11000) resubmitted (array 48603439,
  nice 500).
- **VBFNLO 3.0 (3 Oct evening, built by a sub-agent).** CERN LCG mirror
  tarball, `vbfnlo-3.0-scale20.patch`, `./configure --enable-quad
  --with-LHAPDF=$HOME/.local` (default processes), install in
  `/ptmp/mpp/akarlber/h3j/vbfnlo/install`. Histogram patch
  `/ptmp/mpp/akarlber/h3j/vbfnlo/vbfnlo-3.0-p1506hist.patch`: new
  `utilities/p1506hist.F` books the 23 histograms of `p1506_analysis.f`
  (same names, binning, fill rules) from VBFNLO's own anti-kt 0.4 jets and
  the Higgs momentum; weights of one phase-space point (real + dipoles)
  summed before squaring; last iteration only; no smearing; output
  `p1506_lo.top` / `p1506_nlo.top` (pb) in the run directory. Smoke test:
  the patched σ(≥ 3 jets) equals VBFNLO's integral exactly. VBFNLO
  process 110 only has ≥ 3-jet events, so only the 3- and 4-jet
  histograms are comparable (the "2-jet" ones hold just the ≥ 3-jet part).
  The analysis has no p_T,H for ≥ 3-jet events (would need a new
  histogram in both). LO check and NLO timing test submitted (arrays
  48604035, 48604036; the sub-agent's own sbatch was refused by the
  permission system, so I submitted them after reviewing `run.sh`). The
  pending LO/NLO production arrays got `nice=100` so that the cross-check
  jobs (nice 0) start first; NNLO stays at 500.
- **VBFNLO checks (4 Oct 00:30).** LO check (process 110, ID 20, 2^20 x 4):
  130.596 ± 0.435 fb (expected 130.60 ± 0.44). The patched histogram holds
  the last iteration only (130.264 ± 0.497 fb, exactly VBFNLO's iteration
  4), as VBFNLO's own histograms do; VBFNLO's printed total is the
  weighted mean over iterations. NLO test (2^18 x 2): patched virtual
  151.985 fb = VBFNLO's; patched real −57.88 ± 33.90 fb = VBFNLO's last
  real iteration exactly (the printed −30.23 is again the mean over
  iterations). So the patch is right at NLO too. Comparison point at LO:
  the ≥ 3-jet histograms of our NLO production are tree-level H+3j:
  σ(≥ 3 jets) = 130.349 ± 0.019 fb (1,061 + 140 jobs), VBFNLO 130.60 ±
  0.44 fb. Timing: ~94 s for 2^18 x 2 → about 1 h per job at 2^23 x 5.
- **VBFNLO production (4 Oct 00:40),** dyn cards (μ0, ID 20):
  `/ptmp/mpp/akarlber/h3j/vbfnlo/prod/{lo,nlo}`, LO 50 seeds (SEED 1001-,
  2^24 x 5), NLO 200 seeds (SEED 2001-, 2^23 x 5), arrays 48609527/8,
  run through `run_array.sh` with the wrapper `vbfnlo-run`.
- **Public POWHEG VBF_HJJJ: LO is 6.3x too large (4 Oct 00:50).** The
  sub-agent stopped (usage limit) after building r4135 (only change: the
  imposed W/Z widths in `init_couplings.f`; checked against a fresh svn
  export) with the p1506 analysis, and running LO tests at ptcut 1 GeV,
  μ0 (`runningscales 4` = mthscale = the 1506 μ0, checked in
  `Born_phsp.f`): σ(≥ 3 jets, VBF cuts) = 0.82 pb (with and without Born
  suppression) against 0.130 pb from VBFNLO and proVBFH-cs. Not the
  histogramming (the analysis' "sig incl cuts" equals POWHEG's btilde
  total; the cuts act: M_jj > 600, Δy > 4.5, p_T,j3 > 25), not the event
  record (Higgs at position 3), not the couplings (printed: G_F, M_W,
  M_Z, sin²θ_W 0.2226, widths 2.141/2.4952). Testing the ptcut dependence
  at LO (must be none for ptcut < 25 GeV): runs `lo-ptcut20`, `lo-ptcut5`.
- **LO comparison, bin by bin (4 Oct 02:10).** VBFNLO LO production (50
  seeds) vs the ≥ 3-jet histograms of our NLO production (tree-level
  H+3j): σ(≥ 3 jets) 130.39 vs 130.35 fb; p_T,j3 and y_j3 agree bin by bin
  to 0.1-0.3% (ratio 0.997-1.003). VBFNLO's "2-jet" histograms hold only
  the ≥ 3-jet part, as expected (not comparable).
- **POWHEG VBF_HJJJ LO, ptcut dependence:** σ(≥ 3 jets) = 0.77 ± 0.05 pb
  at ptcut 20 GeV, 0.88 ± 0.09 pb at 5 GeV, 0.82 pb at 1 GeV: no ptcut
  dependence, still ~6x VBFNLO. Not a normalisation either: POWHEG/VBFNLO
  rises with p_T,j3 (5.9 at 30-35 GeV to 10.9 at 90-95 GeV) and towards
  forward y_j3 (4.4 central, 9.4 at 4-4.5). So the public code (r4135, as
  distributed apart from the W/Z widths) has a different tree-level H+3j
  with VBF cuts: extra or mis-weighted sub-processes. Per the
  instructions ("anything else, report to AK before going on") the POWHEG
  NLO production is on hold until AK has seen this. Next diagnostic:
  channel by channel (the code's `channel_type` switch in
  `init_processes.f`: NC/CC, and gluon-initiated vs quark-only), against
  VBFNLO / proVBFH-cs per channel.
- Tried an old-proVBFH Born-only H+3j run as a third LO reference
  (`/ptmp/mpp/akarlber/h3j/old-provbfh/lo-born`): the old code ignores
  `bornonly` ("unused variable"), needs the tags (`withfulltags 1`,
  otherwise it exits in `pwhg_analysis.f:51`), and ran its full exclusive
  part instead. Not pursued (would need code changes); VBFNLO and
  proVBFH-cs already agree with each other at LO.
- **VBFNLO NLO results (4 Oct ~05:00, 200 seeds, scatter errors).**
  Combined in `/ptmp/.../combined/p1506/vbfnlo/` (`vbfnlo-{lo,nlo}.top`).
  σ(≥ 3 jets) at O(α_s²): VBFNLO 125.25 ± 0.34 fb, proVBFH-cs NNLO
  127.16 ± 0.31 fb (1,386 jobs), old 133.24 fb. σ(≥ 4 jets): 16.97 ± 0.02,
  17.07 ± 0.23, 16.88. So the old code is 6.4% above VBFNLO (central y_j3,
  low p_T,j3); proVBFH-cs is 1.5% (4σ) above VBFNLO, at low p_T,j3, and
  agrees above ~50 GeV. thA371a: agreement at μ = m_H, 1.4 fb apart at
  μ0, so most likely the dynamic-scale definition in the real emission
  (VBFNLO ID 20: p_T,H of the actual event). Not settled; for AK.
  χ²/n VBFNLO vs proVBFH-cs: p_T,j3 24/15, y_j3 24/18, y*_j3 45/24,
  min-rap 61/46, p_T,j4 18/15, y_j4 64/18 (proVBFH-cs fluctuations in two
  central bins). Plots (ratio to proVBFH-cs): `plots-h3j-nlo`,
  `plots-h3j-lo`; on the page (version 4).
- **POWHEG LO excess found (4 Oct 11:00): a bug in the public code's
  running scales.** `set_fac_ren_scales` (`Born_phsp.f`) computes `muref`
  in every running-scale branch but sets `muf`, `mur` only in the fixed
  branch, so for `runningscales` > 0 both stay 0 and are floored at
  √2 GeV (instrumented run: mur = 1.414 GeV, α_s = 0.356 at every point).
  The default card (`runningscales 0`, m_H/2) is not affected. Steps that
  excluded the other causes: public `setborn` = proVBFH-cs m1+m2 to all
  digits at fixed points (`/ptmp/.../h3j/pointtest`, `pub.f`, `cs/csp.f90`);
  independent integration of the public `setborn` with its flavour list,
  own phase space, LHAPDF and parton-level cuts (`mc.f`): 128 ± 8 fb;
  pristine `init_couplings.f`, `smartsig 0` and `fullphsp` all unchanged
  (0.77 pb at ptcut 20). Fix (two lines, `muf=muref`, `mur=muref` at the end
  of the running branch) in `VBF_HJJJ/Born_phsp.f` (original kept as
  `.orig`): LO σ(≥ 3 jets) = 121.9 ± 6.5 fb (one seed, ptcut 1). Added as
  problem 4 to the bug report https://claude.ai/artifact/LLJKDy2845Gd7FM8tDwKMP.
- **POWHEG production (4 Oct 11:30, fixed build):** `runs/nlo-pt1`,
  `runs/nlo-pt01` (ptcut 1 and 0.1 GeV): stage 1 20 seeds x 200k, stage 2
  300 seeds x 2M (~2.6 h each), arrays 48623361-4; `runs/lo-prod` (ptcut 1,
  Born suppression): 10 + 50 seeds, arrays 48623365-6.
- **More VBFNLO NLO statistics:** 600 more seeds (2201-2800), array
  48618512.
- **POWHEG LO with the fix, 50 seeds (4 Oct 15:40):** σ(≥ 3 jets) = 124.83 ±
  1.09 fb vs VBFNLO 130.39 ± 0.01: −4.3%, roughly flat in p_T,j3, y_j3.
  Set-up difference, not a bug: the public code has `max_flav = 4` (no b
  quarks), while VBFNLO (VBFHASB), proVBFH-cs and the old proVBFH include
  b quarks in NC. Test with `max_flav = 5` + the scale fix
  (`VBF_HJJJ-nf5fix`, `runs/lo-nf5fix`): 127.7 ± 2.6 fb after 14 of 30 seeds.
  The like-for-like NLO comparison needs the 5-flavour build too:
  `runs/nlo-nf5-pt1` (300 seeds, arrays 48632637/8). The ptcut dependence
  is measured with the 4-flavour pair `nlo-pt1`/`nlo-pt01`.
- POWHEG NLO preview (17:40, 4 flavours): σ(≥ 3 jets) 129.2 ± 1.8 fb (ptcut
  1, 270 seeds), 130.1 ± 2.5 (ptcut 0.1, 71 seeds); VBFNLO 125.73 ± 0.13
  (718 seeds). σ(≥ 4 jets) 15.93 ± 0.11, 16.10 ± 0.28 vs 16.97 ± 0.01 (no b).
- et32 died like et36 (4 Oct ~15:30): all my tasks there cancelled and
  resubmitted with `--exclude=et32,et36`.
- **Negative PDFs (4 Oct ~20:00, AK).** All codes use NNPDF30_nnlo_as_0118
  (261000). POWHEG's `pdfcalls.f` sets negative PDF values to zero
  (counter "negative pdf values": ~9M per run); VBFNLO and proVBFH-cs keep
  them. Variant `VBF_HJJJ-nf5fix-noclip` (5 flavours + scale fix + a local
  `pdfcalls.f` without the zeroing; the shared file untouched). LO test
  `runs/lo-nf5-noclip` (10 + 30 seeds). AK: restart the 5-flavour NLO
  without clipping: the clipping run (48632637/8) cancelled; new run
  `runs/nlo-nf5nc-pt1`, 20 + 1,000 seeds (arrays 48635075/6).
- **LO gap, further checks (4 Oct 20:30).** Negative PDFs kept
  (`lo-nf5-noclip`, same seeds as `lo-nf5fix`): 126.26 ± 1.40 vs 126.62 ±
  1.40 fb, so clipping is −0.3% only. α_s: POWHEG's st_alpha = LHAPDF's
  α_s at the same μ_R (0.12490 at 63.3 GeV). Remaining candidates: the
  Higgs Breit-Wigner of the public phase space (the old agent's narrow-
  width test: +0.9 ± 0.7%) and statistics (2.9σ with 30 heavy-tailed
  seeds). Test: narrow width + 5 flavours + no clipping, 100 seeds
  (`VBF_HJJJ-nf5nc-nw`, `runs/lo-nf5nc-nw`).
- **Narrow width closes the LO gap (5 Oct 03:00).** POWHEG LO with the
  scale fix, 5 flavours, negative PDFs kept and a narrow-width Higgs
  (`parameter (BW=.false.)` in `Born_phsp.f`; VBFNLO and proVBFH-cs use
  narrow width): 131.29 ± 0.90 fb (100 seeds) vs VBFNLO 130.39 ± 0.01
  (+0.7%, 1.0σ). The Breit-Wigner runs give ~126.5 fb (−3.8% ± 1.4%); the
  ±30 Γ window of the public code accounts for only ~1.1% of that, the rest
  is not understood (statistics, or the off-shell Higgs mass against the
  fixed `kn_masses(3)` of the phase-space maps). Not pursued: narrow width
  is the like-for-like set-up. So the like-for-like POWHEG NLO run is
  `runs/nlo-nf5nc-nw-pt1` (20 + 700 seeds, arrays 48641982/3); the
  pending tasks of the Breit-Wigner run `nlo-nf5nc-pt1` were cancelled
  (its ~300 finished + 142 running seeds kept for comparison).

## Smaller queue footprint (5 Oct 08:30)

AK: "There is almost 40k jobs queueing so this is not sustainable." About
22k of them were mine (pending NNLO exclusive tasks). Cancelled the pending
tasks of the four NNLO exclusive arrays (48568055, 48603439, 48568057,
48603568: 20,308 tasks, none started); my queue went from 22,780 to
2,472. From now on `slurm/feed_missing.sh` (run hourly) keeps my total
queue below 2,500 and submits only never-started lines (no `done`, no
`time.log`), recording them in `jobs.list.fed`; hxswg136 first, then
p1506, nice 500. `resubmit.sh` is capped the same way (QCAP, default
2,500). Lines that started but did not finish are left for `resubmit.sh`
at the end.
- **Cluster recovered (5 Oct ~10:00).** alma nodes stuck in completing:
  74 at 08:00, 4-8 at 10:00; ~1,500 of my jobs running. The hourly cap
  went to 20,000 automatically, but the feeder's first submission failed
  ("Pathname ... too long": sbatch rejects an --array list of ~10,900
  single indices); fixed by compressing the indices into ranges
  (`slurm/ranges.py`). Fed 10,887 hxswg136 and 7,586 p1506 NNLO lines
  (queue 20,000). ~2,000 p1506 lines that started earlier but never
  finished (dead nodes, cancellations) are left for `resubmit.sh` at the end.
- **H+3j result, like-for-like (5 Oct 13:00).** POWHEG with the scale fix,
  5 flavours, negative PDFs kept and narrow width: LO 131.29 ± 0.90 fb
  (VBFNLO 130.39), NLO σ(≥ 3 jets) = 133.4 ± 3.0 fb (660 seeds; 40 lost to
  a NODE_FAIL on et13), +6.1% above VBFNLO (125.69 ± 0.12), like the old
  proVBFH (133.24, +6.0%); proVBFH-cs 127.16 ± 0.31 (+1.2%). σ(≥ 4 jets):
  16.84 ± 0.27 vs 16.97 ± 0.01. Breit-Wigner variant (646 seeds): 131.4 ±
  2.9. ptcut 1 / 0.1 GeV (4 flavours, 300 seeds each): 128.3 ± 1.7 / 130.0
  ± 1.3, no significant dependence. Page version 6.
