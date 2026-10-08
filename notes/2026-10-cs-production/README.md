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

## POWHEG with the NLO bugs fixed (5 Oct 17:00)

AK: "set up the fixed POWHEG variants, we need to get to the bottom of
this. Also to make sure that we can actually fix the bugs." On top of the
like-for-like build (`VBF_HJJJ-nf5nc-nw`): `VBF_HJJJ-lfl-fix1` (problem 1,
`compreal_hqqqq.f:529`, `kl = k+4*(2-ftype(7))`), `-fix3` (problem 3,
`real.f:1134` only, `ftype(2)` from `bflav(5)`; the same line text in the
other gg branches is correct and unchanged), `-fix13` (both). Patches
`h3j-crosscheck/powheg/patches/06-*`, `07-*`.

Limit tests of the bug report (its appendix programs, extracted and built
against both builds, `/ptmp/mpp/akarlber/h3j/limittest`): unfixed c =
1.281746 / 0.780187 (problem 1, Q = u / d) and 0.780185 (problem 3), the
report's values; fixed c = 1.000000 for all, k_T → 0. NLO runs (ptcut 1,
20 + 700 seeds each): `runs/nlo-lfl-fix{1,3,13}-pt1`, arrays
48674548-48674553. Problem 2 (missing initial-state FKS region) needs new
flavour-list entries and a split of the four-quark matrix element by graph
class; not done yet.
- **Correction (5 Oct 17:30).** I first estimated problem 2 at ~0.1% of
  σ(≥ 3 jets) per e-fold, so ~2% in total, "too small for the 6%". That
  kept only the logarithmic part. The stage-3 emulation (notes of 2 Oct,
  cs_estimate 6/7) found for the old code: missing region above 1 GeV plus
  the collinear remnant that `btildecoll` adds anyway, 3.15 ± 0.15e-3 pb;
  realised below 1 GeV 1.07 ± 0.30e-3; n_f mismatch 0.87e-3; together
  5.1 of the 6.6e-3 pb excess. So problem 2 is the main driver; fixes 1
  and 3 are expected to move POWHEG little. Pushing ptcut tests only the
  log part (~0.17%/e-fold) and is not worth the CPU.
- `doublefsr` (`find_regions.f:633`) only adds, for final-state regions
  already found, the copy with emitter and radiated parton swapped; it
  cannot create the missing initial-state region.
- Problem 2 fix: delegated to a sub-agent (5 Oct 17:30) working only in a
  new copy `VBF_HJJJ-lfl-fix123`, no job submissions; deliverables and
  validation in `/ptmp/mpp/akarlber/h3j/fix2-tests/`.

## Problem 2 fixed (sub-agent, 5 Oct 17:30-17:50)

Build `VBF_HJJJ-lfl-fix123` (fixes 1 and 3 plus this), patch
`h3j-crosscheck/powheg/patches/08-fix2-nc-isr-region.patch`; test programs
in `tools/fix2/`, outputs in `/ptmp/mpp/akarlber/h3j/fix2-tests/`.
- **Why tags alone are not enough:** entries with the same flavours but
  different tags are separate processes in POWHEG (own matrix element, FKS
  S functions only over their own regions), so a new entry returning the
  full matrix element would double count; and in VBFNLO's `qqh4q` the
  Z-on-line and Z-on-pair graphs interfere coherently.
- **Fix:** 1,000 new NC entries with the pair tag on incoming 1 (from the
  NC gq Borns, tags 5,2,0,1,2,1,5) and on incoming 2 (NC qg Borns, tags
  1,5,0,1,2,2,5), analogous to the CC ones; reals 1,277 -> 2,277 (<
  maxprocreal 2,392). The matrix element is split across the entries with a
  smooth partition of unity built from POWHEG-type distances (d_F of the
  pair, d_I1, d_I2 of the outgoing quark of each line to its beam): the
  outgoing-pair entry gets nc_up d_I1/(d_I1+d_F) + nc_lo d_I2/(d_I2+d_F),
  the new entries nc_up d_F/(d_I1+d_F) and nc_lo d_F/(d_I2+d_F); the
  existing `pair71`/`pair72` routing of `compreal_hqqqq` gives exactly
  nc_up/nc_lo there (checked). No VBFNLO change, no interference dropped;
  btildecoll already adds the qg remnant to the gluon Borns.
- **Validation:** (a) regions 2,432 -> 3,432: each new entry has one ISR
  region (emitters 1 / 2), all old entries and CC unchanged; the bug
  report's region program finds 660 = 160 CC + 500 NC entries with an ISR
  region. (b) Collinear limit (outgoing s || beam 1): the 1/k_T^2 is now in
  the new entry (c_in1 -> 0.99999966), the outgoing-pair entry tends to a
  constant; c_sum equals fix13 digit for digit; summed over entries, every
  NC flavour structure equals the unfixed matrix element to 8.9e-16 at
  random points (checked by me in `compare_sum.log`), all other entries
  bit-identical. (c) POWHEG's own soft/collinear checks: all new regions
  converge, with the same ~0.5% offset as the untouched CC pair-on-beam
  regions (1.0046 / 1.0068); not caused by the fix, not investigated.
- **Side finding:** all like-for-like builds had `st_nlight = 4` with
  `max_flav = 5`: γ_g (sigsoftvirt) misses the g -> b b̄ integrated
  counterterm and the qg remnant misses the b-initiated regions. Patch
  09 sets `st_nlight = 5` (old proVBFH, VBFNLO and proVBFH-cs use 5).
  The like-for-like NLO numbers so far (133.4 ± 3.0 fb) and the fix1/3/13
  runs use 4.
- **Runs (5 Oct 18:00), ptcut 1, 20 + 700 seeds each:**
  `nlo-lfl-fix123-nl5-pt1` (all three fixes, st_nlight 5; arrays
  48674723/4) and the control `nlo-lfl-nl5-pt1` (like-for-like, bugs left
  in, st_nlight 5; 48674725/6).

### 5-6 Oct: node failures and a slow node
- The et nodes keep failing: they go "Not responding" exactly on slurmctld's 1000 s ping cycle, come back on their own (ReturnToService=1), and fail again 1.5-3 h later, killing every running job each time. About 4,000 jobs were lost on 5 Oct; `resubmit_failed.py` reruns them.
- Evidence sent to AK for the admins:
  - The nodes do not reboot, and slurmstepd stays alive.
  - They fail in fixed groups.
  - IPv6 ping fails to 29 of 32 et nodes and works for all other alma nodes.
- ct30 runs jobs about 1.8x slower than other ct nodes, so its 12 h jobs time out. It was added to `bad_nodes`, and the pending arrays were updated with `scontrol update ExcNodeList`.
- Mean wall time of the NNLO exclusive jobs by node family: et 5.8 h, ct 6.3 h, kt 6.9 h, gt 7.5 h.

### 6 Oct 03:50: POWHEG H+3j runs with the fixes (like-for-like, ptcut 1 GeV)
σ(≥3 jets) at O(αs²), μ0, 1506.02660 cuts (`sigtot.py`; errors are the seed scatter):

| run | seeds | σ(≥3 j) [fb] | σ(≥4 j) [fb] |
|---|---|---|---|
| like-for-like, no fixes (`nlo-nf5nc-nw-pt1`) | 660 | 133.4 ± 3.0 | 16.84 ± 0.27 |
| + fix 1 | 698 | 134.0 ± 2.9 | 16.84 ± 0.25 |
| + fix 3 | 668 | 133.3 ± 2.9 | 16.83 ± 0.26 |
| + fixes 1, 3 | 658 | 134.5 ± 3.0 | 16.81 ± 0.27 |
| + fixes 1, 2, 3, `st_nlight 5` | 700 | **127.4 ± 2.9** | 16.79 ± 0.25 |
| proVBFH-cs (NNLO run) | | 127.16 ± 0.31 | 17.06 ± 0.23 |
| VBFNLO 3.0 | 800 | 125.69 ± 0.12 | 16.98 ± 0.01 |

- With all fixes, POWHEG agrees with proVBFH-cs. The χ² against proVBFH-cs is:

  | observable | χ² with all fixes | χ² without fixes |
  |---|---|---|
  | ptj3 | 18.3/15 | 24.9/15 |
  | yj3 | 10.9/18 | 25.7/18 |
  | y*j3 | 12.7/24 | 39.2/24 |
  | ptj4 | 40.2/15 | 30.2/15 |
  | yj4 | 21.8/18 | 27.9/18 |

  So POWHEG with all fixes now reproduces the new result, not the old proVBFH one (133.2 fb).
- Fixes 1 and 3 have no visible effect at this precision (each run's error is ±2.2%). The shift comes from fix 2 and/or `st_nlight 5`.
- To separate the two, the control run `nlo-lfl-nl5-pt1` (`st_nlight 5`, no fixes) is still running. Its stage-1 seed 3 was lost to a NODE_FAIL, which left stage 2 at DependencyNeverSatisfied. Seed 3 was rerun (48685521) and stage 2 (48674726) now depends on it.
- `h3j_update.sh` now plots the run with all fixes instead of the 0.1 GeV ptcut run, because the plotter has only four alternative styles.

### 6 Oct 07:10: control run, and which change causes the shift
- The control run `nlo-lfl-nl5-pt1` (`st_nlight 5`, no fixes; 696 of 700 seeds) gives σ(≥3 j) = 133.2 ± 2.9 fb and σ(≥4 j) = 16.75 ± 0.25 fb.
- All variants use the same `pwgseeds.dat`, so seed i of two runs shares its random numbers. The per-seed values are 92-100% correlated, so differences taken seed by seed are much more precise than the separate errors suggest:

  | change | seeds | Δσ(≥3 j) [fb] | correlation |
  |---|---|---|---|
  | fix 2 (with fixes 1 and 3, `st_nlight 5` in both runs) | 696 | **−5.67 ± 1.14** | 0.92 |
  | `st_nlight` 4 → 5 (no fixes) | 656 | −0.60 ± 0.64 | 0.98 |
  | fixes 1 and 3 (`st_nlight 4`) | 618 | +0.17 ± 0.28 | 1.00 |

- Problem 2, the missing initial-state FKS region for the NC four-quark graphs, accounts for the whole difference between the new and old proVBFH 3-jet rates: 127.16 − 133.24 = −6.1 fb.
- Fixes 1 and 3 and the `st_nlight` choice change σ(≥3 j) by less than 0.5%.
- Paired by seed: fix 1 alone +0.19 ± 0.26 fb on σ(≥3 j) and +0.013 ± 0.001 fb on σ(≥4 j). Fix 3 alone changes nothing: the .top files are bit-identical to the run without fixes, seed by seed.
- Why fix 3 changes nothing (`tools/fix3/gg_all.f`, `gg_alr.f`, linked with `tools/fix2/link.sh`):
  - The 33 gg real entries reach the faulty NC branch of `compreal_hjjj` in the entry order of `init_processes.f`; 25 do.
  - POWHEG calls `setreal` with the flavour order of each FKS region (`flst_alr`) instead.
  - None of the 132 gg regions reaches that branch with pairs of different types.
  - So problem 3 is dormant in POWHEG runs. It is real when the routine is called in entry order, as in our limit test and in proVBFH's copy.
- Bug report updated to version 5:
  - New section "Effect on NLO predictions".
  - Problem 3 marked as dormant in POWHEG runs.
  - Effect sizes for problems 1 and 2 added.

### 6 Oct 11:00: audit of VBFNLO's H+3j NLO (process 110), by a subagent
- No omission, approximation or technical cut that could give the −1.2%:
  - The virtual's gauge-check fallback (qqhqq.F:669-710) changes σ_virt by 2e-7 (tested via LD_PRELOAD).
  - The pair-mass cut 2p_i·p_j < (0.1 GeV)² in the phase space, and the z > 0.999995 cutoff in the K/P terms, are negligible.
  - VBFNLO uses full CS dipoles, with no α_dip parameter and no ycut.
  - The I, K and P operators were checked analytically.
  - All real subprocesses are present, including the NC pair graphs with their initial-state dipoles (so VBFNLO does not have POWHEG's problem 2), gg → H 4q, identical-flavour interference, and b quarks in NC.
  - α_s comes from LHAPDF at μ_R everywhere, with the same PDF set. The μ0 patch is per configuration, and the dipole maps keep p_H.
- Combination:
  - VBFNLO's own printed total (inverse-variance weighting over 5 iterations) is 127.02 fb. It is biased upward by the heavy negative tail of the real weights.
  - Our 125.69 is the plain mean of the last iterations, which is the right estimator.
  - proVBFH-cs is not affected: its histograms weight every point equally (`pwhgaccumup` normalisation over all 3 iterations), and `combine_parts.py` takes a plain mean over jobs.
  - POWHEG: `sigtot.py` takes a plain mean over seeds of the `pwhgaccumup`-normalised .top files.
- The earlier fixed-scale comparison (μ = m_H) agreed: 125.82 ± 0.52 against 125.50 ± 0.67 fb. At μ0 the gap is 1.47 ± 0.33 fb. The two gaps differ by about 2σ only.
- Agent files: `/ptmp/mpp/akarlber/h3j/vbfnlo/audit-2026-10-06/` (jobs.txt, real_iters.txt, preload tests).

### 6 Oct 11:30: fixed-scale (μ = m_H) H+3j comparison, redone with plain means
- **Why:** the stage-3 fixed-scale agreement (VBFNLO 125.82 ± 0.52, proVBFH-cs 125.50 ± 0.67 fb) came from 20 VBFNLO jobs on thA371a. Those runs are not reachable from here. Their μ0 value (126.58) is close to VBFNLO's biased printed total (127.02) and far from the plain mean (125.69), so they probably used the printed, inverse-variance-weighted totals. The agreement at μ = m_H may therefore be an artefact.
- **Updated μ0 gap:** the current production gives σ(≥3 j) = 126.67 ± 0.16 fb (about 4,100 jobs, `combined/p1506/nnlo-W1.top`, 07:50). The 127.16 ± 0.31 in `h3j_update.sh` was from an older combination. The μ0 gap to VBFNLO is now 0.98 ± 0.20 fb (0.8%).
- **New runs** (same seeds as the μ0 runs, so each code's scale effect can be compared seed by seed):
  - **VBFNLO:** 800 jobs in `/ptmp/mpp/akarlber/h3j/vbfnlo/prod/nlo-fixmh` (seeds 2001-2800, ID_MUF = ID_MUR = 0, MUF_USER = MUR_USER = 125). Array 48695214, nice 0, 4 h, 2000 MB. Failed jobs are resubmitted by hand.
  - **proVBFH-cs:** 2,000 NNLO exclusive jobs in `/ptmp/mpp/akarlber/cs-production/fixmh/p1506/excl` (seeds 1000201-1002200, the production card with `runningscales 0`, `cs_scales 1`, `ncall2 4800000`). Array 48695215, nice 50, so they run ahead of the production. `hourly.sh` resubmits failed jobs.
- **Target errors:** about ±0.12 fb for VBFNLO and about ±0.2 fb for proVBFH-cs. This delays the NNLO production by roughly 2-3 h.
- 13:55, VBFNLO at μ = m_H (778 of 800 jobs, plain mean of the last iteration): σ(≥3 j) = **125.24 ± 0.12 fb**.
  - The stage-3 value of 125.82 ± 0.52 was 0.58 fb higher, consistent with that number being a printed, inverse-variance-weighted total.
  - Seed by seed, μ0 minus m_H = +0.45 ± 0.18 fb. The pairing does not help here: the correlation is −0.01, because VBFNLO adapts its grids per run.
  - The 2,000 proVBFH-cs jobs at μ = m_H are all running and should finish tonight.
- 19:10, preliminary proVBFH-cs at μ = m_H (812 of 2,000 jobs, plain mean over jobs of pwg-EXCL.top): σ(≥3 j) = **126.21 ± 0.34 fb**.
  - Against VBFNLO at μ = m_H (125.22 ± 0.12): +0.99 ± 0.36 fb.
  - At μ0, from all 10,598 finished production jobs (plain mean of pwg-EXCL-W1.top): 126.88 ± 0.17 fb, i.e. +1.19 ± 0.21 fb against VBFNLO.
  - So the gap is the same at fixed and dynamic scale. The dynamic-scale treatment is ruled out as its cause.
  - Scale shift μ0 − m_H: proVBFH-cs 0.72 ± 0.28 fb (seed by seed, 793 seeds), VBFNLO 0.47 ± 0.17 fb. These are consistent.
  - σ(≥4 j): μ0 17.13 ± 0.16 fb, m_H 14.88 ± 0.20 fb.

### 6 Oct 23:00: channel split, first look at how to do it
- Neither code has a channel switch: proVBFH-cs has no input option, and VBFNLO's process 110 has none either. Both need a small diagnostic patch, applied to copies and not to the production builds.
- **Initial state** (qq, qg, gq, gg): per-beam flags on the PDFs, set by environment variables.
  - The flags must act on every PDF evaluation: the Born, the reals, the dipoles, and the K/P terms at x/z.
  - proVBFH-cs: `hoppetEval` fills `fB`, `fE` and friends per line.
  - VBFNLO: `pdfproton` is called per beam in m2s_qqh3j.F:166/169 and in the K/P and real drivers. Patch at the call sites, because `pdfproton` does not know which beam it is evaluating.
- **NC/CC:**
  - proVBFH-cs: `compatible(w1,w2)` in cs_exclusive.f90 decides class pairs (NC: w1 = w2 = 0; CC: w1 = −w2 ≠ 0). A switch there selects NC or CC.
  - VBFNLO: the subprocess loop in m2s_qqh3j.F (the `wbf_h3j` calls at lines 316-444), m2s_qqh4j.F and qqh4q_sub.F have to be checked for where ZZ and WW fusion are separated.
- **Plan:** run at μ = m_H. Use 8 runs per code (NC/CC × qq/qg/gq/gg), or 4 with qg+gq together. Size each from the per-channel error of short pilots.

### 7 Oct 01:45: fixed-scale comparison complete; p1506 NNLO production complete
All values are plain means over jobs (`fixcmp.py` in the scratchpad: per-job sig(all VBF cuts N jets) values, seed-scatter errors).

| σ [fb] | proVBFH-cs, μ = m_H (1,998 jobs) | VBFNLO, μ = m_H (800) | proVBFH-cs, μ0 (11,000) | VBFNLO, μ0 (800) |
|---|---|---|---|---|
| ≥3 jets | 126.19 ± 0.18 | 125.22 ± 0.12 | 126.88 ± 0.17 | 125.69 ± 0.12 |
| ≥4 jets | 14.80 ± 0.10 | 14.815 ± 0.007 | 17.12 ± 0.16 | 16.976 ± 0.008 |

- **≥3 jets, proVBFH-cs − VBFNLO:** +0.97 ± 0.21 fb at μ = m_H (+0.8%, 4.5σ) and +1.19 ± 0.20 fb at μ0 (+0.9%, 5.8σ). The offset is the same at both scales, so the dynamic scale is not its cause.
- **Scale shift μ0 − m_H:** proVBFH-cs +0.87 ± 0.14 fb (paired, 1,998 seeds), VBFNLO +0.47 ± 0.17 fb. They differ by 1.8σ, which is not significant.
- **≥4 jets:** both codes agree at both scales. The offset therefore sits in the exactly-3-jet part of the H+3j NLO: virtual, real minus dipoles, and K/P.
- The stage-3 statement "proVBFH-cs agrees with VBFNLO at μ = m_H" is withdrawn. It rested on VBFNLO's printed, inverse-variance-weighted total (125.82) and on a proVBFH-cs value with a ±0.67 fb error.
- The p1506 NNLO production is complete (11,000 of 11,000 exclusive jobs).

### 7 Oct 07:45
- hxswg136 NNLO: 10,858 of 11,000 done. The last 128 failed lines could not be resubmitted under the 2,500 cap: only 52 alma nodes are healthy, and the disorder jobs of the same user fill the queue. They were resubmitted by hand (48723694). `hourly.sh` now uses cap 30,000, since only small resubmissions remain.

### 7 Oct 11:30: hxswg136 NNLO production closed at 10,872 of 11,000 jobs
- The cluster broke down at about 10:20: about 18k jobs of another user, our 128 last hxswg136 resubmissions (after 2h47 of running), and about 7k disorder jobs failed.
- AK: skip the last 128 jobs, they have no impact. The resubmission (48759353) was cancelled by id. The hourly loop is stopped.
- Final sets: p1506 NNLO 11,000 of 11,000; hxswg136 NNLO 10,872 of 11,000; all LO/NLO sets complete (hxswg136 NLO excl 6,540 of 6,600).

### 7 Oct 12:30: spikes in the p1506 NNLO plots (AK)
`tools/spikescan.py` takes the 11,000 per-job W1 histograms. For each bin it finds the job whose removal moves the plain mean most, measured in seed-scatter errors.
- **Single events dominate many tail bins.** In the worst bins one job sits 100-105 standard deviations from the mean (√11,000 ≈ 105). Removing it moves the mean by 1σ, i.e. that single job carries both the bin's value and its error.
- **The same few jobs recur.** job-1005395 (22 bins with a shift above 0.5σ), 1008233 (15), 1001652 (12), 1004554 (10), 1009943 (7), 1010580 (6) and 1008322 (6).
- **Example, job-1005395.** One event of about +3.2 pb/unit appears at the same time in yj1, yj2, yj3, yj4, yH, Δy_jj, R_jj and Δφ_jj, and in σ(≥4 jets): 1.60 pb in that job against a mean of 0.0171. So it is a single four-parton event passing the 4-jet cuts. It alone contributes 0.15 fb (1σ) to σ(≥4 jets).
- **Below the logging threshold.** Only one point in all 11,000 jobs exceeded the logging threshold spike_min = 10 pb (job-1003743). Its real and dipole cancel (+10.40, −10.41). The events behind the visible spikes are therefore below that threshold.
- **Next:** rerun the four worst jobs with the same seed (deterministic) and `cs_spikemin 0.2`. This logs the events with their random numbers, for `cs_replay` (kinematics, which piece, which flavour). Array 48759513, `/ptmp/.../cs-production/spikes/`, about 5.5 h.

## 8 Oct: channel split of NLO H+3j, proVBFH-cs vs VBFNLO (μ = m_H)
Files: `proVBFH-cs/production/h3j-crosscheck/chsplit/`; runs in `/ptmp/mpp/akarlber/chsplit/`.

**Definitions** (same in both codes).
- Initial state: per-beam flags on the parton type of *every* PDF evaluation (Born, virtual and I at the Born x; K/P at x/z; reals and dipoles at the real's x): quark/antiquark or gluon from beam 1, beam 2 → qq, qg, gq, gg. A real and its dipoles share the PDF, so each channel is IR-finite and qq + qg + gq + gg = all exactly. The K/P term P_gq ⊗ f_q of a gluon-initiated Born counts as quark-initiated, like the qq → 4q H real it subtracts.
- Boson: NC (Z) / CC (W). Both codes keep NC and CC incoherent (VBF approximation, no NC–CC interference, also for identical flavours), so this is a clean partition. proVBFH-cs: `compatible()` (class pairs, all Born/virtual/(1,0)/(1,1) terms) and the real groups of `cs_nlo2` (NC iff the Born-level line keeps its flavour; checked to be uniform within each group). VBFNLO: the NC/CC matrix-element arrays scaled on output of every ME routine (`qqHqqj_c_virt`, `qqHqqj_spcor`, `qqh2q2g`, `qqh4q`), which feed Born, virtual, reals and every dipole.
- Switches: environment `CHAN_BOSON=NC|CC|ALL`, `CHAN_INIT=ALL|<b1><b2>` with q, g, a(ny). proVBFH-cs also `CHAN_MULTI=1`: nine weights from the same points (W1 all, W2–5 NC qq/qg/gq/gg, W6–9 CC), through the scale-variation weight machinery with the ME cache.
- Patches: proVBFH-cs `src/cs_chan.f90` + small changes in `cs_nlo2`, `cs_exclusive`, `cs_main`, Makefile (commit of this section); VBFNLO `vbfnlo-3.0-chsplit.patch` (on a copy, `/ptmp/.../chsplit/vbfnlo/{src,install}`).
- Side fix: `p1506_analysis.f` / `hxswg136_analysis.f` had `dsig(7)`; more than 7 weights overflowed it (W8, W9 garbage in the first test). Now `dsig(10)` (= maxmulti). No effect on the production (≤ 7 weights).
- VBFNLO K/P (checked while patching): m2s_qqh3j.F:279 uses `Cx(1)` for antiquarks on beam 2 where the other three cases use `Cx(2)`; harmless, `fincollinear` sets C(2) = C(1).

**Validation (login node, short runs; `topcmp.py`, `vbfnlo_tests.sh`).**
- proVBFH-cs (40k-point run of seed 1000201, `/ptmp/.../chsplit/test-cs`): patched binary with CHAN_MULTI: same VEGAS grid as the production binary (bit-identical `grids-excl.dat`), W1 = production result to 1e-15 in every bin (summation order only); W2 + … + W9 = W1 in all 426 bins to the 8-digit output precision. Env mode `CHAN_BOSON=CC CHAN_INIT=gq` on a fixed grid = W8 of the CHAN_MULTI run on the same grid, bit-identical in every bin. CPU of CHAN_MULTI ≈ the plain run (482 vs 417 s; ME cache).
- VBFNLO (2^17 points, 1 iteration LO and NLO, seed 2001, so all runs see the same points): patched copy with no channel = original install, `p1506_nlo.top` bit-identical; NC + CC = all and qq + qg + gq + gg = all in all 426 bins to the 8-digit output precision.

**Pilots (8 Oct 13:55).** proVBFH-cs CHAN_MULTI, the fixmh card, seeds 1000201-1000300 (the fixmh seeds: W1 must reproduce those jobs), array 48798835; VBFNLO 8 channels (NC/CC × qq/qg/gq/gg) × seeds 2001-2020, full fixmh size (2^23 × 5), array 48798836. Runner: the production `run_array.sh` (unchanged), via `submit.sh`.

### 8 Oct: channel-split pilot results (σ ≥3 jets, all VBF cuts, μ = m_H, fb)
99 of 100 proVBFH-cs jobs (1000201-1000300; job 83 still running, not used) and 8 × 20 VBFNLO seeds. Plain means over seeds, seed-scatter errors. Tools: `chsplit/chdiff.py` (table, output in `chdiff.out`), `chcombine.py` (combined tops `cs-W1..9.top`, `vbfnlo-{NC,CC}-{qq,qg,gq,gg}.top`).

| channel | proVBFH-cs | VBFNLO | cs − VBFNLO | pull |
|---|---|---|---|---|
| NC qq | 29.848 ± 0.105 | 29.163 ± 0.380 | +0.69 ± 0.39 | +1.7 |
| NC qg | 3.894 ± 0.037 | 3.879 ± 0.009 | +0.02 ± 0.04 | +0.4 |
| NC gq | 3.988 ± 0.204 | 3.886 ± 0.015 | +0.10 ± 0.21 | +0.5 |
| NC gg | −0.135 ± 0.003 | −0.136 ± 0.001 | +0.00 ± 0.00 | +0.4 |
| CC qq | 71.929 ± 0.197 | 72.335 ± 0.485 | −0.41 ± 0.52 | −0.8 |
| CC qg | 8.322 ± 0.086 | 8.319 ± 0.020 | +0.00 ± 0.09 | 0.0 |
| CC gq | 8.661 ± 0.545 | 8.326 ± 0.020 | +0.34 ± 0.55 | +0.6 |
| CC gg | −0.233 ± 0.004 | −0.235 ± 0.002 | +0.00 ± 0.01 | +0.3 |
| sum | 126.27 ± 0.98 | 125.54 ± 0.48 | +0.74 ± 1.09 | +0.7 |

- **Checks.** Sum of channels = W1 in every job (max 6e-6 fb); proVBFH-cs W1 = 126.27 ± 0.98 vs the earlier fixmh 126.19 ± 0.18 (1998 jobs), consistent. VBFNLO sum 125.54 ± 0.48 vs 125.22 ± 0.12 (800 jobs), consistent.
- **Cannot localise the offset.** The pilot total difference is +0.74 ± 1.09 fb, so the pilot is ~5x less sensitive than the known +0.97 ± 0.21 fb. No channel deviates by more than 1.7σ. NC qq (+0.69 ± 0.39) and the gq channels (+0.10, +0.34, with errors 0.2-0.55) are the largest, but all are compatible with 0 and with carrying the whole offset. The gg, qg channels agree to <= 0.02 fb and the qg channels are tightly constrained (errors 0.04-0.09): they cannot carry 1 fb. The offset lies in qq (NC/CC) and/or gq, i.e. the channels with large errors.
- **Asymmetry qg vs gq.** In VBFNLO qg = gq as it must be (3.879/3.886, 8.319/8.326). proVBFH-cs gq is noisy (per-job scatter 2.0 / 5.4 fb vs 0.4 / 0.9 for qg): a heavy tail, not a systematic shift.
- **Single job.** job-1000297 has σ(≥3j) = 215 fb (+90 fb above the mean; gq driven: NC gq/CC gq feed it, σ(≥4j) NC gq... 22.6 fb total). Without it the proVBFH-cs W1 mean is 125.37 ± 0.38, i.e. *below* the earlier fixmh value and 0.15 ± 0.4 above VBFNLO; CC gq becomes 8.14 ± 0.15, NC gq 3.80 ± 0.08, NC qq 29.78 ± 0.08, CC qq 71.79 ± 0.15. So one event shifts the 99-job mean by ~0.9 fb (about the size of the offset). Whether the +0.97 fb of the 1998-job mean is carried by such rare events in the gq (and 4-jet-like) tails is not decided here; it is a possibility to test (spikescan on fixmh W1 jobs).
- **σ(≥4 jets)** agrees in every channel (pulls ≤ 0.9); the sum differs by +0.75 ± 0.92, driven by the same job (gq channels). The 2-jet key is not meaningful in these NLO H+3j runs (VBFNLO prints the same value for 2, 3 jets; proVBFH-cs's is a negative real-minus-subtraction piece) and is not reported.
- **Does this contradict earlier conclusions?** No. The 7 Oct statement (offset in the exactly-3-jet part, same at both scales) is untouched; the pilot is just too weak to confirm or localise it.
- **Statistics needed.** Per-job cost is about 1.5 h on 1 core for proVBFH-cs (all nine weights at once), 1.5 h for VBFNLO CC and qg-type channels and 7.7 h for VBFNLO NC qq (one sample). To get errors of ~0.15 fb on the two qq differences: VBFNLO NC qq ~130 seeds (~1000 core-h), CC qq ~210 seeds (~320 core-h), proVBFH-cs ~170 jobs for CC qq (~260 core-h, the rest comes free through CHAN_MULTI); about 1.5-2k core-h in total. The gq channel of proVBFH-cs needs ~3000 jobs (~4.5k core-h) for 0.1 fb because of the heavy tail; a better route is a diagnosis of the tail events (spike replay) rather than brute force. Nothing submitted.
