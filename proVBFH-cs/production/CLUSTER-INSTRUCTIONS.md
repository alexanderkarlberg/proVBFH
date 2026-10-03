# proVBFH-cs production on the cluster: instructions for Claude

Written 2 Oct 2026 (on thA371a) for the Claude session that AK starts on the
cluster. You start inside the proVBFH repository, on `main`. Read this file,
`CLAUDE.md` (repository overview, build, CI) and `proVBFH-cs/docs/DESIGN.md`
(what proVBFH-cs is and how it works) before doing anything else.

## 0. First: ask AK

Do not guess paths or policies; ask AK at the start (one message, numbered):

1. Scratch/work directory for builds and runs, and its quota.
2. Slurm account (if any) besides the partition `alma`; any rules on job
   arrays, core counts per job, memory requests, or I/O.
3. Compilers/modules: gfortran and g++ to use; are hoppet (>= 2.1), LHAPDF 6
   and FastJet installed (where), or should you build them (into the work
   directory)? Where to get the PDF sets `NNPDF30_nnlo_as_0118` (LHAPDF id
   261000) and `PDF4LHC21_40` (id 93100; check that this is the set of the
   HXSWG card, `lhans1 93100`).
4. Where the reference files are on the cluster (AK must copy them; they are
   not in the repository):
   - 1506.02660: `HH.top`, `11.top`, `22.top` (old proVBFH, paper files,
     combined with the old trimming; HH = μ_R = μ_F = μ0/2, 11 = μ0,
     22 = 2μ0). On thA371a: `proVBFH-cs/runs/ref-1506.02660/`.
   - HXSWG 13.6 TeV study: `nnlo-central.top`, `nnlo-min.top`,
     `nnlo-max.top` (and `nlo-*.top`) from
     `LHCHXSWG/vbf-higgs-wg/proVBFH/results/`; the raw per-seed data
     `13.6TeV_NNLO/{HH,11,22}.tgz` (150 MB each) only if needed. Ask AK
     whether the study's min/max are the 3-point (HH, 22) envelope.
5. Where to push results (which repository and branch; see section 7) and
   whether AK wants a page (artifact) with the plots.
6. The targets (section 5): confirm before the full submission, with the
   CPU estimate from the pilot.

Rules from AK that hold everywhere: push only when AK agrees; report at the
checkpoints (pilot result and CPU estimate before the full submission;
results before pushing); log the work in a notes file (section 7); when a
result contradicts an earlier conclusion, say so and correct the record.
Nothing from the private dis-1jet code may be used or committed.

## 1. Physics goal

New CS version of proVBFH (proVBFH-cs: line-by-line projection-to-Born with
Catani–Seymour dipoles, VBF H at NNLO in the factorised approximation),
in the two set-ups of the old proVBFH results, with scale variations, at
statistical errors comparable to the quoted (trimmed) errors of those
studies, and comparison plots new vs old for every histogram:

- **1506.02660** (13 TeV, NNPDF30_nnlo_as_0118, μ0(p_T,H), VBF cuts:
  p_T,j > 25 GeV, |y_j| < 4.5, m_jj > 600 GeV, Δy_jj > 4.5, opposite
  hemispheres, anti-k_t R = 0.4); analysis `p1506`.
- **HXSWG 13.6 TeV study** (PDF set 93100, μ0(p_T,H), p_T,j > 20 GeV,
  m_jj > 300 GeV, STXS bins); analysis `hxswg136`.
- Scale variations: symmetric 3 points, μ_R = μ_F = ξ μ0, ξ = 1, 1/2, 2
  (`cs_scales 3`; output weights W1 = (1,1), W2 = (1/2,1/2), W3 = (2,2)).
  The old results have the same three points (HH, 11, 22).

What is known (thA371a, Sep/Oct 2026; details in
`notes/2026-09-30-cs-p2b-stage3/README.md`, `notes/2026-10-02-scale-variations`):

- proVBFH-cs agrees with the structure functions (inclusive, (1,1)
  validation), with 1506.02660 at NLO, and with VBFNLO 3.0 (H+3j at NLO,
  ≥ 3 jets) at fixed and dynamic scale.
- The old proVBFH is high in ≥ 3 jets (about +5% in σ(≥ 3 jets) and in the
  3-jet distributions): nf = 4/5 mismatch in its `ffunc` and the missing
  initial-state FKS region of the NC pair graphs (unsubtracted, integrated
  down to a sampling-limited cutoff, heavy-tailed, grows like ln N). So
  differences in 3-jet observables are expected and explained; 2-jet
  observables should agree within errors (2-jet σ: −0.56%, −3.4σ in the
  first comparison, of which the ≥ 3-jet region carries all).
- NNLOJET (1802.02445) quotes 844 fb for σ(VBF cuts) at NNLO (no MC error)
  and agrees with the old code; this tension is open.
- The HXSWG comparison so far (8 seeds) agreed within 0.5%.

## 2. Code

Branch `2026-09-cs-p2b` (this file is on it). Build as in CLAUDE.md
(`proVBFH/Makefile.inc` from `./configure` in `proVBFH`, with the paths of
hoppet, LHAPDF, FastJet), then in `proVBFH-cs`:

    make ANALYSIS=p1506        # -> proVBFH-cs-p1506
    make ANALYSIS=hxswg136     # -> proVBFH-cs-hxswg136

`make` without ANALYSIS builds `proVBFH-cs` with proVBFH's own analysis
(not needed here). Before production, run the short regression checks of
section 6 with the cluster build.

## 3. Cards

In `proVBFH-cs/production/<setup>/`: `powheg-excl.input` (exclusive part,
`cs_part 2`, `cs_order 3`: NNLO), `powheg-incl.input` (inclusive part,
`cs_part 1`, `qcd_order 3`), `vbfnlo.input` (couplings, needed by both).
All have `cs_scales 3`; the exclusive cards `cs_hardfrac 0.3`. Every job
needs its own directory with `powheg.input` (a copy of the card with a
unique `iseed`) and `vbfnlo.input`, and runs the binary there without
arguments. A job does its own VEGAS warm-up (`ncall1` x `itmx1`) and then
stage 2 (`ncall2` x `itmx2`); output `pwg-EXCL-W{1,2,3}.top` (exclusive) or
`pwg-LO-W{1,2,3}.top` and `xsct-nnlo.dat` (inclusive), plus `run.log`.
Seeds: all distinct across all jobs of a set-up (e.g. exclusive
1000001..., inclusive 2000001...).

## 4. Cost per job (thserv18, Xeon Silver 4216, central scale only)

- 1506.02660 exclusive, card as given (200k x 2 warm-up, 1.6M x 3 stage 2):
  14,100 CPU-s = 3.9 CPU-h per job. HXSWG exclusive (600k x 3): about
  2.3 CPU-h (8 seeds: 66,580 CPU-s for 2.2M points each, 3.8 ms per point).
- `cs_scales 3` adds about 7% (7 points: 21%).
- Inclusive part: a few CPU-minutes per job; 20–50 jobs per set-up suffice
  (check its error against the exclusive one).
- Cluster cores are faster than thserv18; measure in the pilot.
- Wall time limit 24 h: size `ncall2`/`itmx2` so that a job takes well
  under it (e.g. 6–12 h); more, shorter jobs are better for the seed-scatter
  errors. Keep the warm-up per job (it makes the jobs independent).

## 5. Targets (to confirm with AK)

The per-bin CPU needed to reach the old studies' quoted (trimmed) errors,
estimated on thA371a from 16 (1506) and 8 (HXSWG) seeds, central scale,
thserv18 CPU-h:

| | needed per bin |
|---|---|
| 1506: σ(2 jets) | 1,200 |
| 1506: 2-jet distributions | median 1,000–3,900 per histogram, worst bin 27,000 |
| 1506: σ(≥ 3 jets) | 62,000 |
| 1506: 3-jet distributions | median 50,000–100,000, worst 2·10⁶ |
| 1506: σ(≥ 4 jets), 4-jet distributions | ≥ 5·10⁵ (out of reach) |
| HXSWG: all 946 bins | median 3,600; 16–84%: 860–60,000; total σ 3,000 |

Proposal discussed with AK (2 Oct): with O(10k) cores, all 2-jet bins and
σ(≥ 3 jets) at the quoted errors (about 60–70k CPU-h per set-up with
`cs_scales 3`), possibly the median 3-jet bins (about 100k CPU-h per
set-up). The quoted errors of the old code are artificially small (trimming
of heavy-tailed seed distributions), so matching them is a generous target.

Procedure: (1) pilot, about 200 exclusive jobs and 20 inclusive jobs per
set-up; (2) from the pilot, per bin C = CPU_pilot (err_pilot/err_target)²
(errors from the seed scatter), report the CPU needed for the targets
above to AK; (3) submit the full production after AK's go (Slurm job
arrays; at most 10k running, 30k queued; partition `alma`).

## 6. Checks before production (cluster build)

- Regression: one exclusive job of each set-up with `cs_scales` removed and
  `readingrid 1` on a fixed grid, compared with the same job with
  `cs_scales 3` (W1 must equal the run without variations, histogram by
  histogram; `sig incl cuts` is zero to rounding and differs at 1e-16 pb).
- Monitor the first jobs: running (`squeue`, CPU time growing), no `NaN`
  in run.log, `cs_spikes.dat` small.

## 7. Combination, plots, push

- Combine with `proVBFH-cs/tools/combine_parts.py` (per weight W1, W2, W3:
  exclusive files plus inclusive files as two `--part`s, `--error scatter`
  for the seed-scatter error of the mean; no trimming). Use `--strip -vbf`
  for the 1506 reference names. Scale band: envelope of W1..W3 per bin.
- Plots, one per histogram, both set-ups: upper panel new (central and
  scale band) vs old (central and band: HH/22 for 1506, min/max for
  HXSWG), lower panel new/old (ratio of the centrals with the new band and
  the statistical errors). Expect the 3-jet observables of the old code to
  be high (section 1).
- Push (after AK agrees) to a new branch, e.g. `2026-10-cs-production`:
  the combined .top files per weight and set-up, the plots (PDF/PNG), the
  scripts that made them, the Slurm scripts and cards used, and a notes
  file `notes/2026-10-cs-production/README.md` (what was run, CPU used,
  results, comparison). Not the per-job outputs, grids or logs (keep them
  on the cluster; give AK their location).
- Expected open item for the notes: the NNLOJET 844 fb tension.

## 8. Practicalities learned on thA371a

- A job's numbers depend on its seed only through `iseed`; identical seeds
  give identical results (duplicate work).
- `readingrid 1` skips the warm-up and reads `grids-excl.dat` (exclusive) or
  `grids.dat` (inclusive) from the job directory; only for checks.
- The VBFNLO virtual drops its boxes at near-collinear points and is
  slightly inconsistent for incoming antiquarks (both negligible; see the
  scale-variation notes); nothing to do.
- Kill only your own jobs, by job id; never by pattern.

## 9. Extra task (AK, 3 Oct): NLO VBF H+3j cross-checks with VBFNLO and POWHEG

### Why

proVBFH-cs and the old proVBFH differ significantly in 3-jet observables, but
not in 4-jet ones. Background, all in the notes:
- `notes/2026-09-30-cs-p2b-stage3/README.md`, sections from 1 Oct on;
- `notes/2026-09-29-cs-p2b-stage2/README.md`, "The public POWHEG-BOX-V2 VBF_HJJJ";
- `notes/2026-09-30-hxswg-comparison/README.md`.

State on thA371a: σ(≥ 3 jets, VBF cuts) at O(α_s²) (i.e. NLO H+3j), 1506.02660
set-up, in pb:

| | proVBFH-cs | VBFNLO 3.0 | old proVBFH |
|---|---|---|---|
| μ = m_H | 0.12550 ± 0.00067 | 0.12582 ± 0.00052 | 0.13302 ± 0.00097 |
| μ0(p_T,H) | 0.12795 ± 0.00116 | 0.12658 ± 0.00050 | 0.13324 (paper) |

That is the total rate only, with modest statistics. AK wants to be 100% sure:
1. that proVBFH-cs agrees with VBFNLO, also in **distributions** of 3-jet
   observables;
2. that the public POWHEG-BOX-V2 `VBF_HJJJ` at **fixed order** shows the same
   discrepancy as the old proVBFH. It shares the old code's problems:
   - **(a) unregulated logarithm.** The NC four-quark graphs (Z on the
     q q̄ pair) have no initial-state FKS region; their q || beam
     singularity is cut only by the Born generation cut `ptcut`.
   - **(b) pair-type swap.** In `compreal_hqqqq.f:529`, `kl` swaps u/d pairs.
   - **(c) NC gg pair type** in `real.f`.

   The old proVBFH (current source) has (a)–(c), and in addition an
   n_f = 4/5 mismatch in `ffunc` (+0.87e-3 pb on ≥ 3 jets at μ = m_H). The
   public POWHEG has n_f = 4/4, which is consistent, so it should lie
   slightly below the old code but well above proVBFH-cs and VBFNLO.

Run this after the main production is under way, or interleaved; it is small.

### Set-up: 1506.02660 cuts, fixed scale μ_R = μ_F = m_H first

A fixed scale removes the dynamic-scale ambiguity: VBFNLO's ID 20 is our
own patch, and POWHEG takes one scale per event. Then repeat with
μ0(p_T,H): `ID_MUF = ID_MUR = 20` in VBFNLO, `runningscales 1` in the other
two codes.

The cards are in `proVBFH-cs/production/h3j-crosscheck/` (copied from the
thA371a runs):
- `provbfh-cs/`: `powheg-excl.input` (cs_order 3, runningscales 0) and
  `vbfnlo.input`. The ≥ 3-jet observables at O(α_s²) come from the full
  NNLO run; the inclusive part does not contribute to them, so it is not
  needed here.
- `old-provbfh/`: `powheg.input` (qcd_order 3, testplots 1) and
  `vbfnlo.input`. Optional (AK may already have enough of it); ask AK.
- `vbfnlo/fixmh/`, `vbfnlo/dyn/`: VBFNLO 3.0 cards for process 110 (VBF H+3j
  at NLO, Catani–Seymour).
  - Settings: 13 TeV, NNPDF30_nnlo_as_0118, EWSCHEME 3 with the same G_F,
    M_W, M_Z, VBFHASB (b quarks in NC), anti-kt 0.4, the 1506 cuts.
  - VBFNLO keeps jets with |y| < 4.5 and p_T > 25 GeV (no veto), tags the
    two hardest, and requires ≥ 3 jets. That is the analysis's
    "sig(all VBF cuts 3 jets)".

### VBFNLO 3.0

- **Build.** Take 3.0 final, from the CERN LCG mirror; HepForge is behind an
  anti-bot page. Build the default processes with quad precision; the
  vbf,hjjj-only build does not compile. Apply
  `notes/2026-09-30-cs-p2b-stage3/tools/vbfnlo-3.0-scale20.patch` for
  scale ID 20.
- **Check before production.** At LO (NLO_SWITCH false, 2^20 × 4 points) and
  μ0 (ID 20) it must give 130.60 ± 0.44 fb for ≥ 3 jets with the VBF cuts.
- **Runs on thA371a.** 20 jobs × 2^23 points × 5 iterations per scale,
  different random.dat. Each job uses its own directory, with the cards
  and `vbfnlo --input=.`.
- **Distributions.** VBFNLO only fills its own built-in histograms
  (`histograms.dat`; switch TOP/GNU or data-file output on in `vbfnlo.dat`).
  - First list which of them match observables of
    `proVBFH-cs/analysis/p1506_analysis.f` exactly: same jet definition,
    same ordering in p_T, same binning or re-binnable. The third jet's p_T
    and rapidity and the H p_T in ≥ 3-jet events are the important ones.
  - If too few match, the clean way is to add a call in VBFNLO's
    histogramming routine that fills our observables from VBFNLO's
    momenta and weights. Discuss with AK before patching.
  - The total ≥ 3-jet rate must reproduce the table above within errors.

### POWHEG-BOX-V2 VBF_HJJJ (public, fixed order)

- **Get it.** svn r4135; on thA371a it is exported in
  `~/work/disorder-comparisons/powheg_vbf_hjjj_public`. Do not apply any of
  our fixes: the point is to run it as distributed.
- **Analysis.** Compile `proVBFH-cs/analysis/p1506_analysis.f` into it, as
  its `pwhg_analysis`. Both use POWHEG's `pwhg_bookhist-multi`, so the
  same histograms come out and compare bin by bin.
- **Fixed-order NLO.** Use `testplots 1`; the fixed-order distributions are
  `pwg-*-NLO.top` from the btilde integration (stages 1–2, no event
  generation needed). `bornonly 0`, `withnegweights 1`, a fixed scale via
  `runningscales 0`. Check in `Born_phsp.f`/`init_phys.f` how `muref` is
  set and make it m_H.
- **Parameters.** Use the same values as the old-provbfh card: beams,
  `lhans1/2 261000`, EW inputs in its `vbfnlo.input`, H mass and width,
  NNPDF30_nnlo_as_0118.
  - Note: it ignores ZWIDTH/WWIDTH in vbfnlo.input and computes
    2.5051/2.0950 GeV, where proVBFH reads 2.4952/2.141. This is a small
    effect. Either impose the widths in its `init_couplings.f` (and say so)
    or quantify the difference; `notes/powheg-comparison/powheg.md` has
    the recipe used for VBF_H.
- **Born generation cut.** Because of (a) the result depends on `ptcut`, the
  Born parton p_T generation cut (`#ptcut` in powheg.input, read in
  `Born_phsp.f`), logarithmically.
  - Run at two values, e.g. 1 GeV and 0.1 GeV; both must lie below the
    analysis cuts. Report the difference: it measures the unregulated log.
  - The old proVBFH generates with `phspcuts 1` and a different internal
    cut (`sigreal.f:1067`), so expect the same sign but not the same size.
- **Check before production:**
  - LO (`bornonly 1`): it must agree with proVBFH-cs's and VBFNLO's tree
    level (0.12031 ± 0.00115 / 0.11937 ± 0.00003 pb at μ = m_H).
  - The O(α_s²) ≥ 3-jet rate: expected to be close to the old code, about
    +5% above VBFNLO.

### Comparison and report

- **Totals.** σ(≥ 3 jets) and the exactly-3 / ≥ 4-jet split from all codes,
  at both scales, with seed-scatter errors (heavy tails: show the
  combination method; see `notes/2026-09-30-hxswg-comparison/README.md`).
- **Distributions.** Ratios to proVBFH-cs with errors, for all matching
  3-jet observables. Also the 4-jet ones, which should agree everywhere.
- **Expectation.** proVBFH-cs = VBFNLO within errors; public POWHEG ≈ old
  proVBFH (minus about 0.9e-3 pb from n_f), with a `ptcut` dependence.
  Anything else, report to AK before going on.
- **Tools.** `proVBFH-cs/tools/fixmh_compare.py` is what was used on
  thA371a for the totals.
- **Logging.** Log everything in the notes file of section 7, with a section
  of its own.
