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
