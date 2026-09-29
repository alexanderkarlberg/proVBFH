# proVBFH-cs stage 1: line-by-line P2B at NLO (2026-09-29)

Plan: `proVBFH-cs/docs/DESIGN.md`, and the plan page
https://claude.ai/artifact/RyqC8MbpLzvqYpZ4jHusu9. Code on branch
`2026-09-cs-p2b` (not pushed).

## What was built

`proVBFH-cs/` computes the exclusive part of VBF H at NLO without the
POWHEG-BOX machinery.
- For each VBF Born point (the inclusive code's phase space), radiation
  is generated on each quark line. The line's momentum transfer q_i, the
  Higgs and the other line stay fixed (`src/cs_kinematics.f90`,
  covariant).
- Each three-parton event with weight w has a counterevent −w at its own
  VBF Born point.
- Matrix elements: proVBFH's VBFNLO H+3j Born (`compborn_hjjj`), split
  by the radiating line. It is evaluated once per pair of flavour classes
  (up/down type, quark/antiquark/gluon, NC/CC) and summed with the PDFs.
- The inclusive part is proVBFH's own (`run_inclusive`, `cs_part 1`).
- Inputs are proVBFH's `powheg.input` and `vbfnlo.input`, plus the
  `cs_*` keys.

## Checks

- `tests/test_kinematics`, for power and logarithmic sampling:
  - masses, momentum conservation, q fixed, and x_p and z reproduced, all
    to ≤ 4e-12;
  - the sampled phase-space volume equals the analytic one.
- Flavour classes against an explicit loop over all flavour combinations
  (`cs_flavcheck 1`): 2e-15 at 200 points.
- NLO with VBF cuts. The inputs are those of
  `notes/2026-09-24-validation-against-papers/runs/proVBFH-NLO`: 13 TeV,
  NNPDF30_nnlo_as_0118, M_W = 80.398, Γ_W = 2.141, μ0(pt,H). The run is
  `proVBFH-cs/runs/nlo-vbfcuts-log2`, not committed:

  |                         | proVBFH-cs         | current proVBFH (8 seeds) | paper |
  |-------------------------|--------------------|---------------------------|-------|
  | σ(VBF cuts) [pb]        | 0.8781 ± 0.0031    | 0.8740 ± 0.0014           | 0.876 |
  | distributions, χ²/n     | 0.6–1.7 (11 histograms) | –                    |       |
  | CPU                     | 1.09 h (6 × 6.6M exclusive points + inclusive) | 4.27 h | |

  The inclusive part alone gives the Born-level cross section with VBF
  cuts, 0.9314 ± 0.0012 pb. The reference's "sig incl cuts", 0.9306 pb,
  is the same quantity: `phspcuts 1` means it is not an inclusive cross
  section.

## Speed at NLO

CPU × error², relative to the current code:
- σ(VBF cuts): 1.3 times worse;
- M_jj: 4.6 times better;
- R(j1,j2): 3.4 times better;
- pt,j2: 1.8 times better;
- pt,H: about equal;
- pt,j1: 3.5 times worse.

The reference errors come from `combine_runs` with outlier removal.

So the two are comparable at NLO, as expected. At NLO the current code's
exclusive part is also only the H+3j tree, generated from the VBF Born
with a q-preserving map (`br_real_phsp_isr_new`). Its events and
counterevents already share their Born point. The FKS mismatch, between
the H+4j real and its counterterm, appears only at NNLO, which is where
the line-by-line approach should gain.

## Technical cutoff

Cutoff 1e-6 instead of 1e-8, with the same seeds and statistics
(`runs/nlo-vbfcuts-log2-cut6`):
- The exclusive part agrees with the 1e-8 run in every VBF-cut
  histogram (χ²/n between 0.1 and 1.2): no cutoff dependence.
- The errors are 0.57–0.84 times smaller for the same CPU (about 1.05 h
  exclusive), because fewer points sit deep in the region where event
  and counterevent cancel exactly.
- σ(VBF cuts) = 0.8763 ± 0.0020 pb (current proVBFH 0.8740 ± 0.0014).
- CPU × error² relative to the current code:
  - σ(VBF cuts) 1.8 times better;
  - pt,j1 1.7, yj1 2.9, pt,j2 3, R(j1,j2) 6 and M_jj 7 times better;
  - pt,H about equal.

The default is now 1e-6.

## Problems found on the way (corrections of earlier numbers)

1. **VEGAS adaptation over all of phase space.** The first run adapted
   on |w| everywhere: 0.888 ± 0.006 pb, with errors 4 to 20 times the
   reference's. Now, as with proVBFH's `phspcuts`, the cut decisions for
   the events and the Born come first. The matrix elements are evaluated
   only where they matter, and points where nothing passes are skipped
   (`cs_phspcuts`, on by default). This is exact for all VBF-cut
   histograms, and the exclusive part does not contribute without cuts.
2. **Power sampling (DISENT's, npow = 2) fails at large Q.** At
   Q ≈ 2.5 TeV (a quark going to the opposite hemisphere), emissions
   that are hard for the jet cuts have 1 − x_p ~ W²/Q² ~ 5e-4. Example:
   W = 57 GeV, a 7 GeV gluon at wide angle, Δy_jj shifted below 4.5.
   There the weights grow like 1/√(1−x_p): single points carried
   hundreds of pb, and 1000 points carried 98% of the variance.
   Logarithmic sampling in 1 − x_p and min(z, 1−z) flattens this exactly
   (`cs_npow 0`, now the default).
3. **Histogram normalisation bug (mine).** `pwhg_bookhist-multi`
   normalises by the number of `pwhgaccumup` calls. The integrand did not
   accumulate points that returned early, which inflated the exclusive
   part by 1/(1 − fraction skipped). The factor was about 2 with
   generation cuts, and a few per cent before them. So the intermediate
   results 0.888, 0.858 and 0.809 pb, and the error comparisons between
   them, are not valid. Every point is now accumulated.
