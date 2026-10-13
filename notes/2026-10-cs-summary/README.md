# proVBFH-cs: summary of the study of September-October 2026

proVBFH-cs computes VBF Higgs production at NNLO, fully differential, with a
Catani-Seymour-based exclusive part on top of the structure-function
(inclusive) calculation of proVBFH. This note replaces the six logbooks of the
study (stage 1-3, HXSWG comparison, scale variations, production). They, and
all study tools, are kept in the git tag `cs-study-2026-10`. The code is in
`proVBFH-cs/` (see its `README` and `docs/DESIGN.md`).

## Method

In the factorised approximation each quark line is a DIS process, and every
exclusive piece is a product of per-line DIS pieces contracted through the
VV -> H vertex. The inclusive part is the structure-function calculation of
proVBFH. The exclusive part is built line by line, as in DISENT or disorder:
for a VBF Born point (x1, x2, q1, q2) radiation is generated on one line with
q_i, the Higgs and the other line fixed. The line's initial-final and
final-initial dipole maps keep q_i as well. Every real event and every
counterevent therefore projects to the same VBF Born point, and the projected
terms cancel event by event rather than only statistically, which is what
makes the FKS-based exclusive part of POWHEG so inefficient at NNLO. With
(a, b) for the orders in alpha_s of line 1 and line 2, NLO needs (1,0) + (0,1),
NNLO adds the DIS-like (2,0) + (0,2) (one loop, I, K + P, double real minus
dipoles) and the product (1,1) of the two lines' O(alpha_s) corrections. The
matrix elements are those of VBFNLO 3.0 contained in proVBFH. All scale
variations are computed on the fly (renormalisation-scale dependence of the
virtual by its beta0 shift, factorisation-scale dependence through K + P and
the PDFs), at about +7% CPU for three scale points. One run does the inclusive
part and then the exclusive part, and adds the histograms.

## Sampling problems and fixes

The method is exact, so every problem found was one of sampling, or a bug in
the matrix elements of the old code (see the validation).

| problem | cause | fix |
|---|---|---|
| large errors at large Q (stage 1) | power sampling of 1-x_p (DISENT's) fails when hard emissions have 1-x_p ~ 5e-4; 1000 points carried 98% of the variance | logarithmic sampling in 1-x_p and min(z,1-z), cutoff `cs_cutoff` |
| spikes in the NNLO tails | single hard double emissions from one line at small x_p, sampled with very low density in the four-parton generator | hard channel for the second emission, `cs_hard2mode 1`, `cs_hardfrac2 0.7`: 1.7 times less CPU for sigma(>= 3 jets) and 4.8 times for sigma(>= 4 jets) in a 300-job pilot, no outliers in the 2-jet part |
| heavy tail in the gq/qg channels of sigma(>= 3 jets) | a few jobs with +90 to +190 fb, each from one hard four-parton real at small x of one beam (tau = x1 x2 between 5e-4 and 2.5e-3, below the Born threshold of the cuts, where VEGAS keeps only 1/27 of the uniform density) | small-x channel for the Born phase space, `cs_smallx 0.3`, `cs_smallx_x0 0.03` (exact mixture with the VEGAS map): error at equal CPU 6 times smaller for sigma(>= 3 jets) and 26 times for sigma(>= 4 jets), 33% less CPU per job, means unchanged |

The heavy tail was first suspected to be a missing region in VBFNLO or a bug in
the dipoles. Neither was true: every large point is a finite, hard real with a
resolved dipole map, i.e. variance. Dedicated runs of both codes in x < 0.02
agree (3.06 to 3.13 fb for sigma(>= 3 jets) against 3.07 in VBFNLO).

## Validation

Within the code: unit tests of the kinematics, the small-x channel and the
four-parton generator (`make check`). During the study the real matrix elements
were checked against the dipoles in all singular limits for all 104 flavour
groups, and the beta0 shift of the virtual against evaluating it at every
scale (rms 1.00 of the error over 309 bins). Two bit-level checks in
`proVBFH-cs/validation`: the combined run reproduces the earlier two-job set-up
(V1), and the cleaned-up code reproduces the code before the clean-up (V3).

Against the old proVBFH. At NLO with VBF cuts (13 TeV, 1506.02660 set-up)
proVBFH-cs gives 0.8763 +- 0.0020 pb, proVBFH 0.8740 +- 0.0014 pb, paper
0.876 pb. At NNLO the 2-jet rate of the final production, 0.8405 pb, is 0.4%
below the paper file (0.84383 +- 0.00046 pb); in a 30-seed test the difference
was a normalisation, not a shape. The sigma(>= 3 jets) is different: 126.87 fb against 133.24 fb.
This is the old code, not the new one. The old exclusive part has three problems
in the H+3j region: the pair type is swapped in the NC four-quark real, the
initial-state singularity of the NC pair graphs has no FKS region, and the
H+3j virtual mixes nf = 4 and nf = 5. All three were found with direct tests
(the first also against NNLOJET's matrix elements). POWHEG with the fixes
gives 127.4 +- 2.9 fb, as distributed 133.4 +- 3.0 fb. NNLOJET (1802.02445)
quotes 844 +- 9 fb for the 2-jet rate at NNLO.

Against the HXSWG study (13.6 TeV, PDF4LHC21, STXS cuts). Early runs agree
with the study's results: LO sigma (ptj > 20 GeV) 2.4779 +- 0.0012 pb against
2.4790 +- 0.0001 pb, NLO 2.1826 +- 0.0031 against 2.1794 +- 0.0003 pb, NNLO
2.098 +- 0.010 against 2.090 (trimmed) or 2.096 (plain) pb. The reference
at NNLO depends on how its seeds are trimmed, which is why plain and
symmetric-trim merges are kept in `production/reference/hxswg136`.
The final production (10,872 of 11,000 NNLO jobs) is in
`production/results/final-2026-10-07/hxswg136`.

The final code with its default card (V2: 100 NLO and 100 NNLO jobs, about
200 core-h) reproduces the final p1506 production within errors at both
orders: at NNLO 0.8414 +- 0.0011 against 0.8405 pb for the 2-jet rate, 127.3
+- 0.5 against 126.9 fb for sigma(>= 3 jets), all pulls below 1.5
(`proVBFH-cs/validation/v2_results.txt`).

## Final results for 1506.02660 (p1506)

sigma with all VBF cuts, >= 2 jets, in pb. Scale points (mu_R, mu_F) =
(1,1), (1/2,1/2) and (2,2) times mu_0(p_T,H). Errors from the seed scatter,
no trimming (11,000 exclusive and 1,000 inclusive jobs at NNLO). Files:
`proVBFH-cs/production/results/final-2026-10-07`.

| order | (1,1) | (1/2,1/2) | (2,2) | band |
|---|---|---|---|---|
| LO | 0.95722 +- 0.00001 | 1.02334 | 0.89806 | +6.9%, -6.2% |
| NLO | 0.87690 +- 0.00003 | 0.85959 | 0.88436 | +0.9%, -2.0% |
| NNLO | 0.84047 +- 0.00018 | 0.83007 | 0.85046 | +1.2%, -1.2% |

At NNLO sigma(>= 3 jets) is 126.87 +- 0.17 fb and sigma(>= 4 jets) 17.12 +- 0.15 fb.

## The H+3j comparison with VBFNLO

The sigma(>= 3 jets) of proVBFH-cs is the H+3j cross section at NLO (mu = m_H,
1506.02660 cuts, fb), which VBFNLO 3.0 also computes. The production with the
small-x channel (749 jobs) gave 126.185 +- 0.163 fb against 125.221 +- 0.120
fb, i.e. +0.96 +- 0.20 fb (+0.8%), all of it in the exactly-3-jet part (sigma(>= 4
jets) agrees). To locate it, both codes were run per piece (proVBFH-cs against
VBFNLO 3.0); this larger run supersedes the 749-job number:

| piece | proVBFH-cs | VBFNLO |
|---|---|---|
| LO | 119.349 +- 0.147 | 119.34 |
| V + I + K + P | 30.183 +- 0.072 | 30.21 |
| real - dipoles | -23.713 +- 0.146 | -24.35 +- 0.15 |
| total | 125.82 +- 0.22 | 125.22 |

The total differs by +0.48 to +0.5% (2 to 2.4 sigma), uniform in the
distributions (ptj3, yj3, y*j3, min rap(j,j3): chi2 12.1/15, 14.6/18, 37.3/24,
61.8/46), and independent of the cutoff (`cs_cutoff` 1e-4 and 1e-6 against
1e-5: +0.38 +- 0.47 and -0.21 +- 0.58 fb). Pointwise, all ingredients agree: the
Born (up to the Higgs and boson widths, 0.05 fb), V + I to 2e-14, K + P to 1e-4,
the real matrix elements to 3e-4 or better, and the jet and cut logic in 777 four-parton points, and
VBFNLO's dipoles cancel the real in all limits. VBFNLO's own printed real-dip,
-23.31, is biased by the inverse-variance weighting of each job's VEGAS
iterations: iterations hit by a negative spike also have a large error and get
little weight. The plain mean over 800 jobs is -24.35 +- 0.15. Status: small open discrepancy.

## How to run

```
cd proVBFH && ./configure && make        # proVBFH first: proVBFH-cs uses its Makefile.inc and sources
cd ../proVBFH-cs && make
cd example && ../proVBFH-cs               # reads powheg.input and vbfnlo.input
```

One run gives `pwg-W<k>.top` (all histograms, scale points k) and
`xsct-<order>.dat`, with `cs_order` 1, 2 or 3. For production, `production/slurm/`
has `make_jobs.sh` (one directory per seed), `run_array.sh` (Slurm array) and
`combine_and_plot.sh`. The full list of options is in the README.
