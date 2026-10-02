# On-the-fly scale variations in proVBFH-cs (2 Oct 2026)

AK: "implement the on-the-fly scale variations ... a choice between the
symmetric 3 point or the full 7 point"; then "we need to do the beta0
shift instead" (of evaluating the virtual at every mu_R).

Usage and design: `proVBFH-cs/docs/DESIGN.md`, "On-the-fly scale
variations". Option `cs_scales 3|7` (exclusive and inclusive part), check
options `cs_scalecheck 1|2`.

## Profile (before)

NNLO exclusive part, 1506.02660 set-up, gprof: H+3j one-loop virtual
(`boxline_vg`, box/tensor integrals) about 50% of the CPU, H+3j trees
(`qqhqqj_born_channel`) 25%, H+4-parton reals 6%. Hence a per-point cache
for the matrix elements and no re-evaluation of the loop per scale.

## Validation

All runs on the stored grid of `runs/nnlo-p1506-h03/s1` (`readingrid 1`,
one iteration), so that every run sees the same phase-space points.
(Without `readingrid 1` each run redoes the VEGAS warm-up on its own
weights and the fixed-scale runs see other points: the first attempt was
stopped for this reason.) Scratch directories scaleval2, scaleval3,
scaleval4, inclval; comparison script `scaleval_compare.py` (relative to
the largest |value| of each histogram; `sig incl cuts` is zero to rounding,
1e-16 pb, and left out).

| test | result |
|---|---|
| no `cs_scales`, new vs pre-change code (exclusive, inclusive) | bitwise |
| W1 vs run without variations | identical (histograms), inclusive bitwise |
| inclusive Wk vs fixed-scale runs | bitwise; 1.5e-9, 5.5e-9 at μ_F = 2μ0 |
| exclusive, μ_F-only points (W5, W7) | exact; 8e-8 at μ_F = 2μ0 |
| exclusive, μ_R varied, virtual at each μ_R (first version) | exact; 9e-8 at μ_F = 2μ0 |
| exclusive, μ_R varied, β0 shift, 40k / 400k points | rms 0.06 / 0.06 of the MC error, max 0.3 |
| shift − direct with its own error (`cs_scalecheck 2`, 400k) | rms 1.00 over 309 bins; σ(2 jets) −1.1σ, +1.2σ |

The μ_F = 2μ0 residues: a fixed-scale run with facscfact 2 tabulates the
PDFs up to 2√s (`maxQval`), the on-the-fly run up to √s.

`cs_scalecheck 1` (shifted V + I against a direct evaluation at every
point): slope d(V+I)/d ln μ_R² = β0/(4π) = 0.610094 × Born exactly for
quark- and gluon-initiated radiating lines. Two deviations, both from
VBFNLO's virtual:
- incoming antiquark on the radiating line: 0.6107 (ū), 0.6099 (d̄), i.e.
  1e-3 of the slope, flavour dependent; the virtual routine's own Born is
  exactly 36 × the tree Born (also for antiquarks) and no gauge-check
  fallback occurs, so it is in the loop part; effect on V of order 1e-3
  Born, negligible;
- `qqhqqj-virt.f`: the boxes are switched off when any invariant of the
  line has p_i.p_j < 1 GeV² ("no boxes if there are less than 3 jets,
  these terms do not contribute anyway"); there V is the explicit logs
  plus triangles and its μ_R dependence is not the RG one (up to 20 ×
  Born at such points, which have Born ~ 1e5). In P2B these near-collinear
  emissions mostly cancel against their Born projections. Lowering VBFNLO's
  zero threshold in `C0t` (1e-3 GeV²) changes nothing.
The shift is the RG-consistent μ_R dependence (as DISENT/disorder treat
μ_R: β0 analytically, μ_F through K + P at each μ_F, cf. DISENT's
KPFUNS_SCL_VAR).

## Cost (thA371a, 4 jobs on 6 cores, 40k points)

no variations 44.5 s, 3 points 47.7 s (+7%), 7 points 54.1 s (+21%); with
the virtual evaluated at each μ_R the 7-point run cost +53%.

## Old proVBFH

The shared inclusive files (`proVBFH/src/inclusive/phase_space.f`,
`matrix_element.f90`, `incl_parameters.f90`) get the optional scale
factors and the weight loop; the default path executes the original
statements. The old proVBFH still builds (scratch copy).
