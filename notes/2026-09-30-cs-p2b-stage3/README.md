# proVBFH-cs stage 3: the (1,1) contribution (2026-09-30, in progress)

Design: `proVBFH-cs/docs/DESIGN.md`, section "Stage 3". Code on branch
`2026-09-cs-p2b` (local). Stage 2: `../2026-09-29-cs-p2b-stage2`.

## Building blocks

1. **VBF H + 2 parton Born** (`src/cs_born2.f`). VBFNLO has no H + 2
   parton routine in proVBFH, so this is the amplitude of its H + 3
   parton Born (`qqhqqj-born.f`) without the gluon: couplings `clr`,
   `b(6,V1,V2) xmw`, propagators with `xm2`, `xmg`, and squared currents
   16 (in1.in2)(out1.out2) for equal chiralities and 16 (in1.out2)(out1.in2)
   otherwise (in and out exchanged for antiquark lines).
   **Checked** (`cs_testborn2 1`, `runs/stage3-tests/testborn2.log`)
   against VBFNLO's H+3j in the soft-gluon limit, m_line/(eikonal x B2),
   for all 48 compatible class pairs, both lines, 2 points. The ratio is
   1 + O(delta): at delta = 1e-7 the worst deviation is 5.8e-7. This
   covers NC and CC, quarks and antiquarks, normalisation included.
2. **DIS I and K + P of a line** (`nlo2_ifin_dis`, `nlo2_kp_dis` in
   `cs_nlo2.f90`). `nlo2_kp` was split into a core routine with the
   colour-dependent factors as arguments; the stage-2 histograms are
   bit-for-bit unchanged (`runs/stage3-regress`). For the DIS line
   (T_a.T_b = -CF, one final quark): kqf = 3/2, lsc = ln(muF^2/Q^2).
3. **Line-level dipoles and validation mode** (`cs_order 11`,
   `line_nlo_point`):
   - Born level: (V + I) B f + B (K + P) (x) f, with V = CF (-8 - L^2 - 3L).
   - At the line's three-parton point: R f(xi) - D f(xi).
     - q -> q g: D is the IF + FI dipoles.
     - g -> q qbar: D is the two IF dipoles, with Born -qbar -> q or
       -q -> qbar.

   The dipole map is the stage-1 projection, so every dipole's Born is at
   the VBF Born point.

## Validation against the structure functions (running)

Integrated over everything without cuts, the sum over both lines must equal
sigma_NLO - sigma_LO of the inclusive code (`order_min 2`, `qcd_order 2`).
Same scale, mu0(pt,H) (`runningscales 1`).

**Result** (`runs/stage3-tests/incl-nlo`, `line-nlo`; 13 TeV, NNPDF30,
cutoff 1e-6):
- structure functions: sigma_NLO - sigma_LO = -0.10296 +- 0.00030 pb;
- line pieces from the dipoles: -0.1040 +- 0.0044 pb.

They agree (0.2 sigma). The error is +-0.06 in units of (alpha_s/2 pi)
times the Born coefficient, small against the constants being tested
(CF pi^2/3 = 4.4, 2 gamma_q = 4). This validates the DIS I and K + P and the
line dipoles, and with them the shared structure of the stage-2
`nlo2_ifin`/`nlo2_kp`.

Line-level limit test (`cs_testlines 1`, `runs/stage3-tests/tl/testlines.log`):
- H+3j against the line dipoles: R/D -> 1 as sqrt(lambda) for all classes
  and both lines, in x -> 1 and u -> 0 (and u -> 1 for g -> q qbar).
- Worst deviation at lambda = 1e-7: 8e-4.

## The four-parton (1,1) event and a new matrix-element bug (2026-09-30)

- **Entry map** (`nlo11_init`): every pair of line states (8 quark classes
  q -> q g, 4 gluon classes g -> q qbar, per line) is mapped to a class-12
  real entry (tags 1, 2 on the lines, the extra partons at legs 6, 7 in
  either order). All 54 compatible pairs are found once the NC gg entries
  are accepted: in those, line 1 has the antiquark in slot 4 and the quark
  as its extra parton (swap flag).
- **Limit test** (`cs_testlines11 1`, `runs/stage3-tests/tl11`). One line
  sits at a moderate three-parton point; the other goes to its singular
  limits. The H+4j must approach H3 x that line's dipoles.
  - **Normalisation:** with a factor 2 for two gluons, the ratio tended to
    2. The tagged entries already give the full matrix element for
    gluons on different lines, so the factor is 1.
    **This corrects my first guess** (a symmetry factor 1/2 in setreal).
  - **NC gg -> H q qbar q' qbar' with different pair types failed:** the
    ratio tended to 0.78 or 1.27 = (gL^2 + gR^2)_u/(..)_d or its inverse,
    while same-type pairs passed.
- **Bug (proVBFH `real_vbfnlo.f`, `compreal_hjjj`; identical in the public
  POWHEG-BOX-V2 VBF_HJJJ `real.f:1134`).** In the NC branch
  "g g ---> qb1 q2 q1 qb2 H" (flavours 4 < 0 < 5 with 4 = -6), `ftype(2)`
  is taken from leg 6, which is the first line's quark. Since
  k = -2 ftype(1) - ftype(2) + 7, both q qbar pairs get the first line's Z
  couplings when their types differ. The other gg branches take the
  second line's type correctly.
  - Every NC gg entry of proVBFH's list has this form, and so do the
    public code's NC gg entries 201-216 (listed with
    `~/work/disorder-comparisons/powheg_vbf_hjjj_public/tools/list_gg_public.f`),
    e.g. 201 (cbar d c dbar) and 209 (ubar d u dbar).
  - **Fixed in a marked copy** `proVBFH-cs/src/real_vbfnlo.f` (ftype(2)
    from leg 5). With it, all 216 checks of the limit test pass (worst
    deviation 7.6e-4 at lambda = 1e-7, converging as sqrt(lambda)).
  - Effect: only the gg-initiated NC four-quark reals with different pair
    types (the NLO H+3j of POWHEG; in proVBFH also the old NNLO exclusive
    part). Not quantified yet; it is small (gg luminosity, two g -> q qbar
    splittings).
- **Public build test of the new bug** (tools in
  `../2026-09-29-cs-p2b-stage2/tools/is_limit_gg_public.f` and
  `list_gg_public.f`). The process is gg -> H Qbar d Q dbar with dbar
  collinear to beam 2, which must factorise onto g d -> H Q d Qbar x
  TR [1 - 2x(1-x)].
  - As distributed: c = 0.780185 for Q = u (types differ) and 1.000000
    for Q = d (control).
  - With the one-line fix in the public real.f: 1.000000 for both.
  - Added to the bug report as problem 3.

## Validation of the (1,1) part against the structure functions

`cs_order 13` integrates, without cuts,
sum_{c1,c2} B2 (A1 - J1 D1)(A2 - J2 D2) + E1 + E2 + E3.
- A_i is line i's Born-level O(alpha_s) factor.
- D_i is its dipoles per Born class at its three-parton point.
- The first sum holds the Born-kinematics terms that the projection
  removes from the exclusive part.

The reference is the (1,1) product of the lines' O(alpha_s) structure
functions: the inclusive code with `order_min 3`, `qcd_order 3` and the
new `incl_only11 1` (marked copy `src/matrix_element.f90`, which keeps only
the Fx1(:,2) Fx2(:,2) term; default unchanged).

- First run (`runs/stage3-tests/incl-11`, `line-11`, 1.35M points):
  - inclusive: -0.0017585 +- 0.0000120 pb;
  - dipole pieces: -0.00242 +- 0.00110 pb.

  Consistent, but only at the 60% level: the (1,1) term is a small sum of
  large pieces.
- Large run on thserv18 (`runs/stage3-valid11`, 30 seeds x 3.3M points,
  commit 487a89d) for a test at about 10%: running.

The (1,1) part alone with VBF cuts (`runs/stage3-timing`, cutoff 1e-4,
40k points): 0.19 ms per point with `cs_no20` (the stage-2 evaluations are
skipped), sum |w| stable, no spikes.

## Reference for the full validation (1506.02660, Table I)

13 TeV, NNPDF30_nnlo_as_0118, mu0(pt,H):

| order | sigma (no cuts) [pb] | sigma (VBF cuts) [pb] |
|-------|----------------------|-----------------------|
| LO    | 4.032                | 0.957                 |
| NLO   | 3.929                | 0.876                 |
| NNLO  | 3.888                | 0.844                 |

The NNLO statistical error is about 0.1% with VBF cuts.

Stage 1 gave NLO 0.8763 +- 0.0020 pb. The NNLO correction with VBF cuts
is -0.032 pb, of which (2,0) + (0,2) exclusive is -0.0225 +- 0.0022 pb
(stage 2); the (1,1) exclusive part and the inclusive NNLO come on top.

The old code's known issues shift its result by a few permille:
- nf = 4 in ffunc: -0.18% of sigma(VBF cuts);
- the unsubtracted initial-state region: 1.5e-4 per e-fold of an unknown
  effective cutoff;
- the pair-type swaps: negligible to small.

So agreement is expected at about the 0.3% (0.0025 pb) level.

(1,1) validation result (`runs/stage3-valid11`, `stage3-valid11b`, 2 x 30
seeds x 3.3M points; seed-scatter errors, which agree with VEGAS's):

| sample    | (1,1) from the dipole pieces [pb] | pull vs inclusive |
|-----------|-----------------------------------|-------------------|
| seeds 1-30 (2001-2030) | -0.001652 +- 0.000063 | +1.65 |
| seeds 31-60 (3001-3030) | -0.001668 +- 0.000063 | +1.41 |
| all 60    | -0.001660 +- 0.000044             | +2.14 |

The inclusive reference is -0.0017585 +- 0.0000120 pb.

**Not settled: a 5.6 +- 2.6% (1e-4 pb) difference.** Both halves deviate
in the same direction. Ruled out so far:
- the Born normalisation: LO 4.0354 +- 0.0064 (cs_born2, `cs_order 10`)
  against 4.0329 +- 0.0026 pb (inclusive);
- W flavour counting: Hoppet sums u, c / d, s for W, as my classes do;
- the PDF/alpha_s set-up: LHAPDF and 3-loop alpha_s, independent of
  qcd_order.

The (1,1) term is roughly the product of the two lines' O(alpha_s)
pieces. An error of about 0.004 pb in one line's NLO correction would
shift it by 1e-4, and the first line-NLO validation (+-0.0044) cannot
exclude that. Next: a 10x more precise line-NLO validation
(`runs/stage3-validnlo`, 30 seeds; inclusive reference with 4 more seeds).

**High-precision line-NLO validation** (`runs/stage3-validnlo`, 30 seeds x
3.3M points; `runs/stage3-tests/incl-nlo-seeds`, 4 more inclusive seeds):
- dipole pieces: -0.102618 +- 0.000330 pb (seed scatter; VEGAS 0.000399);
- structure functions: -0.102623 +- 0.000068 pb.

Pull 0.02, precision 0.3% of the NLO correction. This excludes a per-line
error as the source of the (1,1) difference (it would need about 0.004 pb)
at more than 10 sigma. The DIS I and K + P and the line dipoles are right;
stage 2 shares the I and K + P routines. What remains specific to (1,1)
is the H+4j with one extra parton per line (E3), tested so far only in
its singular limits. Otherwise the difference is a fluctuation (p = 0.03).

Size of (1,1) (inclusive code, `runs/stage3-tests/split`):
- without cuts: O(alpha_s^2) = -0.04089 +- 0.00008 pb, of which (1,1) is
  -0.00176, so (2,0) + (0,2) is -0.0391 (ratio 4.5%);
- inclusive part with VBF cuts (Born-level events): -0.00780 pb, of which
  (1,1) is -0.00024 pb (3%).

**E3 away from its limits, against NNLOJET**
(`~/cernbox/disorder-comparisons/vbf_nnlojet/harness_e3.f90`, log
`e3_20.log`). The VBFNLO class-12 H+4j (as used for E3, with the gg fix)
is compared with the one-gluon-per-line ("non-adjacent") piece of
NNLOJET 1.0.2's C2g0VBF. That piece is rebuilt from `C2g0VBFnadj_vbf`,
`coupling_vbf` and `propagator_vbf` exactly as in C2g0VBF, at 20 generic
points (no parton soft or collinear).
- R3/(4 pi R_ctl) is constant to 9 digits, with R_ctl the H+3j control
  ratio:
  - 8/3 for q q (NC s c, u d, d d with identical flavours, CC u d -> d u);
  - 1 for g q (NC u and d pairs);
  - 3/8 for g g (NC same and mixed pair types, CC).
- Each incoming gluon gives a factor 3/8, the quark/gluon colour average;
  the constants are NNLOJET's normalisation conventions. E3's absolute
  normalisation is fixed by the limit tests. So E3 is right everywhere.
- A harness problem on the way, mine: NNLOJET caches amplitudes by their
  six indices, with a flag that flips at every `clearC2gVBF_vbf`. An
  entry not refreshed in the previous epoch becomes valid again two flips
  later and returns an amplitude of an earlier point. NNLOJET's driver
  computes the same combinations at every event, so it is not affected.
  The harness mixed combinations and gave point-dependent ratios (0.03,
  310) until all combinations were computed in every epoch.

So neither the line pieces (0.3% validation) nor E3 explain the 2.1 sigma
of the (1,1) check. The third sample (60 seeds) will show whether it is a
fluctuation.

## To do after validation (agreed with AK, 2026-09-30)

1. **Speed comparison with the old proVBFH:** run its NNLO exclusive part
   for a fixed short time on the same machine and compare CPU x error^2
   per observable. This replaces the estimate of "1e4 one-week runs".
2. **Warm-up strategy:** adapt only the 7 Born dimensions on a cheap
   integrand (LO or inclusive with cuts), freeze the radiation dimensions
   (`jfreeze`) and read that grid in the exclusive run; compare with the
   present adaptation of all 13 dimensions on sum |w|, at the same CPU,
   error^2 per observable. The present warm-up is 2 x 200k of 5.2M points
   per job (8%), with no grid stage or grid combination.
3. **Scale variations** as per-event weight arrays, one run for all scales.

**(1,1) validation: passed** (third sample `runs/stage3-valid11c`, 60 seeds
6001-6060):

| sample | (1,1) from the dipole pieces [pb] | pull vs inclusive |
|--------|-----------------------------------|-------------------|
| seeds 2001-2030 | -0.001652 +- 0.000063 | +1.65 |
| seeds 3001-3030 | -0.001668 +- 0.000063 | +1.41 |
| seeds 6001-6060 | -0.001832 +- 0.000059 | -1.22 |
| all 120         | -0.001746 +- 0.000038 | +0.31 |

The inclusive reference is -0.0017585 +- 0.0000120 pb, so the two agree
to 0.7 +- 2.2%. The earlier 2.1 sigma was a fluctuation: the new seeds
lie on the other side, and the scatter errors of 30-seed samples of this
heavy-tailed integrand are themselves uncertain. Together with the
line-NLO validation (0.3%) and the E3 comparison with NNLOJET, the (1,1)
part is validated. Next: the (1,1) cutoff study with VBF cuts
(`runs/stage3-cut`).

## (1,1) cutoff study with VBF cuts (2026-09-30)

Set-up: `runs/stage3-cut`, commit 487a89d, thserv18, nice 10. (1,1)
exclusive part alone (`cs_order 3`, `cs_only2 1`, `cs_no20 1`). Otherwise
as in the stage-2 study: 13 TeV, VBF cuts, NNPDF30, mu0(pt,H), 10 seeds x
5.2M points per cutoff, 48 CPU-min per job, no spikes.

| cutoff | sigma(VBF cuts) [pb] | error, seed scatter | error, VEGAS |
|--------|----------------------|---------------------|--------------|
| 1e-4   | -0.00461             | 0.00138             | 0.00152      |
| 1e-5   | -0.00335             | 0.00093             | 0.00114      |
| 1e-6   | -0.01159             | 0.01174             | 0.01133      |

- **No cutoff dependence:** over all histograms, chi2/n = 227.8, 253.1 and
  223.9/228. The inclusive pt,H histogram (no VBF cuts) runs a little
  high (60-74/50), as in stage 2; its event/counterevent cancellations
  are the strongest.
- **Size:** the (1,1) exclusive part is -0.0037 +- 0.0008 pb (1e-4 and
  1e-5 combined). With the inclusive (1,1) part (-0.00024) it is -0.0040
  pb with VBF cuts, about 1/7 of (2,0) + (0,2).
- **Cost:** sum |w| is 45, 118 and 257 pb (2.5 x stage 2), but the
  points are cheaper: +-0.0014 pb in 8 CPU-h at 1e-4 (stage 2: +-0.0022
  in 25 CPU-h).
- **Rough NNLO total from the pieces:**
  - inclusive O(alpha_s^2) with cuts -0.0078, (2,0) + (0,2) exclusive
    -0.0225 +- 0.0022, (1,1) exclusive -0.0037 +- 0.0008;
  - total correction -0.0340 +- 0.0023 pb, against -0.032 in 1506.02660;
  - sigma_NNLO(VBF cuts) = 0.842 +- 0.003 pb (with stage 1's NLO
    0.8763 +- 0.0020), against 0.844 in the paper.

## NLO cutoff check and a correction to the cutoff chi2 (2026-09-30)

NLO cutoff check (`runs/nlo-cut-check`): exclusive NLO part with the
stage-1 inputs and seeds (6 seeds, 2M points), thA371a, about 14 CPU-min
per job.

| cutoff | sigma(VBF cuts) [pb] | error, seed scatter | error, VEGAS |
|--------|----------------------|---------------------|--------------|
| 1e-4   | -0.05412             | 0.00141             | 0.00082      |
| 1e-5   | -0.05452             | 0.00166             | 0.00119      |
| 1e-6   | -0.05515             | 0.00155             | 0.00162      |

No cutoff dependence at NLO either.

**Correction.** The Higgs-only histograms (`sig incl cuts`, `ptH-incl`,
`yH-incl`, no VBF cuts) are zero in the exclusive part up to rounding
(|value| < 1e-18 pb): the maps keep q_1 and q_2, so every event and its
Born counterevent have the same Higgs momentum. Their pulls therefore
measure rounding noise only. The "all histograms" chi2 above and in the
stage-2 logbook included their 71 bins. The statement above that the
inclusive pt,H histogram "runs a little high ... its event/counterevent
cancellations are the strongest" is wrong: there is nothing to cancel.
`tools/cutoff_compare.py` now leaves out bins that are zero up to
rounding (below 1e-12 of the largest value). Recomputed chi2/n over the
histograms with VBF cuts:

| study | 1e-4 vs 1e-5 | 1e-4 vs 1e-6 | 1e-5 vs 1e-6 |
|-------|--------------|--------------|--------------|
| NLO (`nlo-cut-check`)      | 117.5/155 | 141.8/155 | 91.5/154  |
| (2,0)+(0,2) (`stage2-cut2`) | 183.5/155 | 149.7/155 | 147.0/155 |
| (1,1) (`stage3-cut`)        | 144.8/154 | 157.1/154 | 133.3/154 |

The conclusion (no cutoff dependence) stands. The runs at different
cutoffs share their seeds, so they are correlated and chi2/n somewhat
below 1 is expected. The largest value, 183.5/155 for stage 2 at 1e-4
vs 1e-5, is 1.6 sigma above the mean.

## Full NNLO against 1506.02660 (2026-09-30)

`runs/nnlo-full` (30 seeds, cs_order 3, cutoff 1e-5, 5.2M points, 3.0
CPU-h per job, thserv18, no spikes) + `runs/nnlo-incl` (4 seeds), against
AK's paper files `runs/ref-1506.02660/11.top` (mu_0; `HH`, `22` = 0.5,
2 mu_0). Combined with `tools/combine_parts.py --strip=-vbf --error max`
(error = larger of VEGAS and seed scatter, bin by bin):

| | proVBFH-cs [pb] | 1506.02660 files [pb] | difference |
|---|---|---|---|
| sig(VBF cuts), NNLO | 0.83854 +- 0.00188 | 0.84383 +- 0.00046 | -0.0053 (-2.7 sigma, -0.63%) |
| inclusive part with VBF cuts | 0.92249 +- 0.00084 | 0.92258 +- 0.00001 (`sig incl cuts`) | +0.1 sigma |
| exclusive part | -0.08395 +- 0.00169 | -0.07875 +- 0.00046 | -0.0052 (-3.0 sigma) |

- The paper's `sig incl cuts` is filled for every event before cuts; in
  a P2B run events and counterevents cancel in it, so it is the inclusive
  part with VBF cuts (phspcuts at the Born). Ours agrees.
- nnlo-full's exclusive part agrees with the sum of the separately run
  pieces at 1e-5 (NLO -0.05452, (2,0)+(0,2) -0.02284, (1,1) -0.00335:
  -0.08071 +- 0.00325; -0.9 sigma).
- Distributions (10 with VBF cuts): chi2 211.6/154; after rescaling ours
  by 1.0063, 141.5/154. The difference is a normalisation, not a shape.
- The O(alpha_s^2) exclusive part is therefore about -0.029 (ours) vs
  -0.024 (paper, using its NLO 0.876 and our NLO inclusive part).
- Known issues of the old code, by the version: the paper files are from
  2018-02-21, after proVBFH 1.1.0 (2018-02-07, fixed H+3j virtual) and
  before 1.1.1 (2018-11-20), which introduced `nf = st_nlight` in the
  explicit logs with `ffunc` still at nf = 4. In 1.1.0 both are 4, and
  the nf terms of `ffunc` and of the explicit gamma_g logs cancel exactly
  when they use the same nf, so the -0.0016 pb of `estimate-nf4` (the
  4/5 mismatch) does not apply to the paper. The missing ISR region
  (issue 2) adds +0.000126 pb per e-fold of the effective collinear
  cutoff (`estimate-coll-e2`): the right sign, and 25-40 e-folds would
  explain the difference, but the old code's effective cutoff (with the
  trimming) is not known. kl swap: 7e-6. NC gg pair type: not yet sized.
- Built `proVBFH-cs-p1506` with the paper's analysis (proVBFH 1.1.2,
  `analysis/p1506_analysis.f`): 23 histograms, including the 3- and 4-jet
  rates. Same seed and card: the 2-jet results are bitwise identical to
  the current analysis. sigma(>= 3 jets) has no counterevents and no
  inclusive part, so it compares the H+3j NLO directly.

## Size of the NC gg pair-type bug (`runs/estimate-gg`, cs_estimate 5)

1506.02660 set-up (nnlo-full card), E3 with the original line minus E3
fixed, as event and Born counterevent, unsubtracted; thserv21, 6 seeds x
5.2M points, about 1.7 CPU-h per job.

| cutoff | sig(VBF cuts 2 jets) [pb] | sig(VBF cuts 3 jets) [pb] |
|---|---|---|
| 1e-4 | +6.1e-6 +- 0.1e-6 | +8.9e-6 |
| 1e-6 | +9.2e-6 +- 0.1e-6 | +1.34e-5 |

(first 3 and 5 seeds). Slope about 7e-7 pb per e-fold of the cutoff:
negligible against the 0.0053 pb difference with 1506.02660 (and against
the missing ISR region's 1.3e-4 per e-fold). Cutoff 1e-8 does not run:
all six jobs stop at the start with "cs_dipoles: no transverse basis"
(the guard for a degenerate transverse basis); the cutoff studies so far
went to 1e-6.

## 60 seeds against 1506.02660, and the 3- and 4-jet rates (19:35)

`runs/nnlo-p1506`: 30 more seeds (iseed 7101-7130; 14 on thserv18, 16
on thserv09) with the paper's analysis (`proVBFH-cs-p1506`), no spikes;
inclusive part rerun with that analysis (`runs/nnlo-incl-p1506`, same
cards and seeds as nnlo-incl). `tools/combine60_1506.py` (seed-scatter
errors):

| sig(VBF cuts) [pb] | proVBFH-cs | 1506.02660 files | difference |
|---|---|---|---|
| 2 jets, p1506 (30 seeds) | 0.84197 +- 0.00255 | 0.84383 +- 0.00046 | -0.7 sigma (-0.22%) |
| 2 jets, nnlo-full (30) | 0.83854 +- 0.00156 | | -3.3 sigma (-0.63%) |
| 2 jets, all 60 | 0.84026 +- 0.00161 | | -2.1 sigma (-0.42%) |
| >= 3 jets (p1506) | 0.12919 +- 0.00157 | 0.13324 +- 0.00006 | -2.6 sigma (-3.0%) |
| >= 4 jets (p1506) | 0.01772 +- 0.00124 | 0.01688 +- 0.00001 | +0.7 sigma |

- The two halves agree with each other (1.2 sigma): the -3.3 sigma of
  the first 30 seeds was partly a fluctuation (correction of the
  -2.7 sigma quoted above with the "max" errors).
- The 3-jet deficit (-0.0040 pb) is as large as the 2-jet one (-0.0036):
  the exactly-2-jet part agrees (+0.0004), the difference sits in the
  H+3j region, where the old code's missing ISR region (issue 2) acts
  (three hard jets and a parton collinear to the beam). At +1.26e-4 pb
  per e-fold (2-jet rate), about 32 e-folds would be needed; POWHEG
  sampling the collinear variable down to about machine precision would
  give about 37. Coefficient for the 3-jet rate: `runs/estimate-coll-p1506`.
- Distributions: 2-jet ones fine; 3-jet ones follow the paper within
  our (much larger) errors, apart from the normalisation. min{rap(j1,j3),
  rap(j3,j2)} has chi2 151.8/45, of which 96 from the last bin, where the
  paper file has 2.4e-17 +- 2.5e-18 (a combine_runs.f zero-bin artefact)
  and ours is 0: 55.8/44 without it.
- Our errors for the 3- and 4-jet observables are 40-190 times the
  paper's (the old code needs no subtraction there at LO and its errors
  are trimmed); the hard-emission channel should help.

### The difference is the shape of the missing ISR region (19:50)

`runs/estimate-coll-p1506` (issue 2, cs_estimate 2, with the paper's
analysis): per e-fold of the collinear cutoff, sig(VBF, 2 jets)
+1.258e-4 +- 0.019e-4 pb (event minus counterevent; reproduces
estimate-coll-e2), >= 3 jets +3.103e-4 +- 0.017e-4 (the events only), 4
jets 0. Fit of (paper - ours) = Delta x (estimate), bin by bin over the
23 histograms of `nnlo-p1506/total.top` (errors of ours and the paper):

| histograms | bins | best Delta [e-folds] | chi2 at Delta = 0 | at best |
|---|---|---|---|---|
| 2-jet (14) | 211 | 30.2 +- 4.7 | 185.5 | 143.8 |
| 3-jet (5) | 102 | 16.0 +- 1.9 | 157.8 | 90.1 |
| all distributions | 344 | 18.7 +- 1.9 | 444.4 | 347.6 |

One parameter, the effective number of e-folds over which the old code
integrates the unsubtracted IS-collinear singularity of the NC pair
graphs, reproduces both the size and the shape of the differences, with
a shape computed independently in proVBFH-cs; the 3-jet chi2 drops by 68.
The 2-jet and 3-jet values of Delta differ by 2.7 sigma; a common
Delta = 19 gives chi2/n about 1. The 4-jet histograms (not affected by
issue 2) have chi2 107/34 from our under-sampling (most bins 10-30% low
with small errors, a few spikes with large errors): the second emission
of gen_four is sampled logarithmically only.

### More seeds (23:50): the difference is established

`nnlo-p1506-h03` (16 seeds, hard channel h = 0.3 on both emissions, iseed
7301-7316; unbiased like the others, see the line-NLO check in the
hxswg notes) pooled with nnlo-p1506 (30) and, for the 2-jet rate,
nnlo-full (30); weighted means of the seed-scatter results:

| sig(VBF cuts) | proVBFH-cs [pb] | 1506.02660 files | difference |
|---|---|---|---|
| 2 jets (76 seeds) | 0.83913 +- 0.00130 | 0.84383 +- 0.00046 | -0.56% (-3.4 sigma) |
| >= 3 jets (46 seeds) | 0.12729 +- 0.00107 | 0.13324 +- 0.00006 | -4.5% (-5.5 sigma) |
| >= 4 jets (46 seeds) | 0.01764 +- 0.00055 | 0.01688 +- 0.00001 | +1.4 sigma |

(Update of the 60-seed -2.1 sigma above.) The 2-jet difference is
-0.0047, the 3-jet one -0.0059: again all of it in the >= 3-jet region
(exactly 2 jets: +0.0012 +- 0.0017). In terms of the missing ISR
region: 19 e-folds from the 3-jet rate, 37 from the 2-jet rate (ratio
of the differences 0.79 +- 0.25 against 0.41 predicted, 1.5 sigma).

### Correction (2026-10-01): NNLOJET agrees with the old proVBFH, not with us

arXiv:1802.02445 (Cruz-Martinez, Gehrmann, Glover, Huss; NNLOJET, antenna
subtraction, factorised VBF) uses exactly the 1506.02660 set-up (NNPDF3.0,
mu_0(pt,H), M_W 80.398, Gamma_W 2.141, anti-kt 0.4, the VBF cuts) and
finds for VBF-2j: LO 957, NLO 877, NNLO 844 +9 -9 (scale) fb (Table 2;
no MC error quoted), and "once these corrections are applied [AK, private
communication] in [proVBFH], we find excellent agreement with our results
for the fiducial cross section, Table 2, and all distributions considered
in [1506.02660]" - i.e. with the 2018 paper files we compare with.
proVBFH-cs: NLO 876.3 +- 2.0, NNLO 839.1 +- 1.3 fb.

So an independent code without POWHEG's missing ISR region agrees with
the old numbers. The interpretation above (the old code high because of
issue 2) is withdrawn: most likely proVBFH-cs is low, by about 0.6% in
the 2-jet rate and 4.5% in the >= 3-jet rate, i.e. in the O(alpha_s^2)
H+3j region. The fit of the differences to the missing-region estimate
does not discriminate: any deficit in the H+3j region has a similar
shape. (Issue 2 itself stands as a code fact, but its numerical effect on
the old results must be small.)

Checked since: the VBFNLO one-loop H+3j used here is in the CDR
(4 pi)^eps/Gamma(1-eps) normalisation assumed by our I operator
(harness_virt: NNLOJET - VBFNLO = (pi^2/6)(2 CF + CA/2) to 1e-13, the
predicted shift from NNLOJET's e^{-eps gamma} normalisation); the colour
factors of nlo2_ifin and the K+P coefficients (kqf, lsc) look right on
inspection. The (2,0)+(0,2) integrated dipoles have no integrated check
(the exclusive part vanishes for Higgs-only observables), so they remain
the prime suspects, together with anything that only shows at >= 3 jets.
The public NNLOJET v1.0.2 cannot run VFH/VFHJ (no driver channels; the
VFH/VFHJ entries are commented out), so a direct NNLOJET run needs its
authors.
