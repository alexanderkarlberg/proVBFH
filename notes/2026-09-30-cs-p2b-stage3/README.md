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
