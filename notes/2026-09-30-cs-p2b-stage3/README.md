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
