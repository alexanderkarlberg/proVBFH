# proVBFH-cs: VBF (H, HH) at NNLO with line-by-line projection-to-Born

Status: design, stage 1 (NLO). Plan agreed with AK on 2026-09-29
(plan page: https://claude.ai/artifact/RyqC8MbpLzvqYpZ4jHusu9).

## Why

The exclusive part of `proVBFH` (VBF H+3j at NLO, POWHEG-BOX, FKS) is very
inefficient.
- POWHEG's FKS counterterms live at mapped kinematics whose q1, q2 differ
  from the real event's. The real event and its counterterm therefore
  project to different VBF Born points, and in every Born-level
  distribution their large weights cancel only statistically.
- There is one scale per event (hence μ0(pt,H) in 1506.02660) and one run
  per scale variation.
- Every point evaluates about 2400 tagged real flavour structures, each
  with its FKS regions.

In the factorised (structure-function) approximation each quark line is a
DIS process. Every exclusive piece is a product of per-line DIS pieces,
contracted through the VV→H (or HH) vertex. So the exclusive part can be
built line by line as in disorder/DISENT, with generation and dipole maps
that keep q_i, the Higgs and the other line fixed. Then every event and
all of its counterevents share the same VBF Born (x1, x2, q1, q2), and the
projected terms cancel event by event.

## Contributions

The notation (a, b) means line 1 at O(αs^a) and line 2 at O(αs^b).

| order | exclusive pieces | per line |
|---|---|---|
| NLO  | (1,0) + (0,1) | three-parton tree (q → q g, g → q q̄) |
| NNLO | (2,0) + (0,2) | DISENT's O(αs²): one loop + I, double real − dipoles, collinear (K, P) terms |
| NNLO | (1,1) | {Born, V + I, R, −D} on line 1 × the same on line 2 |

Each exclusive contribution with weight w is accompanied by a counterevent
with weight −w at the projected VBF Born kinematics:
`p_in,i = x_i P_i` and `p_out,i = x_i P_i − q_i`, where
`x_i = Q_i²/(2 P_i·q_i)`. The inclusive part (structure functions) is
computed and histogrammed at the VBF Born kinematics as in proVBFH.

## Event cycle (stage 1)

1. VBF Born point from the POWHEG VBF_H phase space used by
   `proVBFH-inclusive/src/phase_space.f` (7 variables). From it come x1,
   x2, q1, q2, p_H, and the Born quark momenta.
2. For each line i = 1, 2, a three-parton final state is generated from 3
   more variables (ξ_i ≥ x_i, z_i, φ_i) in that line's Breit frame, keeping
   q_i. The other line and the Higgs are not changed.
   Sampling (`line_radiation`): 1 − x_p and min(z, 1 − z) logarithmic down
   to the cutoff (soft and collinear limits). Optionally (`cs_hardfrac` h,
   2026-09-30) a second channel with probability h: ln x_p uniform in
   [ln x_B, 0] and z uniform, for hard emissions (both partons of the line
   hard, p_T² = Q² z(1 − z)(1 − x_p)/x_p, small x_p), which the
   logarithmic map hardly samples and which dominate the high-p_T tails;
   the weight uses the combined density. The same channel is used in the
   first step of the four-parton generator (`gen_four`, `four_weight`);
   its second emission has its own fraction (`cs_hardfrac2`, default 0:
   at NNLO a hard channel there costs more in the double-unresolved
   corners than it gains in the tails).
3. Weights:
   - w_i = |M_{H+3j}|² for emission from line i × PDFs(ξ_i, x_j) × flux ×
     phase space;
   - the Born-level counterevent is −w_1 − w_2, at the Born point of step 1.
4. The analysis receives the three-parton events (w_1, w_2) and one Born
   event (−w_1 − w_2). The inclusive NLO is a separate VEGAS run filling
   Born-level events. It can later be added to the same cycle.
5. VEGAS adapts on |w_1| + |w_2|. The signed sum per cycle is zero by
   construction.
6. All scale choices are computed at once, as per-event weight arrays
   (μ_i = Q_i per line by default; μ0(pt,H) as an option).

## Matrix elements (option A)

The VBFNLO-3.0 routines already in proVBFH, called directly without the
POWHEG machinery (only a few POWHEG headers are needed for their common
blocks):
- H+3j tree: `qqhqqj_born_channel`, whose `ansc(2:3)` are the gluon
  emitted from the upper and the lower line separately.
- H+4j real: `qqh2q2g_me_new.f`, `qqh4q_me_new.f`.
- H+3j one loop: `qqhqqj-virt.f` (box line).

Couplings come from `convert_coup.f`, and the flavour classes are those of
the existing code.

Later (option B): analytic BDK/MCFM helicity currents per line, contracted
numerically, with the kinematic structures shared across flavours.

## Validation plan (stage 1)

1. **Pointwise.** The per-line H+3j matrix element must equal the sum of
   the tagged pieces of the existing code's Born, flavour by flavour.
2. **Phase space.** The per-line generator must integrate to the analytic
   three-body volume, and its q_i must equal the Born's (to machine
   precision).
3. **Totals.** At NLO the P2B total equals the inclusive NLO total by
   construction; histograms of Born-level observables must agree with the
   inclusive ones when the exclusive part is switched off.
4. **Physics.** NLO distributions under VBF cuts against the current
   proVBFH at NLO and POWHEG VBF H jj at NLO (same scale choice,
   μ0(pt,H)).

## Speed study (stage 1)

Same observables, cuts, PDFs and scales as the current proVBFH at NLO:
- CPU time per point;
- statistical error per CPU hour on σ_VBF and on standard distributions
  (pt,j1, pt,j2, pt,H, Δy_jj), against the current code at equal CPU.

From these, estimate the NNLO cost on the thservs. If it does not fit,
tell AK, who will provide access to a Slurm cluster.

## Stage 3: the (1,1) contribution (design, 2026-09-30)

In the factorised approximation, the (1,1) part is line 1 at O(alpha_s)
times line 2 at O(alpha_s). Per line, the O(alpha_s) pieces are those of
NLO DIS with Catani-Seymour subtraction:
- V + I at the line's Born;
- (K + P) (x) f at the line's Born;
- R at the line's three-parton point;
- -D at the Born.

The line's IF and FI dipoles keep q_i. Their map is exactly the
projection to the VBF Born of stage 1, and the stage-1 radiation variables
(1 - xp, z) are the dipole variables (x, u). So every D_i lies at the
Born point of the line.

Products whose kinematics are the VBF Born on both lines cancel with their
own projection and are left out. That leaves, per point:

| event | kinematics | weight |
|---|---|---|
| E1 | line 1 radiated, line 2 Born | H3(R1) f1(xi1) [(V2 + I2) f2(x2) + (K+P)_2 (x) f2 (x2) - K2 f2(xi2) J2] |
| E2 | line 1 Born, line 2 radiated | the same with 1 <-> 2 |
| E3 | both radiated | H4(R1, R2) f1(xi1) f2(xi2) J1 J2 |
| Born | VBF Born | -(E1 + E2 + E3) |

Notation:
- H3(R1) is the H+3j tree with the extra parton on line 1 (stage 1).
- V2 is the vertex correction of line 2 in the H+3j virtual,
  CF (-8 - L^2 - 3 L), L = ln(mu^2/Q_2^2), which stage 2 removes from the
  (2,0) virtual.
- I2 and (K+P)_2 are the DIS I operator and K + P of line 2.
- K2 J2 is line 2's dipole (IF + FI kernels over 2 p.p x) times the
  radiation Jacobian, at line 2's three-parton point.
- H4(R1, R2) is the H+4j tree with one extra parton on each line: the
  class-12 entries of the real flavour list, via their tags.

The four-parton event E3 uses the three-parton points of both lines from
stage 1 (random numbers 8:13): no new phase space. In the single limits
(line 2 unresolved), E3 cancels the K2 term of E1. In the double limit,
E1 + E2 + E3 -> -H2 K1 K2 at Born-like kinematics, which cancels against
the projection event by event.

Validation, before the (1,1) runs:
1. Line-level dipoles against H3 in the singular limits of each line
   (q -> q g: x -> 1, u -> 0, u -> 1; g -> q qbar: u -> 0, 1).
2. **Structure functions.** With the dipole pieces of one line and the
   other line at Born level,
   int [(V + I) B f + B (K+P) (x) f + (R - D)]
   summed over both lines must equal sigma_NLO - sigma_LO of the
   inclusive (structure-function) code without cuts. This tests the DIS I
   and K + P constants. They share their structure with the stage-2
   `nlo2_ifin`/`nlo2_kp`, whose constants are otherwise tested only by the
   comparison with the old code.
3. E3 against H3 x line-2 dipole in line 2's singular limits, including
   the normalisation for two gluons on different lines (no symmetry
   factor, distinguishable lines).
