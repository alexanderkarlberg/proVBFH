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
