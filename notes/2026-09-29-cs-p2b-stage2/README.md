# proVBFH-cs stage 2: the (2,0) and (0,2) contributions (2026-09-29, in progress)

Plan: `proVBFH-cs/docs/DESIGN.md`; stage 1 in `../2026-09-29-cs-p2b-stage1`.
Code on branch `2026-09-cs-p2b` (not pushed, not yet committed for stage 2).

Stage 2 is the NLO correction of the line that radiated (DISENT's NLO
2+1-jet structure per VBF line): the H+4j real with both extra partons on
the line, its Catani-Seymour dipoles (FF, FI, IF; one incoming parton), the
one-loop H+3j with the loop on the line, the I operator and K + P.

## What is built

- `src/cs_dipoles.f90`: CS maps and their inverses for a line (q fixed),
  the kernels (spin correlations with POWHEG's B^{mu nu}), and the
  four-parton phase space: the three-parton configuration of stage 1 plus
  one FF or FI splitting, the outputs in random order, with the
  multichannel weight from the 24 ways of reaching a point. Its seven
  random numbers are VEGAS dimensions that are not adapted (the density is
  computed analytically; `integration.f` gained `jfreeze`).
- `src/cs_glue2.f`: proVBFH's real flavour list with its line tags
  (`init_processes`), `setreal` for one entry, the H+3j Born with the extra
  parton on one line and its spin correlations, the one-loop H+3j per line.
- `src/cs_nlo2.f90`: the real entries of each line (1018 per line, class 11
  or 22) in four line structures (S1 g -> Q Qbar g; S2 q -> Q g g; S3
  q -> Q + pair from a gluon; S4 q' -> Q Qbar q' with the pair tag on the
  incoming quark), grouped by the value of their matrix element and
  dipoles (checked numerically at start-up: 52 groups per line); the
  dipole tables; I (finite part in the virtual's normalisation) and K + P
  (DISENT's KPFUNS, MSbar, convolution by Gauss quadrature).
- `src/cs_exclusive.f90`: `cs_order 2` adds, per line, V + I + K + P at the
  three-parton event and the real event with its six dipole counterevents,
  all projected to the VBF Born point.

## Checks

- `tests/test_four`: maps and inverses to 1e-13; weights consistent;
  integrals of test functions with gen_four agree with flat (RAMBO)
  sampling of the line's three-body phase space (2% statistical errors,
  all pulls < 2.1 in 24 comparisons).
- Real against dipoles in every singular limit (`cs_testlimits 1`), for all
  104 groups: R/D - 1 below 1.5e-3 at lambda = 1e-6 for S1, S2, S4 and CC
  S3; g -> q qbar collinear limits approach 1 more slowly (0.994 at
  lambda = 1e-7). This includes the spin correlations (g -> gg, g -> q qbar,
  q -> g initial state) and required no colour-average factor for the
  flavour-changing initial-state dipoles (as expected with averaged matrix
  elements and CS's averaged kernels).
  - The first versions of the test were misleading in two ways, both test
    artefacts: its three-parton seed points were themselves near a
    singular limit (logarithmic sampling, cutoff 1e-12), and it kept
    dipoles whose own three-parton Born is singular (dropped by the
    three-parton cutoff in the integration).
- Direct factorisation checks: the H+3j tree (VBFNLO) in the
  initial-state collinear limit against the W-fusion Born, and the H+4j
  real against the H+3j, both constant to 1e-4 in x.

## Problems found in the matrix elements proVBFH uses

1. **NC four-quark real: pair type swapped** (`../proVBFH/src/exclusive/
   compreal_hqqqq_new.f:559`, `kl = k+4*ftype(7)-4`). The list selected by
   kl has up-type pairs in 1..4 and down-type pairs in 5..8, but ftype = 2
   is up-type. This affects the graphs in which the Z couples to the pair
   (the incoming quark emits a t-channel gluon that fuses with the Z into
   the pair), which VBFNLO includes in the NC q -> Q + pair entries.
   Evidence: their initial-state collinear limit, averaged over the
   azimuth, is 0.79217 (d-type pair) and 1.26237 (u-type pair) times
   P_gq(x) times the gluon-initiated Born, constant in x; the two are exact
   reciprocals. With `kl = k+4*(2-ftype(7))` the ratio is 1.00000, as for
   the CC entries (S4). proVBFH-cs uses a marked copy
   (`src/compreal_hqqqq_new.f`); proVBFH is not changed.
   **Confirmed independently with NNLOJET** (v1.0.2 core library, harness
   `~/cernbox/disorder-comparisons/vbf_nnlojet/harness_4q.f90`, Z fusion,
   sW^2 = 1 - MW^2/MZ^2, fixed seed; logs `4q_fixed.log`, `4q_orig.log`):
   - control, q-initiated H+3j: VBFNLO/NNLOJET = 3.586228e-3 for the line
     flavours (s,c), (d,u), (u,d), (c,s) at three points, identical to all
     digits;
   - g-initiated H+3j, g -> Q Qbar: 1.344836e-3 (3/8 of the control, the
     colour average) for Q = d and Q = u at all points, with NNLOJET's slot
     i1 the outgoing antiquark (the other assignment varies);
   - four quarks s c -> s c Q Qbar H with the outgoing s collinear to the
     beam (pT 1e-3 GeV, 1 - x = 0.09), where the graphs with the Z on the
     pair dominate (NNLOJET's E0g0VBF with the incoming pair quark):
     corrected VBFNLO/NNLOJET = 4.50671e-2 (Q = d) and 4.50663e-2 (Q = u),
     i.e. 4 pi times the control (the extra g_s^2); the original gives
     3.514e-2 and 5.780e-2 (u/d = 1.645).

2. **No FKS region for that singularity in the POWHEG part.** The same
   graphs are singular (R ~ 1/lambda) when the outgoing quark of the
   incoming flavour becomes collinear to the beam. The merged Born would
   need pair tags (11) on the q qbar pair, which the Born list does not
   have (tags 0, 1, 2 only), so POWHEG's region finder cannot produce the
   region; the old exclusive NNLO part integrates this singularity
   unsubtracted. **Confirmed with POWHEG's own region finder**
   (`tools/list_regions.f`: proVBFH's init_processes and finalize_tags, then
   genflavreglist from find_regions.f, output `tools/list_regions.log`): all
   500 NC entries q(1) -> q(1) + pair(11) on line 1 have exactly one
   region, emitter 6 radiating 7 (the pair, g -> q qbar), and no
   initial-state region; all 160 entries with the pair tag on the incoming
   quark (S4, CC only) have their initial-state region. So the CC version of
   these graphs is subtracted in the old code and the NC one is not. In
   proVBFH-cs these NC entries get the two IF dipoles q -> q with the
   gluon-initiated Born (dipole structure 5); their integrated form is the
   QG term of K + P.
3. **One-loop H+3j: scale of the ERT function — withdrawn (my test was
   wrong).** I first concluded that the radiating line's finite virtual in
   `qqhqqj-virt.f` has no mu_R dependence, because the logs of `ffunc` (FHZ's
   F, arXiv:0710.5621 eq. 2.24, called with `st_muren2`) cancel the explicit
   logs exactly, and `cs_testvirt` showed V + I violating the RG constraint
   by (CF + CA/2) ln^2 mu_R terms. That test used a three-parton point that
   was itself nearly singular (logarithmic sampling with cutoff 1e-9 gave
   1 - x_p ~ 1.4e-6, the same mistake as in the first limit test). There,
   VBFNLO's box-line gauge check fails and it falls back to a constant
   (`cvirtH3j`), which is where the double logs of the radiating line
   actually sit; the explicit logs and F cancel between themselves as they
   should. At moderate points (linear sampling, and generic points in the
   NNLOJET harness) V + I is linear in ln mu_R^2, with slope 4.167 instead
   of b0 = 3.833: F has nf = 4 hard-coded while the explicit logs (and
   alpha_s) use nf = 5. That small inconsistency, (alpha_s/2 pi) B (1/3)
   times a log, is all that remains of problem 3. The proposed "fix"
   (F at Q_i^2) was wrong and is not used; the estimate runs based on it
   (`runs/estimate-ffunc-*`) were stopped and are not valid.
   - **Comparison with NNLOJET's one-loop H+3j** (v1.0.2 core library;
     harness `~/cernbox/disorder-comparisons/vbf_nnlojet/harness_virt.f90`,
     logs `virt*.log`). NNLOJET splits it into leading colour (`C1g1VBF`),
     subleading colour (`Ct1g1VBF`, with the q qbar dipole of the radiating
     line) and the nf part (`Ch1g1VBF`); the full result is
     N C1g1 - Ct1g1/N + nf Ch1g1. (My first comparison used C1g1 alone with a
     normalisation of 2 CF, which is why its double logs came out as -3 CF
     instead of -(2 CF + CA/2).) With all three:
     - NNLOJET's 1/eps^2 coefficient is Catani's, -(4 CF + CA) = -8.3333, at
       every point; its nf part has zero finite part and the 1/eps pole 1/3
       per flavour (the nf part of -gamma_g);
     - CC u d -> d u g H (6 points, mu_R = 25 ... 400 GeV and Q_1) and NC with
       the same quark type on both lines, q- and g-initiated:
       NNLOJET - VBFNLO = (pi^2/6)(2 CF + CA/2) + (1/6)[ln(|t|/mu^2) +
       ln(|u|/mu^2)], averaged over the two gluon placements with their Borns,
       to 1e-11 to 1e-13. The constant is the change of normalisation
       (e^{-eps gamma} against 1/Gamma(1-eps), times the 1/eps^2
       coefficients); the logs are VBFNLO's nf = 4 in `ffunc` against nf = 5
       in its explicit logs: with a consistent nf the nf-dependent finite part
       is zero, as in NNLOJET, and VBFNLO has a spurious
       (alpha_s/2 pi) B (1/6)[ln(mu^2/|t|) + ln(mu^2/|u|)], t, u the gluon's
       invariants with the quarks of its line. This is the slope 1/3 of the
       mu_R test above. **Size** on the old results (`cs_estimate 4`: the
       (1,0) events and their projections weighted by that term, at
       mu0(pt,H) as in the old runs; `proVBFH-cs/runs/estimate-nf4`, 5 seeds x
       6.6M points): delta sigma(VBF cuts) = -(1.55 +- 0.07) e-3 pb, -0.18% of
       sigma_NLO(VBF cuts) (about 4% of the NNLO correction); in the
       distributions -0.1% to -0.4% where significant (p_t,j2 110-120 GeV:
       -0.35%, M_jj 3.8-4 TeV: -0.4%).
     - NC with different quark types on the two lines: a further difference of
       up to 0.2 in units of alpha_s/(2 pi) B, independent of mu_R, zero for
       equal types, antisymmetric in the types ((s,c): -0.191, (c,s):
       +0.191; g -> d dbar with c: +0.108, g -> u ubar with s: -0.108), and
       odd under the reflection p_y -> -p_y of the event (checked). So it
       multiplies the parity-violating coupling combination
       g_L1^2 g_R2^2 - g_R1^2 g_L2^2, i.e. the LR and RL helicities, which are
       equal at tree level: a T-odd absorptive term, in which the two codes
       differ (both have T-odd parts, and they agree for LL, RR and for
       same-type lines). It integrates to zero for any reflection-symmetric
       observable and cuts; all of proVBFH's histograms and the VBF cuts are
       (phi(j1,j2) is |Delta phi|). Which code is right in this term is not
       settled: VBFNLO's one-mass box `D01m_fin` (the "corrected" file) agrees
       with the standard analytic continuation (Ellis, Zanderighi) in every
       sign region; interposing NNLOJET's complex log with its conjugate
       also changes its T-even parts and is inconclusive.
       **Corrected (later on 2026-09-29, see "The T-odd term settled"
       below): the codes do not differ. My harness gave NNLOJET the line
       flavours in the wrong order (`nfC1` must be the flavour of the line
       that does not radiate); with the right order NNLOJET and VBFNLO
       agree to 2e-12 for all processes, and an independent unitarity
       calculation reproduces VBFNLO's T-odd term.**
     - **Correction on the origin of the nf = 4 (later on 2026-09-29):** the
       mismatch is proVBFH's own. The public POWHEG-BOX-V2 VBF_HJJJ has
       nf = 4 in both places (the explicit gamma_g logs, line 641 of its
       `vbfnlo-files/qqhqqj-virt.f`, and `ffunc`), and with a consistent nf
       the finite part does not depend on nf. proVBFH's copy
       (`../proVBFH/src/exclusive/vbfnlo-files/qqhqqj-virt.f:692-693`)
       comments out `nf = 4d0` in the explicit logs and uses
       `nf = dble(st_nlight)` there, but leaves `ffunc` at 4. Checked: a
       build with nf = 4 in both places (`qqhqqj-virt-nf44.f`,
       `harness_virt_nf44`, log `virt14_nf44_swap.log`) agrees with NNLOJET
       to 2e-12, as the nf = 5/5 build does. So the -0.18% above is an
       effect in proVBFH's results only. The fix in proVBFH-cs (`st_nlight`
       in `ffunc` as well) is equivalent to the public code.
     - VBFNLO's gauge-check fallback (box line replaced by a constant when the
       gauge amplitude exceeds 0.1 of the Born one) never triggered at these
       points; it does at nearly singular points, which misled the first
       mu_R test.

## Size of problems 1 and 2 on the old results

In the old code the NC graphs with the Z on the pair are integrated with
their initial-state collinear singularity unsubtracted, so their
contribution grows like ln(1/cutoff), the cutoff being whatever limits the
sampling near the beam. Its coefficient is the collinear limit itself:
alpha_s/(2 pi) [P_gq (x) sum_q f_q] times the NC gluon-initiated H+3j
Born, at the three-parton events of stage 1 with their Born projection
(`cs_estimate 2`, swapped pair type as in proVBFH; `cs_estimate 3`, correct
type; `proVBFH-cs/runs/estimate-coll-e2`, `-e3`, 5 seeds x 6.6M points each,
13 TeV, VBF cuts, mu0(pt,H), NNPDF30_nnlo_as_0118).

| per e-fold of the cutoff        | swapped (proVBFH)       | correct                 |
|---------------------------------|-------------------------|-------------------------|
| sigma(VBF cuts) [pb]            | (1.258 +- 0.025) e-4    | (1.318 +- 0.027) e-4    |
| relative to sigma_NLO(VBF cuts) | 1.44e-4                 | 1.51e-4                 |

- Problem 1 alone (the u/d swap) changes this by -7e-6 of sigma_NLO per
  e-fold: negligible.
- Problem 2 is small but formally divergent: 1.5e-4 of sigma per e-fold, so
  0.15-0.3% for an effective ln(1/cutoff) of 10-20, compared with an NNLO
  correction of about 4% after VBF cuts (abstract of 1506.02660). In the
  distributions it is at most about 1e-3 per e-fold (the central y_j1 bin,
  M_jj above 5 TeV), elsewhere 1e-4 to 4e-4. In a Monte Carlo run it would
  show up as a slowly growing, spiky contribution rather than as a clean
  shift. The effective cutoff of the old runs is not known.

## Checking the nf fix against NNLOJET (2026-09-29)

With `nf = dble(st_nlight)` in `ffunc` of proVBFH-cs's marked copy of
`qqhqqj-virt.f`, NNLOJET - VBFNLO is the normalisation constant
(pi^2/6)(2 CF + CA/2) = 6.853892 at every mu_R, for CC and for NC with the
same quark type on both lines, to 1e-13. V + I then has the RG slope b0
to 4e-11.

## The T-odd term settled (2026-09-29, late)

1. **Helicities isolated** (`harness_virt.f90`, s c -> s c g H, Z fusion;
   VBFNLO's `clr` and NNLOJET's `enL`/`enR` set to keep one helicity per
   line; dumps `heli_points.dat`, `heli_vbfnlo.dat`). LL and RR agree
   between the codes. LR and RL differed by +-0.850 (point 1) and +-0.279
   (point 2) in units of alpha_s/(2 pi), and the magnitudes were those of
   the other helicity order.
2. **Independent calculation from unitarity**
   (`~/cernbox/disorder-comparisons/vbf_nnlojet/unitarity_todd.py`, log
   `unitarity_todd.log`, copy in `tools/`). The only timelike channel of the
   one-loop line amplitude q + J -> Q + g is s_Qg. Its absorptive part is
   (1/2) sum int dPhi_2 A(q J -> Q' g') A(Q' g' -> Q g) from tree amplitudes
   only. These are written from scratch in Python: explicit Weyl spinors and
   SU(3) matrices, and the other line as the current J = ubar gamma u.
   Checks:
   - both trees vanish with eps -> k (1e-12);
   - q g -> q g reproduces the textbook |M|^2 to 1e-14.

   T = 8 pi^2 * 2 Re(A0* i Abs)/|A0|^2 is independent of the t-channel
   cutoff (1e-5 against 1e-7). The collinear log is proportional to A0 and
   drops out. The cut integrand also has an integrable point singularity
   where the cut gluon is collinear to the incoming quark. A tensor Gauss
   rule converged badly there (-0.433, -0.435, -0.429). Splitting the sphere
   into two charts, with a partition of unity, gives -0.411094, stable to
   1e-6. **Result:** for both points, both gluon placements and all four
   helicities, the unitarity T-odd term equals VBFNLO's [V(h) - V(-h)]/2 to
   five digits:
   - point 1, gluon on line 1: LL -0.41109, LR -0.41072;
   - point 1, gluon on line 2: LL +0.44548, LR -0.43355;
   - point 2, gluon on line 1: LL -0.09527, LR -0.09505;
   - point 2, gluon on line 2: LL +0.15015, LR -0.14864.

   VBFNLO is right.
3. **Why NNLOJET seemed to differ: my harness.** In NNLOJET's
   `C1g1VBF(i1,i5,i3,i2,i4,...)`, the gluon is on the i1-i4 line. The
   couplings are `coupling_vbf(nfC1, nfC2, ...)` with glr(2) =
   enL(nfC1) enR(nfC2), labelled "LR". The amplitude in helicity slot 2 is
   slot 1 (++++, RR) with i2 and i3 exchanged, i.e. with the non-radiating
   line flipped. The two agree only if nfC1 is the flavour of the i2-i3 line
   (the one that does not radiate) and nfC2 that of the radiating line. The
   same holds for the second placement (`coupling_vbf(nfC2, nfC1)`). The
   harness had set nfC1 to the radiating line's flavour. That swaps LR and
   RL, which changes only the T-odd term (tree level and the P-even parts
   are symmetric under the swap). It also explains why only NC with mixed
   line types showed a difference.

   With `NFSWAP=1` (nfC1 = the non-radiating line's flavour; logs
   `virt13_swap.log`, and `virt13_noswap.log` for the old order):
   NNLOJET - VBFNLO - 6.853892 is at most 2.2e-12 for all eight processes
   (q- and g-initiated, same and mixed line types, reflected points). For
   the isolated helicities it is at most 7e-10.

   The public NNLOJET release has no VBF driver, so how NNLOJET itself sets
   nfC1 cannot be checked. There is no evidence of an NNLOJET problem.

**Corrects** my earlier statements, the same evening, that "the codes differ
in a T-odd term" and that NNLOJET's sign "follows the other line's
helicity". Both were artefacts of the harness.

## The public POWHEG-BOX-V2 VBF_HJJJ (svn r4135, last changed r3728, 2020-04-23)

Exported to `~/work/disorder-comparisons/powheg_vbf_hjjj_public` for the
report to its authors. Status of the issues there:

1. **Pair-type swap: present.** `compreal_hqqqq.f:529` has the same
   `kl = k+4*ftype(7)-4`, with the same `ftype(7) = 2-mod(abs(bflav(7)),2)`
   (2 = up) and the same `NCmatrix_r` list (1..4 u-type pair, 5..8 d-type).
   It calls the same VBFNLO routine `qqh4q`, without proVBFH's tag
   filtering, so it includes the graphs with the Z on the pair.
2. **Missing initial-state FKS region: present.** Built with its own
   Makefile, then POWHEG's `genflavreglist` was run on its flavour list
   (`tools/list_regions_public.f`, log `list_regions_public.log`; its tags:
   lines 1, 2, the four-quark pair 5; 192 Borns, 856 reals, 1696 regions).
   - All 256 NC entries q(1) -> q(1) + pair(5) have one region, emitter 6
     radiating 7 (g -> q qbar), and no initial-state region.
   - All 128 CC entries with the pair tag on the incoming quark have theirs.

   In that code the singularity is cut off by the Born generation cut
   (`ptcut` on the Born partons), so the NC four-quark channel depends
   logarithmically on it.
3. **nf in `ffunc`: not a problem there** (nf = 4 in both places, see the
   correction above).
4. **T-odd term: no problem** (VBFNLO confirmed by unitarity).

### Direct tests on the public build (`~/work/disorder-comparisons/powheg_vbf_hjjj_public/tools`, copies in `tools/`)

- **Problem 1, initial-state collinear limit** (`is_limit_public.f`, logs
  `is_limit_public*.log`). Uses the public `init_phys`, `setreal` and
  `setborn`. The real is s c -> H s c Q Qbar with the outgoing s || beam,
  x = 0.91 and k_T = 1 ... 1e-3 GeV, from the exact CS initial-initial map.
  The test quantity is c = <R>_phi x 2 p1.k4/(16 pi^2 P_gq(x) B_g), which
  must go to 1. As distributed: c = 1.281755 (Q = u) and 0.780179 (Q = d),
  reciprocals, and <R>(u) equals the fixed <R>(d) to all digits. With
  `kl = k+4*(2-ftype(7))`: c = 1.000000 for both. <R> grows like 1/k_T^2,
  so the unsubtracted singularity of problem 2 is present in the public
  matrix element.
- **One-loop H+3j of the public build** (`virt_heli_public.f`). At the
  helicity-isolated points, against NNLOJET with the correct flavour order:
  - with the public code's own Z/W widths (it ignores ZWIDTH/WWIDTH in
    vbfnlo.input and computes 2.5051/2.0950 GeV; proVBFH reads
    2.4952/2.141), agreement to 4e-6, about 3e-7 relative;
  - with the harness widths imposed, 2e-10, within the helicity-isolation
    precision (7e-10).

  The widths matter because V/B is a Born-weighted average over the two
  gluon placements, whose boson momenta differ.

  The public build links `brakets.f` rather than `new-vbfnlo/brakets_new.f`
  (both are in libvbfnlo.a with the same symbols, and `brakets.o` comes
  first). This makes no difference here (tested). The public
  `hjjj_amp_aux.f` is identical to proVBFH's `hjjj_amp_aux_corrected.F`.

## First integration of the O(alpha_s^2) part (2026-09-30, night)

Timing test (`runs/stage2-timing`, 13 TeV, VBF cuts, cutoff 1e-6, 40k
points, `cs_order 2` with the new `cs_only2 1`, which keeps only the (2,0) +
(0,2) weights): 0.7 ms per point on thA371a, with 70% of the points failing
all cuts before any matrix element.

1. **Bug in my integration driver, fixed.** `nlo2_real_kin` keeps the
   four-parton point in module state. `excl_point` set up both lines first
   and only then called `nlo2_real_me` per line. So line 1's real and
   dipoles were evaluated at line 2's point (with line 2's groups and the
   wrong PDFs) and attached to line 1's events. It now sets up each line's
   point again right before its matrix elements, and checks that the
   dipole flags agree. Effect of the fix: sum |w| went from 494 +- 242 to
   134 +- 13 pb. The limit tests (kinematics and matrix elements line by
   line) could not see it, and no earlier result used this path.
2. **Remaining heavy tail: double-unresolved corners.** New diagnostics:
   `cs_dump2` (per-point contributions), `cs_replay` (per-group,
   per-dipole printout of dumped points) and per-path densities of
   `four_weight`.
   - The largest weights come from the real minus dipoles; the virtual
     part is smooth.
   - Worst point: a soft (2% and 0.4% of the beam energy), collinear
     (s_23/Q^2 = 3e-6) pair of final-state partons next to the hard quark.
     The FF dipole of that pair has y = 0.027, which is not small (the
     three partons form a narrow jet). Its map subtracts y/(1-y) k_1 from
     the pair and leaves the merged gluon at 1e-4 of the beam energy
     (mapped z3 = 1.3e-5, 1 - xp3 = 1.2e-4, just above the cutoff).
   - The dipole's Born is therefore far more singular than the real at
     that point: D = 5e5 against R = 1.6e3, in weight units.
   - The multichannel density of the path through that mapped
     configuration dominates, as it should. The weight equals what one
     expects for the log-sampled double-unresolved measure
     (ln^4(1/cutoff), about 4e4).
   - So this is the known behaviour of CS dipoles in double-unresolved
     regions, not a bug. The observable cancels against the Born
     projection only where event, counterevents and Born pass the cuts
     together, so the spikes in sigma(VBF cuts) come from such corners
     near cut boundaries.
   - **Correction of my first guess:** I suspected the spin-correlated
     g -> gg / g -> q qbar term (soft mapped gluon nearly collinear to the
     beam). Replacing it by its azimuthal average (`cs_spinavg`, diagnostic
     only) changes D by 2%, so that was wrong.
   - Possible remedies if the variance is too large: a larger cutoff (if
     the results are cutoff-independent), Nagy's alpha restriction of the
     dipoles (with its analytic terms in I, K + P), or averaging over
     related points.
3. **Cutoff study** running on thserv18 (`runs/stage2-cut`, commit a552b94,
   30 jobs, nice 10): the O(alpha_s^2) part alone at cutoffs 1e-4, 1e-5
   and 1e-6, 10 seeds x 5.2M points each. Comparison with
   `tools/cutoff_compare.py` (seed-scatter errors besides VEGAS's).
4. **The first cutoff study was aborted** (`runs/stage2-cut-aborted`,
   commit a552b94). This corrects item 3, which described it as running.
   - After one iteration of VEGAS adaptation, 5 of the 10 seeds at 1e-6
     and 2 at 1e-5 hit single points with weights up to 1e13, and their
     grids collapsed. The one seed that finished (c6/s3) is garbage:
     sigma(VBF cuts) = -8e-6 pb, with 88% of the points failing all cuts.
   - Reproduced locally (`runs/stage2-spike`, spike dump `cs_spikes.dat`)
     and replayed. At the spike point, k1 is collinear to the incoming
     parton to s_a1/Q^2 = 2.5e-11, far below the cutoff, and k3 is
     extremely soft (1.4e-6 of the beam energy).
   - The FI paths that cover the initial-state collinear region exclude
     the point (their z is below the cutoff). Only the FF path with y ~ 1
     reaches it, with a tiny density (weight 283).
   - The IF dipoles, still active, cancel only half of R: the soft
     spectator is softer than the collinear recoil it absorbs.
   - Cause (mine): the generator's support is "all six maps at least the
     cutoff away from their singular limits", but the real and dipoles
     were evaluated wherever other paths placed a point.
   - **Fix (9310f3e): a consistent technical cut.** The whole four-parton
     point (real and all counterevents) is dropped if any map has FF y, z,
     1-z or FI 1-x, z, 1-z below the cutoff. R - sum D is integrable, so
     this costs O(c) up to logs; the cutoff study tests it.
   - With the cut, the c6/s1 warm-up is stable (136 +- 6, 143 +- 2) with
     no spike. 43% of the four-parton points are dropped: for narrow
     final-state jets, the FI variable 1-x = s_ij/(2 pa.(pi+pj)) acts
     like a cut on s_ij/Q^2.
5. **Cutoff study, second attempt** (`runs/stage2-cut2`, commit 9310f3e,
   thserv18, 30 jobs, nice 10, started 2026-09-30 03:06): same set-up as
   item 3.

## Cutoff study with the technical cut (2026-09-30)

Report page: https://claude.ai/artifact/K9dmJhRDs5Sd97fb49CevZ (source `report.html`; tables `compare.txt`, `cost.txt`).

Set-up: `runs/stage2-cut2`, commit 9310f3e, thserv18, nice 10. 13 TeV,
VBF cuts, NNPDF30_nnlo_as_0118, mu = mu0(pt,H) of 1506.02660 (`runningscales 1`, i.e. scale_choice 3; corrected, I first wrote mu = Q_i per line). O(alpha_s^2) (2,0) +
(0,2) exclusive part alone (`cs_only2 1`). Cutoffs 1e-4, 1e-5, 1e-6, with
10 seeds x 5.2M points each (warm-up 2 x 200k, production 3 x 1.6M).

Running:
- All 30 jobs finished, with no spikes (no point above 10 pb in
  `cs_spikes.dat`).
- 8881 s CPU per job (1.7 ms per point on thserv18; 0.9 ms on thA371a).
- The technical cut drops 4.6-4.8M four-parton points per job (of
  10.4M line evaluations).

| cutoff | sigma(VBF cuts) [pb] | error, seed scatter | error, VEGAS |
|--------|----------------------|---------------------|--------------|
| 1e-4   | -0.02253             | 0.00219             | 0.00203      |
| 1e-5   | -0.02284             | 0.00228             | 0.00255      |
| 1e-6   | -0.01649             | 0.00632             | 0.00559      |

- **No cutoff dependence.** sigma agrees within 1 sigma. Over all
  histograms (226 bins), chi2/n with seed-scatter errors is 267, 262 and
  254/226 for 1e-4 vs 1e-5, 1e-4 vs 1e-6 and 1e-5 vs 1e-6. With variances
  from 10 seeds, about 1.1-1.3 is expected. The VEGAS errors agree with
  the seed scatter.
- **Size:** -2.6% of sigma_NLO(VBF cuts) = 0.878 pb (stage 1). In the VBF
  distributions it is typically -1% (pt,H) to -6% (pt,j2) of NLO per bin.
- **Cost** (`cost.txt`; median bin with at least 2% of the peak, error
  relative to the NLO bin):
  - At 1e-4, 1% per bin needs 14 (pt,j2), 16 (pt,H), 20 (phi_jj), 37
    (y_j1), 41 (pt,j1) and 51 (M_jj) CPU-h on thserv18. sigma(VBF cuts)
    reaches 0.25% of NLO with 25 CPU-h.
  - At 1e-6 the same costs 3-6 times more, because the variance grows
    with the logs of the cutoff.
  - So the (2,0) + (0,2) part fits comfortably on the thservs, far from
    the "10k one-week runs" of the old code. The (1,1) part and scale
    variations come on top. Tails (high pt, M_jj) will need more.
- **Default cutoff for stage 2: 1e-4** (cutoff-independent at this
  precision and cheapest). Recheck when the precision is pushed further.
