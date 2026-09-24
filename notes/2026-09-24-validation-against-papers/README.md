# Validation against published results (2026-09-24)

Before setting up CI, we checked that the current code, built against
the latest dependency releases, reproduces the numbers in the papers
based on it:

- [1506.02660] Cacciari, Dreyer, Karlberg, Salam, Zanderighi: VBF H at NNLO, fully differential
- [1606.00840] Dreyer, Karlberg: VBF H at N3LO (inclusive)
- [1811.07906] Dreyer, Karlberg: VBF HH at N3LO (inclusive)
- [1811.07918] Dreyer, Karlberg: VBF HH at NNLO, fully differential

PDFs of these papers, plus [2005.11334] (Dreyer, Karlberg, Tancredi:
non-factorisable corrections), are in `papers/`.

Versions: hoppet 2.3.0, LHAPDF 6.5.6, FastJet 3.5.1, gfortran 15.2 (Linux).
Code at commit b7da4bf plus an rpath fix in `proVBFH-inclusive/Makefile`
(no physics change) and the coupling fix in `proVBFHH` described below.

## Inclusive (proVBFH-inclusive)

Common flags: `-ncall1 1000000 -ncall2 4000000 -itmx1 3 -itmx2 3 -7scaleuncert`.
Each order is a separate run (`-lo`, `-nlo`, `-nnlo`, `-n3lo`), because
only the highest order requested is written out.

**H, 13 TeV, 1606.00840 table 1** (μ = Q_i):

```
provbfh_incl -n3lo -sqrts 13000 -pdf PDF4LHC15_nnlo_mc -mh 125.09 -mw 80.385 -wwidth 2.085 -scale-choice 1
```

| order | this run [pb]          | paper [pb]                 |
|-------|------------------------|----------------------------|
| LO    | 4.0962 +0.0543 −0.0643 | 4.099 +0.051 −0.067        |
| NLO   | 3.9677 +0.0274 −0.0201 | 3.970 +0.025 −0.023        |
| NNLO  | 3.9294 +0.0178 −0.0079 | 3.932 +0.015 −0.010        |
| N3LO  | 3.9270 +0.0071 −0.0000 | 3.928 +0.005 −0.001        |

MC error ±0.0013 pb on each.

**HH, 14 TeV, 1811.07906 table 1** (μ = Q_i):

```
provbfhh_incl -n3lo -sqrts 14000 -pdf PDF4LHC15_nnlo_mc -mh 125 -mw 80.379 -wwidth 2.085 -hwidth 4.030e-3 -scale-choice 1
```

| order | this run [fb]          | paper [fb]                 |
|-------|------------------------|----------------------------|
| LO    | 2.0749 +0.1800 −0.1489 | 2.079 +0.177 −0.152        |
| NLO   | 2.0616 +0.0216 −0.0147 | 2.065 +0.022 −0.018        |
| NNLO  | 2.0519 +0.0046 −0.0025 | 2.056 +0.003 −0.005        |
| N3LO  | 2.0514 +0.0034 −0.0000 | 2.055 +0.001 −0.001        |

MC error ±0.0011 fb on each. At 27 TeV, N3LO: 8.3945 ± 0.0050 vs 8.407 fb.

**No-cuts column of the differential papers** (μ₀(p_t) scale, `-scale-choice 3`):

| process, order | this run | paper |
|---|---|---|
| H 13 TeV NNPDF30_nnlo_as_0118, LO  | 4.0276 pb | 4.032 pb |
| H 13 TeV NNPDF30_nnlo_as_0118, NLO | 3.9250 pb | 3.929 pb |
| HH 14 TeV PDF4LHC15_nnlo_mc, LO    | 2.0135 fb | 2.016 fb |
| HH 14 TeV PDF4LHC15_nnlo_mc, NLO   | 2.0470 fb | 2.049 fb |

(H: `-mh 125 -mw 80.398 -wwidth 2.141`; HH: `-mh 125 -mw 80.379 -wwidth 2.141 -hwidth 4.030e-3`.)

**Open point.** All inclusive numbers come out slightly low: by
0.03–0.1% for H and about 0.18% for HH (roughly 3σ of MC error for HH).
The offset is the same at every order, including LO, so it comes from
the inputs or PDF handling rather than the QCD corrections. Halving
the hoppet grid spacing (`-dy 0.025`) changes the HH LO result by only
2e-6 relative, so the y-grid is not the cause.

The papers predate bfdb73e (2025-01-30), which raised the phase-space
cut on min(Q1, Q2) from `Qmin = 1 GeV` to `max(sqrt(2), sqrt(Q2minPDF))`
to protect downward μ_F variations. Rerunning LO with `Qmin = 1.0_dp`
(same seed and statistics, so the MC fluctuations are largely
correlated) gives:

| LO setup                  | Qmin = √2 (current) | Qmin = 1 GeV | shift  | paper |
|---------------------------|---------------------|--------------|--------|-------|
| H 13 TeV, μ = Q           | 4.09624             | 4.09776      | +0.037% | 4.099 |
| H 13 TeV, μ₀(p_{t,H})     | 4.02760             | 4.02920      | +0.040% | 4.032 |
| HH 14 TeV, μ = Q          | 2.07489             | 2.07594      | +0.051% | 2.079 |

So the Qmin change accounts for about half of the H offset; with it
reverted, H at μ = Q agrees with the paper to about 1σ.

For HH, two more changes since the paper contribute (LO, 14 TeV, same
seed, applied cumulatively on top of Qmin = 1):

| change relative to the paper-era code                          | shift   |
|----------------------------------------------------------------|---------|
| Qmin √2 → 1 GeV (bfdb73e)                                       | +0.051% |
| EW masses in hoppet's `StartStrFct` (sin²θ_W in the Z couplings of the structure functions) | +0.038% |
| old tensor matrix element (`-tensorME`, before 03adc62/edc7a24) instead of the non-tensor default | +0.023% |

On the EW masses: until the hoppet v2 port the code did not pass
`mw, mz` to `StartStrFct`, so hoppet used its built-in M_W = 80.398,
M_Z = 91.187 (sin²θ_W = 0.22265) for the Z couplings of the structure
functions, while the matrix element used M_W = 80.379. The code now
passes the same M_W to both (sin²θ_W = 0.22301). That makes it
consistent, so this is a deliberate change, not a regression.

On the matrix element: 2.0.0 always called
`eval_matrix_element_tensor`. Since 7d5e5fa the faster non-tensor
(analytic) `eval_matrix_element` is the default. The two differed by
about 2e-4 at LO (the +0.023% above).

**Correction (later the same day).** A first analysis of this, with a
debug build evaluating both implementations on the same phase-space
points, took the tensor version as the reference and concluded that
the F3 × F3 interference terms of the analytic `F3F3` expression in
`ME_expressions.f` were wrong. That was backwards: the tensor version
had the bug (fixed in 03adc62). It stored the contravariant components
of q, P and the Higgs momenta in tensors declared with lower indices,
i.e. it worked with the parity-flipped vectors, but built
iε_{μνρσ}P^ρq^σ from the unflipped ones. Terms without ε are parity
invariant and came out right, as did F3 × F3 for the g^{μν} part of the
current (A) alone; the F3 × F3 interference of A with the t/u-channel
parts (B, C = (2k+q)^μ(k'−k−q)^ν) did not.

An independent calculation with explicit index handling at a physical
phase-space point (`scripts/f3-check/`) settles it: the analytic
`F3F3` coefficients AA, AB, AC and BC agree with it to all printed
digits (BB and CC vanish), while the old tensor convention reproduces
AA only (AB is off by a factor 39, AC by 4.4, BC by 22). The F1 × F1
coefficients agree in both.

After the fix, and with the order combination below, same-point
comparisons in provbfhh_incl (HH, 14 TeV, PDF4LHC21_40, 0.6M points,
debug build from `scripts/dbg_patch.py`) give:

| comparison                                            | analytic/tensor − 1 (integral) | pointwise |
|-------------------------------------------------------|---------------------|-----------|
| zero widths, LO                                       | −2e-16 ± 9e-16      | ≤ 1.4e-5; points above 1e-8 carry 4e-9 of σ |
| zero widths, N3LO                                     | +4e-16 ± 9e-16      | ≤ 5.9e-6 |
| physical widths, LO                                   | (−2.2 ± 1.2)e-4     | up to 46% of the scale |
| physical widths, tensor averaged over the point and its mirror image (y → −y), LO | +7e-16 | ≤ 1.4e-9 of the scale |
| same, N3LO                                            | −1e-15              | ≤ 3.9e-9 of the scale |
| without F3, physical widths, LO / N3LO                | −4e-14 / −1e-13     | 92% of points within 1e-12 (99.7% of σ) |

(The pointwise scale is the sum of the absolute values of all terms.
The few points with larger relative deviations have tiny matrix
elements at extreme kinematics, where both implementations lose
digits; old and new tensor code differ there by up to 1e-9 as well,
although they are algebraically identical at LO.)

So with finite widths the only difference left is the F1·F3 and F2·F3
terms, Im(A B*) etc. times the iε part of the hadronic tensor, which
the analytic version does not have. Each contains a single ε tensor, so
they are parity odd: large pointwise, but they integrate to zero for
the cross section and for any parity-even distribution. The analytic
version is correct for all parity-even observables.

The old tensor code (the 2.0.0 default, used for the 2018 HH papers and
always used for the NF corrections) had wrong F3 interference terms.
Their integrated effect is at the 1e-4 level: analytic/old tensor − 1 =
(−1.5 ± 0.9)e-4 in the first analysis (6.5M points, physical widths),
(−1.9 ± 3.2)e-4 at zero width here.

Decision (user, 2026-09-24): the analytic matrix element stays the
default; the tensor one (`-tensorME` / `tensorME 1`) is needed only for
parity-odd observables, the diagram switches and the non-factorisable
corrections.

Also found on the way (af68ff5): proVBFHH's matrix elements used
`cmplx()` without `kind=dp`, i.e. single-precision propagators (a
1.5e-8 effect on the LO cross section).

**Order combination.** The two versions also combined perturbative
orders differently, in proVBFH-inclusive and in proVBFHH alike. The
non-tensor version expands strictly, summing Fx1(i)·Fx2(j) with
i + j = n + 1. The old tensor version summed the structure functions
over all orders up to `order_stop` first and contracted the product, so
it included products beyond the working order (for example NLO×NLO at
NLO), and it ignored `order_start`. In proVBFHH this mattered only with
`tensorME 1` beyond LO: `order_min` is hard-wired to 1, the non-tensor
version is the default, and the non-factorisable path, which always
uses the tensor version, is LO-only. Fixed in edc7a24 (see below).

Both (old) versions evaluated on the same 6.5M points (HH, 14 TeV;
`-ncall1 1000000 -ncall2 2000000`, 3+3 iterations). Removing F3 isolates
the order-combination effect:

| order | non-tensor/old tensor − 1, without F3 (order combination) | extra terms in the old tensor version | F3 share of σ (old tensor / non-tensor) | non-tensor/old tensor − 1, total (±0.9e-4) |
|-------|---------|--------------------------|-----------------|---------|
| LO    | 2e-15   | none                     | 9.1e-4 / 7.6e-4 | −1.5e-4 |
| NLO   | +2.7e-4 | NLO×NLO                  | 8.5e-4 / 7.0e-4 | +1.2e-4 |
| NNLO  | +0.65e-4| products up to NNLO×NNLO | 8.5e-4 / 7.0e-4 | −0.9e-4 |
| N3LO  | +0.06e-4| products up to N3LO×N3LO | 8.5e-4 / 7.0e-4 | −1.5e-4 |

The order-combination terms fall off by about a factor 4–10 per order,
as expected for terms beyond the working accuracy. The
order-independent remainder of about −1.5e-4 in the last column is the
old tensor code's F3 error. The 2018 papers used the old tensor
version, so their HH inclusive numbers include both (at the 1e-4
level). The tensor version now contracts per order pair: since the
hadronic tensor is linear in the structure functions,
W_b = Σ_k F_k G_b(k), the contractions T(k,l) = Re Tr[(G_1(k) M)(M* G_2(l))]
of the three basis tensors are computed once per point and boson, and
the orders are summed with the structure functions. The debug build
reproduces the table above with the new code: old/new tensor − 1
without F3 = −2.8e-4 (NLO), −0.70e-4 (NNLO), −0.057e-4 (N3LO).

Separate same-seed runs of `provbfhh_incl` with the old code showed the
same pattern within their ±1e-4 resolution: tensor/non-tensor − 1 =
+2.3e-4 (LO), −0.2e-4 (NLO), +1.7e-4 (NNLO), +2.2e-4 (N3LO). With the
new code and zero widths, same-seed runs with and without `-tensorME`
give identical results and MC errors at every order (LO
2.081409 ± 0.005302 fb; N3LO 2.059101 ± 0.005539 fb), also with
coupling modifiers (cVVHfact 1.1, cVVHHfact 0.7, lambdafact 2.5: LO
7.393876, NNLO 7.467001 fb), which the provbfhh_incl tensor version
ignored before 44d0d79. proVBFHH `inclusive_only 1`, NNLO, zero
widths: 1.729837 ± 0.008886 fb for both (old tensor code 1.729126).

Where the non-tensor version is used: the proVBFHH-inclusive default
since 7d5e5fa (2025-01), and the projected inclusive part of the full
proVBFHH code (`src/exclusive/incl2pwhg.f:51`, which ignores
`tensorME`) since proVBFHH 1.0.0. Since the analytic `F3F3` is correct,
the 1811.07918 differential results are not affected by the F3 issue
(they lack only the parity-odd width terms).

Timing (provbfhh_incl, 1M points, one core): the old tensor version
took 7.4 s at LO and 9.1 s at N3LO, the new one 4.4 s and 6.1 s, the
analytic one 2.7 s and 4.4 s. (Earlier measurement with the old code:
12.3 s vs 5.7 s at LO for 1.2M points, on a busier machine.) The
remaining cost is the 4×4 complex algebra in `trace_matrix`; the
structure functions dominate the analytic version.

Averaging four further seeds (11–14) to reduce MC error:

| configuration                                    | HH LO 14 TeV [fb]  | vs paper 2.079 |
|--------------------------------------------------|--------------------|----------------|
| current code                                     | 2.07745 ± 0.00054  | −0.075%        |
| Qmin = 1, 2018 hoppet EW defaults, old `-tensorME` | 2.07983 ± 0.00055  | +0.040% (1.5σ) |

(Seed 10, used in the tables above, fluctuates low by about 0.1%.) The
paper-era configuration reproduces 1811.07906, so the HH offset is fully
explained.

## Differential (proVBFH, proVBFHH), with VBF cuts

These start from the `example/` inputs with `runningscales 1` and the
standard VBF cuts. The PDF, beam energy and EW parameters (in
`vbfnlo.input`) are set to each paper's values:

- H (1506.02660): 13 TeV, `lhans 261000` (NNPDF30_nnlo_as_0118), M_W = 80.398, Γ_W = 2.141
- HH (1811.07918): 14 TeV, `lhans 91500` (PDF4LHC15_nnlo_mc), M_W = 80.379, Γ_W = 2.141, Γ_H = 4.030e-3

LO (`qcd_order 1`): 8 seeds with `ncall1 1000000`, `ncall2 10000000`,
merged with `aux/combine_runs pwg-LO-*.top`. NLO (`qcd_order 2`):
`aux/runpar.sh` on 8 cores (3 grid iterations, then stage 2), merged with
`aux/combine_runs pwg-*-NLO.top`. The H NLO run used the example
statistics (`ncall2 5000000`); the HH NLO run used `ncall2 1000000`
because proVBFHH is slower. The number quoted is the
`sig(all VBF cuts 2 jets)` bin of `total_distrib.top`.

| process, order | this run             | paper    |
|----------------|----------------------|----------|
| H LO           | 0.95722 ± 0.00025 pb | 0.957 pb |
| H NLO          | 0.8740 ± 0.0014 pb   | 0.876 pb |
| HH LO          | 0.79879 ± 0.00018 fb | 0.799 fb |
| HH NLO         | 0.7310 ± 0.0024 fb   | 0.726 fb |

H NLO wall time: about 3 min for the grids and 29 min for stage 2 on 8 cores.
HH NLO (ncall2 1M): about 1 min for the grids and 11 min for stage 2.

**Bug found: proVBFHH NLO/NNLO gave about 1e20 fb.** Since 828ff81
(2025-03-24), `src/exclusive/convert_coup.f` has read `cVVHHfact`,
`cVVHfact` and `lambdafact` with `powheginput` and no default. When the
keys are missing from `powheg.input` (as in `example/`), the value is
−1e6, which scales the HHjj/HHjjj couplings on the exclusive side.
Inclusive-only and LO runs were unaffected, because
`incl_parameters.f90` already defaults these to 1. With the fix (default
1 unless set, as on the inclusive side), all 8 seeds agree and the
result is the number in the table above.

## Non-factorisable corrections (2005.11334)

Dreyer, Karlberg, Tancredi: non-factorisable (NF) NNLO corrections in
the eikonal approximation, for VBF H and HH at 13 TeV with VBF cuts.
Setup as in 1506.02660: NNPDF30_nnlo_as_0118, M_W = 80.398,
Γ_W = 2.141, Γ_Z = 2.4952, M_H = 125, Γ_H = 4.030e-3 (HH), central
scale μ₀(p_t) (`runningscales 1`). NF runs use `qcd_order 1`,
`nonfact 1`. For HH, `tri_off` / `box_t_off` + `box_u_off` select the
diagram classes (TT = triangle topologies, BB = boxes, TB = their
interference = all − TT − BB), and `tri1_off`/`box1_off` or
`tri2_off`/`box2_off` keep only the 2-loop or 1-loop pieces. The gluon
mass is hard-wired to λ = M_V (`matrix_element.f90:959,1045`), as in the
paper's table. The diagram switches exist only in the tensor code, so
the HH Born breakdown uses `tensorME 1`. The NF path always uses the
tensor code.

Statistics: H NF 6.6M calls; HH NF 2.4M calls per run; HH Born 16.5M.
Scripts: `scripts/nf-setup.sh`, `scripts/nf-runs.sh`; results in
`runs/nonfact/` (`sig(all VBF cuts 2 jets)` bin of `pwg-LO-0001.top`).
Only central scales were checked; scale variations were not.

**Single Higgs:** δσ(NF) = −0.003025 ± 0.000010 pb, vs −0.0030 pb in
the paper's table 3. That is −0.316% of LO (0.95722 pb), vs −0.32%
quoted in the text.

**HH, the paper's table 2:** our number first, the paper's in brackets.

|              | σ_TT              | σ_BB              | σ_TB                | Σ                |
|--------------|-------------------|-------------------|---------------------|------------------|
| Born [fb]    | 10.368 (10.393)   | 14.128 (14.172)   | −23.839 (−23.904)   | 0.6577 (0.662)   |
| 1-loop NF    | 0.342% (0.339%)   | 0.303% (0.300%)   | 0.321% (0.318%)     | 0.274% (0.286%)  |
| 2-loop NF    | −0.666% (−0.667%) | −0.620% (−0.621%) | −0.643% (−0.644%)   | −0.513% (−0.516%)|
| full NF      | −0.324% (−0.327%) | −0.318% (−0.320%) | −0.323% (−0.326%)   | −0.239% (−0.230%)|

The full-NF row is the sum of the 1-loop and 2-loop runs. A direct
full-NF run with all diagrams gives −0.00152 ± 0.00004 fb = −0.231%
(paper: −0.0015 fb, −0.230%). MC errors on the NF pieces are 0.2–0.5%
for TT and BB, 0.6% for the 2-loop Σ and 3.4% for the 1-loop Σ, which
has a strong cancellation. All NF percentages agree with the paper
within about 1% relative, or within the MC error.

**HH after the tensor fixes (03adc62, af68ff5, edc7a24).** The NF path
and the Born breakdown always use the tensor matrix element, so the HH
runs were repeated with the fixed code (same setup, statistics and
seeds; hoppet 2.3.0, LHAPDF 6.5.6; `scripts/nf-table.py`, runs in
`runs/nonfact-after-tensor-fix/`):

|              | σ_TT              | σ_BB              | σ_TB                | Σ                |
|--------------|-------------------|-------------------|---------------------|------------------|
| Born [fb]    | 10.368 (10.393)   | 14.128 (14.172)   | −23.838 (−23.904)   | 0.6577 (0.662)   |
| 1-loop NF    | 0.342% (0.339%)   | 0.303% (0.300%)   | 0.320% (0.318%)     | 0.283% (0.286%)  |
| 2-loop NF    | −0.666% (−0.667%) | −0.620% (−0.621%) | −0.643% (−0.644%)   | −0.516% (−0.516%)|
| full NF      | −0.324% (−0.327%) | −0.318% (−0.320%) | −0.323% (−0.326%)   | −0.233% (−0.230%)|

Direct full-NF run: −0.00153 ± 0.00005 fb = −0.233% (paper −0.0015 fb,
−0.230%). The Born TT result is identical to the last digit (the
triangle diagrams have only the g^{μν} structure, so the F3 bug, an
A × B interference, cannot affect them); BB and Σ change by 1–2e-5.
TT, BB and TB NF percentages are unchanged within 0.001; the Σ column,
which has the strongest cancellation, moves towards the paper
(1-loop 0.274% → 0.283%, 2-loop −0.513% → −0.516%, full −0.239% →
−0.233%), within its MC error (±0.010, ±0.004, ±0.011) in both cases.
The paper was computed with the old tensor code. H NF is unaffected
(proVBFH uses the analytic matrix element).

**HH Born normalisation (closed, accepted as is).** The HH Born numbers are
0.25–0.3% below the paper (Σ −0.65%), well outside the 0.05% MC error.
The paper does not state the PDF used for HH; it only says the EW
parameters are "set identically" to the single-Higgs study. With
PDF4LHC15_nnlo_mc (used in the earlier HH papers) instead, the results
come out high by about the same amount: TT 10.422 (+0.28%), BB 14.213
(+0.29%), TB −23.972 (+0.28%), Σ 0.6630 (+0.15%). Neither PDF matches
exactly. This does not affect the NF validation, since the NF
percentages are ratios to the Born. For comparison, the same code
reproduces the 1811.07918 HH Born (14 TeV, PDF4LHC15, M_W = 80.379) to
0.03%. Reviewed with the user on 2026-09-24 and accepted as is for now,
since the NF validation does not depend on it.

## Practical notes

- At NLO, the `sig incl cuts` bin is not the no-cuts cross section,
  because `phspcuts 1` applies generation cuts. Take no-cuts numbers
  from the inclusive programs.
- Also at NLO, the per-seed `pwg-st2-*-stat.dat` totals are not
  physical cross sections. Use the combined histograms.
- If `LHAPDF_DATA_PATH` is set in the environment, it overrides the
  data directory of the LHAPDF installation the program is linked to.

## Files in this directory

- `papers/`: the papers listed at the top.
- `runs/`: inputs, logs and results of every run quoted above. VEGAS
  and POWHEG grids and the per-point debug dumps (about 110 MB) were
  not kept; they are regenerated by rerunning.
  - `incl/`: inclusive runs, all orders (tables at the top)
  - `proVBFH-LO/`, `proVBFH-NLO/`, `proVBFHH-LO/`, `proVBFHH-NLO/`:
    differential runs; `total_distrib.top` holds the combined result
  - `qmin/`: current code vs Qmin = 1 GeV, vs paper-era structure-function
    EW parameters, each with and without `-tensorME`
  - `seeds/`: seeds 11–14 of the current and paper-era HH configurations
  - `width/`: tensor vs non-tensor with widths set to zero
  - `timing/`: tensor vs non-tensor timing
  - `nonfact-after-tensor-fix/`: the HH NF and Born runs repeated
    with the fixed tensor code
  - `tensor-fix/`: logs of the same-point debug runs of the correction
    above (`runs-dbg4-*`: zero widths, physical widths, no F3;
    `runs-dbg5-*`: mirror-averaged), same-seed tensor vs analytic runs
    at all orders with zero widths (`orders-zero-width/`) and with
    coupling modifiers (`kappa/`)
  - `orders/`, `orders-same/`, `meint/`, `mecmp*/`: tensor vs non-tensor
    comparisons (separate runs and same-point harness)
- `patches/`: the modified builds of proVBFH-inclusive used above, as
  diffs against `proVBFH-inclusive/src` at the `ci` branch commit
  1326b19. Apply with `patch -p0` from a copy of `proVBFH-inclusive/`.
  - `qmin-1GeV.patch`: `Qmin = 1.0_dp`
  - `paper-era-config.patch`: Qmin = 1 GeV plus hoppet's pre-v2
    default M_W = 80.398, M_Z = 91.187 in `StartStrFct` (use with
    `-tensorME` to reproduce the 2018 configuration)
  - `me-debug-harness.patch`: switches `tensor_F_mode`, `fA/fB/fC` and
    `keep1/keep2` in `matrix_element_dihiggs.f90`, and a
    `debug_compare` routine in `phase_space_dihiggs.f`, called when
    `-tensorME` is set. In this final version it accumulates
    VEGAS-weighted integrals of both implementations on the same points
    and writes `me_integrals.dat` every 500k points. Earlier stages of
    the investigation used variants of `debug_compare` that wrote
    per-point dumps (`me_compare.dat`, `me_terms.dat`, `me_F.dat`,
    `me_F3.dat`, analysed with `scripts/terms.py`-style scripts).
    This harness belongs to the first analysis, which took the old
    tensor code as the reference (see the correction above); it applies
    to the code before the tensor fixes.
- `scripts/`: the driver scripts for the differential runs
  (`diff-runs.sh`, `hh-nlo.sh`), the NF runs (`nf-setup.sh`,
  `nf-runs.sh`, `nf-table.py` for table 2) and `terms.py` (per-term analysis of `me_terms.dat`).
  Added with the tensor fixes:
  - `f3-check/`: the independent, index-correct evaluation of the F3
    terms at a physical point, compared with the analytic `F3F3` and
    with the old tensor convention (see its README).
  - `dbg_patch.py`: same-point debug harness for the current code
    (analytic, new tensor and old tensor matrix element on the same
    points; `NOF3`, `MIRROR` switches), used for the tables in the
    correction above. They contain hard-coded scratch paths (`SP=...`) and
  expect a dependency install made with `ci/install-deps.sh`, with
  `LHAPDF_DATA_PATH` pointing at it.
