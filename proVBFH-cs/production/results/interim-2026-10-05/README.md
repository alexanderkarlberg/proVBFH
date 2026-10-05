# Interim combined results (5 Oct 2026, 08:00)

Not final: the production is still running (see
`notes/2026-10-cs-production/README.md`). Errors from the seed scatter of
the finished jobs, no trimming. Weights W1 = (1,1), W2 = (1/2,1/2),
W3 = (2,2) of μ0(p_T,H). Exclusive + inclusive parts added.

| set-up | order | exclusive jobs | inclusive jobs |
|---|---|---|---|
| p1506 | NNLO | 1,386 | 1,000 |
| p1506 | NLO | 1,650 | 150 |
| p1506 | LO | - | 150 |
| hxswg136 | NNLO | **pilot only**: 198 | 20 |
| hxswg136 | NLO | 3,779 | 150 |
| hxswg136 | LO | - | 150 |

`h3j-crosscheck/`: ≥ 3-/4-jet histograms of the 1506.02660 analysis at
μ0 (pb): `vbfnlo-lo-3j` (50 seeds), `vbfnlo-nlo-3j` (800 seeds),
`powheg-lo` (4 flavours, scale fix), `powheg-lo-nf5` (5 flavours; the
like-for-like LO with 5 flavours, negative PDFs kept and a narrow-width
Higgs gives 131.29 ± 0.90 fb, see the notes),
`powheg-nlo-pt1`, `powheg-nlo-pt01` (4 flavours, ptcut 1 / 0.1 GeV, 300
seeds each). `powheg-lo-lfl`, `powheg-nlo-lfl-pt1`: like-for-like (scale fix, 5 flavours, negative PDFs kept, narrow-width Higgs), 100 and 660 seeds.
