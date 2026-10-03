# HXSWG 13.6 TeV study: direct (untrimmed) merge of the NNLO seeds

For estimating the bias of the study's trimming (`combine_runs.f`, which
drops seeds outside median ± 5 × (84%−16% range) per bin). Made on 3 Oct
2026 from the study's raw per-seed files
(`LHCHXSWG/vbf-higgs-wg/proVBFH/13.6TeV_NNLO/{HH,11,22}.tgz`) with
`notes/2026-09-30-hxswg-comparison/tools/merge_plain.py`:

- `nnlo-{HH,11,22}-plain.top`: per scale (μ0/2, μ0, 2μ0) the plain mean
  over all seeds (9973, 9940, 9965; one empty HH file skipped), error =
  seed scatter/√N. Note: the seed distributions have power-law tails
  (Hill index 1.0–1.6, see the notes), so these errors are themselves
  uncertain and dominated by a few seeds.
- `nnlo-{central,min,max}-plain.top`: central = 11, min/max = the per-bin
  envelope of HH, 11, 22 (the errors those of the extremal scale). This is
  the study's definition: the same construction with the trimmed
  combinations reproduces the study's `nnlo-{central,min,max}.top` to
  5e-8 in all 946 non-empty bins.

Format as the study's files: xlow xhigh value error fraction (values per
bin width; fraction = 1).

| σ(ptj > 20) [pb] | plain | trimmed (study) |
|---|---|---|
| μ0/2 (HH) | 2.08247 ± 0.00847 | 2.07425 ± 0.00095 |
| μ0 (11) | 2.09614 ± 0.00649 | 2.08980 ± 0.00075 |
| 2μ0 (22) | 2.11588 ± 0.00506 | 2.11109 ± 0.00060 |

Central, all bins: plain − trimmed has median +0.52% (+0.43 plain errors),
positive in 59% of the 946 bins.
