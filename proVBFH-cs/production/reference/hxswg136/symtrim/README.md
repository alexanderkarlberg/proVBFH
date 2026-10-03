# HXSWG 13.6 TeV study: robust (symmetric-trim) merge of the NNLO seeds

A third reference next to the study's trimmed files (`../nnlo-*.top`) and
the plain merge (`../plain/`). Made on 3 Oct 2026 from the same raw
per-seed files with `notes/2026-09-30-hxswg-comparison/tools/merge_plain.py`.

Method, per bin and per scale: sort the seeds, drop the lowest and highest
0.5% (49 of ~9950 on each side), average the rest. The error is the
standard deviation of this estimator over 400 bootstrap resamplings of the
seeds. Compared with the study's `combine_runs`:
- the cut is a fixed fraction on both sides, the same in every bin (the
  study's cut, median ± 5 × (q84 − q16), removes a bin-dependent and
  usually one-sided number of seeds, and has an off-by-one);
- the error is that of the estimator actually used (the study combines the
  per-seed VEGAS errors of the kept seeds, which underestimates it).
The trimmed mean is biased low by an amount set by the upper tail (the
seed distributions have power-law tails, Hill index 1.0–1.6); the plain
mean is unbiased but its error is unreliable for such tails. The two
bracket the true value.

Checks (central scale, all 946 non-empty bins, scratch scripts):
- trim fraction 0.25%, 1%, 2% instead of 0.5%: median shift 0.27, 0.28,
  0.49 bootstrap errors;
- winsorising at 0.5% instead of trimming: same values, errors 1.26× larger;
- median-of-means agrees in the total cross sections, but drifts with the
  number of groups in bins with skewed tails, so not used;
- Peng's tail-corrected mean is unstable for Hill index ≈ 1, not used.

Files: `nnlo-{HH,11,22}-symtrim.top` per scale; `nnlo-{central,min,max}-symtrim.top`
with central = 11 and min/max the per-bin envelope of HH, 11, 22 (the
study's definition). Format: xlow xhigh value error fraction (values per
bin width; fraction of seeds used).

| σ(ptj > 20) [pb] | symmetric trim | study (trimmed) | plain |
|---|---|---|---|
| μ0/2 (HH) | 2.07465 ± 0.00111 | 2.07425 ± 0.00095 | 2.08247 ± 0.00847 |
| μ0 (11) | 2.09022 ± 0.00088 | 2.08980 ± 0.00075 | 2.09614 ± 0.00649 |
| 2μ0 (22) | 2.11129 ± 0.00072 | 2.11109 ± 0.00060 | 2.11588 ± 0.00506 |

All bins (central; min and max are similar):
- symmetric trim − study: median +0.05% (+0.15 symmetric-trim errors);
  more than 3 errors apart in 81 of 946 bins;
- symmetric-trim errors are a median 1.20× the study's;
- plain − symmetric trim: median +0.34 plain errors.
