# Reference results of the old proVBFH (for the comparison plots)

- `p1506/`: the 1506.02660 results (old proVBFH, AK's paper files of 21 Feb
  2018, combined with the old trimming): `11.top` μ_R = μ_F = μ0(p_T,H),
  `HH.top` μ0/2, `22.top` 2μ0. NNLO, VBF cuts, histogram names with the
  suffix `-vbf` (strip it to match the p1506 analysis: `combine_parts.py
  --strip -vbf`). From thA371a `proVBFH-cs/runs/ref-1506.02660/`.
- `hxswg136/`: the HXSWG 13.6 TeV study (old proVBFH v2.1.0), from
  `LHCHXSWG/vbf-higgs-wg/proVBFH/results/`: `{lo,nlo,nnlo}-central.top` and
  the scale bands `-min.top`, `-max.top` (see `README-study.txt`; values are
  per bin width). The raw per-seed NNLO data (`13.6TeV_NNLO/{HH,11,22}.tgz`,
  150 MB each) are not included.
- `hxswg136/plain/`: the direct (untrimmed) merge of the study's NNLO
  seeds per scale, and the central/min/max band built the same way as the
  study's (see its README), for estimating the bias of the trimming.
- `hxswg136/symtrim/`: a third reference, the symmetric 0.5% trimmed mean
  per bin with bootstrap errors (see its README), same band construction.
