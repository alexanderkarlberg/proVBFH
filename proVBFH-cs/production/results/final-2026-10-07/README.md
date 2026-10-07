# Final combined results of the proVBFH-cs production (7 Oct 2026)

Errors come from the seed scatter of the finished jobs: plain means, no trimming. Weights
W1 = (1,1), W2 = (1/2,1/2), W3 = (2,2) of μ0(p_T,H). The exclusive and inclusive parts are
added. The full log is in `notes/2026-10-cs-production/README.md`. The raw per-job outputs
stay in `/ptmp/mpp/akarlber/cs-production`.

| set-up | order | exclusive jobs | inclusive jobs |
|---|---|---|---|
| p1506 | NNLO | 11,000 | 1,000 |
| p1506 | NLO | 1,650 | 150 |
| p1506 | LO | - | 150 |
| hxswg136 | NNLO | 10,872 of 11,000 (*) | 1,000 |
| hxswg136 | NLO | 6,540 of 6,600 | 150 |
| hxswg136 | LO | - | 150 |

(*) The last 128 jobs were killed by a cluster breakdown on 7 Oct, and the remaining ones
by node failures. AK decided to leave them out (no impact).

Some tail bins of the NNLO distributions are dominated by single events, i.e. one job out of
11,000 carries the bin's value and error (`tools/spikescan.py`; see the notes, 7 Oct).
These are under investigation and are kept in the plain mean.

## `h3j-crosscheck/`

≥ 3- and 4-jet histograms of the 1506.02660 analysis, in pb, at O(α_s²) for ≥ 3 jets unless
marked LO:

| file | content |
|---|---|
| `vbfnlo-nlo-3j.top` | VBFNLO 3.0, process 110, NLO H+3j at μ0, 800 seeds (plain mean of the last iteration) |
| `vbfnlo-lo-3j.top` | the same at LO, 50 seeds |
| `vbfnlo-nlo-3j-fixmh.top` | VBFNLO at μ_R = μ_F = m_H, 800 seeds |
| `cs-nnlo-3j-fixmh.top` | proVBFH-cs NNLO exclusive at μ = m_H, 1,998 jobs (≥ 3-jet observables only; the inclusive part does not contribute to them) |
| `powheg-lo.top`, `powheg-lo-nf5.top` | POWHEG-BOX-V2 VBF_HJJJ LO with the running-scale fix: 4 and 5 flavours |
| `powheg-lo-nf5nc.top` | POWHEG LO, like-for-like with VBFNLO and proVBFH-cs: 5 flavours, negative PDFs kept, narrow-width Higgs |
| `powheg-nlo-pt1.top` | POWHEG NLO, 4 flavours, ptcut 1 GeV |
| `powheg-nlo-nf5-pt1.top` | POWHEG NLO, like-for-like, NLO problems 1-3 present, 660 seeds |
| `powheg-nlo-fix123.top` | POWHEG NLO, like-for-like with problems 1-3 fixed and `st_nlight 5`, 700 seeds |

σ(≥ 3 jets) in fb:

| | μ = m_H | μ0 |
|---|---|---|
| proVBFH-cs | 126.19 ± 0.18 | 126.87 ± 0.17 |
| VBFNLO 3.0 | 125.22 ± 0.12 | 125.69 ± 0.12 |
| POWHEG, fixes 1-3 | - | 127.4 ± 2.9 |
| POWHEG, as distributed (like-for-like) | - | 133.4 ± 3.0 |
| old proVBFH (paper) | - | 133.24 |

proVBFH-cs lies 0.8-0.9% (about 1 fb, 4.5-6σ) above VBFNLO at both scales; the 4-jet rates
agree. The cause is open (see the notes).
