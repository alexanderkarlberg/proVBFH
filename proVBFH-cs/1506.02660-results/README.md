# 1506.02660 set-up: plots of the proVBFH-cs output

Inputs: `{lo,nlo,nnlo}-{central,min,max}[-nlo-ratio].top` (new proVBFH-cs, cards `powheg.input-*`),
`old/` (copies of `11.top` central, `HH.top` (1/2,1/2), `22.top` (2,2), and their envelopes
`max.top`/`min.top` from the old proVBFH results, `1506.02660_bug_fix_results/`).

The `.top` inputs come from runs with the cards `powheg.input-{lo,nlo,nnlo}` (scripts `run-*.sh`; the
outputs `provbfh-<seed>-<order>-W<k>.top` of the seeds are merged per scale point, W1 central, W2
(1/2,1/2), W3 (2,2)). The cards were used with a version of the code that gave the inclusive part
`ncall/10` for `1e4 < ncall <= 1e7`; both parts now use `ncall` as given.

Reproduce (needs numpy, matplotlib):

    python3 plot_hists.py new       # plots-new/      all 14 histograms, LO/NLO/NNLO, ratio to NLO
    python3 plot_hists.py overlap   # plots-overlap/  the 11 histograms also in old/, adds NNLO (old);
                                    #                 ratio to the NEW NLO central

Style follows Fig. 2 of 1506.02660 (LO yellow, NLO blue, NNLO red hatched bands = scale envelope).
Each directory has one PDF/PNG per observable and `all-<mode>.pdf`.

Notes
- Old NNLO band: envelope of 11, HH, 22. The old files have +-Infinity in some empty bins
  (mergedata artefact); these are read as 0. Column 5 of the old files (fraction of untrimmed
  seeds) is not used.
- `sig incl cuts`, `ptH-incl`, `yH-incl` have no counterpart in `old/` and appear only in `plots-new`.
  (`sig incl cuts` differs by a factor of ~4 between the codes (3.89 vs 0.92 pb at NNLO), so it is
  presumably defined differently; the VBF-cut cross section agrees to 0.07%.)
- `lo-*-nlo-ratio.top` differ from LO/NLO recomputed from the main files by up to 6e-3; the NLO and
  NNLO ratio files agree to 5e-8. The plots recompute the ratios.
