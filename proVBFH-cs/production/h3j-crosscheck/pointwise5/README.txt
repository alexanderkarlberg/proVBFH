VBFNLO 3.0 NLO H+3j (process 110): do real + sum(dipoles) vanish in the soft / collinear limits?  (13 Oct)
Verdict: YES.  The 10-90% residuals of pointwise4 (PW_LIM) were harness artefacts (two bugs, below).

Files
  m2s_qqh4j_lim5.patch   against /ptmp/mpp/akarlber/chsplit/vbfnlo/src/amplitudes/hjjj/m2s_qqh4j.F (supersedes pointwise4's
                         m2s_qqh4j_dump.patch; same LD_PRELOAD wrapper m2s_wbfh4j -> m2s_wbfh4j_in, PW_LIM=n limit sequences)
  qqh4q_sub_dump.patch   unchanged from pointwise4 (pw_rec labels the 4-quark subprocesses)
  build_run.sh           builds libt.so (flags incl. -DWITH_NLO, which the first attempt forgot: R=0 without it) and runs
                         PW_LIM=N points; needs a vbfnlo.dat (PROCESS 110, NLO only, NLO_POINTS 16) in WORKDIR/run
  agg_limits.py          summary: only singular points are used (R(1e-6)/R(1e-4) > 20 and R(1e-4)/R(1e-2) > 5)
  results.txt            summary table with cuts (60 points); results_nocuts.txt: same with PW_NOCUTS=1 (old behaviour)
Output files in run/: fort.99 (LIM: lambda, R, sum D, ratio), fort.98 (LK per-subprocess R and sum D; LE per-dipole values for
  one flavour key; SIJ invariants), fort.96 (JET: njets and pT of the jets of the real and each of the 24 mapped configs).

Bugs of the pointwise4 harness
 1. The mapped (dipole) configurations were only checked for on-shellness/momentum conservation, not for VBFNLO's jet cuts.
    In the integrand lok(nd) = Cuts(...) with njets >= 3 drops every dipole whose Born configuration does not have 3 jets.
    In a soft limit of parton k, the dipoles that do not emit k keep k as a (soft) Born parton: Born -> infinity, Born config
    has 2 hard jets, fails the cut in VBFNLO (D = 0) but was kept in the harness: e.g. soft slot 5, 4q entries R~1e-7,
    D~-6.7 (dipole (1,6,3), Born g u -> u d dbar with soft dbar).  Fix: defineJets/Cuts on pl,vl for all 25 configs
    (same calls as reals_phasespacepoint in vbfnlo_main.F90).
 2. The final-state collinear sequences scaled the transverse momenta w.r.t. the pair 3-momentum axis, which makes
    back-to-back (not collinear) pairs when one parton has negative longitudinal momentum: s_35 stayed at 1.5e5 GeV^2.  Fix
    (pw_collpair): same pair 3-momentum, pair mass^2 -> lambda*s, same rest-frame direction, recoil as before.
 The "ISR a1 || s3" sequence is not singular when p3 points to the beam of incoming 2 (R flat ~1e-15; slot 3 only belongs to
 incoming 1 in create_ikj: dipoles (a,3,.) exist, (b,3,.) do not, by design); the summary keeps only singular points.

Result (60 points, singular ones only; median / max of |R+sum D|/R)
 see results.txt.  All limits converge: soft 5e-5 / 6e-4 (flat plateau for lambda<1e-4), coll (5,6) 4e-5 / 7e-4,
 coll (3,5), (4,6) 6e-5 / 8e-4 at lambda=1e-8, ISR 1e-4 median at 1e-8 (sqrt(lambda)-like convergence, max 1e-2 at 1e-8 is numerics).
 Remaining feature: in the soft limit (R+D)/R is a lambda-independent -1e-7...-7e-4 (D over-subtracts), flavour-independent at
 a given point, larger for identical-flavour keys.  Most likely the colour-suppressed inter-line (non-factorisable) terms of the
 real that dipsub/dipole Born (msq_b(1|2), line-factorised) do not contain.  Bound: <= 7e-4 of the singular soft real.

Notes: scratch runs were in the login-node /tmp scratchpad (not /ptmp); with PW_NOCUTS=1 (old behaviour, but fixed collinear) the residuals are 0.2-0.8 median (results_nocuts.txt), i.e. the missing cuts are the dominant bug.
 variable declared double (gfortran warning only, works).
