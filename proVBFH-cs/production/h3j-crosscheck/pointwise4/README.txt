Pointwise comparison of the 4-parton (real-emission) ingredients of NLO H+3j, proVBFH-cs vs VBFNLO 3.0 (12 Oct).
Follow-up of ../pointwise/ (Born, V+I, K+P agreed).  Results in results.txt.

Scratch (this session): /tmp/claude-*/.../scratchpad/pw4 (login-node local disk; /ptmp was at its hard file quota,
1,572,025 of 1,572,864 files, so nothing could be created there).

vbfnlo/   LD_PRELOAD instrumentation of the (chsplit-patched) VBFNLO 3.0 copy (patches against
          /ptmp/mpp/akarlber/chsplit/vbfnlo/src/amplitudes/hjjj/{m2s_qqh4j.F,qqh4q_sub.F}):
          m2s_wbfh4j is renamed m2s_wbfh4j_in and wrapped.  For every 3rd call with nlo=1 and some lok(nd) true
          (up to 800 points) it writes pw4_vbfnlo.dat: EV (xi, alpha_s, lok(1:25)), P (all 25 momentum sets with their
          lok flag), D (per subprocess: real and the 24 dipole terms), and for the first 250 points with lok(1) the
          sums 'A' of the real matrix elements (without PDFs) over all 13 classes x 24 assignments of the four outgoing
          momenta to the outgoing slots, keyed by the flavours at (p1, p2 | k1..k4).  Env PW_VAR=1/2/3 integrates
          the real only (all / pure quark / gluon channels) instead of real+dipoles; PW_LIM=n runs limit sequences
          (soft, collinear, ISR) for the first n points (harness inconclusive, see below).
          Build: gfortran -fPIC -shared <flags of libhjjj> -I.../src/{utilities,include,amplitudes/hjjj} m2s_qqh4j.F
          qqh4q_sub.F -lVBFNLOUtilities -lVBFNLOLoops -lVBFNLOAmplitudes, run stock vbfnlo (PROCESS 110, NLO only,
          NLO_POINTS 14) with LD_PRELOAD.  Variants used: linterference=.false. and the CC colour flows 3,4 of the
          4-quark real switched off in qqh4q_sub.F (sed on the parameter / the two CCmatrix_r(..,3|4,..) terms).
cs/       drivers linked against the objects of a scratch build of proVBFH-cs (buildpw.sh; buildpw2.sh links the
          p1506 analysis instead of user_analysis.o):
          pw4_real2.f90  real MEs: for each point all 2333 real entries x 24 leg assignments through cs_real_me (setreal),
                         keyed like VBFNLO; writes per (a,b) sums (VBFNLO, CS all entries, CS entries used by the
                         (2,0)/(0,2)/(1,1) classes); full key set for CS (restricting CS to VBFNLO's keys loses
                         CS-only flavour assignments and gives spurious deficits)
          pw4_jets.f90   cs_passes (fastjet anti-kt, buildjets, vbfcuts of analysis/p1506_analysis.f) and njets >= 3
                         against VBFNLO's lok(nd) for the real point and the 24 dipole-mapped configurations
          pw4_lim.f90    proVBFH-cs nlo2_test_limits (line-level R/D limits)
          pw4_dbg.f90    contributions of the individual real entries to one flavour key
py/       jetcmp.py, jetcmp3.py: independent python anti-kt (E-scheme) + cuts (reproduces VBFNLO's lok exactly with
          njets >= 3); ana9.py: real-ME summary; ana8.py: variants.  Inputs: powheg.input/vbfnlo.input of csrun
          (p1506 production, WWIDTH 2.141, ZWIDTH 2.4952).

Things to know when reading the outputs
- VBFNLO's 4-quark CC branch (qqh4q_sub.F) sums four colour/flavour columns into one subprocess keyed by the flavours of
  column 1 and with column 1's PDFs.  Per-flavour-assignment comparisons with proVBFH-cs therefore differ in the labelling
  (s<->d and the like), not in physics.  The invariant comparison is the flavour-blind sum per incoming pair (a,b)
  (the cuts do not depend on flavour), done in ana9.py.
- The hjjj cuts of VBFNLO include njets >= 3 (min_jets); the dijet VBF cuts alone are looser (without njets >= 3: 39 extra real
  and 2664 extra mapped configurations on top of VBFNLO's 294 / 4559 accepted).
- The PW_LIM limit sequences (rescaled momenta + boost/scale recoil, VBFNLO's own ptilde for the dipole kinematics)
  gave residuals (R+sum D)/R of O(10-90%) that do not vanish for lambda -> 0 and depend on the point; a natural VBFNLO
  point at 2 p.p ~ 2 GeV^2 (lchsb debug print) gives 2.5e-3.  Not understood (harness artefact or a VBFNLO subtlety),
  not used.  The proVBFH-cs line-level test (pw4_lim) converges (R/D -> 1 like sqrt(lambda), max |R/D-1| 2.0% at 1e-6).
- The two codes use different dipole sets (VBFNLO: 20-24 dipoles per point, spectators include the other extra parton and
  both incoming partons, mapped configurations differ; proVBFH-cs: 6 maps per line, spectator the third parton of the line
  or the incoming one), so individual dipoles and mapped momenta cannot be compared one to one; only R - sum D at the
  same point (finite, not equal) and its singular limits are meaningful.
