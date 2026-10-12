Per-piece integration of NLO VBF H+3j (>=3 jets, p1506 set-up, fixed mu = m_H): proVBFH-cs vs VBFNLO 3.0
12 Oct 2026. Sizing round; nothing but pilots was run (2.6 core-h in /ptmp/mpp/akarlber/pieces/pilot).

SPLIT
 VBFNLO (native, no patch): one run prints LO, "NLO virtual" = Born + V + I + K+P, and real - dipoles
   (QCD real emission). The p1506hist patch writes the same three as "integration nlo = -1, 1, 2".
 proVBFH-cs (patch cs_piece.patch, bitmask cs_piece in powheg.input, default all = 127; build in
   /ptmp/mpp/akarlber/pieces/proVBFH-cs with make ANALYSIS=p1506):
   bit 1 w (three-parton tree, (1,0)+(0,1) = LO H+3j; equals cs_order 1),
   bit 2 wv (one-loop + I + K+P of the radiating line), bit 4 e12 virtual-type (other line's
   (V+I)f + (K+P)xf at Born times the 3-parton event), bit 8 e12 dipole-type (minus other line's
   dipoles at its 3-parton point), bit 16 wr (4-parton real, one line), bit 32 wd (its six dipoles),
   bit 64 e3 ((1,1) 4-parton, both lines).  Comparable groups (all with cs_order 3, >=3j histogram of
   pwg-EXCL.top; the Born counterevent has 2 jets and never enters):
     LO              = 1    (cs_order 1)        <-> VBFNLO LO
     V+I+KP (corr.)  = 2+4  = 6                 <-> VBFNLO (Born+V+I+KP) - LO
     real - dipoles  = 8+16+32+64 = 120         <-> VBFNLO real emission
     all             = 127  (= production)
   The masks are applied after the matrix elements, so the cost of a piece is that of the full run
   (except LO). VEGAS adapts to the masked piece. e12 dipole-type belongs to "real - dipoles" because
   it is minus the dipole counterterm of a real-like configuration.

TECHNICAL CUTS
 VBFNLO: no dipole alpha, no parton pT cut in the generator (pTjmin = PT_JET_MIN only sets a split of the
   1/pT and a flat soft sampling, 20% of points with pT in [0, max(pTjmin,10)]); the one technical cut
   is m2min: weight 0 if 2 p_i.p_j < (0.1 GeV)^2 for any parton pair (phasespace/ps_vnj.F). Variation:
   VBFNLO_M2MIN in piece.env (GeV), implemented by an LD_PRELOAD library vbfnlo-ps/libps_m2min.so
   (ps_vnj.F with the env var; build: gfortran -cpp -O2 -fPIC -shared -ffixed-line-length-132
   -I<vbfnlo>/include -I<vbfnlo>/install/include ps_vnj.F -o libps_m2min.so). Propose 1.0 (and 0.01) GeV.
 proVBFH-cs: cs_cutoff (default 1d-5): four-parton points with 1-xp, z, 1-z of a dipole below it are
   dropped (about 2% of the four-parton points). Variation: cs_cutoff 1d-4 and 1d-6.
   Paired runs (same seeds as the baseline jobs, the points differ only where dropped) should give the
   difference with a much smaller error than two independent runs.

PILOT NUMBERS (>=3j, fb; 1.1M points per cs job, 2^20 points per VBFNLO iteration)
 cs   LO  118.51 +- 2.14 (0.03 core-h)   V+I+KP corr. 30.33 +- 0.79 (0.52 h)   real-dip -21.9 +- 5.7 (0.51 h)
      all 118.5 +- 5.9 (0.57 h)   real-dip with cs_cutoff 1d-4: -30.6 +- 5.9 (0.60 h, unpaired)
 VBFNLO (2^20: LO 118.9 +- 0.4, B+V+I+KP 148.9 +- 0.5, real-dip -17.0 +- 4.2; 0.19 h);
      with M2MIN = 1 GeV: 119.6 +- 0.4, 150.3 +- 0.6, real-dip -10.7 +- 6.8 (0.17 h)
 NOTE: the pieces agree where measured (cs V+I+KP 30.3 +- 0.8 vs VBFNLO 149.55 - 119.34 = 30.2);
 the real-dip pilots are too noisy to say more.

EXISTING VBFNLO PRODUCTION (800 jobs, 1.46 core-h each) ALREADY SPLITS INTO PIECES
 (vbfnlo_prod_pieces.py, vbfnlo_prod_iterations.py), mean over jobs +- error of the mean:
   p1506hist estimator (last iteration, plain sum):  LO 119.34, B+V+I+KP 149.552 +- 0.007,
        real-dip -24.33 +- 0.12   -> total 125.22
   VBFNLO's own result ("QCD real emission"): real-dip -23.31 +- 0.06 -> total 126.24 +- 0.06
   real-dip means per VEGAS iteration: -25.0(4.0), -22.2(1.0), -23.9(0.4), -24.9(0.3), -24.3(0.12)
   The own number is the inverse-variance weighted combination of the 5 iterations per job; it is
   pulled to -23.3, outside iterations 3-5, i.e. biased by the weighting (heavy-tailed real). The
   hist estimator (iteration 5) is the unbiased one, so the 125.22 stands, and the offset to cs is
   unchanged. Born+V+I+KP is identical in both (149.552).

SIZING (target <= 0.15 fb on each piece; per-job error scaled from the pilots; the production jobs
 (600 checked) have a scatter of the >=3j value 1.7x larger than the quoted error because of spikes,
 a factor 1.7 is included for real-dip)
 cs: LO 20 jobs x 0.4 h = 8 core-h;  V+I+KP 20 jobs (ncall2 1.2M) x 1.6 h = 32 core-h;
     real-dip 330 jobs (4.8M, 6.3 h) = 2100 core-h (730 core-h without the spike factor).
     Cutoff variations, paired: 40 jobs each for 1d-4 and 1d-6 = 500 core-h.
 VBFNLO: the 800 existing jobs already give every piece to <= 0.12 fb: 0 core-h. Variation M2MIN =
     1.0 and 0.01 GeV, paired with the seeds 2001-2100 of the existing jobs: 2 x 100 x 1.46 = 292 core-h.
 Total about 2900 core-h.

COMMANDS (not submitted). P=/u/akarlber/work/proVBFH/proVBFH-cs/production/h3j-crosscheck/pieces; B=/ptmp/mpp/akarlber/pieces/prod; mkdir -p $B/{lo,v,rd,rd-c4,rd-c6,vb-m1,vb-m001}/logs
 $P/mkjobs.sh cs $B/lo 7200 20 1 1 4800000 200000;       $P/submit.sh $B/lo/jobs.list 1-20 2:00:00 cs
 $P/mkjobs.sh cs $B/v 7300 20 3 6 1200000 100000;        $P/submit.sh $B/v/jobs.list 1-20 4:00:00 cs
 $P/mkjobs.sh cs $B/rd 7400 330 3 120 4800000 200000;    $P/submit.sh $B/rd/jobs.list 1-330 10:00:00 cs
 $P/mkjobs.sh cs $B/rd-c4 7400 40 3 120 4800000 200000 1d-4;  $P/submit.sh $B/rd-c4/jobs.list 1-40 10:00:00 cs
 $P/mkjobs.sh cs $B/rd-c6 7400 40 3 120 4800000 200000 1d-6;  $P/submit.sh $B/rd-c6/jobs.list 1-40 10:00:00 cs
 $P/mkjobs.sh vbfnlo $B/vb-m1 2001 100 23 1.0;    $P/submit.sh $B/vb-m1/jobs.list 1-100 3:00:00 vbfnlo
 $P/mkjobs.sh vbfnlo $B/vb-m001 2001 100 23 0.01; $P/submit.sh $B/vb-m001/jobs.list 1-100 3:00:00 vbfnlo
 (summarise with: summ.py <jobs.list>)
