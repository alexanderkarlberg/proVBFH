Pointwise comparison with VBFNLO 3.0 (12 Oct).
1. vbfnlo_m2s_qqh3j_dump.patch: dump patch for m2s_qqh3j.F of the (chsplit-patched) VBFNLO copy. Build as a
   scratch shared lib (gfortran -fPIC -shared, flags of libhjjj, link -lVBFNLOUtilities -lVBFNLOLoops
   -lVBFNLOAmplitudes) and run the stock vbfnlo with LD_PRELOAD=libdump.so in a copy of the fixmh run dir
   (LO_POINTS 12, NLO_POINTS 16): writes pw_vbfnlo.dat (momenta, Born, finite virtual, V+I, K+P convolutions).
2. pw_cs.f90: reads pw_vbfnlo.dat, evaluates proVBFH-cs Born, virtual, I, K+P at the same momenta.
   Build in a scratch copy of proVBFH-cs (objects of all but cs_main.o), run in a dir with the production
   powheg.input / vbfnlo.input of fixmh p1506.
