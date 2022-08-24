c -*- Fortran -*-

c The user must set nlegborn to the appropriate value for his process
c in file nlegborn.h in the user directory
c
c maxprocborn and maxprocreal must be greater or equal to the number of
c independent flavour structures for the born and real process.
c
c flst_nborn and flst_nreal should be set to the  number of
c independent flavour structures for the born and real process.
c
c flst_born and flst_real should be filled with the flavour structure
c of the born and real subprocesses.
c The meaning is: flst_real(j,k) is the flavour of leg j in subprocess k.
c The flavour is taken incoming for the two incoming particles and outgoing
c for the outgoing particles.
c The flavour number is according to PDG conventions, EXCEPT for gluons,
c where we take 0 instead of 21.
c As an example, for p p -> (Z->e+e-)+2 j,
c [1,0,11,-11,1,0] represents d g -> e- e+ d g
c Notice that only one ordering of final state flavour can appear for
c a given final state.
c It is required that legs in the Born and real processes
c should be ordered as follows
c leg 1: incoming parton with positive rapidity
c leg 2: incoming parton with negative rapidity
c from leg 3 onward: final state particles, ordered as follows:
c Colorless particles first, massive coloured particles,
c massless coloured particles.
c The flavour of colored neutral and massive colored particles should
c be the same for all real and born subprocesses.
c This is the case for most QCD NLO calculation.


      integer flst_alr_tag(nlegreal,2*maxalr),flst_real_tag(nlegreal,maxprocreal)
      integer alr_tag, realequiv_tag
      common/pwhg_flst2/flst_alr_tag, alr_tag, flst_real_tag, realequiv_tag
      save /pwhg_flst2/
