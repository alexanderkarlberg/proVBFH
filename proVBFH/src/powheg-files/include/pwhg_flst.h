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

      integer flst_nborn, flst_born(nlegborn,maxprocborn),
     1     flst_borntags(nlegborn,maxprocborn),
     2     flst_bornres(nlegborn,maxprocborn),
     3     flst_nreson,flst_reslist(nlegborn)
c res arrays contain pointers to the mother resonance for
c resonace's decay products, zero for the other particles.
      integer flst_nreal, flst_real(nlegreal,maxprocreal),
     1     flst_realtags(nlegreal,maxprocreal),
     2     flst_realres(nlegreal,maxprocreal)
      integer flst_nalr,flst_alr(nlegreal,maxalr),
     1     flst_alrtags(nlegreal,maxalr),
     2     flst_alrres(nlegreal,maxalr)
      integer flst_nregular,flst_regular(nlegreal,maxprocreal),
     1     flst_regulartags(nlegreal,maxprocreal),
     2     flst_regularres(nlegreal,maxprocreal)
      integer flst_uborn(nlegborn,maxalr),
     1     flst_uborntags(nlegborn,maxalr),
     2     flst_ubornres(nlegborn,maxalr)
      integer flst_mult(maxalr), flst_ubmult(maxalr),
     1        flst_emitter(maxalr)
      integer maxregions
      parameter (maxregions=nlegreal*(nlegreal-1)/2)
c     list of all regions      
      integer flst_allreg(2,0:maxregions,maxalr)
c     pointer from flst_alr to flst_born
      integer flst_alr2born(maxalr)
c     pointer from flst_born to flst_alr
c     born2alr(0,k) = number of alpha_r with underlying born k
c     born2alr(1.. born2alr(0,k),k)= pointer to corresponding alpha_r
      integer flst_born2alr(0:maxalr,maxprocborn)
      integer flst_lightpart
      logical flst_sonof,flst_isfs,flst_isres(nlegborn)
      integer flst_cur_iborn,flst_cur_alr 
      common/pwhg_flst/
     1     flst_nborn, flst_born,flst_borntags,flst_bornres,
     2     flst_nreal,flst_real,flst_realtags,flst_realres,
     3     flst_nalr,flst_alr,flst_alrtags,flst_alrres,
     4     flst_nregular,flst_regular,flst_regulartags,flst_regularres,
     5     flst_uborn,flst_uborntags,flst_ubornres,flst_mult,
     7     flst_ubmult,flst_emitter,flst_allreg,flst_alr2born,
     8     flst_born2alr,flst_lightpart,
     9     flst_nreson,flst_reslist,
     1     flst_isres,flst_cur_iborn,flst_cur_alr
      save /pwhg_flst/
