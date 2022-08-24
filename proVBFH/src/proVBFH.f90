!----------------------------------------------------------------------
!----------------------------------------------------------------------
!---- proVBFH, Version 1.2.1 ------------------------------------------
!     
!  Program for the computation of inclusive cross sections in VBFH up
!  to N3LO in QCD and fully differential distributions up to NNLO.
!
!  Written by
!  Matteo Cacciari cacciari@lpthe.jussieu.fr 
!  Frederic Dreyer, frederic.dreyer@physics.ox.ac.uk
!  Alexander Karlberg, karlberg@physik.uzh.ch
!  Gavin Salam, gavin.salam@cern.ch
!  Giulia Zanderighi, g.zanderighi1@physics.ox.ac.uk
!
!  based on 
!  Phys.Rev.Lett. 115 (2015) no.8, 082002 [arXiv:1506.02660]
!  Phys.Rev.Lett. 117 (2016) no.7, 072001 [arXiv:1606.00840]
!
!  First version: 28/06/2017
!  Last edited: 28/06/2017
!     
!----------------------------------------------------------------------
!----------------------------------------------------------------------

program proVBFH
  use incl_vbfh
  real * 8 powheginput

  if (powheginput('#inclusive_only').eq.1d0 &
       & .or. powheginput('#qcd_order').eq.1d0 &
       & .or. powheginput('#qcd_order').eq.4d0 &
       &     ) then
     ! run the inclusive part only for total cross sections up to N3LO
     call run_inclusive()
  else
     ! run exclusive + inclusive for fully differential results up to NNLO
     call inclusive_init()
     call excl_vbfh()
  endif
end program proVBFH
