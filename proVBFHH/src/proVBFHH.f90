!----------------------------------------------------------------------
!----------------------------------------------------------------------
!---- proVBFHH, Version 1.2.0 ------------------------------------------
!     
!  Program for the computation of inclusive cross sections in VBFHH up
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
!  [arXiv:1811.07918], submitted to Phys.Rev.Lett.
!  Phys.Rev.Lett. 115 (2015) no.8, 082002 [arXiv:1506.02660]
!
!  as well as for the N3LO capabilities, on
!  [arXiv:1811.07906], submitted to Phys.Rev.D
!  Phys.Rev.Lett. 117 (2016) no.7, 072001 [arXiv:1606.00840]
!
!  First version: 21/11/2018
!  Last edited: 21/11/2018
!     
!----------------------------------------------------------------------
!----------------------------------------------------------------------

program proVBFHH
  use incl_vbfhh
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
     call excl_vbfhh()
  endif
end program proVBFHH
