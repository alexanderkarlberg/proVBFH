! proVBFH-cs line-level R/D limit test (cs_excl_testlimits(1) = nlo2_test_limits) for the report
program pw4_lim
  use incl_vbfh, only: inclusive_init
  use incl_parameters, only: sqrts
  use phase_space, only: set_beams
  use cs_exclusive, only: cs_excl_setup, cs_excl_setup2, cs_excl_testlimits
  implicit none
  call inclusive_init()
  call cs_init_vbfnlo()
  call cs_excl_setup()
  call set_beams(sqrts)
  call cs_excl_setup2(.false.)
  call cs_excl_testlimits(1)
end program
