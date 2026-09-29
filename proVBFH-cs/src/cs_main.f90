!----------------------------------------------------------------------
! proVBFH-cs: VBF H with a line-by-line projection-to-Born
! (docs/DESIGN.md). Stage 1 (NLO): run with the same powheg.input and
! vbfnlo.input as proVBFH, plus
!   cs_part   1: inclusive part (structure functions, Born-level events),
!             2: exclusive part, (1,0) + (0,1) with counterevents
!                (default)
!   cs_npow   sampling power for 1-xp and z; 0 (default): logarithmic
!   cs_cutoff invariant cutoff on 1-xp, z, 1-z (default 1d-6)
!   cs_order  1: the O(alpha_s) exclusive part (NLO, default);
!             2: also the O(alpha_s^2) (2,0) + (0,2) (stage 2)
!   cs_testlimits 1, cs_testvirt 1: stage-2 unit tests (real against
!             dipoles in the singular limits; scale dependence of V + I)
! The two parts write separate histogram files, to be added.
!----------------------------------------------------------------------
program provbfh_cs
  use types, only: dp
  use incl_parameters
  use incl_vbfh, only: inclusive_init, run_inclusive
  use phase_space, only: set_beams
  use integration
  use cs_exclusive
  use cs_nlo2, only: nlo2_ncount
  use cs_dipoles, only: spin_avg
  implicit none
  integer, parameter :: maxdim = 20
  integer :: ndim
  real(dp) :: region(2*maxdim), integ, err, chi2, t0, t1
  real(dp) :: powheginput
  external powheginput
  integer :: part, ilast
  common/last_integ/ilast
  character(len=20) :: pwgprefix
  integer :: lprefix
  common/cpwgprefix/pwgprefix,lprefix
  character(len=60) :: histname

  call cpu_time(t0)
  part = 2
  if (powheginput('#cs_part') > 0) part = nint(powheginput('#cs_part'))

  if (part == 1) then
     call run_inclusive()
     call cpu_time(t1)
     write(6,'(a,f12.1,a)') ' proVBFH-cs inclusive part: CPU ', t1 - t0, ' s'
     stop
  endif

  ! exclusive part
  call inclusive_init()          ! parameters, PDFs, alpha_s, scales
  call cs_init_vbfnlo()          ! VBFNLO couplings from vbfnlo.input
  call cs_excl_setup()
  call set_beams(sqrts)
  if (powheginput('#cs_npow') >= 0) excl_npow = nint(powheginput('#cs_npow'))
  if (powheginput('#cs_cutoff') > 0) excl_cutoff = powheginput('#cs_cutoff')
  if (powheginput('#cs_flavcheck') > 0) excl_flavcheck = 200
  if (powheginput('#cs_phspcuts') == 0) excl_phspcuts = .false.
  if (powheginput('#cs_order') > 0) excl_order = nint(powheginput('#cs_order'))
  if (powheginput('#cs_only2') == 1) excl_only2 = .true.
  if (powheginput('#cs_dump2') == 1) then
     excl_dump2 = .true.
     if (powheginput('#cs_dump2min') > 0) dump2_min = powheginput('#cs_dump2min')
     open(79, file='cs_dump2.dat', status='replace')
  endif
  if (powheginput('#cs_estimate') > 0) excl_estimate = nint(powheginput('#cs_estimate'))
  if (powheginput('#cs_estimu') > 0) excl_estimu = nint(powheginput('#cs_estimu'))
  if (excl_order >= 2 .or. powheginput('#cs_testlimits') >= 1 .or. powheginput('#cs_testvirt') == 1) then
     call cs_excl_setup2(.true.)
     call cpu_time(t1)
     write(6,'(a,f10.2,a)') ' proVBFH-cs stage-2 set-up: CPU ', t1 - t0, ' s'
  endif
  if (powheginput('#cs_testlimits') >= 1) then
     call cs_excl_testlimits(nint(powheginput('#cs_testlimits')))
     stop
  endif
  if (powheginput('#cs_testvirt') == 1) then
     call cs_excl_testvirt()
     stop
  endif
  if (powheginput('#cs_spinavg') == 1) spin_avg = .true.
  if (powheginput('#cs_replay') == 1) then
     call cs_excl_replay()
     stop
  endif
  if (powheginput('#cs_dump') == 1) then
     excl_dump = .true.
     open(77, file='cs_dump.dat', status='replace')
  endif

  ! the seven dimensions of the four-parton events are not adapted (their
  ! density enters the multichannel weight)
  ndim = 13
  if (excl_order >= 2) then
     ndim = 20
     jfreeze = 14
  endif
  region(1:ndim) = 0
  region(ndim+1:2*ndim) = 1
  ilast = 0
  outgridfile = 'grids-excl'//trim(seedstr)//'.dat'
  outgridtopfile = 'grids-excl'//trim(seedstr)//'.top'
  if (.not. readin) then
     writeout = .true.
     excl_fill = .false.
     call vegas(region, ndim, cs_excl_dsigma, 0, ncall1, itmx1, 0, integ, err, chi2)
     writeout = .false.
  endif
  ingridfile = 'grids-excl'//trim(seedstr)//'.dat'
  excl_fill = do_analysis_call
  if (excl_fill) call init_hist
  call vegas(region, ndim, cs_excl_dsigma, 1, ncall2, itmx2, 0, integ, err, chi2)
  if (excl_fill) then
     call pwhgsetout
     histname = pwgprefix(1:lprefix)//'-EXCL'//trim(seedstr)
     call pwhgtopout(histname)
  endif
  call cpu_time(t1)
  write(6,'(a,es13.5,a,es10.3,a)') ' proVBFH-cs exclusive part: sum |w1|+|w2| = ', integ, ' +- ', err, ' pb'
  write(6,'(a,4i14)') ' points, line radiations cut off, NaN, points failing all cuts: ', excl_stats
  if (excl_order >= 2) write(6,'(a,2i14)') ' stage 2: real and Born matrix-element calls: ', nlo2_ncount
  write(6,'(a,f12.1,a)') ' proVBFH-cs exclusive part: CPU ', t1 - t0, ' s'
end program provbfh_cs
