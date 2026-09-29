!----------------------------------------------------------------------
! proVBFH-cs: VBF H with a line-by-line projection-to-Born
! (docs/DESIGN.md). Stage 1 (NLO): run with the same powheg.input and
! vbfnlo.input as proVBFH, plus
!   cs_part   1: inclusive part (structure functions, Born-level events),
!             2: exclusive part, (1,0) + (0,1) with counterevents
!                (default)
!   cs_npow   sampling power for 1-xp and z (default 2)
!   cs_cutoff invariant cutoff on 1-xp, z, 1-z (default 1d-8)
! The two parts write separate histogram files, to be added.
!----------------------------------------------------------------------
program provbfh_cs
  use types, only: dp
  use incl_parameters
  use incl_vbfh, only: inclusive_init, run_inclusive
  use phase_space, only: set_beams
  use integration
  use cs_exclusive
  implicit none
  integer, parameter :: ndim = 13
  real(dp) :: region(2*ndim), integ, err, chi2, t0, t1
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
  if (powheginput('#cs_npow') > 0) excl_npow = nint(powheginput('#cs_npow'))
  if (powheginput('#cs_cutoff') > 0) excl_cutoff = powheginput('#cs_cutoff')
  if (powheginput('#cs_flavcheck') > 0) excl_flavcheck = 200
  call set_beams(sqrts)

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
  write(6,'(a,3i14)') ' points, rejected by cutoff (per line), NaN: ', excl_stats
  write(6,'(a,f12.1,a)') ' proVBFH-cs exclusive part: CPU ', t1 - t0, ' s'
end program provbfh_cs
