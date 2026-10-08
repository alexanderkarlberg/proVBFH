!----------------------------------------------------------------------
! proVBFH-cs: VBF H with a line-by-line projection-to-Born
! (docs/DESIGN.md). Stage 1 (NLO): run with the same powheg.input and
! vbfnlo.input as proVBFH, plus
!   cs_part   1: inclusive part (structure functions, Born-level events),
!             2: exclusive part, (1,0) + (0,1) with counterevents
!                (default)
!   cs_npow   sampling power for 1-xp and z; 0 (default): logarithmic
!   cs_hardfrac  fraction of the line radiations in the hard channel (small
!             xp, z ~ 1/2; for the high-pT tails), also the first step of
!             the four-parton generator; default 0
!   cs_hardfrac2  the same for the second emission of the four-parton
!             generator; default 0
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
  use cs_nlo2, only: nlo2_ncount, nlo2_ncut, nlo2_emul_kappa
  use cs_memo, only: memo_hits, memo_misses, memo_on
  use cs_chan, only: chan_init, chan_multi, chan_trivial
  use incl_parameters, only: incl_nscale, incl_scr, incl_scf
  use cs_dipoles, only: spin_avg, four_hard, four_hard2
  use matrix_element, only: incl_only11
  implicit none
  integer, parameter :: maxdim = 20
  integer :: ndim
  real(dp) :: region(2*maxdim), integ, err, chi2, t0, t1
  real(dp) :: powheginput
  external powheginput
  integer :: part, ilast, i
  common/last_integ/ilast
  character(len=20) :: pwgprefix
  integer :: lprefix
  common/cpwgprefix/pwgprefix,lprefix
  character(len=60) :: histname

  call cpu_time(t0)
  part = 2
  if (powheginput('#cs_part') > 0) part = nint(powheginput('#cs_part'))

  if (part == 1) then
     if (powheginput('#incl_only11') == 1) incl_only11 = .true.
     ! on-the-fly scale variations (the same points as the exclusive part)
     if (powheginput('#cs_scales') > 1) then
        incl_nscale = nint(powheginput('#cs_scales'))
        if (incl_nscale /= 3 .and. incl_nscale /= 7) stop 'cs_scales must be 1, 3 or 7'
        incl_scr(1:incl_nscale) = sc_r(1:incl_nscale)
        incl_scf(1:incl_nscale) = sc_f(1:incl_nscale)
        write(6,'(a,i2,a)') ' proVBFH-cs inclusive part: ', incl_nscale, ' scale points, weights W1, W2, ...:'
        write(6,'(7(a,f4.2,a,f4.2,a))') (' (', sc_r(i), ',', sc_f(i), ')', i = 1, incl_nscale)
     endif
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
  if (powheginput('#cs_hardfrac') > 0) excl_hardfrac = powheginput('#cs_hardfrac')
  four_hard = excl_hardfrac
  if (powheginput('#cs_hardfrac2') > 0) four_hard2 = powheginput('#cs_hardfrac2')
  if (powheginput('#cs_cutoff') > 0) excl_cutoff = powheginput('#cs_cutoff')
  if (powheginput('#cs_flavcheck') > 0) excl_flavcheck = 200
  if (powheginput('#cs_phspcuts') == 0) excl_phspcuts = .false.
  if (powheginput('#cs_order') > 0) excl_order = nint(powheginput('#cs_order'))
  if (powheginput('#cs_only2') == 1) excl_only2 = .true.
  if (powheginput('#cs_no20') == 1) excl_no20 = .true.
  if (powheginput('#cs_dump2') == 1) then
     excl_dump2 = .true.
     if (powheginput('#cs_dump2min') > 0) dump2_min = powheginput('#cs_dump2min')
     open(79, file='cs_dump2.dat', status='replace')
  endif
  if (powheginput('#cs_estimate') > 0) excl_estimate = nint(powheginput('#cs_estimate'))
  call cs_set_ggbug(.false.)
  if (powheginput('#cs_estimu') > 0) excl_estimu = nint(powheginput('#cs_estimu'))
  if (excl_order >= 2 .or. powheginput('#cs_testlimits') >= 1 .or. powheginput('#cs_testvirt') == 1 &
       & .or. powheginput('#cs_testborn2') == 1 .or. powheginput('#cs_testlines') == 1 &
       & .or. powheginput('#cs_testlines11') == 1) then
     call cs_excl_setup2(.true.)
     call cpu_time(t1)
     write(6,'(a,f10.2,a)') ' proVBFH-cs stage-2 set-up: CPU ', t1 - t0, ' s'
  endif
  ! cs_estimate 6: k_T cut [GeV] of the emulated old-code treatment (set after
  ! the stage-2 set-up, whose grouping of the real entries uses the full matrix elements)
  if (excl_estimate == 6) then
     nlo2_emul_kappa = 0.1_dp
     if (powheginput('#cs_kappa') > 0) nlo2_emul_kappa = powheginput('#cs_kappa')
  endif
  ! cs_estimate 7 (Born-type weights only): no four-parton matrix elements
  if (excl_estimate == 7) nlo2_emul_kappa = huge(1.0_dp)
  ! on-the-fly scale variations: cs_scales 3 (symmetric) or 7 points
  if (powheginput('#cs_scales') > 0) then
     call cs_scales_setup(nint(powheginput('#cs_scales')))
     if (powheginput('#cs_scalecheck') == 1) excl_scalecheck = .true.
     if (powheginput('#cs_scalecheck') == 2) then
        if (excl_nscale /= 7) stop 'cs_scalecheck 2 needs cs_scales 7'
        excl_scalecheck2 = .true.
        sc_r(1:7) = [1.0_dp, 0.5_dp, 2.0_dp, 0.5_dp, 2.0_dp, 1.0_dp, 1.0_dp]
        sc_f(1:7) = [1.0_dp, 0.5_dp, 2.0_dp, 0.5_dp, 2.0_dp, 1.0_dp, 1.0_dp]
        write(6,'(a)') ' cs_scalecheck 2: W1-W3 (1,1), (1/2,1/2), (2,2) with the beta0 shift, W4, W5 (1/2,1/2), (2,2)' &
             & //' with the virtual evaluated directly, W6 = W2 - W4, W7 = W3 - W5 (W6, W7: scale labels not meaningful)'
     endif
     write(6,'(a,i2,a)') ' proVBFH-cs: ', excl_nscale, ' scale points (mu_R/mu_R0, mu_F/mu_F0), weights W1, W2, ...:'
     write(6,'(7(a,f4.2,a,f4.2,a))') (' (', sc_r(i), ',', sc_f(i), ')', i = 1, excl_nscale)
  endif
  ! diagnostic channel split (cs_chan; environment CHAN_BOSON, CHAN_INIT, CHAN_MULTI)
  call chan_init()
  if (excl_estimate /= 0 .and. (chan_multi .or. .not. chan_trivial())) stop 'cs_chan: not with cs_estimate'
  if (chan_multi) then
     if (excl_nscale /= 1) stop 'cs_chan: CHAN_MULTI needs cs_scales 1'
     excl_nscale = 9
     memo_on = .true.
     sc_r = 1; sc_f = 1
  endif
  if (powheginput('#cs_testlimits') >= 1) then
     call cs_excl_testlimits(nint(powheginput('#cs_testlimits')))
     stop
  endif
  if (powheginput('#cs_testborn2') == 1) then
     call cs_excl_testborn2()
     stop
  endif
  if (powheginput('#cs_testlines') == 1) then
     call cs_excl_testlines()
     stop
  endif
  if (powheginput('#cs_testlines11') == 1) then
     call cs_excl_testlines11()
     stop
  endif
  if (powheginput('#cs_testvirt') == 1) then
     call cs_excl_testvirt()
     stop
  endif
  if (powheginput('#cs_spinavg') == 1) spin_avg = .true.
  if (excl_order >= 2) then
     if (powheginput('#cs_spikemin') > 0) spike_min = powheginput('#cs_spikemin')
     open(80, file='cs_spikes.dat', status='replace')
  endif
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
  if (excl_order >= 2 .and. excl_order /= 10 .and. excl_order /= 11 .and. excl_order /= 13) then
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
  if (excl_order >= 2) write(6,'(a,i14)') ' stage 2: four-parton points dropped by the technical cut: ', nlo2_ncut
  if (excl_nscale > 1) write(6,'(a,2i14)') ' scale variations: matrix-element cache hits, misses: ', memo_hits, memo_misses
  if (excl_scalecheck) write(6,'(a,i12,a,2es10.2)') ' scale check: ', scalecheck_n, &
       & ' shifted V+I against direct, max deviation relative to |V+I|, to |born|:', scalecheck_dev
  write(6,'(a,f12.1,a)') ' proVBFH-cs exclusive part: CPU ', t1 - t0, ' s'
end program provbfh_cs
