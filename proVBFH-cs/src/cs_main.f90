!----------------------------------------------------------------------
! proVBFH-cs: VBF H with a line-by-line projection-to-Born
! (docs/DESIGN.md). One run gives the inclusive part (structure
! functions, Born-level events; 7 dimensions) and the exclusive part
! (O(alpha_s) and O(alpha_s^2) real and virtual corrections with
! counterevents), each with its own VEGAS integration, and adds their
! histograms. The output is the sum only: provbfh-<iseed>-<order>-W<k>.top
! for the k scale points (without -W<k> for one) and the cross section with
! the cuts of the analysis in provbfh-<iseed>-<order>.dat; <iseed> has four
! digits, <order> is LO, NLO or NNLO. With cs_write_parts 1 the parts are
! also written, provbfh-<iseed>-<order>-incl-W<k>.top and ...-excl-W<k>.top.
! Same powheg.input and vbfnlo.input as proVBFH, plus
!   inclusive_only 1: only the inclusive part, at the order qcd_order (any order)
!   qcd_order 1: LO (inclusive part only), 2: NLO, 3: NNLO (default); the
!             order of the inclusive part and of the run (as in the other
!             proVBFH codes; cs_order is accepted as a deprecated alias)
!   cs_scales number of on-the-fly scale points: 1, 3 (default; (1,1),
!             (1/2,1/2), (2,2)) or 7 (needs runningscales 1 for > 1)
!   ncall1, ncall2, itmx1, itmx2  VEGAS set-up (grid and main pass), used as
!             given by the inclusive and the exclusive part
!   cs_npow   sampling power for 1-xp and z; 0 (default): logarithmic
!   cs_hardfrac  fraction of the line radiations in the hard channel (small
!             xp, z ~ 1/2; for the high-pT tails), also the first step of
!             the four-parton generator; default 0.3
!   cs_hardfrac2  the same for the second emission of the four-parton
!             generator; default 0.7
!   cs_hard2mode  0: cs_hardfrac2 after every first step;
!             1 (default): only after the first step's hard channel
!             (targeted at hard double emissions from one line, small xp)
!   cs_smallx  fraction h of the Born phase-space points in the small-x
!             channels (cs_smallx.f90; default 0.3, 0 = off), cs_smallx_x0
!             their x limit (default 0.03)
!   cs_cutoff invariant cutoff on 1-xp, z, 1-z (default 1d-5)
!   cs_4qint  1 (default): identical-quark interference in the four-quark
!             reals (compreal_hqqqq_new.f); 0: off; 2: the interference
!             alone (exclusive part, provbfh-<iseed>-NNLO-excl-W<k>.top, which
!             is then written; needs qcd_order 3)
!   cs_write_parts  1: also write the inclusive and the exclusive part
!             (-incl-, -excl-) besides the sum; default 0
!----------------------------------------------------------------------
program provbfh_cs
  use types, only: dp
  use incl_parameters
  use incl_vbfh, only: inclusive_init, run_inclusive
  use phase_space, only: set_beams
  use integration
  use cs_exclusive
  use cs_nlo2, only: nlo2_ncount, nlo2_ncut, nlo2_intonly
  use cs_memo, only: memo_hits, memo_misses, memo_on
  use incl_parameters, only: incl_nscale, incl_scr, incl_scf
  use cs_smallx, only: sx_h, sx_x0, sx_dim0
  use cs_dipoles, only: four_hard, four_hard2, four_hard2mode
  implicit none
  integer, parameter :: maxdim = 25
  integer :: ndim
  real(dp) :: region(2*maxdim), integ, err, chi2, t0, t1
  real(dp) :: powheginput
  external powheginput
  integer :: part, ilast, i, nscales
  common/last_integ/ilast
  character(len=60) :: histname
  logical :: readin_job, write_parts
  integer :: idum_job, corder
  integer :: idum
  common/ranno/idum
  real(dp) :: sig_i, err_i, sig_k, err_k
  character(len=60) :: runbase    ! provbfh-<iseed>-<order>: stem of the output files
  character(len=20) :: seedtxt
  character(len=4), parameter :: oname(3) = ['LO  ', 'NLO ', 'NNLO']

  call cpu_time(t0)
  corder = 3
  if (powheginput('#qcd_order') > 0) corder = nint(powheginput('#qcd_order'))
  if (powheginput('#cs_order') > 0) then
     ! deprecated alias of qcd_order
     if (powheginput('#qcd_order') > 0 .and. nint(powheginput('#cs_order')) /= corder) &
          stop 'qcd_order and cs_order (deprecated alias of qcd_order) differ'
     if (powheginput('#qcd_order') <= 0) corder = nint(powheginput('#cs_order'))
     write(6,'(a)') ' proVBFH-cs: WARNING: cs_order is deprecated, use qcd_order'
  endif
  if (corder < 1 .or. corder > 3) stop 'qcd_order must be 1 (LO), 2 (NLO) or 3 (NNLO)'
  nscales = 3
  if (powheginput('#cs_scales') > 0) nscales = nint(powheginput('#cs_scales'))
  if (nscales /= 1 .and. nscales /= 3 .and. nscales /= 7) stop 'cs_scales must be 1, 3 or 7'

  ! inclusive part, at the order of the run
  incl_order_override = corder
  call inclusive_init()          ! parameters, PDFs, alpha_s, scales
  idum_job = idum                ! -iseed, for the exclusive part
  ! stem of the output files: always the seed (at least four digits), not
  ! the seedstr of the inclusive part, which is empty without a seeds file
  write(seedtxt,'(i0)') -idum_job
  if (len_trim(seedtxt) < 4) seedtxt = repeat('0', 4 - len_trim(seedtxt))//trim(seedtxt)
  runbase = 'provbfh-'//trim(seedtxt)//'-'//trim(oname(corder))
  ! the parts are also written on request, and for cs_4qint 2, whose
  ! sum is meaningless
  write_parts = nint(powheginput('#cs_write_parts')) == 1 .or. nint(powheginput('#cs_4qint')) == 2
  ! on-the-fly scale variations (the same points as the exclusive part)
  if (nscales > 1) then
     incl_nscale = nscales
     incl_scr(1:incl_nscale) = sc_r(1:incl_nscale)
     incl_scf(1:incl_nscale) = sc_f(1:incl_nscale)
  endif
  write(6,'(a,i2,a,i2,a)') ' proVBFH-cs: qcd_order ', corder, ', ', nscales, ' scale points (mu_R/mu_R0, mu_F/mu_F0), weights W1, W2, ...:'
  write(6,'(7(a,f4.2,a,f4.2,a))') (' (', sc_r(i), ',', sc_f(i), ')', i = 1, nscales)
  write(6,'(a,2i12)') ' proVBFH-cs: ncall1, ncall2 (both parts) = ', ncall1, ncall2
  ! the grid file of the inclusive part (grids<seed>.dat; readingrid) is
  ! its own, the exclusive part has grids-excl<seed>.dat
  readin_job = readin
  call run_inclusive(.true., sig_i, err_i)
  readin = readin_job
  call cpu_time(t1)
  write(6,'(a,es13.5,a,es10.3,a)') ' proVBFH-cs inclusive part: total (no cuts) ', sig_i, ' +- ', err_i, ' pb'
  write(6,'(a,f12.1,a)') ' proVBFH-cs inclusive part: CPU ', t1 - t0, ' s'
  call pwhgsetout
  call cs_hist_save()
  ! (at qcd_order 1 and with inclusive_only the sum is this part)
  if (write_parts .and. corder > 1 .and. nint(powheginput('#inclusive_only')) /= 1) then
     histname = trim(runbase)//'-incl'
     call pwhgtopout(histname)
  endif
  ! qcd_order 1, or inclusive_only 1 (as in the other proVBFH codes): the
  ! inclusive part alone, at the order qcd_order
  if (corder == 1 .or. nint(powheginput('#inclusive_only')) == 1) then
     call write_total()
     stop
  endif

  ! exclusive part: qcd_order 2 is the (1,0) + (0,1) of excl_order 1, qcd_order 3
  ! also the (2,0) + (0,2) and (1,1)
  if (corder == 2) excl_order = 1
  if (corder == 3) excl_order = 3
  call cs_init_vbfnlo()          ! VBFNLO couplings from vbfnlo.input
  call cs_excl_setup()
  call set_beams(sqrts)
  if (powheginput('#cs_npow') >= 0) excl_npow = nint(powheginput('#cs_npow'))
  if (powheginput('#cs_hardfrac') >= 0) excl_hardfrac = powheginput('#cs_hardfrac')
  four_hard = excl_hardfrac
  if (powheginput('#cs_hardfrac2') >= 0) four_hard2 = powheginput('#cs_hardfrac2')
  if (powheginput('#cs_hard2mode') >= 0) four_hard2mode = nint(powheginput('#cs_hard2mode'))
  if (powheginput('#cs_cutoff') > 0) excl_cutoff = powheginput('#cs_cutoff')
  ! cs_4qint 2: only the identical-quark interference of the four-quark reals
  ! (no dipoles: it is integrable) and its counterevents
  if (nint(powheginput('#cs_4qint')) == 2) then
     if (excl_order /= 3) stop 'cs_4qint 2 needs qcd_order 3'
     excl_intonly = .true.
     nlo2_intonly = .true.
     write(6,'(a)') ' proVBFH-cs: cs_4qint 2, exclusive part = the four-quark interference alone ('//trim(runbase)//'-excl-W<k>.top)'
  endif
  if (excl_order >= 2) then
     call cs_excl_setup2(.true.)   ! stage 2: real flavour groups, dipole tables
     call cpu_time(t1)
     write(6,'(a,f10.2,a)') ' proVBFH-cs stage-2 set-up: CPU ', t1 - t0, ' s'
  endif
  ! on-the-fly scale variations
  if (nscales > 0) then
     call cs_scales_setup(nscales)
     write(6,'(a,i2,a)') ' proVBFH-cs: ', excl_nscale, ' scale points (mu_R/mu_R0, mu_F/mu_F0), weights W1, W2, ...:'
     write(6,'(7(a,f4.2,a,f4.2,a))') (' (', sc_r(i), ',', sc_f(i), ')', i = 1, excl_nscale)
  endif
  ! the seven dimensions of the four-parton events are not adapted (their
  ! density enters the multichannel weight)
  ndim = 13
  if (excl_order >= 2) then
     ndim = 20
     jfreeze = 14
  endif
  ! small-x channel for the Born phase space (cs_smallx, cs_smallx_x0): five
  ! extra, non-adapted dimensions (selector, four numbers)
  if (powheginput('#cs_smallx') >= 0) sx_h = powheginput('#cs_smallx')
  if (sx_h > 0) then
     if (powheginput('#cs_smallx_x0') > 0) sx_x0 = powheginput('#cs_smallx_x0')
     sx_dim0 = ndim
     ndim = ndim + 5
     write(6,'(a,f6.3,a,f7.4,a,i3)') ' proVBFH-cs: small-x channel, fraction ', sx_h, ', x0 ', sx_x0, ', dimensions ', ndim
  endif
  idum = idum_job   ! random numbers as in a run of the exclusive part alone
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
     if (write_parts) then
        histname = trim(runbase)//'-excl'
        call pwhgtopout(histname)
     endif
     call cs_hist_add()
     call write_total()
  endif
  call cpu_time(t1)
  write(6,'(a,es13.5,a,es10.3,a)') ' proVBFH-cs exclusive part: sum |w1|+|w2| = ', integ, ' +- ', err, ' pb'
  write(6,'(a,4i14)') ' points, line radiations cut off, NaN, points failing all cuts: ', excl_stats
  if (excl_order >= 2) write(6,'(a,2i14)') ' stage 2: real and Born matrix-element calls: ', nlo2_ncount
  if (excl_order >= 2) write(6,'(a,i14)') ' stage 2: four-parton points dropped by the technical cut: ', nlo2_ncut
  if (excl_nscale > 1) write(6,'(a,2i14)') ' scale variations: matrix-element cache hits, misses: ', memo_hits, memo_misses
  write(6,'(a,f12.1,a)') ' proVBFH-cs exclusive part: CPU ', t1 - t0, ' s'

contains

  ! total histograms provbfh-<iseed>-<order>-W<k>.top and the cross section with
  ! the cuts of the analysis (histogram 1) for every scale point,
  ! provbfh-<iseed>-<order>.dat
  subroutine write_total()
    character(len=60) :: xsname
    integer :: k
    call pwhgtopout(runbase)
    xsname = trim(runbase)//'.dat'
    open(unit=11, file=trim(xsname), action='write')
    write(11,'(a)') '# proVBFH-cs: total '//trim(oname(corder))//' cross section with the cuts of the analysis'
    if (corder == 1 .or. nint(powheginput('#inclusive_only')) == 1) then
       write(11,'(a)') '# (inclusive part only), one line per scale point (mu_R/mu_R0, mu_F/mu_F0)'
    else
       write(11,'(a)') '# (inclusive + exclusive part), one line per scale point (mu_R/mu_R0, mu_F/mu_F0)'
    endif
    do k = 1, nscales
       call cs_hist_total(k, sig_k, err_k)
       write(11,'(es13.6,a,es13.6,a,2f5.2)') sig_k, ' +/- ', err_k, ' (pb)  # W', sc_r(k), sc_f(k)
       write(6,'(a,i2,a,es13.6,a,es13.6,a)') ' proVBFH-cs total, W', k, ': ', sig_k, ' +- ', err_k, ' pb'
    enddo
    close(11)
  end subroutine write_total
end program provbfh_cs
