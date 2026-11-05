!----------------------------------------------------------------------
! Pointwise test of the subtraction of the (2,0)/(0,2) real (stage 2):
! in a final-final collinear limit k_i || k_j the real of every group
! with that pair and its dipoles for the pair cancel,
!     M_real - D_ij -> 0   (relative to M_real; group_me returns the dipoles
!                           with the sign of the real, nlo2_real_me subtracts them),
! at every azimuth of the splitting (the spin-correlated kernels of the
! dipoles must match the real pointwise, not only after averaging).
!
! Points: Born kinematics from the inclusive phase space (as in
! cs_excl_setup2), real configurations from nlo2_real_kin (without the
! technical cut). For a pair (i,j) of a real configuration, map_ff gives the
! massless parents (ptij, ptk) of the dipole, and split_ff, the exact inverse
! of the map, splits ptij at y = s_ij/s_ijk from 1e-2 to 1e-8, azimuth phi,
! with z fixed: the total momentum is unchanged and s_ij = y s_ijk up to
! rounding (checked, to max(1d-6, 1d-12/y)).
! nlo2_group_test recomputes the six maps and clears the Born cache for the
! momenta it is given (group_me reads the maps from the module state that
! nlo2_real_kin sets). For each group whose real has a 1/y
! singularity in the pair (|M| grows by more than SING from y = 1e-6 to
! 1e-8; 1/y gives 100, the 1/sqrt(y) of the interference terms of the
! identical-quark groups, which needs no subtraction, gives 10) the residual
!     r(y) = |M - D| / |M|
! must be below TOL at the smallest y and smaller than at larger y.
! Corrections to the leading power are O(sqrt(y)), with coefficients up to
! ~30, hence the smallest y of 1e-8.
! The collinear singularity of a final pair is cancelled by all dipoles
! with that pair: the final-final one (spectator the third parton) and the
! final-initial one (spectator the incoming parton), slots ip and 3+ip of
! group_me. Not covered: final-initial collinear limits and soft limits (the
! four-quark flip terms have their own limit test, in
! validation/v4_4qint/limits_4qint.f90).
! Run from the example/ directory, which holds the cards it reads:
!   cd example && ../tests/test_limits
! Prints 'all N singular cases cancel' on success (the caller checks for
! this line: a stop in the library exits with status 0), and exits with
! status 1 if a check fails.
!----------------------------------------------------------------------
program test_limits
  use cs_nlo2, only: nlo2_real_kin, nlo2_group_test, nlo2_ngrp, nlo2_group_line
  use cs_exclusive, only: cs_excl_setup, cs_excl_setup2
  use cs_dipoles, only: map_ff, split_ff
  use cs_kinematics, only: mdot
  use incl_parameters, only: sqrts, set_parameters
  use phase_space, only: gen_phsp, set_phsp, set_beams
  implicit none
  integer, parameter :: dp = kind(1d0)
  real(dp), parameter :: tol = 1d-2, sing = 3d1, cutoff = 1d-5
  integer, parameter :: nborn = 3, nphi = 3, ny = 7, nreal_try = 50
  real(dp), parameter :: ys(ny) = [1d-2, 1d-3, 1d-4, 1d-5, 1d-6, 1d-7, 1d-8]
  real(dp), parameter :: phis(nphi) = [0.3_dp, 1.9_dp, 4.2_dp]
  real(dp), parameter :: zsplit = 0.37_dp
  integer :: ipair(2,3), lpair(3)
  real(dp) :: pb(0:3,5), xb(2), xr(7), r(7), pr(0:3,7), pd(0:3,6,6)
  real(dp) :: xb1, xb2, jacb, sbeams, y0, z0
  real(dp) :: kk(0:3,3), kn(0:3,3,ny), pa(0:3), mt(0:3), mk(0:3)
  real(dp) :: pi_(0:3), pj(0:3), pk(0:3)
  real(dp) :: mev(ny), dv(ny), rat(ny), sij, sref
  real(dp) :: dtmp(6), mtmp, worst_r, worst_split
  integer :: ib, it, line, ig, ip, ii, jj, ll, ph, iy, itry, ntest, nsing, nfail, npass
  integer :: nbadmom, nbadsplit
  logical :: dok(6), ok, failed
  external cs_get_born, cs_init_vbfnlo, hoppetStart, hoppetSetCoupling
  integer :: seed(8)

  ! reproducible random numbers (the real configurations are drawn at random)
  seed = 20261105
  call random_seed(put=seed)
  ! the parameters of the inclusive part (sqrts, Higgs mass, ...), without its PDFs
  call set_parameters()
  ! the running coupling of hoppet (as in validation/v4_4qint/limits_4qint.f90)
  call hoppetStart(0.5d0, 2)
  call hoppetSetCoupling(0.118d0, 91.1876d0, 2)
  call cs_init_vbfnlo()          ! VBFNLO couplings from vbfnlo.input (as cs_main)
  ipair = reshape([1,2, 1,3, 2,3], [2,3])
  lpair = [3, 2, 1]
  call cs_excl_setup()
  call set_beams(sqrts)
  call cs_excl_setup2(.false.)

  ntest = 0; nsing = 0; nfail = 0; npass = 0; nbadmom = 0; nbadsplit = 0
  worst_r = 0; worst_split = 0
  write(*,'(a)') 'test_limits: final-final collinear limits of the (2,0)/(0,2) real'
  write(*,'(a)') '  line  pair  phi  group    |M|(1e-2)      |M|(1e-8)      r(1e-2)      r(1e-4)      r(1e-8)   result'

  ! Born points, from the same deterministic sequence and cuts as cs_excl_setup2
  ! (random Born points can be degenerate, where line_radiation stops)
  ib = 0
  do it = 1, 1000
     xr = modulo(0.1234567_dp*it*[1.0_dp, 1.7_dp, 2.3_dp, 3.1_dp, 3.7_dp, 4.3_dp, 5.9_dp], 1.0_dp)
     call gen_phsp(xr)
     call set_phsp()
     call cs_get_born(pb, xb1, xb2, jacb, sbeams)
     if (jacb == 0) cycle
     if (max(xb1, xb2) > 0.5_dp) cycle
     xb = [xb1, xb2]
     ib = ib + 1
     if (ib > nborn) exit
     do line = 1, 2
        ! a real configuration of this line, from random numbers
        ok = .false.
        do itry = 1, nreal_try
           call random_number(r)
           call nlo2_real_kin(line, pb, xb, r, cutoff, pr, pd, dok, ok, nocut=.true.)
           if (ok) exit
        enddo
        if (.not. ok) cycle
        pa = pr(:, line)
        kk = pr(:, 4:6)
        do ip = 1, 3
           ii = ipair(1, ip); jj = ipair(2, ip); ll = lpair(ip)
           ! the massless parents of the dipole of this pair
           call map_ff(kk(:,ii), kk(:,jj), kk(:,ll), mt, mk, y0, z0)
           sref = 2*mdot(mt, mk)
           do ph = 1, nphi
              ! the splittings of the parents at the ny values of y
              do iy = 1, ny
                 call split_ff(mt, mk, ys(iy), zsplit, phis(ph), pi_, pj, pk)
                 kn(:,:,iy) = kk
                 kn(:,ii,iy) = pi_; kn(:,jj,iy) = pj; kn(:,ll,iy) = pk
                 ! s_ij = y s_ijk exactly, and the total momentum is unchanged
                 sij = 2*mdot(pi_, pj)
                 worst_split = max(worst_split, abs(sij/(ys(iy)*sref) - 1))
                 ! (to 1d-6 down to y = 1d-6, as before; below, rounding limits s_ij to an
                 ! absolute accuracy of ~1d-12 s_ijk: observed 2e-5 at y = 1e-8)
                 if (abs(sij/(ys(iy)*sref) - 1) > max(1d-6, 1d-12/ys(iy))) &
                      nbadsplit = nbadsplit + 1
                 if (maxval(abs(kn(:,ii,iy) + kn(:,jj,iy) + kn(:,ll,iy) - kk(:,ii) - kk(:,jj) - kk(:,ll))) &
                      > 1d-8*sqrts) nbadmom = nbadmom + 1
              enddo
              ! each group of this line: its real and dipole at the ny points
              do ig = 1, nlo2_ngrp()
                 if (nlo2_group_line(ig) /= line) cycle
                 do iy = 1, ny
                    call nlo2_group_test(ig, pa, kn(:,:,iy), mtmp, dtmp)
                    mev(iy) = mtmp
                    ! the dipoles singular in i || j: final-final with the third parton as
                    ! spectator (slot ip) and final-initial with the incoming parton (slot 3+ip;
                    ! that slot also holds the initial-final dipoles of the map, regular here)
                    dv(iy) = dtmp(ip) + dtmp(3+ip)
                    rat(iy) = abs(mtmp - dv(iy)) / max(abs(mtmp), tiny(1d0))
                 enddo
                 ntest = ntest + 1
                 ! no 1/y singularity in this pair: nothing to check
                 if (abs(mev(ny)) <= sing * max(abs(mev(ny-2)), tiny(1d0))) cycle
                 nsing = nsing + 1
                 failed = .false.
                 if (.not. (rat(ny) <= tol)) failed = .true.          ! also catches NaN
                 ! must decrease towards the limit (against the largest residual at
                 ! y >= 1e-4: the smaller y can be at the rounding floor)
                 if (.not. (rat(ny) < maxval(rat(1:3)))) failed = .true.
                 if (dv(ny) == 0d0) failed = .true.                   ! a singular real needs its dipole
                 if (failed) then
                    nfail = nfail + 1
                 else
                    npass = npass + 1
                 endif
                 worst_r = max(worst_r, rat(ny))
                 write(*,'(2x,i3,2x,i1,a1,i1,2x,f5.2,2x,i5,2x,es13.4,2x,es13.4,2x,es11.3,2x,es11.3,2x,es11.3,2x,a)') &
                      line, ii, '-', jj, phis(ph), ig, mev(1), mev(ny), rat(1), rat(3), rat(ny), &
                      merge('FAIL', 'ok  ', failed)
              enddo
           enddo
        enddo
     enddo
  enddo

  write(*,'(a,i0,a,i0,a)') ' singular (group, pair, azimuth) cases: ', nsing, ' of ', ntest, ' group evaluations'
  write(*,'(a,es10.3,a,es10.3)') ' worst residual at y = 1e-8: ', worst_r, ' (tolerance ', tol
  write(*,'(a,es10.3)') ' worst |s_ij/(y s_ijk) - 1| of the splits: ', worst_split
  write(*,'(a,i0,a,i0)') ' split failures (invariant mass, momentum): ', nbadsplit, ', ', nbadmom
  if (nsing == 0) then
     write(*,'(a)') 'test_limits: FAILED, no singular case found (the test is not testing anything)'
     stop 1
  endif
  if (nfail > 0 .or. nbadsplit > 0 .or. nbadmom > 0) then
     write(*,'(a,i0,a)') 'test_limits: FAILED with ', nfail, ' cases'
     stop 1
  endif
  write(*,'(a,i0,a)') 'test_limits: all ', npass, ' singular cases cancel'
end program test_limits
