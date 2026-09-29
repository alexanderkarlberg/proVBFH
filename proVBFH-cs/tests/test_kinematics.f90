!----------------------------------------------------------------------
! Unit test of cs_kinematics (line_radiation), at random VBF-like line
! kinematics:
!   - a, b massless, pin along the beam with pin = pB/xp,
!   - momentum conservation pin + (pOB - pB) = a + b, i.e. the momentum
!     transfer q = pOB - pB (and Q^2) is unchanged,
!   - xp = Q^2/(2 pin.qin) and z = pin.a/pin.qin are reproduced,
!   - the sampled phase-space volume, sum of wrad over the unit cube,
!     equals Q^2/(16 pi^2) ln(1/xB) (the 2-body phase space of the
!     line's final state with f = |M|^2 = 1, cutoff -> 0).
! Both samplings (npow = 2 and logarithmic, npow = 0) are tested; with the
! logarithmic one the volume is that inside the cutoff.
! Exits with status 1 if a check fails.
!----------------------------------------------------------------------
program test_kinematics
  use cs_kinematics
  implicit none
  integer, parameter :: dp = kind(1d0)
  real(dp), parameter :: pi = 3.141592653589793238462643383279502884197_dp
  real(dp) :: pB(0:3), pOB(0:3), pin(0:3), a(0:3), b(0:3), qin(0:3), r(3)
  real(dp) :: xB, xp, z, wrad, Q2, E, th, ph, errmax(5), vol, vol2, exact, sig
  logical :: ok, fail
  integer :: ipt, k, n, ncut, npw, ipw
  errmax = 0; fail = .false.
  call random_seed()
  do ipw = 1, 2
  npw = merge(2, 0, ipw == 1)
  do ipt = 1, 20
     ! incoming Born parton along +z (or -z), outgoing at a random angle
     call random_number(r)
     xB = 0.01_dp + 0.8_dp*r(1)
     E = 1000*xB
     pB = [E, 0.0_dp, 0.0_dp, merge(E, -E, ipt <= 10)]
     call random_number(r)
     th = pi*r(1); ph = 2*pi*r(2)
     pOB = 300*r(3)*[1.0_dp, sin(th)*cos(ph), sin(th)*sin(ph), cos(th)] + 1d-3
     pOB(0) = sqrt(sum(pOB(1:3)**2))
     Q2 = 2*mdot(pB, pOB)
     qin = pOB - pB
     ! pointwise checks
     do k = 1, 200
        call random_number(r)
        call line_radiation(pB, pOB, xB, r, npw, 1d-12, pin, a, b, xp, z, wrad, ok)
        if (.not. ok) cycle
        errmax(1) = max(errmax(1), abs(mdot(a,a))/Q2, abs(mdot(b,b))/Q2)
        errmax(2) = max(errmax(2), maxval(abs(pin + qin - a - b))/sqrt(Q2))
        errmax(3) = max(errmax(3), maxval(abs(pin*xp - pB))/pB(0))
        errmax(4) = max(errmax(4), abs(Q2/(2*mdot(pin, qin)) - xp))
        errmax(5) = max(errmax(5), abs(mdot(pin, a)/mdot(pin, qin) - z))
     enddo
     ! phase-space volume (plain MC over the unit cube)
     n = 400000; vol = 0; vol2 = 0; ncut = 0
     do k = 1, n
        call random_number(r)
        call line_radiation(pB, pOB, xB, r, npw, 1d-12, pin, a, b, xp, z, wrad, ok)
        if (.not. ok) then
           ncut = ncut + 1; cycle
        endif
        vol = vol + wrad; vol2 = vol2 + wrad**2
     enddo
     vol = vol/n; sig = sqrt((vol2/n - vol**2)/n)
     exact = Q2/(16*pi**2)*log(1/xB)
     if (npw == 0) exact = volcut(xB, 1d-12)*Q2/(16*pi**2)
     if (abs(vol - exact) > 4*sig + 1d-6*exact) then
        fail = .true.
        write(*,'(a,i3,a,4es14.6)') ' FAIL volume, point', ipt, ': MC, error, exact, pull ', &
             & vol, sig, exact, (vol-exact)/sig
     endif
  enddo
  enddo
  write(*,'(a,5es10.2)') ' max errors (masses, momentum, pin, xp, z): ', errmax
  if (maxval(errmax) > 1d-10) fail = .true.
  if (fail) then
     write(*,*) 'test_kinematics: FAILED'
     call exit(1)
  endif
  write(*,*) 'test_kinematics: passed (20 line kinematics, pointwise and phase-space volume)'
contains
  ! integral of dxp/xp dz over 1-xp in [c, 1-xB], z in [c, 1-c]
  real(dp) function volcut(xB, c)
    real(dp), intent(in) :: xB, c
    volcut = (log((1-c)/xB))*(1 - 2*c)
  end function volcut
end program test_kinematics
