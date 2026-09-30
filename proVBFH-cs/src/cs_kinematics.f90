!----------------------------------------------------------------------
! Three-parton kinematics for one VBF quark line with the momentum
! transfer kept fixed (DIS-like, as DISENT's GENTHR), in covariant form.
!
! Input: the line's Born momenta, pB (incoming, = xB * P_beam) and pOB
! (outgoing), so that q = pOB - pB and Q^2 = -q^2 = 2 pB.pOB, and three
! random numbers r(1:3). Output: the incoming parton pin = pB/xp (along
! the same beam, momentum fraction xB/xp) and two outgoing partons a, b
! with pin + q = a + b... i.e. pin - (pB - pOB) = a + b, parametrised by
!   xp = Q^2/(2 pin.q_in) in [xB, 1],  z = pin.a/pin.q_in in (0,1),
!   phi = azimuth of a around the (pB, pOB) plane,
! where q_in = pOB - pB. Explicitly
!   a = alpha pB + z pOB + pT (cos(phi) e1 + sin(phi) e2),
!   alpha = (1-xp)(1-z)/xp,  pT^2 = Q^2 z (1-xp)(1-z)/xp,
! with e1, e2 unit vectors orthogonal to pB and pOB.
!
! Weight: relative to the Born line factor 2 pi f(xB)/Q^2, the line's
! phase space and flux are
!   Q^2/(16 pi^2) dxp/xp f(xB/xp) dz dphi/(2 pi),
! so that (docs/DESIGN.md) the (1,0) weight is
!   kn_jacborn/(2 x1 x2 S) f(xi1) f(x2) |M(H+3j), line 1|^2 * wrad,
!   wrad = Q^2/(16 pi^2) (1/xp) dxp/dr1 dz/dr2
! (dphi/(2 pi) = dr3).
! Sampling (npow > 0): 1 - xp = (1 - xB) r1^npow;  z = (2 r2)^npow/2 for
! r2 < 1/2, 1 - z = (2 (1-r2))^npow/2 otherwise. npow = 0: logarithmic in
! 1 - xp and min(z, 1-z) down to the cutoff (the default in proVBFH-cs).
! With npow = 0 and hard = h > 0 (optional), a second channel with
! probability h samples the hard region: ln(xp) uniform in [ln xB, 0] and
! z uniform; the weight uses the combined density
!   g = (1-h) g_log(xp) g_log(z) + h g_hard(xp) g_hard(z),
! g_log(xp) = 1/((1-xp) ln((1-xB)/c)), g_log(z) = 1/(2 min(z,1-z) ln(1/(2c))),
! g_hard(xp) = 1/(xp ln(1/xB)), g_hard(z) = 1. The logarithmic channel puts
! few points at small xp, where both partons of the line are hard
! (pT^2 = Q^2 z (1-z) (1-xp)/xp), which dominates the high-pT tails.
! Events with 1-xp, z or 1-z below cutoff are rejected (ok = .false.), as
! DISENT's invariant cutoff.
! Momenta are (E, px, py, pz) with index 0:3.
!----------------------------------------------------------------------
module cs_kinematics
  implicit none
  private
  public :: line_radiation, mdot

  integer, parameter :: dp = kind(1d0)
  real(dp), parameter :: pi = 3.141592653589793238462643383279502884197_dp

contains

  real(dp) function mdot(p, q)
    real(dp), intent(in) :: p(0:3), q(0:3)
    mdot = p(0)*q(0) - p(1)*q(1) - p(2)*q(2) - p(3)*q(3)
  end function mdot

  subroutine line_radiation(pB, pOB, xB, r, npow, cutoff, pin, a, b, xp, z, wrad, ok, hard)
    real(dp), intent(in)  :: pB(0:3), pOB(0:3), xB, r(3), cutoff
    integer,  intent(in)  :: npow
    real(dp), intent(out) :: pin(0:3), a(0:3), b(0:3), xp, z, wrad
    logical,  intent(out) :: ok
    real(dp), intent(in), optional :: hard
    real(dp) :: Q2, jx, jz, u, alpha, pT, phi, e1(0:3), e2(0:3), qin(0:3)
    real(dp) :: h, s, lz, lx
    ok = .false.
    pin = 0; a = 0; b = 0; xp = 1; z = 0; wrad = 0
    Q2 = 2*mdot(pB, pOB)
    if (Q2 <= 0 .or. xB >= 1) return
    h = 0
    if (present(hard)) h = hard
    if (npow == 0 .and. h > 0) then
       ! logarithmic channel (probability 1-h) and hard channel (h)
       if (1 - xB <= cutoff) return
       u = log((1 - xB)/cutoff)
       lz = log(0.5_dp/cutoff)
       lx = log(1/xB)
       if (r(1) < 1 - h) then
          s = r(1)/(1 - h)
          xp = 1 - cutoff*exp(u*s)
          if (r(2) < 0.5_dp) then
             z = cutoff*exp(lz*2*r(2))
          else
             z = 1 - cutoff*exp(lz*2*(1 - r(2)))
          endif
       else
          s = (r(1) - (1 - h))/h
          xp = xB**(1 - s)
          z = r(2)
       endif
       if (1 - xp < cutoff .or. z < cutoff .or. 1 - z < cutoff) return
       ! jx*jz = 1/g(xp, z)
       jx = 1/((1 - h)/((1 - xp)*u)/(2*min(z, 1 - z)*lz) + h/(xp*lx))
       jz = 1
    elseif (npow > 0) then
       ! power sampling: xp in [xB, 1], z symmetric about 1/2
       xp = 1 - (1 - xB)*r(1)**npow
       jx = (1 - xB)*npow*r(1)**(npow-1)
       if (r(2) < 0.5_dp) then
          u = 2*r(2)
          z = 0.5_dp*u**npow
       else
          u = 2*(1 - r(2))
          z = 1 - 0.5_dp*u**npow
       endif
       jz = npow*u**(npow-1)
    else
       ! logarithmic sampling of 1-xp in [cutoff, 1-xB] and of min(z,1-z)
       ! in [cutoff, 1/2]: flattens the 1/((1-xp) z (1-z)) behaviour of the
       ! matrix element exactly. (With power sampling the weights grow like
       ! 1/sqrt(1-xp) and dominate the variance at large Q, where emissions
       ! that are hard for the jet cuts have 1-xp ~ W^2/Q^2 << 1.)
       if (1 - xB <= cutoff) return
       u = log((1 - xB)/cutoff)
       xp = 1 - cutoff*exp(u*r(1))
       jx = (1 - xp)*u
       if (r(2) < 0.5_dp) then
          z = cutoff*exp(log(0.5_dp/cutoff)*2*r(2))
          jz = 2*z*log(0.5_dp/cutoff)
       else
          z = 1 - cutoff*exp(log(0.5_dp/cutoff)*2*(1 - r(2)))
          jz = 2*(1 - z)*log(0.5_dp/cutoff)
       endif
    endif
    if (1 - xp < cutoff .or. z < cutoff .or. 1 - z < cutoff) return
    phi = 2*pi*r(3)
    call transverse_basis(pB, pOB, e1, e2)
    alpha = (1 - xp)*(1 - z)/xp
    pT = sqrt(Q2*z*(1 - xp)*(1 - z)/xp)
    pin = pB/xp
    a = alpha*pB + z*pOB + pT*(cos(phi)*e1 + sin(phi)*e2)
    qin = pOB - pB
    b = pin + qin - a
    wrad = Q2/(16*pi**2)/xp*jx*jz
    ok = .true.
  end subroutine line_radiation

  ! unit spacelike vectors e1, e2 (e.e = -1) orthogonal to the massless
  ! pB, pOB and to each other, from the lab x (or y) and y (or z) axes
  subroutine transverse_basis(pB, pOB, e1, e2)
    real(dp), intent(in)  :: pB(0:3), pOB(0:3)
    real(dp), intent(out) :: e1(0:3), e2(0:3)
    real(dp) :: ref(0:3,3), n2
    integer :: i, got
    ref = 0
    ref(1,1) = 1; ref(2,2) = 1; ref(3,3) = 1
    got = 0
    do i = 1, 3
       if (got == 0) then
          e1 = perp(ref(:,i))
          n2 = -mdot(e1, e1)
          if (n2 > 1d-6*maxval(abs(ref(:,i)))**2) then
             e1 = e1/sqrt(n2)
             got = i
          endif
       elseif (got > 0) then
          e2 = perp(ref(:,i))
          e2 = e2 + mdot(e2, e1)*e1        ! e1.e1 = -1
          n2 = -mdot(e2, e2)
          if (n2 > 1d-6) then
             e2 = e2/sqrt(n2)
             return
          endif
       endif
    enddo
    stop 'cs_kinematics: no transverse basis'
  contains
    function perp(v) result(w)
      real(dp), intent(in) :: v(0:3)
      real(dp) :: w(0:3), d
      d = mdot(pB, pOB)
      w = v - mdot(v, pOB)/d*pB - mdot(v, pB)/d*pOB
    end function perp
  end subroutine transverse_basis

end module cs_kinematics
