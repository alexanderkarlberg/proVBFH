!----------------------------------------------------------------------
! Small-x channel for the Born phase space of the exclusive part
! (cs_smallx, 11 Oct 2026; docs in notes/2026-10-cs-production/README.md).
!
! Why: the VEGAS grids of the Born variables r2 (tau = x1 x2), r3 (y) and
! r6 (the dijet angle) fall to 1/25 of uniform below the Born threshold of
! the VBF cuts (tau < 3e-3) and in the corners of y, while hard real
! emissions off a line with a small Bjorken x do pass the cuts there. The
! rare points that land in these holes carry weights up to 1e3-1e4 and make
! sigma(>=3 j) heavy-tailed.
!
! Channels (an exact mixture in the unit cube b = (r2, r3, r4, r5, r6) of
! the five Born variables of gen_phsp, mapped to the kinematics as before):
!   A    (probability 1-h): b from the VEGAS variables, density g_A(b) = the
!        product of the current VEGAS grid densities,
!   B_k  (probability h/2 each, k = 1, 2 the beam with the small x): (ln x_k,
!        ln x_o) uniform in R_k = {x_k < x0, x1 x2 >= taumin} of the plane
!        (ln x1, ln x2), r4, r5, r6 uniform; the variables come from extra,
!        non-adapted dimensions. Both beams are treated alike.
! The density of the mixture at any b is
!   rho(b) = (1-h) g_A(b) + h rho_B(b),
!   rho_B(b) = sum_k (1/2) 1_{R_k}/Area * |d(ln x1, ln x2)/d(r2, r3)|,
!   Area = (ln x0 - ln taumin)^2/2,
!   |d(ln x1, ln x2)/d(r2, r3)| = 2 r2 |ln taumin| |ln tau|
! (tau = exp(ln taumin (1 - r2^2)), y = -(1 - 2 r3) ln(tau)/2 in gen_phsp),
! and the integrand is multiplied by g_D/rho, g_D = the VEGAS density of the
! VEGAS variables r2..r6 (the product cancels the VEGAS weight, so channel
! B points do not depend on their unused VEGAS draws).
!----------------------------------------------------------------------
module cs_smallx
  implicit none
  private
  public :: sx_h, sx_x0, sx_dim0, sx_choose, sx_gen, sx_rho

  integer, parameter :: dp = kind(1d0)
  ! fraction of the points in channels B (default 0.3; 0: off), x0, and the number of
  ! the VEGAS dimensions before the extra ones (selector and four numbers)
  real(dp), save :: sx_h = 0.3_dp, sx_x0 = 0.03_dp
  integer, save :: sx_dim0 = 20

contains

  ! channel from the selector s in [0,1): ch = 0 (A), 1 or 2 (B_k); u: s
  ! rescaled to a uniform number within the channel
  subroutine sx_choose(s, h, ch, u)
    real(dp), intent(in) :: s, h
    integer, intent(out) :: ch
    real(dp), intent(out) :: u
    if (s < 0.5_dp*h) then
       ch = 1; u = s/(0.5_dp*h)
    elseif (s < h) then
       ch = 2; u = (s - 0.5_dp*h)/(0.5_dp*h)
    else
       ch = 0; u = 0
    endif
  end subroutine sx_choose

  ! channel B_k: b = (r2, r3, r4, r5, r6) from five uniform numbers
  subroutine sx_gen(k, u, lntmin, x0, b)
    integer, intent(in) :: k
    real(dp), intent(in) :: u(5), lntmin, x0
    real(dp), intent(out) :: b(5)
    real(dp) :: L, uk, v, l1, l2, lt, y
    L = log(x0) - lntmin
    uk = lntmin + L*sqrt(u(1))              ! ln x_k, density ~ (ln x_k - ln taumin)
    v = (lntmin - uk)*(1 - u(2))            ! ln x_o uniform in [ln taumin - ln x_k, 0]
    if (k == 1) then
       l1 = uk; l2 = v
    else
       l1 = v; l2 = uk
    endif
    lt = l1 + l2
    y = 0.5_dp*(l1 - l2)
    b(1) = sqrt(max(1 - lt/lntmin, 0.0_dp))
    b(2) = 0.5_dp
    if (lt < 0) b(2) = 0.5_dp + y/lt
    b(3:5) = u(3:5)
  end subroutine sx_gen

  ! density of channel B (both beams) at b, with respect to d r2 ... d r6
  real(dp) function sx_rho(b, lntmin, x0) result(rho)
    real(dp), intent(in) :: b(5), lntmin, x0
    real(dp) :: lt, y, l1, l2, L, area, lx0
    lx0 = log(x0)
    L = lx0 - lntmin
    rho = 0
    if (L <= 0) return
    area = 0.5_dp*L*L
    lt = lntmin*(1 - b(1)**2)
    y = -(1 - 2*b(2))*lt/2
    l1 = 0.5_dp*lt + y
    l2 = 0.5_dp*lt - y
    if (l1 < lx0) rho = rho + 0.5_dp/area
    if (l2 < lx0) rho = rho + 0.5_dp/area
    rho = rho*2*b(1)*abs(lntmin)*abs(lt)
  end function sx_rho

end module cs_smallx
