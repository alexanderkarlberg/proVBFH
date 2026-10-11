!----------------------------------------------------------------------
! Unit test of cs_smallx (the small-x channel of the Born phase space).
!   1. rho_B (sx_rho) integrates to 1 over the unit cube b = (r2..r6);
!   2. the points of sx_gen follow sx_rho: E_gen[f] = int f rho_B db
!      (against flat sampling of b), for several test functions;
!   3. the mixture rho = (1-h) g_A + h rho_B with a stand-in for the
!      VEGAS density g_A (a product of piecewise-constant densities with a
!      hole of 1/25 below r = 0.4, like the grids of r2, r3, r6): the
!      estimator F/rho of sx_choose/sx_gen/g_A points agrees with the flat
!      integral of F, for several F, including ones peaked in the hole;
!   4. negative control: with a wrong rho_B (the Jacobian without the
!      factor r2) test 3 must fail.
! lntmin = ln(mH^2/s) at 13 TeV, x0 = 0.03 and 0.1. Exit status 1 if a
! check fails or the negative control does not.
!----------------------------------------------------------------------
program test_smallx
  use cs_smallx
  implicit none
  integer, parameter :: dp = kind(1d0)
  integer, parameter :: nf = 5
  real(dp) :: lnt, x0, h, pull(nf), worst
  integer :: ix0, nfail, nneg
  logical :: wrong
  call random_seed()
  lnt = log((125.0_dp/13000.0_dp)**2)
  h = 0.3_dp
  nfail = 0
  do ix0 = 1, 2
     x0 = merge(0.03_dp, 0.1_dp, ix0 == 1)
     write(*,'(a,f6.3)') ' x0 =', x0
     call check_norm(lnt, x0, nfail)
     call check_gen(lnt, x0, nfail)
     wrong = .false.
     call check_mix(lnt, x0, h, wrong, pull)
     worst = maxval(abs(pull))
     write(*,'(a,es10.2)') '   mixture, worst pull: ', worst
     if (worst > 4) nfail = nfail + 1
  enddo
  ! negative control
  wrong = .true.
  call check_mix(lnt, 0.03_dp, h, wrong, pull)
  nneg = count(abs(pull) > 4)
  write(*,'(a,i2,a,es10.2)') ' negative control (wrong Jacobian): ', nneg, ' of 5 test functions fail, worst pull ', maxval(abs(pull))
  if (nneg == 0) then
     write(*,*) 'test_smallx: negative control did not fail'
     nfail = nfail + 1
  endif
  if (nfail > 0) then
     write(*,*) 'test_smallx: FAILED'
     call exit(1)
  endif
  write(*,*) 'test_smallx: passed'
contains
  real(dp) function ranu()
    call random_number(ranu)
  end function ranu

  ! test functions of b
  real(dp) function tf(i, b, lnt)
    integer, intent(in) :: i
    real(dp), intent(in) :: b(5), lnt
    real(dp) :: lt, y
    lt = lnt*(1 - b(1)**2); y = -(1 - 2*b(2))*lt/2
    select case (i)
    case (1); tf = 1
    case (2); tf = b(1)
    case (3); tf = b(2)**2 + 0.5_dp*b(5)
    case (4); tf = exp(-((b(1) - 0.2_dp)/0.1_dp)**2)*merge(3.0_dp, 1.0_dp, b(2) < 0.1_dp)
    case default; tf = merge(1.0_dp, 0.0_dp, 0.5_dp*lt + y < log(0.01_dp)) + 0.3_dp*b(3)
    end select
  end function tf

  subroutine check_norm(lnt, x0, nfail)
    real(dp), intent(in) :: lnt, x0
    integer, intent(inout) :: nfail
    integer, parameter :: n = 4000000
    real(dp) :: b(5), s, s2, v, m, e
    integer :: i
    s = 0; s2 = 0
    do i = 1, n
       call random_number(b)
       v = sx_rho(b, lnt, x0)
       s = s + v; s2 = s2 + v*v
    enddo
    m = s/n; e = sqrt((s2/n - m*m)/n)
    write(*,'(a,f8.4,a,f8.4)') '   int rho_B db = ', m, ' +-', e
    if (abs(m - 1) > 4*e) nfail = nfail + 1
  end subroutine check_norm

  subroutine check_gen(lnt, x0, nfail)
    real(dp), intent(in) :: lnt, x0
    integer, intent(inout) :: nfail
    integer, parameter :: n = 4000000
    real(dp) :: b(5), u(5), s1(nf), s2(nf), q1(nf), q2(nf), f, w, pl, e1, e2
    integer :: i, k, j
    s1 = 0; s2 = 0; q1 = 0; q2 = 0
    do i = 1, n
       call random_number(u)
       k = merge(1, 2, ranu() < 0.5_dp)
       call sx_gen(k, u, lnt, x0, b)
       do j = 1, nf
          f = tf(j, b, lnt); s1(j) = s1(j) + f; q1(j) = q1(j) + f*f
       enddo
       call random_number(b)
       w = sx_rho(b, lnt, x0)
       do j = 1, nf
          f = tf(j, b, lnt)*w; s2(j) = s2(j) + f; q2(j) = q2(j) + f*f
       enddo
    enddo
    pl = 0
    do j = 1, nf
       s1(j) = s1(j)/n; s2(j) = s2(j)/n
       e1 = sqrt(max(q1(j)/n - s1(j)**2, 0.0_dp)/n); e2 = sqrt(max(q2(j)/n - s2(j)**2, 0.0_dp)/n)
       pl = max(pl, abs(s1(j) - s2(j))/sqrt(e1**2 + e2**2))
    enddo
    write(*,'(a,es10.2)') '   generated points vs rho_B, worst pull: ', pl
    if (pl > 4) nfail = nfail + 1
  end subroutine check_gen

  ! density and inverse-CDF sample of the stand-in grid: g = 1/25 below 0.4
  ! dimensions 1, 2 (r2, r3): the hole; 3 (r4): a mild one; 4, 5: flat
  real(dp) function gA1(j, r)
    integer, intent(in) :: j
    real(dp), intent(in) :: r
    select case (j)
    case (1, 2); gA1 = merge(0.04_dp, (1 - 0.016_dp)/0.6_dp, r < 0.4_dp)
    case (3); gA1 = merge(0.25_dp, (1 - 0.075_dp)/0.7_dp, r < 0.3_dp)
    case default; gA1 = 1
    end select
  end function gA1
  real(dp) function gAinv(j, u)
    integer, intent(in) :: j
    real(dp), intent(in) :: u
    select case (j)
    case (1, 2)
       if (u < 0.016_dp) then
          gAinv = u/0.04_dp
       else
          gAinv = 0.4_dp + (u - 0.016_dp)/((1 - 0.016_dp)/0.6_dp)
       endif
    case (3)
       if (u < 0.075_dp) then
          gAinv = u/0.25_dp
       else
          gAinv = 0.3_dp + (u - 0.075_dp)/((1 - 0.075_dp)/0.7_dp)
       endif
    case default; gAinv = u
    end select
  end function gAinv

  subroutine check_mix(lnt, x0, h, wrong, pull)
    real(dp), intent(in) :: lnt, x0, h
    logical, intent(in) :: wrong
    real(dp), intent(out) :: pull(nf)
    integer, parameter :: n = 6000000
    real(dp) :: b(5), u(5), s, rho, gA, rb, f, s1(nf), q1(nf), s2(nf), q2(nf), e1, e2
    integer :: i, j, ch, jj
    s1 = 0; q1 = 0; s2 = 0; q2 = 0
    do i = 1, n
       call sx_choose(ranu(), h, ch, u(1))
       call random_number(u(2:5))
       if (ch == 0) then
          do jj = 1, 5
             b(jj) = gAinv(jj, ranu())
          enddo
       else
          call sx_gen(ch, u, lnt, x0, b)
       endif
       gA = 1
       do jj = 1, 5
          gA = gA*gA1(jj, b(jj))
       enddo
       rb = sx_rho(b, lnt, x0)
       if (wrong) rb = rb/(2*max(b(1), 1d-12))
       rho = (1 - h)*gA + h*rb
       do j = 1, nf
          f = tf(j, b, lnt)/rho; s1(j) = s1(j) + f; q1(j) = q1(j) + f*f
       enddo
       call random_number(b)
       do j = 1, nf
          f = tf(j, b, lnt); s2(j) = s2(j) + f; q2(j) = q2(j) + f*f
       enddo
    enddo
    do j = 1, nf
       s1(j) = s1(j)/n; s2(j) = s2(j)/n
       e1 = sqrt(max(q1(j)/n - s1(j)**2, 0.0_dp)/n); e2 = sqrt(max(q2(j)/n - s2(j)**2, 0.0_dp)/n)
       pull(j) = (s1(j) - s2(j))/sqrt(e1**2 + e2**2)
       write(*,'(a,i2,a,es12.5,a,es9.2,a,es12.5,a,es9.2,a,f7.2)') '   F', j, ': mixture ', s1(j), ' +-', e1, '  flat ', s2(j), ' +-', e2, '  pull ', pull(j)
    enddo
  end subroutine check_mix
end program test_smallx
