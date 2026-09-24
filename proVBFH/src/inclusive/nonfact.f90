module nonfact
  use helper
  use constants
  use incl_parameters
!  use grids
  implicit none

  public f1
!  public f2
  public chinf

  private f1_integrand
!  private f2_integrand
  private adaptive_integral, gauss_kronrod_15, scalar_integral, tri_integrand

  real(dp),save :: q1sq_rk,q2sq_rk,qHsq_rk
  real(dp),save :: p1x_rk, p2x_rk, p2y_rk, p3x_rk, p3y_rk
  real(dp), save :: MVsq ! Mass of vector boson, ie either MW or MZ
contains

  function f1(q1sq,q2sq,qHsq,MV) result(res)
    real(dp), intent(in) :: q1sq,q2sq,qHsq,MV
    real(dp) :: res

    ! Save input for runge kutta routine
    q1sq_rk = q1sq
    q2sq_rk = q2sq
    qHsq_rk = qHsq
    MVsq    = MV**2

    res = scalar_integral(f1_integrand,zero,one)
  end function f1

  function f1_analytic(q1sq,q2sq,qHsq,MV) result(res)
    real(dp), intent(in) :: q1sq,q2sq,qHsq,MV
    real(dp) :: res, delta1, delta2

    MVsq    = MV**2
    delta1 = q1sq + MVsq
    delta2 = q2sq + MVsq

    res = & 
     & -(((-MVsq - q1sq)* &
     &       (q2sq*(q1sq - q2sq + qHsq) + & 
     &         MVsq*(-q1sq + q2sq + qHsq))* &
     &       Log(delta2/delta1))/ &
     &     (-(MVsq*(q1sq - q2sq)**2) + delta1*delta2*qHsq) &
     &     ) + 2*Log(1 + q1sq/Mv**2) + &
     &  (delta2*(-MVsq - q1sq)* &
     &     (-2*MVsq + q1sq + q2sq - qHsq)* &
     &     Sqrt(qHsq/(4*MVsq + qHsq))* &
     &     Log((-Sqrt(qHsq) + Sqrt(4*MVsq + qHsq))/ &
     &       (Sqrt(qHsq) + Sqrt(4*MVsq + qHsq))))/ &
     &   (-(MVsq*(q1sq - q2sq)**2) + delta1*delta2*qHsq)
  end function f1_analytic
  ! Eq. 9 + 10 of 1906.10899 
  subroutine f1_integrand(x, res)
    real(dp), intent(in) :: x
    real(dp), intent(out) :: res(:)
    ! Internal
    real(dp) :: delta1,delta2,r1,r2,r12

    r1 = q1sq_rk * x + q2sq_rk * (one - x) - qHsq_rk * x *(one - x)
    r2 = qHsq_rk * x * (one - x) + MVsq
    r12 = r1 + r2
    delta1 = q1sq_rk + MVsq
    delta2 = q2sq_rk + MVsq

    res(1) = delta1*delta2/r12**2 * (log(r12**2/(r2*MVsq)) + (r1 - r2)/r2)
  end subroutine f1_integrand

!  function f2(q1sq,q2sq,qHsq,MV) result(res)
!    real(dp), intent(in) :: q1sq,q2sq,qHsq,MV
!    real(dp) :: res
!!
!    ! Save input for runge kutta routine
!    q1sq_rk = q1sq
!    q2sq_rk = q2sq
!    qHsq_rk = qHsq
!    MVsq    = MV**2
!!
!    res = scalar_integral(f2_integrand,zero,one)
!  end function f2
!  ! Eq. 9 + 10 of 1906.10899 
!  subroutine f2_integrand(x, res)
!    real(dp), intent(in) :: x
!    real(dp) :: res
!    ! Internal
!    real(dp) :: delta1,delta2,r1,r2,r12
!    complex(dp) :: HPL2
!    complex(dp) :: z
!    integer :: n1,n2
!    real(dp) :: Li2
!
!    ! This returns Li2
!    n1 = 0
!    n2 = 1
!!    print*, 'HPL2', HPL2(n1,n2,z)
!    
!    r1 = q1sq_rk * x + q2sq_rk * (one - x) - qHsq_rk * x *(one - x)
!    r2 = qHsq_rk * x * (one - x) + MVsq
!    r12 = r1 + r2
!    delta1 = q1sq_rk + MVsq
!    delta2 = q2sq_rk + MVsq
!
!    z = r1/r12
!    Li2 = real(HPL2(n1,n2,z))
!    res = delta1*delta2/r12**2 * ( &
!         & (log(r12**2/(r2*MVsq)) + (r1 - r2)/r2)**2 &
!         & - log(r12/r2)**2 - two*r12/r2*log(r12/r2) &
!         & - two*Li2 - ((r1-r2)/r2)**2 + two*zeta2)
!  end subroutine f2_integrand

    function tri_2loop(MV, p1x, p2x, p2y, lambda) result(res)
    use nonfact_expressions
    use incl_parameters
    real(dp), intent(in) :: MV, p1x, p2x, p2y, lambda
    real(dp) :: TT012, TT022, TT12, TT22, logl, TT(2)
    real(dp) :: res
    p1x_rk = p1x
    p2x_rk = p2x
    p2y_rk = p2y
    MVsq = MV**2
    logl = log(lambda/MVsq)
    ! Tri2lKirill = T012l + Integrate[  2*Re[ T022l]  ,{\[Xi],0,2 \[Pi] }] +  Integrate[ 2*Re[T12l]  ,{\[Xi],0,2 \[Pi] }]*Log[\[Lambda]/Mv^2]  +  T22l*Log[\[Lambda]/Mv^2]^2

    TT012 = t012(MVsq, pi, p1x, p2x, p2y)
    TT22 = t22(MVsq, pi, p1x, p2x, p2y)
    ! t01 and t022 in one adaptive integration (they share the roots
    ! and logs); t12 = -2 t01 pointwise
    TT = adaptive_integral(tri_integrand, 2, zero, 2.0_dp*pi)
    TT022 = TT(2)
    TT12 = -two*TT(1)
    
    res = (TT012 + two*TT022 + two*TT12*logl + TT22*logl**2)/TT22
  end function tri_2loop

  function chinf(qT1,qT2,MV) result(res)
    real(dp), intent(in) :: qT1(1:2), qT2(1:2),MV
    real(dp) :: qTH(1:2), ptH, pt1, pt2, res, lambda
    real(dp) :: q1rot(1:3), q2rot(1:3), cosphi, sinphi
    double precision, parameter :: z(1:3) = (/zero, zero, one/)

    q1rot = zero
    q2rot = zero

    q1rot(1:2) = qT1(1:2)
    q2rot(1:2) = qT2(1:2)
    if(abs(qT1(1)).gt.0d0) then
       cosphi = qT1(1)/sqrt(qT1(1)**2 + qT1(2)**2)
       sinphi = sqrt(one - cosphi**2)
       if(qT1(2).gt.zero) then
          sinphi = -sinphi
       endif
    else
       cosphi = one
       sinphi = zero
    endif
    call mrotate(z, sinphi, cosphi, q1rot(1:3))
    call mrotate(z, sinphi, cosphi, q2rot(1:3))

    qTH = - qT1 - qT2
    ptH = qTH(1)**2 + qTH(2)**2
    pt1 = qT1(1)**2 + qT1(2)**2
    pt2 = qT2(1)**2 + qT2(2)**2

    !    res = f1(pt1,pt2,ptH,MV)**2 - f2(pt1,pt2,ptH,MV)

    ! Gluon-mass regulator lambda = nf_regfact * MV^2; the result does
    ! not depend on it (the log(lambda/MV^2) terms of the 1-loop square
    ! and the 2-loop term cancel)
    lambda = nf_regfact * MV**2
    res = zero
    if(oneloop_on) res = res + (f1_analytic(pt1,pt2,ptH,MV) - log(lambda/MV**2))**2
!    if(twoloop_on) res = res - f2(pt1,pt2,ptH,MV)
    if(twoloop_on) res = res - tri_2loop(MV,q1rot(1),q2rot(1),q2rot(2),lambda)
  end function chinf

  subroutine tri_integrand(x, res)
    use nonfact_expressions
    use incl_parameters, only: pi
    real(dp), intent(in) :: x
    real(dp), intent(out) :: res(:)
    call tri_integrands(MVsq, pi, p1x_rk, p2x_rk, p2y_rk, x, res(1), res(2))
  end subroutine tri_integrand

  ! Integral of a one-component integrand
  function scalar_integral(f, x0, x1) result(res)
    real(dp), intent(in) :: x0, x1
    real(dp) :: res, r(1)
    interface
       subroutine f(x, res)
         use helper
         implicit none
         real(dp), intent(in) :: x
         real(dp), intent(out) :: res(:)
       end subroutine f
    end interface
    r = adaptive_integral(f, 1, x0, x1)
    res = r(1)
  end function scalar_integral

  !----------------------------------------------------------------------
  ! Integrals of the n components of f over [x0,x1] by adaptive
  ! Gauss-Kronrod (7-point Gauss, 15-point Kronrod) quadrature: the range
  ! is split into 4 intervals, and the interval with the largest error
  ! estimate (relative to the integral of |f| of each component) is
  ! bisected until, for every component, the summed error estimate is
  ! below nf_epsrel times the integral of |f| (or maxint intervals are
  ! reached). The angular integrands are smooth for typical kinematics,
  ! but have peaks of width ~ MV/pT at large transverse momenta, which a
  ! fixed rule does not resolve. Components that share the expensive
  ! parts (roots, logs) are integrated together on the same nodes.
  function adaptive_integral(f, n, x0, x1) result(res)
    use incl_parameters, only: nf_epsrel
    integer, intent(in) :: n
    real(dp), intent(in) :: x0, x1
    real(dp) :: res(n)
    interface
       subroutine f(x, res)
         use helper
         implicit none
         real(dp), intent(in) :: x
         real(dp), intent(out) :: res(:)
       end subroutine f
    end interface
    integer, parameter :: ninit = 4, maxint = 200
    real(dp) :: lo(maxint), hi(maxint), r(n,maxint), e(n,maxint), ra(n,maxint)
    real(dp) :: tol(n), worst(maxint)
    integer :: nint, i, k

    nint = ninit
    do i = 1, nint
       lo(i) = x0 + (x1 - x0)*(i-1)/nint
       hi(i) = x0 + (x1 - x0)*i/nint
       call gauss_kronrod_15(f, n, lo(i), hi(i), r(:,i), e(:,i), ra(:,i))
    enddo
    do
       tol = nf_epsrel*sum(ra(:,1:nint), dim=2)
       if (all(sum(e(:,1:nint), dim=2) <= tol) .or. nint >= maxint) exit
       do i = 1, nint
          worst(i) = maxval(e(:,i)/max(tol, tiny(one)))
       enddo
       k = maxloc(worst(1:nint), 1)
       nint = nint + 1
       lo(nint) = half*(lo(k) + hi(k))
       hi(nint) = hi(k)
       hi(k) = lo(nint)
       call gauss_kronrod_15(f, n, lo(k), hi(k), r(:,k), e(:,k), ra(:,k))
       call gauss_kronrod_15(f, n, lo(nint), hi(nint), r(:,nint), e(:,nint), ra(:,nint))
    enddo
    res = sum(r(:,1:nint), dim=2)
  end function adaptive_integral

  ! 15-point Kronrod estimates of the integrals of the n components of f
  ! over [a,b], the differences to the embedded 7-point Gauss rule as
  ! error estimates, and the integrals of |f| (QUADPACK's qk15 nodes and
  ! weights)
  subroutine gauss_kronrod_15(f, n, a, b, res, err, resabs)
    integer, intent(in) :: n
    real(dp), intent(in) :: a, b
    real(dp), intent(out) :: res(n), err(n), resabs(n)
    interface
       subroutine f(x, res)
         use helper
         implicit none
         real(dp), intent(in) :: x
         real(dp), intent(out) :: res(:)
       end subroutine f
    end interface
    real(dp), parameter :: xgk(8) = (/ 0.991455371120812639206854697526329_dp, &
         & 0.949107912342758524526189684047851_dp, 0.864864423359769072789712788640926_dp, &
         & 0.741531185599394439863864773280788_dp, 0.586087235467691130294144845693013_dp, &
         & 0.405845151377397166906606412076961_dp, 0.207784955007898467600689403773245_dp, &
         & 0.0_dp /)
    real(dp), parameter :: wgk(8) = (/ 0.022935322010529224963732008058970_dp, &
         & 0.063092092629978553290700663189204_dp, 0.104790010322250183839876322541518_dp, &
         & 0.140653259715525918745189590510238_dp, 0.169004726639267902826583426598550_dp, &
         & 0.190350578064785409913256402421014_dp, 0.204432940075298892414161999234649_dp, &
         & 0.209482141084727828012999174891714_dp /)
    real(dp), parameter :: wg(4) = (/ 0.129484966168869693270611432679082_dp, &
         & 0.279705391489276667901467771423780_dp, 0.381830050505118944950369775488975_dp, &
         & 0.417959183673469387755102040816327_dp /)
    real(dp) :: c, h, fc(n), f1(n), f2(n), rg(n), rk(n)
    integer :: j

    c = half*(a + b)
    h = half*(b - a)
    call f(c, fc)
    rg = fc*wg(4)
    rk = fc*wgk(8)
    resabs = abs(fc)*wgk(8)
    do j = 1, 7
       call f(c - h*xgk(j), f1)
       call f(c + h*xgk(j), f2)
       rk = rk + wgk(j)*(f1 + f2)
       resabs = resabs + wgk(j)*(abs(f1) + abs(f2))
       if (mod(j,2) == 0) rg = rg + wg(j/2)*(f1 + f2)
    enddo
    res = rk*h
    resabs = resabs*abs(h)
    err = abs((rk - rg)*h)
  end subroutine gauss_kronrod_15

end module nonfact
