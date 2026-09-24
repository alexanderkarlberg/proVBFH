module nonfact
  use helper
  use constants, pi_const => pi
!  use incl_parameters
!  use grids
  implicit none

  public f1
!  public f2
  public chi_tri1
!  public chi_tri2

!  private f2_integrand
  private adaptive_integral, gauss_kronrod_15, scalar_integral
  private box_angular_integrals, box_pair_integrand

  real(dp),save :: q1sq_rk,q2sq_rk,qHsq_rk, s_sk, t_sk
  real(dp),save :: p1x_rk, p2x_rk, p2y_rk, p3x_rk, p3y_rk
  real(dp), save :: MVsq, MVHsq ! Mass of vector boson, ie either MW or MZ
contains
  ! Expressions from Kirill below 
  ! Eq. 9 + 10 of 1906.10899 
  function f1(q1sq,q2sq,qHsq,MV) result(res)
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
  end function f1

  
  ! Eq. 8 below of 1906.10899 
  function chi_tri1(q1sq,q2sq,qHsq,MV,lambda) result(res)
    real(dp), intent(in) :: q1sq,q2sq,qHsq,MV,lambda
    real(dp) :: res

    res = - log(lambda/MV**2) + f1(q1sq,q2sq,qHsq,MV)
  end function chi_tri1

  ! Lorenzo's expressions below
  
  function box_2loop(MV, MVH2, p1x, p2x, p2y, p3x, p3y, lambda) result(res)
    use nonfact_expressions
    use incl_parameters
    real(dp), intent(in) :: MV, MVH2, p1x, p2x, p2y, p3x, p3y, lambda
    real(dp) :: BB012, BB022, BB12, BB22, BB01, logl
    real(dp) :: res
    p1x_rk = p1x
    p2x_rk = p2x
    p2y_rk = p2y
    p3x_rk = p3x
    p3y_rk = p3y
    MVsq = MV**2
    MVHsq = MVH2 
    logl = log(lambda/MVsq)
    ! Box2l = B012l + Integrate[  2*Re[ B022l]  ,{\[Xi],0,2 \[Pi] }] +  Integrate[ 2*Re[B12l]  ,{\[Xi],0,2 \[Pi] }]*Log[\[Lambda]/Mv^2]  +  B22l*Log[\[Lambda]/Mv^2]^2

    BB012 = b012(MVsq, MVH2, pi, p1x, p2x, p2y, p3x, p3y)
    BB22 = b22(MVsq, MVH2, pi, p1x, p2x, p2y, p3x, p3y)

    call box_angular_integrals(MV, MVH2, p1x, p2x, p2y, p3x, p3y, BB01, BB022)
    ! only needed for lambda /= MV^2
    BB12 = zero
    if (logl /= zero) BB12 = scalar_integral(b12_integrand,zero,2.0_dp*pi)
    
    res = (BB012 + two*BB022 + two*BB12*logl + BB22*logl**2)/BB22
  end function box_2loop
  
  function box_1loop_new(MV, MVH2, p1x, p2x, p2y, p3x, p3y, lambda) result(res)
    use nonfact_expressions
    use incl_parameters
    real(dp), intent(in) :: MV, MVH2, p1x, p2x, p2y, p3x, p3y, lambda
    real(dp) :: BB01, BB11, BB022, logl
    real(dp) :: res
    p1x_rk = p1x
    p2x_rk = p2x
    p2y_rk = p2y
    p3x_rk = p3x
    p3y_rk = p3y
    MVsq = MV**2
    MVHsq = MVH2
    logl = log(lambda/MVsq)
    ! Box1l = Integrate[ 2*Re[B01l]  ,{\[Xi],0,2 \[Pi] }]+  B11l*Log[\[Lambda]/Mv^2]

    BB11 = b11(MVsq, MVH2, pi, p1x, p2x, p2y, p3x, p3y)

    call box_angular_integrals(MV, MVH2, p1x, p2x, p2y, p3x, p3y, BB01, BB022)
    ! Minus sign because of factored out 1/i
    res = -(two*BB01 + BB11*logl)/BB11
  end function box_1loop_new
  
  function tri_2loop(MV, p1x, p2x, p2y, lambda) result(res)
    use nonfact_expressions
    use incl_parameters
    real(dp), intent(in) :: MV, p1x, p2x, p2y, lambda
    real(dp) :: TT012, TT022, TT12, TT22,logl
    real(dp) :: res
    p1x_rk = p1x
    p2x_rk = p2x
    p2y_rk = p2y
    MVsq = MV**2
    logl = log(lambda/MVsq)
    ! Tri2lKirill = T012l + Integrate[  2*Re[ T022l]  ,{\[Xi],0,2 \[Pi] }] +  Integrate[ 2*Re[T12l]  ,{\[Xi],0,2 \[Pi] }]*Log[\[Lambda]/Mv^2]  +  T22l*Log[\[Lambda]/Mv^2]^2

    TT012 = t012(MVsq, pi, p1x, p2x, p2y)
    TT22 = t22(MVsq, pi, p1x, p2x, p2y)

    TT022 = scalar_integral(t022_integrand,zero,2.0_dp*pi)
    ! only needed for lambda /= MV^2
    TT12 = zero
    if (logl /= zero) TT12 = scalar_integral(t12_integrand,zero,2.0_dp*pi)
    
    res = (TT012 + two*TT022 + two*TT12*logl + TT22*logl**2)/TT22
  end function tri_2loop
  
  function tri_1loop(MV, p1x, p2x, p2y, lambda) result(res)
    use nonfact_expressions
    use incl_parameters
    real(dp), intent(in) :: MV, p1x, p2x, p2y, lambda

    real(dp) :: TT01, TT11,logl
    real(dp) :: res
    p1x_rk = p1x
    p2x_rk = p2x
    p2y_rk = p2y
    MVsq = MV**2
    logl = log(lambda/MVsq)
    ! Tri1lKirill = Integrate[ 2*Re[T01l]  ,{\[Xi],0,2 \[Pi] }]+  T11l*Log[\[Lambda]/Mv^2]

    TT11 = t11(MVsq, pi, p1x, p2x, p2y)

    TT01 = scalar_integral(t01_integrand,zero,2.0_dp*pi)
    
    ! Minus sign because of factored out 1/i
    res = -(two*TT01 + TT11*logl)/TT11
  end function tri_1loop

  !----------------------------------------------------------------------
  ! The azimuthal integrals of b01 (1-loop box) and b022 (2-loop box),
  ! computed together (they share the roots and logs, and the peaks of
  ! the integrands). box_1loop_new and box_2loop are called with the
  ! same arguments for each boson and t/u channel, so the last few
  ! results are cached (keyed on the exact arguments, so a new
  ! phase-space point or nf_epsrel always recomputes).
  subroutine box_angular_integrals(MV, MVH2, p1x, p2x, p2y, p3x, p3y, BB01, BB022)
    use incl_parameters, only: nf_epsrel, pi
    real(dp), intent(in) :: MV, MVH2, p1x, p2x, p2y, p3x, p3y
    real(dp), intent(out) :: BB01, BB022
    integer, parameter :: ncache = 4
    real(dp), save :: key(8,ncache) = -one, val(2,ncache) = zero
    integer, save :: next = 1
    real(dp) :: args(8), res(2)
    integer :: i

    args = (/ MV, MVH2, p1x, p2x, p2y, p3x, p3y, nf_epsrel /)
    do i = 1, ncache
       if (all(key(:,i) == args)) then
          BB01 = val(1,i)
          BB022 = val(2,i)
          return
       endif
    enddo
    p1x_rk = p1x
    p2x_rk = p2x
    p2y_rk = p2y
    p3x_rk = p3x
    p3y_rk = p3y
    MVsq = MV**2
    MVHsq = MVH2
    res = adaptive_integral(box_pair_integrand, 2, zero, 2.0_dp*pi)
    BB01 = res(1)
    BB022 = res(2)
    key(:,next) = args
    val(:,next) = res
    next = mod(next, ncache) + 1
  end subroutine box_angular_integrals

  subroutine box_pair_integrand(x, res)
    use nonfact_expressions
    use incl_parameters
    real(dp), intent(in) :: x
    real(dp), intent(out) :: res(:)
    call box_integrands(MVsq, MVHsq, pi, p1x_rk, p2x_rk, p2y_rk, p3x_rk, p3y_rk, x, &
         & res(1), res(2))
  end subroutine box_pair_integrand

  subroutine b12_integrand(x, res)
    use nonfact_expressions
    use incl_parameters
    real(dp), intent(in) :: x
    real(dp), intent(out) :: res(:)
    res(1) = b12(MVsq, MVHsq, pi, p1x_rk, p2x_rk, p2y_rk, p3x_rk, p3y_rk, x)
  end subroutine b12_integrand

  subroutine t022_integrand(x, res)
    use nonfact_expressions
    use incl_parameters
    real(dp), intent(in) :: x
    real(dp), intent(out) :: res(:)
    res(1) = t022(MVsq, pi, p1x_rk, p2x_rk, p2y_rk, x)
  end subroutine t022_integrand

  subroutine t12_integrand(x, res)
    use nonfact_expressions
    use incl_parameters
    real(dp), intent(in) :: x
    real(dp), intent(out) :: res(:)
    res(1) = t12(MVsq, pi, p1x_rk, p2x_rk, p2y_rk, x)
  end subroutine t12_integrand

  subroutine t01_integrand(x, res)
    use nonfact_expressions
    use incl_parameters
    real(dp), intent(in) :: x
    real(dp), intent(out) :: res(:)
    res(1) = t01(MVsq, pi, p1x_rk, p2x_rk, p2y_rk, x)
  end subroutine t01_integrand

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
