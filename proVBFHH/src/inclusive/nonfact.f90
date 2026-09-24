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
  private adaptive_integral, gauss_kronrod_15

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
    real(dp) :: BB012, BB022, BB12, BB22,logl
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

    BB022 = adaptive_integral(b022_integrand,zero,2.0_dp*pi)
    ! only needed for lambda /= MV^2
    BB12 = zero
    if (logl /= zero) BB12 = adaptive_integral(b12_integrand,zero,2.0_dp*pi)
    
    res = (BB012 + two*BB022 + two*BB12*logl + BB22*logl**2)/BB22
  end function box_2loop
  
  function box_1loop_new(MV, MVH2, p1x, p2x, p2y, p3x, p3y, lambda) result(res)
    use nonfact_expressions
    use incl_parameters
    real(dp), intent(in) :: MV, MVH2, p1x, p2x, p2y, p3x, p3y, lambda
    real(dp) :: BB01, BB11,logl
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

    BB01 = adaptive_integral(b01_integrand,zero,2.0_dp*pi)
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

    TT022 = adaptive_integral(t022_integrand,zero,2.0_dp*pi)
    ! only needed for lambda /= MV^2
    TT12 = zero
    if (logl /= zero) TT12 = adaptive_integral(t12_integrand,zero,2.0_dp*pi)
    
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

    TT01 = adaptive_integral(t01_integrand,zero,2.0_dp*pi)
    
    ! Minus sign because of factored out 1/i
    res = -(two*TT01 + TT11*logl)/TT11
  end function tri_1loop

  function b022_integrand(x) result(res)
    use nonfact_expressions
    use incl_parameters
    real(dp), intent(in) :: x
    real(dp) :: res
    res = b022(MVsq, MVHsq, pi, p1x_rk, p2x_rk, p2y_rk, p3x_rk, p3y_rk, x)
  end function b022_integrand

  function b12_integrand(x) result(res)
    use nonfact_expressions
    use incl_parameters
    real(dp), intent(in) :: x
    real(dp) :: res
    res = b12(MVsq, MVHsq, pi, p1x_rk, p2x_rk, p2y_rk, p3x_rk, p3y_rk, x)
  end function b12_integrand

  function b01_integrand(x) result(res)
    use nonfact_expressions
    use incl_parameters
    real(dp), intent(in) :: x
    real(dp) :: res
    res = b01(MVsq, MVHsq, pi, p1x_rk, p2x_rk, p2y_rk, p3x_rk, p3y_rk, x)
  end function b01_integrand

  function t022_integrand(x) result(res)
    use nonfact_expressions
    use incl_parameters
    real(dp), intent(in) :: x
    real(dp) :: res
    res = t022(MVsq, pi, p1x_rk, p2x_rk, p2y_rk, x)
  end function t022_integrand

  function t12_integrand(x) result(res)
    use nonfact_expressions
    use incl_parameters
    real(dp), intent(in) :: x
    real(dp) :: res
    res = t12(MVsq, pi, p1x_rk, p2x_rk, p2y_rk, x)
  end function t12_integrand

  function t01_integrand(x) result(res)
    use nonfact_expressions
    use incl_parameters
    real(dp), intent(in) :: x
    real(dp) :: res
    res = t01(MVsq, pi, p1x_rk, p2x_rk, p2y_rk, x)
  end function t01_integrand

  !----------------------------------------------------------------------
  ! Integral of f over [x0,x1] by adaptive Gauss-Kronrod (7-point Gauss,
  ! 15-point Kronrod) quadrature: the range is split into 4 intervals,
  ! and the interval with the largest error estimate is bisected until
  ! the summed error estimate is below nf_epsrel times the integral of
  ! |f| (or maxint intervals are reached). The angular integrands are
  ! smooth for typical kinematics, but have peaks of width ~ MV/pT at
  ! large transverse momenta, which a fixed rule does not resolve.
  function adaptive_integral(f, x0, x1) result(res)
    use incl_parameters, only: nf_epsrel
    real(dp), intent(in) :: x0, x1
    real(dp) :: res
    interface
       function f(x) result(res)
         use helper
         implicit none
         real(dp), intent(in) :: x
         real(dp) :: res
       end function f
    end interface
    integer, parameter :: ninit = 4, maxint = 200
    real(dp) :: lo(maxint), hi(maxint), r(maxint), e(maxint), ra(maxint)
    integer :: n, i, k

    n = ninit
    do i = 1, n
       lo(i) = x0 + (x1 - x0)*(i-1)/n
       hi(i) = x0 + (x1 - x0)*i/n
       call gauss_kronrod_15(f, lo(i), hi(i), r(i), e(i), ra(i))
    enddo
    do while (sum(e(1:n)) > nf_epsrel*sum(ra(1:n)) .and. n < maxint)
       k = maxloc(e(1:n), 1)
       n = n + 1
       lo(n) = half*(lo(k) + hi(k))
       hi(n) = hi(k)
       hi(k) = lo(n)
       call gauss_kronrod_15(f, lo(k), hi(k), r(k), e(k), ra(k))
       call gauss_kronrod_15(f, lo(n), hi(n), r(n), e(n), ra(n))
    enddo
    res = sum(r(1:n))
  end function adaptive_integral

  ! 15-point Kronrod estimate of the integral of f over [a,b], the
  ! difference to the embedded 7-point Gauss rule as error estimate,
  ! and the integral of |f| (QUADPACK's qk15 nodes and weights)
  subroutine gauss_kronrod_15(f, a, b, res, err, resabs)
    real(dp), intent(in) :: a, b
    real(dp), intent(out) :: res, err, resabs
    interface
       function f(x) result(res)
         use helper
         implicit none
         real(dp), intent(in) :: x
         real(dp) :: res
       end function f
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
    real(dp) :: c, h, fc, f1, f2, rg, rk
    integer :: j

    c = half*(a + b)
    h = half*(b - a)
    fc = f(c)
    rg = fc*wg(4)
    rk = fc*wgk(8)
    resabs = abs(fc)*wgk(8)
    do j = 1, 7
       f1 = f(c - h*xgk(j))
       f2 = f(c + h*xgk(j))
       rk = rk + wgk(j)*(f1 + f2)
       resabs = resabs + wgk(j)*(abs(f1) + abs(f2))
       if (mod(j,2) == 0) rg = rg + wg(j/2)*(f1 + f2)
    enddo
    res = rk*h
    resabs = resabs*abs(h)
    err = abs((rk - rg)*h)
  end subroutine gauss_kronrod_15

end module nonfact
