module nonfact
  use helper
  use constants, pi_const => pi
!  use runge_kutta
  use incl_parameters
!  use grids
  implicit none

  public f1
!  public f2
  public chi_tri1
!  public chi_tri2

!  private f2_integrand
  private runge_kutta_dx

  real(dp),save :: q1sq_rk,q2sq_rk,qHsq_rk, s_sk, t_sk
  real(dp),save :: p1x_rk, p2x_rk, p2y_rk, p3x_rk, p3y_rk
  real(dp), save :: MVsq, MVHsq ! Mass of vector boson, ie either MW or MZ
contains
  ! Expressions from Kirill below 
  ! Eq. 9 + 10 of 1906.10899 
  function f1(q1sq,q2sq,qHsq,MV) result(res)
    real(dp), intent(in) :: q1sq,q2sq,qHsq,MV
    real(dp) :: res, delta1, delta2
    integer :: iterations

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
    integer :: iterations
    real(dp) :: res
    p1x_rk = p1x
    p2x_rk = p2x
    p2y_rk = p2y
    p3x_rk = p3x
    p3y_rk = p3y
    iterations = niter
    MVsq = MV**2
    MVHsq = MVH2 
    logl = log(lambda/MVsq)
    ! Box2l = B012l + Integrate[  2*Re[ B022l]  ,{\[Xi],0,2 \[Pi] }] +  Integrate[ 2*Re[B12l]  ,{\[Xi],0,2 \[Pi] }]*Log[\[Lambda]/Mv^2]  +  B22l*Log[\[Lambda]/Mv^2]^2

    BB012 = b012(MVsq, MVH2, pi, p1x, p2x, p2y, p3x, p3y)
    BB22 = b22(MVsq, MVH2, pi, p1x, p2x, p2y, p3x, p3y)

    call runge_kutta_dx(b022_integrand,zero,2.0_dp*pi,zero,BB022,iterations)
    call runge_kutta_dx(b12_integrand,zero,2.0_dp*pi,zero,BB12,iterations)
    
    res = (BB012 + two*BB022 + two*BB12*logl + BB22*logl**2)/BB22
  end function box_2loop
  
  function box_1loop_new(MV, MVH2, p1x, p2x, p2y, p3x, p3y, lambda) result(res)
    use nonfact_expressions
    use incl_parameters
    real(dp), intent(in) :: MV, MVH2, p1x, p2x, p2y, p3x, p3y, lambda
    real(dp) :: BB01, BB11,logl
    integer :: iterations
    real(dp) :: res
    p1x_rk = p1x
    p2x_rk = p2x
    p2y_rk = p2y
    p3x_rk = p3x
    p3y_rk = p3y
    iterations = niter
    MVsq = MV**2
    MVHsq = MVH2
    logl = log(lambda/MVsq)
    ! Box1l = Integrate[ 2*Re[B01l]  ,{\[Xi],0,2 \[Pi] }]+  B11l*Log[\[Lambda]/Mv^2]

    BB11 = b11(MVsq, MVH2, pi, p1x, p2x, p2y, p3x, p3y)

    call runge_kutta_dx(b01_integrand,zero,2.0_dp*pi,zero,BB01,iterations)
    ! Minus sign because of factored out 1/i
    res = -(two*BB01 + BB11*logl)/BB11
  end function box_1loop_new
  
  function tri_2loop(MV, p1x, p2x, p2y, lambda) result(res)
    use nonfact_expressions
    use incl_parameters
    real(dp), intent(in) :: MV, p1x, p2x, p2y, lambda
    real(dp) :: TT012, TT022, TT12, TT22,logl
    integer :: iterations
    real(dp) :: res
    p1x_rk = p1x
    p2x_rk = p2x
    p2y_rk = p2y
    iterations = niter
    MVsq = MV**2
    logl = log(lambda/MVsq)
    ! Tri2lKirill = T012l + Integrate[  2*Re[ T022l]  ,{\[Xi],0,2 \[Pi] }] +  Integrate[ 2*Re[T12l]  ,{\[Xi],0,2 \[Pi] }]*Log[\[Lambda]/Mv^2]  +  T22l*Log[\[Lambda]/Mv^2]^2

    TT012 = t012(MVsq, pi, p1x, p2x, p2y)
    TT22 = t22(MVsq, pi, p1x, p2x, p2y)

    call runge_kutta_dx(t022_integrand,zero,2.0_dp*pi,zero,TT022,iterations)
    call runge_kutta_dx(t12_integrand,zero,2.0_dp*pi,zero,TT12,iterations)
    
    res = (TT012 + two*TT022 + two*TT12*logl + TT22*logl**2)/TT22
  end function tri_2loop
  
  function tri_1loop(MV, p1x, p2x, p2y, lambda) result(res)
    use nonfact_expressions
    use incl_parameters
    real(dp), intent(in) :: MV, p1x, p2x, p2y, lambda

    real(dp) :: TT01, TT11,logl
    integer :: iterations
    real(dp) :: res
    p1x_rk = p1x
    p2x_rk = p2x
    p2y_rk = p2y
    iterations = niter
    MVsq = MV**2
    logl = log(lambda/MVsq)
    ! Tri1lKirill = Integrate[ 2*Re[T01l]  ,{\[Xi],0,2 \[Pi] }]+  T11l*Log[\[Lambda]/Mv^2]

    TT11 = t11(MVsq, pi, p1x, p2x, p2y)

    call runge_kutta_dx(t01_integrand,zero,2.0_dp*pi,zero,TT01,iterations)
    
    ! Minus sign because of factored out 1/i
    res = -(two*TT01 + TT11*logl)/TT11
  end function tri_1loop

  function b022_integrand(x,y) result(res)
    use nonfact_expressions
    use incl_parameters
    real(dp), intent(in) :: x,y
    real(dp) :: res
    res = b022(MVsq, MVHsq, pi, p1x_rk, p2x_rk, p2y_rk, p3x_rk, p3y_rk, x)
  end function b022_integrand

  function b12_integrand(x,y) result(res)
    use nonfact_expressions
    use incl_parameters
    real(dp), intent(in) :: x,y
    real(dp) :: res
    res = b12(MVsq, MVHsq, pi, p1x_rk, p2x_rk, p2y_rk, p3x_rk, p3y_rk, x)
  end function b12_integrand

  function b01_integrand(x,y) result(res)
    use nonfact_expressions
    use incl_parameters
    real(dp), intent(in) :: x,y
    real(dp) :: res
    res = b01(MVsq, MVHsq, pi, p1x_rk, p2x_rk, p2y_rk, p3x_rk, p3y_rk, x)
  end function b01_integrand

  function t022_integrand(x,y) result(res)
    use nonfact_expressions
    use incl_parameters
    real(dp), intent(in) :: x,y
    real(dp) :: res
    res = t022(MVsq, pi, p1x_rk, p2x_rk, p2y_rk, x)
  end function t022_integrand

  function t12_integrand(x,y) result(res)
    use nonfact_expressions
    use incl_parameters
    real(dp), intent(in) :: x,y
    real(dp) :: res
    res = t12(MVsq, pi, p1x_rk, p2x_rk, p2y_rk, x)
  end function t12_integrand

  function t01_integrand(x,y) result(res)
    use nonfact_expressions
    use incl_parameters
    real(dp), intent(in) :: x,y
    real(dp) :: res
    res = t01(MVsq, pi, p1x_rk, p2x_rk, p2y_rk, x)
  end function t01_integrand

  subroutine runge_kutta_dx(f,x0,x1,y0,integral,n_iter)

    real(dp), intent(in) :: x0,x1,y0
    integer, intent(in) :: n_iter
    real(dp), intent(out) :: integral
    real(dp) :: k(1:4), delx_6, delx_2
    real(dp) :: x(1:4), delx
    integer :: iter
    logical :: debug = .false.
    interface
       function f(x,y) result(res)
         use helper
         implicit none
         real(dp), intent(in) :: x, y
         real(dp) :: res
       end function f
    end interface

    if(n_iter.lt.1) then
       print*, n_iter
       stop 'n_iter should be positive integer'
    endif

    integral = y0 ! Initial value
    delx = (x1 - x0) / n_iter ! Step size
    delx_6 = delx/six
    delx_2 = delx/two
    k = zero
    x(1) = x0
    x(2) = x0 + half * delx
    !      x(3) = x0 + half * delx
    x(3) = x(2)
    x(4) = x0 + delx

    if(debug) then
       print*, 'x0        ', x0
       print*, 'x1        ', x1
       print*, 'y0        ', y0
       print*, 'n_iter    ', n_iter
       print*, 'x(1:4)    ', x
       print*, 'delx      ', delx
       print*, 'delx_2    ', delx_2
       print*, 'delx_6    ', delx_6
       print*, ''
    endif

    do iter = 1, n_iter
       k(1) = f(x(1),integral)
       k(2) = f(x(2),integral + delx_2 * k(1))
       k(3) = f(x(3),integral + delx_2 * k(2))
       k(4) = f(x(4),integral + delx * k(3))

       integral = integral + &
            & delx_6 * (k(1) + two*k(2) + two*k(3) + k(4))
       x(1:4) = x(1:4) + delx
       if(debug) then
          print*, 'x0        ', x0
          print*, 'x1        ', x1
          print*, 'y0        ', y0
          print*, 'iter      ', iter
          print*, 'x(1:4)    ', x
          print*, 'delx      ', delx
          print*, 'delx_2    ', delx_2
          print*, 'delx_6    ', delx_6
          print*, 'k(1:4)    ', k
          print*, 'integral  ', integral
          if(integral.ne.integral) stop
       endif
    enddo

  end subroutine runge_kutta_dx

end module nonfact
