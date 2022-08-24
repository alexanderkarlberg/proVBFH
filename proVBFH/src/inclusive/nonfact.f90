module nonfact
  use helper
  use constants
  use runge_kutta
  use incl_parameters
!  use grids
  implicit none

  public f1
!  public f2
  public chinf

  private f1_integrand
!  private f2_integrand
  private runge_kutta_dx

  real(dp),save :: q1sq_rk,q2sq_rk,qHsq_rk
  real(dp),save :: p1x_rk, p2x_rk, p2y_rk, p3x_rk, p3y_rk
  real(dp), save :: MVsq ! Mass of vector boson, ie either MW or MZ
contains

  function f1(q1sq,q2sq,qHsq,MV) result(res)
    real(dp), intent(in) :: q1sq,q2sq,qHsq,MV
    real(dp) :: res
    integer :: iterations

    ! Save input for runge kutta routine
    q1sq_rk = q1sq
    q2sq_rk = q2sq
    qHsq_rk = qHsq
    MVsq    = MV**2
    iterations = niter

    call runge_kutta_dx(f1_integrand,zero,one,zero,res,iterations)
  end function f1

  function f1_analytic(q1sq,q2sq,qHsq,MV) result(res)
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
  end function f1_analytic
  ! Eq. 9 + 10 of 1906.10899 
  function f1_integrand(x,y) result(res)
    real(dp), intent(in) :: x,y
    real(dp) :: res
    ! Internal
    real(dp) :: delta1,delta2,r1,r2,r12

    r1 = q1sq_rk * x + q2sq_rk * (one - x) - qHsq_rk * x *(one - x)
    r2 = qHsq_rk * x * (one - x) + MVsq
    r12 = r1 + r2
    delta1 = q1sq_rk + MVsq
    delta2 = q2sq_rk + MVsq

    res = delta1*delta2/r12**2 * (log(r12**2/(r2*MVsq)) + (r1 - r2)/r2)
  end function f1_integrand

!  function f2(q1sq,q2sq,qHsq,MV) result(res)
!    real(dp), intent(in) :: q1sq,q2sq,qHsq,MV
!    real(dp) :: res
!    integer :: iterations
!
!    ! Save input for runge kutta routine
!    q1sq_rk = q1sq
!    q2sq_rk = q2sq
!    qHsq_rk = qHsq
!    MVsq    = MV**2
!    iterations = niter
!
!    call runge_kutta_dx(f2_integrand,zero,one,zero,res,iterations)
!  end function f2
!  ! Eq. 9 + 10 of 1906.10899 
!  function f2_integrand(x,y) result(res)
!    real(dp), intent(in) :: x,y
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
!  end function f2_integrand

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

  function chinf(qT1,qT2,MV) result(res)
    real(dp), intent(in) :: qT1(1:2), qT2(1:2),MV
    real(dp) :: qTH(1:2), ptH, pt1, pt2, res
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

    res = zero
    if(oneloop_on) res = res + f1_analytic(pt1,pt2,ptH,MV)**2 
!    if(twoloop_on) res = res - f2(pt1,pt2,ptH,MV)
    if(twoloop_on) res = res - tri_2loop(MV,q1rot(1),q2rot(1),q2rot(2),MV**2)
  end function chinf

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
end module nonfact
