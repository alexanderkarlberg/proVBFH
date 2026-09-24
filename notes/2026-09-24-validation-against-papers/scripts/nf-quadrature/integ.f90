! Convergence study of the angular integrals of the NF corrections
module integrands
  use nonfact_expressions, only: sb01 => b01, sb022 => b022, st022 => t022
  use nonfact_expressions_dp, only: db01 => b01, db022 => b022, dt022 => t022
  implicit none
  integer, parameter :: dp = kind(1d0)
  real(dp), parameter :: pi = 3.14159265358979323846264338327950288_dp
  real(dp) :: a(7)       ! MV, MVH2, p1x, p2x, p2y, p3x, p3y
  integer  :: kind_int   ! 1: b01, 2: b022, 3: t022
  logical  :: single
contains
  real(dp) function f(xi)
    real(dp), intent(in) :: xi
    real(dp) :: MV2
    MV2 = a(1)**2
    select case (kind_int)
    case (1)
       if (single) then; f = sb01(MV2,a(2),pi,a(3),a(4),a(5),a(6),a(7),xi)
       else;             f = db01(MV2,a(2),pi,a(3),a(4),a(5),a(6),a(7),xi); end if
    case (2)
       if (single) then; f = sb022(MV2,a(2),pi,a(3),a(4),a(5),a(6),a(7),xi)
       else;             f = db022(MV2,a(2),pi,a(3),a(4),a(5),a(6),a(7),xi); end if
    case default
       if (single) then; f = st022(MV2,pi,a(3),a(4),a(5),xi)
       else;             f = dt022(MV2,pi,a(3),a(4),a(5),xi); end if
    end select
  end function f
  ! the current scheme: RK4 on dy/dx = f(x), i.e. Simpson with n panels
  real(dp) function rk(n)
    integer, intent(in) :: n
    real(dp) :: h, x, y
    integer :: i
    h = 2*pi/n; y = 0; x = 0
    do i = 1, n
       y = y + h/6*(f(x) + 2*f(x+h/2) + 2*f(x+h/2) + f(x+h))
       x = x + h
    end do
    rk = y
  end function rk
  real(dp) function trap(n)
    integer, intent(in) :: n
    integer :: i
    trap = 0
    do i = 0, n-1
       trap = trap + f(2*pi*i/n)
    end do
    trap = trap * 2*pi/n
  end function trap
  real(dp) function gauss(n)
    integer, intent(in) :: n
    real(dp) :: x(n), w(n)
    integer :: i
    call gauleg(n, x, w)
    gauss = 0
    do i = 1, n
       gauss = gauss + w(i)*f(pi*(x(i)+1))
    end do
    gauss = gauss*pi
  end function gauss
  ! Gauss-Legendre nodes and weights on [-1,1] (Newton on P_n)
  subroutine gauleg(n, x, w)
    integer, intent(in) :: n
    real(dp), intent(out) :: x(n), w(n)
    integer :: i, j, it
    real(dp) :: z, p1, p2, p3, pp
    do i = 1, (n+1)/2
       z = cos(pi*(i-0.25_dp)/(n+0.5_dp))
       do it = 1, 100
          p1 = 1; p2 = 0
          do j = 1, n
             p3 = p2; p2 = p1
             p1 = ((2*j-1)*z*p2 - (j-1)*p3)/j
          end do
          pp = n*(z*p1-p2)/(z*z-1)
          z = z - p1/pp
          if (abs(p1/pp) < 1e-16_dp) exit
       end do
       x(i) = -z; x(n+1-i) = z
       w(i) = 2/((1-z*z)*pp*pp); w(n+1-i) = w(i)
    end do
  end subroutine gauleg
end module integrands

