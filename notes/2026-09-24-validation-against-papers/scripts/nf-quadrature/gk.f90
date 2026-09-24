module gk
  use integrands
  implicit none
  integer :: ncalls
  ! Gauss-Kronrod 7-15 nodes/weights on [-1,1] (QUADPACK qk15)
  real(dp), parameter :: xgk(8) = (/ 0.991455371120812639206854697526329_dp, &
       0.949107912342758524526189684047851_dp, 0.864864423359769072789712788640926_dp, &
       0.741531185599394439863864773280788_dp, 0.586087235467691130294144845693013_dp, &
       0.405845151377397166906606412076961_dp, 0.207784955007898467600689403773245_dp, 0.0_dp /)
  real(dp), parameter :: wgk(8) = (/ 0.022935322010529224963732008058970_dp, &
       0.063092092629978553290700663189204_dp, 0.104790010322250183839876322541518_dp, &
       0.140653259715525918745189590510238_dp, 0.169004726639267902826583426598550_dp, &
       0.190350578064785409913256402421014_dp, 0.204432940075298892414161999234649_dp, &
       0.209482141084727828012999174891714_dp /)
  real(dp), parameter :: wg(4) = (/ 0.129484966168869693270611432679082_dp, &
       0.279705391489276667901467771423780_dp, 0.381830050505118944950369775488975_dp, &
       0.417959183673469387755102040816327_dp /)
contains
  subroutine qk15(a, b, res, err, resabs)
    real(dp), intent(in) :: a, b
    real(dp), intent(out) :: res, err, resabs
    real(dp) :: c, h, fc, f1, f2, rg, rk
    integer :: j
    c = 0.5_dp*(a+b); h = 0.5_dp*(b-a)
    fc = f(c); ncalls = ncalls + 1
    rg = fc*wg(4); rk = fc*wgk(8); resabs = abs(fc)*wgk(8)
    do j = 1, 7
       f1 = f(c - h*xgk(j)); f2 = f(c + h*xgk(j)); ncalls = ncalls + 2
       rk = rk + wgk(j)*(f1+f2); resabs = resabs + wgk(j)*(abs(f1)+abs(f2))
       if (mod(j,2) == 0) rg = rg + wg(j/2)*(f1+f2)
    end do
    res = rk*h; resabs = resabs*abs(h)
    err = abs((rk-rg)*h)
  end subroutine qk15
  ! adaptive bisection of the interval with the largest error, until the
  ! summed error is below epsrel * integral of |f|
  real(dp) function adapt(epsrel, ninit)
    real(dp), intent(in) :: epsrel
    integer, intent(in) :: ninit
    integer, parameter :: nmax = 200
    real(dp) :: lo(nmax), hi(nmax), r(nmax), e(nmax), ra(nmax), tot, etot, atot
    integer :: n, i, k
    n = ninit
    do i = 1, n
       lo(i) = 2*pi*(i-1)/n; hi(i) = 2*pi*i/n
       call qk15(lo(i), hi(i), r(i), e(i), ra(i))
    end do
    do
       tot = sum(r(1:n)); etot = sum(e(1:n)); atot = sum(ra(1:n))
       if (etot <= epsrel*atot .or. n >= nmax-1) exit
       k = maxloc(e(1:n), 1)
       n = n + 1
       lo(n) = 0.5_dp*(lo(k)+hi(k)); hi(n) = hi(k); hi(k) = lo(n)
       call qk15(lo(k), hi(k), r(k), e(k), ra(k))
       call qk15(lo(n), hi(n), r(n), e(n), ra(n))
    end do
    adapt = tot
  end function adapt
end module gk

program testgk
  use gk
  implicit none
  character(len=3) :: tag
  integer :: ios, i, k, n, ie, ni
  real(dp) :: ref, scale, epsv(3) = (/ 1e-6_dp, 1e-8_dp, 1e-10_dp /)
  integer :: ninits(2) = (/ 1, 4 /)
  real(dp), allocatable :: errs(:), calls(:)
  character(len=*), parameter :: names(3) = (/ 'b01 ', 'b022', 't022' /)
  allocate(errs(10000), calls(10000))
  single = .false.
  do k = 1, 3
    do ni = 1, 2
     do ie = 1, 3
        n = 0
        open(10, file='args.dat')
        do
           read(10,*,iostat=ios) tag, a
           if (ios /= 0) exit
           if ((k==1 .and. tag/='B1') .or. (k==2 .and. tag/='B2') .or. (k==3 .and. tag/='T2')) cycle
           kind_int = k; n = n + 1
           ref = trap(4096)
           scale = 0
           do i = 0, 4095; scale = scale + abs(f(2*pi*i/4096)); end do
           scale = scale*2*pi/4096
           ncalls = 0
           errs(n) = abs(adapt(epsv(ie), ninits(ni)) - ref)/scale
           calls(n) = ncalls
        end do
        close(10)
        call sort(errs(1:n))
        write(*,'(a5,a,i2,a,es7.0,a,f7.1,a,i5,a,3es9.2)') names(k), ' ninit', ninits(ni), ' eps', epsv(ie), &
             ': mean calls', sum(calls(1:n))/n, ' max', int(maxval(calls(1:n))), &
             '  err 50%/99%/max', errs(n/2), errs(99*n/100), errs(n)
     end do
    end do
  end do
contains
  subroutine sort(x)
    real(dp), intent(inout) :: x(:)
    integer :: i, j; real(dp) :: t
    do i = 2, size(x)
       t = x(i); j = i - 1
       do while (j >= 1)
          if (x(j) <= t) exit
          x(j+1) = x(j); j = j - 1
       end do
       x(j+1) = t
    end do
  end subroutine sort
end program testgk
