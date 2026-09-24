program worst
  use integrands
  implicit none
  character(len=3) :: tag
  integer :: ios, i, n, iw
  real(dp) :: ref, scale, e, ew, aw(7), errs(10000), s
  real(dp) :: v
  kind_int = 3; single = .false.; n = 0; ew = 0
  open(10, file='args.dat')
  do
     read(10,*,iostat=ios) tag, a
     if (ios /= 0) exit
     if (tag /= 'T2') cycle
     n = n + 1
     ref = trap(4096)
     scale = 0
     do i = 0, 4095; scale = scale + abs(f(2*pi*i/4096)); end do
     scale = scale*2*pi/4096
     e = abs(rk(100)-ref)/scale; errs(n) = e
     if (e > ew) then; ew = e; aw = a; end if
  end do
  call sort(errs(1:n))
  write(*,'(a,5es10.2)') ' RK100 t022 error quantiles 50/90/99/99.9/max: ', errs(n/2), errs(9*n/10), errs(99*n/100), errs(999*n/1000), errs(n)
  write(*,'(a,7es12.4)') ' worst point args: ', aw
  a = aw
  do i = 0, 64
     v = f(2*pi*i/64)
     write(*,'(f8.4,es14.5)') 2*pi*i/64, v
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
end program worst
