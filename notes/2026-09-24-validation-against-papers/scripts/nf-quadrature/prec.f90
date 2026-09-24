program prec
  use integrands
  implicit none
  character(len=3) :: tag
  integer :: ios, i, k, n
  real(dp) :: rs, rd, scale, d(10000)
  character(len=*), parameter :: names(3) = (/ 'b01 ', 'b022', 't022' /)
  do k = 1, 3
     n = 0
     open(10, file='args.dat')
     do
        read(10,*,iostat=ios) tag, a
        if (ios /= 0) exit
        if ((k==1 .and. tag/='B1') .or. (k==2 .and. tag/='B2') .or. (k==3 .and. tag/='T2')) cycle
        kind_int = k; n = n + 1
        single = .false.; rd = trap(1024)
        scale = 0
        do i = 0, 1023; scale = scale + abs(f(2*pi*i/1024)); end do
        scale = scale*2*pi/1024
        single = .true.; rs = trap(1024)
        d(n) = abs(rs-rd)/abs(rd)
     end do
     close(10)
     call sort(d(1:n))
     write(*,'(a5,a,3es10.2)') names(k), ': |single/double - 1| of the integral, 50%/99%/max:', d(n/2), d(99*n/100), d(n)
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
end program prec
