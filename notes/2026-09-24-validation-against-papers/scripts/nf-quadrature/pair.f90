program pair
  use nonfact_expressions
  implicit none
  real*8 :: a(7), xi, b1, b2, c1, c2, d1, d2, pi, sc
  complex*16 :: r(6), dd
  integer :: k, j
  character(len=3) :: tag
  integer :: ios, i, n
  pi = 4d0*atan(1d0); d1 = 0; d2 = 0; n = 0
  open(10, file='args.dat')
  do
     read(10,*,iostat=ios) tag, a
     if (ios /= 0) exit
     if (tag /= 'B1' .and. tag /= 'B2') cycle
     do i = 0, 63
        xi = 2*pi*i/64 + 0.01d0
        b1 = b01(a(1)**2,a(2),pi,a(3),a(4),a(5),a(6),a(7),xi)
        b2 = dble(b022(a(1)**2,a(2),pi,a(3),a(4),a(5),a(6),a(7),xi))
        call box_integrands(a(1)**2,a(2),pi,a(3),a(4),a(5),a(6),a(7),xi,c1,c2)
        r(1)=r1(a(1)**2,a(3),xi); r(2)=r2(a(1)**2,a(3),xi); r(3)=r3(a(1)**2,a(4),a(5),xi); r(4)=r4(a(1)**2,a(4),a(5),xi)
        r(5)=r5(a(2),a(3),a(6),a(7),xi); r(6)=r6(a(2),a(3),a(6),a(7),xi); sc = 0
        do k = 1, 5, 2
           dd = pi*r(k)
           do j = 1, 6
              if (j /= k) dd = dd*(r(k)-r(j))
           end do
           sc = sc + abs(log(-r(k)/a(1))/dd)
        end do
        d1 = max(d1, abs(c1-b1)/sc); d2 = max(d2, abs(c2-b2)/max(abs(b2),1d-300)); n = n+1
     end do
  end do
  print '(a,i7,a,2es10.2)', ' points', n, '  max rel diff b01, b022:', d1, d2
end program pair
