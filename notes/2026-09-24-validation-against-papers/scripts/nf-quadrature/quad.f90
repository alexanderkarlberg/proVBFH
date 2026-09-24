program quad
  use integrands
  implicit none
  character(len=3) :: tag
  integer :: ios, npt, i, k
  integer, parameter :: nt = 6, ng = 5
  integer :: ntr(nt) = (/ 8, 16, 24, 32, 48, 64 /), ngl(ng) = (/ 16, 32, 48, 64, 128 /)
  real(dp) :: ref, ref2, e_rk_sp, e_rk_dp, e_tr(nt), e_gl(ng), refchk
  real(dp) :: m_rk_sp, m_rk_dp, m_tr(nt), m_gl(ng), m_chk, scale
  character(len=*), parameter :: names(3) = (/ 'b01 ', 'b022', 't022' /)
  do k = 1, 3
     m_rk_sp = 0; m_rk_dp = 0; m_tr = 0; m_gl = 0; m_chk = 0; npt = 0
     open(10, file='args.dat')
     do
        read(10, *, iostat=ios) tag, a
        if (ios /= 0) exit
        if (k == 3 .and. tag /= 'T2') cycle
        if (k == 1 .and. tag /= 'B1') cycle
        if (k == 2 .and. tag /= 'B2') cycle
        kind_int = k; npt = npt + 1
        single = .false.
        ref = trap(4096); ref2 = trap(2048)
        ! scale: integral of |f|, to measure errors without cancellations
        scale = 0
        do i = 0, 4095
           scale = scale + abs(f(2*pi*i/4096))
        end do
        scale = scale*2*pi/4096
        m_chk = max(m_chk, abs(ref2-ref)/scale)
        m_rk_dp = max(m_rk_dp, abs(rk(100)-ref)/scale)
        do i = 1, nt; m_tr(i) = max(m_tr(i), abs(trap(ntr(i))-ref)/scale); end do
        do i = 1, ng; m_gl(i) = max(m_gl(i), abs(gauss(ngl(i))-ref)/scale); end do
        single = .true.
        m_rk_sp = max(m_rk_sp, abs(rk(100)-ref)/scale)
     end do
     close(10)
     write(*,'(a,a,i6,a)') '== ', trim(names(k)), npt, ' points; max |error| / int|f|'
     write(*,'(a,es9.2)') '   reference check (trap 2048 vs 4096):  ', m_chk
     write(*,'(a,es9.2)') '   current: RK niter=100 (400 calls), complex*8 : ', m_rk_sp
     write(*,'(a,es9.2)') '   RK niter=100 (400 calls), complex*16        : ', m_rk_dp
     do i = 1, nt; write(*,'(a,i4,a,es9.2)') '   trapezoid ', ntr(i), ' calls, complex*16 : ', m_tr(i); end do
     do i = 1, ng; write(*,'(a,i4,a,es9.2)') '   Gauss-Legendre ', ngl(i), ' calls        : ', m_gl(i); end do
  end do
end program quad
