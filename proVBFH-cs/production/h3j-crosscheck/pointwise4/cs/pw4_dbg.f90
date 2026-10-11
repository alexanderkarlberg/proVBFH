! debug: contributions of the individual real entries to one flavour key at one dumped event
program pw4_dbg
  use types, only: dp
  use incl_vbfh, only: inclusive_init
  use incl_parameters, only: sqrts
  use phase_space, only: set_beams
  use cs_exclusive, only: cs_excl_setup, cs_excl_setup2
  implicit none
  integer, external :: cs_nreal
  character(len=2000) :: line
  character(len=100) :: arg
  integer :: ios, iev, nkin, nd, i, j, l, m, ip, nreal, tev, ka, kb, kout(4), perms(4,24), flav(7), tags(7), fk(6), ktag
  real(dp) :: xi(2), als, mur2, p(0:3,6), pk(0:3,6), p7(0:3,7), h(0:3), me, f
  logical :: lok(25), found
  integer :: nsel(4), nmatch
  call inclusive_init()
  call cs_init_vbfnlo()
  call cs_excl_setup()
  call set_beams(sqrts)
  call cs_excl_setup2(.false.)
  nreal = cs_nreal()
  ip = 0
  do i = 1,4; do j = 1,4; do l = 1,4; do m = 1,4
     if (i/=j .and. i/=l .and. i/=m .and. j/=l .and. j/=m .and. l/=m) then
        ip = ip + 1; perms(:,ip) = [i,j,l,m]
     endif
  enddo; enddo; enddo; enddo
  call get_command_argument(1,arg); read(arg,*) tev
  call get_command_argument(2,arg); read(arg,*) ka
  call get_command_argument(3,arg); read(arg,*) kb
  do i = 1,4
     call get_command_argument(3+i,arg); read(arg,*) kout(i)
  enddo
  open(10, file='pw4_vbfnlo.dat', status='old')
  found = .false.
  do
     read(10,'(a)',iostat=ios) line
     if (ios /= 0) exit
     if (line(1:2) == 'EV') then
        read(line(3:),*) iev, nkin, xi, als, mur2
        found = (iev == tev)
     elseif (line(1:2) == 'P ' .and. found) then
        read(line(3:),*) nd, lok(1), pk
        if (nd == 1) then
           p = pk
           exit
        endif
     endif
  enddo
  h = p(:,1) + p(:,2) - p(:,3) - p(:,4) - p(:,5) - p(:,6)
  f = als**2/(xi(1)*xi(2))
  write(6,'(a,i5,a,2es12.4)') 'event', tev, ' xi', xi
  do l = 1,6
     write(6,'(a,i2,4es14.5)') ' p', l, p(:,l)
  enddo
  nmatch = 0
  do j = 1, nreal
     call cs_get_real(j, flav, tags)
     if (flav(1) /= ka .or. flav(2) /= kb) cycle
     do ip = 1, 24
        fk(1) = flav(1); fk(2) = flav(2)
        do l = 1, 4
           fk(2+perms(l,ip)) = flav(3+l)
        enddo
        if (.not. same_multiset(fk(3:6), kout)) cycle
        p7(:,1) = p(:,1); p7(:,2) = p(:,2); p7(:,3) = h
        do l = 1, 4
           p7(:,3+l) = p(:,2+perms(l,ip))
        enddo
        call cs_real_me(j, p7, me)
        write(6,'(a,i5,a,7i4,a,7i3,a,4i3,a,es13.5)') ' j', j, ' flav', flav, ' tags', tags, ' key', fk(3:6), ' CS*f', me*f
     enddo
  enddo
contains
  logical function same_multiset(a, b)
    integer, intent(in) :: a(4), b(4)
    integer :: x(4), y(4), i, j, t
    x = a; y = b
    do i = 1,3; do j = 1,4-i
       if (x(j) > x(j+1)) then; t = x(j); x(j) = x(j+1); x(j+1) = t; endif
       if (y(j) > y(j+1)) then; t = y(j); y(j) = y(j+1); y(j+1) = t; endif
    enddo; enddo
    same_multiset = all(x == y)
  end function
end program
