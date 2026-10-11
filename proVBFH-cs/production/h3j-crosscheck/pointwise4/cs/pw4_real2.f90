! Step 1 of the pointwise 4-parton comparison with VBFNLO 3.0 (12 Oct): real-emission matrix elements.
! Reads pw4_vbfnlo.dat (written by the instrumented m2s_qqh4j.F: for each 4-parton point the sum of
! the VBFNLO real matrix elements over all subprocess classes and all 24 assignments of the 4 outgoing
! momenta to the outgoing slots, keyed by the flavours (a, b | f_1..f_4) attached to the incoming
! partons and to the 4 outgoing momenta, without PDFs), and evaluates proVBFH-cs's real (setreal,
! all entries of the flavour list, same sum over the 24 leg assignments) at the same momenta.
program pw4_real2
  use types, only: dp
  use incl_vbfh, only: inclusive_init
  use incl_parameters, only: sqrts
  use phase_space, only: set_beams
  use cs_exclusive, only: cs_excl_setup, cs_excl_setup2
  implicit none
  integer, external :: cs_nreal
  integer, parameter :: nkey = 11**6
  real(dp), allocatable :: accv(:), accc(:)
  integer, allocatable :: tv(:), tc(:)
  character(len=2000) :: line
  character(len=100) :: arg
  integer :: ios, iev, nkin, nd, i, j, l, m, ip, nreal, ntv, ntc, nfull, nev, maxev, ndone
  integer :: perms(4,24), flav(7), tags(7), fk(6), idx, k
  real(dp) :: xi(2), als, mur2, p(0:3,6), pk(0:3,6), p7(0:3,7), h(0:3), me, f, tmp(4)
  logical :: lok(25), full, havep, haveA
  real(dp) :: sv, sc, rv, rc, r0, ppx
  real(dp), allocatable :: accu(:), pdfa(:,:)
  logical :: used, lk
  integer :: ia, ib
  real(dp) :: vab(-6:6,-6:6), cab(-6:6,-6:6), uab(-6:6,-6:6), pab(-6:6,-6:6)
  integer :: ncat, cat
  real(dp) :: catv(0:80), catc(0:80), catn(0:80)
  integer :: catcnt(0:80)
  real(dp) :: rmin, rmax
  integer :: unit_det
  real :: t0, t1
  call inclusive_init()
  write(0,*) 'init1'
  call cs_init_vbfnlo()
  write(0,*) 'init2'
  call cs_excl_setup()
  call set_beams(sqrts)
  call cs_excl_setup2(.false.)
  write(0,*) 'setup done'
  nreal = cs_nreal()
  write(0,*) 'nreal', nreal
  allocate(accv(0:nkey-1), accc(0:nkey-1), tv(nkey), tc(nkey))
  accv = 0; accc = 0
  allocate(accu(0:nkey-1)); accu = 0
  open(12, file='pw4_ab.txt', status='replace')
  ip = 0
  do i = 1,4; do j = 1,4; do l = 1,4; do m = 1,4
     if (i/=j .and. i/=l .and. i/=m .and. j/=l .and. j/=m .and. l/=m) then
        ip = ip + 1; perms(:,ip) = [i,j,l,m]
     endif
  enddo; enddo; enddo; enddo
  maxev = 1000; nfull = 3
  if (command_argument_count() >= 1) then; call get_command_argument(1,arg); read(arg,*) maxev; endif
  if (command_argument_count() >= 2) then; call get_command_argument(2,arg); read(arg,*) nfull; endif
  open(10, file='pw4_vbfnlo.dat', status='old')
  open(11, file='pw4_real_details.txt', status='replace')
  catv = 0; catc = 0; catn = 0; catcnt = 0
  ndone = 0; nev = 0; ntv = 0; haveA = .false.; havep = .false.
  call cpu_time(t0)
  do
     read(10,'(a)',iostat=ios) line
     if (ios /= 0) exit
     if (line(1:2) == 'EV') then
        read(line(3:),*) iev, nkin, xi, als, mur2, lok(1:nkin)
        ntv = 0; haveA = .false.; havep = .false.; pab = 0
     elseif (line(1:2) == 'P ') then
        read(line(3:),*) nd, lk, pk
        if (nd == 1) then; p = pk; havep = .true.; endif
     elseif (line(1:2) == 'A ') then
        read(line(3:),*) fk, f, ppx
        idx = (((((fk(1)+5)*11+(fk(2)+5))*11+(fk(3)+5))*11+(fk(4)+5))*11+(fk(5)+5))*11+(fk(6)+5)
        if (accv(idx) == 0) then; ntv = ntv + 1; tv(ntv) = idx; endif
        accv(idx) = accv(idx) + f
        pab(fk(1),fk(2)) = ppx
        haveA = .true.
     elseif (line(1:3) == 'END') then
        if (haveA .and. havep .and. nev < maxev .and. ntv > 5000) then
           nev = nev + 1
           full = nev <= nfull
           h = p(:,1) + p(:,2) - p(:,3) - p(:,4) - p(:,5) - p(:,6)
           ntc = 0
           do j = 1, nreal
              call cs_get_real(j, flav, tags)
              do ip = 1, 24
                 fk(1) = flav(1); fk(2) = flav(2)
                 do l = 1, 4
                    fk(2+perms(l,ip)) = flav(3+l)
                 enddo
                 idx = (((((fk(1)+5)*11+(fk(2)+5))*11+(fk(3)+5))*11+(fk(4)+5))*11+(fk(5)+5))*11+(fk(6)+5)
                 if (.not. full .and. accv(idx) == 0) cycle
                 p7(:,1) = p(:,1); p7(:,2) = p(:,2); p7(:,3) = h
                 do l = 1, 4
                    p7(:,3+l) = p(:,2+perms(l,ip))
                 enddo
                 call cs_real_me(j, p7, me)
                 if (me == 0) cycle
                 used = entry_used(tags)
                 if (used) accu(idx) = accu(idx) + me
                 if (accc(idx) == 0) then; ntc = ntc + 1; tc(ntc) = idx; endif
                 accc(idx) = accc(idx) + me
              enddo
           enddo
           ! compare
           f = als**2/(xi(1)*xi(2))
           sv = 0; sc = 0
           do j = 1, ntv; sv = sv + accv(tv(j)); enddo
           do j = 1, ntc; sc = sc + accc(tc(j)); enddo
           write(6,'(a,i4,a,i6,a,es12.4,a,es12.4,a,es12.4,a,i6,a,i6)') 'event', iev, ' keys vb', ntv, ' sumvb', sv, ' sumcs*f', sc*f, &
                & ' ratio', sv/(sc*f), ' keys cs', ntc, ' nreal', nreal
           call flush(6)
           ! per key details
           do j = 1, ntv
              idx = tv(j)
              rv = accv(idx); rc = accc(idx)*f
              call decode(idx, fk)
              write(11,'(i5,6i3,2es16.8,es14.6)') iev, fk, rv, rc, merge(rv/rc, 0.0_dp, rc /= 0)
           enddo
           do j = 1, ntc
              idx = tc(j)
              if (accv(idx) == 0) then
                 call decode(idx, fk)
                 write(11,'(i5,6i3,2es16.8,a)') iev, fk, 0.0_dp, accc(idx)*f, '  CS-only'
              endif
           enddo
           vab = 0; cab = 0; uab = 0
           do j = 1, ntv
              call decode(tv(j), fk); vab(fk(1),fk(2)) = vab(fk(1),fk(2)) + accv(tv(j))
           enddo
           do j = 1, ntc
              call decode(tc(j), fk); cab(fk(1),fk(2)) = cab(fk(1),fk(2)) + accc(tc(j))*f
              uab(fk(1),fk(2)) = uab(fk(1),fk(2)) + accu(tc(j))*f
           enddo
           do ia = -6,6; do ib = -6,6
              if (vab(ia,ib) /= 0 .or. cab(ia,ib) /= 0) write(12,'(i5,2i3,4es16.8)') iev, ia, ib, vab(ia,ib), cab(ia,ib), uab(ia,ib), pab(ia,ib)
           enddo; enddo
           call flush(12)
           ! reset
           do j = 1, ntc; accc(tc(j)) = 0; accu(tc(j)) = 0; enddo
        endif
        do j = 1, ntv; accv(tv(j)) = 0; enddo
        ntv = 0
     endif
  enddo
  call cpu_time(t1)
  write(6,'(a,i5,a,f9.1,a)') 'done events', nev, '  cpu', t1 - t0, ' s'
contains
  logical function entry_used(tg)
    integer, intent(in) :: tg(7)
    integer :: st
    st = sum(tg)
    entry_used = (st == 8 .or. st == 28 .or. st == 10 .or. st == 30)
    if (tg(1) == 1 .and. tg(2) == 2 .and. tg(4) == 1 .and. tg(5) == 2) then
       if ((tg(6) == 1 .and. tg(7) == 2) .or. (tg(6) == 2 .and. tg(7) == 1)) entry_used = .true.
    endif
  end function entry_used
  subroutine decode(idx, fk)
    integer, intent(in) :: idx
    integer, intent(out) :: fk(6)
    integer :: r, i
    r = idx
    do i = 6, 1, -1
       fk(i) = mod(r,11) - 5
       r = r/11
    enddo
  end subroutine decode
end program
