! V4 (cs_4qint): pointwise check of the four-quark real matrix element of
! proVBFH-cs (compreal_hqqqq via setreal, i.e. what the program uses) against
! the independent Dirac-matrix reference ref.f90 (Weyl representation; only the
! couplings, masses and widths are read from the VBFNLO commons).
! For every four-quark entry of the real list (CC and NC, all tag patterns, all
! flavours incl. u/c mixing and crossings) and every phase-space point (x 24
! label permutations of the final-state partons):
!   cs_4qint 1: caller = N + 1/2 I  if the flipped diagram is physical, else N
!   cs_4qint 0: caller = N
!   cs_4qint 2: caller = 1/2 I  if the flipped diagram is physical, else 0
! N = |normal|^2, I = 2 colint Re(M_N M_F^*) (Fermi sign in M_F), colour-summed,
! times couplings, V propagators and the spin-colour average. The flag is read
! from powheg.input in the run directory.
! Usage: check_4qint [npt per mode (default 20)] [permutations (default 24)]
program check_4qint
  use ref
  implicit none
  double precision :: clr(4,5,-1:1), xm2(6), xmg(6), bb(6,6,6), vv(4,5), aa(4,5)
  common /bkopou/ clr, xm2, xmg, bb, vv, aa
  double complex :: ahvv(3,4,4), ahvvL(3,4,4)
  common /tensorhvv/ ahvv, ahvvL
  double precision :: mufsq(2,3), als(2,3), xir, xif, qmin(3)
  integer :: id_mur, id_muf
  common /cscales/ mufsq, als, xir, xif, qmin, id_mur, id_muf
  integer :: tag_factor
  common /cdoubletag/ tag_factor
  double precision, parameter :: pi = 3.141592653589793238d0
  integer, parameter :: maxe = 2400
  double precision :: powheginput
  external powheginput

  integer :: nreal, j, fl(7), tg(7), ipt, mode, pat, icc, i, k, npt, ne, e, cl
  integer :: efl(7,maxe), etg(7,maxe), ereal(maxe), epat(maxe), ecc(maxe), ecl(maxe)
  integer :: sigN(3,maxe), sigF(3,maxe), sn(3), sf(3)
  logical :: eflip(maxe), eident(maxe), cc, flipvalid, ident, ok, wpair, found
  integer :: mode4q
  double precision :: p(0:3,7), pp(0:3,7), amp2, Nn, Ff, Ii, pc2, ps2, norm, pred, scl, d
  double precision :: dmax(0:1,4,3), rmin(0:1,4,3), rmax(0:1,4,3), isum(0:1,3), nsum(0:1,3)
  integer :: cnt(0:1,4,3), nzero
  integer :: perm(4,24), iperm, ip, a, b, c, dd, nperm, nuse
  character(len=32) :: arg
  integer :: seed(8)

  seed = 12345
  call random_seed(put=seed)
  npt = 20
  call get_command_argument(1, arg)
  if (len_trim(arg) > 0) read(arg,*) npt
  nuse = 24
  call get_command_argument(2, arg)
  if (len_trim(arg) > 0) read(arg,*) nuse

  call hoppetStart(0.5d0,2)
  call hoppetSetCoupling(0.118d0,91.1876d0,2)
  call cs_init_vbfnlo()
  call cs_init_reals(nreal)
  call ref_init()
  do i = 3, 4
     do k = 2, 3
        cpl(i,k,-1) = clr(i,k,-1); cpl(i,k,1) = clr(i,k,1)
     end do
  end do
  mode4q = nint(powheginput('#cs_4qint'))
  if (mode4q < 0) mode4q = 1      ! absent: default 1
  write(*,'(a,i2)') ' V4: cs_4qint', mode4q

  ! the four-quark entries and their class (as in the w4q-spinor check)
  call gen(1, p)
  ne = 0
  do j = 1, nreal
     call cs_get_real(j, fl, tg)
     if (any(fl([1,2,4,5,6,7]) == 0)) cycle
     if (fl(3) == 240025) then
        cc = .true.
     else if (fl(3) == 230025) then
        cc = .false.
     else
        cycle
     end if
     ne = ne + 1
     efl(:,ne) = fl; etg(:,ne) = tg; ereal(ne) = j
     ecc(ne) = 0; if (cc) ecc(ne) = 1
     if (tg(1) >= 10) then
        epat(ne) = 3
     else if (tg(2) >= 10) then
        epat(ne) = 4
     else if (mod(tg(7),10) == 1) then
        epat(ne) = 1
     else
        epat(ne) = 2
     end if
     call ref_entry(fl, tg, tag_factor, p, cc, Nn, Ff, Ii, pc2, ps2, flipvalid, ident, ok, sn, sf, wpair)
     if (.not. ok) stop 'entry not classified'
     sigN(:,ne) = sn; sigF(:,ne) = sf; eflip(ne) = flipvalid; eident(ne) = ident
  end do
  ! class: 1 no physical flip, 2 flip, identical particles, 3 flip = normal of another entry, 4 flip not in the list
  do e = 1, ne
     if (.not. eflip(e)) then
        ecl(e) = 1
     else if (eident(e)) then
        ecl(e) = 2
     else
        found = .false.
        do k = 1, ne
           if (ecc(k) == ecc(e) .and. all(sigN(:,k) == sigF(:,e))) found = .true.
        end do
        ecl(e) = 4
        if (found) ecl(e) = 3
     end if
  end do
  write(*,'(a,i5)') ' four-quark entries (CC+NC):', ne
  do icc = 0, 1
     write(*,'(a,i2,a,4i6)') ' cc=', icc, '  entries in class 1..4:', (count(ecc(1:ne)==icc .and. ecl(1:ne)==cl), cl=1,4)
  end do
  if (count(ecl(1:ne) == 4) > 0) stop 'flip not covered by the list: 1/2 I is not the right weight'

  nperm = 0
  do a = 1, 4
     do b = 1, 4
        do c = 1, 4
           do dd = 1, 4
              if (a==b .or. a==c .or. a==dd .or. b==c .or. b==dd .or. c==dd) cycle
              nperm = nperm + 1
              perm(:,nperm) = [a,b,c,dd]
           end do
        end do
     end do
  end do
  nperm = min(nperm, nuse)

  dmax = 0; rmin = 1d300; rmax = -1d300; cnt = 0; nzero = 0; isum = 0; nsum = 0
  do mode = 1, 3
     do ipt = 1, npt
        call gen(mode, p)
        do iperm = 1, nperm
           pp = p
           do ip = 1, 4
              pp(:,3+ip) = p(:,3+perm(ip,iperm))
           end do
           do e = 1, ne
              fl = efl(:,e); tg = etg(:,e); cc = (ecc(e) == 1)
              call cs_real_me(ereal(e), pp, amp2)
              call ref_entry(fl, tg, tag_factor, pp, cc, Nn, Ff, Ii, pc2, ps2, flipvalid, ident, ok, sn, sf, wpair)
              k = 2; if (cc) k = 3
              norm = (16d0*pi**2)*als(1,1)*als(2,1)*abs(ahvv(1,k,k+merge(1,0,cc)))**2 / (36d0 * &
                   abs(cmplx(pc2-xm2(k), xmg(k), kind=8))**2 * abs(cmplx(ps2-xm2(k), xmg(k), kind=8))**2)
              pred = Nn
              if (mode4q == 2) pred = 0
              scl = Nn
              if (mode4q /= 0 .and. ecl(e) >= 2) then
                 pred = pred + 0.5d0*Ii
                 scl = Nn + 0.5d0*abs(Ii)
              end if
              pred = norm*pred; scl = norm*scl
              icc = ecc(e); pat = epat(e); cl = min(ecl(e),3)
              if (scl == 0d0) then
                 if (amp2 /= 0d0) write(*,'(a,i5,es14.5)') ' ZERO MISMATCH entry', ereal(e), amp2
                 nzero = nzero + 1
                 cycle
              end if
              d = abs(amp2 - pred)/scl
              dmax(icc,pat,cl) = max(dmax(icc,pat,cl), d)
              if (pred /= 0d0) then
                 rmin(icc,pat,cl) = min(rmin(icc,pat,cl), amp2/pred)
                 rmax(icc,pat,cl) = max(rmax(icc,pat,cl), amp2/pred)
              end if
              cnt(icc,pat,cl) = cnt(icc,pat,cl) + 1
              if (ecl(e) >= 2) then
                 isum(icc,mode) = isum(icc,mode) + abs(0.5d0*Ii*norm)
                 nsum(icc,mode) = nsum(icc,mode) + Nn*norm
              end if
           end do
        end do
     end do
  end do
  write(*,'(a)') ' cc pat class  evaluations   max |caller-pred|/(N+|I|/2)   min caller/pred       max caller/pred'
  do icc = 0, 1
     do pat = 1, 4
        do cl = 1, 3
           if (cnt(icc,pat,cl) > 0) write(*,'(i3,i4,i6,i12,3es22.14)') icc, pat, cl, cnt(icc,pat,cl), &
                dmax(icc,pat,cl), rmin(icc,pat,cl), rmax(icc,pat,cl)
        end do
     end do
  end do
  write(*,'(a,i8)') ' evaluations with N = I = 0 (and caller = 0):', nzero
  write(*,'(a)') ' size of the term: sum |I/2| / sum N over the flip entries, per mode (1 VBF-like, 2 generic, 3 p1506-like)'
  do mode = 1, 3
     write(*,'(a,i2,a,2es12.4)') '  mode', mode, '  NC, CC:', isum(0,mode)/nsum(0,mode), isum(1,mode)/nsum(1,mode)
  end do
  write(*,'(a,es12.4)') ' V4 overall max |caller-pred|/(N+|I|/2):', maxval(dmax)
  if (maxval(dmax) < 1d-10) then
     write(*,'(a)') ' V4 PASS'
  else
     write(*,'(a)') ' V4 FAIL'
  end if

contains

  subroutine gen(mode, pin)
    integer, intent(in) :: mode
    double precision, intent(out) :: pin(0:3,7)
    double precision :: r(20), pt(7), y(7), ph(7), mt, e, pz, sx, sy, ea, eb, mh
    integer :: k
    mh = 125d0
    do
       call random_number(r)
       if (mode == 1 .or. mode == 3) then
          pt(4) = 30d0 + 170d0*r(1)**2;  y(4) = 2.5d0 + 2d0*r(2)
          pt(5) = 30d0 + 170d0*r(3)**2;  y(5) = -2.5d0 - 2d0*r(4)
          pt(6) = 15d0 + 65d0*r(5)**2;   y(6) = -2.5d0 + 5d0*r(6)
          pt(3) = 10d0 + 140d0*r(7);      y(3) = -2d0 + 4d0*r(8)
          y(7) = -3d0 + 6d0*r(9)
       else
          pt(4) = 20d0 + 180d0*r(1);  y(4) = -3d0 + 6d0*r(2)
          pt(5) = 20d0 + 180d0*r(3);  y(5) = -3d0 + 6d0*r(4)
          pt(6) = 20d0 + 180d0*r(5);  y(6) = -3d0 + 6d0*r(6)
          pt(3) = 10d0 + 140d0*r(7);  y(3) = -2d0 + 4d0*r(8)
          y(7) = -3d0 + 6d0*r(9)
       end if
       ph(3) = 2d0*pi*r(10); ph(4) = 2d0*pi*r(11); ph(5) = 2d0*pi*r(12); ph(6) = 2d0*pi*r(13)
       sx = 0; sy = 0
       do k = 3, 6
          sx = sx + pt(k)*cos(ph(k)); sy = sy + pt(k)*sin(ph(k))
       end do
       pt(7) = sqrt(sx**2+sy**2)
       if (pt(7) < 12d0) cycle
       if (mode == 3) then
          if (pt(6) < 25d0 .or. pt(7) < 25d0 .or. abs(y(4)) > 4.5d0 .or. abs(y(5)) > 4.5d0) cycle
          if (abs(y(7)) > 4.5d0 .or. abs(y(6)) > 4.5d0) cycle
       end if
       ph(7) = atan2(-sy,-sx)
       do k = 3, 7
          if (k == 3) then
             mt = sqrt(mh**2+pt(3)**2)
          else
             mt = pt(k)
          end if
          pin(0,k) = mt*cosh(y(k)); pin(3,k) = mt*sinh(y(k))
          pin(1,k) = pt(k)*cos(ph(k)); pin(2,k) = pt(k)*sin(ph(k))
       end do
       if (mode == 3) then
          if (abs(y(4)-y(5)) < 4.5d0) cycle
          if ((pin(0,4)+pin(0,5))**2 - (pin(1,4)+pin(1,5))**2 - (pin(2,4)+pin(2,5))**2 - (pin(3,4)+pin(3,5))**2 < 600d0**2) cycle
       end if
       e = sum(pin(0,3:7)); pz = sum(pin(3,3:7))
       ea = 0.5d0*(e+pz); eb = 0.5d0*(e-pz)
       if (ea <= 0d0 .or. eb <= 0d0) cycle
       if (ea > 6500d0 .or. eb > 6500d0) cycle
       pin(:,1) = [ea, 0d0, 0d0, ea]
       pin(:,2) = [eb, 0d0, 0d0, -eb]
       exit
    end do
  end subroutine gen
end program check_4qint
