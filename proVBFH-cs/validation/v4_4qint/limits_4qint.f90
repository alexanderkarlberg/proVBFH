! V4 (cs_4qint), single-unresolved limits of the identical-quark interference.
! For every four-quark entry with a physical flipped diagram (classes 2, 3 of
! check_4qint) and every single-unresolved limit of the H+4-parton (6-quark)
! real (final-final collinear i || j, initial-final collinear, soft final quark),
! the scaling of N = |normal|^2 and I/2 (the term added by cs_4qint) is measured
! with the reference ref.f90:  |X| ~ s^(-a) (collinear, s = 2 p_i.p_j) or
! |X| ~ E^(-a) (soft quark of energy E). Max over 4 azimuths at each step.
! Integrable: a < 1 (collinear, dPS ~ ds), a < 2 (soft, dPS ~ E dE).
! The caller (cs_4qint 1, from powheg.input) is checked against N + I/2 along
! the way. a_F: the flipped square F (= the normal square of the partner entry).
! |I/2| <= |colint|/colsq sqrt(N F) (Cauchy-Schwarz) is checked as well.
program limits_4qint
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
  integer, parameter :: maxe = 2400, nst = 4, nphi = 4, nbase = 2
  double precision, parameter :: th(nst) = [1d-2, 1d-3, 1d-4, 1d-5]
  integer :: nreal, j, fl(7), tg(7), i, k, ne, e, ib, ist, iph, ltype, li, lj, icc
  integer :: efl(7,maxe), etg(7,maxe), ereal(maxe), ecc(maxe), sn(3), sf(3)
  logical :: cc, flipvalid, ident, ok, wpair, sing
  double precision :: p0(0:3,7), p(0:3,7), Nn, Ff, Ii, pc2, ps2, amp2, norm, x
  double precision :: xN(nst), xI(nst), xF(nst), aN, aI, aF, prec, csmax
  ! [cc, limit type 1 FF coll, 2 IF coll, 3 soft]: max a_I over all, max a_I where N is singular, counts
  double precision :: aImax(0:1,3), aNmax(0:1,3), aFmax(0:1,3), precmax
  integer :: nlim(0:1,3), nsing(0:1,3)
  double precision :: powheginput
  external powheginput
  integer :: seed(8)

  seed = 777
  call random_seed(put=seed)
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
  if (nint(powheginput('#cs_4qint')) == 0) stop 'run with cs_4qint 1'

  call gen(p0)
  ne = 0
  do j = 1, nreal
     call cs_get_real(j, fl, tg)
     if (any(fl([1,2,4,5,6,7]) == 0)) cycle
     if (fl(3) /= 240025 .and. fl(3) /= 230025) cycle
     cc = (fl(3) == 240025)
     call ref_entry(fl, tg, tag_factor, p0, cc, Nn, Ff, Ii, pc2, ps2, flipvalid, ident, ok, sn, sf, wpair)
     if (.not. flipvalid) cycle
     ne = ne + 1
     efl(:,ne) = fl; etg(:,ne) = tg; ereal(ne) = j; ecc(ne) = merge(1,0,cc)
  end do
  write(*,'(a,i5)') ' four-quark entries with a physical flip:', ne

  aImax = -1d300; aNmax = -1d300; aFmax = -1d300; csmax = 0; nlim = 0; nsing = 0; precmax = 0
  write(*,'(a)') ' entry  cc lim  i  j     a_N      a_I      a_F   (singular N: a_N > 0.9)'
  do ib = 1, nbase
     call gen(p0)
     do e = 1, ne
        fl = efl(:,e); tg = etg(:,e); cc = (ecc(e) == 1); icc = ecc(e)
        do ltype = 1, 3
           do lj = 4, 7
              do li = 1, 7
                 if (ltype == 3 .and. li /= 1) cycle      ! soft: one entry per final quark
                 if (ltype == 1 .and. (li < 4 .or. li >= lj)) cycle
                 if (ltype == 2 .and. li > 2) cycle
                 xN = 0; xI = 0; xF = 0
                 do ist = 1, nst
                    do iph = 1, nphi
                       call limitpoint(ltype, li, lj, th(ist), 2d0*pi*(iph-0.5d0)/nphi + 0.3d0*ib, p0, p, ok)
                       if (.not. ok) exit
                       ! (cs_real_me first: it sets the couplings als used in norm)
                       call cs_real_me(ereal(e), p, amp2)
                       call ref_entry(fl, tg, tag_factor, p, cc, Nn, Ff, Ii, pc2, ps2, flipvalid, ident, ok, sn, sf, wpair)
                       k = 2; if (cc) k = 3
                       norm = (16d0*pi**2)*als(1,1)*als(2,1)*abs(ahvv(1,k,k+merge(1,0,cc)))**2 / (36d0 * &
                            abs(cmplx(pc2-xm2(k), xmg(k), kind=8))**2 * abs(cmplx(ps2-xm2(k), xmg(k), kind=8))**2)
                       xN(ist) = max(xN(ist), abs(norm*Nn))
                       xI(ist) = max(xI(ist), abs(norm*0.5d0*Ii))
                       xF(ist) = max(xF(ist), abs(norm*Ff))
                       ! Cauchy-Schwarz: |I/2| <= |colint|/colsq sqrt(N F)
                       if (Nn*Ff > 0d0) csmax = max(csmax, abs(0.5d0*Ii)/(abs(colint)/colsq*sqrt(Nn*Ff)))
                       if (Nn > 0d0) then
                          prec = abs(amp2 - norm*(Nn+0.5d0*Ii))/(norm*(Nn+0.5d0*abs(Ii)))
                          precmax = max(precmax, prec)
                       end if
                    end do
                    if (.not. ok) exit
                 end do
                 if (.not. ok) cycle
                 if (xN(nst) <= 0d0 .or. xN(nst-1) <= 0d0) cycle
                 if (ltype == 3) then
                    x = log(th(nst-1)/th(nst))            ! E ~ th
                 else
                    x = 2d0*log(th(nst-1)/th(nst))        ! s ~ th^2
                 end if
                 aN = log(xN(nst)/xN(nst-1))/x
                 aI = -1d0
                 if (xI(nst) > 0d0 .and. xI(nst-1) > 0d0) aI = log(xI(nst)/xI(nst-1))/x
                 aF = -1d0
                 if (xF(nst) > 0d0 .and. xF(nst-1) > 0d0) aF = log(xF(nst)/xF(nst-1))/x
                 aFmax(icc,ltype) = max(aFmax(icc,ltype), aF)
                 nlim(icc,ltype) = nlim(icc,ltype) + 1
                 aImax(icc,ltype) = max(aImax(icc,ltype), aI)
                 aNmax(icc,ltype) = max(aNmax(icc,ltype), aN)
                 sing = aN > 0.9d0
                 if (sing) nsing(icc,ltype) = nsing(icc,ltype) + 1
                 if (ib == 1 .and. (aI > 0.6d0 .or. (sing .and. e <= 400 .and. mod(e,37) == 0))) &
                      write(*,'(i6,i4,i4,2i3,3f9.4)') ereal(e), icc, ltype, li, lj, aN, aI, aF
              end do
           end do
        end do
     end do
  end do
  write(*,'(a)') ' ---- summary: exponents a, |X| ~ s^-a (collinear) or E^-a (soft), from the last two steps'
  write(*,'(a)') ' cc  limit           cases  N singular   max a_N   max a_I   max a_F   (integrable: a < 1 coll., < 2 soft)'
  do icc = 0, 1
     write(*,'(i3,a,2i10,3f10.4)') icc, '  FF collinear', nlim(icc,1), nsing(icc,1), aNmax(icc,1), aImax(icc,1), aFmax(icc,1)
     write(*,'(i3,a,2i10,3f10.4)') icc, '  IF collinear', nlim(icc,2), nsing(icc,2), aNmax(icc,2), aImax(icc,2), aFmax(icc,2)
     write(*,'(i3,a,2i10,3f10.4)') icc, '  soft quark  ', nlim(icc,3), nsing(icc,3), aNmax(icc,3), aImax(icc,3), aFmax(icc,3)
  end do
  write(*,'(a,f8.4)') ' max |I/2| / (|colint|/colsq sqrt(N F)) (Cauchy-Schwarz bound, must be <= 1):', csmax
  write(*,'(a,es10.2)') ' caller (cs_4qint 1) vs N + I/2 along the limits, max rel. difference:', precmax

contains

  ! final momenta of a point: legs 4..7 massless (pt, y, phi), Higgs (3) takes the
  ! transverse recoil at fixed rapidity, incoming from E +- pz
  subroutine gen(pin)
    double precision, intent(out) :: pin(0:3,7)
    double precision :: r(12), pt, y, ph
    integer :: k
    logical :: okk
    do
       call random_number(r)
       do k = 4, 7
          pt = 20d0 + 100d0*r(k-3); y = -2.5d0 + 5d0*r(k+1); ph = 2d0*pi*r(k+5)
          pin(:,k) = [pt*cosh(y), pt*cos(ph), pt*sin(ph), pt*sinh(y)]
       end do
       pin(3,3) = 0d0
       call balance(pin, okk)
       if (okk) exit
    end do
  end subroutine gen

  subroutine balance(pin, okk)
    double precision, intent(inout) :: pin(0:3,7)
    logical, intent(out) :: okk
    double precision :: px, py, mt, e, pz, ea, eb, yh
    yh = 0.2d0
    px = -sum(pin(1,4:7)); py = -sum(pin(2,4:7))
    mt = sqrt(125d0**2 + px**2 + py**2)
    pin(:,3) = [mt*cosh(yh), px, py, mt*sinh(yh)]
    e = sum(pin(0,3:7)); pz = sum(pin(3,3:7))
    ea = 0.5d0*(e+pz); eb = 0.5d0*(e-pz)
    okk = ea > 0d0 .and. eb > 0d0 .and. ea < 6500d0 .and. eb < 6500d0
    pin(:,1) = [ea, 0d0, 0d0, ea]
    pin(:,2) = [eb, 0d0, 0d0, -eb]
  end subroutine balance

  ! the limit: 1 final j collinear to final i (angle t), 2 final j collinear to
  ! incoming i (angle t), 3 final j soft (energy x t); phi: azimuth of the approach
  subroutine limitpoint(ltype, i, j, t, phi, pb, pout, okk)
    integer, intent(in) :: ltype, i, j
    double precision, intent(in) :: t, phi, pb(0:3,7)
    double precision, intent(out) :: pout(0:3,7)
    logical, intent(out) :: okk
    double precision :: n(3), u(3), v(3), d(3), ej
    pout = pb
    ej = pb(0,j)
    if (ltype == 3) then
       pout(:,j) = t*pb(:,j)
    else
       if (ltype == 1) then
          n = pb(1:3,i)/pb(0,i)
       else
          n = [0d0, 0d0, merge(1d0,-1d0,i==1)]
       end if
       ! orthonormal u, v perpendicular to n
       u = [n(2), -n(1), 0d0]
       if (sqrt(sum(u**2)) < 1d-3) u = [1d0, 0d0, 0d0]
       u = u - sum(u*n)*n; u = u/sqrt(sum(u**2))
       v = [n(2)*u(3)-n(3)*u(2), n(3)*u(1)-n(1)*u(3), n(1)*u(2)-n(2)*u(1)]
       d = cos(t)*n + sin(t)*(cos(phi)*u + sin(phi)*v)
       ej = min(ej, 200d0)
       pout(:,j) = [ej, ej*d(1), ej*d(2), ej*d(3)]
    end if
    call balance(pout, okk)
  end subroutine limitpoint
end program limits_4qint
