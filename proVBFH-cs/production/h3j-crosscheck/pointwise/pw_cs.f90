! Pointwise comparison with VBFNLO 3.0: reads pw_vbfnlo.dat (dumped by the instrumented
! m2s_qqh3j.F, see README), evaluates the proVBFH-cs Born, virtual and I operator at the same
! momenta and prints the ratios.  Scratch driver; links the proVBFH-cs objects except cs_main.
program pw_cs
  use types, only: dp
  use incl_vbfh, only: inclusive_init
  use incl_parameters, only: sqrts
  use phase_space, only: set_beams
  use cs_nlo2, only: nlo2_ifin, nlo2_ifin_dis, nlo2_kp, nlo2_kp_dis
  use cs_chan, only: kp_beam
  use cs_exclusive, only: cs_excl_setup, cs_excl_setup2
  implicit none
  real(dp), parameter :: twopi = 6.283185307179586476925286766559005768394_dp
  character(len=400) :: line
  integer :: ios, gs, nlo, fs(4), p2d(5), i, ch, jc, ld, lo, k, d, ip, io
  real(dp) :: p(0:3,5), h(0:3), xi(2), als1, als2, mu2, fac(3), arr(6,6), p6(0:3,6)
  real(dp) :: m1, m2, bcs, virt, ifl, ifo, rcs, rvb, rvbn, born_vb, fl1(6), fl2(6)
  integer :: bflav(6), a(2), b(2), sig(2), lines_of(2), f1, f2, nrow
  character(len=4) :: tag
  integer :: npt, jcl(2), idx, jx
  integer :: ngp
  real(dp) :: skp(-5:5,4), f0(-5:5,2), fq(-6:6), fg, fgq, fqo(-6:6), kpv, kpc, kpvo, kpco, asotpi, ef, muf
  real(dp), parameter :: CF = 4.0_dp/3, CA = 3
  call inclusive_init()
  call cs_init_vbfnlo()
  call cs_excl_setup()
  call set_beams(sqrts)
  call cs_excl_setup2(.false.)
  open(10, file='pw_vbfnlo.dat', status='old')
  npt = 0
  do
     read(10,'(a)',iostat=ios) line
     if (ios /= 0) exit
     if (line(1:2) /= 'PT') cycle
     read(line(3:),*) gs, nlo, fs, p2d
     read(10,'(a)') line; read(line(3:),*) xi, als1, als2, mu2
     do i = 1, 5
        read(10,'(a)') line; read(line(3:),*) p(:,i)
     enddo
     read(10,'(a)') line; read(line(2:),*) h
     read(10,'(a)') line; read(line(4:),*) fac
     do i = 1, 6
        read(10,'(a)') line; read(line(5:),*) arr(:,i)
     enddo
     do i = 1, 4
        read(10,'(a)') line; read(line(3:),*) idx, skp(:,i)
     enddo
     do i = 1, 2
        read(10,'(a)') line; read(line(3:),*) idx, f0(:,i)
     enddo
     if (gs == -1) then
        call gluon_point()
        cycle
     endif
     if (gs /= 1) cycle
     npt = npt + 1
     ! cs momenta: 1,2 incoming, 3 H, 4,5 outgoing of lines 1,2, 6 gluon
     ! line 21 = diag labels 1,2 ; line 43 = diag 3,4.  Outgoing partner of beam l:
     ! the outgoing parton (phys 3 or 4) with the same line as beam l.
     do d = 1, 2   ! beam d
        ip = d
        do io = 3, 4
           if ((p2d(ip) <= 2) .eqv. (p2d(io) <= 2)) exit
        enddo
        lines_of(d) = io
        jcl(d) = merge(2, 3, p2d(ip) <= 2)   ! jc of beam d's line (2: line 21, 3: line 43)
     enddo
     h = p(:,1) + p(:,2) - p(:,3) - p(:,4) - p(:,5)   ! the dumped v is not the Higgs
     p6(:,1) = p(:,1); p6(:,2) = p(:,2); p6(:,3) = h
     p6(:,4) = p(:,lines_of(1)); p6(:,5) = p(:,lines_of(2)); p6(:,6) = p(:,5)
     sig(1) = fs(p2d(1)); sig(2) = fs(p2d(2))
     do ch = 1, 6
        ! ch 1..4 NC: uucc uuss ddcc ddss (line-1 type, line-2 type); 5,6 CC udsc ducs
        select case (ch)
        case (1); f1 = 2; f2 = 4
        case (2); f1 = 2; f2 = 3
        case (3); f1 = 1; f2 = 4
        case (4); f1 = 1; f2 = 3
        case (5); f1 = 2; f2 = 3
        case (6); f1 = 1
        end select
        if (ch <= 4) then
           a = [sig(1)*f1, sig(2)*f2]; b = a
        else
           if (ch == 5) then
              f1 = 2; f2 = merge(3, 4, sig(1)*sig(2) == 1)
              a = [sig(1)*f1, sig(2)*f2]; b = [sig(1)*1, sig(2)*merge(4, 3, sig(1)*sig(2) == 1)]
              ! line 1 u -> d, line 2 d -> u (quarks) or ubar-type -> dbar-type: s->c / c->s
              if (sig(1)*sig(2) == -1) then
                 a(2) = sig(2)*4; b(2) = sig(2)*3
              else
                 a(2) = sig(2)*3; b(2) = sig(2)*4
              endif
           else
              f1 = 1
              a = [sig(1)*1, 0]; b = [sig(1)*2, 0]
              if (sig(1)*sig(2) == -1) then
                 a(2) = sig(2)*3; b(2) = sig(2)*4
              else
                 a(2) = sig(2)*4; b(2) = sig(2)*3
              endif
           endif
        endif
        bflav = [a(1), a(2), 25, b(1), b(2), 0]
        call cs_hjjj_lines(p6, bflav, m1, m2)
        do ld = 1, 2   ! cs line radiating
           if (ch == 1 .or. ch == 4) then
              lo = 3 - ld
              kp_beam = ld
              call nlo2_kp(xi(ld), p6(:,ld), p6(:,3+ld), p6(:,6), 125.0_dp, fq, fg, fgq)
              kp_beam = lo
              call nlo2_kp_dis(xi(lo), p6(:,lo), p6(:,3+lo), 125.0_dp, fqo)
              asotpi = als1/twopi
              if (ld == 1) then
                 kpv = ((CF - CA/2)/CF*skp(a(1),1) + CA/(2*CF)*skp(a(1),2))/asotpi
                 kpvo = skp(a(2),3)/asotpi
              else
                 kpv = ((CF - CA/2)/CF*skp(a(2),3) + CA/(2*CF)*skp(a(2),4))/asotpi
                 kpvo = skp(a(1),1)/asotpi
              endif
              write(6,'(a,i4,2i2,2x,a,3es13.5,2x,a,3es13.5)') 'KPR', npt, ch, ld, 'L: vb cs ratio', kpv, fq(a(ld)), kpv/fq(a(ld)), &
                   & 'O: vb cs ratio', kpvo, fqo(a(lo)), kpvo/fqo(a(lo))
           endif
           jc = jcl(ld)
           bcs = merge(m1, m2, ld == 1)
           call cs_hjjj_virt_line(p6, bflav, ld, mu2, virt)
           lo = 3 - ld
           ifl = nlo2_ifin(1, p6(:,ld), p6(:,3+ld), p6(:,6), mu2)
           ifo = nlo2_ifin_dis(p6(:,lo), p6(:,3+lo), mu2)
           rcs = virt/bcs + (ifl + ifo)/twopi
           ! VBFNLO arrays: columns (b2,b3,v2,v3,u2,u3)
           born_vb = arr(merge(1, 2, jc == 2), ch)
           rvbn = (arr(merge(5, 6, jc == 2), ch)/born_vb - 1)/als1
           write(6,'(a,i3,i2,i2,2es14.6,3es13.5,es14.6,2x,es13.5)') 'PW', npt, ch, ld, bcs, born_vb, &
                & virt/bcs, (ifl + ifo)/twopi, rcs, rvbn, rvbn - rcs
        enddo
     enddo
  enddo
contains
  ! gluon-initiated points (NC only)
  subroutine gluon_point()
    integer :: jcg, pd1, pd2, od1, od2, ig, ib, dpp(5), iq, iqb, oo, tp, to, fp, fo, so, ch2, ap_
    real(dp) :: pq(0:3), pqb(0:3), pO(0:3), r1, r2, kpl, kpo, fgv
    do i = 1, 5
       dpp(p2d(i)) = i
    enddo
    if (fs(1) == -fs(2) .and. fs(3) == fs(4)) then
       jcg = 2; pd1 = 1; pd2 = 2; od1 = 3; od2 = 4
    elseif (fs(1) == fs(2) .and. fs(3) == -fs(4)) then
       jcg = 3; pd1 = 3; pd2 = 4; od1 = 1; od2 = 2
    else
       return
    endif
    ig = merge(1, 2, p2d(1) == 5)
    ib = 3 - ig
    h = p(:,1) + p(:,2) - p(:,3) - p(:,4) - p(:,5)
    if (fs(pd1) == 1) then
       iq = dpp(pd1); iqb = dpp(pd2)
    else
       iq = dpp(pd2); iqb = dpp(pd1)
    endif
    ! outgoing parton of the other line
    oo = dpp(od1); if (oo == ib) oo = dpp(od2)
    so = fs(p2d(ib))
    p6(:,ig) = p(:,ig); p6(:,ib) = p(:,ib); p6(:,3) = h
    p6(:,3+ig) = p(:,iq); p6(:,6) = p(:,iqb); p6(:,3+ib) = p(:,oo)
    do ch2 = 1, 4
       select case (ch2)
       case (1); tp = 2; to = 2    ! (type of line 21, line 43)
       case (2); tp = 2; to = 1
       case (3); tp = 1; to = 2
       case (4); tp = 1; to = 1
       end select
       ! ch2: (t21, t43) = up/up, up/down, down/up, down/down; pair line is 21 for jc=2
       if (jcg == 2) then
          fp = merge(2, 1, ch2 <= 2); fo = merge(4, 3, ch2 == 1 .or. ch2 == 3)
       else
          fp = merge(2, 1, ch2 == 1 .or. ch2 == 3); fo = merge(4, 3, ch2 <= 2)
       endif
       a(ig) = 0; a(ib) = so*fo
       b(ig) = fp; b(ib) = so*fo
       if (ig == 1) then
          bflav = [0, a(2), 25, fp, b(2), -fp]
       else
          bflav = [a(1), 0, 25, b(1), fp, -fp]
       endif
       call cs_hjjj_lines(p6, bflav, m1, m2)
       bcs = merge(m1, m2, ig == 1)
       call cs_hjjj_virt_line(p6, bflav, ig, mu2, virt)
       ifl = nlo2_ifin(2, p6(:,ig), p6(:,3+ig), p6(:,6), mu2)
       ifo = nlo2_ifin_dis(p6(:,ib), p6(:,3+ib), mu2)
       rcs = virt/bcs + (ifl + ifo)/twopi
       born_vb = arr(merge(1, 2, jcg == 2), ch2)
       rvbn = (arr(merge(5, 6, jcg == 2), ch2)/born_vb - 1)/als1
       kp_beam = ig
       call nlo2_kp(xi(ig), p6(:,ig), p6(:,3+ig), p6(:,6), 125.0_dp, fq, fg, fgq)
       kp_beam = ib
       call nlo2_kp_dis(xi(ib), p6(:,ib), p6(:,3+ib), 125.0_dp, fqo)
       asotpi = als1/twopi
       if (jcg == 2) then
          kpl = 0.5_dp*(skp(0,1) + skp(0,2))/asotpi; kpo = skp(a(ib),3)/asotpi
       else
          kpl = 0.5_dp*(skp(0,3) + skp(0,4))/asotpi; kpo = skp(a(ib),1)/asotpi
       endif
       write(6,'(a,i4,3i2,2es13.5,3es13.5,es10.2,2x,2es13.5,2x,2es13.5)') 'GPW', npt, ch2, ig, jcg, bcs, born_vb, virt/bcs, &
            & (ifl + ifo)/twopi, rcs, rcs - rvbn, rvbn, kpl/(fg*xi(ig)), kpo/(fqo(a(ib))*xi(ib)), kpl, fg*xi(ig)
    enddo
  end subroutine
end program
