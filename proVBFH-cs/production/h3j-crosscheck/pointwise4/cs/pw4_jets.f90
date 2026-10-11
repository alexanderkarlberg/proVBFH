! Step 3 of the pointwise 4-parton comparison (12 Oct): jet algorithm and cuts.
! For every configuration of pw4_vbfnlo.dat (the real event, nd = 1, and the 24 dipole-mapped
! three-parton configurations nd = 2..25 that VBFNLO's phase space provides) compare VBFNLO's
! accept/reject decision (lok(nd), Cuts of VBFNLO: anti-kt R = 0.4 with E-scheme, pt, y, VBF cuts)
! with proVBFH-cs's (cs_passes: fastjet anti-kt, buildjets, vbfcuts of analysis/p1506_analysis.f).
program pw4_jets
  implicit none
  logical, external :: cs_passes
  character(len=2000) :: line
  integer :: ios, iev, nkin, nd, i, n, ndis(25)
  double precision :: xi(2), als, mur2, p(0:3,6), pk(0:3,6), pn(0:3,7), h(0:3), mh2, e
  logical :: lok(25), lk, pass, valid, okp
  double precision :: ptj(4), yj(4), pj(0:3,4)
  integer :: njets
  integer :: cnt(0:1,0:1,2), nval(2), ndet, nprint
  integer :: cntnd(25,0:1,0:1)
  cnt = 0; cntnd = 0; nprint = 0
  open(10, file='pw4_vbfnlo.dat', status='old')
  open(11, file='pw4_jets_details.txt', status='replace')
  iev = 0
  do
     read(10,'(a)',iostat=ios) line
     if (ios /= 0) exit
     if (line(1:2) == 'EV') then
        read(line(3:),*) iev, nkin, xi, als, mur2, lok(1:nkin)
     elseif (line(1:2) == 'P ') then
        read(line(3:),*) nd, lk, pk
        p = pk
        if (nd == 1) then
           n = 7
           h = p(:,1) + p(:,2) - p(:,3) - p(:,4) - p(:,5) - p(:,6)
           pn(:,1) = p(:,1); pn(:,2) = p(:,2); pn(:,3) = h
           do i = 1, 4; pn(:,3+i) = p(:,2+i); enddo
           valid = .true.
        else
           n = 6
           h = p(:,1) + p(:,2) - p(:,3) - p(:,4) - p(:,5)
           pn(:,1) = p(:,1); pn(:,2) = p(:,2); pn(:,3) = h
           do i = 1, 3; pn(:,3+i) = p(:,2+i); enddo
           valid = .true.
           do i = 3, 5
              e = p(0,i)
              if (e <= 0d0) valid = .false.
              if (abs(p(0,i)**2 - p(1,i)**2 - p(2,i)**2 - p(3,i)**2) > 1d-6*max(e*e,1d0)) valid = .false.
           enddo
           mh2 = h(0)**2 - h(1)**2 - h(2)**2 - h(3)**2
           if (abs(mh2 - 125d0**2) > 1d0) valid = .false.
        endif
        if (.not. valid) cycle
        pass = cs_passes(n, pn, .false.)
        call cs_fill_phep(n, pn)
        call buildjets(pj, njets, ptj, yj)
        pass = pass .and. njets >= 3
        cnt(merge(1,0,lok(nd)), merge(1,0,pass), merge(1,2,nd == 1)) = cnt(merge(1,0,lok(nd)), merge(1,0,pass), merge(1,2,nd == 1)) + 1
        cntnd(nd, merge(1,0,lok(nd)), merge(1,0,pass)) = cntnd(nd, merge(1,0,lok(nd)), merge(1,0,pass)) + 1
        if (lok(nd) .neqv. pass) then
           write(11,'(a,i5,a,i3,a,l2,a,l2)') 'DISAGREE event', iev, ' nd', nd, ' VBFNLO', lok(nd), ' CS', pass
           call cs_debug_jets(n, pn)
           call flush(11)
        endif
     endif
  enddo
  write(6,'(a)') 'accept(VBFNLO, CS) counts: [VB=0,CS=0] [VB=0,CS=1] [VB=1,CS=0] [VB=1,CS=1]'
  write(6,'(a,4i8)') 'real event (nd=1):       ', cnt(0,0,1), cnt(0,1,1), cnt(1,0,1), cnt(1,1,1)
  write(6,'(a,4i8)') 'mapped configs (nd>=2):  ', cnt(0,0,2), cnt(0,1,2), cnt(1,0,2), cnt(1,1,2)
  do nd = 1, 25
     if (sum(cntnd(nd,:,:)) > 0) write(6,'(a,i3,4i8)') ' nd', nd, cntnd(nd,0,0), cntnd(nd,0,1), cntnd(nd,1,0), cntnd(nd,1,1)
  enddo
end program
