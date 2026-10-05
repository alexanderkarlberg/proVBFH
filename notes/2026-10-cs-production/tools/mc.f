c Independent LO integration of the public VBF_HJJJ tree-level H+3j
c (setborn) with its own phase space, LHAPDF PDFs and alpha_s, the
c Born flavour list of the public init_processes, parton-level
c 1506.02660 cuts (all partons pT > 25, |y| < 4.5, pairwise dR > 0.4;
c tagging = two hardest: mjj > 600, dy > 4.5, opposite hemispheres),
c mu_R = mu_F = mu0(pT,H) = sqrt(mH/2 sqrt(mH^2/4 + pT,H^2)).
c Phase space: x1, x2 from tau, y; s -> H + Q (Q^2 flat), Q -> k6 + R
c (R^2 flat), R -> k4 + k5, isotropic two-body decays in rest frames.
      program mc
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      include 'PhysPars.h'
      integer*8 npts, ip, nacc
      integer j, k, mu, iflv
      real*8 p(0:3,nlegborn), born, bornjk(nlegborn,nlegborn),
     $     bmunu(0:3,0:3,nlegborn), bornsub(2)
      real*8 sb, rs, mh, tau, taumin, y, x1, x2, s, jac, w, sum, sum2
      real*8 q2, r2, q2max, r2max, pv, pq(0:3), pr(0:3), ph(0:3),
     $     k4(0:3), k5(0:3), k6(0:3), pd(0:3)
      real*8 f1(-6:6), f2(-6:6), mu2, alphasPDF, ptH, lum, rnd(9)
      real*8 pt(3), yy(3), mjj, dy, phi(3)
      integer iord(3)
      logical ok
      real*8 hbarc2
      parameter (hbarc2 = 0.3893793d9) ! GeV^2 pb
      external alphasPDF
      call init_flsttag
      call init_phys
      call initPDFsetByName('NNPDF30_nnlo_as_0118')
      call initPDF(0)
      mh = ph_Hmass
      sb = 13000d0**2
      npts = 300000000
      sum = 0; sum2 = 0; nacc = 0
      call srand(12345)
      do ip = 1, npts
         do k = 1, 9
            rnd(k) = rand()
         enddo
         jac = 1
c        tau, y
         taumin = mh**2/sb
         tau = taumin**rnd(1)
         jac = jac*tau*log(1/taumin)
         y = 0.5d0*log(tau)*(1-2*rnd(2))
         jac = jac*log(1/tau)
         x1 = sqrt(tau)*exp(y)
         x2 = sqrt(tau)*exp(-y)
         s = tau*sb
         rs = sqrt(s)
c        s -> H + Q, Q^2 in (0, (rs-mh)^2)
         q2max = (rs-mh)**2
         q2 = q2max*rnd(3)
         jac = jac*q2max/(2*pi)
         call twobody(rs, mh**2, q2, rnd(4), rnd(5), ph, pq, pv)
         jac = jac*pv/(16*pi**2*rs)*4*pi
c        Q -> k6 + R, R^2 in (0, Q^2)
         r2max = q2
         r2 = r2max*rnd(6)
         jac = jac*r2max/(2*pi)
         call twobody(sqrt(q2), 0d0, r2, rnd(7), rnd(8), k6, pr, pv)
         jac = jac*pv/(16*pi**2*sqrt(q2))*4*pi
         call boost(pq, k6)
         call boost(pq, pr)
c        R -> k4 + k5
         call twobody(sqrt(r2), 0d0, 0d0, rnd(9), rand(), k4, k5, pv)
         jac = jac*pv/(16*pi**2*sqrt(r2))*4*pi
         call boost(pr, k4)
         call boost(pr, k5)
c        lab frame: boost along z with the rapidity y
         call zboost(y, ph); call zboost(y, k4); call zboost(y, k5)
         call zboost(y, k6)
c        cuts on the three partons
         call ptyphi(k4, pt(1), yy(1), phi(1))
         call ptyphi(k5, pt(2), yy(2), phi(2))
         call ptyphi(k6, pt(3), yy(3), phi(3))
         ok = .true.
         do k = 1, 3
            if (pt(k).lt.25d0 .or. abs(yy(k)).gt.4.5d0) ok = .false.
         enddo
         if (.not.ok) cycle
         if (dr(yy(1),phi(1),yy(2),phi(2)).lt.0.4d0) cycle
         if (dr(yy(1),phi(1),yy(3),phi(3)).lt.0.4d0) cycle
         if (dr(yy(2),phi(2),yy(3),phi(3)).lt.0.4d0) cycle
c        tagging jets: the two hardest
         iord = (/1,2,3/)
         call sort3(pt, iord)
         if (iord(1).eq.1) then
            pd = k4
         elseif (iord(1).eq.2) then
            pd = k5
         else
            pd = k6
         endif
         call tagcuts(iord, k4, k5, k6, yy, ok)
         if (.not.ok) cycle
c        scales and PDFs
         ptH = sqrt(ph(1)**2+ph(2)**2)
         mu2 = mh/2*sqrt(mh**2/4+ptH**2)
         st_muren2 = mu2
         st_mufact2 = mu2
         st_alpha = alphasPDF(sqrt(mu2))
         call evolvePDF(x1, sqrt(mu2), f1)
         call evolvePDF(x2, sqrt(mu2), f2)
c        momenta in POWHEG order: 1, 2 incoming, 3 H, 4, 5, 6
         p(:,1) = 0; p(0,1) = x1*sqrt(sb)/2; p(3,1) = p(0,1)
         p(:,2) = 0; p(0,2) = x2*sqrt(sb)/2; p(3,2) = -p(0,2)
         p(:,3) = ph; p(:,4) = k4; p(:,5) = k5; p(:,6) = k6
         w = 0
         do iflv = 1, flst_nborn
            call setborn(p, flst_born(1,iflv), born, bornjk, bmunu,
     $           bornsub)
            lum = f1(flst_born(1,iflv))*f2(flst_born(2,iflv))/(x1*x2)
            w = w + born*lum
         enddo
c        flux and units; the generated final state is an ordered
c        triple (k4, k5, k6) matched to the flavour list's legs 4, 5, 6
         w = w*jac/(2*s)*hbarc2
         sum = sum + w
         sum2 = sum2 + w**2
         nacc = nacc + 1
      enddo
      sum = sum/npts
      sum2 = sum2/npts
      write(*,'(a,i10,a,es14.6,a,es12.4,a)') ' accepted', nacc,
     $     '  sigma(>=3 jets, VBF cuts) =', sum*1d3, ' +-',
     $     sqrt((sum2-sum**2)/npts)*1d3, ' fb'
      contains
      real*8 function dr(y1, p1, y2, p2)
      real*8 y1, p1, y2, p2, dp
      dp = abs(p1-p2)
      if (dp.gt.pi) dp = 2*pi-dp
      dr = sqrt((y1-y2)**2+dp**2)
      end function
      end

      subroutine twobody(m, m1s, m2s, r1, r2, p1, p2, pmod)
c     decay of a particle of mass m at rest into masses^2 m1s, m2s
      implicit none
      real*8 m, m1s, m2s, r1, r2, p1(0:3), p2(0:3), pmod, ct, st, ph
      real*8 lam
      lam = (m**2-m1s-m2s)**2-4*m1s*m2s
      pmod = sqrt(max(lam,0d0))/(2*m)
      ct = 1-2*r1
      st = sqrt(max(1-ct**2,0d0))
      ph = 2*3.141592653589793d0*r2
      p1(1) = pmod*st*cos(ph); p1(2) = pmod*st*sin(ph); p1(3) = pmod*ct
      p1(0) = sqrt(pmod**2+m1s)
      p2(1:3) = -p1(1:3)
      p2(0) = sqrt(pmod**2+m2s)
      end

      subroutine boost(q, p)
c     boost p from the rest frame of q to the frame where q is given
      implicit none
      real*8 q(0:3), p(0:3), m, e, pq, f
      m = sqrt(max(q(0)**2-q(1)**2-q(2)**2-q(3)**2, 1d-30))
      pq = p(1)*q(1)+p(2)*q(2)+p(3)*q(3)
      e = (q(0)*p(0)+pq)/m
      f = (p(0)+e)/(q(0)+m)
      p(1:3) = p(1:3)+f*q(1:3)
      p(0) = e
      end

      subroutine zboost(y, p)
      implicit none
      real*8 y, p(0:3), e, z
      e = cosh(y)*p(0)+sinh(y)*p(3)
      z = sinh(y)*p(0)+cosh(y)*p(3)
      p(0) = e; p(3) = z
      end

      subroutine ptyphi(p, pt, y, phi)
      implicit none
      real*8 p(0:3), pt, y, phi
      pt = sqrt(p(1)**2+p(2)**2)
      y = 0.5d0*log((p(0)+p(3))/(p(0)-p(3)))
      phi = atan2(p(2), p(1))
      end

      subroutine sort3(pt, iord)
c     indices of pt in decreasing order
      implicit none
      real*8 pt(3)
      integer iord(3), i, j, t
      do i = 1, 2
         do j = i+1, 3
            if (pt(iord(j)).gt.pt(iord(i))) then
               t = iord(i); iord(i) = iord(j); iord(j) = t
            endif
         enddo
      enddo
      end

      subroutine tagcuts(iord, k4, k5, k6, yy, ok)
      implicit none
      integer iord(3)
      real*8 k4(0:3), k5(0:3), k6(0:3), yy(3), a(0:3), b(0:3), m2
      logical ok
      call pick(iord(1), k4, k5, k6, a)
      call pick(iord(2), k4, k5, k6, b)
      m2 = (a(0)+b(0))**2-(a(1)+b(1))**2-(a(2)+b(2))**2
     $     -(a(3)+b(3))**2
      ok = m2.gt.600d0**2 .and. abs(yy(iord(1))-yy(iord(2))).gt.4.5d0
     $     .and. yy(iord(1))*yy(iord(2)).lt.0d0
      end

      subroutine pick(i, k4, k5, k6, a)
      implicit none
      integer i
      real*8 k4(0:3), k5(0:3), k6(0:3), a(0:3)
      if (i.eq.1) a = k4
      if (i.eq.2) a = k5
      if (i.eq.3) a = k6
      end
