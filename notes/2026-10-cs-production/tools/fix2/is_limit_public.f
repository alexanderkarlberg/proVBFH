c Initial-state collinear limit of the NC four-quark real of the PUBLIC
c POWHEG-BOX-V2 VBF_HJJJ (svn r4135), for the bug report to its authors.
c
c s(1) c(2) -> H s(4) c(5) Q(6) Qbar(7), outgoing s || incoming s: the
c graphs with the Z on the Q Qbar pair factorise onto the gluon-initiated
c Born g(x p1) c(p2) -> H Q c Qbar with the splitting q -> q + g (P_gq):
c   R -> 8 pi as / (x 2 p1.k4) P_gq(x) B_g   (azimuthal average),
c i.e. c = <amp2>_phi 2 p1.k4 x / (16 pi^2 P_gq(x) B_g) -> 1, where amp2
c is setreal's |M|^2/(as/2pi) and B_g setborn's Born, both spin/colour
c averaged. Printed for Q = u and Q = d and a sequence of k_T.
c
c Kinematics: a fixed Born point; the real point from the Catani-Seymour
c initial-initial map (k4 = (1-x-b) p1 + b p2 + kT, the final state
c Lorentz-transformed from K~ = x p1 + p2 to K = p1 + p2 - k4), exact.
c Link against the public build's objects (without pwhg_main.o); run in
c a directory with its testrun powheg.input and vbfnlo.input.
      program is_limit_public
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      include 'pwhg_flst_2.h'
      include 'PhysPars.h'
      integer HWW, HZZ, iq, iphi, nphi, ikt, mu, j
      parameter (nphi = 16)
      real * 8 pb(0:3,nlegborn), pr(0:3,nlegreal), born, bornjk(nlegborn,
     $     nlegborn), bmunu(0:3,0:3,nlegborn), bornsub(2), amp2, avg
      real * 8 x, kt, phi, p1(0:3), p2(0:3), a(0:3), b(0:3), k4(0:3),
     $     kk(0:3), kb(0:3), ftil(0:3,4), s12, beta, alpha, pgq, c,
     $     kts(4), dot
      integer bflav(nlegborn), rflav(nlegreal)
      data kts /1d0, 1d-1, 1d-2, 1d-3/
      external dot

c     the public init_phys (PDF set, couplings, flavour list, regions);
c     run with softtest 0 and colltest 0 in powheg.input
      call init_flsttag
      call init_phys
      call particle_identif(HWW, HZZ)
      st_alpha = 0.118d0
      st_muren2 = 100d0**2
      st_mufact2 = 100d0**2
      alr_tag = 0
      realequiv_tag = 0

c the Born point: final Q(6), c(5), Qbar(7) massless, H on shell
      call mom(60d0, 1.2d0, 0.3d0, ftil(0,2))
      call mom(55d0, -2.6d0, 2.9d0, ftil(0,3))
      call mom(40d0, 0.4d0, 4.6d0, ftil(0,4))
      ftil(1,1) = -(ftil(1,2) + ftil(1,3) + ftil(1,4))
      ftil(2,1) = -(ftil(2,2) + ftil(2,3) + ftil(2,4))
      ftil(3,1) = 30d0
      ftil(0,1) = sqrt(ph_Hmass**2 + ftil(1,1)**2 + ftil(2,1)**2
     $     + ftil(3,1)**2)
      do mu = 0, 3
         kb(mu) = ftil(mu,1) + ftil(mu,2) + ftil(mu,3) + ftil(mu,4)
      enddo
      a = 0; b = 0
      a(0) = (kb(0) + kb(3))/2; a(3) = a(0)
      b(0) = (kb(0) - kb(3))/2; b(3) = -b(0)
      x = 0.91d0
      p1 = a/x
      p2 = b
      s12 = 2*dot(p1, p2)
      pgq = 4d0/3d0*(1 + (1 - x)**2)/x
      write(*,'(a,f8.4,a,3f10.3)') ' x =', x, '  Higgs mass, K~^2 ^1/2:',
     $     ph_Hmass, sqrt(dot(kb, kb))

      do iq = 2, 1, -1
c        Born g(a) c(b) -> H Q c Qbar, legs 4 = Q, 5 = c, 6 = Qbar
         bflav = (/ 0, 4, HZZ, iq, 4, -iq /)
         pb(:,1) = a; pb(:,2) = b; pb(:,3) = ftil(:,1)
         pb(:,4) = ftil(:,2); pb(:,5) = ftil(:,3); pb(:,6) = ftil(:,4)
         call setborn(pb, bflav, born, bornjk, bmunu, bornsub)
c        real s c -> H s c Q Qbar
         rflav = (/ 3, 4, HZZ, 3, 4, iq, -iq /)
         do ikt = 1, 4
            kt = kts(ikt)
            avg = 0
            do iphi = 1, nphi
               phi = 2*pi*(iphi - 0.5d0)/nphi
c              k4 = alpha p1 + beta p2 + kT, alpha = 1 - x - beta
               beta = (1 - x)/2*(1 - sqrt(1 - 4*kt**2/(s12*(1 - x)**2)))
               alpha = 1 - x - beta
               k4 = alpha*p1 + beta*p2
               k4(1) = kt*cos(phi); k4(2) = kt*sin(phi)
               kk = p1 + p2 - k4
               pr(:,1) = p1; pr(:,2) = p2; pr(:,4) = k4
               call lmap(kb, kk, ftil(0,1), pr(0,3))
               call lmap(kb, kk, ftil(0,3), pr(0,5))
               call lmap(kb, kk, ftil(0,2), pr(0,6))
               call lmap(kb, kk, ftil(0,4), pr(0,7))
               call setreal(pr, rflav, amp2)
               avg = avg + amp2/nphi
            enddo
            c = avg*2*dot(p1, k4)*x/(16*pi**2*pgq*born)
            write(*,'(a,i2,a,es9.1,a,es14.6,a,es14.6,a,f10.6)')
     $           ' Q =', iq, '  kT =', kt, '  <R> =', avg, '  B_g =',
     $           born, '  c =', c
         enddo
      enddo
      end

      subroutine mom(pt, y, phi, p)
      implicit none
      real * 8 pt, y, phi, p(0:3)
      p(0) = pt*cosh(y); p(1) = pt*cos(phi); p(2) = pt*sin(phi)
      p(3) = pt*sinh(y)
      end

      real * 8 function dot(p, q)
      implicit none
      real * 8 p(0:3), q(0:3)
      dot = p(0)*q(0) - p(1)*q(1) - p(2)*q(2) - p(3)*q(3)
      end

c q = Lambda p, Lambda: K~ -> K (K^2 = K~^2), the inverse of the CS map
c Lambda^mu_nu = g - 2 (K+K~)^mu (K+K~)_nu/(K+K~)^2 + 2 K^mu K~_nu/K~^2
      subroutine lmap(kt, k, p, q)
      implicit none
      real * 8 kt(0:3), k(0:3), p(0:3), q(0:3), s(0:3), dot, a1, a2
      external dot
      s = k + kt
      a1 = 2*dot(s, p)/dot(s, s)
      a2 = 2*dot(kt, p)/dot(kt, kt)
      q = p - a1*s + a2*k
      end
