c Problem 2 (bug report): initial-state collinear limits of the NC
c four-quark real s c -> H s c Q Qbar, per tagged real entry, along
c POWHEG's evaluation path (alr_tag set). Kinematics as in
c is_limit_public.f (exact Catani-Seymour initial-initial map, x = 0.91,
c 16 azimuths).
c  ilim = 1: outgoing s || beam 1, Born g(x p1) c -> H Q c Qbar
c  ilim = 2: outgoing c || beam 2, Born s g(x p2) -> H s Q Qbar
c For each entry type (1: pair outgoing, 2: pair tag on incoming 1,
c 3: on incoming 2) c_e = <R_e>_phi x 2 p_beam.k / (16 pi^2 P_gq(x) B_g).
c Expected after the fix: the entry with the initial-state region of
c this limit (2 for ilim = 1, 3 for ilim = 2) -> 1, the others -> 0
c (R_e finite, c_e ~ kT^2), and sum = the unfixed matrix element.
      program fix2_limit
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      include 'pwhg_flst_2.h'
      include 'PhysPars.h'
      integer HWW, HZZ, iq, iphi, nphi, ikt, mu, ilim, it, ient(3),
     $     findent, ialr(3), nkt
      parameter (nphi = 16, nkt = 6)
      real * 8 pb(0:3,nlegborn), pr(0:3,nlegreal), born, bornjk(nlegborn,
     $     nlegborn), bmunu(0:3,0:3,nlegborn), bornsub(2), amp2, avg(3)
      real * 8 x, kt, phi, p1(0:3), p2(0:3), a(0:3), b(0:3), k4(0:3),
     $     kk(0:3), kb(0:3), ftil(0:3,4), s12, beta, alpha, pgq, c(3),
     $     kts(nkt), dot, pbeam(0:3), csum, rmax(3), rmin(3)
      integer bflav(nlegborn)
      data kts /1d0, 1d-1, 1d-2, 1d-3, 1d-4, 1d-5/
      external dot, findent

      call init_flsttag
      call init_phys
      call particle_identif(HWW, HZZ)
      st_alpha = 0.118d0
      st_muren2 = 100d0**2
      st_mufact2 = 100d0**2
      alr_tag = 0
      realequiv_tag = 0
      write(*,*) 'nreal = ', flst_nreal

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
      pgq = 4d0/3d0*(1 + (1 - x)**2)/x

      do ilim = 1, 2
      do iq = 2, 1, -1
         do it = 1, 3
            ient(it) = findent(it, 3, 4, iq, HZZ)
         enddo
         if (ilim.eq.1) then
            p1 = a/x
            p2 = b
            pbeam = p1
c           Born g(a) c(b) -> H Q c Qbar
            bflav = (/ 0, 4, HZZ, iq, 4, -iq /)
            pb(:,4) = ftil(:,2); pb(:,5) = ftil(:,3); pb(:,6) = ftil(:,4)
         else
            p1 = a
            p2 = b/x
            pbeam = p2
c           Born s(a) g(b) -> H s Q Qbar
            bflav = (/ 3, 0, HZZ, 3, iq, -iq /)
            pb(:,4) = ftil(:,3); pb(:,5) = ftil(:,2); pb(:,6) = ftil(:,4)
         endif
         pb(:,1) = a; pb(:,2) = b; pb(:,3) = ftil(:,1)
         call setborn(pb, bflav, born, bornjk, bmunu, bornsub)
         s12 = 2*dot(p1, p2)
         write(*,'(/,a,i2,a,i2,a,3i6)') ' limit', ilim, '  Q =', iq,
     $        '  real entries (out, in1, in2):', ient
         write(*,'(a,es14.6)') ' Born B_g =', born
         write(*,'(a)') '   kT [GeV]    c_out        c_in1        c_in2'
     $        //'        c_sum        <R_out>      <R_in1>      <R_in2>'
     $        //'      max/min_phi R_in(lim)'
         do ikt = 1, nkt
            kt = kts(ikt)
            avg = 0
            rmax = -1d300; rmin = 1d300
            do iphi = 1, nphi
               phi = 2*pi*(iphi - 0.5d0)/nphi
               beta = (1 - x)/2*(1 - sqrt(1 - 4*kt**2/(s12*(1 - x)**2)))
               alpha = 1 - x - beta
               if (ilim.eq.1) then
                  k4 = alpha*p1 + beta*p2
               else
                  k4 = alpha*p2 + beta*p1
               endif
               k4(1) = kt*cos(phi); k4(2) = kt*sin(phi)
               kk = p1 + p2 - k4
c              roles: S_in, C_in, H, S_out, C_out, U, Ub
               pr(:,1) = p1; pr(:,2) = p2
               call lmap(kb, kk, ftil(0,1), pr(0,3))
               call lmap(kb, kk, ftil(0,2), pr(0,6))
               call lmap(kb, kk, ftil(0,4), pr(0,7))
               if (ilim.eq.1) then
                  pr(:,4) = k4
                  call lmap(kb, kk, ftil(0,3), pr(0,5))
               else
                  pr(:,5) = k4
                  call lmap(kb, kk, ftil(0,3), pr(0,4))
               endif
               do it = 1, 3
                  call evalent(ient(it), it, pr, amp2, ialr(it))
                  avg(it) = avg(it) + amp2/nphi
                  rmax(it) = max(rmax(it), amp2)
                  rmin(it) = min(rmin(it), amp2)
               enddo
            enddo
            c = avg*2*dot(pbeam, k4)*x/(16*pi**2*pgq*born)
            csum = c(1) + c(2) + c(3)
            write(*,'(es9.1,4f13.8,3es13.5,f9.4)') kt, c, csum, avg,
     $           rmax(ilim+1)/max(rmin(ilim+1),1d-300)
         enddo
      enddo
      enddo
      end
