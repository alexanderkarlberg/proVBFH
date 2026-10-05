c Problem 2: the total real matrix element of every flavour structure,
c summed over the real entries that carry it, at random real points.
c Writes (unit 22) one line per point and per entry 1..nold (the
c entries of the unfixed list): for every entry setreal along POWHEG's
c path (alr_tag set) with the momenta assigned to its legs in order
c [roles of the NC outgoing-pair entries: S_in, C_in, H, S_out, C_out,
c U, Ub]; for the NC outgoing-pair entries (i,j,H,i,j,q,-q) also the
c new entries of the same flavour structure (pair tag on incoming 1 or
c 2), evaluated at the same physical point, and the sum.
      program fix2_sum
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      include 'pwhg_flst_2.h'
      include 'PhysPars.h'
      integer nold, np
      parameter (nold = 1277, np = 5)
      integer HWW, HZZ, ip, j, k, mu, ient2, ient3, findent, ialr, nnc
      real * 8 pr(0:3,7), pt, y, phi, ptot(0:3), r1, r2, r3, rnd
      real * 8 seed
      common/cseed/seed
      external findent, rnd
      call init_flsttag
      call init_phys
      call particle_identif(HWW, HZZ)
      st_alpha = 0.118d0
      st_muren2 = 100d0**2
      st_mufact2 = 100d0**2
      alr_tag = 0
      realequiv_tag = 0
      seed = 12345d0
      write(*,*) 'nreal = ', flst_nreal
      do ip = 1, np
c        four massless partons (legs 4..7), H (leg 3), beams from
c        momentum conservation
         ptot = 0
         do k = 4, 7
            pt = 20 + 100*rnd()
            y = -3 + 6*rnd()
            phi = 2*pi*rnd()
            call mom(pt, y, phi, pr(0,k))
            ptot = ptot + pr(:,k)
         enddo
         pr(1,3) = -ptot(1); pr(2,3) = -ptot(2)
         pr(3,3) = 200*(rnd() - 0.5d0)
         pr(0,3) = sqrt(ph_Hmass**2 + pr(1,3)**2 + pr(2,3)**2
     $        + pr(3,3)**2)
         ptot = ptot + pr(:,3)
         pr(:,1) = 0; pr(:,2) = 0
         pr(0,1) = (ptot(0) + ptot(3))/2; pr(3,1) = pr(0,1)
         pr(0,2) = (ptot(0) - ptot(3))/2; pr(3,2) = -pr(0,2)
         nnc = 0
         do j = 1, nold
            call evalent(j, 1, pr, r1, ialr)
            r2 = 0; r3 = 0; ient2 = 0; ient3 = 0
            if (flst_real(3,j).eq.HZZ .and. flst_realtags(6,j).eq.5
     $           .and. flst_realtags(7,j).eq.5) then
               nnc = nnc + 1
               ient2 = findent(2, flst_real(1,j), flst_real(2,j),
     $              flst_real(6,j), HZZ)
               ient3 = findent(3, flst_real(1,j), flst_real(2,j),
     $              flst_real(6,j), HZZ)
               call evalent(ient2, 2, pr, r2, ialr)
               call evalent(ient3, 3, pr, r3, ialr)
            endif
            write(22,'(i2,i5,7i4,2i5,4es24.15)') ip, j,
     $           (flst_real(k,j),k=1,7), ient2, ient3, r1, r2, r3,
     $           r1 + r2 + r3
         enddo
      enddo
      write(*,*) 'NC outgoing-pair entries per point: ', nnc
      end

      real * 8 function rnd()
      implicit none
      real * 8 seed
      common/cseed/seed
      seed = mod(seed*16807d0, 2147483647d0)
      rnd = seed/2147483647d0
      end
