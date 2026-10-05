c Helpers for the problem-2 tests. A physical NC four-quark point is
c given by role: q(1) = S_in, q(2) = C_in, 3 = H, 4 = S_out, 5 = C_out,
c 6 = U (quark), 7 = Ub (antiquark), S (C) the line from beam 1 (2),
c U Ub the q qbar pair. Entries with these flavours:
c   ityp 1: outgoing pair  (i, j, H, i, j, q, -q),   tags 1 2 0 1 2 5 5
c   ityp 2: tag on beam 1  (i, j, H, q, j, -q, i),   tags 5 2 0 1 2 1 5
c   ityp 3: tag on beam 2  (i, j, H, i, q, -q, j),   tags 1 5 0 1 2 2 5
c findent returns the real entry (0 if absent); evalent evaluates the
c entry with setreal along POWHEG's path (alr_tag = its first alr, the
c momenta mapped to the alr's leg order by (flavour, tag)).
      integer function findent(ityp,i,j,q,HZZ)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      integer ityp,i,j,q,HZZ,f(7),t(7),k,l
      if (ityp.eq.1) then
         f = (/ i, j, HZZ, i, j, q, -q /); t = (/ 1,2,0,1,2,5,5 /)
      elseif (ityp.eq.2) then
         f = (/ i, j, HZZ, q, j, -q, i /); t = (/ 5,2,0,1,2,1,5 /)
      else
         f = (/ i, j, HZZ, i, q, -q, j /); t = (/ 1,5,0,1,2,2,5 /)
      endif
      findent = 0
      do k=1,flst_nreal
         do l=1,7
            if (flst_real(l,k).ne.f(l)) goto 10
            if (l.ne.3.and.flst_realtags(l,k).ne.t(l)) goto 10
         enddo
         findent = k
         return
 10      continue
      enddo
      end

c role of each real leg, per entry type
      subroutine roles(ityp,r)
      implicit none
      integer ityp,r(7)
      if (ityp.eq.1) then
         r = (/ 1,2,3,4,5,6,7 /)
      elseif (ityp.eq.2) then
         r = (/ 1,2,3,6,5,7,4 /)
      else
         r = (/ 1,2,3,4,6,7,5 /)
      endif
      end

      subroutine evalent(ient,ityp,prole,amp2,ialr)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_flst_2.h'
      integer ient,ityp,ialr,r(7),alr,l,m,used(7)
      real*8 prole(0:3,7),amp2,pa(0:3,7)
      logical same
      amp2 = 0
      ialr = 0
      if (ient.eq.0) return
      call roles(ityp,r)
      do alr=1,flst_nalr
         if (same(flst_real(1,ient),flst_realtags(1,ient),
     $        flst_alr(1,alr),flst_alrtags(1,alr))) goto 20
      enddo
      stop 'evalent: no alr'
 20   ialr = alr
      used = 0
      do l=1,7
         do m=1,7
            if (used(m).eq.0.and.flst_alr(l,alr).eq.flst_real(m,ient)
     $           .and.flst_alrtags(l,alr).eq.flst_realtags(m,ient)) then
               used(m) = 1
               pa(:,l) = prole(:,r(m))
               goto 30
            endif
         enddo
         stop 'evalent: leg map'
 30      continue
      enddo
      alr_tag = alr
      realequiv_tag = 0
      call setreal(pa,flst_alr(1,alr),amp2)
      alr_tag = 0
      end

      logical function same(f1,t1,f2,t2)
      implicit none
      integer f1(7),t1(7),f2(7),t2(7),i,k,used(7)
      same = .false.
      if (f1(1).ne.f2(1) .or. t1(1).ne.t2(1)) return
      if (f1(2).ne.f2(2) .or. t1(2).ne.t2(2)) return
      used = 0
      do i=3,7
         do k=3,7
            if (used(k).eq.0 .and. f1(i).eq.f2(k)
     $           .and. t1(i).eq.t2(k)) then
               used(k) = 1
               goto 10
            endif
         enddo
         return
 10      continue
      enddo
      same = .true.
      end

      real * 8 function dot(p, q)
      implicit none
      real * 8 p(0:3), q(0:3)
      dot = p(0)*q(0) - p(1)*q(1) - p(2)*q(2) - p(3)*q(3)
      end

      subroutine mom(pt, y, phi, p)
      implicit none
      real * 8 pt, y, phi, p(0:3)
      p(0) = pt*cosh(y); p(1) = pt*cos(phi); p(2) = pt*sin(phi)
      p(3) = pt*sinh(y)
      end

c q = Lambda p, Lambda: K~ -> K (K^2 = K~^2), the inverse of the CS map
      subroutine lmap(kt, k, p, q)
      implicit none
      real * 8 kt(0:3), k(0:3), p(0:3), q(0:3), s(0:3), dot, a1, a2
      external dot
      s = k + kt
      a1 = 2*dot(s, p)/dot(s, s)
      a2 = 2*dot(kt, p)/dot(kt, kt)
      q = p - a1*s + a2*k
      end
