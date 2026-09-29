c List the FKS singular regions (alr) that POWHEG's genflavreglist builds
c from the flavour list of the PUBLIC POWHEG-BOX-V2 VBF_HJJJ (svn r4135),
c with its tags, for the real entries of the NC type q(1) -> q(1) + q'qbar'
c (pair tag 5 in that code): does any region have the outgoing quark of the
c incoming flavour collinear to the beam? Uses the tags of this code
c (lines tagged 1, 2; the q'qbar' pair of the four-quark reals tagged 5).
c Initialisation as in the public init_phys.f: init_processes,
c setup_reson, genflavreglist (no finalize_tags).
      program list_regions
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_flg.h'
      integer j,k,alr,nfound,nmatch,m1(7),m2(7),t1(7),t2(7)
      logical match, same
      call init_processes
      call setup_reson
      call genflavreglist
      write(*,*) 'nborn, nreal, nalr: ',flst_nborn,flst_nreal,flst_nalr
      nfound = 0
      do j=1,flst_nreal
c        NC entries of type q(1) -> q(1) + pair(5) on line 1: legs
c        1 and 4 the same flavour with tag 1, legs 6, 7 the pair (tag 5)
         if (flst_realtags(1,j).ne.1 .or. flst_realtags(4,j).ne.1) cycle
         if (flst_realtags(6,j).ne.5 .or. flst_realtags(7,j).ne.5)
     $        cycle
         if (flst_real(1,j).ne.flst_real(4,j)) cycle
         nfound = nfound + 1
         nmatch = 0
         do alr=1,flst_nalr
c           the same flavour structure (legs possibly reordered)?
            if (same(flst_real(1,j),flst_realtags(1,j),
     $           flst_alr(1,alr),flst_alrtags(1,alr))) then
               nmatch = nmatch + 1
               if (nfound.le.12) write(*,'(a,i5,a,7i4,a,7i4,a,i3,a,i3)')
     $              ' real',j,' alr flav',(flst_alr(k,alr),k=1,7),
     $              ' tags',(flst_alrtags(k,alr),k=1,7),
     $              '  emitter',flst_emitter(alr),
     $              '  radiated leg',nlegreal
            endif
         enddo
         if (nfound.le.12 .or. nmatch.eq.0) write(*,'(a,i5,a,7i4,a,i3)')
     $        ' real',j,' flav',(flst_real(k,j),k=1,7),
     $        '  number of regions:',nmatch
      enddo
      write(*,*) 'NC q -> q + pair(5) entries on line 1: ',nfound
c     summary over all of them: number of regions and emitters
      call summary(1)
c     and the CC entries with the pair tag on the incoming quark,
c     which need the initial-state region
      call summary(2)
      end

c isel = 1: NC q(1) -> q(1) + pair(5); isel = 2: the pair tag 5 on the
c incoming quark of line 1 (only CC entries exist). Count the entries by the
c set of their emitters (0: both initial, 1/2: one initial, >2: final)
      subroutine summary(isel)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      integer isel,j,alr,nisr,nfsr,nent,hist(0:3,0:3),a,b
      logical same
      hist = 0; nent = 0
      do j=1,flst_nreal
         if (isel.eq.1) then
            if (flst_realtags(1,j).ne.1 .or. flst_realtags(4,j).ne.1)
     $           cycle
            if (flst_realtags(6,j).ne.5.or.flst_realtags(7,j).ne.5)
     $           cycle
            if (flst_real(1,j).ne.flst_real(4,j)) cycle
         else
            if (flst_realtags(1,j).ne.5) cycle
         endif
         nent = nent + 1
         nisr = 0; nfsr = 0
         do alr=1,flst_nalr
            if (same(flst_real(1,j),flst_realtags(1,j),
     $           flst_alr(1,alr),flst_alrtags(1,alr))) then
               if (flst_emitter(alr).le.2) then
                  nisr = nisr + 1
               else
                  nfsr = nfsr + 1
               endif
            endif
         enddo
         a = min(nisr,3); b = min(nfsr,3)
         hist(a,b) = hist(a,b) + 1
      enddo
      write(*,'(a,i2,a,i5)') ' selection',isel,': entries ',nent
      do a=0,3
         do b=0,3
            if (hist(a,b).gt.0) write(*,'(a,i2,a,i2,a,i5)')
     $           '   initial-state regions',a,', final-state regions',b,
     $           ': entries',hist(a,b)
         enddo
      enddo
      end

c the same flavours and tags on the incoming legs 1, 2 and the same
c multiset of (flavour, tag) on the outgoing legs 3..7
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

