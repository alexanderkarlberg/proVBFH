c Adapted from list_regions_public.f (bug report, problem 2): FKS regions
c that POWHEG's genflavreglist builds from init_processes' flavour list
c and tags. For every real entry: number of initial- and final-state
c regions (dumped to fort.21 for a diff between builds), and a summary
c per class of four-quark entries:
c   NC/CC, pair tag 5 on the outgoing pair (legs 6,7), on incoming 1,
c   on incoming 2. Initialisation as in init_phys: init_processes,
c   setup_reson, genflavreglist.
      program list_regions_fix2
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      integer HWW,HZZ,j,k,alr,ic,nisr,nfsr,ncls,hist(0:3,0:3,8),a,b,
     $     nent(8),firstex(8)
      character*40 cname(8)
      logical same
      data cname /'NC pair out (6,7), q(1)->q(1) line 1',
     $     'NC pair out (6,7), all','NC pair tag on incoming 1',
     $     'NC pair tag on incoming 2','CC pair out (6,7)',
     $     'CC pair tag on incoming 1','CC pair tag on incoming 2',
     $     'other (gluon) reals'/
      call init_processes
      call setup_reson
      call genflavreglist
      call particle_identif(HWW,HZZ)
      write(*,*) 'nborn, nreal, nalr, nregular: ',flst_nborn,flst_nreal,
     $     flst_nalr,flst_nregular
      hist = 0; nent = 0; firstex = 0
      do j=1,flst_nreal
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
         write(21,'(i5,7i4,a,7i3,a,2i3)') j,(flst_real(k,j),k=1,7),
     $        ' |',(flst_realtags(k,j),k=1,7),' | isr fsr',nisr,nfsr
         ic = 8
         if (flst_real(1,j)*flst_real(2,j)*flst_real(4,j)*flst_real(5,j)
     $        *flst_real(6,j)*flst_real(7,j).ne.0) then
            if (flst_realtags(6,j).eq.5.and.flst_realtags(7,j).eq.5)
     $           ic = 2
            if (flst_realtags(1,j).eq.5) ic = 3
            if (flst_realtags(2,j).eq.5) ic = 4
            if (flst_real(3,j).eq.HWW) ic = ic + 3
         endif
         call addto(ic)
         if (ic.eq.2.and.flst_realtags(1,j).eq.1.and.
     $        flst_realtags(4,j).eq.1.and.
     $        flst_real(1,j).eq.flst_real(4,j)) call addto(1)
      enddo
      do ic=1,8
         write(*,'(a,a40,a,i5)') ' class ',cname(ic),': entries',
     $        nent(ic)
         do a=0,3
            do b=0,3
               if (hist(a,b,ic).gt.0) write(*,'(a,i2,a,i2,a,i5)')
     $              '     initial-state regions',a,
     $              ', final-state regions',b,': entries',hist(a,b,ic)
            enddo
         enddo
         j = firstex(ic)
         if (j.gt.0) then
            write(*,'(a,i5,7i4,a,7i3)') '     e.g. real',j,
     $           (flst_real(k,j),k=1,7),' tags',(flst_realtags(k,j),k=1,7)
            do alr=1,flst_nalr
               if (same(flst_real(1,j),flst_realtags(1,j),
     $              flst_alr(1,alr),flst_alrtags(1,alr)))
     $              write(*,'(a,i5,7i4,a,7i3,a,i2,a,i4)') '       alr',
     $              alr,(flst_alr(k,alr),k=1,7),' tags',
     $              (flst_alrtags(k,alr),k=1,7),' emitter',
     $              flst_emitter(alr),' uborn',flst_alr2born(alr)
            enddo
         endif
      enddo
      contains
      subroutine addto(icl)
      integer icl
      a = min(nisr,3); b = min(nfsr,3)
      hist(a,b,icl) = hist(a,b,icl) + 1
      nent(icl) = nent(icl) + 1
      if (firstex(icl).eq.0) firstex(icl) = j
      end subroutine
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
