c List proVBFH's real flavour structures with their tags, as built by its
c init_processes (run in a directory with powheg.input), and the line
c assignment the tags encode (11, 22: both extra partons on line 1/2;
c 12, 21: one on each line). For developing stage 2 of proVBFH-cs.
      program list_reals
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      integer j,k,n(0:4),ig,tf
      integer tag_factor
      common/cdoubletag/tag_factor
      call init_processes
      call finalize_tags
      write(*,*) 'nborn, nreal, tag_factor: ',flst_nborn,flst_nreal,
     $     tag_factor
      n = 0
      do j=1,flst_nreal
         call classify(flst_realtags(1,j),ig)
         n(ig) = n(ig)+1
         write(*,'(i5,a,7i4,a,7i6,a,i2)')
     $        j,' flav',(flst_real(k,j),k=1,nlegreal),
     $        ' tags',(flst_realtags(k,j),k=1,nlegreal),'  class',ig
      enddo
      write(*,*) 'entries per class (0 = other, 1..4 = 11,22,12,21): ',n
      do j=1,flst_nborn
         write(*,'(i5,a,6i4,a,6i6)') j,' born',
     $        (flst_born(k,j),k=1,nlegborn),
     $        ' tags',(flst_borntags(k,j),k=1,nlegborn)
      enddo
      end

      subroutine classify(tags,ig)
c     as proVBFH's pwhg_analysis.f
      implicit none
      include 'nlegborn.h'
      include 'tags.h'
      integer tags(nlegreal),ig,s,tag_factor
      common/cdoubletag/tag_factor
      s = sum(tags)
      ig = 0
      if (s.eq.8+2*iqpairtag*tag_factor .or. s.eq.8) ig = 1
      if (s.eq.10+2*iqpairtag*tag_factor .or. s.eq.10) ig = 2
      if (s.eq.9 .and. tags(6).eq.1 .and. tags(7).eq.2) ig = 3
      if (s.eq.9 .and. tags(6).eq.2 .and. tags(7).eq.1) ig = 4
      end
