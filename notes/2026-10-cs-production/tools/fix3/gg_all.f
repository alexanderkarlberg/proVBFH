c Count the gg-initiated real entries of VBF_HJJJ before and after
c genflavreglist, and the FKS regions (alr) with a gg-initiated real.
      program gg_entries
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      integer j,ngg,nggalr
      call init_processes
      ngg = 0
      do j=1,flst_nreal
         if (flst_real(1,j).eq.0 .and. flst_real(2,j).eq.0) then
            ngg = ngg + 1
            write(*,'(a,7i4,a,l2)') ' gg real ',
     $           flst_real(:,j),' fix3 branch:',
     $           flst_real(4,j).eq.-flst_real(6,j)
         endif
      enddo
      write(*,*) 'before genflavreglist: nreal, gg reals: ',
     $     flst_nreal, ngg
      call setup_reson
      call genflavreglist
      nggalr = 0
      do j=1,flst_nalr
         if (flst_alr(1,j).eq.0 .and. flst_alr(2,j).eq.0)
     $        nggalr = nggalr + 1
      enddo
      write(*,*) 'after: nreal, nalr, gg regions: ',
     $     flst_nreal, flst_nalr, nggalr
      end
