      program pub
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_st.h'
      include 'pwhg_flst_2.h'
      real*8 p(0:3,nlegborn), born, bornjk(nlegborn,nlegborn),
     $     bmunu(0:3,0:3,nlegborn), bornsub(2)
      integer j, k, n
      call init_flsttag
      call init_phys
      st_alpha = 0.118d0
      st_muren2 = 100d0**2
      call testpoint(p)
      write(*,*) 'nborn', flst_nborn
      n = 0
      do j=1,flst_nborn
         if (flst_born(6,j).ne.0) cycle
         if (flst_born(1,j).le.0 .or. flst_born(2,j).le.0) cycle
         n = n + 1
         if (n.gt.12) exit
         call setborn(p, flst_born(1,j), born, bornjk, bmunu, bornsub)
         write(*,'(a,6i4,a,6i3,a,es16.8)') ' flav',(flst_born(k,j),k=1,6),
     $        ' tags',(flst_borntags(k,j),k=1,6),' born/as ',born/st_alpha
      enddo
      end
