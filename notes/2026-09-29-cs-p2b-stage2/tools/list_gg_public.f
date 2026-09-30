c print the gg-initiated real entries of the public VBF_HJJJ flavour list
      program list_gg_public
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      integer j,k,n
      call init_processes
      n = 0
      do j=1,flst_nreal
         if (flst_real(1,j).eq.0 .and. flst_real(2,j).eq.0) then
            n = n + 1
            write(*,'(i5,7i4,a,7i3)') j,(flst_real(k,j),k=1,7),
     $           '  tags',(flst_realtags(k,j),k=1,7)
         endif
      enddo
      write(*,*) 'gg entries:',n
      end
