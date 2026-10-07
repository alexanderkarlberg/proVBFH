c For the gg-initiated FKS regions (flst_alr, the flavour order passed to
c setreal), count those that reach the NC branch of compreal_hjjj changed
c by fix 3 (leg4 < 0 < leg5, leg4 = -leg6) with pairs of different type.
      program gg_alr
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      integer j,n,nbr,ndiff
      call init_processes
      call setup_reson
      call genflavreglist
      n = 0
      nbr = 0
      ndiff = 0
      do j=1,flst_nalr
         if (flst_alr(1,j).ne.0 .or. flst_alr(2,j).ne.0) cycle
         n = n + 1
         if (n.le.6) write(*,'(a,7i4)') ' gg alr ',flst_alr(:,j)
         if (flst_alr(4,j).lt.0 .and. flst_alr(5,j).gt.0 .and.
     $        flst_alr(4,j).eq.-flst_alr(6,j)) then
            nbr = nbr + 1
            if (mod(abs(flst_alr(4,j)),2).ne.mod(abs(flst_alr(5,j)),2))
     $           ndiff = ndiff + 1
         endif
      enddo
      write(*,*) 'gg regions, in fix-3 branch, with types differing: ',
     $     n, nbr, ndiff
      end
