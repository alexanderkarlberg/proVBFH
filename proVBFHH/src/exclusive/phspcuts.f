      subroutine phspcuts_analysis(dsig0)
      implicit none
      real * 8 dsig0

      integer maxjet
      parameter (maxjet=4)

      real *8 ptj(maxjet),yj(maxjet)
      real * 8 pj(0:3,4)
      integer njets
      double precision invmjj, mjj
      external mjj
      logical ini
      data ini/.true./
      save ini

c     COMMON block to cut on phasespace
      include 'pwhg_flg.h'
      include 'phspcuts.h'
      include 'nlegborn.h'
      include 'pwhg_kn.h'
c===============================================
      if (ini) then

         call setup_vbf_cuts

         write(*,*) '**************************************************'
         write(*,*) '**************************************************'
         write(*,*) '             USING PHASE SPACE CUTS               '
         write(*,*) '**************************************************'
         write(*,*) '**************************************************'
         ini = .false.
      endif

      if(kn_jacborn.eq.0d0) then ! Bad phase space point. No reason to continue.
         passed_cuts = .false.
         return
      endif
      
      if(.not.phspcuts) then
         print*, 'We are running without phase space cuts.'
         stop 'Should not have entered here.'
      endif

c     Build jets
      call buildjets(pj,njets,ptj,yj)
c     Check if jets satisfy VBF cuts
      call vbfcuts(pj,njets,ptj,yj,passed_cuts)
      
c     Ak We only cut on the pt of the jets and the invariant mass
c     for the phasespace. This should give better convergence than
c     cutting on everything.
      if(.not.flg_nlotest) then
         invmjj = mjj(pj(0,1),pj(0,2))
         passed_cuts = (min(pTj(1),pTj(2)).gt.ptjetmin) .and.
     $        (invmjj.gt.mjjmin) ! Only pt and mjj cuts for grids
      endif

      end


