
      subroutine incl2pwhg(www)
      use matrix_element
      implicit none
      include 'hepevt.h'
!      include 'nlegborn.h'
      include 'brinclude.h'
      include 'incl2pwhg.h'
      include 'pwhg_anexinf.h'
      include 'pwhg_kn.h'
      include 'phspcuts.h'
      
      double precision www

      double precision psave(5,br_nlegborn),psub(5,br_nlegborn,4)
      real * 8 p1(1:5),p2(1:5),p3(1:5),p4(1:5),ptmp(1:5)
      integer ig,j
      logical ini
      data ini/.true./
      save ini

      if(ini) then
         call getQ2min(0,Qmin)
         Qmin = sqrt(Qmin)
         if(Qmin.lt.1d-1) then
            print*, 'WARNING :Qmin in PDF is', Qmin
            Qmin = 1d0          ! The PDF should return a Qmin value
                                ! which is around 1 GeV
            print*, 'Setting Qmin =', Qmin
         endif
         ini = .false.
      endif

      incl_tot = 0d0
      psave = 0d0

      if(.not.btildennloon) return 
      if(dummy_analysis) return
      
      psave(1:3,1:br_nlegborn)=brkn_pborn(1:3,1:br_nlegborn)
      psave(4,1:br_nlegborn)=brkn_pborn(0,1:br_nlegborn)

      do j=1,br_nlegborn
         psave(5,j) = sqrt(abs(psave(4,j)**2 - psave(1,j)**2 - psave(2
     $        ,j)**2 - psave(3,j)**2))
      enddo
      
      call phep2struct_funct(psave(:,1:5))

      if((brkn_jacborn.ne.0d0).and.(min(Q1_sq,Q2_sq).gt.Qmin**2)) then
         incl_tot= eval_matrix_element(1,order, x1, x2, kn_beams(:,1),
     $        kn_beams(:,2), vq1, vq2, ptH)
         incl_tot = incl_tot * Q1_sq * Q2_sq / (x1 * x2) 

         incl_tot = incl_tot * www * brkn_jacborn
      else
         incl_tot = 0d0
      endif
      
!      print*, 'incl_tot', incl_tot, www, brkn_jacborn
         
      end subroutine


      subroutine phep2struct_funct(p)
      implicit none
      include 'incl2pwhg.h'
      include 'nlegborn.h'
      include 'pwhg_kn.h'
      double precision invmass
      external invmass
      
      double precision p(1:5,1:nlegborn-1)

      x1 = 2d0*p(4,1)/(sqrt(kn_sbeams))
      x2 = 2d0*p(4,2)/(sqrt(kn_sbeams))
      
      ptH = sqrt(p(1,3)**2 + p(2,3)**2) 

      vq1(1:3) = p(1:3,4) - p(1:3,1) 
      vq2(1:3) = p(1:3,5) - p(1:3,2) 
      vq1(0) = p(4,4) - p(4,1) 
      vq2(0) = p(4,5) - p(4,2) 

      Q1_sq = invmass(vq1)**2
      Q2_sq = invmass(vq2)**2
      end

      function invmass(p)
      implicit none
      real * 8 p(0:3),invmass
      invmass = sqrt(abs(p(0)**2-p(1)**2-p(2)**2-p(3)**2))
      end function
