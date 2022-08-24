      subroutine analysis(dsig0)
      implicit none
      real * 8 dsig(7),dsig0
      include 'hepevt.h'
      include 'pwhg_math.h'  
      include 'brinclude.h'
!     include 'nlegborn.h'
      include 'pwhg_weights.h'
!      include 'pwhg_flst.h'
      include 'pwhg_flg.h'
      include 'phspcuts.h'
      include 'incl2pwhg.h'
c     
C     =====================================================================
C     extra info needed for particular analysis 
      include 'pwhg_anexinf.h' 
      include 'tags.h' 
      integer tag_factor
      common/cdoubletag/tag_factor
      character * 6 WHCPRG
      common/cWHCPRG/WHCPRG
      data WHCPRG/'NLO   '/
      integer ig,igmin,igmax,j,iuborn,alr  
      double precision tot(0:4), psave(5,nlegreal),psub(5,nlegreal,4)
      real * 8 p1(1:5),p2(1:5),p3(1:5),p4(1:5),ptmp(1:5)

      integer entered_analysis, passed_analysis
      data entered_analysis,passed_analysis/0,0/
      save entered_analysis, passed_analysis


c      print*, 'Entered with ', anexinf_type, 'current xsec ', vbftot
      tot = 0d0 
      vbftot_save = vbftot ! Save for later
c      print*, anexinf_sigarr
c      stop
      if (anexinf_type == 'born' .or. anexinf_type == 'virt' .or. 
     C     anexinf_type == 'colr' ) then 

C     sanity check 
         if (anexinf_narr /= flst_nborn ) then 
            write(*,*) 'anexinf_narr, flst_nborn', anexinf_narr, flst_nborn
            write(*,*) 'analysis: anexinf_narr /= flst_nborn' 
            call exit(-1) 
         endif

         ig=0; igmin = 0; igmax = 2 
         do j=1,flst_nborn
            ig=flst_borntags(nlegborn,j)
            if (ig == 0) then 
               write(*,*) 'running with no tags, 
     C but using analysis with subtraciton' 
               call exit(-1) 
            endif
            if (ig > igmax) then 
               write(*,*) 'analysis: ig out of range' 
            endif 
            tot(0)=tot(0)+anexinf_sigarr(j)
!            print*, 'flst_nborn, anexinf_sigarr', j, anexinf_sigarr(j)
            if(ig.ne.0) then
               tot(ig)=tot(ig)+anexinf_sigarr(j)
            endif
         enddo
      
      elseif (anexinf_type .eq. 'realct') then 
         ig = 0; igmin = 0; igmax = 2 
         tot = 0d0 
         do alr=1,flst_nalr
            iuborn = flst_alr2born(alr) 
            ig=flst_borntags(nlegborn,iuborn)
            tot(0)=tot(0)+anexinf_sigarr(alr)
            if (ig > igmax .or. ig <1 ) then 
               write(*,*) 'analysis: ig out of range' 
            endif 
            if(ig.ne.0) then
               tot(ig)=tot(ig)+anexinf_sigarr(alr)
            endif
         enddo

            
      elseif (anexinf_type .eq. 'real') then 
         ig = 0; igmin = 0; igmax = 4 
         do j=1,flst_nalr
c     use tag to split into 11, 22, 12 or 21
            if (any(flst_alrtags(1:2,j) == 0) .or. any(flst_alrtags(5:nlegreal,j) == 0)) then
c     do nothing, since we are running without tags, only set igmax
               igmax = 0 
            elseif (any(flst_alrtags(1:2,j) == iqpairtag*tag_factor) .or. 
     $              any(flst_alrtags(5:nlegreal,j) == iqpairtag*tag_factor)) then
c     do nothing, since we are running without tags and it's the qqqqqq channel
               igmax = 0 
            else
               if (sum(flst_alrtags(1:nlegreal,j)) == 8+2*iqpairtag*tag_factor) then
                  ig=1          ! 11 case and qqqqqq
               elseif (sum(flst_alrtags(1:nlegreal,j)) == 10+2*iqpairtag*tag_factor) then
                  ig=2          ! 22 case and qqqqqq
               elseif (sum(flst_alrtags(1:nlegreal,j)) == 8) then
                  ig=1          ! 11 case
               elseif (sum(flst_alrtags(1:nlegreal,j)) == 10) then
                  ig=2          ! 22 case
               elseif (sum(flst_alrtags(1:nlegreal,j)) == 9 .and. 
     $                 flst_alrtags(nlegborn,j) == 1 .and. flst_alrtags(nlegreal,j) == 2) then ! case 1 2 -> 1 2 1 2 
                  ig=3          ! 12 case
               elseif (sum(flst_alrtags(1:nlegreal,j)) == 9 .and. 
     $                 flst_alrtags(nlegborn,j) == 2 .and. flst_alrtags(nlegreal,j) == 1) then ! case 1 2 -> 1 2 2 1 
                  ig=4          ! 21 case
               else
                  write(*,*) 'flst_alrtags', flst_alrtags 
                  stop 'pwhg_analysis_new: invalid tags'
                  call exit(-1) 
               endif
C     reconstruct the total 
               tot(0) = tot(0) + anexinf_sigarr(j) 
C     and the various contributions 
               if(ig.ne.0) then
                  tot(ig) = tot(ig) + anexinf_sigarr(j)
               endif
            endif
         enddo
         
CCC     -- SANITY CHECKS.. 
C         if (abs(tot(0)/dsig0-1d0) .gt. 1d-2) then 
C            write(*,*) 'tot, weight',anexinf_type, tot(0),dsig0,
C     X           sum(anexinf_sigarr(1:flst_nalr)),tot(0)/dsig0, flst_nalr
C            write(*,*) 'pwhg_analysis_new: total not
C     X reconstructed properly'
C         endif
C         if (abs(tot(0)-sum(tot(1:4))) .gt. 1d-6) then 
C            write(*,*) 'tot', tot 
C            write(*,*) 'pwhg_analysis_new: total not
C     X split properly'
C         endif

      elseif(anexinf_type == 'nnlo') then
         igmin = 1
         igmax = 1
         tot = 0d0
         tot(1) = - incl_tot ! Minus due to convention of subtracting below
      else
         write(*,*) 'anexinf_type',anexinf_type
         write(*,*) 'analysis: type not allowed' 
         call exit(-1) 
      endif

c     store the momenta in a temporary array
      psave(1:5,1:nlegreal)=phep(1:5,1:nlegreal)
      do j=1,4
         psub(1:5,1:nlegreal,j)= psave(1:5,1:nlegreal) 
      enddo

      if (igmax .eq. 1) then    ! NNLO case
         psub = 0d0
         psub(1:3,1:br_nlegborn,1)=brkn_pborn(1:3,1:br_nlegborn)
         psub(4,1:br_nlegborn,1)=brkn_pborn(0,1:br_nlegborn)

         do j=1,br_nlegborn
            psub(5,j,1) = sqrt(abs(psub(4,j,1)**2 - psub(1,j,1)**2 -
     $           psub(2,j,1)**2 - psub(3,j,1)**2))
         enddo

      endif
      
c     Mapping for one or two gluons on upper line
      if (igmax .ge. 2) then 
         call mapping_1_1(psave(:,1),psave(:,5),psave(:,nlegborn),psave(:,nlegreal)
     $        ,p1,p2)
         psub(:,1,1)=p1
         psub(:,5,1)=p2
         psub(:,nlegborn,1)=0d0
         psub(:,nlegreal,1)=0d0
         
c     Mapping for one or two  gluons on lower line
         call mapping_2_2(psave(:,2),psave(:,6),psave(:,nlegborn),psave(:,nlegreal)
     $        ,p1,p2)
         psub(:,2,2)=p1
         psub(:,6,2)=p2
         psub(:,nlegborn,2)=0d0
         psub(:,nlegreal,2)=0d0            
      endif

      if (igmax .eq.4) then 
c     Mapping for Born gluon on upper line and real gluon on lower
         call mapping_1_2(psave(:,1),psave(:,5),psave(:,nlegborn)
     $        ,psave(:,2),psave(:,6),psave(:,nlegreal)
     $        ,p1,p2,p3,p4)
         psub(:,1,3)=p1
         psub(:,5,3)=p2
         psub(:,2,3)=p3
         psub(:,6,3)=p4
         psub(:,nlegborn,3)=0d0
         psub(:,nlegreal,3)=0d0

c     Mapping for Born gluon on lower line and real gluon on upper
         call mapping_1_2(psave(:,1),psave(:,5),psave(:,nlegreal)
     $        ,psave(:,2),psave(:,6),psave(:,nlegborn)
     $        ,p1,p2,p3,p4)
         psub(:,1,4)=p1
         psub(:,5,4)=p2
         psub(:,2,4)=p3
         psub(:,6,4)=p4
         psub(:,nlegborn,4)=0d0
         psub(:,nlegreal,4)=0d0
      endif

!      print*, tot(1:2), incl_tot(1:2)

c      print*, anexinf_type, tot(0:igmax), dummy_analysis
C     call now analysis several times to include subtraction terms ig =
C     0 denotes positive weight, ig > 0 denote subtraction terms
      passed_cuts_ig = .false. ! Initialise
      do ig = igmin,igmax 
         if (ig ==0) then 
            WHCPRG='NLO   '
            dsig(1) = tot(0) 
            phep(:5,:nlegreal) = psave(:5,:nlegreal) 
         elseif (ig > 0) then 
c     set WHCPRG to say that we are doing projection to Born (later) 
            WHCPRG='PROJEC'
            dsig(1) = - tot(ig) ! change sign since this is a subtraction term
            phep(:5,:nlegreal) = psub(:5,:nlegreal,ig) 
         endif
         if(phspcuts.and.(dummy_analysis.or..not.flg_nlotest)) then
c     We only enter here if we have turned on phase space cuts. If
c     flg_nlotest is false, it means that we are computing grids. Then
c     we should never enter below, as we don't want to fill
c     histograms. If dummy_analysis is true, it means we have yet to
c     compute the matrix element, but we want to check whether or not
c     this phase space point will pass the cuts.
            call phspcuts_analysis(dsig0)
         else
            call user_analysis(dsig(1))
         endif

         passed_cuts_ig(ig) = passed_cuts
         if(passed_cuts.and..not.dummy_analysis) then
            vbftot = vbftot + dsig(1) ! This is the xsec passing the VBF
                                      ! cuts in the analysis
         endif 
      enddo

      if(all(passed_cuts_ig(0:igmax)).and..not.dummy_analysis) then 
! returns true if all elements up to 2 (born) or 4 (real) are true. In
! this case the sum of the parts have to be zero
         vbftot = vbftot_save   ! Numerically more stable as dsig can
                                ! become huge
      endif
!      print*, 'passed_cuts_ig', passed_cuts_ig(0:igmax), dsig(1), dsig0
!      print*, 'tot', tot
!      print*, 'vbftot', vbftot
      end
