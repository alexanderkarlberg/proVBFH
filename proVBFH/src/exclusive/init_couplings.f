      subroutine init_couplings
      implicit none
      include 'PhysPars.h'
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      include 'pwhg_physpar.h'
      include 'PhysPars_Higgs.h'
      include 'coupl.inc'
      real * 8 masswindow
      logical verbose
      parameter(verbose=.true.)
      integer i,j

      real*8 sthw,cthw,g2

      real *8 powheginput
      external powheginput 

      integer ckm_offdiag
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc   INDEPENDENT QUANTITIES       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c fermion masses:
      physpar_ml(1)=0d0 ! e
      physpar_ml(2)=0d0 ! mu
      physpar_ml(3)=0d0  ! tau
      physpar_mq(1)=0d0   ! up
      physpar_mq(2)=0d0   ! down
      physpar_mq(3)=0d0   ! strange
      physpar_mq(4)=0d0   ! charm
      physpar_mq(5)=0d0   ! bottom
c
c     number of light flavors
c     AK this number is not used anywhere. Should change flst_lightflav in init_processes.f
      st_nlight = 5
c
c diagonal or non-diag. CKM matrix at event generation level? 
c (note: fixed order MEs are always with diag. CKM)
c     
      ckm_offdiag = powheginput("#ckm_offdiag")
      if (ckm_offdiag.ne.1) ckm_offdiag = 0

c
      if (ckm_offdiag.eq.0) then 
c     DIAGONAL CKM 
         ph_CKM(1,1)=1d0 
         ph_CKM(1,2)=0d0 
         ph_CKM(1,3)=0d0
         ph_CKM(2,1)=0d0 
         ph_CKM(2,2)=1d0 
         ph_CKM(2,3)=0d0
         ph_CKM(3,1)=0d0
         ph_CKM(3,2)=0d0
         ph_CKM(3,3)=1d0
      else
         ph_CKM(1,1)=0.9748 	
         ph_CKM(1,2)=0.2225  	 
         ph_CKM(1,3)=0.0036  	
         ph_CKM(2,1)=0.2225  	
         ph_CKM(2,2)=0.9740 	
         ph_CKM(2,3)=0.041	
         ph_CKM(3,1)=0.009    
         ph_CKM(3,2)=0.0405   
         ph_CKM(3,3)=0.9992
      endif

c     initialize CKM with flavor indexes
      call inizialize_ph_CKM_matrix

      call coup_powheg_to_vbfnlo
c fermion masses:
c      physpar_ml(1)=0.511d-3 ! e
c      physpar_ml(2)=0.1057d0 ! mu
c      physpar_ml(3)=1.777d0  ! tau
c      physpar_mq(1)=0.33d0   ! up
c      physpar_mq(2)=0.33d0   ! down
c      physpar_mq(3)=0.50d0   ! strange
c      physpar_mq(4)=1.50d0   ! charm
c      physpar_mq(5)=4.80d0   ! bottom
      
c Higgs parameters: 
      ph_Hmass  = hmass
      ph_Hwidth = hwidth !4.02964352284941d-3 !hwidth
C
C gauge boson parameters:
      ph_Zmass   = zmass   
      ph_Zwidth  =  zwidth !2.4952d0 !zwidth
      ph_Wmass   = wmass     
      ph_Wwidth  = wwidth! 2.141d0 !wwidth

c EW parameters:
      ph_sthw2 = 1d0 - (ph_Wmass/ph_Zmass)**2
      sthw = SQRT(ph_sthw2)
      cthw = SQRT(1.d0 -ph_sthw2 )
      G2 = SQRT(8.d0*GFERMI/SQRT(2.d0))*ph_Zmass*cthw
      ph_alphaem = g2**2*ph_sthw2/(4.d0*PI)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc   DEPENDENT QUANTITIES       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      ph_sthw = sqrt(ph_sthw2)
      ph_cthw = sqrt(1-ph_sthw2)
      ph_Zmass2 = ph_Zmass**2
      ph_Wmass2 = ph_Wmass**2
      ph_Hmass2 = ph_Hmass**2

c     set mass windows around Z-mass peak in units of ph_Zwidth
c     It is used in the generation of the Born phase space
      masswindow = 30
      ph_Zmass2low=(ph_Zmass-masswindow*ph_Zwidth)**2
      ph_Zmass2high=(ph_Zmass+masswindow*ph_Zwidth)**2
      ph_ZmZw = ph_Zmass * ph_Zwidth

c     set mass window around W-mass peak in units of ph_Wwidth
c     It is used in the generation of the Born phase space
      masswindow = 30
      ph_Wmass2low=(ph_Wmass-masswindow*ph_Wwidth)**2
      ph_Wmass2high=(ph_Wmass+masswindow*ph_Wwidth)**2
      ph_WmWw = ph_Wmass * ph_Wwidth

c     Higgs parameters (not used)
      masswindow = 30
      if (powheginput("#higgsmasswindow").gt.0d0) then
         masswindow = powheginput("#higgsmasswindow")
      endif
      ph_Hmass2low=max(0d0,ph_Hmass-masswindow*ph_Hwidth)
      ph_Hmass2low=ph_Hmass2low**2
      ph_Hmass2high=(ph_Hmass+masswindow*ph_Hwidth)**2
      ph_HmHw = ph_Hmass * ph_Hwidth
      ph_Hmass2 = ph_Hmass**2    
      ph_unit_e = sqrt(4*pi*ph_alphaem)

      if(verbose) then
      write(*,*) '*************************************'
      write(*,*) 'Z mass = ',ph_Zmass
      write(*,*) 'Z width = ',ph_Zwidth
      write(*,*) 'W mass = ',ph_Wmass
      write(*,*) 'W width = ',ph_Wwidth
      write(*,*) 'H mass = ',ph_Hmass
      write(*,*) 'H width = ',ph_Hwidth
      write(*,*) '1/alphaem = ',1d0/ph_alphaem
      write(*,*) 'sthw2 = ',ph_sthw2
      write(*,*) '(unit_e)^2 = ',ph_unit_e**2   
      write(*,*) '(g_w)^2 = ',ph_unit_e*ph_unit_e/ph_sthw2   
      write(*,*) 'CKM matrix for events' 
      do i=1,3
         write(*,*) (ph_CKM(i,j),j=1,3)
      enddo
      write(*,*) '*************************************'
      endif      

c convert couplings into format needed by EW matrixelements:
      
      end


      

      subroutine inizialize_ph_CKM_matrix
      implicit none     
      include 'PhysPars.h'  
      integer i,j
      do i=1,6
         do j=1,6
            ph_CKM_matrix(i,j) = 0d0
         enddo
      enddo
      ph_CKM_matrix(1,2) = ph_CKM(1,1)
      ph_CKM_matrix(1,4) = ph_CKM(2,1)
      ph_CKM_matrix(1,6) = ph_CKM(3,1)
      ph_CKM_matrix(2,3) = ph_CKM(1,2)
      ph_CKM_matrix(2,5) = ph_CKM(1,3)
      ph_CKM_matrix(3,4) = ph_CKM(2,2)
      ph_CKM_matrix(3,6) = ph_CKM(3,2)
      ph_CKM_matrix(4,5) = ph_CKM(2,3)
      ph_CKM_matrix(5,6) = ph_CKM(3,3)
      do i=1,6
         do j=i+1,6
            ph_CKM_matrix(j,i) = ph_CKM_matrix(i,j)
         enddo
      enddo
      end


