c     HXSWG VBF 13.6 TeV study (STXS and fiducial set-ups): copied without
c     changes from the study's proVBFH copy, analysis/13.6TeV_analysis.f
c     (the analysis of its results/*.top), with its bins.h (copied to
c     analysis/bins.h, unchanged). Build: make ANALYSIS=hxswg136
c  The next subroutines open some histograms and prepare them 
c      to receive data 
c  You can substitute these with your favourite ones

      subroutine init_hist
      implicit none
      include 'pwhg_bookhist-multi.h'
      include 'PhysPars.h'
      include 'bins.h'
      real *8, parameter :: pi = 4.0d0*atan(1.0)
      integer i,j,k,l
      character * 100 histostr
      real *8, parameter ::  mjjbins(1:mjjnbins+1) =
     $     (/300d0, 500d0, 700d0, 900d0, 1100d0, 1d50/)
      real *8, parameter ::  ptHbins(1:pthnbins+1) =
     $     (/0d0, 80d0, 120d0, 260d0, 500d0, 850d0, 1d50/)
      real *8, parameter ::  delyjjbins(1:delyjjnbins+1) =
     $     (/2d0, 4d0, 5d0, 6d0, 7d0, 1d50/)
      real *8, parameter ::  delphijjbins(1:delphijjnbins+1) =
     $     (/0d0, pi/4d0, pi/2d0, 3d0*pi/4d0, pi/)
      real *8, parameter ::  njetbins(1:njetnbins+1) =
     $     (/2d0, 3d0, 4d0, 5d0, 6d0/)
      call inihists
      call setup_bin_names()
!     Inclusive histos
      call bookupeqbins('sig incl cuts (ptj > 20 GeV)',1d0, 0d0, 1d0)
      call bookupeqbins('sig incl cuts (ptj > 30 GeV)',1d0, 0d0, 1d0)
!     STXS bins. 
      do i = 1, stxsmjjnbins
         do j = 1,stxspthnbins
            do k = 1, stxspthjjnbins
               do l = 1, stxsdelphijjnbins
                  histostr = "STXS"//trim(stxsmjjstr(i))//trim(stxspthstr(j))
     $                 //trim(stxspthjjstr(k))//trim(stxsdelphijjstr(l))
                  print*, 'Booking ', trim(histostr)
                  call bookupeqbins(trim(histostr),1d0, 0d0, 1d0)
               enddo
            enddo
         enddo
      enddo
!     Our own bins. The two bins in each histograms correspond to
!     30>ptj>20 and ptj>30 respectively

!     Mjj vs pth
      i = 0
      do j = 1,pthnbins
         do k = 1,ptjnbins
            histostr = "HISTO"//trim(mjjstr(i))//trim(pthstr(j))//trim(ptjstr(k))
            print*, 'Booking ', trim(histostr)
            call bookup(trim(histostr),mjjnbins,mjjbins)
         enddo
      enddo

!     Mjj vs delta phi jj
      do j = 1,delphijjnbins
         do k = 1,ptjnbins
            histostr = "HISTO"//trim(mjjstr(i))//trim(delphijjstr(j))//trim(ptjstr(k))
            print*, 'Booking ', trim(histostr)
            call bookup(trim(histostr),mjjnbins,mjjbins)
         enddo
      enddo
!     Mjj vs delta y jj
      do j = 1,delyjjnbins
         do k = 1,ptjnbins
            histostr = "HISTO"//trim(mjjstr(i))//trim(delyjjstr(j))//trim(ptjstr(k))
            print*, 'Booking ', trim(histostr)
            call bookup(trim(histostr),mjjnbins,mjjbins)
         enddo
      enddo

!     ptH vs delta y jj
      do j = 1,delyjjnbins
         do k = 1,ptjnbins
            histostr = "HISTO"//trim(pthstr(i))//trim(delyjjstr(j))//trim(ptjstr(k))
            print*, 'Booking ', trim(histostr)
            call bookup(trim(histostr),pthnbins,pthbins)
         enddo
      enddo

!     Extra distribution for the PS simulations

!     ptH vs n jets
      do j = 1,njetnbins
         do k = 1,ptjnbins
            histostr = "HISTO"//trim(pthstr(i))//trim(njetstr(j))//trim(ptjstr(k))
            print*, 'Booking ', trim(histostr)
            call bookup(trim(histostr),pthnbins,pthbins)
         enddo
      enddo

      call bookupeqbins("ptHjj-ptj20", 20d0, 0d0, 1000d0)
      call bookupeqbins("ptHjj-ptj20-overflow", 1d0, 0d0, 1d0)
      call bookupeqbins("ptHjj-ptj30", 20d0, 0d0, 1000d0)
      call bookupeqbins("ptHjj-ptj30-overflow", 1d0, 0d0, 1d0)

!     Extra histograms for comparisons with Mathieu. These use the STXS
!     cuts.

      call bookupeqbins("ptj1-STXS", 20d0, 0d0, 2000d0)
      call bookupeqbins("ptj2-STXS", 20d0, 0d0, 2000d0)
      call bookupeqbins("mjj-STXS", 50d0, 0d0, 5000d0)
      call bookupeqbins("yjj-STXS", 0.25d0, -4.5d0, 4.5d0)
      call bookupeqbins("ptH-STXS", 20d0, 0d0, 2000d0)
      call bookupeqbins("ptHjj-STXS", 20d0, 0d0, 1000d0)
      call bookupeqbins("yj1-STXS", 0.25d0, -4.5d0, 4.5d0)
      call bookupeqbins("yj2-STXS", 0.25d0, -4.5d0, 4.5d0)
      call bookupeqbins("yH-STXS", 0.25d0, -4.5d0, 4.5d0)

      end

      subroutine setup_bin_names()
      implicit none
      include 'bins.h'

      stxsmjjstr(1) = "-MJJ-350-700"
      stxsmjjstr(2) = "-MJJ-700-1000"
      stxsmjjstr(3) = "-MJJ-1000-1500"
      stxsmjjstr(4) = "-MJJ-1500-INFTY"

      stxspthstr(1) = "-PTH-0-200"
      stxspthstr(2) = "-PTH-200-INFTY"

      stxspthjjstr(1) = "-PTHJJ-0-25"
      stxspthjjstr(2) = "-PTHJJ-25-INFTY"
      
      stxsdelphijjstr(1) = "-DPHIJJ-0-PIov2"
      stxsdelphijjstr(2) = "-DPHIJJ-PIov2-PI"

      mjjstr(0) = "-MJJ"
      mjjstr(1) = "-MJJ-300-500"
      mjjstr(2) = "-MJJ-500-700"
      mjjstr(3) = "-MJJ-700-900"
      mjjstr(4) = "-MJJ-900-1100"
      mjjstr(5) = "-MJJ-1100-INFTY"

      pthstr(0) = "-PTH"
      pthstr(1) = "-PTH-0-80"
      pthstr(2) = "-PTH-80-120"
      pthstr(3) = "-PTH-120-260"
      pthstr(4) = "-PTH-260-500"
      pthstr(5) = "-PTH-500-850"
      pthstr(6) = "-PTH-850-INFTY"

      delphijjstr(0) = "-DPHIJJ"
      delphijjstr(1) = "-DPHIJJ-0-PIov4"
      delphijjstr(2) = "-DPHIJJ-PIov4-PIov2"
      delphijjstr(3) = "-DPHIJJ-PIov2-3PIov4"
      delphijjstr(4) = "-DPHIJJ-3PIov4-PI"

      delyjjstr(0) = "-DYJJ"
      delyjjstr(1) = "-DYJJ-2-4"
      delyjjstr(2) = "-DYJJ-4-5"
      delyjjstr(3) = "-DYJJ-5-6"
      delyjjstr(4) = "-DYJJ-6-7"
      delyjjstr(5) = "-DYJJ-7-INFTY"

      njetstr(0) = "-NJETS"
      njetstr(1) = "-NJETS-2"
      njetstr(2) = "-NJETS-3"
      njetstr(3) = "-NJETS-4"
      njetstr(4) = "-NJETS-5"

      ptjstr(0) = "-PTJ"
      ptjstr(1) = "-PTJ-20"
      ptjstr(2) = "-PTJ-30"
      end
      
      subroutine user_analysis(dsig0)
      implicit none
      real * 8 dsig(7),dsig0
      include 'hepevt.h'
      include 'pwhg_math.h'  
      include 'nlegborn.h'
      include 'pwhg_kn.h'
      include 'LesHouches.h'
      include 'pwhg_weights.h'
      include 'bins.h'
c
c arrays to reconstruct jets
      integer maxtrack,maxjet
      parameter (maxtrack=4,maxjet=4)
      real *8 ptj(maxjet),yj(maxjet)
      integer mu,njets,ijet
      real * 8 R,ptmin_fastkt,palg,ptjbin
      integer  k, j, i
      real * 8 pH(0:3),ptH,yH,inv_mH,pj(0:3,4),pHjj(0:3), mHjj
      real * 8 ptHjj
      real * 8 rsepn,getrapidity0,mjj,azi,mjj2,rsepn4
      external rsepn,getrapidity0,mjj,azi,mjj2,rsepn4
      real * 8 Rmin,delphijj,invmjj,delyjj,rapj1j3,rapj2j3
      logical pass_cuts, passed_cuts_jet
c     we need to tell to the this analysis file which program is running it
      character * 6 WHCPRG
      common/cWHCPRG/WHCPRG
C      data WHCPRG/'NLO   '/
      double precision getdelphiv4,powheginput
      external getdelphiv4,powheginput
      logical passed_cuts_vbf
      real * 8 ptjmjjbinsize,bin,s
      common/cptjmjjbinsize/ptjmjjbinsize
      save /cptjmjjbinsize/
      logical ini,ptjsqrtsrwgt
      data ini/.true./
      save ini,ptjsqrtsrwgt
      character * 100 histostr(5)

c     COMMON block to cut on phasespace. Contains the cuts.
      include 'pwhg_flg.h'
      include 'phspcuts.h'
      real * 8 parallelstage
c===============================================
      if(ini) then
         if(WHCPRG.eq.'NLO   '.or.WHCPRG.eq.'LHE   ') then
            weights_num=0
         endif
         if(weights_num.eq.0) then
            call setupmulti(1)
         else
            call setupmulti(weights_num)
         endif
      endif !ini

      dsig(:)=0d0
      if(weights_num.eq.0) then
         dsig(1)=dsig0
      else
         dsig(1:weights_num)=weights_val(1:weights_num)
      endif

      if(sum(abs(dsig)).eq.0) return
      
      if (ini) then
C     The routines "setup_vbf_cuts", "buildjets" and "vbfcuts" are also
C     called inside the phspcuts analysis. Hence it is important that
C     changes take place inside these routines rather than outside, as
C     otherwise the phase space cuts will no longer work.
         call setup_vbf_cuts()
C     Any additional cut parameters can be added here, like a central
C     jet veto. VBF like cuts should enter in the above subroutine, as
C     these are also used by the phase space cuts routine.
         write(*,*) '**************************************************'
         write(*,*) '**************************************************'
         write(*,*) '                ANALYSIS CUTS                     '
         write(*,*) '**************************************************'
         write(*,*) '**************************************************'
         write(*,*) 'ptalljetmin = ',ptalljetmin
         write(*,*) 'ptjetmin = ',ptjetmin
         write(*,*) 'yjetmax = ',yjetmax
         write(*,*) 'mjjmin = ',mjjmin 
         write(*,*) 'deltay_jjmin = ',deltay_jjmin
         write(*,*) 'Rsep_jjmin = ',Rsep_jjmin
         write(*,*) 'jet_opphem = ',jet_opphem 
         write(*,*) '**************************************************'
         write(*,*) '**************************************************'
         ptjsqrtsrwgt = powheginput('#ptjsqrtsrwgt').eq.1d0
         ini = .false.
      endif

c     Build jets passing ptalljetmin cut and yjetmax cut.
      call buildjets(pj,njets,ptj,yj)
c     Check if jets satisfy VBF cuts
      call vbfcuts(pj,njets,ptj,yj,passed_cuts)
      
c     At this point we have at least two jets with ptj > 20 and the two hardest jets satisfying mjj>300
c     Return already here if VBF cuts are not passed. 
      if(.not.passed_cuts) return

c     Higgs momentum
      do mu=1,3
         pH(mu) = phep(mu,3)
      enddo
      pH(0) = phep(4,3)

      ptH = sqrt(pH(1)**2+pH(2)**2)
      yH = getrapidity0(pH) 
      
c     Now compute some jet variables
      invmjj       = mjj(pj(0,1),pj(0,2)) !Invariant dijet mass
      delphijj     = getdelphiv4(pj(:,1),pj(:,2))
      delyjj       = abs(yj(1) - yj(2)) ! Rapidity separations

c     ptHjj
      pHjj(:) = pH(:) + pj(:,1) + pj(:,2)
      ptHjj = sqrt(pHjj(1)**2+pHjj(2)**2)

      call select_stxs_bin(ptj,invmjj,pth,pthjj,delphijj,ptjbin,histostr)
      ! Fill STXS bins if we pass the cut
      if(invmjj.gt.350d0.and.abs(yH).lt.2.5d0.and.ptj(2).gt.30d0) then
         call filld(trim(histostr(1)),ptjbin,dsig)
         call filld("ptj1-STXS", ptj(1)        ,dsig)
         call filld("ptj2-STXS", ptj(2)        ,dsig)
         call filld("mjj-STXS",  invmjj        ,dsig)
         call filld("yjj-STXS",  yj(1) - yj(2) ,dsig)
         call filld("ptH-STXS",  ptH           ,dsig)
         call filld("ptHjj-STXS",ptHjj         ,dsig)
         call filld("yj1-STXS",  yj(1)         ,dsig)
         call filld("yj2-STXS",  yj(2)         ,dsig)
         call filld("yH-STXS",   yH            ,dsig)

      endif

c     Additional cuts needed for this setup. We rebuild the jets, which
c     is not the most efficient, but perhaps the cleanest?
      deltay_jjmin = 2d0
      yjetmax = 4.7d0

c     Build jets passing ptalljetmin cut and yjetmax cut.
      call buildjets(pj,njets,ptj,yj)
c     Check if jets satisfy VBF cuts
      call vbfcuts(pj,njets,ptj,yj,passed_cuts)

      deltay_jjmin = 0d0 ! Reset cut
      yjetmax = 4.7d10 ! Reset cut

      if(.not.passed_cuts) return

c     Now compute some jet variables
      invmjj       = mjj(pj(0,1),pj(0,2)) !Invariant dijet mass
      delphijj     = getdelphiv4(pj(:,1),pj(:,2))
      delyjj       = abs(yj(1) - yj(2)) ! Rapidity separations

c     ptHjj
      pHjj(:) = pH(:) + pj(:,1) + pj(:,2)
      ptHjj = sqrt(pHjj(1)**2+pHjj(2)**2)

c     Always fill inclusive cross section         
      call filld('sig incl cuts (ptj > 20 GeV)',0.5d0, dsig)
      if(ptj(2).gt.30d0) call filld('sig incl cuts (ptj > 30 GeV)',0.5d0, dsig)
      
      call select_histo_bin(ptj,invmjj,pth,delphijj,delyjj,njets
     $     ,histostr)

      call filld(trim(histostr(1))//trim(ptjstr(1)),invmjj,dsig)
      call filld(trim(histostr(2))//trim(ptjstr(1)),invmjj,dsig)
      call filld(trim(histostr(3))//trim(ptjstr(1)),invmjj,dsig)
      call filld(trim(histostr(4))//trim(ptjstr(1)),pth,dsig)
      call filld(trim(histostr(5))//trim(ptjstr(1)),pth,dsig)

      if(ptj(2).gt.30d0) then
         call filld(trim(histostr(1))//trim(ptjstr(2)),invmjj,dsig)
         call filld(trim(histostr(2))//trim(ptjstr(2)),invmjj,dsig)
         call filld(trim(histostr(3))//trim(ptjstr(2)),invmjj,dsig)
         call filld(trim(histostr(4))//trim(ptjstr(2)),pth,dsig)
         call filld(trim(histostr(5))//trim(ptjstr(2)),pth,dsig)
      endif
      
      if(ptHjj.lt.1000d0) then
         call filld('ptHjj-ptj20',ptHjj,dsig)
         if(ptj(2).gt.30d0) call filld('ptHjj-ptj30',ptHjj,dsig)
      else
         call filld('ptHjj-ptj20-overflow',0.5d0,dsig)
         if(ptj(2).gt.30d0) call filld('ptHjj-ptj30-overflow',0.5d0,dsig)
      endif
      end

      subroutine select_stxs_bin(ptj,invmjj,pth,pthjj,delphijj,ptjbin
     $     ,histostr)
      implicit none
      include 'bins.h'
      real * 8 pi
      parameter (pi = 4d0*atan(1d0))
      real * 8 ptj(2),invmjj,pth,pthjj,ptjbin,delphijj
      integer invmjjbin,pthbin,pthjjbin,delphijjbin
      character * 100 histostr(1)

      ! Select which bin to fill - always do 30 GeV, so just one bin
      ptjbin = 0.5d0
      
      if(delphijj.lt.0.5d0*pi) then
         delphijjbin = 1
      else 
         delphijjbin = 2
      endif
      
      if(ptHjj.lt.25d0) then
         pthjjbin = 1
      else
         pthjjbin = 2
      endif

      if(ptH.lt.200d0) then
         pthbin = 1
      else
         pthbin = 2
      endif

      if(invmjj.gt.350d0.and.invmjj.lt.700d0) then
         invmjjbin = 1
      elseif(invmjj.gt.700d0.and.invmjj.lt.1000d0) then
         invmjjbin = 2
      elseif(invmjj.gt.1000d0.and.invmjj.lt.1500d0) then
         invmjjbin = 3
      elseif(invmjj.gt.1500d0) then
         invmjjbin = 4
      endif

      histostr(1) = "STXS"//trim(stxsmjjstr(invmjjbin))
     $     //trim(stxspthstr(pthbin))
     $     //trim(stxspthjjstr(pthjjbin))
     $     //trim(stxsdelphijjstr(delphijjbin))
      end
      
      subroutine select_histo_bin(ptj,invmjj,pth,delphijj,delyjj
     $     ,njets,histostr)
      implicit none
      include 'bins.h'
      real * 8 pi
      parameter (pi = 4d0*atan(1d0))
      real * 8 ptj(2),invmjj,pth,delphijj,delyjj
      integer invmjjbin,pthbin,pthjjbin,delphijjbin,delyjjbin,njets,
     $     njetsbin, ptjbin
      character * 100 histostr(5)

      if(delphijj.lt.0.25d0*pi) then
         delphijjbin = 1
      elseif(delphijj.lt.0.5d0*pi) then
         delphijjbin = 2
      elseif(delphijj.lt.0.75d0*pi) then
         delphijjbin = 3
      else 
         delphijjbin = 4
      endif
      
      if(ptH.lt.80d0) then
         pthbin = 1
      elseif(ptH.lt.120d0) then
         pthbin = 2
      elseif(ptH.lt.260d0) then
         pthbin = 3
      elseif(ptH.lt.500d0) then
         pthbin = 4
      elseif(ptH.lt.850d0) then
         pthbin = 5
      else
         pthbin = 6
      endif

!      if(invmjj.gt.300d0.and.invmjj.lt.500d0) then
!         invmjjbin = 1
!      elseif(invmjj.lt.700d0) then
!         invmjjbin = 2
!      elseif(invmjj.lt.900d0) then
!         invmjjbin = 3
!      elseif(invmjj.lt.1100d0) then
!         invmjjbin = 4
!      elseif(invmjj.gt.1100d0) then
!         invmjjbin = 5
!      endif

      if(delyjj.gt.2d0.and.delyjj.lt.4d0) then
         delyjjbin = 1
      elseif(delyjj.lt.5d0) then
         delyjjbin = 2
      elseif(delyjj.lt.6d0) then
         delyjjbin = 3
      elseif(delyjj.lt.7d0) then
         delyjjbin = 4
      elseif(delyjj.gt.7d0) then
         delyjjbin = 5
      endif

      if(njets.eq.2) then
         njetsbin = 1
      elseif(njets.eq.3) then
         njetsbin = 2
      elseif(njets.eq.4) then
         njetsbin = 3
      elseif(njets.eq.5) then
         njetsbin = 4
      endif
      
      histostr(1) = "HISTO"//trim(mjjstr(0))
     $     //trim(pthstr(pthbin))

      histostr(2) = "HISTO"//trim(mjjstr(0))
     $     //trim(delphijjstr(delphijjbin))

      histostr(3) = "HISTO"//trim(mjjstr(0))
     $     //trim(delyjjstr(delyjjbin))

      histostr(4) = "HISTO"//trim(pthstr(0))
     $     //trim(delyjjstr(delyjjbin))

      histostr(5) = "HISTO"//trim(pthstr(0))
     $     //trim(njetstr(njetsbin))


      end
      
      subroutine particle_identif(HWW,HZZ)
      implicit none
      integer pdg_Higgs,pdg_Z,pdg_W,HZZ,HWW
      pdg_Higgs = 25
      pdg_Z = 23
      pdg_W = 24      
c     build an identifier for Higgs production in WW and ZZ fusion 
      HWW = 10000*pdg_W + pdg_Higgs
      HZZ = 10000*pdg_Z + pdg_Higgs
      end

      subroutine getrapidity(p,y)
      implicit none
      real * 8 p(4),y
      y=0.5d0*log((p(4)+p(3))/(p(4)-p(3)))
      end

      function getrapidity0(p)
      implicit none
      real * 8 p(0:3),getrapidity0
      getrapidity0=0.5d0*log((p(0)+p(3))/(p(0)-p(3)))
      end

      subroutine getinvmass(p,m)
      implicit none
      real * 8 p(4),m
      m=dsqrt(abs(p(4)**2-p(1)**2-p(2)**2-p(3)**2))
      end
      
      function azi(p)
      implicit none
      include 'pwhg_math.h'  
      real * 8 azi,p(0:3)
      azi = atan(p(2)/p(1))
      if (p(1).lt.0d0) then
         if (azi.gt.0d0) then               
            azi = azi - pi
         else
            azi = azi + pi
         endif
      endif    
      end

c     calculate the separation in the lego plot between the two momenta
c     p1 and p2
      function rsepn(p1,p2)
      implicit none
      include 'pwhg_math.h'  
      real * 8 rsepn,p1(0:3),p2(0:3)
      real * 8 y1,phi1,y2,phi2
      real * 8 delphi
      real * 8 getrapidity0,azi
      external getrapidity0,azi

      phi1 = azi(p1)   
      phi2 = azi(p2)
      y1 = getrapidity0(p1)
      y2 = getrapidity0(p2)

      delphi = abs(phi1-phi2)
      if (delphi.gt.pi) then
         delphi = 2*pi-delphi
      endif
      if (delphi.lt.0 .or. delphi.gt.pi) then
         print*,' problem in rsepn. delphi = ',delphi
      endif
      rsepn = sqrt( (y1-y2)**2 + delphi**2 )
      end

c     calculate the separation in the lego plot between the two momenta
c     p1 and p2 for pj(1:4)
      function rsepn4(p1,p2)
      implicit none
      include 'pwhg_math.h'  
      real * 8 rsepn4,p1(4),p2(4)
      real * 8 y1,phi1,y2,phi2
      real * 8 delphi
      real * 8 azi4
      external azi4

      phi1 = azi4(p1)   
      phi2 = azi4(p2)
      call getrapidity(p1,y1)
      call getrapidity(p2,y2)

      delphi = abs(phi1-phi2)
      if (delphi.gt.pi) then
         delphi = 2*pi-delphi
      endif
      if (delphi.lt.0 .or. delphi.gt.pi) then
         print*,' problem in rsepn. delphi = ',delphi
      endif
      rsepn4 = sqrt( (y1-y2)**2 + delphi**2 )
      end
      
      function azi4(p)
      implicit none
      include 'pwhg_math.h'  
      real * 8 azi4,p(4)
      azi4 = atan(p(2)/p(1))
      if (p(1).lt.0d0) then
         if (azi4.gt.0d0) then               
            azi4 = azi4 - pi
         else
            azi4 = azi4 + pi
         endif
      endif    
      end

c mjj^2 = (p1+p2)^2 = p1^2 + p2^2 + 2*dotp(p1,p2)
      function mjj(p1,p2)
      implicit none
      real * 8 mjj,p1(0:3),p2(0:3)
      real * 8 p(0:3)
      integer mu
      do mu=0,3
         p(mu)=p1(mu)+p2(mu)
      enddo
      mjj = sqrt(abs(p(0)**2-p(1)**2-p(2)**2-p(3)**2))
      end

c     find the first "nhardjets" hardest jets in pjet (that contains njets)
c     and return their position.
c     foundhardjets is the number of found hard jets (.le.nhardjets)
      subroutine find_hardest_jets(njets,pjet,nhardjets,
     #     foundhardjets,jj)
      implicit none
      integer njets
      real *8 pjet(4,njets) 
      integer nhardjets,jj(nhardjets)
      real * 8 ptj(nhardjets),pt
      integer ijet,hjet,foundhardjets,i
      logical is_i_in_array
      external is_i_in_array

      if (njets.eq.0) then
         write(*,*) 'WARNING!!!!!!!!!!!  EMPTY  PJET ARRAY'
         nhardjets=0
         return
      endif

      do hjet=1,nhardjets
         jj(hjet)=0d0
         ptj(hjet)=0d0
      enddo
      foundhardjets=1
      do ijet=1,njets   
         pt=sqrt(pjet(1,ijet)**2 + pjet(2,ijet)**2)
         do hjet=1,min(foundhardjets,nhardjets)
            if (pt.gt.ptj(hjet).and.
     $           .not.is_i_in_array(nhardjets,ijet,jj)) then
               foundhardjets = foundhardjets + 1
               do i=nhardjets,hjet+1,-1
                  ptj(i)=ptj(i-1)
                  jj(i)=jj(i-1)
               enddo
               ptj(hjet)=pt
               jj(hjet)=ijet
            endif
         enddo
      enddo
c     set number of jets found
      foundhardjets = min(foundhardjets-1,nhardjets)
      end

      function is_i_in_array(nhardjets,i,jj)
      implicit none
      logical is_i_in_array
      integer nhardjets,i,jj(nhardjets)
      integer j
      is_i_in_array = .false.
      do j=1,nhardjets
         if (i.eq.jj(j)) then
            is_i_in_array = .true.
            return
         endif
      enddo
      end


c----- ------ ----- ----- ----- ----- ----- ----- ----- -----                                           
      function getdelphi(ptx1,pty1,ptx2,pty2)
      implicit none
      real*8 getdelphi,ptx1,pty1,ptx2,pty2,tiny,pt1,pt2,tmp
      parameter (tiny=1.d-5)
c                                                                                                       
      pt1=sqrt(ptx1**2+pty1**2)
      pt2=sqrt(ptx2**2+pty2**2)
      if(pt1.ne.0.d0.and.pt2.ne.0.d0)then
        tmp=ptx1*ptx2+pty1*pty2
	tmp=tmp/(pt1*pt2)
	if(abs(tmp).gt.1.d0+tiny)then
          write(*,*)'Cosine larger than 1'
          stop
        elseif(abs(tmp).ge.1.d0)then
          tmp=sign(1.d0,tmp)
        endif
        tmp=acos(tmp)
      else
        tmp=1.d8
      endif
      getdelphi=tmp
      return
      end

c----- ------ ----- ----- ----- ----- ----- ----- ----- -----                                           
      function getdelphiv4(p1,p2)
      implicit none
      real*8 getdelphiv4,p1(0:3),p2(0:3)
      real*8 getdelphi
c                                                                                                       
      getdelphiv4=getdelphi(p1(1),p1(2),
     #                      p2(1),p2(2))
      return
      end



      subroutine setup_vbf_cuts()
      implicit none
      double precision powheginput
      external powheginput

      include 'phspcuts.h'

      jet_opphem = .false.
      
c replace default cuts by input values (if active):
      ptalljetmin=powheginput('#ptalljetmin')
      ptjetmin=powheginput('#ptjetmin')
      mjjmin=powheginput('#mjjmin')  
      deltay_jjmin=powheginput('#deltay_jjmin')
      yjetmax=powheginput('#yjetmax') 
      ptjsqrtsmax=powheginput('#ptjsqrtsmax')
      Rsep_jjmin=powheginput('#Rsep_jjmin')
      if(powheginput('#jet_opphem').eq.1d0) jet_opphem=.true.

      if(ptalljetmin.gt.ptjetmin) then
         print*, 'Inconsistent jet cuts.', ptalljetmin, ptjetmin
         stop
      endif
      end

      subroutine buildjets(pj,njets,ptj,yj)
      implicit none
      include 'hepevt.h'
      include 'phspcuts.h'
      integer maxtrack,maxjet
      parameter (maxtrack=4,maxjet=4)
      character * 6 WHCPRG
      common/cWHCPRG/WHCPRG

!     Output
      double precision pj(0:3,maxjet)
!     Internal
      integer mu, njets, njets_all, ntracks, ijet, j
      integer jetvec(maxtrack)
      double precision pjet(4,maxjet), pjet_all(4,maxjet)
      double precision ptrack(4,maxtrack)
      double precision ptj(maxjet),yj(maxjet)
      double precision ptj_all(maxjet),yj_all(maxjet)
      double precision R, ptmin_fastkt, palg


      ptrack = 0d0
      jetvec = 0
      pjet = 0d0
      pjet_all = 0d0
      njets=0
      njets_all = 0
      ntracks = nhep - 3 ! 2 initial states and one Higgs 
      if(WHCPRG.eq.'PROJEC') ntracks = 2

      do mu=1,4
         ptrack(mu,1:ntracks)=phep(mu,4:nhep)
      enddo
      
      
************************************************************************
*     siscone algorithm
**********************************************************************
c     R =  radius parameter
c     f =  overlapping fraction
c.....run the clustering        
c      call fastjetsiscone(ptrack,ntracks,R,f,pjet,njets) 
************************************************************************
*     fastkt algorithm
**********************************************************************
      R = Rsep_jjmin          
      ptmin_fastkt = 0d0
      palg = -1d0
c     -1 is anti_kt, +1 is kt

      call fastjetppgenkt(ptrack,ntracks,r,palg,ptmin_fastkt,pjet,njets,
     $                        jetvec)

      if(njets.gt.4) then
         print*, 'njets out of bounds!!', njets
         stop
      endif
c     Find the number of jets inside the detector and passing the jet pt
c     cut.
      j=0
      pjet_all=pjet
      njets_all=njets
      pjet=0d0
      njets=0d0
      ptj = 0d0
      yj = 0d0
      do ijet=1,njets_all
         ptj_all(ijet) = sqrt(pjet_all(1,ijet)**2 + pjet_all(2,ijet)**2)
         call getrapidity(pjet_all(:,ijet),yj_all(ijet))
         if(ptj_all(ijet).gt.ptalljetmin.and.
     $        abs(yj_all(ijet)).lt.yjetmax) then
            j=j+1
            pjet(:,j)=pjet_all(:,ijet)
            ptj(j) = ptj_all(ijet)
            yj(j) = yj_all(ijet)
         endif
      enddo
      njets=j

      pj = 0d0

      do ijet=1,njets
         do mu=1,3
            pj(mu,ijet)=pjet(mu,ijet)
         enddo
         pj(0,ijet)=pjet(4,ijet)
      enddo

      end

      subroutine vbfcuts(pj,njets,ptj,yj,passed)
      implicit none
      include 'phspcuts.h'
      logical passed
      double precision pj(0:3,4), ptj(4),yj(4)
      integer njets
      double precision invmjj, delyjj,yj1dotyj2, mjj
      external mjj
      passed = .false.
      
      if(njets.lt.2) return

      invmjj = mjj(pj(0,1),pj(0,2))
      delyjj = abs(yj(1) - yj(2))
      yj1dotyj2 = yj(1) * yj(2)

      passed = (min(pTj(1),pTj(2)).gt.ptjetmin) .and.
     $        (max(abs(yj(1)),abs(yj(2))).lt.yjetmax) .and.
     $        (invmjj.gt.mjjmin) .and.
     $        (delyjj.gt.deltay_jjmin)

      if (jet_opphem) then
         passed = passed .and. 
     $        (yj1dotyj2.lt.0)
      endif
      
      end

      subroutine jetcuts(pj,njets,ptj,yj,passed)
      implicit none
      include 'phspcuts.h'
      logical passed
      double precision pj(0:3,4), ptj(4),yj(4)
      integer njets
      passed = .false.
      
      if(njets.lt.2) return

      passed = (min(pTj(1),pTj(2)).gt.ptjetmin) .and.
     $        (max(abs(yj(1)),abs(yj(2))).lt.yjetmax)
      
      end

