c  The next subroutines open some histograms and prepare them 
c      to receive data 
c  You can substitute these with your favourite ones

      subroutine init_hist
      implicit none
      include 'pwhg_bookhist-multi.h'
      include 'PhysPars.h'
      real * 8 pi
      parameter(pi = 3.141592653589793D0)

      call inihists

      call bookupeqbins('sig incl cuts',1d0, 0d0, 1d0)

      call bookupeqbins('sig(all VBF cuts 2 jets)',1d0, 0d0, 1d0)
      call bookupeqbins('sig(all VBF cuts 3 jets)',1d0, 0d0, 1d0)
      call bookupeqbins('sig(all VBF cuts 4 jets)',1d0, 0d0, 1d0)

      call bookupeqbins('ptj1',20d0, 0d0, 600d0)
      call bookupeqbins('ptj2',10d0, 0d0, 320d0)
      call bookupeqbins('ptH1',30d0, 0d0, 600d0)
      call bookupeqbins('ptH2',30d0, 0d0, 600d0)
      call bookupeqbins('ptH_h',30d0, 0d0, 800d0)
      call bookupeqbins('ptH_s',30d0, 0d0, 600d0)

      call bookupeqbins('yj1',5d-1, -6d0, 6d0)
      call bookupeqbins('yj2',5d-1, -6d0, 6d0)
      call bookupeqbins('yH1',5d-1, -5d0, 5d0)
      call bookupeqbins('yH2',5d-1, -5d0, 5d0)
      call bookupeqbins('yH_h',5d-1, -5d0, 5d0)
      call bookupeqbins('yH_s',5d-1, -5d0, 5d0)
      
!      call bookupeqbins('mH1',200d0, 0d0, 2000d0)
!      call bookupeqbins('mH2',200d0, 0d0, 2000d0)
      
      call bookupeqbins('phi(H1,H2)',pi/10d0, 0d0, pi)

      call bookupeqbins('yHH',5d-1, -5d0, 5d0)
      call bookupeqbins('mHH',200d0, 0d0, 7000d0)
      call bookupeqbins('ptHH',30d0, 0d0, 600d0)

      call bookupeqbins('M_jj',200d0, 0d0, 7000d0)
      call bookupeqbins('rap(j1,j2)',5d-1, 0d0, 10d0)
      call bookupeqbins('phi(j1,j2)',pi/10d0, 0d0, pi)
      call bookupeqbins('R(j1,j2)',5d-1, 0d0, 10d0)

      call bookupeqbins('ptj3',5d0, 0d0, 100d0)
      call bookupeqbins('yj3',5d-1, -6d0, 6d0)
      call bookupeqbins('y*j3',5d-1, -6d0, 6d0)
      call bookupeqbins('min{rap(j1,j3),rap(j3,j2)}',2d-1, -5d0, 5d0)
      call bookupeqbins('rapj1j3',5d-1, 0d0, 10d0)
      call bookupeqbins('rapj2j3',5d-1, 0d0, 10d0)

      call bookupeqbins('ptj4',5d0, 0d0, 100d0)
      call bookupeqbins('yj4',5d-1, -6d0, 6d0)
      call bookupeqbins('rapj1j4',5d-1, 0d0, 10d0)
      call bookupeqbins('rapj2j4',5d-1, 0d0, 10d0)
      call bookupeqbins('rapj3j4',5d-1, 0d0, 10d0)
      
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
c
c arrays to reconstruct jets
      integer maxtrack,maxjet
      parameter (maxtrack=4,maxjet=4)
      real *8 ptj(maxjet),yj(maxjet)
      integer mu,njets,ijet
      real * 8 R,ptmin_fastkt,palg,rsepj1j2
      integer  k, j
      real * 8 pH1(0:3),pH2(0:3),yH1,yH2,mH1,mH2,mHH,pj(0:3,4)
      real * 8  mHHjj,pHHjj(0:3),ptHH,yHH,ptH1,ptH2
      real * 8 ptHs, ptHd, yHs, yHd, phiHH
      real * 8 y3star, min_yj1yj2yj3, ht_jets, ht_all
      real * 8 rsepn,getrapidity0,mjj,azi,mjj2,rsepn4
      external rsepn,getrapidity0,mjj,azi,mjj2,rsepn4
      real * 8 Rmin,delphi_jj,invmjj,rapj1j2,rapj1j3,rapj2j3,rapj1j4
     $     ,rapj2j4,rapj3j4
      logical pass_cuts
      real * 8 tchan_vec(1:4), uchan_vec(1:4), t, u, Qtilde,s
      real * 8 ptjsqrtsmax, powheginput
      save ptjsqrtsmax
c     we need to tell to the this analysis file which program is running it
      character * 6 WHCPRG
      common/cWHCPRG/WHCPRG
C      data WHCPRG/'NLO   '/
      double precision getdelphiv4
      external getdelphiv4

      logical ini
      data ini/.true./
      save ini

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
C     othersie the phase space cuts will no longer work.
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
         ini = .false.
         ptjsqrtsmax = 1d0
         if(powheginput('#ptjsqrtsmax').gt.0d0) then
            ptjsqrtsmax = powheginput('#ptjsqrtsmax')
         endif
      endif

c     Build jets passing ptalljetmin cut and yjetmax cut.
      call buildjets(pj,njets,ptj,yj)
!      print*, 'njets', njets
c     Check if jets satisfy VBF cuts
      call vbfcuts(pj,njets,ptj,yj,passed_cuts)
      
      tchan_vec = phep(1:4,5) + phep(1:4,3)
      uchan_vec = phep(1:4,6) + phep(1:4,3)
      
      t = sqrt(tchan_vec(1)**2 + tchan_vec(2)**2)
      u = sqrt(uchan_vec(1)**2 + uchan_vec(2)**2)
      Qtilde = max(t,u)
      Qtilde = max(ptj(1),Qtilde)
!     Fill inclusive distributions
      s = 2d0 * (phep(4,1)*phep(4,2) - phep(3,1)*phep(3,2))
      if(Qtilde/sqrt(s).gt.ptjsqrtsmax) return
c     Always fill inclusive cross section
      call filld('sig incl cuts',0.5d0, dsig)
c     Return already here if VBF cuts are not passed. 
      if(.not.passed_cuts) return

c     Higgs momentum
      do mu=1,3
         pH1(mu) = phep(mu,3)
         pH2(mu) = phep(mu,4)
      enddo
      pH1(0) = phep(4,3)
      pH2(0) = phep(4,4)

      ptH1 = sqrt(pH1(1)**2+pH1(2)**2)
      ptH2 = sqrt(pH2(1)**2+pH2(2)**2)
      ptHH = sqrt((pH1(1)+pH2(1))**2+(pH1(2)+pH2(2))**2)
      call getrapidity(phep(1,3),yH1)
      call getinvmass(phep(1,3),mH1)
      call getrapidity(phep(1,4),yH2)
      call getinvmass(phep(1,4),mH2)
      call getrapidity(phep(:,3)+phep(:,4),yHH)
      call getinvmass(phep(:,3)+phep(:,4),mHH)
      
c     Now compute some jet variables
      invmjj        = mjj(pj(0,1),pj(0,2)) !Invariant dijet mass
      rapj1j2       = yj(1) - yj(2) ! Rapidity separations
      Rsepj1j2      = rsepn(pj(0,1),pj(0,2)) !R-separation between tagging jets
      delphi_jj     = getdelphiv4(pj(:,1),pj(:,2))

      phiHH = getdelphiv4(pH1,pH2)
      if (ptH1.gt.ptH2) then
         ptHd = ptH1
         yHd  = yH1
         ptHs = ptH2
         yHs  = yH2
      else
         ptHd = ptH2
         yHd  = yH2
         ptHs = ptH1
         yHs  = yH1
      endif
     
c     compute invariant mass of the Hjj system
      do mu=0,3
         pHHjj(mu) = pH1(mu)+pH2(mu)+pj(mu,1)+pj(mu,2)
      enddo
      mHHjj = sqrt(abs(phhjj(0)**2-phhjj(1)**2-phhjj(2)**2-phhjj(3)**2))
c     
      
c     compute HT of jet system and of entire system
      ht_jets = sum(ptj(:))
      ht_all  = ht_jets + ptH1 + ptH2

      if(njets.gt.2) then
         rapj1j3       = yj(1) - yj(3)
         rapj2j3       = yj(2) - yj(3)
         y3star        = yj(3) - (yj(1)+yj(2))*0.5d0
         
         min_yj1yj2yj3 = min(abs(rapj1j3),abs(rapj2j3))
c     If the third jet is in the j1,j2 window the two rapidity separations
c     will have opposite sign. Otherwise they will have the same sign.
c     We assign plus if the third jet is in the window and minus otherwise.
         min_yj1yj2yj3 = sign(min_yj1yj2yj3,-rapj1j3*rapj2j3) 
      endif
      if(njets.gt.3) then
         rapj1j4       = yj(1) - yj(4)
         rapj2j4       = yj(2) - yj(4)
         rapj3j4       = yj(3) - yj(4)
      endif
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCAPPLY CUTS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c     
c     If 3/4 jets present we have a pt cut.
      pass_cuts=passed_cuts
       
      if (pass_cuts) then
         call filld('sig(all VBF cuts 2 jets)',0.5d0, dsig)
         call filld('ptj1',ptj(1), dsig)
         call filld('ptj2',ptj(2), dsig)
         call filld('ptH1',ptH1, dsig)
         call filld('ptH2',ptH2, dsig)
         call filld('ptH_h',ptHd,  dsig)
         call filld('ptH_s',ptHs,  dsig)
         
         call filld('yj1',yj(1),dsig)
         call filld('yj2',yj(2),dsig)
         call filld('yH1',yH1,dsig)
         call filld('yH2',yH2,dsig)
         call filld('yH_h',yHd,dsig)
         call filld('yH_s',yHs,dsig)

!         call filld('mH1',mH1,dsig)
!         call filld('mH2',mH2,dsig)
         
         call filld('yHH',yHH,dsig)
         call filld('mHH',mHH,dsig)
         call filld('ptHH',ptHH,dsig)
         call filld('phi(H1,H2)',phiHH,dsig)

         call filld('M_jj',invmjj,dsig)
         call filld('rap(j1,j2)',abs(rapj1j2),dsig)
         call filld('phi(j1,j2)',delphi_jj,dsig)
         call filld('R(j1,j2)',Rsepj1j2,dsig)
      endif
      
      if(njets.ge.3) then
         pass_cuts=pass_cuts.and.(ptj(3).gt.ptalljetmin)
      endif         
      
      
      if((njets.ge.3).and.pass_cuts) then 
         call filld('sig(all VBF cuts 3 jets)',0.5d0, dsig)
         call filld('ptj3',ptj(3),dsig)
         call filld('yj3',yj(3),dsig)
         call filld('y*j3',y3star,dsig)
         call filld('min{rap(j1,j3),rap(j3,j2)}',min_yj1yj2yj3,dsig)
         call filld('rapj1j3',abs(rapj1j3),dsig)
         call filld('rapj2j3',abs(rapj2j3),dsig)
      endif
      
      if(njets.ge.4) then
         pass_cuts=pass_cuts.and.(ptj(4).gt.ptalljetmin)
      endif
      
      if((njets.ge.4).and.pass_cuts) then 
         call filld('sig(all VBF cuts 4 jets)',0.5d0, dsig)
         call filld('ptj4',ptj(4),dsig)
         call filld('yj4',yj(4),dsig)
         call filld('rapj1j4',abs(rapj1j4),dsig)
         call filld('rapj2j4',abs(rapj2j4),dsig)
         call filld('rapj3j4',abs(rapj3j4),dsig)
      endif
      
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
!     VBFHHMOD
      ntracks = nhep - 4 ! 2 initial states and two Higgs 
!      print*, 'ntracks', ntracks
      if(WHCPRG.eq.'PROJEC') ntracks = 2
!      VBFHHMOD
      do mu=1,4
         ptrack(mu,1:ntracks)=phep(mu,5:nhep)
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
      double precision invmjj, rapj1j2,yj1dotyj2, mjj
      external mjj
      passed = .false.
      
      if(njets.lt.2) return

      invmjj = mjj(pj(0,1),pj(0,2))
      rapj1j2 = abs(yj(1) - yj(2))
      yj1dotyj2 = yj(1) * yj(2)

      passed = (min(pTj(1),pTj(2)).gt.ptjetmin) .and.
     $        (max(abs(yj(1)),abs(yj(2))).lt.yjetmax) .and.
     $        (invmjj.gt.mjjmin) .and.
     $        (rapj1j2.gt.deltay_jjmin)

      if (jet_opphem) then
         passed = passed .and. 
     $        (yj1dotyj2.lt.0)
      endif
      
      end
    
