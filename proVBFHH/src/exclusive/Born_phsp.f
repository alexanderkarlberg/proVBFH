      subroutine born_phsp(xborn)
      implicit none
      include 'brinclude.h'
      include 'pwhg_kn.h'
      include 'coupl.inc'
      include 'PhysPars.h'
      real * 8 xborn(ndiminteg-3),jac,tmpvec(0:3),bwcutoff,mass
      external mass
      integer k, mu, j, i
      logical ini,fullphsp,check
      data check/.false./     
      data ini/.true./
      save ini,fullphsp
      real * 8 powheginput
      logical higgsfixedwidth
      real * 8 pt1cut,pt2cut,pt1,pt2,m2jjmin,m2jj,pt3,pt3cut
      real * 8 mjj1, mjj2, mjj3, dotp, m2, m2min
      real * 8 mjj_max, mjj_maxcut
      save pt1cut, pt2cut, pt3cut, m2jjmin, mjj_maxcut
      external dotp
      logical higgs_use_BW
      common/useBW/higgs_use_BW

C -   N.B: To check the phase space volume is correct,                                     
C -   replace res(j) in sigborn.f by res(j)=kn_jacborn                                     
C -   then check that the cross section is equal to,                                       
C -   hc^2/(8*pi) = 1.54929*10^7.            

      if(ini) then
c     set initial- and final-state masses for Born and real
         do k=1,nlegborn
            kn_masses(k)=0d0
         enddo
         kn_masses(nlegreal)=0d0
         kn_masses(3)=hmass
         kn_masses(4)=hmass

c     set initial- and final-state masses for Born and real
         higgsfixedwidth=powheginput("#higgsfixedwidth").gt.0
         ini=.false.


c     determine whether the Higgs boson propagator is
C     narrow width or breit wigner
         higgs_use_BW = .false.
         if (powheginput("#higgsbreitwigner").eq.1) higgs_use_BW=.true.
         
      endif      

      brkn_ktmin=0
      kn_ktmin=0

      !VBFHHMOD: Changed to hhjj phase space
      call born_phsp_hhjj(xborn)
      if(brkn_jacborn.eq.0d0) then
         kn_jacborn = 0d0
         return
      endif
      
      brkn_emitter=0
      call br_real_phsp_isr_new(xborn(ndiminteg-5),jac)
      kn_cmpborn=brkn_cmpreal
      kn_pborn=brkn_preal
      kn_xb1=brkn_x1
      kn_xb2=brkn_x2
      jac=jac

      kn_jacborn=brkn_jacborn*jac

c     set the CMS energy 
      kn_sborn=brkn_sreal

c     minimal final state mass 
      kn_minmass=kn_ktmin + sqrt(kn_ktmin**2 + ph_Hmass2low)
      
      if (check) then
         write(*,*) ''
         print*, 'kn_jacborn', kn_jacborn
         write(*,*) 'new set CM _ALL'  
         write(*,*) 'random no', xborn(ndiminteg-3)
         do j=1,nlegborn
            write(*,*) 'mom ',j,(kn_cmpborn(mu,j),mu=0,3)
         enddo      
         do j=1,nlegborn
            write(*,*) 'mom ',j,(kn_pborn(mu,j),mu=0,3)
            write(*,*) 'mass ',j,mass(kn_pborn(0,j))
         enddo
         call checkmomzero(nlegborn,kn_pborn)
      endif

      end


      subroutine born_phsp_hhjj(xborn)
      implicit none
      include 'brinclude.h'
      include 'pwhg_kn.h'
      include 'pwhg_par.h'
      include 'pwhg_math.h'
      include 'PhysPars.h'
      include 'PhysPars_Higgs.h'
      include 'incl2pwhg.h'
      real * 8 xborn(ndiminteg-6)
      real * 8 xx(ndiminteg-6)
      real * 8 xjac,taumin,tau,y,beta,betaCM,vec(3),cth,s,
     #     z,zhigh,zlow,n,m2,m2H1,m2H2
      integer mu,k,j
      parameter (n = 2d0)
      logical ini
      data ini/.true./
      save ini
      real * 8 Vmass2,Vmass2low,Vmass2high,VmVw  
      real * 8 m2jj,pV(0:3),pVmod,pVmod2,pJ(0:3,2),cthj,phij,pcmjj(0:3),
     #     pcmjjmod,ptmp(0:3,2), weight, m22,pt3cut
      real * 8 mass
      real * 8 pH1(0:3), pH2(0:3)
      real * 8 pH13m(1:4), pH23m(1:4), pV3m(1:4)
      external mass      
      logical check
      parameter(check=.false.)
      logical higgs_use_BW
      common/useBW/higgs_use_BW
c      parameter (BW=.true.)
      real * 8 epsilon
      parameter (epsilon=1d-10)
      real * 8 pt1cut,pt2cut,pt1,pt2,m2jjmin
      real * 8 BW_fixed,BW_running
      logical higgsfixedwidth
      save higgsfixedwidth
      real * 8 powheginput
      external powheginput
!     VBFNLO variables
      DOUBLE PRECISION RM2(0:2),RMG(0:2),RM2MIN(0:2),RM2MAX(0:2),ecm,pTjmin
      save RM2,RMG,RM2MIN,RM2MAX,ecm, pTjmin
      double precision q(0:4), k1(0:3), k2(0:3)
      logical Resonance, TwoToJetsPlusX, TwoBodyDecay
      external Resonance, TwoToJetsPlusX, TwoBodyDecay

      xx = xborn
!     Start by initialising everything we need
      brkn_masses = 0d0
      brkn_sborn = 0d0 ! Isn't needed for the inclusive part of the code..
      brkn_xb1 = 1d0
      brkn_xb2 = 1d0
      brkn_pborn = 0d0
      brkn_jacborn = 1d0
      q = 0d0
      if(ini) then
         rm2(0) = (2d0*ph_Hmass)**2 ! Pseudo resonance with 2*Higgs mass
         rmg(0) = ph_Hmass*500d0 ! Pseudo resonance width times mass (broad resonance)
         rm2min(0) = 1d-3       ! This value taken from vbfnlo. Should test
         rm2max(0) = kn_sbeams*0.5d0 ! Maximum allowed value for intermediate mass
         ecm = sqrt(kn_sbeams)  ! Centre-of-mass energy
         pTjmin = 0d0 ! minimum jet pt

         ! Currently we assume on-shell Higgs. Should be changed later for breit-wigner perhaps?
         rm2(1) = ph_Hmass**2 ! Mass squared of first Higgs
         rm2(2) = ph_Hmass**2 ! Mass squared of first Higgs
         ini = .false.
      endif
      
!     We borrow the phase space for HHjj from vbfnlo.  First generate
!     q^2 of intermediate pseudo-state.

!     AK: Understanding here is that we generate a fake resonance with a
!     large width to perform pseudo-decay into two on-shell (or
!     off-shell) Higgs. The routine returns brkn_jacborn (0d0 if the
!     function is false) and the q^2 as q(4).
      brkn_jacborn = 3d0*xx(1)**(3d0-1d0)*brkn_jacborn
      xx(1) = 1d0 - xx(1)**3
      if (.not. Resonance(rm2(0), rmg(0), rm2min(0), rm2max(0),
     $     xx(1), brkn_jacborn, q(4))) return !continue
      
      brkn_jacborn = 2d0*xx(2)**(2d0-1d0)*brkn_jacborn
      xx(2) = 1d0 - xx(2)**2
      if(xx(4).lt.0.5d0) then
         xx(4) = 2d0*xx(4)
         xx(4) = xx(4)**n
         brkn_jacborn =  brkn_jacborn * n * xx(4)**(1d0-1d0/n)
         xx(4) = (1d0 - xx(4))*0.5d0
      else
         xx(4) = 2d0*(1d0 - xx(4))
         xx(4) = xx(4)**n
         brkn_jacborn =  brkn_jacborn * n * xx(4)**(1d0-1d0/n)
         xx(4) = (1d0 + xx(4))*0.5d0
      endif
      brkn_jacborn = 2d0*xx(5)**(2d0-1d0)*brkn_jacborn
      xx(5) = xx(5)**2
!     Now we generate the two-jet system from the pseudo-resonance q^2
!      if (.not.TwoToJetsPlusX(2, xx(2:7), 0d0, ecm, pTjmin, q(4),
      if (.not.TwoToJetsPlusX(2, xx(2:7), xx(10), ecm, pTjmin, q(4),
     $     brkn_pborn(:,1), brkn_pborn(:,2), brkn_xb1, brkn_xb2, q(0),
     $     brkn_pborn(:,5:6), brkn_jacborn)) return! continue

      brkn_sborn = brkn_xb1*brkn_xb2*kn_sbeams ! AK I think this is right...
!      print*, brkn_xb1, brkn_xb2, kn_sbeams, brkn_sborn,brkn_jacborn
!      stop

!     Now decay the pseudo-resonance
      if(xx(8).lt.0.5d0) then
         xx(8) = 2d0*xx(8)
         xx(8) = xx(8)**n
         brkn_jacborn =  brkn_jacborn * n * xx(8)**(1d0-1d0/n)
         xx(8) = 0.5d0*xx(8)
      else
         xx(8) = 2d0*(1d0 - xx(8))
         xx(8) = xx(8)**n
         brkn_jacborn =  brkn_jacborn * n * xx(8)**(1d0-1d0/n)
         xx(8) = (2d0 - xx(8))*0.5d0
      endif
      if (.not. TwoBodyDecay(xx(8), xx(9), q(0), q(4), rm2(1),
     $     rm2(2), brkn_pborn(:,3), brkn_pborn(:,4), brkn_jacborn))
     $     return !continue


      brkn_masses(3:4) = sqrt(rm2(1:2))

      brkn_jacborn = brkn_jacborn * 0.5d0 ! symmetry factor from two identical particles in final state
!      print*, 'Born momenta'
!      print*, '1', brkn_pborn(:,1), brkn_masses(1), invmass(brkn_pborn(:,1))
!      print*, '2', brkn_pborn(:,2), brkn_masses(2), invmass(brkn_pborn(:,2))
!      print*, '3', brkn_pborn(:,3), brkn_masses(3), invmass(brkn_pborn(:,3))
!      print*, '4', brkn_pborn(:,4), brkn_masses(4), invmass(brkn_pborn(:,4))
!      print*, '5', brkn_pborn(:,5), brkn_masses(5), invmass(brkn_pborn(:,5))
!      print*, '6', brkn_pborn(:,6), brkn_masses(6), invmass(brkn_pborn(:,6))
!     stop
      
c     boost in the CM frame
      betaCM=(brkn_xb1-brkn_xb2)/(brkn_xb1+brkn_xb2)
      vec(1)=0
      vec(2)=0
      vec(3)=1
      call mboost(br_nlegborn,vec,-betaCM,brkn_pborn(0,1),brkn_cmpborn(0,1))      
      if (check) then
         write(*,*) ''
         write(*,*) 'new set'     
         do j=1,nlegborn
            write(*,*) 'mom ',j,(brkn_cmpborn(mu,j),mu=0,3)
         enddo
      endif

      call br_compute_csimax_fsr
      
      end

c Mappings of the underlying born configuration in
c brkn_cmpborn(0:3,br_nlegborn), and the xrad(1:3) variables
c in the unit cube, into brkn_real(0:3,nlegreal).
c The factor jac_over_csi*csi*brkn_csimax, multiplied
c by the Born phase space jacobian, yields the real phase
c space jacobian.
c More explicitly:
c d Phi_n = d^3 xrad jac_over_csi csi csimax d Phi_{n-1}
c Since
c  d Phi_n = d phi d y d csi Jrad d Phi_{n-1}
c (where Jrad is given in FNO2006) we get
c                                  d phi d y d csi
c csimax csi jac_over_csi = Jrad  ----------------
c                                    d^3 xrad
c Notice that using d csi=d csitilde csimax the csimax
c factor cancels, and jac_over_csi is as given in the
c code below (see notes on xscaled.tm).
c br_real_phsp_fsr: provides the mapping for the final state
c radiation, assuming that the emitter is the brkn_emitter-th
c particle, and the emitted particle is the nlegreal-th particle
c br_real_phsp_isr: mapping for the initial state radiation
      subroutine br_real_phsp_fsr(xrad,jac)
      implicit none
      real * 8 xrad(3),jac
      include 'pwhg_math.h'
      include 'brinclude.h'
      real * 8 q0,q2,xjac
c Boost the underlying Born variables to their cm frame
      q0=2*brkn_cmpborn(0,1)
      q2=brkn_sborn
      brkn_csitilde=xrad(1)
      xjac=1
      brkn_y=1-2*xrad(2)
      xjac=xjac*2
c importance sampling for brkn_y
      xjac=xjac*1.5d0*(1-brkn_y**2)
      brkn_y=1.5d0*(brkn_y-brkn_y**3/3)
      brkn_azi=2*pi*xrad(3)
      xjac=xjac*2*pi
      brkn_csimax=brkn_csimax_arr(brkn_emitter)
      brkn_csi=brkn_csitilde*brkn_csimax     
c remember: no csimax in the jacobian factor, we are integrating in csitilde
      call br_real_phsp_fsr_rad
      jac=xjac*brkn_jacreal*brkn_csimax
      end

      subroutine br_real_phsp_fsr_new(xrad,jac)
      implicit none
      real * 8 xrad(3),jac, one_minus_y_cut
      real * 8 n,nn
      parameter (n = 8d0)
      parameter (nn = 2d0)
      include 'pwhg_math.h'
      include 'pwhg_par.h'
      include 'brinclude.h'
      real * 8 q0,q2,xjac
c Boost the underlying Born variables to their cm frame
!      print*, 'xrad', xrad
      q0=2*brkn_cmpborn(0,1)
      q2=brkn_sborn

      brkn_csitilde=xrad(1)
      xjac=1
c      brkn_csitilde = xrad(1)**nn
c      xjac = nn * brkn_csitilde**(1d0-1d0/nn)

c      if(xrad(1).lt.0.5d0) then
c         brkn_csitilde = 2d0*xrad(1)
c         brkn_csitilde = brkn_csitilde**nn
c         xjac = nn * brkn_csitilde**(1d0-1d0/nn)
c         brkn_csitilde = 0.5d0 * brkn_csitilde
c      else
c         brkn_csitilde = 2d0*(1d0 - xrad(1))
c         brkn_csitilde = brkn_csitilde**nn
c         xjac = nn * brkn_csitilde**(1d0-1d0/nn)
c         brkn_csitilde = 1d0 - 0.5d0 * brkn_csitilde
c      endif

      if(xrad(2).lt.0.5d0) then
         brkn_y = 2d0*xrad(2)
         brkn_y = brkn_y**n
         xjac = 2d0 * xjac * n * brkn_y**(1d0-1d0/n)
         brkn_y = 1d0 - brkn_y
      else
         brkn_y = 2d0*(1d0 - xrad(2))
         brkn_y = brkn_y**n
         xjac = 2d0 * xjac * n * brkn_y**(1d0-1d0/n)
         brkn_y = - 1d0 + brkn_y
      endif
      brkn_csitilde = brkn_csitilde * (1 - 2 * par_isrtinycsi) + par_isrtinycsi
      brkn_y = brkn_y * (1 - par_fsrtinyy)
c      one_minus_y_cut = 1d-4 ! Any higher than 1d-10 and things start crashing...
c      if(abs(brkn_y).gt.1d0-one_minus_y_cut) brkn_y = sign(brkn_y,1d0) *
c     $     (1d0 - one_minus_y_cut)

!      brkn_azi=2*pi*xrad(3)
!      xjac=xjac*2*pi

      if(xrad(3).gt.0.5d0) then
         brkn_azi=(xrad(3)-0.5d0)*2d0
         brkn_azi = pi*brkn_azi**nn 
         xjac = 2d0*xjac*nn*pi*brkn_azi**(1d0-1d0/nn)
         brkn_azi = brkn_azi + pi
      else
         brkn_azi=(0.5d0 - xrad(3))*2d0
         brkn_azi = pi*brkn_azi**nn 
         xjac = 2d0*xjac*nn*pi*brkn_azi**(1d0-1d0/nn)
         brkn_azi = pi - brkn_azi 
      endif

      brkn_csimax=brkn_csimax_arr(brkn_emitter)
      brkn_csi=brkn_csitilde*brkn_csimax     
c remember: no csimax in the jacobian factor, we are integrating in csitilde 
!      print*, brkn_csitilde, brkn_csimax, xjac, brkn_azi
      call br_real_phsp_fsr_rad_new
      jac=xjac*brkn_jacreal*brkn_csimax
      end

c br_real_phsp_fsr_rad: provides the mapping for the final state
c radiation, assuming that we are considering the region rad_kinreg
c and the emitted particle is the nlegreal-th particle,
c for given brkn_csi, brkn_y, brkn_azi. Sets the jacobian
c brkn_jacreal so that brkn_jacreal d brkn_csi d brkn_y d brkn_azi times
c the underlying Born jacobian is the phase space volume
      subroutine br_real_phsp_fsr_rad
      implicit none
      include 'pwhg_math.h'
      include 'brinclude.h'
      real * 8 vec(3),q0,beta
      integer i
      data vec/0d0,0d0,1d0/
      save vec
      q0=2*brkn_cmpborn(0,1)
c remember: no csimax factor, we are integrating in csitilde 
c      print*, 'brkn_csi,brkn_y,brkn_azi', brkn_csi,brkn_y,brkn_azi
      call barradmap(br_nlegborn-2,brkn_emitter-2,q0,brkn_cmpborn(0,3),
     1    brkn_csi,brkn_y,brkn_azi,brkn_preal(0,3),brkn_jacreal)
      beta=(brkn_xb1-brkn_xb2)/(brkn_xb1+brkn_xb2)
      call mboost(br_nlegreal-2,vec,beta,
     1    brkn_preal(0,3),brkn_preal(0,3))
      do i=0,3
         brkn_preal(i,1)=brkn_pborn(i,1)
         brkn_preal(i,2)=brkn_pborn(i,2)
      enddo
      brkn_x1=brkn_xb1
      brkn_x2=brkn_xb2
      brkn_sreal=brkn_sborn
c      call checkmomzero(br_nlegreal,brkn_preal)
      call br_compcmkin
      call br_compdij
      end

c br_real_phsp_fsr_rad: provides the mapping for the final state
c radiation, assuming that we are considering the region rad_kinreg
c and the emitted particle is the nlegreal-th particle,
c for given brkn_csi, brkn_y, brkn_azi. Sets the jacobian
c brkn_jacreal so that brkn_jacreal d brkn_csi d brkn_y d brkn_azi times
c the underlying Born jacobian is the phase space volume
      subroutine br_real_phsp_fsr_rad_new
      implicit none
      include 'pwhg_math.h'
      include 'brinclude.h'
      real * 8 vec(3),q0,beta
      integer i
      data vec/0d0,0d0,1d0/
      save vec
      q0=2*brkn_cmpborn(0,1)
c remember: no csimax factor, we are integrating in csitilde 
c      print*, 'brkn_csi,brkn_y,brkn_azi', brkn_csi,brkn_y,brkn_azi
      call barradmap_new(br_nlegborn-2,brkn_emitter-2,q0,brkn_cmpborn(0,3),
     1    brkn_csi,brkn_y,brkn_azi,brkn_preal(0,3),brkn_jacreal)
      beta=(brkn_xb1-brkn_xb2)/(brkn_xb1+brkn_xb2)
      call mboost(br_nlegreal-2,vec,beta,
     1    brkn_preal(0,3),brkn_preal(0,3))
      do i=0,3
         brkn_preal(i,1)=brkn_pborn(i,1)
         brkn_preal(i,2)=brkn_pborn(i,2)
      enddo
      brkn_x1=brkn_xb1
      brkn_x2=brkn_xb2
      brkn_sreal=brkn_sborn
c      call checkmomzero(br_nlegreal,brkn_preal)
      call br_compcmkin
      call br_compdij
      end

      subroutine br_real_phsp_isr(xrad,jac)
      implicit none
      real * 8 xrad(3),jac
      include 'pwhg_math.h'
      include 'brinclude.h'
      real * 8 xjac
      brkn_csitilde=(3-2*xrad(1))*xrad(1)**2
      xjac=6*(1-xrad(1))*xrad(1)
c      print*, 'brkn_csitilde', brkn_csitilde
c      print*, 'xjac', xjac
      brkn_y=1-2*xrad(2)
      xjac=xjac*2
      xjac=xjac*1.5d0*(1-brkn_y**2)
      brkn_y=1.5d0*(brkn_y-brkn_y**3/3)
      brkn_azi=2*pi*xrad(3)
      xjac=xjac*2*pi
      call br_compcsimax
      brkn_csi=brkn_csitilde*brkn_csimax
      call br_real_phsp_isr_rad
c      print*, 'brkn_y', brkn_y
c      print*, 'brkn_jacreal', brkn_jacreal
c      print*, 'brkn_csimax', brkn_csimax
c      print*, 'brkn_csitilde', brkn_csitilde
c      print*, 'brkn_csi', brkn_csi
c      print*, 'xjac real', xjac
      jac=xjac*brkn_jacreal*brkn_csimax
      end
      
      subroutine br_real_phsp_isr_new(xrad,jac)
      implicit none
      real * 8 xrad(3),jac, one_minus_y_cut
      real * 8 n, nn
      parameter (n = 8d0)
      parameter (nn = 6d0)
      include 'pwhg_math.h'
      include 'brinclude.h'
      include 'pwhg_par.h'
      real * 8 xjac

      brkn_csitilde=(3-2*xrad(1))*xrad(1)**nn
      xjac=(3*nn-2*(nn+1)*xrad(1))*xrad(1)**(nn-1)

!      brkn_csitilde=(3-2*xrad(1))*xrad(1)**2
!      xjac=6*(1-xrad(1))*xrad(1)
      if(xrad(2).lt.0.5d0) then
         brkn_y = 2d0*xrad(2)
         brkn_y = brkn_y**n
         xjac = 2d0 * xjac * n * brkn_y**(1d0-1d0/n)
         brkn_y = 1d0 - brkn_y
      else
         brkn_y = 2d0*(1d0 - xrad(2))
         brkn_y = brkn_y**n
         xjac = 2d0 * xjac * n * brkn_y**(1d0-1d0/n)
         brkn_y = - 1d0 + brkn_y
      endif
      brkn_y = brkn_y * (1 - par_isrtinyy)
      brkn_csitilde = brkn_csitilde * (1 - 2 * par_isrtinycsi) + par_isrtinycsi

      brkn_azi=2*pi*xrad(3)
      xjac=xjac*2*pi
      call br_compcsimax
      brkn_csi=brkn_csitilde*brkn_csimax
      call br_real_phsp_isr_rad
c      print*, 'brkn_y', brkn_y
c      print*, 'brkn_jacreal', brkn_jacreal
c      print*, 'brkn_csimax', brkn_csimax
c      print*, 'brkn_csitilde', brkn_csitilde
c      print*, 'brkn_csi', brkn_csi
c      print*, 'xjac real', xjac
c      stop
      jac=xjac*brkn_jacreal*brkn_csimax
      end
      
      subroutine br_compcsimax
      implicit none
      include 'brinclude.h'
      real * 8 y,xb1,xb2
      xb1=brkn_xb1
      xb2=brkn_xb2
      y=brkn_y
      brkn_csimax=1-max(2*(1+y)*xb1**2/
     1    (sqrt((1+xb1**2)**2*(1-y)**2+16*y*xb1**2)+(1-y)*(1-xb1**2)),
     1            2*(1-y)*xb2**2/
     1    (sqrt((1+xb2**2)**2*(1+y)**2-16*y*xb2**2)+(1+y)*(1-xb2**2)))
      end

      subroutine br_real_phsp_isr_rad
      implicit none
      include 'pwhg_math.h'
      include 'brinclude.h'
      include 'pwhg_kn.h'
      real * 8 y,xb1,xb2,x1,x2,betal,betat,vecl(3),vect(3),
     1         cth,sth,cph,sph,csi,pt2
      integer i,mu
      real * 8 dotp
      external dotp
c the following call sets brkn_csimax, brkn_csimaxp, brkn_csimaxm
c also when br_real_phsp_isr_rad is called directly
c (i.e. not through br_real_phsp_isr_rad0)
      call br_compcsimax
      y=brkn_y
      xb1=brkn_xb1
      xb2=brkn_xb2
      csi=brkn_csi
      cth=y
      sth=sqrt(1-cth**2)
      cph=cos(brkn_azi)
      sph=sin(brkn_azi)
      x1=xb1/sqrt(1-csi)*sqrt((2-csi*(1-y))/(2-csi*(1+y)))
      x2=xb2/sqrt(1-csi)*sqrt((2-csi*(1+y))/(2-csi*(1-y)))
      brkn_x1=x1
      brkn_x2=x2
      do mu=0,3
         brkn_preal(mu,1)=kn_beams(mu,1)*x1
         brkn_preal(mu,2)=kn_beams(mu,2)*x2
      enddo
      brkn_sreal=brkn_sborn/(1-csi)
c Build k_n+1 in the rest frame of brkn_preal
c      write(*,*) ' br_nlegreal ',br_nlegreal
      brkn_preal(0,br_nlegreal)=sqrt(brkn_sreal)*csi/2
      brkn_preal(1,br_nlegreal)=brkn_preal(0,br_nlegreal)*sth*sph
      brkn_preal(2,br_nlegreal)=brkn_preal(0,br_nlegreal)*sth*cph
      brkn_preal(3,br_nlegreal)=brkn_preal(0,br_nlegreal)*cth
c boost it to the frame of brkn_preal
      do i=1,3
         vecl(i)=(brkn_preal(i,1)+brkn_preal(i,2))
     1          /(brkn_preal(0,1)+brkn_preal(0,2))
      enddo      
      betal=sqrt(vecl(1)**2+vecl(2)**2+vecl(3)**2)
      if(betal.gt.0) then
         do i=1,3
            vecl(i)=vecl(i)/betal
         enddo
      else
         vecl(1)=1
         vecl(2)=0
         vecl(3)=0
      endif
      call mboost(1,vecl,betal,
     1    brkn_preal(0,br_nlegreal),brkn_preal(0,br_nlegreal))
c longitudinal boost of underlying Born to zero rapidity frame
      do i=1,3
         vecl(i)=(brkn_pborn(i,1)+brkn_pborn(i,2))
     1          /(brkn_pborn(0,1)+brkn_pborn(0,2))
      enddo
      betal=sqrt(vecl(1)**2+vecl(2)**2+vecl(3)**2)
      if(betal.gt.0) then
         do i=1,3
            vecl(i)=vecl(i)/betal
         enddo
      else
         vecl(1)=1
         vecl(2)=0
         vecl(3)=0
      endif
      call mboost(br_nlegborn-2,vecl,-betal,
     1 brkn_pborn(0,3),brkn_preal(0,3))
c      call printtot(br_nlegborn,brkn_preal(0,1))
c construct transverse boost velocity
      vect(3)=0
      vect(1)=brkn_preal(1,br_nlegreal)
      vect(2)=brkn_preal(2,br_nlegreal)
      pt2=vect(1)**2+vect(2)**2
c      betat=1/sqrt(1+(brkn_sreal*(1-csi))/pt2)
      betat=sqrt(pt2/(pt2+brkn_sreal*(1-csi)))
      if(pt2.eq.0) then
         vect(1)=1
         vect(2)=0
      else
         vect(1)=vect(1)/sqrt(pt2)
         vect(2)=vect(2)/sqrt(pt2)
      endif
c     write(*,*) ' k+1: ',(brkn_preal(mu,br_nlegreal),mu=0,3)
      call mboost(br_nlegborn-2,vect,-betat,
     1     brkn_preal(0,3),brkn_preal(0,3))
c     call printtot(nlegborn,brkn_preal(0,1))
c     longitudinal boost in opposite direction
      call mboost(br_nlegborn-2,vecl,betal,
     1     brkn_preal(0,3),brkn_preal(0,3))
c     call printtot(br_nlegreal,brkn_preal(0,1))
c      print*, 'brkn_sreal', brkn_sreal
c      print*, 'csi', csi
c      print*, 'csi/(1-csi)',csi/(1-csi) 
      brkn_jacreal=brkn_sreal/(4*pi)**3*csi/(1-csi)
      call br_compcmkin
      call br_compdij
      end


      subroutine br_compcmkin
      implicit none
      include 'brinclude.h'
      real * 8 vecl(3),betal
      data vecl/0d0,0d0,1d0/
      save vecl
      betal=-(brkn_preal(3,1)+brkn_preal(3,2))/
     1 (brkn_preal(0,1)+brkn_preal(0,2))
      call mboost(br_nlegreal,vecl,betal,brkn_preal,brkn_cmpreal)
      end

      subroutine br_compdij
      implicit none
      include 'brinclude.h'
      integer j,k
      real * 8 y
      real * 8 crossp,dotp
      external crossp,dotp
!     VBFHHMOD
!     do j=4,br_nlegreal
      do j=5,br_nlegreal
         y=1-dotp(brkn_cmpreal(0,1),brkn_cmpreal(0,j))
     1 /(brkn_cmpreal(0,1)*brkn_cmpreal(0,j))
         brkn_dijterm(0,j)=(brkn_cmpreal(0,j)**2
     1 *(1-y**2))**brpar_diexp
         brkn_dijterm(1,j)=(brkn_cmpreal(0,j)**2
     1 *2*(1-y))**brpar_diexp
         brkn_dijterm(2,j)=(brkn_cmpreal(0,j)**2
     1 *2*(1+y))**brpar_diexp
      enddo
!      do j=4,br_nlegreal
!     VBFHHMOD
      do j=5,br_nlegreal
         do k=j+1,br_nlegreal
            brkn_dijterm(j,k)=
     1(2*dotp(brkn_cmpreal(0,k),brkn_cmpreal(0,j))*
     1       brkn_cmpreal(0,k)*brkn_cmpreal(0,j)
     2    /  (brkn_cmpreal(0,k)+brkn_cmpreal(0,j))**2)**brpar_dijexp
c     2    /  ((brkn_cmpreal(1,k)+brkn_cmpreal(1,j))**2+
c     3        (brkn_cmpreal(2,k)+brkn_cmpreal(2,j))**2+
c     4        (brkn_cmpreal(3,k)+brkn_cmpreal(3,j))**2))**brpar_dijexp
         enddo
      enddo
      end

      subroutine br_compute_csimax_fsr
      implicit none
c Compute csimax for all possible final state emitters;
c for initial state emitters it is not possible, since
c csimax depends upon y in this case.
      include 'brinclude.h'
      integer j
      real * 8 q0,mrec2
      logical valid_emitter
      external valid_emitter
!       
!      do j=4,5
!     VBFHHMOD
      do j=5,6
          q0=2*brkn_cmpborn(0,1)
          mrec2=(q0-brkn_cmpborn(0,j))**2
     1           -brkn_cmpborn(1,j)**2-brkn_cmpborn(2,j)**2
     1           -brkn_cmpborn(3,j)**2
          brkn_csimax_arr(j)=1-mrec2/q0**2
      enddo
      end




      function mass(p)
      implicit none
      real * 8 p(0:3),mass
      mass = sqrt(abs(p(0)**2-p(1)**2-p(2)**2-p(3)**2))
      end



      subroutine born_suppression(fact)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_flg.h'

      integer i,j
      real*8 pt1,pt2,pt3,Ht
      real * 8 fact,ptmin, mjjmin, mj1j2min
      real * 8 dotp,pt4sq,pt5sq,pt6sq,pt7sq
      real * 8 kp, km, mjjminsq,ptminsq
      parameter (kp = 2d0)
      parameter (km = 2d0)
      save mjjminsq,ptminsq
      real * 8 m45sq, m46sq, m56sq
      logical, save :: ini = .true. 
      integer, save :: btype
      real*8 powheginput
      external powheginput

      real*8 ptsqrel(nlegborn,nlegborn)
      real*8 fact1,fact2
      real*8 p,scale1,scale2
      real*8 pt1sq,pt2sq,pt3sq
      real*8 msq12,msq23,msq31

      if (ini) then 
         if (flg_weightedev) then 
            btype = powheginput('#bornsuppfact')
            ptmin=powheginput('#ptsuppfact')
            mjjmin=powheginput('#mjjsuppfact')
            if (btype.eq.1) then 
               write(*,*) 'Using exponential Born suppression factor' 
               write(*,*) 'with scale1[GeV]=',ptmin ,"^(2)"
               write(*,*) 'with scale2[GeV]=',mjjmin ,"^(2)"
            elseif (btype.eq.2) then 
               write(*,*) 'Using multiplicative Born suppression factor' 
               write(*,*) 'with scalep[GeV]=',ptmin ,"^(2)"
               write(*,*) 'with scalem[GeV]=',mjjmin ,"^(2)"
               ptminsq  = ptmin**2
               mjjminsq = mjjmin**2
            else
               write(*,*) 'this value of bornsuppfact is not supported'
               write(*,*) 'please change your input'
               stop
            endif   
         else
            write(*,*) 'Using no Born suppression' 
         endif
         ini = .false. 
      endif


      if(flg_weightedev) then
         if (btype.eq.1) then 
c     expon. damping (c.f. trijets, arXiv:1402.4001):
            scale1 = ptmin**2
            scale2 = mjjmin**2
            p = 2d0
            
            pt1sq = (kn_cmpborn(1,4)**2+kn_cmpborn(2,4)**2)
            pt2sq = (kn_cmpborn(1,5)**2+kn_cmpborn(2,5)**2)
            pt3sq = (kn_cmpborn(1,6)**2+kn_cmpborn(2,6)**2)

            do i = 1, nlegborn
               do j = 1, nlegborn 
                  ptsqrel(i,j)=2*dotp(kn_cmpborn(0,i),kn_cmpborn(0,j))
     1                 *  kn_cmpborn(0,i)*kn_cmpborn(0,j)
     2                 / (kn_cmpborn(0,i)**2+kn_cmpborn(0,j)**2)
               enddo
            enddo

            msq12=ptsqrel(4,5)
            msq31=ptsqrel(4,6)
            msq23=ptsqrel(5,6)
            
            fact1 = exp(-scale1**p * (
     1           1/pt1sq**p
     2           +1/pt2sq**p
     3           +1/pt3sq**p
     4           +1/msq12**p
     5           +1/msq23**p
     6           +1/msq31**p ))
            
            pt1 = sqrt(pt1sq)
            pt2 = sqrt(pt2sq)
            pt3 = sqrt(pt3sq)
            Ht  = pt1 + pt2 + pt3
            
            fact = fact1/scale2**p/(1d0/scale2 + 1d0/Ht**2)**p
            
         elseif (btype.eq.2) then 
            
            m45sq=2d0*dotp(kn_cmpborn(0,4),kn_cmpborn(0,5))
            m46sq=2d0*dotp(kn_cmpborn(0,4),kn_cmpborn(0,6))
            m56sq=2d0*dotp(kn_cmpborn(0,5),kn_cmpborn(0,6))
            
            pt4sq=(kn_cmpborn(1,4)**2+kn_cmpborn(2,4)**2)
            pt5sq=(kn_cmpborn(1,5)**2+kn_cmpborn(2,5)**2)
            pt6sq=(kn_cmpborn(1,6)**2+kn_cmpborn(2,6)**2)
c            pt7sq=(kn_cmpreal(1,7)**2+kn_cmpreal(2,7)**2)
            
            fact=
     &           (pt4sq/(pt4sq+ptminsq))**kp * 
     &           (pt5sq/(pt5sq+ptminsq))**kp *  
     &           (pt6sq/(pt6sq+ptminsq))**kp * 
c     &           (pt7sq/(pt7sq+ptminsq))**kp * 
     &           (m45sq/(m45sq+mjjminsq))**km * 
     &           (m46sq/(m46sq+mjjminsq))**km * 
     &           (m56sq/(m56sq+mjjminsq))**km
            
         endif                  !btype
         
      else                      ! no Born suppression
         
         fact=1
         
      endif
c      print*, 'fact', fact
      end


      subroutine set_fac_ren_scales(muf,mur)
      implicit none
      include 'PhysPars.h'
      include 'PhysPars_Higgs.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_st.h'
      include 'pwhg_flg.h'
      real*8 pwhg_alphas
      external pwhg_alphas
      real * 8 powheginput
      external powheginput
      real * 8 muf,mur, muref,renscfact,facscfact
      save renscfact,facscfact
      logical runningscales
      real*8 ptj1, ptj2, ptj3, ptj4, eth, pth
      
      logical ini
      data ini/.true./
      save ini

      if (ini) then
         if(powheginput('#runningscales').ne.0d0) then
            runningscales =.true.
            if (powheginput('#runningscales').eq.2d0) 
     $           stop 'runningscales = 2 unavailable for exclusive runs'
         else
            runningscales = .false.
         endif
         renscfact=powheginput("#renscfact")
         facscfact=powheginput("#facscfact")
      endif
      
      if (renscfact .eq. 0d0) stop 'renscale = 0 not allowed'
      if (facscfact .eq. 0d0) stop 'facscale = 0 not allowed'

      if(runningscales) then
         if (ini) then
            write(*,*) '*************************************'
            write(*,*) 'Factorization and renormalization '
            write(*,*) 'scales (mur, muf) set to '
            write(*,*) '((MH/2)**4+(MH*ptH/2)**2)**1/4'
            if (renscfact .gt. 0d0) 
     &           write(*,*) 'Renormalization scale rescaled by', renscfact
            if (facscfact .gt. 0d0) 
     &           write(*,*) 'Factorization scale rescaled by  ', facscfact
            write(*,*) '*************************************'
            ini=.false.
         endif
         
c     default is Born kinematics:
         ptH = sqrt((kn_pborn(1,3)+kn_pborn(1,4))**2
     $        +(kn_pborn(2,3)+kn_pborn(2,4))**2) ! ptHH
         
         
         muref = ( ( ph_Hmass * 0.5d0 ) ** 4d0
     $        + ( ph_Hmass * ptH * 0.5d0 ) ** 2d0 
     $        ) ** 0.25d0
         
      else
         muref=ph_Hmass
         if (ini) then
            st_mufact2= muref**2*st_facfact**2
            st_muren2 = muref**2*st_renfact**2
            st_alpha  = pwhg_alphas(st_muren2,st_lambda5MSB,st_nlight)
c     print*, st_alpha,st_muren2,st_lambda5MSB,st_nlight
c     stop
            write(*,*) '**********************************'
            write(*,*) 'RENORMALIZATION SCALE = ',sqrt(st_muren2)
            write(*,*) 'FACTORIZATION   SCALE = ',sqrt(st_mufact2)
            write(*,*) 'alfa_s =        ', st_alpha !pwhg_alphas(muf**2,st_lambda5MSB,st_nlight)
            write(*,*) '**********************************'
            ini=.false.
         endif      
      endif

      muf=muref
      mur=muref
      
c     make sure muf never falls below min. cutoff value:
      muf = max(muf,dsqrt(2d0))       
c     make sure mur never falls below min. cutoff value:
      mur = max(mur,dsqrt(2d0)) 
      
      end

c This routine performs the inverse mapping from barred and radiation
c variables to the n+1 momenta, as in Sec. 5.2.1 in fno2006.
c All particle can have masses, except for the n+1-th and j-th.
c conventions: vector(4)=(x,y,z,t)
c Input:
c n           : number of final state barred momenta
c j           : the emitter
c q0          : CM energy
c barredk(4,n): the n barred-k 4-vectors
c csi,y,phi   : the radiation variables
c Output:
c xk(4,n+1)   : the n+1 real final state momenta
c jac         : jacobian factor on phirad
      subroutine barradmap_new(n,j,q0,barredk,csi,y,phi,xk,jac)
      implicit none
c parameters
      include 'pwhg_math.h'
      integer n,j
      real * 8 q0,barredk(0:3,n),csi,y,phi,xk(0:3,n+1),jac
C Local variables
      real * 8 q2,mrec2,k0np1,uknp1,ukj,uk,cpsi,cpsi1,ubkj,vec(3),
     #     norm,k0rec,ukrec,beta,k2
      integer i
!      print*, 'Entered with', csi, y, phi
c     according to fno2006: by k0 we mean the 0 component in the CM, by
c     uk (underlined k) we mean the modulus of its 3-momentum n and np1
c     in a variable name suggests n and n+1, etc.
      q2=q0**2
c (5.42) of fnw2006
      k0np1=csi*q0/2
      uknp1=k0np1
c compute Mrec^2 (5.45)
      mrec2=(q0-barredk(0,j))**2
     #     -barredk(1,j)**2-barredk(2,j)**2-barredk(3,j)**2
      ukj=(q2-mrec2-2*q0*uknp1)/(2*(q0-uknp1*(1-y)))
c compute the length of k (5.44)
      uk=sqrt(ukj**2+uknp1**2+2*ukj*uknp1*y)
c compute cos psi (angle between knp1 and k)
      cpsi=(uk**2+uknp1**2-ukj**2)/(2*uk*uknp1)
c get the cosine of the angle between kn and k
      cpsi1=(uk**2+ukj**2-uknp1**2)/(2*uk*ukj)
c Set k_j and k_n+1 parallel to kbar_n
      ubkj=barredk(0,j)
      do i=0,3
         xk(i,j)=ukj*barredk(i,j)/ubkj
         xk(i,n+1)=uknp1*barredk(i,j)/ubkj
      enddo
c Set up a unit vector orthogonal to kbar_n and to the z axis
      vec(3)=0
      norm=sqrt(barredk(1,j)**2+barredk(2,j)**2)
      vec(1)=barredk(2,j)/norm
      vec(2)=-barredk(1,j)/norm
c Rotate k_n+1 around vec of an amount psi
      call mrotate(vec,sqrt(abs(1-cpsi**2)),cpsi,xk(1,n+1))
c      print*, 1-cpsi1**2,cpsi1, csi
c Rotate k_j around vec of an amount psi1 in opposite direction
      call mrotate(vec,-sqrt(abs(1-cpsi1**2)),cpsi1,xk(1,j))
c set up a unit vector parallel to kbar_j
      do i=1,3
         vec(i)=barredk(i,j)/ubkj
      enddo
c Rotate k_j and k_n+1 around this vector of an amount phi
      call mrotate(vec,sin(phi),cos(phi),xk(1,n+1))
      call mrotate(vec,sin(phi),cos(phi),xk(1,j))
c compute the boost velocity
      k0rec=q0-ukj-uknp1
c use abs to fix tiny negative root FPE
      ukrec=sqrt(abs(k0rec**2-mrec2))
      beta=(q2-(k0rec+ukrec)**2)/(q2+(k0rec+ukrec)**2)
c     Boost all other barred k (i.e. 1 to j-1,j+1 to n) along vec with velocity
c     beta in the k direction (same as barred k_j)
      do i=1,3
         vec(i)=barredk(i,j)/ubkj
      enddo
      call mboost(j-1,vec,beta,barredk(0,1),xk(0,1))
      if(n-j.gt.0) call mboost(n-j,vec,beta,barredk(0,j+1),xk(0,j+1))
      k2=2*ukj*uknp1*(1-y)
c returns jacobian of FNO 5.40 (i.e. jac*d csi * d y * d phi is phase space)
      jac=q2*csi/(4*pi)**3*ukj**2/ubkj/(ukj-k2/(2*q0))
      end

      subroutine swap_momenta(p1,p2)
      implicit none
      real * 8 p1(0:3),p2(0:3),tmp(0:3)
      tmp=p1
      p1=p2
      p2=tmp
      end
! All routines below this line have been taken from the phasespace of
! vbfnlo (version 3.0.0 beta 5). Thye shoudl be credited whenever using
! the code.
****************************************************************************
      LOGICAL FUNCTION Resonance(xm2, xmg, rm2min, rm2max, R, weight, qsq)
****************************************************************************
*     Calculates the Q-square for a resonance
*
*    INPUT
*       xm2       : mass**2 of the resonance
*       xmg       : mass*width of the resonance
*       rm2min    : lowerBound**2 of the resonance
*       rm2max    : upperBound**2 of the resonance
*       R         : a random number between 0 and 1
*       weight    : phasespace weight
*
*    OUTPUT
*       weight    : new phasespace weight (multiplied by Jacobi factor)
*       qsq       : the Q**2 
*       Resonance : .true.,  if everything is ok 
*                   .false., if kinematics not ok -> weight = 0
****************************************************************************
      IMPLICIT NONE

      DOUBLE PRECISION xm2, xmg, rm2min, rm2max, R, weight, qsq, tanx
      DOUBLE PRECISION twopi
      PARAMETER (twopi = 2d0*3.141592653589793d0)
      DOUBLE PRECISION xmin, xmax, deltax, x 

      IF (rm2max.LT.rm2min) THEN
         weight = 0d0 
         Resonance = .false.
         RETURN
      ENDIF
      
      xmin = ATAN( (rm2min - xm2)/xmg )
      xmax = ATAN( (rm2max - xm2)/xmg )
      deltax = xmax - xmin
      x = xmin + deltax * R
      tanx = tan(x)
      qsq = xm2 + tanx * xmg

      
      IF (qsq.LE.0 .or. deltax.le.0) THEN
         weight = 0d0 
         Resonance = .false.
         RETURN 
      ENDIF

      weight = weight * deltax / twopi
      weight = weight * xmg * (1 + tanx**2 )

      Resonance = .true.

      RETURN 
      END

****************************************************************************
      LOGICAL FUNCTION TwoBodyDecay(R1, R2, P_in, M2_in, M2_out1, M2_out2, 
     &                              P_out1, P_out2, weight)
****************************************************************************
*     Calculates the two-body decay P_in -> P_out1 + P_out2
*
*    INPUT
*       R1, R2         : two random numbers between 0 and 1
*       P_in           : 4-momentum of mother particle
*       M_in           : mass**2 of mother particle
*       M_out1, M_out2 : mass**2 of daughter particles
*       weight         : phasespace weight
*
*    OUTPUT
*       weight         : new phasespace weight (multiplied by Jacobi factor)
*       P_out1, P_out2 : 4-momenta of daughter particles
*       TwoBodyDecay   : .true.,  if everything is ok 
*                        .false., if kinematics not ok -> weight = 0
****************************************************************************
      IMPLICIT NONE

      DOUBLE PRECISION R1, R2, P_in(0:3), M2_in, M2_out1, M2_out2
      DOUBLE PRECISION P_out1(0:3), P_out2(0:3), weight
      DOUBLE PRECISION twopi
      PARAMETER (twopi = 2d0*3.141592653589793d0)
      INTEGER*4 mu
      double precision Q(0:3), CTH, STH, PHID, M_in, M_out1, M_out2, E1, P1

      M_in = sqrt(M2_in)
      M_out1 = sqrt(M2_out1)
      M_out2 = sqrt(M2_out2)

      IF (M_in .LT. (M_out1 + M_out2) ) THEN
         weight = 0d0 
         TwoBodyDecay = .false.
         RETURN 
      ENDIF

c generate outgoing momenta in CMS and boost back
      CTH = -1d0 + 2d0 * R1
      if ((1D0 - CTH**2) .le. 1D-10) then
         weight = 0D0
         return
      end if
      STH = sqrt(1d0 - CTH**2 )
      PHID = twopi * R2 - twopi/2.0d0
      
      E1 = (M2_in - M2_out2 + M2_out1) / (2d0*M_in)
      if ((E1**2 - M2_out1) .le. 1D-10) then
         weight = 0d0
         TwoBodyDecay = .false.
         return
      end if
      P1 = sqrt(E1**2 - M2_out1)
      Q(0) = E1
      Q(1) = P1 * STH*COS(PHID)
      Q(2) = P1 * STH*SIN(PHID)
      Q(3) = P1 * CTH

      call boostn(Q,P_in,P_out1)

      do mu = 0,3
         P_out2(mu) = P_in(mu) - P_out1(mu)
      enddo

      weight = weight * P1/(M_in * 2d0 * twopi) 
      
      TwoBodyDecay = .true.
      RETURN
      END

****************************************************************************
      LOGICAL FUNCTION TwoToJetsPlusX(N,R,RN,ECM,PTmin,XQ2,
     1                                K1,K2,X1,X2,Q,P,weight)
****************************************************************************
*     Calculates the phasespace for the 2 -> N massless jets + X
*     X is a pseudo-particle with arbitrary mass>0
*
*    INPUT
*       N      : number of massless jets      
*       R      : array with 3*N random numbers between 0 and 1
*       RN     : one additional random number 
*       ECM    : center-of-mass energy square
*       PTmin  : minimum PT of jets
*       XQ2    : Q**2 of X-particle
*       weight : phasespace weight
*
*    OUTPUT
*       K1,K2  : 4-vectors for 2 incoming partons
*       X1,X2  : corresponding values of Feynman x
*       Q      : 4-vector of X-particle
*       P      : 4-vectors for the N massless jets
*       weight : new phasespace weight (multiplied by Jacobi factor)
****************************************************************************
      IMPLICIT NONE
      INTEGER*4 N, N1
      PARAMETER (N1 = 2)
      DOUBLE PRECISION TPI,PI
      PARAMETER (PI=3.141592653589793d0)
      PARAMETER (TPI=2D0*PI )
      DOUBLE PRECISION R(*),RN,ECM,PTmin,XQ2,
     1       K1(0:3),K2(0:3),Q(0:3),P(0:3,N1),
     2       X1,X2,weight
      
      INTEGER*4 I,J
      DOUBLE PRECISION S,FAC,DELY,XMIN,XMAX,X,DELX,YCM,YMIN,YMAX,YBAR
      DOUBLE PRECISION PCM(0:3), PNJ(0:4), Y(N1), PT(N1), PHI(N1)
      DOUBLE PRECISION SINHY, COSHY, M2, M2min, QV, XMJ, PJT2, SHA, ACM
      DOUBLE PRECISION SUMPST, Q0ST, RSHAT
      PARAMETER (m2min = 0.1d0**2)

      double precision dotrr
      external dotrr

      S = ECM**2
      FAC = 1d0/(16d0*PI**3)
      weight = weight*TPI

C  Reset local variables and incoming momenta      
      DO I=0,3
         PCM(I) = 0d0
         PNJ(I) = 0d0
         K1(I) = 0d0
         K2(I) = 0d0
      ENDDO
      PNJ(4) = 0d0

C  Get 4-momenta of the N massless jets
      DO I = 1,N
         weight = weight*FAC
         J = 3*(I-1)
         
         if (R(J+2).gt.0.2) then
            XMIN = 2D0/ECM
            XMAX = 1/max(PTmin,10d0)
            DELX = XMAX-XMIN
            X = XMIN + (R(J+2)-0.2)/0.8d0*DELX
            PT(I) = 1/X
            weight = weight/0.8d0*DELX*PT(I)**3
         else
            xmin = 0
            xmax = max(PTmin,10d0)
            DELX = XMAX-XMIN
            PT(i) = XMIN + R(J+2)/0.2d0*DELX
            weight = weight/0.2d0*DELX*PT(I)
         endif
         YMAX = ECM/(2d0*PT(I))
         YMAX = LOG(YMAX+SQRT(YMAX**2-1d0))

         DELY = 2*YMAX

         IF (I.LE.3) THEN
            Y(I) = -YMAX+R(J+1)*DELY
         ELSE
C  fix random number to 0 or 1 for delta-y = 0 with previous jet: R=0 or 1 
C  corresponds to lowest dijet invariant mass
            IF (Y(I-1).LT.0D0) THEN
               IF (Y(I-1).LT.-YMAX) THEN
                  Y(I) = -YMAX + R(J+1)*DELY
               ELSEIF (Y(I-1).LT.0d0) THEN
                  Y(I) = Y(I-1) + R(J+1)*DELY
                  IF (Y(I).GT.YMAX) Y(I) = -2*YMAX + Y(I)
               ELSE
                  Y(I) =  - R(J+1)*DELY
                  IF (Y(I).LT.-YMAX) Y(I) = 2*YMAX + Y(I)
               ENDIF
            ELSE
               IF (Y(I-1).GT.YMAX) THEN
                  Y(I) = YMAX - R(J+1)*DELY
               ELSEIF (Y(I-1).GT.0D0) THEN
                  Y(I) = Y(I-1) - R(J+1)*DELY
                  IF (Y(I).LT.-YMAX) Y(I) = 2*YMAX + Y(I)
               ELSE
                  Y(I) = R(J+1)*DELY
                  IF (Y(I).GT.YMAX) Y(I) = -2*YMAX + Y(I)
               ENDIF
            ENDIF
         ENDIF

         SINHY = SINH(Y(I))
         COSHY = SQRT(1+SINHY**2)
         weight = weight*DELY

         P(0,I) = PT(I)*COSHY

         IF (I.EQ.1) THEN
            PHI(I) = PI*(2*RN-1)
         ELSE
            PHI(I) = PHI(I-1) + TPI*R(J+3)
            if (phi(i).gt.pi) phi(i) = phi(i)-tpi
         ENDIF

         weight = weight*TPI

         P(1,I) = PT(I)*COS(PHI(I))
         P(2,I) = PT(I)*SIN(PHI(I))
         P(3,I) = PT(I)*SINHY

         DO J = 0,3
            PNJ(J) = PNJ(J)+P(J,I)
         ENDDO
      ENDDO

C  determine N-jet invariant mass
      IF (N.GT.1) THEN
         PNJ(4) = PNJ(0)**2-PNJ(1)**2-PNJ(2)**2-PNJ(3)**2
      ELSE
         PNJ(4) = 1D-18
      ENDIF
      IF (PNJ(4).GE.S) THEN
         weight=0d0
         TwoToJetsPlusX = .false.
         RETURN
      ENDIF

c  check invariant masses of all parton pairs: eliminate events which 
c  may cause numerical problems 
      do I = 1,N-1
         do J = I+1,N
            m2 = 2*dotrr(p(0,I),p(0,J))
            if (m2.lt.m2min) then
               weight = 0d0
               TwoToJetsPlusX = .false.
               RETURN
            endif
         enddo
      enddo
      
C  Next get the center of mass rapidity YCM
      QV = SQRT(XQ2)
      XMJ = SQRT(PNJ(4))
      if (PNJ(0).le.PNJ(3)) then
        ybar = 100
      elseif (PNJ(0).le.-PNJ(3)) then
        ybar = -100
      else
        YBAR = 0.5D0*LOG( (PNJ(0)+PNJ(3))/(PNJ(0)-PNJ(3)) )
      endif
      PJT2 = PNJ(1)**2+PNJ(2)**2
      SHA = (S-(QV+XMJ)**2-4*PJT2) / (4*(PNJ(4)+PJT2)) 
      IF (SHA.LE.0D0) THEN
        weight=0d0
        TwoToJetsPlusX = .false.
        RETURN
      ENDIF
      SHA = SQRT(SHA)
      ACM = LOG( SHA+SQRT(SHA**2+1) ) 
      YMIN = YBAR - ACM
      YMAX = YBAR + ACM

      DELY = YMAX-YMIN
      YCM = YMIN + R(3)*DELY
      SINHY = SINH(YCM)
      COSHY = SQRT(1+SINHY**2)
      weight = weight*DELY

C  now get Q and SHAT etc.
      SUMPST = PNJ(1)**2+PNJ(2)**2+(PNJ(3)*COSHY-PNJ(0)*SINHY)**2
      Q0ST = SQRT(XQ2+SUMPST)
      RSHAT = Q0ST + SQRT(PNJ(4)+SUMPST)
      PCM(0) = RSHAT*COSHY
      PCM(3) = RSHAT*SINHY

      DO J = 0,3
        Q(J) = PCM(J)-PNJ(J)
      ENDDO

      X1 = (PCM(0)+PCM(3))/ECM
      X2 = (PCM(0)-PCM(3))/ECM
      IF (X1.GE.1.0D0 .OR. X2.GE.1.0D0) THEN
         weight=0d0
         TwoToJetsPlusX = .false.
         RETURN
      ENDIF
      K1(0) =  X1*ECM/2
      K1(3) =  K1(0)
      K2(0) =  X2*ECM/2
      K2(3) = -K2(0)

c  check initial-final momentum transfer for all partons: eliminate events 
c  which may cause numerical problems 
      do I = 1,N         
         m2 = 2*min( p(0,I)*k1(0)-p(3,I)*k1(3),
     1               p(0,I)*k2(0)-p(3,I)*k2(3) )
         if (m2.lt.m2min) then
            weight = 0d0
            TwoToJetsPlusX = .false.
            RETURN
         endif
      enddo

C  insert Jacobian for dtau * delta(PCM(0)-Ql(0,0)-PNJ(0))/(2 Ql(0,0))
      weight = weight * RSHAT/(S*Q0ST)
      
C  finally put in flux factor for the parton collisions and conversion to fb
!      weight = weight / (2*RSHAT**2) * 3.89379304D+11

      TwoToJetsPlusX = .true.
      RETURN
      END


!      subroutine born_phsp_hhjj(xborn)
!      implicit none
!!       include 'nlegborn.h'
!!       include 'pwhg_flst.h'
!      include 'brinclude.h'
!      include 'pwhg_kn.h'
!      include 'pwhg_par.h'
!      include 'pwhg_math.h'
!      include 'PhysPars.h'
!      include 'PhysPars_Higgs.h'
!      include 'incl2pwhg.h'
!      real * 8 xborn(ndiminteg-6)
!      real * 8 xjac,taumin,tau,y,beta,betaCM,vec(3),cth,s,
!     #     z,zhigh,zlow,n,m2,m2H1,m2H2
!      integer mu,k,j
!      parameter (n = 8d0)
!      logical ini
!      data ini/.true./
!      save ini
!      real * 8 Vmass2,Vmass2low,Vmass2high,VmVw  
!      real * 8 m2jj,pV(0:3),pVmod,pVmod2,pJ(0:3,2),cthj,phij,pcmjj(0:3),
!     #     pcmjjmod,ptmp(0:3,2), weight, m22,pt3cut
!      real * 8 mass
!      real * 8 pH1(0:3), pH2(0:3), mH1, mH2, wt
!      real * 8 pH13m(1:4), pH23m(1:4), pV3m(1:4)
!      external mass      
!      logical check
!      parameter(check=.false.)
!      logical higgs_use_BW
!      common/useBW/higgs_use_BW
!c      parameter (BW=.true.)
!      real * 8 epsilon
!      parameter (epsilon=1d-10)
!      real * 8 pt1cut,pt2cut,pt1,pt2,m2jjmin
!      real * 8 BW_fixed,BW_running
!      logical higgsfixedwidth
!      save higgsfixedwidth
!      real * 8 powheginput
!      external powheginput
!
!
!      Vmass2 = ph_Hmass2
!      Vmass2low = ph_Hmass2low
!      Vmass2high = ph_Hmass2high
!      VmVw = ph_HmHw
!      
!!     VBFHHMOD: Removed narrow-width option for internal Higgs
!      zlow=atan((Vmass2low  - Vmass2)/VmVw)
!      zhigh=atan((min(Vmass2high,kn_sbeams)  - Vmass2)/VmVw)
!      z=zlow+(zhigh-zlow)*xborn(1)
!      xjac=zhigh-zlow
!      m2=VmVw*tan(z)+Vmass2
!c     The BW integrates to Pi ==> divide by Pi
!      xjac=xjac/pi
!
!
!
!c     d x1 d x2 = d tau d y;
!      taumin=m2/kn_sbeams
!!      tau=exp(log(taumin)*(1-xborn(2)**2))
!!      xjac=xjac*tau*abs(log(taumin))*2*xborn(2)
!      tau=exp(log(taumin)*(1-xborn(2)**(1.5d0)))
!      xjac=xjac*tau*abs(log(taumin))*1.5d0*xborn(2)**(0.5d0)
!      s=kn_sbeams*tau
!      brkn_sborn=s
!c     ymax=|log(tau)|/2
!      y=-(1-2*xborn(3))*log(tau)/2
!      xjac=-xjac*log(tau)
!
!c     generate dijet squared mass m2jj
!      m2jj = xborn(4)*(sqrt(s)-sqrt(m2))**2
!      xjac = xjac * (sqrt(s)-sqrt(m2))**2
!      
!      pV(0)=(s+m2-m2jj)/(2*sqrt(s))
!      pVmod2 = pV(0)**2-m2
!      
!      if (pVmod2.lt.epsilon) then
!         pVmod = 0d0
!      else
!         pVmod = sqrt(pVmod2)
!      endif
!      
!      z=1-2*xborn(5)
!      xjac=xjac*2
!      cth=1.5d0*(z-z**3/3)
!      xjac=xjac*1.5d0*(1-z**2)
!
!      pV(1) = pVmod*sqrt(1-cth**2)
!      pV(2) = 0d0
!      pV(3) = pVmod*cth
!
!      
!!     VBFHHMOD: mass of the two final state higgses
!!     todo: implement higgs BW and decays
!!     AK: I think this is how to do it.
!      if (higgs_use_BW) then
!         zlow=atan((Vmass2low  - Vmass2)/VmVw)
!         zhigh=atan((min(Vmass2high,kn_sbeams)  - Vmass2)/VmVw)
!         z=zlow+(zhigh-zlow)*xborn(ndiminteg-7)
!         xjac=zhigh-zlow
!         mH1=VmVw*tan(z)+Vmass2
!c     The BW integrates to Pi ==> divide by Pi
!         xjac=xjac/pi
!
!         zlow=atan((Vmass2low  - Vmass2)/VmVw)
!         zhigh=atan((min(Vmass2high,kn_sbeams)  - Vmass2)/VmVw)
!         z=zlow+(zhigh-zlow)*xborn(ndiminteg-6)
!         xjac=zhigh-zlow
!         mH2=VmVw*tan(z)+Vmass2
!c     The BW integrates to Pi ==> divide by Pi
!         xjac=xjac/pi
!
!         write(*,*) ' error in Born phase space: BW not implemented'
!         call pwhg_exit(-1)
!      else
!         xjac = 1
!         mH1 = ph_Hmass
!         mH2 = ph_Hmass
!      endif
!
!!     The phi3m routine assumes momenta of the form p(1,2,3,4)
!      pV3m(1:3) = pV(1:3)
!      pV3m(4) = pV(0) 
!     
!!     AK+FD: Schematic of decay
!!     xborn(ndiminteg-7) and xborn(ndiminteg-6) are for BW higgs masses
!      call phi3m(xborn(ndiminteg-9),xborn(ndiminteg-8),pV3m,pH13m,pH23m
!     $     ,mH1,mH2,wt)
!      xjac=xjac*wt
!
!      pH1(1:3) = pH13m(1:3)
!      pH1(0) = pH13m(4)
!      pH2(1:3) = pH23m(1:3)
!      pH2(0) = pH23m(4)
!
!     
!c     supply 2 pi for azimuthal integration (not performed)
!      xjac=xjac*2*pi
!c     supply the other factors to the jacobian
!c     factor for the two-body jet phase space
!      xjac=xjac/(8*(2*pi)**2)
!c     factor for V production
!      xjac=xjac/(4*(2*pi)**3)*pVmod/sqrt(s)
!
!c     Build kinematics
!      brkn_xb1=sqrt(tau)*exp(y)
!      brkn_xb2=tau/brkn_xb1
!c     boost back in the lab frame
!c     now boost everything along 3rd axis
!      betaCM=(brkn_xb1-brkn_xb2)/(brkn_xb1+brkn_xb2)
!      vec(1)=0
!      vec(2)=0
!      vec(3)=1
!      
!!     VBFHHMOD: This needs to be changed to pH1, pH2 going to
!!     brkn_pborn(:,3) & prkn_pborn(:,4)
!!     call mboost(1,vec,betaCM,pV,brkn_pborn(0,3))      
!
!      call mboost(1,vec,betaCM,pH1,brkn_pborn(0,3))      
!      call mboost(1,vec,betaCM,pH2,brkn_pborn(0,4))      
!
!c     build jet momenta in the jet CM frame
!      pJ(0,1) = sqrt(m2jj)/2
!      pJ(0,2) = pJ(0,1)
!c     azimuth and polar angle of a jet
!c      z=1-2*xborn(6)**2
!c      xjac=xjac*4*xborn(6)
!      z=1-2*xborn(6)**4
!      xjac=xjac*8*xborn(6)**3
!      cthj=1.5d0*(z-z**3/3)
!      xjac=xjac*1.5d0*(1-z**2)
!
!      phij = 2*pi*xborn(7)
!      xjac=xjac*2*pi
!
!      brkn_jacborn = xjac
!      
!      pJ(1,1) = pJ(0,1)*sqrt(1-cthj**2)*sin(phij)
!      pJ(2,1) = pJ(0,1)*sqrt(1-cthj**2)*cos(phij)
!      pJ(3,1) = pJ(0,1)*cthj
!      do mu=1,3
!         pJ(mu,2)=-pJ(mu,1)
!      enddo
!
!      do mu=0,3
!         brkn_pborn(mu,1)=brkn_xb1*kn_beams(mu,1)
!         brkn_pborn(mu,2)=brkn_xb2*kn_beams(mu,2)
!      enddo
!
!      if (check) then
!         call mboost(2,vec,-betaCM,brkn_pborn(0,1),ptmp(0,1))      
!         write(*,*) 'CM vec1',(ptmp(mu,1),mu=0,3)
!         write(*,*) 'CM vec2',(ptmp(mu,2),mu=0,3)
!      endif
!
!c     boost in the lab frame
!c     compute first p_plus+p_minus-pV
!      do mu=0,3
!      !VBFHHMOD: need to add second higgs 
!         pcmjj(mu)= brkn_pborn(mu,1) + brkn_pborn(mu,2) - brkn_pborn(mu
!     $        ,3) - brkn_pborn(mu,4)
!      enddo
!      pcmjjmod = sqrt(pcmjj(1)**2+pcmjj(2)**2+pcmjj(3)**2)
!c     recompute pcmjj(0) from m2jj, otherwise there are points where 
!c     beta > 1 or beta < 0
!      pcmjj(0) = sqrt(m2jj+pcmjjmod**2)
!
!      beta=pcmjjmod/pcmjj(0)
!
!      do mu=1,3
!         vec(mu)=pcmjj(mu)/pcmjjmod
!      enddo
!
!      !VBFHHMOD: need to switch to index 5
!!      call mboost(2,vec,beta,pJ(0,1),brkn_pborn(0,4))
!      call mboost(2,vec,beta,pJ(0,1),brkn_pborn(0,5))
!
!      if (check) then
!         call mboost(1,vec,-beta,pcmjj(0),ptmp(0,1))
!         write(*,*) 'only time component ==> ',(ptmp(mu,1),mu=0,3)
!         write(*,*) ''
!         write(*,*) 'new set'
!         do j=1,br_nlegborn
!            write(*,*) 'mom ',j,(brkn_pborn(mu,j),mu=0,3)
!            write(*,*) 'mass ',j,mass(brkn_pborn(0,j))
!         enddo
!         call checkmomzero(br_nlegborn,brkn_pborn)
!      endif
!      
!c     boost in the CM frame
!      vec(1)=0
!      vec(2)=0
!      vec(3)=1
!      call mboost(br_nlegborn,vec,-betaCM,brkn_pborn(0,1),brkn_cmpborn(0,1))      
!      if (check) then
!         write(*,*) ''
!         write(*,*) 'new set'     
!         do j=1,nlegborn
!            write(*,*) 'mom ',j,(brkn_cmpborn(mu,j),mu=0,3)
!         enddo
!      endif
!
!      call br_compute_csimax_fsr
!      
!      end

