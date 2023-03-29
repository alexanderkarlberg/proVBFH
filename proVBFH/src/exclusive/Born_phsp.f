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
      logical higgs_use_BW, generation_cuts
      save generation_cuts
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

c     set initial- and final-state masses for Born and real
         higgsfixedwidth=powheginput("#higgsfixedwidth").gt.0
         ini=.false.

         if(powheginput("#ptcut").gt.0d0) then
            pt1cut = powheginput("#ptcut")
            pt2cut = powheginput("#ptcut")
            pt3cut = powheginput("#ptcut")
         else
            pt1cut = 0d0
            pt2cut = 0d0
            pt3cut = 0d0
         endif

c     determine whether the Higgs boson propagator is
C     narrow width or breit wigner
         higgs_use_BW = .false.
         if (powheginput("#higgsbreitwigner").eq.1) higgs_use_BW=.true.
         
         if(powheginput("#mjjmaxcut").gt.0d0) then
            mjj_maxcut = powheginput("#mjjmaxcut")**2
         else
            mjj_maxcut=0d0
         endif
         if(powheginput("#mjjcut").gt.0d0) then
            m2jjmin = powheginput("#mjjcut")**2
         else
            m2jjmin = 0d0
         endif
         
         m2min=0d0
c         generation_cuts = (m2min.ne.0d0)
         generation_cuts = .false.

         if ((pt1cut.ne.0d0).or.(pt2cut.ne.0d0).or.(pt3cut.ne.0d0)
     &                      .or.(m2jjmin.ne.0)) then
            write(*,*) '***************************************'
            write(*,*) '****  GENERATION CUTS ACTIVATED  ! ****' 
            write(*,*) '***************************************'
            write(*,*) '****  ptj>',pt1cut,            '   ****' 
            write(*,*) '****  mj1j2>',sqrt(mjj_maxcut),'   ****' 
            write(*,*) '****  mjj>',sqrt(m2jjmin),     '   ****' 
            write(*,*) '****  m2min>',sqrt(m2min),     '   ****'             
            write(*,*) '***************************************' 
            generation_cuts = .true.
         endif
      endif      

      brkn_ktmin=0
      kn_ktmin=0

      call born_phsp_hjj(xborn)

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
      
      if(generation_cuts) then

         pt1 = sqrt(kn_pborn(1,4)**2+kn_pborn(2,4)**2)
         pt2 = sqrt(kn_pborn(1,5)**2+kn_pborn(2,5)**2)  
         pt3 = sqrt(kn_pborn(1,6)**2+kn_pborn(2,6)**2)  
         mjj1=dotp(kn_pborn(0,4),kn_pborn(0,6))
         mjj2=dotp(kn_pborn(0,5),kn_pborn(0,6))
         mjj3=dotp(kn_pborn(0,4),kn_pborn(0,5))
         mjj_max = max(mjj1,mjj2)
         mjj_max = max(mjj_max,mjj3)
      
         if(mjj_max.lt.mjj_maxcut) then
            kn_jacborn=0d0
         endif
         
         if ((pt1.lt.pt1cut).or.(pt2.lt.pt2cut).or.(pt3.lt.pt3cut)) then
            kn_jacborn=0d0
         endif      
         
         if ((mjj1.lt.m2jjmin).or.(mjj2.lt.m2jjmin).or.(mjj3.lt.m2jjmin)) then
            kn_jacborn=0d0
         endif  
         
         do i=1,2
            do j=i+1,3
               m2=2d0*abs(dotp(kn_pborn(0,i+3),kn_pborn(0,j+3)))
               if(m2.lt.m2min) kn_jacborn=0d0
            enddo
         enddo
         
         do i=1,3
            m2=min(abs(2d0*dotp(kn_pborn(0,i+3),kn_pborn(0,1))),
     &           abs(2d0*dotp(kn_pborn(0,i+3),kn_pborn(0,2)))) 
            if(m2.lt.m2min) kn_jacborn=0d0
            
         enddo
      endif ! Generation cuts

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


      subroutine born_phsp_hjj(xborn)
      implicit none
!       include 'nlegborn.h'
!       include 'pwhg_flst.h'
      include 'brinclude.h'
      include 'pwhg_kn.h'
      include 'pwhg_par.h'
      include 'pwhg_math.h'
      include 'PhysPars.h'
      include 'PhysPars_Higgs.h'
      include 'incl2pwhg.h'
      real * 8 xborn(ndiminteg-6)
      real * 8 m2,xjac,taumin,tau,y,beta,betaCM,vec(3),cth,s,
     #     z,zhigh,zlow,n
      integer mu,k,j
      parameter (n = 8d0)
      logical ini
      data ini/.true./
      save ini
      real * 8 Vmass2,Vmass2low,Vmass2high,VmVw  
      real * 8 m2jj,pV(0:3),pVmod,pVmod2,pJ(0:3,2),cthj,phij,pcmjj(0:3),
     #     pcmjjmod,ptmp(0:3,2), weight, m22,pt3cut
      real * 8 mass      
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


      Vmass2 = ph_Hmass2
      Vmass2low = ph_Hmass2low
      Vmass2high = ph_Hmass2high
      VmVw = ph_HmHw

      if (higgs_use_BW) then
         zlow=atan((Vmass2low  - Vmass2)/VmVw)
         zhigh=atan((min(Vmass2high,kn_sbeams)  - Vmass2)/VmVw)
         z=zlow+(zhigh-zlow)*xborn(1)
         xjac=zhigh-zlow
         m2=VmVw*tan(z)+Vmass2
c     The BW integrates to Pi ==> divide by Pi
         xjac=xjac/pi
      else
         xjac = 1
         m2 = Vmass2
      endif



c     d x1 d x2 = d tau d y;
      taumin=m2/kn_sbeams
!      tau=exp(log(taumin)*(1-xborn(2)**2))
!      xjac=xjac*tau*abs(log(taumin))*2*xborn(2)
      tau=exp(log(taumin)*(1-xborn(2)**(1.5d0)))
      xjac=xjac*tau*abs(log(taumin))*1.5d0*xborn(2)**(0.5d0)
      s=kn_sbeams*tau
      brkn_sborn=s
c     ymax=|log(tau)|/2
      y=-(1-2*xborn(3))*log(tau)/2
      xjac=-xjac*log(tau)

c     generate dijet squared mass m2jj
      m2jj = xborn(4)*(sqrt(s)-sqrt(m2))**2
      xjac = xjac * (sqrt(s)-sqrt(m2))**2
      
      pV(0)=(s+m2-m2jj)/(2*sqrt(s))
      pVmod2 = pV(0)**2-m2
      
      if (pVmod2.lt.epsilon) then
         pVmod = 0d0
      else
         pVmod = sqrt(pVmod2)
      endif
      
      z=1-2*xborn(5)
      xjac=xjac*2
      cth=1.5d0*(z-z**3/3)
      xjac=xjac*1.5d0*(1-z**2)

      pV(1) = pVmod*sqrt(1-cth**2)
      pV(2) = 0d0
      pV(3) = pVmod*cth

c     supply 2 pi for azimuthal integration (not performed)
      xjac=xjac*2*pi
c     supply the other factors to the jacobian
c     factor for the two-body jet phase space
      xjac=xjac/(8*(2*pi)**2)
c     factor for V production
      xjac=xjac/(4*(2*pi)**3)*pVmod/sqrt(s)

c     Build kinematics
      brkn_xb1=sqrt(tau)*exp(y)
      brkn_xb2=tau/brkn_xb1
c     boost back in the lab frame
c     now boost everything along 3rd axis
      betaCM=(brkn_xb1-brkn_xb2)/(brkn_xb1+brkn_xb2)
      vec(1)=0
      vec(2)=0
      vec(3)=1
      call mboost(1,vec,betaCM,pV,brkn_pborn(0,3))      

c     build jet momenta in the jet CM frame
      pJ(0,1) = sqrt(m2jj)/2
      pJ(0,2) = pJ(0,1)
c     azimuth and polar angle of a jet
c      z=1-2*xborn(6)**2
c      xjac=xjac*4*xborn(6)
      z=1-2*xborn(6)**4
      xjac=xjac*8*xborn(6)**3
      cthj=1.5d0*(z-z**3/3)
      xjac=xjac*1.5d0*(1-z**2)

      phij = 2*pi*xborn(7)
      xjac=xjac*2*pi

      brkn_jacborn = xjac
      
      pJ(1,1) = pJ(0,1)*sqrt(1-cthj**2)*sin(phij)
      pJ(2,1) = pJ(0,1)*sqrt(1-cthj**2)*cos(phij)
      pJ(3,1) = pJ(0,1)*cthj
      do mu=1,3
         pJ(mu,2)=-pJ(mu,1)
      enddo

      do mu=0,3
         brkn_pborn(mu,1)=brkn_xb1*kn_beams(mu,1)
         brkn_pborn(mu,2)=brkn_xb2*kn_beams(mu,2)
      enddo

      if (check) then
         call mboost(2,vec,-betaCM,brkn_pborn(0,1),ptmp(0,1))      
         write(*,*) 'CM vec1',(ptmp(mu,1),mu=0,3)
         write(*,*) 'CM vec2',(ptmp(mu,2),mu=0,3)
      endif

c     boost in the lab frame
c     compute first p_plus+p_minus-pV
      do mu=0,3
         pcmjj(mu)= brkn_pborn(mu,1) + brkn_pborn(mu,2) - brkn_pborn(mu,3)
      enddo
      pcmjjmod = sqrt(pcmjj(1)**2+pcmjj(2)**2+pcmjj(3)**2)
c     recompute pcmjj(0) from m2jj, otherwise there are points where 
c     beta > 1 or beta < 0
      pcmjj(0) = sqrt(m2jj+pcmjjmod**2)

      beta=pcmjjmod/pcmjj(0)

      do mu=1,3
         vec(mu)=pcmjj(mu)/pcmjjmod
      enddo

      call mboost(2,vec,beta,pJ(0,1),brkn_pborn(0,4))

      if (check) then
         call mboost(1,vec,-beta,pcmjj(0),ptmp(0,1))
         write(*,*) 'only time component ==> ',(ptmp(mu,1),mu=0,3)
         write(*,*) ''
         write(*,*) 'new set'
         do j=1,br_nlegborn
            write(*,*) 'mom ',j,(brkn_pborn(mu,j),mu=0,3)
            write(*,*) 'mass ',j,mass(brkn_pborn(0,j))
         enddo
         call checkmomzero(br_nlegborn,brkn_pborn)
      endif
      
c     boost in the CM frame
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
      do j=4,br_nlegreal
         y=1-dotp(brkn_cmpreal(0,1),brkn_cmpreal(0,j))
     1 /(brkn_cmpreal(0,1)*brkn_cmpreal(0,j))
         brkn_dijterm(0,j)=(brkn_cmpreal(0,j)**2
     1 *(1-y**2))**brpar_diexp
         brkn_dijterm(1,j)=(brkn_cmpreal(0,j)**2
     1 *2*(1-y))**brpar_diexp
         brkn_dijterm(2,j)=(brkn_cmpreal(0,j)**2
     1 *2*(1+y))**brpar_diexp
      enddo
      do j=4,br_nlegreal
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
      do j=4,5
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
            if (powheginput('#runningscales').eq.2d0) then 
               print*, 'Using Q1 and Q2 as scale'
               flg_minlo = .true.
               flg_minlo_real = .false.
               runningscales = .false.
            endif
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
         if(flg_btildepart.eq.'r') then ! Use real momenta
            ptH = sqrt(kn_preal(1,3)**2+kn_preal(2,3)**2) ! ptH
         else ! Otherwise, always use born
            ptH = sqrt(kn_pborn(1,3)**2+kn_pborn(2,3)**2) ! ptH
         endif
         
         muref = ( ( ph_Hmass * 0.5d0 ) ** 4d0
     $        + ( ph_Hmass * ptH * 0.5d0 ) ** 2d0 
     $        ) ** 0.25d0
         
      else
!         muref=ph_Hmass
         muref=ph_Wmass
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
