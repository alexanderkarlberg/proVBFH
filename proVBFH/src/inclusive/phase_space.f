C     Module phase_space is for the main part taken from
C       POWHEG-BOX/VBF_H/Born_phsp.f
C     Usage:
C       1) Set up the beams with set_beams(sqrts). [only once]
C       2) Generate phase space with gen_phsp(xborn) using an array
C          of random numbers. This fills the internal variables
C          contained in pwhg_kn.h with new values
C       3) Map internal variables to the ones of the module with
C          set_phsp().
C       4) If using a powheg analysis, map module variables to
C          phep(:,:), the array in hepevt.h, used in powheg.
      module phase_space

      logical fill_plots
      double precision x1, x2, vq1(0:3), vq2(0:3), phi, jacobian
      double precision Q1_sq, q1_perp, q1_perp_sq, Q2_sq, q2_perp
      
      integer nlegborn, nlegreal
      parameter (nlegborn = 5)
      parameter (nlegreal = nlegborn + 1)
      
      contains

  
C------------------------------------------------------------
C     dsigma function
C     calculate the cross section, using the phase space from
C     POWHEG-BOX/VBF_H/Born_phsp.f and matrix element computed with Hoppet
      double precision function dsigma(xrand, vegas_weight)
      use matrix_element
      use incl_parameters
      implicit none
      include 'pwhg_kn.h'
C     xrand contains a vector of random numbers in [0,1]
      double precision xrand(7), vegas_weight
      double precision ptH
      integer vegas_ncall
      common/vegas_ncall/vegas_ncall

      dsigma = 0d0
C     generate phase space using phase_space module
      call gen_phsp(xrand(1:7))
C     set phase space variables x1, x2, vq1, vq2, phi, Q1_sq,
C     q1_perp, q1_perp_sq, Q2_sq, q2_perp using generated phase space
      call set_phsp()
      
C     skip phase space points with vanishing jacobian or with minQ1,Q2) < Qmin
      if ((jacobian.ne.0d0).and.(min(Q1_sq,Q2_sq).gt.(Qmin**2))) then
C     For scale_choice=3 we need ptH. We have in the transverse plane
C     q1(1:2) = (q1_perp         ,    0            )
C     q2(1:2) = (q2_perp*cos(phi), q2_perp*sin(phi))
C     pH(:)   = - q1(:) - q2(:)
C     Hence ptH is given by
         ptH = sqrt((q1_perp + q2_perp*cos(phi))**2
     $        + (q2_perp*sin(phi))**2)
C     compute dsigma using the squared hadronic tensor
!         dsigma = eval_matrix_element_tensor(order_min,order_max, x1, x2, 
!     $        kn_beams(:,1), kn_beams(:,2), vq1, vq2, ptH)
         dsigma = eval_matrix_element(order_min,order_max, x1, x2, 
     $        kn_beams(:,1), kn_beams(:,2), vq1, vq2, ptH)
C     convert to [pb] and add in jacobian
         dsigma = dsigma * gev2pb * jacobian

      endif
      
C     map to powheg variables for use by analysis
      call set_phep()
C     remove the rare outliers where we get dsigma = NaN
      if (dsigma.ne.dsigma) dsigma = 0d0
      
C     fill histograms using powheg analysis if requested
C     (only at LO when running exclusive)
      if (fill_plots) then
         call user_analysis(dsigma*vegas_ncall*vegas_weight)
         call pwhgaccumup
      endif

      end function dsigma

C----------------------------------------------------------------------
C     set_beams(sqrts)
C     set kn_beams(:,:) and beams(:,:) to required value
C     NB: This only needs to be called once
      subroutine set_beams(sqrts)
      implicit none
      include 'pwhg_kn.h'
      double precision sqrts
      kn_sbeams = sqrts**2
      kn_beams = 0d0
      kn_beams(0,1) = sqrts/2d0
      kn_beams(0,2) = sqrts/2d0
      kn_beams(3,1) = sqrts/2d0
      kn_beams(3,2) = -sqrts/2d0
      end subroutine

C----------------------------------------------------------------------
C     set_phsp()
C     set phase space variables
C       x1, x2, vq1, vq2, phi, Q1_sq, q1_perp, q1_perp_sq, Q2_sq, q2_perp
C     to values generated with gen_phsp
      subroutine set_phsp()
      implicit none
      include 'pwhg_kn.h'
      double precision acos_arg
      x1 = kn_xb1
      x2 = kn_xb2
      vq1(:) = kn_pborn(:,4) - kn_pborn(:,1)
      vq2(:) = kn_pborn(:,5) - kn_pborn(:,2)
      Q1_sq = invmass(vq1)**2
      Q2_sq = invmass(vq2)**2
      q1_perp_sq = vq1(1)**2 + vq1(2)**2
      q1_perp = sqrt(q1_perp_sq)
      q2_perp = sqrt(vq2(1)**2 + vq2(2)**2)
C     phi = arccos(vq1 . vq2 / q1_T q2_T)
      acos_arg = (vq1(1)*vq2(1) + vq1(2)*vq2(2))/(q1_perp*q2_perp)
C     protection against acos_arg being larger than one
      if (abs(acos_arg) > 1d0) then
C     write(6,*)  "WARNING: |acos_arg| > 1; acos_arg = ", acos_arg
         acos_arg = sign(1d0, acos_arg)
      endif
      phi = acos(acos_arg)
      jacobian = kn_jacborn * Q1_sq * Q2_sq / (x1 * x2) 
      end subroutine

C----------------------------------------------------------------------
C     gen_phsp(xborn(1:7))
C     generated phase space using an array of [0,1] random numbers
C     this code is taken from POWHEG-BOX/VBF_H
      subroutine gen_phsp(xborn)
      ! rename s, pi to avoid conflict with internal variables
      use incl_parameters, sval => s, pival => pi 
      implicit none
      include 'pwhg_kn.h'
      include 'pwhg_math.h'
      include 'PhysPars.h'
      real * 8 xborn(7)
      real * 8 m2,xjac,taumin,tau,y,beta,betaCM,vec(3),cth,s,
     #     z,zhigh,zlow,m
      integer mu,k,j
      logical ini
      data ini/.true./
      save ini
      real * 8 Vmass2,Vmass2low,Vmass2high,VmVw  
      real * 8 m2jj,pV(0:3),pVmod,pVmod2,pJ(0:3,2),cthj,phij,pcmjj(0:3),
     #     pcmjjmod,ptmp(0:3,2)
      real * 8 mass      
      external mass      
      logical check
      parameter(check=.false.)
      real * 8 epsilon
      parameter (epsilon=1d-10)
      real * 8 pt1cut,pt2cut,pt1,pt2,m2jjmin
      real * 8 BW_fixed,BW_running
      real * 8 weight
      integer BWflag
      integer nlegborn, nlegreal
      parameter (nlegborn = 5)
      parameter (nlegreal = nlegborn + 1)
      if(ini) then
c     set initial- and final-state masses for Born and real
         do k=1,nlegborn
            kn_masses(k)=0
         enddo
         kn_masses(nlegreal)=0
         if (complexpole) then
            if (higgsfixwdth) then
               write(*,*) ' '
               write(*,*) '*******************************************'
               write(*,*) 'The higgsfixedwidth and complexpolescheme '//
     $              'flags are both set to true.'
               write(*,*) 'The two flags are incompatible.'
               write(*,*) 'The POWHEG BOX exits.'
               write(*,*) '*******************************************'
               call exit(1)
            endif
            write(*,*) '*******************************************'
            write(*,*) '*******************************************'
            write(*,*) '****        COMPLEX-POLE SCHEME        ****'
            write(*,*) '****          Passarino et al          ****'
            write(*,*) '*******************************************'
            write(*,*) '*******************************************'
         endif
         ini=.false.
         pt1cut = 0d0
         pt2cut = 0d0
c         m2jjmin = 0.1d0**2
         m2jjmin = 0d0
         if ((pt1cut.ne.0d0).or.(pt2cut.ne.0d0).or.(m2jjmin.ne.0)) then
            write(*,*) '*************************************'
            write(*,*) '****       CUTS IN PLACE!!!      ****' 
            write(*,*) '*************************************'
         endif
      endif

      Vmass2 = mh_sq
      Vmass2low = (max(0d0,mh-hmasswindow*hwidth))**2
      Vmass2high = (mh+hmasswindow*hwidth)**2
      VmVw = mh*hwidth

      if (higgs_use_BW) then
         zlow=atan((Vmass2low  - Vmass2)/VmVw)
         zhigh=atan((min(Vmass2high,kn_sbeams)  - Vmass2)/VmVw)
         z=zlow+(zhigh-zlow)*xborn(1)
         xjac=zhigh-zlow
         m2=VmVw*tan(z)+Vmass2
c         print*, 'zlow,zhigh,z,xjac,m2', zlow,zhigh,z,xjac,m2
c     The BW integrates to Pi ==> divide by Pi
         xjac=xjac/pi
c         if(.not.higgsfixedwidth.and..not.complexpole) then
c     running width
c            BW_fixed=ph_HmHw/((m2-ph_Hmass2)**2 + ph_HmHw**2)
c            BW_running= (m2*ph_Hwidth/ph_Hmass) /
c     $           ((m2-ph_Hmass2)**2+(m2*ph_Hwidth/ph_Hmass)**2)
c            xjac = xjac * BW_running/BW_fixed
c         endif
      else
         xjac = 1
         m2 = Vmass2
      endif

      kn_masses(3)=sqrt(m2)

c     d x1 d x2 = d tau d y;
      taumin=m2/kn_sbeams
      tau=exp(log(taumin)*(1-xborn(2)**2))
      xjac=xjac*tau*abs(log(taumin))*2*xborn(2)
      s=kn_sbeams*tau
      kn_sborn=s
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
c         write(*,*) '============================>',pVmod2
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
      kn_xb1=sqrt(tau)*exp(y)
      kn_xb2=tau/kn_xb1
c     boost back in the lab frame
c     now boost everything along 3rd axis
      betaCM=(kn_xb1-kn_xb2)/(kn_xb1+kn_xb2)
      vec(1)=0
      vec(2)=0
      vec(3)=1
      call mboost(1,vec,betaCM,pV,kn_pborn(0,3))      

c     build jet momenta in the jet CM frame
      pJ(0,1) = sqrt(m2jj)/2
      pJ(0,2) = pJ(0,1)
c     azimuth and polar angle of a jet
      z=1-2*xborn(6)**2
      xjac=xjac*4*xborn(6)
      cthj=1.5d0*(z-z**3/3)
      xjac=xjac*1.5d0*(1-z**2)

      phij = 2*pi*xborn(7)
      xjac=xjac*2*pi

      kn_jacborn = xjac
      
      pJ(1,1) = pJ(0,1)*sqrt(1-cthj**2)*sin(phij)
      pJ(2,1) = pJ(0,1)*sqrt(1-cthj**2)*cos(phij)
      pJ(3,1) = pJ(0,1)*cthj
      do mu=1,3
         pJ(mu,2)=-pJ(mu,1)
      enddo

      do mu=0,3
         kn_pborn(mu,1)=kn_xb1*kn_beams(mu,1)
         kn_pborn(mu,2)=kn_xb2*kn_beams(mu,2)
      enddo

      if (check) then
         call mboost(2,vec,-betaCM,kn_pborn(0,1),ptmp(0,1))      
         write(*,*) 'CM vec1',(ptmp(mu,1),mu=0,3)
         write(*,*) 'CM vec2',(ptmp(mu,2),mu=0,3)
      endif

c     boost in the lab frame
c     compute first p_plus+p_minus-pV
      do mu=0,3
         pcmjj(mu)= kn_pborn(mu,1) + kn_pborn(mu,2) - kn_pborn(mu,3)
      enddo
      pcmjjmod = sqrt(pcmjj(1)**2+pcmjj(2)**2+pcmjj(3)**2)
c     recompute pcmjj(0) from m2jj, otherwise there are points where 
c     beta > 1 or beta < 0
      pcmjj(0) = sqrt(m2jj+pcmjjmod**2)

c      write(*,*) '1 ===========> ',(pcmjj(0)**2-pcmjjmod**2)/m2jj

      beta=pcmjjmod/pcmjj(0)

      do mu=1,3
         vec(mu)=pcmjj(mu)/pcmjjmod
      enddo

      call mboost(2,vec,beta,pJ(0,1),kn_pborn(0,4))

c     CUTS!!!!!!!!
      pt1 = sqrt(kn_pborn(1,4)**2+kn_pborn(2,4)**2)
      pt2 = sqrt(kn_pborn(1,5)**2+kn_pborn(2,5)**2)           
      if ((pt1.lt.pt1cut).or.(pt2.lt.pt2cut).or.m2jj.lt.m2jjmin) then
         write(*,*) 'setting jacobian to 0'
         kn_jacborn=0d0
      endif


      if (check) then
         call mboost(1,vec,-beta,pcmjj(0),ptmp(0,1))
         write(*,*) 'only time component ==> ',(ptmp(mu,1),mu=0,3)
         write(*,*) ''
         write(*,*) 'new set'
         do j=1,nlegborn
            write(*,*) 'mom ',j,(kn_pborn(mu,j),mu=0,3)
            write(*,*) 'mass ',j,mass(kn_pborn(0,j))
         enddo
         call checkmomzero(nlegborn,kn_pborn)
      endif
      
c     boost in the CM frame
      vec(1)=0
      vec(2)=0
      vec(3)=1
      call mboost(nlegborn,vec,-betaCM,kn_pborn(0,1),kn_cmpborn(0,1))      
      if (check) then
         write(*,*) ''
         write(*,*) 'new set'     
         do j=1,nlegborn
            write(*,*) 'mom ',j,(kn_cmpborn(mu,j),mu=0,3)
         enddo
      endif

      kn_minmass=sqrt(Vmass2low)
      
      end subroutine


      function invmass(p)
      implicit none
      real * 8 p(0:3),invmass
      invmass = sqrt(abs(p(0)**2-p(1)**2-p(2)**2-p(3)**2))
      end function

      
      subroutine mboost(m,vec,beta,vin,vout)
c     boosts the m vectors vin(0:3,m) into the vectors vout(0:3,m) (that can
c     be the same) in the direction of vec(3) (|vec|=1) with velocity
c     beta.  Lorents convention: (t,x,y,z).
      implicit none
      integer m
      real * 8 vec(3),beta,vin(0:3,m),vout(0:3,m)
      real * 8 betav,gamma
      real * 8 vdotb
      integer ipart,idim
      gamma=1/sqrt(1-beta**2)
      do ipart=1,m
         vdotb=vin(1,ipart)*vec(1)
     #         +vin(2,ipart)*vec(2)+vin(3,ipart)*vec(3)
         do idim=1,3
            vout(idim,ipart)=vin(idim,ipart)
     #           +vec(idim)*((gamma-1)*vdotb
     #           +gamma*beta*vin(0,ipart))
         enddo
         vout(0,ipart)=gamma*(vin(0,ipart)+vdotb*beta)
      enddo
      end subroutine

C------------------------------------------------------------
C     set_phep()
C     Set the momentum vectors phep(1:5,1:5) using generated phase space
C     This is needed to be able to call the powheg analysis, which uses phep(:,:)
      subroutine set_phep()
      use incl_parameters
      include 'hepevt.h'
      include 'pwhg_kn.h'
      double precision q1(1:4), q2(1:4)
      double precision p1q1, p2q1, p1q2, p2q2, psum, msq(1:5)
      integer i


      phep = 0d0
      phep(4,1:5) = kn_pborn(0,1:5) ! Energy component
      phep(1:3,1:5) = kn_pborn(1:3,1:5) ! Vector components
      msq(1:5) = phep(4,1:5)**2 - phep(1,1:5)**2 -
     $     phep(2,1:5)**2 - phep(3,1:5)**2
      ! for values that are slightly negative (eg. -10^(-14))
      ! round to zero to avoid error when taking the square root
      if (msq(1).lt.0d0) msq(1) = 0d0
      if (msq(2).lt.0d0) msq(2) = 0d0
      if (msq(4).lt.0d0) msq(4) = 0d0
      if (msq(5).lt.0d0) msq(5) = 0d0
      phep(5,1:5) = sqrt(msq(1:5)) 
C     assign particle ID (idhep) and information on in/out (isthep)      
C     This belongs in a common block....
      idhep(1:5) = 0
      idhep(3) = 25 ! = Higgs boson
      isthep(1) =-1 ! ingoing
      isthep(2) =-1 ! ingoing
      isthep(3) = 1 ! outgoing
      isthep(4) = 1 ! outgoing
      isthep(5) = 1 ! outgoing
      nhep=6
      return

      end subroutine

      end module phase_space
