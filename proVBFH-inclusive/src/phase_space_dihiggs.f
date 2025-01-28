C     Module phase_space is for the main part taken from VBFNLO
C     Usage:
C       1) Set up the beams with set_beams(sqrts). [only once]
C       2) Generate phase space with gen_phsp(xborn) using an array
C          of random numbers. This fills the internal variables
C          contained in pwhg_kn.h with new values
C       3) Map internal variables to the ones of the module with
C          set_phsp().
C       4) If using a powheg analysis, map module variables to
C          phep(:,:), the array in hepevt.h, used in powheg.
      module phase_space_dihiggs

      double precision x1, x2, vq1(0:3), vq2(0:3), pH1(0:3), pH2(0:3)
      double precision phi, jacobian
      double precision Q1_sq, q1_perp, q1_perp_sq, Q2_sq, q2_perp
      
      integer nlegborn, nlegreal
      parameter (nlegborn = 6) ! VBFHHMOD: Changed to 6 born legs
      parameter (nlegreal = nlegborn + 1)
      
      contains

  
C------------------------------------------------------------
C     dsigma function
C     calculate the cross section, using the phase space from VBFNLO
C     and matrix element computed with Hoppet
      double precision function dsigma(xrand, vegas_weight)
      use matrix_element_dihiggs
      use parameters
      implicit none
      include 'pwhg_kn.h'
C     xrand contains a vector of random numbers in [0,1]
      double precision xrand(11), vegas_weight
      double precision ptH1H2, dsig_temp,ds2
      integer vegas_ncall
      common/vegas_ncall/vegas_ncall

      dsigma = 0d0
C     generate phase space using phase_space module
      call gen_phsp(xrand(1:11))
C     set phase space variables x1, x2, vq1, vq2, phi, Q1_sq,
C     q1_perp, q1_perp_sq, Q2_sq, q2_perp using generated phase space
      call set_phsp()

C     skip phase space points with vanishing jacobian or with minQ1,Q2) < Qmin
      if ((jacobian.ne.0d0).and.(min(Q1_sq,Q2_sq).gt.(Qmin**2))) then
C     For scale_choice=3 we need ptH. We have in the transverse plane
C     q1(1:2) = (q1_perp         ,    0            )
C     q2(1:2) = (q2_perp*cos(phi), q2_perp*sin(phi))
C     pH1H2(:) = pH1(:) + pH2(:) = - q1(:) - q2(:)
C     Hence ptH1H2 = ptH1 + ptH2 is given by
         ptH1H2 = sqrt((q1_perp + q2_perp*cos(phi))**2
     $        + (q2_perp*sin(phi))**2)
C     compute dsigma using the squared hadronic tensor
C         dsigma = eval_matrix_element_vbfnlo(order_min,order_max, x1, x2,
C     $        kn_beams(:,1), kn_beams(:,2), vq1, vq2, pH1, pH2, ptH1H2)
C     print*, 'vbfnlo', dsigma
         if(tensorME) then
            dsigma = eval_matrix_element_tensor(order_min,order_max, x1, x2,
     $           kn_beams(:,1), kn_beams(:,2), vq1, vq2, pH1, pH2, ptH1H2)
!     print*, 'tensor', dsigma
         else
            dsigma = eval_matrix_element(order_min,order_max, x1, x2,
     $           kn_beams(:,1), kn_beams(:,2), vq1, vq2, pH1, pH2, ptH1H2)
         endif
!     print*, 'default', dsigma
!         stop
C     convert to [pb] and add in jacobian
         dsigma = dsigma * gev2pb * jacobian
      endif

      
C     map to powheg variables for use by analysis
      call set_phep()
C     remove the rare outliers where we get dsigma = NaN
      if (dsigma.ne.dsigma) dsigma = 0d0
      
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
      vq1(:) = kn_pborn(:,5) - kn_pborn(:,1)
      vq2(:) = kn_pborn(:,6) - kn_pborn(:,2)
      pH1(:) = kn_pborn(:,3)
      pH2(:) = kn_pborn(:,4)
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
C     gen_phsp(xborn(1:11))
C     generated phase space using an array of [0,1] random numbers
C     this code is taken from POWHEG-BOX/VBF_H
      subroutine gen_phsp(xborn)
      ! rename s, pi to avoid conflict with internal variables
      use parameters, sval => s, pival => pi 
      implicit none
      include 'pwhg_kn.h'
      include 'pwhg_math.h'
      include 'PhysPars.h'
      real * 8 xborn(11),xx(11)
      real * 8 m2,xjac,taumin,tau,y,beta,betaCM,vec(3),cth,s,
     #     z,zhigh,zlow,m
      integer mu,k,j
      logical ini
      data ini/.true./
      save ini
      real * 8 Vmass2,Vmass2low,Vmass2high,VmVw  
      real * 8 m2jj,pV(0:3),pVmod,pVmod2,pJ(0:3,2),cthj,phij,pcmjj(0:3),
     #     pcmjjmod,ptmp(0:3,2)
      real * 8 pH1(0:3), pH2(0:3), mH1, mH2, wt
      real * 8 pH13m(1:4), pH23m(1:4), pV3m(1:4)
      real * 8 mass      
      external mass      
      logical check
      parameter(check=.false.)
      real * 8 epsilon
      parameter (epsilon=1d-10)
      real * 8 m2jjmin
      real * 8 BW_fixed,BW_running
      real * 8 weight
      integer BWflag
      integer nlegborn, nlegreal
      parameter (nlegborn = 6)
      parameter (nlegreal = nlegborn + 1)
!     VBFNLO variables
      DOUBLE PRECISION RM2(0:2),RMG(0:2),RM2MIN(0:2),RM2MAX(0:2),ecm,pTjmin
      save RM2,RMG,RM2MIN,RM2MAX,ecm, pTjmin
      double precision q(0:4), k1(0:3), k2(0:3)
      double precision n
      parameter (n = 2d0)

      xx = xborn
!     Start by initialising everything we need
      kn_masses = 0d0
      kn_sborn = 0d0 ! Isn't needed for the inclusive part of the code..
      kn_xb1 = 1d0
      kn_xb2 = 1d0
      kn_pborn = 0d0
      kn_jacborn = 1d0
      q = 0d0
      if(ini) then
         rm2(0) = (2d0*mh)**2 ! Pseudo resonance with 2*Higgs mass
         rmg(0) = mh*500d0 ! Pseudo resonance width times mass (broad resonance)
         rm2min(0) = 1d-3       ! This value taken from vbfnlo. Should test
         rm2max(0) = kn_sbeams*0.5d0 ! Maximum allowed value for intermediate mass
         ecm = sqrt(kn_sbeams)  ! Centre-of-mass energy
         pTjmin = 0d0 ! minimum jet pt

         ! Currently we assume on-shell Higgs. Should be changed later for breit-wigner perhaps?
         rm2(1) = mh**2 ! Mass squared of first Higgs
         rm2(2) = mh**2 ! Mass squared of first Higgs
         ini = .false.
      endif
      
!     We borrow the phase space for HHjj from vbfnlo.  First generate
!     q^2 of intermediate pseudo-state.

!     AK: Understanding here is that we generate a fake resonance with a
!     large width to perform pseudo-decay into two on-shell (or
!     off-shell) Higgs. The routine returns kn_jacborn (0d0 if the
!     function is false) and the q^2 as q(4).
      kn_jacborn = 3d0*xx(1)**(3d0-1d0)*kn_jacborn
      xx(1) = 1d0 - xx(1)**3
      if (.not. Resonance(rm2(0), rmg(0), rm2min(0), rm2max(0),
     $     xx(1), kn_jacborn, q(4))) return

      kn_jacborn = 2d0*xx(2)**(2d0-1d0)*kn_jacborn
      xx(2) = 1d0 - xx(2)**2
      if(xx(4).lt.0.5d0) then
         xx(4) = 2d0*xx(4)
         xx(4) = xx(4)**n
         kn_jacborn =  kn_jacborn * n * xx(4)**(1d0-1d0/n)
         xx(4) = (1d0 - xx(4))*0.5d0
      else
         xx(4) = 2d0*(1d0 - xx(4))
         xx(4) = xx(4)**n
         kn_jacborn =  kn_jacborn * n * xx(4)**(1d0-1d0/n)
         xx(4) = (1d0 + xx(4))*0.5d0
      endif

      kn_jacborn = 2d0*xx(5)**(2d0-1d0)*kn_jacborn
      xx(5) = xx(5)**2
!     Now we generate the two-jet system from the pseudo-resonance q^2
      if (.not.TwoToJetsPlusX(2, xx(2:7), 0d0, ecm, pTjmin, q(4),
     $     kn_pborn(:,1), kn_pborn(:,2), kn_xb1, kn_xb2, q(0),
     $     kn_pborn(:,5:6), kn_jacborn)) return

      if(xx(8).lt.0.5d0) then
         xx(8) = 2d0*xx(8)
         xx(8) = xx(8)**n
         kn_jacborn =  kn_jacborn * n * xx(8)**(1d0-1d0/n)
         xx(8) = 0.5d0*xx(8)
      else
         xx(8) = 2d0*(1d0 - xx(8))
         xx(8) = xx(8)**n
         kn_jacborn =  kn_jacborn * n * xx(8)**(1d0-1d0/n)
         xx(8) = (2d0 - xx(8))*0.5d0
      endif
!     Now decay the pseudo-resonance
      if (.not. TwoBodyDecay(xx(8), xx(9), q(0), q(4), rm2(1),
     $     rm2(2), kn_pborn(:,3), kn_pborn(:,4), kn_jacborn)) return

      kn_masses(3:4) = sqrt(rm2(1:2))

      kn_jacborn = kn_jacborn * 0.5d0 ! symmetry factor from two identical particles in final state

!     print*, kn_xb1, kn_xb2, kn_jacborn

!      print*, 'Born momenta'
!      print*, '1', kn_pborn(:,1), kn_masses(1), invmass(kn_pborn(:,1))
!      print*, '2', kn_pborn(:,2), kn_masses(2), invmass(kn_pborn(:,2))
!      print*, '3', kn_pborn(:,3), kn_masses(3), invmass(kn_pborn(:,3))
!      print*, '4', kn_pborn(:,4), kn_masses(4), invmass(kn_pborn(:,4))
!      print*, '5', kn_pborn(:,5), kn_masses(5), invmass(kn_pborn(:,5))
!      print*, '6', kn_pborn(:,6), kn_masses(6), invmass(kn_pborn(:,6))
!     stop
      
c      boost in the CM frame 
      betaCM=(kn_xb1-kn_xb2)/(kn_xb1+kn_xb2)
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
C     Set the momentum vectors phep(1:6,1:5) using generated phase space
C     This is needed to be able to call the powheg analysis, which uses phep(:,:)
      subroutine set_phep()
      use parameters
      include 'hepevt.h'
      include 'pwhg_kn.h'
      double precision q1(1:4), q2(1:4)
      double precision p1q1, p2q1, p1q2, p2q2, psum, msq(1:nlegborn)
      integer i


      phep = 0d0
      phep(4,1:6) = kn_pborn(0,1:6) ! Energy component
      phep(1:3,1:6) = kn_pborn(1:3,1:6) ! Vector components
      msq(1:6) = phep(4,1:6)**2 - phep(1,1:6)**2 -
     $     phep(2,1:6)**2 - phep(3,1:6)**2
      ! for values that are slightly negative (eg. -10^(-14))
      ! round to zero to avoid error when taking the square root
      if (msq(1).lt.0d0) msq(1) = 0d0
      if (msq(2).lt.0d0) msq(2) = 0d0
      if (msq(3).lt.0d0) msq(3) = 0d0
      if (msq(4).lt.0d0) msq(4) = 0d0
      if (msq(5).lt.0d0) msq(5) = 0d0
      if (msq(6).lt.0d0) msq(6) = 0d0
      phep(5,1:6) = sqrt(msq(1:6)) 
!      write(*,'(5(f13.5))') phep(:,1)
!      write(*,'(5(f13.5))') phep(:,2)
!      write(*,'(5(f13.5))') phep(:,3)
!      write(*,'(5(f13.5))') phep(:,4)
!      write(*,'(5(f13.5))') phep(:,5)
!      write(*,'(5(f13.5))') phep(:,6)
!      print*, ''
C     assign particle ID (idhep) and information on in/out (isthep)      
C     This belongs in a common block....
      idhep(1:6) = 0
      idhep(3) = 25 ! = Higgs boson
      idhep(4) = 25 ! = Higgs boson
      isthep(1) =-1 ! ingoing
      isthep(2) =-1 ! ingoing
      isthep(3) = 1 ! outgoing
      isthep(4) = 1 ! outgoing
      isthep(5) = 1 ! outgoing
      isthep(6) = 1 ! outgoing
      nhep=nlegborn
      return

      end subroutine
      
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
      
C********** BOOSTN ****************************************************
C
      SUBROUTINE BOOSTN(P,R,Q)
C
C
C     The four vector P is assumed to be given in the rest frame of R,
C     which must be a timelike vector.
C     output Q is the vector P boosted to the frame in which R is given.
C                                              Compare Jackson, p.517
C                                              D. Zeppenfeld (28.6.1985)
C     New version with energy stored in zeroth component of four vector
C     arrays. Checked on May 27, 1988.
C
      REAL*8 P(0:3),R(0:3),Q(0:3)
      REAL*8 BETA(3), X, Y, GAMMA
      INTEGER I
      X = 0D0
      Y = 0D0
      DO I = 1,3
         BETA(I) = R(I)/R(0)
         X = X + BETA(I)**2
         Y = Y + BETA(I)*P(I)
      ENDDO
      IF (X.GT.1D-16.AND.X.LT.(1D0-1D-12)) THEN
         GAMMA = 1D0/DSQRT(1D0-X)
         DO I = 1,3
            Q(I) = P(I)+BETA(I)*(Y*(GAMMA-1D0)/X + GAMMA*P(0))
         ENDDO
         Q(0) = GAMMA*(P(0) + Y)
      ELSE
         DO I = 0,3
            Q(I) = P(I)
         ENDDO
         IF(X.GE.(1D0-1D-12)) 
     *      WRITE(6,1000) R,R(0)**2-R(1)**2-R(2)**2-R(3)**2
      ENDIF
 1000 FORMAT (" The reference vector ",4G12.3," is not timelike."/
     1        " R**2 = ",G12.3)
      END
      
c***************************************************************************
      function dotrr(p1,p2)
c***************************************************************************
c
c     dotrr(p1,p2) = p1.p2
c
c***************************************************************************
      implicit none

      double precision dotrr,p1(0:3),p2(0:3)

      dotrr = p1(0)*p2(0) - p1(1)*p2(1) - p1(2)*p2(2) - p1(3)*p2(3)

      end

      end module phase_space_dihiggs
