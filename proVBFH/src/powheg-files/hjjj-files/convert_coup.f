c*************************************************************************
      SUBROUTINE COUP_POWHEG_TO_VBFNLO
c*************************************************************************
c     This is the interface to the koppln-routine for the initialization
c     of coupling constants. 
c*************************************************************************
      implicit none
      include "koppln_ew.inc"
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      include 'global.inc'

** Used for determination of slha input file
      character*100 GetInputPath
      external GetInputPath
      character*100 path

** Function which determines whether a parameter has already been set with a 
** reasonable value
      logical replace
      external replace

** 1/Alfa_QED: this is the default read-in value, but is immediately converted
** to alfa
      double precision invalfa
      double precision vbfnloinput, pwhg_alphas
      external vbfnloinput, pwhg_alphas
** Coupling parameters
      double precision e, g2, s, c, z, w, q, g

** Branching ratios and widths
      double precision BWNE,BWUD,BWTB,BZNN,BZEE,BZUU,BZDD,BZTT,
     &                BHWW,BHZZ,BHGG,BHTT,BHBB,BHCC,BHTAU,BHMU,
     &                BHGAM, BHGAMZ, XGW, XGZ, XGH, GAMT
      COMMON /BRANCH/ BWNE,BWUD,BWTB,BZNN,BZEE,BZUU,BZDD,BZTT,
     &                BHWW,BHZZ,BHGG,BHTT,BHBB,BHCC,BHTAU,BHMU,
     &                BHGAM, BHGAMZ, XGW, XGZ, XGH, GAMT


      call clearwidths
      model=1
      higgstype=0

c strong coupling from powheg:
!       call setscalesbtilde
      print *," "
      print *,"              Physics parameters"
      print *,"-----------------------------------------------"
      xmh=vbfnloinput('#HMASS')!call read_Real("HMASS",xmh,120d0)
      xgh=-999d0
      XGH=vbfnloinput('#HWIDTH')
      if (replace(xmt,0)) xmt=vbfnloinput('#TOPMASS')!call read_ReaL("TOPMASS",xmt,172.4d0)
      if (replace(xmb,0)) xmb= vbfnloinput("#BOTTOMMASS")
      if (replace(xmc,0)) xmc=vbfnloinput("#CHARMMASS")
      if (replace(xmtau,0)) xmtau=vbfnloinput("#TAU_MASS")
      

c     ewscheme=int(vbfnloinput("#EWSCHEME"))
c only scheme=3 implemented consistently:
      ewscheme = 3
      
      Select Case (EWSCHEME)
      Case(1)
         if (replace(gf,0)) gf= vbfnloinput("#FERMI_CONST")
         if (replace(alfa,0)) then 
            alfa= vbfnloinput("#INVALFA")
            alfa=1/alfa
         endif
         if (replace(xmz,0)) xmz= vbfnloinput("#ZMASS")
         xmw = -1.d0
         sin2w = -1.d0
      CASE(2)
         if (replace(gf,0)) gf= vbfnloinput("#FERMI_CONST")
         if (replace(sin2w,0)) sin2w= vbfnloinput("#SIN2W")
         if (replace(xmz,0)) xmz= vbfnloinput("#ZMASS")
         xmw = -1.d0
         alfa = -1.d0
      CASE(3)
         if (replace(gf,0)) gf= vbfnloinput("FERMI_CONST")
         if (replace(xmw,0)) xmw= vbfnloinput("WMASS")
         if (replace(xmz,0)) xmz= vbfnloinput("ZMASS")
         alfa = -1.d0
         sin2w = -1.d0
      CASE(4)
         if (replace(gf,0)) gf= vbfnloinput("FERMI_CONST")
         if (replace(alfa,0)) then 
            alfa= vbfnloinput("#INVALFA")
            alfa=1/alfa
         endif
         sin2w= vbfnloinput("SIN2W")
         if (replace(xmw,0)) xmw= vbfnloinput("WMASS")
         if (replace(xmz,0)) xmz= vbfnloinput("ZMASS")
      Case Default
         print*,"Invalid choice for EWSCHEME"
      End Select

      alfas= pwhg_alphas(xmz**2,st_lambda5MSB,st_nlight)

** Setting the electroweak parameters
      call setEWpara(e,g2,s,c,z,w,q,g)

      CALL KOPPLN(2,e,g2,s,c,z,w,q,g)
      call ctrans(xmb)

      RETURN
      END 
! 
! 
C========= subroutine KOPPLN  ==============================================
C
C  Adapted from RKWW by R.Kleiss
C
C  H --> gamma gamma decay added by D. Rainwater,
C        from HHG, Mar 97
c  #included D5 operators by TFigy, June 04
C  H --> b bbar decay modified by D. Rainwater,
C        from HHG, Mar 98
C
c  H --> V*V* -> all added by T.Figy, Oct 2003
c
C  Calculate the fermion-fermion-boson and three-boson couplings
C  in the standard model, without colour factors.
C  also calculate mass**2 and mass*width of the bosons.
C  V   = FFB vector coupling
C  A   = FFB axial-vector coupling
C  CLR = FFB left- and righthanded couplings
C  B   = BBB coupling
C  The fermion indices are:  1: neutrino of electron (muon,tau)
C                            2: electron (muon,tau)
C                            3: up quark (charm,top)
C                            4: down quark (strange,bottom)
C  The boson indices are     1: photon
C                            2: Z0 boson
C                            3: W+ boson
C                            4: W- boson
C                            5: gluon
C                            6: Higgs
C  The helicity indices are  -1: lefthanded
C                            +1: righthanded
C  No Kobayashi-Maskawa mixing is implemented here.
C  The W+- are identified as outgoing from the BBB vertex.
C  XM2 = boson mass**2
C  XMG = boson mass*width
C  Also the branching ratios of W and Z decay are calculated in the SM
C
C  The fundamental parameters:
C  ALFA = The Q.E.D. fine structure constant
C  SIN2W = sin**2 of the weak mixing angle
C  ALFAS = the QCD fine structure constant
C  XMZ = mass of the Z0 boson in GeV
C  XMW = mass of the W+- bosons in GeV
C  XMH = mass of the Higgs boson in GeV
C  XMT = mass of the top quark in GeV
C
C  printout is steered by INFO:
C   <=0 :   no printout
C   <=1 :   printout of couplings only
C    =2 :   printout of input parameters and branching ratios
C   else:   printout of couplings and W/Z branching ratios
C
C===========================================================================
C

*** IMPORTANT NOTE for VBFNLO!
*** In order to avoid hard-wiring values of the electroweak parameters when
*** using anomalous couplings, the setting of these parameters has been moved to
*** a new subroutine, 'setEWpara', which can be found in utilities/parameters.F


      SUBROUTINE KOPPLN(INFO,e,g2,s,c,z,w,q,g)

      implicit none
      double complex  f_tau_t, f_tau_b, f_tau_W, F_t, F_b, F_W, dum
      double precision  beta0, beta1, gamma0, gamma1
      double precision  mllb
      double precision  zero
      parameter ( zero = 0.0d0 )
      integer  info, i, j, k
      logical ldebug
c
c INPUT
      include 'pwhg_st.h'
      include "koppln_ew.inc"
      include 'pwhg_math.h'
c
c OUTPUT
c
      double precision XM2s(6),XMGs(6),
     & BWNE,BWUD,BWTB,BZNN,BZEE,BZUU,BZDD,BZTT,
     & BHWW,BHZZ,BHGG,BHTT,BHBB,BHCC,BHTAU,BHMU,BHGAM
      COMMON /BKOPOUshort/ XM2s,XMGs
 

** clrCT is the counterterm for qqV
      double complex clrCT(3:4,2:4,-1:1)
      common /bkopouCT/ clrCT
      DOUBLE PRECISION CLR,XM2,XMG,B,V,A
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),V(4,5),A(4,5)

c     ahvv(i,v1,v2) i = 1,2,3 are coefficients of eq 1 of 
c     Physics Letters B 591,297
      double complex ahvv(3,4,4), ahvvL(3,4,4)
      common/tensorhvv/ ahvv, ahvvL
      real*8 treefacW, treefacZ, loopfac 
      common/lhcoup/ treefacW, treefacZ, loopfac
      double precision g5hvv(2,1:4,1:4),lambda5
      common/hcoupl/g5hvv,lambda5
c     form factor for hvv couplings
      integer ffac
      double precision m2ff, mff 
      logical lff
      common/ formfacmass/ m2ff,mff,ffac,lff
C
C local variable
C

** coupling parameters, widths, ..
      double precision e,g,s,c,g2,z,w,q, 
     & gwh, gzh, decw, gwne, gwud, gwtb, decz, gznn, gzee, gzuu, gzdd,
     & gztt, xgz, fh, y1,y2,x,xgw,xlw,xlz,dhzz,dhwwst,
     & ghww,ghzz,ghtt, betac,betab,betat, qcdc,qcdb,qcdt,
     & ghcc,ghbb,ghtau,ghgg,ghgam,tau_t,eta_p,eta_m,tau_b,tau_w,xgh
      double precision xmMU,betaMU,ghmu,gh4f,gamt,alphas5,alphas
      external alphas5, alphas
      double precision BHGAMZ,GHGAMZ
      COMMON /BRANCH/ BWNE,BWUD,BWTB,BZNN,BZEE,BZUU,BZDD,BZTT,
     &                BHWW,BHZZ,BHGG,BHTT,BHBB,BHCC,BHTAU,BHMU,
     &                BHGAM, BHGAMZ, XGW, XGZ, XGH, GAMT
      double precision xmax,xmin
      double complex II_1,II_2,A_sum,A_w,A_t,A_b
      double precision cos2w,tan2w,cotw,lamb_b,lamb_t,lamb_w

c parameters needed for calculation of H -> B BBAR
      double precision alphas_H, alphas_B, BMASS_POLE, BMASS_B
      double precision alphas_C, CMASS_POLE, CMASS_C
      double precision DELTA_qq, DELTA_H, C1, C2
      double precision A_12_t, A_12_b      
      integer Nf
      logical replace
      double precision pwhg_alphas, heavy
      external replace, pwhg_alphas, heavy
c
c mass dependence for H--->W*W,Z*Z (From Gunion et al., Higgs Hunters Guide)
c
c      FH(X) = - (1.d0 - X**2)*(23.5d0*X**2 - 6.5d0 + 1.d0/X**2) 
c     &        - 3.d0*(1.d0 - 6.d0*X**2 + 4.d0*X**4)*LOG(X)
c     &        + 3.d0*(1.d0 - 8.d0*X**2 + 20.d0*X**4)/
c     &                                   SQRT(4.d0*X**2 - 1.d0)
c     &        * ACOS( (3.d0*X**2 - 1.d0)/(2.d0*X**3) )
c      PI = 4.d0*ATAN(1.d0)

      ldebug = .false.   ! alfa.eq.0.122431


c     initialize ahvv to zero
      A_sum=(0.0d0,0.0d0)
      dum =(0.0d0,0.0d0)
         do i = 1,4
            do j = 1,4
               ahvv(1,i,j) = (0.0d0,0.0d0)
               ahvv(2,i,j) = (0.0d0,0.0d0)
               ahvv(3,i,j) = (0.0d0,0.0d0)
               g5hvv(2,i,j) =0.0d0
               g5hvv(1,i,j) =0.0d0
            enddo
         enddo
         treefacW = 1.0d0
         treefacZ = 1.0d0
         loopfac = 1.0d0
         lambda5 = 1D0
         lff=.false.


c THE FFB VECTOR COUPLING
c
* photon-fermion coupling
      V(1,1) = zero
      V(2,1) = -E
      V(3,1) = 2.d0*E/3.d0
      V(4,1) = -E/3.d0
* Z-fermion coupling
      V(1,2) =  Z
      V(2,2) = -Z*(1.d0 - 4.d0*SIN2W)
      V(3,2) =  Z*(1.d0 - 8.d0/3.d0*SIN2W)
      V(4,2) = -Z*(1.d0 - 4.d0/3.d0*SIN2W)

* W-fermion coupling
      DO 1 K = 1,4
         DO 1 J = 3,4
            V(K,J) = W
 1    continue
* gluon-fermion coupling
      V(1,5) = zero
      V(2,5) = zero
      V(3,5) = G
      V(4,5) = G

c
c THE FFB AXIAL-VECTOR COUPLING
c
* photon/gluon - fermion coupling
      DO 2 K = 1,4
         DO 2 J = 1,5,4
            A(K,J) = zero
 2    continue
* Z-fermion coupling
      A(1,2) = -Z
      A(2,2) =  Z
      A(3,2) = -Z
      A(4,2) =  Z
* W-fermion coupling
      do K = 1,4
         do J = 3,4
            A(K,J) = -W
         enddo
      enddo

* switching sin-theta-w conventions for the mssm
      if (model .eq. 2) then
         do k = 1, 4
            do j = 2, 4
               V(k,j) = -V(k,j)
               A(k,j) = -A(k,j)
            end do
         end do
      end if

c
c THE FFB LEFT- AND RIGHTHANDED COUPLINGS
c
      DO K = 1,4
         DO J = 1,5
            CLR(K,J,0) = zero
            CLR(K,J,-1)= V(K,J)-A(K,J)
            CLR(K,J,1) = V(K,J)+A(K,J)
         end do
      end do
      if (ldebug) then
         print*," electron L coupling "
         print*, (clr(2,i,-1),i=1,5)
         print*," electron R coupling "
         print*, (clr(2,i,+1),i=1,5)
      endif

c
c THE BBB COUPLINGS
c
      DO I = 1,6
         DO J = 1,6
            DO K = 1,6
               B(I,J,K) = zero
            end do
         end do
      end do
* W-W-photon coupling
      B(3,4,1) = -E
      B(4,1,3) = -E
      B(1,3,4) = -E
      B(1,4,3) =  E
      B(4,3,1) =  E
      B(3,1,4) =  E
* W-W-Z coupling
      B(3,4,2) = -Q
      B(4,2,3) = -Q
      B(2,3,4) = -Q
      B(2,4,3) =  Q
      B(4,3,2) =  Q
      B(3,2,4) =  Q
* altering sin-theta-w conventions for mssm
      if (model .eq. 2) then
         B(3,4,2) = Q
         B(4,2,3) = Q
         B(2,3,4) = Q
         B(2,4,3) = -Q
         B(4,3,2) = -Q
         B(3,2,4) = -Q
      end if         
* gluon coupling
      B(5,5,5) =  G
c
c Higgs couplings to W and Z
c
      GWH = G2
      GZH = GWH/C**2

* Z-Z-Higgs coupling
      B(2,2,6) = GZH
      B(2,6,2) = GZH
      B(6,2,2) = GZH
* W-W-Higgs coupling
      B(3,4,6) = GWH
      B(4,3,6) = GWH
      B(3,6,4) = GWH
      B(4,6,3) = GWH
      B(6,3,4) = GWH
      B(6,4,3) = GWH
      
* a1 are sm hzz and hww couplings
      ahvv(1,2,2) = treefacZ*B(6,2,2)*xmw 
      do i=3,4
         do j=3,4
            ahvv(1,i,j) = treefacW*B(6,i,j)*xmw 
         enddo
      enddo
c
c Calculate the total decay widths and branching ratios of the heavy bosons if 
c we're not using a slha file, or if they're not included in the slha file
c
c W first
c
      DECW = XMW/12.d0/PI

      if (replace(BWNE,0)) then
         GWNE = DECW*( V(1,3)**2 + A(1,3)**2 )
      else
         GWNE = BWNE*XGW
      end if

      if (replace(BWUD,0)) then
         GWUD = DECW*( V(3,3)**2 + A(3,3)**2 )*3.d0*(1.d0 + ALFAS/PI)
      else
         GWUD = BWUD*XGW
      end if

      if (replace(BWTB,0)) then
         GWTB = zero
         IF(XMT.GT.XMW) GOTO 31
         GWTB = GWUD*HEAVY( V(3,3)/A(3,3), XMT/XMW, zero )
      else
         GWTB = BWTB*XGW
      end if

 31   if (replace(XGW,0)) then
         XGW  = 3.d0*GWNE + 2.d0*GWUD + GWTB
      end if

      BWNE = GWNE/XGW
      BWUD = GWUD/XGW
      BWTB = GWTB/XGW
c
c Z next
c
      DECZ = XMZ/12.d0/PI

      if (replace(BZNN,0)) then
         GZNN = DECZ*( V(1,2)**2 + A(1,2)**2 )
      else
         GZNN = BZNN*XGZ
      end if

      if (replace(BZEE,0)) then
         GZEE = DECZ*( V(2,2)**2 + A(2,2)**2 )
      else
         GZEE = BZEE*XGZ
      end if

      if (replace(BZUU,0)) then
         GZUU = DECZ*( V(3,2)**2 + A(3,2)**2 )*3.d0*(1.d0 + ALFAS/PI)
      else
         GZUU = BZUU*XGZ
      end if

      if (replace(BZDD,0)) then
         GZDD = DECZ*( V(4,2)**2 + A(4,2)**2 )*3.d0*(1.d0 + ALFAS/PI)
      else
         GZDD = BZDD*XGZ
      end if

      if (replace(BZTT,0)) then
         GZTT = zero
         IF(XMT.GT.XMZ/2.d0) GOTO 32
         GZTT = GZUU*HEAVY( V(3,2)/A(3,2) , XMT/XMZ , XMT/XMZ )
      else
         GZTT = BZTT*XGZ
      end if

 32   if (replace(XGZ,0)) then
         XGZ  = 3.d0*GZNN + 3.d0*GZEE + 2.d0*GZUU + 3.d0*GZDD + GZTT
      end if
      BZNN = GZNN/XGZ
      BZEE = GZEE/XGZ
      BZUU = GZUU/XGZ
      BZDD = GZDD/XGZ
      BZTT = GZTT/XGZ

c
c     mass and width info 
c
      XM2(1) = zero
      XM2(2) = XMZ**2
      XM2(3) = XMW**2
      XM2(4) = XM2(3)
      XM2(5) = zero
      xm2(6) = xmh**2
c
      XMG(1) = zero
      XMG(2) = XMZ*XGZ
      XMG(3) = XMW*XGW
      XMG(4) = XMG(3)
      XMG(5) = zero
c     cannot compute xmg(6) at this point
c
c dummy xm2 and xmg for common block BKOPOUshort 
c
      do i = 1,6
         xm2s(i) = xm2(i)
      end do
      do i = 1,5
         xmgs(i) = xmg(i)
      end do     


c******************************************************************
c Next the Higgs decays to WW, ZZ, bb, cc, tau tau and gamma gamma.
c These are either taken from the slha file, from FH, or are 
c calculated below
c******************************************************************

** Kill the code if we have a very low higgs mass and are not working in the
** complex MSSM
c      IF ((XMH.LT.60.d0) .and. (higgsmix .ne. 3)) THEN
c         XMH = 60.d0
c         print*," Higgs mass below 60 GeV is not valid: m_H = ",xmh
c      ENDIF

      XLZ = (XMZ/XMH)**2
      XLW = (XMW/XMH)**2


c     Higgs decays to WW and ZZ
c     subroutines below use  BKOPOUshort common block
c **********  First h -> Z*Z* **********************
      if (replace(bhzz,3) .or. replace(xgh,3)) then

         call flavor_sum(ghzz)
         if (replace(xgh,3) .and. (.not. replace(bhzz,3))) then
            xgh = bhzz*ghzz
         end if
      else
         ghzz = bhzz*xgh
      end if


c*********** Next h -> W*W* *************************
      if (replace(bhww,3) .or. replace(xgh,3)) then
         call initbosid(4)
         call calZ(0.0d0,xmin)      
         call calZ(xmh**2,xmax)
         if (ldebug) print*," xmin,xmax = ",xmin,xmax
        
         call quad2d(xmin,xmax,ghww)
         if(ghww.lt.0.0d0) ghww = 0.0d0
         if (replace(xgh,3) .and. (.not. replace(bhww,3))) then
            xgh = bhww*ghww
         end if
      else 
         ghww = bhww*xgh
      end if


** Calculating running bottom mass
      BMASS_POLE = xmb

c the coefficient fuction c is evaluated within the five flavour 
c approximation for the whole mass range. For m_Higgs > m_top,
c the deviation to the six flavour scheme is less the 1%
      Nf=5

      alphas_B=pwhg_alphas(BMASS_POLE**2,st_lambda5MSB,-1)!alphas5(BMASS_POLE**2,1)
      alphas_H=pwhg_alphas(xmh**2,st_lambda5MSB,st_nlight)!alphas5(xmh**2,1)


c calculate the running MSbar bottom mass mb(m_b) = BMASS_B
      BMASS_B=BMASS_POLE*(1.0d0-4.d0/3.d0*alphas_B/PI 
     &        +(1.0414d0*Nf-14.3323d0)*alphas_B**2/PI**2)

c calculate the running MSbar bottom mass mb(m_H) = BMASS_H
      alphas_B = pwhg_alphas(BMASS_B**2,st_lambda5MSB,-1)!alphas5(BMASS_B**2,1)

      C1=(23.0d0/6.0d0*alphas_H/PI)**(12.d0/23.d0)*(1.d0+1.17549d0
     &     *alphas_H/PI+1.50071d0*alphas_H**2/PI**2+0.172478d0
     &     *alphas_H**3/PI**3)
      C2=(23d0/6.d0*alphas_B/PI)**(12.d0/23d0)*(1.d0+1.17549d0
     &     *alphas_B/PI+1.50071d0*alphas_B**2/PI**2+0.172478d0
     &     *alphas_B**3/PI**3)

      BMASS_H=BMASS_B*C1/C2


c*********** Next h -> b bbar ************************
      if (replace(bhbb,3) .or. replace(xgh,3)) then
c Now calculate the H -> bbar width
         DELTA_qq = 5.67d0*alphas_H/PI + (35.94-1.36*Nf)*
     &        alphas_H**2/PI**2 + (164.14-25.77*Nf+0.26*Nf**2)*
     &        alphas_H**3/PI**3

         DELTA_H = alphas_H**2/PI**2 *(1.57-2.0d0/3.0d0*
     -        LOG(xmh**2/xmt**2) + 1.0d0/9.0d0*
     -        (LOG(BMASS_H**2/xmh**2))**2)

         GHBB = 3.0d0*GF/(4.0d0*PI*SQRT(2.0d0)) * XMH * BMASS_H**2
     -        *( 1 + DELTA_qq + DELTA_H )
         if (replace(xgh,3) .and. (.not. replace(bhbb,3))) then
            xgh = bhbb*ghbb
         end if
      else 
         ghbb = bhbb*xgh
      end if  ! end of Higgs to b bar width calculation


c*********** Next h -> c cbar ************************
      if (replace(bhcc,3) .or. replace(xgh,3)) then
         CMASS_POLE = xmc
c determine number of flavors with mass smaller than the Higgs mass
c and the running alphas according to that number
         if (xmh.gt.xmt) then
            Nf = 6
         else
            Nf = 5
         endif
         alphas_H = pwhg_alphas(xmh**2,st_lambda5MSB,-1)!alphas5(xmh**2,1)
         alphas_C = pwhg_alphas(cmass_pole**2,st_lambda5MSB,-1)!alphas5(CMASS_POLE**2,1)

c calculate the running MSbar charm mass mc(m_c) = CMASS_C
         CMASS_C = CMASS_POLE*(1.0d0-4.d0/3.d0*alphas_C/PI 
     1        + (1.0414*Nf-14.3323)*alphas_C**2/PI**2)

c calculate the running MSbar charm mass mc(m_H) = CMASS_H
         C1=(23.0d0/6.0d0*alphas_H/PI)**(12.d0/23.d0)*(1.d0+1.175d0*
     -        alphas_H/PI+1.501d0*alphas_H**2/PI**2+0.1725d0*
     &        alphas_H**3/PI**3)
         C2=(23d0/6.d0*alphas_C/PI)**(12.d0/23d0)*(1.d0+1.175d0*
     &        alphas_C/PI+1.501d0*alphas_C**2/PI**2+0.1725d0*
     &        alphas_C**3/PI**3)

         CMASS_H = CMASS_C*C1/C2

c Now calculate the H -> ccbar width
         DELTA_qq = 5.67d0*alphas_H/PI + (35.94-1.36*Nf)*
     &        alphas_H**2/PI**2 + (164.14-25.77*Nf+0.26*Nf**2)*
     &        alphas_H**3/PI**3

         DELTA_H = alphas_H**2/PI**2 *(1.57-2.0d0/3.0d0*
     &        LOG(xmh**2/xmt**2)+1.0d0/9.0d0*
     &        (LOG(CMASS_H**2/xmh**2))**2)

         GHCC = 3.0d0*GF/(4.0d0*PI*SQRT(2.0d0)) * XMH * CMASS_H**2
     &        *( 1 + DELTA_qq + DELTA_H )
         if (replace(xgh,3) .and. (.not. replace(bhcc,3))) then
            xgh = bhcc*ghcc
         end if
      else 
         ghcc = bhcc*xgh
      end if  ! end of higgs to c cbar width calculation


c*********** Next h -> t tbar ************************
      if (replace(bhtt,3) .or. replace(xgh,3)) then
         IF ( XMH.GT.2.d0*XMT ) THEN
            BETAT = SQRT(1.d0 - 4.d0*(XMT/XMH)**2)
            beta0 = 7.d0        !11.d0 - 12.d0/3.d0
            QCDT = ( LOG(2.d0*XMT/0.2d0) / 
     &           LOG(XMH/0.2d0) )**(8.d0/beta0)*( 1.d0 + 3.d0*ALFAS/PI )
            GHTT = 3.d0*GF*XMT**2/(4.d0*PI*SQRT(2.d0)) * BETAT**3 * 
     &           XMH * QCDT
         ELSE
            GHTT = zero
         ENDIF
         if (replace(xgh,3) .and. (.not. replace(bhtt,3))) then
            xgh = bhtt*ghtt
         end if
      else 
         ghtt = bhtt*xgh
      end if  ! end of higgs to t tbar 


** Higgs to tau-plus tau-minus decay
      if (replace(bhtau,0) .or. replace(xgh,3)) then
*         XMtau = 1.77684d0      !  and  H ----> tau+ tau-
*         XMtau = ML
         BETAt = SQRT(1.d0 - 4.d0*(XMtau/XMH)**2)
         GHtau = GF*XMtau**2/(4.d0*PI*SQRT(2.d0)) * BETAt**3 * XMH 
         if (replace(xgh,3) .and. (.not. replace(bhtau,3))) then
            xgh = bhtau*ghtau
         end if
      else 
         ghtau = bhtau*xgh
      end if


** Higgs to gg 
** Leading-order is exact result -- see e.g. Djouadi, hep-ph/0503172
      tau_t = xmh**2/(4d0*xmt**2)
      if (tau_t .le. 1d0) then
        A_12_t = asin(sqrt(tau_t))**2
      else
        A_12_t = -( 
     &    log( (1+sqrt(1-1/tau_t))/(1-sqrt(1-1/tau_t)) ) - (0d0,pi)
     &    )**2/4d0
      endif
      A_12_t = 2*(tau_t+(tau_t-1)*A_12_t)/tau_t**2

      tau_b = xmh**2/(4d0*xmb**2)
      if (tau_b .le. 1d0) then
        A_12_b = asin(tau_b)**2
      else
        A_12_b = -( 
     &    log( (1+sqrt(1-1/tau_b))/(1-sqrt(1-1/tau_b)) ) - (0d0,pi)
     &    )**2/4d0
      endif
      A_12_b = 2*(tau_b+(tau_b-1)*A_12_b)/tau_b**2

      GHgg = GF*alphas_H**2*XMH**3/(36d0*sqrt(2d0)*pi**3)
      GHgg = GHgg*abs(3d0/4d0*(A_12_t+A_12_b))**2

      if (tau_t .le. 1d0) then
c higher-order corrections -- see eg. Schreck, Steinhauser, arXiv:0708.0916
      GHgg = GHgg * (1d0 + 
     &  alphas_H/pi*(17.9167d0+9.574d0*tau_t+5.571d0*tau_t**2
     &              +3.533d0*tau_t**3+2.395d0*tau_t**4+1.708d0*tau_t**4)
     & +(alphas_H/pi)**2*(156.808d0+109.365d0*tau_t+74.434d0*tau_t**2)
     & +(alphas_H/pi)**3*(467.684d0)
     & )
      endif


c for  H --> gamma gamma  I use p. 25 of Higgs Hunter"s Guide,
c using t,b and W loops only
c Physics Letters B 318,1 eq 12 with cpodd d5 operator
c      GHGAM = (alfa*G2/xmw)**2 * (xmh/pi)**3 / 1024.d0
      if (replace(bhgam,3) .or. replace(xgh,3)) then
         ghgam = alfa*g2/(4.0d0*pi*xmw)
         tau_t = 4.d0*(xmt/xmh)**2
         if ( tau_t.lt.1.d0 ) then
            eta_p = 1.d0 + sqrt(1.d0 - tau_t)
            eta_m = 1.d0 - sqrt(1.d0 - tau_t)
            f_tau_t = -0.25d0 * ( log(eta_p/eta_m) + pi*(0d0,-1.d0) )**2
         else
            f_tau_t = ( asin(sqrt(1.d0/tau_t)) )**2
         end if
         F_t = -2.d0*tau_t* (1.d0 + (1.d0 - tau_t)*f_tau_t)

         tau_b = 4.d0*(xmb/xmh)**2
         if ( tau_b.lt.1.d0 ) then
            eta_p = 1.d0 + sqrt(1.d0 - tau_b)
            eta_m = 1.d0 - sqrt(1.d0 - tau_b)
            f_tau_b = -0.25d0 * ( log(eta_p/eta_m) + pi*(0d0,-1.d0) )**2
         else
            f_tau_b = ( asin(sqrt(1.d0/tau_b)) )**2
         end if
         F_b = -2.d0*tau_b* (1.d0 + (1.d0 - tau_b)*f_tau_b)
         
         tau_W = 4.d0*(xmw/xmh)**2
         if ( tau_W.lt.1.d0 ) then
            eta_p = 1.d0 + sqrt(1.d0 - tau_W)
            eta_m = 1.d0 - sqrt(1.d0 - tau_W)
            f_tau_W = -0.25d0 * ( log(eta_p/eta_m) + pi*(0d0,-1.d0) )**2
         else
            f_tau_W = ( asin(sqrt(1.d0/tau_W)) )**2
         end if
         F_W = 2.d0 + 3.d0*tau_W*(1.d0 + (2.d0 - tau_W)*f_tau_W)
      
         dum = 3.d0*(4.d0*F_t + F_b)/9.d0 + F_W
         
         dum = ghgam*dum
         ahvv(2,1,1) = dum

         ghgam = xmh**3/(64.0d0*pi) * 
     &        ( dble(ahvv(2,1,1))**2+dimag(ahvv(2,1,1))**2+
     &        dble(ahvv(3,1,1))**2+dimag(ahvv(3,1,1))**2 )
c      GHGAM = GHGAM * ( dble(dum)**2 + dimag(dum)**2 ) 
         if (replace(xgh,3) .and. (.not. replace(bhgam,3))) then
            xgh = bhgam*ghgam
         end if
      else if (xgh .gt. 0d0) then
         ghgam = bhgam*xgh
      end if


c for H -> gamma Z (on-shell)
c cp even and odd d5 operators
c     Using Higgs Hunters guide
c     t,b ,W loops are include in SM calculation
c     compute At,Ab,Aw
      if (replace(bhgamz,3) .or. replace(xgh,3)) then
         cos2w = 1.0d0 - sin2w
         tan2w = sin2w/cos2w
         cotw = sqrt(1.0d0/tan2w)
c     
c     b quark loop
         tau_b = 4.d0*(xmb/xmh)**2
         lamb_b = 4.0d0*(xmb/xmz)**2
         call compints(tau_b,lamb_b,II_1,II_2)
c   
         A_b = -3.0d0 *2.0d0*V(4,1)*V(4,2)*(2.0d0*c/(g2*E))*(II_1 -II_2)
         A_b = A_b/(s*c) 
c     top quark loop
         tau_t = 4.d0*(xmt/xmh)**2
         lamb_t = 4.0d0*(xmt/xmz)**2
         call compints(tau_t,lamb_t,II_1,II_2)
c   
         A_t = -3.0d0 *2.0d0*V(3,1)*V(3,2)*(2.0d0*c/(g2*E))*(II_1 -II_2)
         A_t = A_t/(s*c) 
c     W loop
         tau_w = 4.0d0*(xmw/xmh)**2
         lamb_w = 4.0d0*(xmw/xmz)**2
         call compints(tau_w,lamb_w,II_1,II_2)
         A_w = 4.0d0*(3.0d0-tan2w)*II_2
         A_w = A_w + ((1.0d0+2.0d0/tau_w)*tan2w
     $        -(5.0d0+2.0d0/tau_w))*II_1
         A_w = -cotw*A_w
         A_sum = A_w + A_t + A_b
         A_sum = A_sum * alfa*g2/(4.0d0*pi*xmw)
         ahvv(2,1,2) = -A_sum 
         ahvv(2,2,1) = -A_sum
c     
         if(xmh.gt.xmz) then
            call width_hgamff(ghgamz)
         else
            ghgamz = 0.0d0
         endif

         if (replace(xgh,3) .and. (.not. replace(bhgamz,3))) then
            xgh = bhgamz*ghgamz
         end if
      else 
         ghgamz = bhgamz*xgh
      end if  ! end of Higgs decay to gamma Z

c H total decay width and branching ratios
 50   if (replace(xgh,3)) then
         XGH = GHWW + GHZZ + GHBB + GHTT + GHCC + GHTAU + GHMU 
     $        + GHGG + GHGAM + GHGAMZ
      end if
c
      BHWW  = GHWW/XGH
      BHZZ  = GHZZ/XGH
      BHGG  = GHGG/XGH
      BHTT  = GHTT/XGH
      BHBB  = GHBB/XGH
      BHTAU = GHTAU/XGH
      BHMU  = GHMU/XGH
      BHCC  = GHCC/XGH
      BHGAM = GHGAM/XGH
      BHGAMZ = GHGAMZ/XGH      

** Checking that the relevant branching ratio is less than one
      call checkBR


** top width: this is calculated if we're not using a slha file or if it's not 
** present in the slha file
      if (replace(gamt,0)) then      
         gamt = gf*xmt**3/(8d0*pi*sqrt(2d0))*(1-(xmw/xmt)**2)**2 *
     -        (1+2d0*(xmw/xmt)**2)*
     -        (1 - 2d0*pwhg_alphas(xmt**2,st_lambda5MSB,st_nlight)/(3d0*pi)*(2*pi**2/3d0 - 2.5d0))
      end if


c width info for higgs
c
      XMG(6) = XMH*XGH 
c
c dummy xmg for commbon block BKOPOUshort for higgs
c
      xmgs(6) = xmg(6)      
c
      print *," "
      print *,"  Information on couplings, masses and widths"
      print *,"==============================================="
      print *," "
      print *,"  SM coupling parameters "
      print*,"---------------------------"
      write(*,"(T4,A,T36,A,F14.12)") "G_Fermi",": ",GF
      write(*,"(T4,A,T36,A,F10.8)") "SIN**2 of the weak mixing angle",
     &     ": ",SIN2W
      write(*,"(T4,A,T36,A,G13.8)") "Q.E.D. fine structure constant",
     &     ": 1/",1.D0/ALFA

!       write(*,"(T4,A,T36,A,F10.8)") "Q.C.D. LO alfas(MZ)",
!      &     ": ",ALFAS_LO
      write(*,"(T4,A,T36,A,F10.8)") "Q.C.D. NLO alfas(MZ)",
     &     ": ",ALFAS 
c output, which is not so relevant to the user
      if (.false.) then
         write(*,*) " "
         write(*,12) C,S,Z,W,Q
   12 FORMAT(/" C,S,Z,W,Q =",5F10.6)
      write(6,*) " "
      write(6,13) ((V(J,K),K=1,5),J=1,4)
      write(6,14) ((A(J,K),K=1,5),J=1,4)
      write(6,15) ((CLR(J,K,-1),K=1,5),J=1,4)
      write(6,16) ((CLR(J,K,1),K=1,5),J=1,4)
      write(6,17) (((B(I,J,K),K=1,6),J=1,6),I=1,6)
   13 FORMAT(" FFB VECTOR COUPLING ( ROW=F, COL= B )"/,4(5F10.6/))
   14 FORMAT(" FFB AXIAL-VECTOR COUPLING ( ROW=F, COL= B )"/,4(5F10.6/))
   15 FORMAT(" FFB LEFTHANDED COUPLING ( ROW=F, COL= B )"/,4(5F10.6/))
   16 FORMAT(" FFB RIGHTHANDED COUPLING ( ROW=F, COL= B )"/,4(5F10.6/))
   17 FORMAT(" BBB NON-ABELIAN COUPLINGS"/,6( 6(6F10.6/)/) )
      write(6,*) " "
      write(6,19) (XM2(J),XMG(J),J=1,6)
   19 FORMAT(/," MASS**2 AND MASS*WIDTH VALUES"/,6(2F10.2/))
      endif

      print *," "
      print *,"     particle masses"
      print*,"---------------------------"
      write(*,"(T4,A,T36,A,F9.4)") "mass of the Z0 boson",": ",XMZ
      write(*,"(T4,A,T36,A,F9.4)") "mass of the W+/- boson",": ",XMW
      write(*,"(T4,A,T36,A,F9.4)") "mass of the Higgs boson",": ",XMH
      write(*,"(T4,A,T36,A,F9.4)") "mass of the top quark",": ",XMT
      write(*,"(T4,A,T36,A,F9.4)") "mass of the bottom quark at MH",
     &     ": ",BMASS_H
      write(*,"(T3,A)") "(All other fermions have mass equal ZERO)"
      IF (INFO.LE.1) RETURN
      print *," "
      print *,"     particle widths"
      print*,"---------------------------"
      write(*,"(T4,A,T36,A,F10.6,A)") "total widths of Z0",": ",
     &     XGZ," GeV"
      write(*,"(T4,A,T36,A,F10.6,A)") "total widths of W+ or W-",": ",
     &     XGW," GeV"
      write(*,"(T4,A,T36,A,F10.6,A)") "total widths of Higgs",": ",
     &     XGH," GeV"
      write(*,"(T4,A,T36,A,F10.6,A)") "total widths of top quark",": ",
     &     gamt," GeV"

      IF (.true.) then 
         print*," "         
         PRINT*," partial Higgs decay width "
         print*,"---------------------------"
         write(*,"(T4,A,T36,A,D12.6,A)") "H -> WW",": ",ghww," GeV"
         write(*,"(T4,A,T36,A,D12.6,A)") "H -> ZZ",": ",ghzz," GeV"
         write(*,"(T4,A,T36,A,D12.6,A)") "H -> bottom anti-bottom",
     &        ": ", ghbb," GeV"
         write(*,"(T4,A,T36,A,D12.6,A)") "H -> top anti-top",
     &        ": ",ghtt," GeV"
         write(*,"(T4,A,T36,A,D12.6,A)") "H -> charm anti-charm",
     &        ": ",  ghcc," GeV"
         write(*,"(T4,A,T36,A,D12.6,A)") "H -> tau anti-tau",": ",
     &        ghtau," GeV"
         write(*,"(T4,A,T36,A,D12.6,A)") "H -> mu anti-mu",": ",
     &        ghmu," GeV"
         write(*,"(T4,A,T36,A,D12.6,A)") "H -> gluon gluon",": ",
     &        ghgg," GeV"
         write(*,"(T4,A,T36,A,D12.6,A)") "H -> gamma gamma",": ",
     &        ghgam," GeV"
         write(*,"(T4,A,T36,A,D12.6,A)") "H -> Z gamma",": ",
     &        ghgamz," GeV"
      ENDIF

      print *," "
      print *,"    branching ratios"
      print*,"---------------------------"
      write(*,"(T4,A,T36,A,F10.7,A)") "W -> lepton anti-lepton",": ",
     &     BWNE*100.d0," %"
      write(*,"(T4,A,T36,A,F10.7,A)") 
     &     "W -> massless quark anti-quark",": ", BWUD*100.d0," %"
      write(*,"(T4,A,T36,A,F10.7,A)") "W -> top bottom",": ",
     &     BWTB*100.d0," %"
      print *," "
      write(*,"(T4,A,T36,A,F10.7,A)") 
     &     "Z -> neutrino anti-neutrino",": ",BZNN*100.d0," %"
      write(*,"(T4,A,T36,A,F10.7,A)") 
     &     "Z -> charged lepton anti-lepton",": ",BZEE*100.d0," %"
      write(*,"(T4,A,T36,A,F10.7,A)") 
     &     "Z -> up-type quark anti-quark",": ",BZUU*100.d0," %"
      write(*,"(T4,A,T36,A,F10.7,A)") 
     &     "Z -> down-type quark anti-quark",": ",BZDD*100.d0," %"
      write(*,"(T4,A,T36,A,F10.7,A)") "Z -> top anti-top",": ",
     &     BZTT*100.d0," %"
      print *," "
      write(*,"(T4,A,T36,A,F10.7,A)") "H -> WW",": ",bhww*100.d0," %"
      write(*,"(T4,A,T36,A,F10.7,A)") "H -> ZZ",": ",bhzz*100.d0," %"
      write(*,"(T4,A,T36,A,F10.7,A)") "H -> bottom anti-bottom",": ",
     &     BHBB*100.d0," %"
      write(*,"(T4,A,T36,A,F10.7,A)") "H -> top anti-top",": ",
     &     BHTT*100.d0," %"
      write(*,"(T4,A,T36,A,F10.7,A)") "H -> charm anti-charm",": ",
     &     BHcc*100.d0," %"
      write(*,"(T4,A,T36,A,F10.7,A)") "H -> tau anti-tau",": ",
     &     BHTAU*100.d0," %"
      write(*,"(T4,A,T36,A,F10.7,A)") "H -> mu anti-mu",": ",
     &     BHMU*100.d0," %"
      write(*,"(T4,A,T36,A,F10.7,A)") "H -> gluon gluon",": ",
     &     BHGG*100.d0," %"
      write(*,"(T4,A,T36,A,F10.7,A)") "H -> gamma gamma",": ",
     &     BHGAM*100.d0," %"
      write(*,"(T4,A,T36,A,F10.7,A)") "H -> Z gamma",": ",
     &     BHGAMZ*100.d0," %"
      print *," "

      if (ldebug) stop

      RETURN

      END
      

********************************************************************************
********************************************************************************
c                         
c     The subroutine below calculate the partial width for H-> V*V* semi-
c     analytically.
      
      function fun(xx,yy)
      implicit none

      double precision lamb1,lamb2,xx,yy,fun,
     $     fun5de,fun5do,funsm5de,funsm
c     
      include 'pwhg_math.h'

      double precision g5hvv(2,1:4,1:4),lambda5
      common/hcoupl/g5hvv,lambda5
      double precision XM2s(6),XMGs(6)
      COMMON /BKOPOUshort/ XM2s,XMGs
      double precision xmv,xmw, xmh
      double precision lambda,constants
      double precision Q2(2),prop1,prop2,DV(2:4,2)
      integer id,bos,bos1,bos2,i,j
      common /partid/ id
      double complex ahvv(3,4,4), ahvvL(3,4,4)
      common/tensorhvv/ ahvv, ahvvL

      bos = id
      if(id.eq.4) then
         bos1 = 3
         bos2 = 4
      else
         bos1 = 2
         bos2 = 2
      endif
      xmv = sqrt(xm2s(bos1))
      xmh = sqrt(xm2s(6))
      xmw = sqrt(xm2s(4))
c      print*,"mh=",xmh
      call calQ2(xx,Q2(1))
c      Q2(1) = xx 
      lamb1 = Q2(1)/xmh**2
c      
      call calQ2(yy,Q2(2))
c      Q2(2) = yy
      lamb2 = Q2(2)/xmh**2
c
      call lambda1(1.0d0,lamb1,lamb2,lambda) !call 1 time for all
      funsm = lambda * (1.0d0 + lamb1**2 + lamb2**2 + 
     1       10.0d0*lamb1*lamb2 -  2.0d0 * lamb1 -  2.0d0 * lamb2)
     
      constants = dble(ahvv(1,bos1,bos2))**2 + 
     1     dimag(ahvv(1,bos1,bos2))**2
      constants = constants * xmh**3/((PI**3) * 64.0d0 * xmv**4)

      funsm = funsm * constants

  
      funsm5de = 0d0
      fun5do = 0d0
      fun5de = 0d0


c     cp odd contribution
      fun5do = lambda * (1.0d0 + lamb1**2 + lamb2**2 - 
     1     2.0d0*lamb1*lamb2 -  2.0d0 * lamb1 -  2.0d0 * lamb2)
     2     *lamb1*lamb2
      constants = g5hvv(2,bos1,bos2)**2 * 
     $     xmh**7/( 8.0d0 * xmv**4 * (PI**3) * (lambda5**2))
      fun5do = fun5do * constants
      
c     cp even contribution 
      fun5de = lambda * (1.0d0 + lamb1**2 + lamb2**2 + 
     1     4.0d0*lamb1*lamb2 -  2.0d0 * lamb1 -  2.0d0 * lamb2)
     2     *lamb1*lamb2
      constants = g5hvv(1,bos1,bos2)**2 * 
     $     xmh**7/( 8.0d0 * xmv**4 * (PI**3) * (lambda5**2))
      fun5de = fun5de*constants


c     interefence of SM and CP even 
      funsm5de = lambda * (1.0d0 - lamb1 - lamb2) 
     2     *lamb1*lamb2
c
c id = 1 for photon, 2 for Z , 3 for W , 4 for W
c
      constants = -g5hvv(1,bos1,bos2)*dble(ahvv(1,bos1,bos2))* 
     $     3.0d0*xmh**5/( 8.0d0 * xmv**4 * (PI**3) * lambda5)
      funsm5de = funsm5de * constants

      fun = funsm + funsm5de + fun5do + fun5de !sum 


      return 

      end


********************************************************************************
********************************************************************************
c
      subroutine calQ2(x,Q2)
      implicit none
      double precision XM2s(6),XMGs(6)
      COMMON /BKOPOUshort/ XM2s,XMGs
      double precision x,Q2
      integer id                !id = 3,4 W"s 2 is Z
      common /partid/ id
c      for W bosons
      Q2 = xm2s(id) + xmgs(id) * dtan(x)
      return

      end


********************************************************************************
********************************************************************************
c
      subroutine calZ(QQ2,Z)
      implicit none
      double precision XM2s(6),XMGs(6)
      COMMON /BKOPOUshort/ XM2s,XMGs
      double precision QQ2,Z
      integer id
      common /partid/ id

c      print*,"mv=",sqrt(xm2s(id))
      Z = datan((QQ2 - xm2s(id))/xmgs(id))
c      
      return

      end


********************************************************************************
********************************************************************************

      subroutine lambda1(x,y,z,lambda)
      implicit none
      double precision x,y,z,lambda2, lambda
c
      lambda2 = x**2 + y**2 + z**2 - 2.0d0*x*y - 2.0d0*y*z - 2.0d0*x*z
      if (lambda2.gt.0d0) then
         lambda = sqrt(lambda2)
      else
         lambda = 0d0
      endif
      return

      end


********************************************************************************
********************************************************************************

      subroutine quad2d(x1,x2,ss)
      implicit none
      double precision eps
      parameter(eps = -1.0d-5)
      double precision ss,x1,x2,h,gaus2
      external h
      external gaus2
c
      ss = gaus2(h,x1,x2,eps)   !partial width for h-> V*V* -> all
      return

      end


********************************************************************************
********************************************************************************
c
      function ff2(yy)
      implicit none
      double precision ff2,yy,fun,x,y
      common/xy/ x,y

      y = yy
      
      ff2 = fun(x,y)
      return
 
      end


********************************************************************************
********************************************************************************
c
      function h(xx)
      implicit none
      double precision eps
      parameter(eps = 1.0d-7)
      double precision h,xx,y1,y2,x,y,ff2
      common /xy/ x,y
      double precision ss,gaus
      external ff2
      external gaus
      double precision Q2,q2max
      double precision  xmh
      double precision XM2s(6),XMGs(6)
      COMMON /BKOPOUshort/ XM2s,XMGs
      integer id,boson
      common /partid/ id
c     
      xmh = sqrt(xm2s(6))
      x = xx
      call calQ2(x,Q2)
c      Q2 = x
      q2max = (xmh - sqrt(Q2))**2
      call calZ(q2max,y2) 
c      y2 = q2max                !no BW map
      call calZ(0.0d0,y1)
c      y1 = 0.0d0

      ss = gaus(ff2,y1,y2,eps)
      if(id.eq.2) then          !symmetry factor of h->ZZ
         h =ss/2.0d0
      else
         h = ss
      endif

      return

      end


********************************************************************************
********************************************************************************

      subroutine initbosid(boson)
      implicit none
      integer id,boson
      common /partid/ id

      id = boson

      end


********************************************************************************
********************************************************************************
c
* THIS IS ITERATIVE INTEGRATION PROCEDURE                               
* ORIGINATES  PROBABLY FROM CERN LIBRARY                                
* IT SUBDIVIDES INTEGRATION RANGE UNTIL REQUIRED PRECISION IS REACHED    
* PRECISION IS A DIFFERENCE FROM 8 AND 16 POINT GAUSS ITEGR. RESULT     
* EEPS POSITIVE TREATED AS ABSOLUTE PRECISION                           
* EEPS NEGATIVE TREATED AS RELATIVE PRECISION                           
      DOUBLE PRECISION FUNCTION GAUS(F,A,B,EEPS)                        
*     *************************             
      implicit none

      double precision A, B, AA, BB, C1, C2, CONST, DELTA, eeps, EPS
      double precision F, S8, S16, U, W, X, Y
      integer i
      integer count
      EXTERNAL F                                                        
      DIMENSION W(12),X(12)                                             
      DATA CONST /1.0D-19/                                              
      DATA W                                                            
     1/0.10122 85362 90376, 0.22238 10344 53374, 0.31370 66458 77887,   
     2 0.36268 37833 78362, 0.02715 24594 11754, 0.06225 35239 38648,   
     3 0.09515 85116 82493, 0.12462 89712 55534, 0.14959 59888 16577,   
     4 0.16915 65193 95003, 0.18260 34150 44924, 0.18945 06104 55069/   
      DATA X                                                            
     1/0.96028 98564 97536, 0.79666 64774 13627, 0.52553 24099 16329,   
     2 0.18343 46424 95650, 0.98940 09349 91650, 0.94457 50230 73233,   
     3 0.86563 12023 87832, 0.75540 44083 55003, 0.61787 62444 02644,   
     4 0.45801 67776 57227, 0.28160 35507 79259, 0.09501 25098 37637/   

      count = 0

      EPS=dABS(EEPS)                                                    
      DELTA=CONST*dABS(A-B)                                             
      GAUS=0.d0                                                         
      AA=A                                                              
    5 Y=B-AA             
      IF(dABS(Y) .LE. DELTA) RETURN                                     
    2 BB=AA+Y              
      count = count + 1
      if (count .ge. 150000) then
         write(*,*)'Sorry!  We are having problems calculating a '
         write(*,*)'decay width.  Try altering the anomalous couplings'
         write(*,*)'or putting the width in by-hand using a SLHA file'
         stop
      end if      
      C1=0.5d0*(AA+BB)                                                  
      C2=C1-AA                                                          
      S8=0.d0                                                           
      S16=0.d0                                                          
      DO I=1,4                                                        
         U=X(I)*C2                                                         
         S8=S8+W(I)*(F(C1+U)+F(C1-U))   
      end do
      DO I=5,12                                                       
         U=X(I)*C2                                                         
         S16=S16+W(I)*(F(C1+U)+F(C1-U)) 
      end do
      S8=S8*C2                                                          
      S16=S16*C2           
      IF(EEPS.LT.0D0) THEN       
        IF(dABS(S16-S8) .GT. EPS*dABS(S16)) GO TO 4                     
      ELSE                
        IF(dABS(S16-S8) .GT. EPS) goto 4
      ENDIF                                                             
      GAUS=GAUS+S16 
      AA=BB                                       
      GO TO 5                                                          
    4 Y=0.5d0*Y                        
      IF(dABS(Y) .GT. DELTA) GOTO 2                                     
      PRINT 7                                                           
      GAUS=0.d0                   
      RETURN                                                            
    7 FORMAT(1X,36HGAUS  ... TOO HIGH ACCURACY REQUIRED)                

      END                                                               


********************************************************************************
********************************************************************************
                                                                        
      DOUBLE PRECISION FUNCTION GAUS2(F,A,B,EEPS)                       
*     *************************       
      implicit none

      integer count

      double precision A, B, AA, BB, F, EPS, EEPS, CONST
      double precision C1, C2, DELTA, S8, S16, U, W, X, Y
      integer i
      EXTERNAL F                                                        
      DIMENSION W(12),X(12)                                             
      DATA CONST /1.0D-19/                                              
      DATA W                                                            
     1/0.10122 85362 90376, 0.22238 10344 53374, 0.31370 66458 77887,   
     2 0.36268 37833 78362, 0.02715 24594 11754, 0.06225 35239 38648,   
     3 0.09515 85116 82493, 0.12462 89712 55534, 0.14959 59888 16577,   
     4 0.16915 65193 95003, 0.18260 34150 44924, 0.18945 06104 55069/   
      DATA X                                                            
     1/0.96028 98564 97536, 0.79666 64774 13627, 0.52553 24099 16329,   
     2 0.18343 46424 95650, 0.98940 09349 91650, 0.94457 50230 73233,   
     3 0.86563 12023 87832, 0.75540 44083 55003, 0.61787 62444 02644,   
     4 0.45801 67776 57227, 0.28160 35507 79259, 0.09501 25098 37637/   

      count = 0

      EPS=dABS(EEPS)                                                    
      DELTA=CONST*dABS(A-B)                                             
      GAUS2=0.D0                                                        
      AA=A                                                              
    5 Y=B-AA                                                          
      IF(dABS(Y) .LE. DELTA) RETURN                                     
    2 BB=AA+Y                          
      count = count + 1
      if (count .ge. 150000) then
         write(*,*)'Sorry!  We are having problems calculating a '
         write(*,*)'decay width.  Try altering the anomalous couplings'
         write(*,*)'or putting the width in by-hand using a SLHA file'
         stop
      end if      
      C1=0.5d0*(AA+BB)                                                  
      C2=C1-AA                                                          
      S8=0.d0                                                           
      S16=0.d0                                                        
      DO 1 I=1,4             
      U=X(I)*C2         
    1 S8=S8+W(I)*(F(C1+U)+F(C1-U))                                      

      DO 3 I=5,12                                                       
      U=X(I)*C2                                                         
    3 S16=S16+W(I)*(F(C1+U)+F(C1-U))                                    
      S8=S8*C2                                                          
      S16=S16*C2                                                     
      IF(EEPS.LT.0D0) THEN                                              
        IF(dABS(S16-S8) .GT. EPS*dABS(S16)) GO TO 4 
      ELSE                                                              
        IF(dABS(S16-S8) .GT. EPS) GO TO 4 
      ENDIF                                                             
      GAUS2=GAUS2+S16                                                   
      AA=BB                                                             
      GO TO 5                                                           
    4 Y=0.5d0*Y                                                         
      IF(dABS(Y) .GT. DELTA) GOTO 2                                     
      PRINT 7                                                           
      GAUS2=0.D0                                                        
      RETURN                                                            
    7 FORMAT(1X,36HGAUS2 ... TOO HIGH ACCURACY REQUIRED)                

      END                                                               


********************************************************************************
********************************************************************************
C This subroutine computes the parametric integrals in the 
c Higgs Hunters Guide on pg 29 (eq. 2.24)

      subroutine compints(a,b,I1,I2)
      implicit none
      double precision a,b
      double complex I1,I2
      double complex fci,g
      external fci
      external g

      I1 = a*b/(2.0d0*(a-b)) + a**2 * b**2/(2.0*(a-b)**2) * (fci(a)-fci(b))
     $     + a**2 * b/((a-b)**2) *(g(a) - g(b))
      I2 = -a*b/(2.0d0*(a -b))* (fci(a) - fci(b))
      return

      end


********************************************************************************
********************************************************************************
c
** Function g from eq. 2.24 of the Higgs Hunters' Guide (p29), used in the
** calculation of H -> Z photon

      double complex function g(aa)
      implicit none
      double precision aa,nup,num
      double precision pi
      parameter(pi =3.141592653589793d0) 

      if(aa.ge.1.0d0) then
         g = sqrt(aa -1.0d0)*ASIN(1.0d0/sqrt(aa))
      else
         nup = 1.0d0 + sqrt(1.0d0 - aa)
         num = 1.0d0 - sqrt(1.0d0 - aa)
         g = 0.5d0* sqrt(1.0d0 - aa)*(log(nup/num)*(1.d0,0.d0)
     $        + pi*(0.d0,-1.d0))
      endif
      return

      end


********************************************************************************
********************************************************************************
c
** Function from eq. 2.19 of Higgs Hunters' Guide (p26), used in calculation of
** H -> Z photon

      double complex function fci(aa)
      implicit none
      double precision aa,nup,num
      double precision pi
      parameter(pi =3.141592653589793d0) 
      
      if(aa.ge.1.0d0) then
         fci = (ASIN(sqrt(1.0d0/aa)))**2
      
      else
         nup = 1.0d0 + sqrt(1.0d0 - aa)
         num = 1.0d0 - sqrt(1.0d0 - aa)
         fci = -0.25d0*(log(nup/num)
     $        + pi*(0.d0,-1.d0))**2
      endif
      return

      end


********************************************************************************
********************************************************************************
c	all couplings, widths, and masses are computed from 
c	values stored in common/ bkopout/
c
	subroutine ctrans(xmb)
c     some coupling parameters are modified here. 
c***********************************************************************
c 	This subroutine sets up the couplings, masses and widths 
c	needed by HELAS as they are given in koppln.f
c***********************************************************************

	implicit none      
 
c common blocks which need to be filled for HELAS (masses, widths, couplings):	
	
        include "coupl.inc"
        include "scales.inc"
			
      double precision g2,ee2,sw,cw,vv,sc2,pi		
      parameter ( pi = 3.141592653589793d0 )
      logical loutput
      parameter (loutput = .false.)
	
c     for linit=0, alfas is filled according to the input value
c     for linit=1, alfas is filled according to scales.inc input, als(1,1).
      integer linit
      data linit /0/
      save linit
c common block with input values (only needed for top mass xmt):
      double precision alphas, xmt, alfa, xmz, xmw,s2w, xmh, gf
      common /bkopin/ alphas, xmt, alfa, xmz, xmw, s2w, xmh, gf
      double precision xmb
	
c common block which has been filled in koppln.f and contains 
c FFB, BBB, and Higgs couplings (but not in form needed by HELAS):
      double complex clrCT
      common /bkopouCT/ clrCT(3:4,2:4,-1:1)
      DOUBLE PRECISION CLR,XM2,XMG,B,V,A
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),V(4,5),A(4,5)
      

      if (linit.eq.0) then
c fill basic parameters:
	G      =  B(5,5,5)
	GG(1)  = -G
	GG(2)  = -G
	alfas  =  G**2/(pi*4.d0)  ! strong coupling at fixed ren. scale MZ
        linit = linit+1
      else if (linit.eq.1) then
	G      =  sqrt(als(1,1)*4d0*pi)
	GG(1)  = -G
	GG(2)  = -G
	alfas  =  G**2/(pi*4.d0)  ! strong coupling at selected ren. scale

        return
      endif

      gfermi = gf

	ee     = -CLR(2,1,-1)
	ee2    =  ee**2
	alpha  =  ee2/(pi*4.d0)
	
	g2     = sqrt(2.d0)*CLR(2,3,-1)
	cw     = B(2,4,3)/g2
	sin2w  = 1.d0 - cw**2
	sw     = sqrt(sin2w)
	
	
c fill masses:

	hmass = sqrt(xm2(6))
	wmass = sqrt(xm2(3))
	zmass = sqrt(xm2(2))
	amass = 0.d0
	tmass = xmt
        bmass = xmb
	cmass = 0.d0
	lmass = 0.d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
	
c  fill widths:			

	hwidth = xmg(6)/hmass 
	wwidth = xmg(3)/wmass
	zwidth = xmg(2)/zmass
	awidth = 0.d0	
	lwidth =  0d0
	twidth = 1.6d0	! no top emerges in WW-prod. -> value here "arbitrary"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
        	
c  auxiliary variables:	
	
	vv  = 2.d0*wmass*sw/ee
	sc2 = sin2w*(1.d0-sin2w)
	
c  fill gauge couplings: 
c
c	fermion-fermion-vector couplings: 

	gal(1) = -dcmplx(clr(2,1,-1))	! photon
	gal(2) = -dcmplx(clr(2,1, 1))
	gau(1) = -dcmplx(clr(3,1,-1))	
	gau(2) = -dcmplx(clr(3,1, 1))	
	gad(1) = -dcmplx(clr(4,1,-1))	
	gad(2) = -dcmplx(clr(4,1, 1))	
	
	gwf(1) = -dcmplx(clr(1,3,-1))	! W+- boson (same for all fermions)
	gwf(2) = -dcmplx(clr(1,3, 1))	! zero
	
	gzn(1) = -dcmplx(clr(1,2,-1))	! Z boson
	gzn(2) = -dcmplx(clr(1,2, 1))	! zero	
	gzl(1) = -dcmplx(clr(2,2,-1))
	gzl(2) = -dcmplx(clr(2,2, 1))	
	gzu(1) = -dcmplx(clr(3,2,-1))
	gzu(2) = -dcmplx(clr(3,2, 1))	
	gzd(1) = -dcmplx(clr(4,2,-1))
	gzd(2) = -dcmplx(clr(4,2, 1))
	
c  	vector boson couplings:

	gw   = ee/sw
	gwwa = -dcmplx(b(3,4,1))
	gwwz = -dcmplx(b(3,4,2))
	
c	gauge - higgs boson couplings:

	gwwh  = dcmplx(b(3,4,6)*wmass)	
	gzzh  = dcmplx(b(2,2,6)*wmass)	
	ghhh  = dcmplx(-hmass**2/vv*3d0)
	gwwhh = dcmplx( ee2/sin2w*0.5d0 )
        gzzhh = dcmplx( ee2/sc2*0.5d0)
        ghhhh = ghhh/vv
	

c-----------------------------------------------

c for reference:
c	info on couplings returned by ctrans:

      if (loutput) then
	write(6,*) "couplings returned by ctrans:"
	write(6,*) "alpha,ee, sin2w, alfas,g,gg:"
	write(6,*)  alpha,ee, sin2w, alfas,g,gg
	write(6,*) 
	write(6,*) "hmass, wmass, zmass, amass:"
	write(6,*)  hmass, wmass, zmass, amass
	write(6,*) 
	write(6,*) "tmass, bmass, lmass, cmass:"
	write(6,*)  tmass, bmass, lmass, cmass
	write(6,*) 
	write(6,*) "hwidth, wwidth, zwidth:" 
	write(6,*)  hwidth, wwidth, zwidth
	write(6,*) 
	write(6,*) "twidth, lwidth, awidth:"
	write(6,*)  twidth, lwidth, awidth
	write(6,*) 
	write(6,*) "gal   , gad   , gau   , gwf:"
	write(6,*)  gal   , gad   , gau   , gwf
	write(6,*) 
	write(6,*) "gzn   , gzl   , gzd   , gzu:"
	write(6,*)  gzn   , gzl   , gzd   , gzu
	write(6,*) 
	write(6,*) "gw, gwwa, gwwz:"
	write(6,*)  gw, gwwa, gwwz
	write(6,*) 
	write(6,*) "gwwh, gzzh, gwwhh, gzzhh, ghhh, ghhhh:"
	write(6,*)  gwwh, gzzh, gwwhh, gzzhh, ghhh, ghhhh
	write(6,*) 
      endif
	
c--------------------------------------------------------------------------	

	return
	end


********************************************************************************
********************************************************************************
c
c PHASE SPACE REDUCTION FACTOR FOR HEAVY FERMIONIC DECAYS
c
      double precision function heavy(X,Y1,Y2)

      implicit none

      double precision X, Y1, Y2
      
      HEAVY = ( 1.d0 - .5d0*(Y1**2 + Y2**2) - 
     &           0.5d0*(Y1**2 - Y2**2)**2 +
     &           3.d0*Y1*Y2*((X**2 - 1.d0)/(X**2 + 1.d0)) )*
     &           SQRT((1.d0 - Y1**2 - Y2**2)**2 - 4.d0*Y1**2*Y2**2)

      end 


********************************************************************************

