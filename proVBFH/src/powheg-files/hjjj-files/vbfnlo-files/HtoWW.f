********************************************************************************
********************************************************************************
*** HtoWW.F                                                                  ***
*** 27 July 2011                                                             ***
*** sophy@particle.uni-karlsruhe.de, rauch@particle.uni-karlsruhe.de         ***
***                                                                          ***
*** This file contains the code needed to compute the H -> VV - 4 lepton     ***
*** matrix squared elements.  These were originally duplicated in the        ***
*** folders amplitudes/vvjj and amplitudes/hjjj, and have now been moved     ***
*** here.  Used for Hjj via VBF and gluon fusion, Hjjj and HAjj.             ***
***                                                                          ***
********************************************************************************
********************************************************************************

c     F1,F2,F3
      subroutine computeFs(qsq1,qsq2,tau1,tau3,flav1,flav3,
     $     bos1,bos2,F1,F2,F3)
      implicit none
      real*8 qsq1,qsq2          !q1^2 and q2^2
      integer tau1,tau3         !helicities of f1 and f3
      integer flav1,flav3       !flavor of fermions
      integer bos1,bos2         !boson type : W+ 3
c                                             W- 4
c                                             Z  2 
c                                         photon 1
c
      double complex prop1(4),prop2(4) !1 = photon propagator
c     2 = Z propagator
      double complex F1,F2,F3
c     koppln commom blocks
      double complex ahvv(3,4,4), ahvvL(3,4,4)
      common/tensorhvv/ ahvv, ahvvL
      DOUBLE PRECISION CLR,XM2,XMG,B,V,A
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),V(4,5),A(4,5)


** propagator factors:
      prop1(1) = 1.d0/dcmplx(qsq1,0.0d0)
      prop2(1) = 1.d0/dcmplx(qsq2,0.0d0)
      prop1(2) = 1.d0/dcmplx(qsq1-xm2(2),xmg(2))
      prop2(2) = 1.d0/dcmplx(qsq2-xm2(2),xmg(2))
      prop1(3) = 1.d0/dcmplx(qsq1-xm2(3),xmg(3))
      prop2(3) = 1.d0/dcmplx(qsq2-xm2(3),xmg(3))
      prop1(4) = prop1(3)
      prop2(4) = prop2(3)

      if(bos1.eq.2) then  ! H -> ZZ -> llll
         F1 = CLR(flav1,2,tau1)*CLR(flav3,2,tau3)*ahvv(1,2,2)*
     1                              prop1(2)*prop2(2)
         
** NOTE: this includes effective H-gamma-Z and H-gamma-gamma vertices in the SM
         F2 = CLR(flav1,2,tau1)*CLR(flav3,2,tau3)*ahvv(2,2,2)*
     1                              prop1(2)*prop2(2)+
     $       V(flav1,1)*CLR(flav3,2,tau3)*ahvv(2,2,1)*prop1(1)*prop2(2)+
     $       V(flav3,1)*CLR(flav1,2,tau1)*ahvv(2,2,1)*prop1(2)*prop2(1)+
     $       V(flav1,1)*V(flav3,1)*ahvv(2,1,1)*prop1(1)*prop2(1)
c
         F3 = CLR(flav1,2,tau1)*CLR(flav3,2,tau3)*ahvv(3,2,2)*
     1                              prop1(2)*prop2(2)+
     $       V(flav1,1)*CLR(flav3,2,tau3)*ahvv(3,2,1)*prop1(1)*prop2(2)+
     $       V(flav3,1)*CLR(flav1,2,tau1)*ahvv(3,2,1)*prop1(2)*prop2(1)+
     $       V(flav1,1)*V(flav3,1)*ahvv(3,1,1)*prop1(1)*prop2(1)

      else  ! H -> WW -> lvlv

         F1 = CLR(flav1,3,tau1)*CLR(flav3,4,tau3)*ahvv(1,3,4)*
     1                              prop1(3)*prop2(4)
c         
         F2 = CLR(flav1,3,tau1)*CLR(flav3,4,tau3)*ahvv(2,3,4)*
     1                              prop1(3)*prop2(4)
c     
         F3 = CLR(flav1,3,tau1)*CLR(flav3,4,tau3)*ahvv(3,3,4)*
     1                              prop1(3)*prop2(4)

      endif

      return

      end
c

