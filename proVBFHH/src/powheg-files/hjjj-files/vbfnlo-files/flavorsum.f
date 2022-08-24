********************************************************************************
********************************************************************************

c     This subroutine computes the decay width of 
c     h -> Z/gamma Z/gamma -> f1 f1bar f3 f3bar for 
c     flav1,flav3
c     Final state fermions are treated as non-identical particles.
c
      subroutine flavor_sum(width)
c     only works for h->Z*Z* decay

      implicit none

      include "koppln_ew.inc"
      include "para_blocks.inc"

      real*8 width,dum
      real*8 matrix(11,11),mass(11),Nc(4),gm(4),matrix1(4,4)
      real*8 Nf(4),factor
c     nu_e,e,nu_mu,mu,nu_tau,tau,u,d,c, s, b
c     1    2   3   4    5    6   7 8 9 10 11
      real*8 flavor(11)
      integer i,j,k,l
      real*8 pi
      parameter(pi=3.141592653589793d0)
      logical lcomputeALL
      data flavor/1,2,1,2,1,2,3,4,3,4,4/
      parameter(lcomputeALL=.false.)

c     color factors 
      Nc(1) = 1.d0
      Nc(2) = 1.d0
      Nc(3) = 3.d0*(1.d0 + ALFAS/PI)   
      Nc(4) = 3.d0*(1.d0 + ALFAS/PI)

c     Number of fermion in generation
      Nf(1) = 3.0d0
      Nf(2) = 3.0d0 
      Nf(3) = 2.0d0             !no top
      Nf(4) = 3.0d0

c     compute geometric means of masses 
      gm(1) = 0d0     ! NEUTRINO MASS
      gm(2) = (Mf(2,1)*Mf(2,2)*Mf(2,3))**(1.0d0/3.d0)
      gm(3) = (Mf(3,1)*Mf(3,2))**(1.0d0/2.d0)
      gm(4) = (Mf(4,1)*Mf(4,2)*Mf(4,3))**(1.0d0/3.d0)

c    setting mass
      mass(1) = Mf(1,1)
      mass(2) = Mf(2,1)
      mass(3) = Mf(1,2)
      mass(4) = Mf(2,2)
      mass(5) = Mf(1,3)
      mass(6) = Mf(2,3)
      mass(7) = Mf(3,1)
      mass(8) = Mf(4,1)
      mass(9) = Mf(3,2)
      mass(10) = Mf(4,2)
      mass(11) = Mf(4,3)

      if(lcomputeALL) then
         dum = 0.0d0
         do k=1,11
            do l = k,11
               if (l.eq.k) then
                  factor = 2D0
               else
                  factor = 1D0
               end if
               i = flavor(k)
               j = flavor(l)
               call hzgammawidth(mass(k),mass(l),i,j,matrix(k,l))
               dum = matrix(k,l)*Nc(i)*Nc(j)/factor + dum
            enddo
         enddo
      else
         dum =0.0d0
         factor = 2.0d0
         do i=1,4
            do j = 1,4         
               call hzgammawidth(gm(i),gm(j),i,j,matrix1(i,j))
               dum = matrix1(i,j)*Nc(i)*Nc(j)*Nf(i)*Nf(j)/factor + dum
            enddo
         enddo
      endif
      width = dum

      return
      end


********************************************************************************
c
*** Subroutine called in the calculation of h -> ZZ -> 4f, from flavor_sum
*** mf1 and mf3 and the fermion masses 
*** iflav1 and iflav3 are the fermion flavours 
*** flavor_sum runs over all of the end-state fermions

      subroutine hzgammawidth(mf1,mf3,iflav1,iflav3,width)

      implicit none

      real*8 width,mf1,mf3,mh,qsq1,qsq2,ss1,ss2
      integer flav1,flav3
      integer iflav1,iflav3      
      common/flavors/flav1,flav3
      real*8 mu1,mu3,x1,x2
      common/masses/mu1,mu3
      DOUBLE PRECISION CLR,XM2,XMG,B,V,A
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),V(4,5),A(4,5)
      integer id                !id = 3,4 W"s 2 is Z
      common /partid/ id 
         
      id = 2                    !always -- i.e. Z
      flav1 = iflav1
      flav3 = iflav3
      
      mu1 = mf1
      mu3 = mf3
      

      x1 = 0.0d0
      x2 = 1.0d0

      call quad2dim(x1,x2,ss2)
      width = ss2 

      return
      end
c    
c
********************************************************************************
*** This is the function that's (eventually) called to calculate the 
*** H -> ZZ -> 4f width
*** It's (I think) the equivalent of the function 'fun' in koppln.F that's 
*** used to calculated H -> ZZ in the SM
      
      real*8 function dgam(xx,yy)

      implicit none

      integer flav1,flav3,i,j
      common/flavors/flav1,flav3
      real*8 qsq1,qsq2,xx,yy
      double complex F1,F2,F3,F1c,F2c,F3c
      real*8 mu1,mu3
      common/masses/mu1,mu3
      real*8 dum1
      real*8 ints(3,3),lamb1,lamb2,factor,pi
      DOUBLE PRECISION CLR,XM2,XMG,B,V,A
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),V(4,5),A(4,5)
      parameter(pi=3.141592653589793d0)
      real*8 Jac2,Jac1
      real*8 lowq2(2),hiq2(2)
 
      
c     convention here is 2 is yy
c     1 is xx 
c      qsq1 = yy
c      qsq2 = xx
      lowq2(1) = 4.0d0*mu1**2
      hiq2(1) = xm2(6)
         
      call com_Jq2(1,lowq2,hiq2,xx,qsq1,Jac1) !outer integral
      
c     compute bound on inner integration dy or 2
      lowq2(2) = 4.0d0*mu3**2
      hiq2(2) = (sqrt(xm2(6))-sqrt(qsq1))**2

      call com_Jq2(2,lowq2,hiq2,yy,qsq2,Jac2) !inner integral

      dum1 = 0.0d0
      call computeIs(qsq1,qsq2,ints)

      do i=-1,1,2
         do j = -1,1,2

** this subroutine is in HtoWW.F.  It calculates the combinations of H-V and
** V-f couplings
            call computeFs(qsq1,qsq2,i,j,flav1,flav3,2,2,
     $           F1,F2,F3)
  

            F1c = conjg(F1)
            F2c = conjg(F2)
            F3c = conjg(F3)

            dum1 = (dble(f1)**2+dimag(f1)**2)*ints(1,1)+
     $           (dble(f2)**2+dimag(f2)**2)*ints(2,2)+
     $           (dble(f3)**2+dimag(f3)**2)*ints(3,3)+
     $           2.0d0*dble(f1*f2c)*ints(2,1) + dum1
         enddo
      enddo

      lamb1 = qsq1/xm2(6)
      lamb2 = qsq2/xm2(6)
      call lambda1(1.0d0,lamb1,lamb2,factor)
      dgam = dum1*factor*(1.0d0/(8.0d0*pi))*
     $     (1.0d0/(32.d0*pi**2))**2 
     $     *(1.d0/(2.0d0*pi))**2 * 1.d0/(2.0d0*sqrt(xm2(6)))

      dgam = dgam*Jac1*Jac2     !mult by jacobians

      return
      end
c      


********************************************************************************
c
*** Function used in calculation of H -> ZZ -> 4f

      real*8 function ddwidth(yy)
      implicit none
      real*8 yy,x,y,dgam
      common/xy/ x,y
      external dgam

      y = yy
      ddwidth = dgam(x,y)

      return
      end


********************************************************************************

*** Function used in calculation of H -> ZZ -> 4f

      real*8 function dwidth(xx)
      implicit none
      real*8 eps
      parameter(eps = 1.0d-5)
      real*8 xx,y1,y2,x,y
      common /xy/ x,y
      real*8 ss,gaus,ddwidth      
      
      DOUBLE PRECISION CLR,XM2,XMG,B,V,A
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),V(4,5),A(4,5)
     
      integer flav1,flav3      
      common/flavors/flav1,flav3
      
     
      external ddwidth
      external gaus     
c
      x = xx
      
      y1 = 0.0d0
      y2 = 1.0d0
      ss = gaus(ddwidth,y1,y2,eps)
        
      dwidth = ss
      return
      end


********************************************************************************

*** Function used in calculation of H -> ZZ -> 4f

      subroutine quad2dim(x1,x2,ss)
      implicit none
      integer i
      real*8 eps
      parameter(eps = -1.0d-2)
      real*8 ss,x1,x2,dwidth,gaus2
      external dwidth
      external gaus2
     
      ss = gaus2(dwidth,x1,x2,eps)  

      return
      end



********************************************************************************

*** Function used in calculation of H -> ZZ -> 4f

      subroutine computeIs(qsq1,qsq2,arrayOFints)

      implicit none

      integer i,j
      real*8 qsq1,qsq2
      real*8 arrayOFints(3,3)
      real*8 q1q2,pi               !q1.q2
      parameter(pi =3.141592653589793d0)
      DOUBLE PRECISION CLR,XM2,XMG,B,V,A
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),V(4,5),A(4,5)

c     initialize
      do i=1,3
         do j=1,3
            arrayOFints(i,j) = 0.0d0
         enddo
      enddo
c     compute q1.q2
      q1q2 = (xm2(6) - qsq1 - qsq2)/2.d0
c
c      Fill arrayOFints
      arrayOFints(1,1) = (64.d0/9.d0)*pi**2 *
     1     (q1q2**2 + 2.0d0*qsq1*qsq2)
      arrayOFints(1,2) = (64.d0/3.d0)*pi**2* qsq1*qsq2* q1q2
      arrayOFints(3,3) = (128.d0/9.d0)*pi**2* 
     $     qsq1*qsq2*(q1q2**2 - qsq1*qsq2)
      arrayOFints(2,1) = arrayOFints(1,2)
      arrayOFints(2,2) = (64.d0/9.d0)*pi**2 *qsq1*qsq2 *
     1     (2.0d0*q1q2**2 + qsq1*qsq2)
      return
      end


********************************************************************************
c
c     compute jacobian and q2.  USed in H -> ZZ -> 4f and H -> Zgamma -> 2f gam

      subroutine com_Jq2(ich,lowq2,hiq2,r,qsq,J)

      implicit none

      integer id                !id = 3,4 W"s 2 is Z
      common /partid/ id
      double precision XM2s(6),XMGs(6)
      COMMON /BKOPOUshort/ XM2s,XMGs
      real*8 r,qsq,J             ! J is the jacobian and r is the variable 
                                 ! being integrated over
      real*8 rmid,rmidl,a,c,b         !c is the midpoint on [a,b]
      real*8 x,xmin,q2cut,xmax
      integer i,ich
 
      real*8 lowq2(2),hiq2(2)
     
      logical lmap,logmapon
      parameter(logmapon=.true.)
c
      q2cut = 60.d0**2
      rmid = 0.2d0
      
      i = ich                   !1 or 2 = inner and outer bounds
c     if hiq2(i) is less than q2cut then let q2cut = hiq2(i)
      if(hiq2(i).lt.q2cut) then ! whole is integrated using ln map
         a = lowq2(i)
         b = hiq2(i)
         q2cut = hiq2(i)
         c = q2cut
         rmidl = 1.0d0
      else 
         a = lowq2(i)
         b = hiq2(i)
         c = q2cut
         rmidl = rmid
      endif
         
      lmap = (r.lt.rmidl).and.(a.gt.0.0d0).and.logmapon
      if(lmap) then        !ln map

         qsq = a * exp(R*(dlog(c/a))/rmidl)
         J = qsq *(dlog(c/a))/rmidl
      else
         if(a.eq.0.0d0) then ! for a =0 always use bw map
            c = a
            rmidl =0.0d0
         else 
            c= q2cut
            rmidl = rmid
         endif
         call calZ(c,xmin)
         call calZ(b,xmax)
         x = xmin +(r - rmidl)/(1.0d0 - rmidl)*(xmax - xmin)
         call calQ2(x,qsq)
         
         J = (xmax - xmin)/(1.0d0 - rmidl) 
     $        * xmgs(id)*(dtan(x)**2 + 1.0d0)
      endif
      return
      end
c


********************************************************************************

*** This is the subroutine that calculates the width 
***          h -> gamma Z -> gamma 2f 
*** It's called from width_hgamff (which sums over fermion type)
*** idf is fermion type, mf is fermion mass

      subroutine widthhgz(idf,mf,q2,width)
      implicit none
      integer idf,i,j               !1 for neutrino,2 for electron,3 for up, 4
                                    !for down
      real*8 cc(3),width(4),ratio(3) !1 is gg* 2 is gz* and 3 is the 
                                     ! interference term
      common/constants/cc
      real*8 q2,FF2(3),psfactor ! q squared of decaying boson
      real*8 delta,mg,mz,mh,mf,Pi,Lamb1,c1,c2,c3
      parameter(Pi= 3.141592653589793d0)
     
      double complex ahvv(3,4,4), ahvvL(3,4,4)
      common/tensorhvv/ ahvv, ahvvL

      DOUBLE PRECISION CLR,XM2,XMG,B,V,A
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),V(4,5),A(4,5)
      real*8 mf1
*      parameter(mf1=0d0)
      mf1 = mf

      mh = Sqrt(xm2(6))
c      write(6,*) mh
      mz = Sqrt(xm2(2))
      mg = xmg(2)
c      mg =0.001d0
c     Put coupling factors here:
      do i =1,3
         FF2(i) = 0.0d0
         cc(i) = 0.0d0
      enddo

* H -> gamma gamma -> ff gamma
      cc(1) = 2.0d0*v(idf,2)**2 * 
     1     ( dble(ahvv(2,1,1))**2+dimag(ahvv(2,1,1))**2 +
     2       dble(ahvv(3,1,1))**2+dimag(ahvv(3,1,1))**2 ) 
* H -> Z gamma -> ff gamma
      cc(2) = (CLR(idf,2,-1)**2 + CLR(idf,2,+1)**2) * 
     1     ( dble(ahvv(2,2,1))**2+dimag(ahvv(2,2,1))**2 +
     2       dble(ahvv(3,2,1))**2 + dimag(ahvv(3,2,1))**2 ) 
* interference: H -> Z gamma -> ff gamma WITH H -> gamma gamma -> ff gamma
      cc(3) = 2.d0*(CLR(idf,2,-1) + CLR(idf,2,+1))*v(idf,1) * 
     1     ( dble(ahvv(2,1,1)*conjg(ahvv(2,1,2)))+
     2       dble(ahvv(3,1,1)*conjg(ahvv(3,1,2))) )
c
      FF2(1) = cc(1)/q2**2 
      FF2(2) = cc(2)/((q2 - mz**2)**2 + mg**2)
      FF2(3) = cc(3)*(q2-mz**2)/((q2 - mz**2)**2 + mg**2)/q2
c
    
       psfactor =(1.0d0/(mh**3))*(1.0d0/(32.0d0*Pi**2))**2 * 
     1           (mh**2 - q2)*sqrt((q2 - 4.0d0*mf1**2)/q2) 

c     We use a fermion mass to regulate photon pole.
      do i = 1,3

         width(i) = FF2(i) * psfactor * (4.0d0/3.0d0)* Pi *
     1             (q2-mf1**2)*(mh**2 - q2)**2 ! m2s is computed with 
                                               ! massive fermions
c         write(6,*)"q2 =",q2
c         write(6,*)"dw/dq2 =",width(i),"for",i
      enddo
      
      width(4) = width(1)+width(2)+width(3)

      end


********************************************************************************

*** Function used in calculation of H -> gamma Z -> gamma ff

      function dw(xx)  ! interface with subroutine dwidth

      implicit none

      integer idf,iflav
      common/flav/ idf
      integer id                !id = 3,4 W"s 2 is Z
      common /partid/ id
      real*8 dw,xx,q2
      integer jj
      common/choice/jj
      real*8 mass(2),mf
* mass1 is set by the arguments of quad1d, mass(1) = mfermion, mass(2) = mhiggs
      common/mass1/mass
      real*8 ss,gaus,dpw(4)
      real*8 jac,lowq2(2),hiq2(2)
      external gaus
      
      id = 2
c     
      lowq2(1) = 4.0d0*mass(1)**2
      hiq2(1) = mass(2)**2
      call com_Jq2(1,lowq2,hiq2,xx,q2,Jac)
      iflav = idf
      mf = mass(1)
      call widthhgz(iflav,mf,q2,dpw)

c      dw  = dpw(2)+dpw(3)+dpw(1)
      dw = dpw(jj)*jac          !mult by jacobian

      return
      end


********************************************************************************
c
** USed in calculation of H -> Z gamma -> ff gamma

      subroutine quad1d(ii,iflav,mf,hmass,ss)
 
      implicit none

      integer idf,iflav
      common/flav/ idf
      integer ii,jj
      common/choice/jj
      real*8 mass(2),mf,hmass
      common/mass1/mass
      real*8 eps
      parameter(eps = 1.0d-7)
      real*8 ss,x1,x2,dw,gaus
      external dw
      external gaus

      jj = ii
      idf = iflav
      mass(1) = mf
      mass(2) = hmass
      x1 = 0.0d0
      x2 = 1.0d0
      ss = gaus(dw,x1,x2,eps)  

      return

      end


********************************************************************************

*** Subroutine called from koppln.F to calculate 
***          h-> gamma/Z* gamma -> f fbar gamma
*** for anomalous higgs-V couplings

      subroutine width_hgamff(width)

      implicit none

      include "koppln_ew.inc"
      include "para_blocks.inc"

      real*8 width,dum
      real*8 Nc(4),gm(4),matrix1(4),matrix(11)
      real*8 Nf(4)
c     nu_e,e,nu_mu,mu,nu_tau,tau,u,d,c, s, b
c     1    2   3   4    5    6   7 8 9 10 11
      integer i,j,k,l,flavor(11)
       real*8 pi,q2min,q2max
      parameter(pi=3.141592653589793d0)
      logical ldoALL
      data flavor/1,2,1,2,1,2,3,4,3,4,4/
      data ldoALL/.false./

c
c     color factors
      Nc(1) = 1.d0
      Nc(2) = 1.d0
      Nc(3) = 3.d0*(1.d0 + ALFAS/PI)
      Nc(4) = 3.d0*(1.d0 + ALFAS/PI)

c     Number of fermion in generation
      Nf(1) = 3.0d0
      Nf(2) = 3.0d0
      Nf(3) = 2.0d0             !no top
      Nf(4) = 3.0d0

c     compute geometric means of masses
      gm(1) = 0d0  ! NEUTRINO MASS
      gm(2) = (Mf(2,1)*Mf(2,2)*Mf(2,3))**(1.0d0/3.d0)
      gm(3) = (Mf(3,1)*Mf(3,2))**(1.0d0/2.d0)
      gm(4) = (Mf(4,1)*Mf(4,2)*Mf(4,3))**(1.0d0/3.d0)

      dum =0.0d0
      if(ldoALL) then
          do k=1,11
             j = flavor(k)
             call quad1d(4,j,gm(j),xmh,matrix(k))
             dum = matrix(k)*Nc(j) + dum
          enddo
       else
          do j = 1,4
             call quad1d(4,j,gm(j),xmh,matrix1(j))
             dum = matrix1(j)*Nc(j)*Nf(j) + dum
          enddo
       endif
       
      width = dum

      return
      end



********************************************************************************
********************************************************************************
