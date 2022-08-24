c     taken from VBFNLO
c     last modified by Barbara Jaeger: Feb. 2013
C
      subroutine qqHqqj_vonly_channel(nlo,pbar,sign,qbar,gsign,k,ans,ansc)


c >0  includes matrix elements and coupl
      implicit none 
      real * 8 pi,pi2
      parameter (pi=3.141592653589793238462643383279502884197D0,
     1           pi2=pi**2)   
      include 'pwhg_st.h' 
      include 'nlegborn.h' 
      include 'pwhg_flst.h' 
      include 'global_col.inc'

      double precision  fpi
      parameter ( fpi=4d0*pi)
      
      double precision crealH3j, cvirtH3j
      parameter (cvirtH3j = pi**2/3d0-8d0, crealH3j = 2d0*pi**2/3d0-13d0/2d0)
        
c    This subroutine computes the born graphs and finite virtual graphs.
c
C     nlo =  0  |M_born|^2
c     nlo =  1  2Re[M_born conj(M_virt)]
c     nlo =  4  2Re[M_born conj(M_virt-M_virt,boxline)]
c     nlo = -4  2Re[M_born conj(M_virt,boxline)]
c
C  QQHQQ calculates the matrix elements**2 for light Higgs production by
C  electroweak boson fusion in quark quark scattering
C
C        q1 q3    ---->   q2 q4 g H
C
C  QQHQQI must be called first to initialize some couplings
C
C  The main task of QQHQQI is to set up the products of coupling constants
C  needed in Feynman graphs a, ..., g which are stored in
      real*8 fcpl(4,6)
C  and which are calculated from output of KOPPLN
      real*8 clr,xm2,xmg,b,v,a
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),V(4,5),A(4,5)
      include 'koppln_ew.inc'  
C
C  Here fcpl(sig,i) contains the coupling constant factors. 
C     sig = 1,..,4    identifies the 4 different helicity combinations of the
C                     quarks
C     i   = 1,..,6    numbers the possible subprocesses (4 NC and 2 CC)
c
C  
C  The various processes are identified by the following
      integer v2ini(6), v3ini(6)
      integer fl1(6), fl3(6)
      integer sig1(4,2),sig3(4,2)

      data fl1 /3,3,4,4,3,4/
      data fl3 /3,4,3,4,4,3/

      save fl1,fl3
      common /cqqhqqj/ fcpl, sig1,sig3, v2ini,v3ini
c
c alfas, scales etc
      include 'scales.inc'
C  Other local variables for QQBQQI
      integer i,mu,nu,v2,v3, isig,isig1,isig3
      integer imin,nlo
C
C  Variables for the main part of the program
      real*8 ans,ansc(3)
      integer k
      real*8 pbar(0:3,4+nv),qbar(0:4),!uucc,uuss,ddcc,ddss,udsc,ducs,
     1       p(0:3,4+nv),q(0:4),p21(0:4),p43(0:4),p65(0:4),
     2       eps(0:3,0:2),p21g(0:4),p43g(0:4), res(6,2), resv(6,2),
     3       al21, al43
      real*8 p87(0:4),ph(0:4)
c      real*8 uucc_c(3),uuss_c(3),ddcc_c(3),ddss_c(3),udsc_c(3),ducs_c(3)
      integer sign(8), gsign
      complex*16 psi(2,-1:1,4), braket(2,-1:1,4,0:2), 
     1        j21(0:3,-1:1), j43(0:3,-1:1), j21e43(-1:1,-1:1,0:2),
     2        jh1(0:3,-1:1), jh2(0:3,-1:1), e21(0:3,-1:1,0:2),
     3        e43(0:3,-1:1,0:2), e21j43(-1:1,-1:1,0:2)
      complex*16 mm21(6,4,0:2), mm43(6,4,0:2),mv21(6,4,0:2),
     1     mv43(6,4,0:2)
      complex*16 virt(4),born(2),mbox21(-1:1,-1:1,0:2),
     1     mbox43(-1:1,-1:1,0:2)
      complex*16 eps1(0:3,0:2)
      real*8 pk(0:4,4)
      complex*16 prop21(2:4), prop43(2:4), prop21g(2:4), prop43g(2:4)
      real*8 betah, lamb, fac,colfac,ratio
      real*8 den,p57,p68
      real*8 qvec(0:3)          !q_mu from dips subr
      real*8 gm(0:3,0:3)
c
      double precision c2,c2o4pi
      parameter (c2=4d0/3d0, c2o4pi=c2/4d0/pi)
c
c variables for powheg:
      double precision s13,sa3,sa1,sb2,s23,sb3    
      double precision r13,ra3,ra1,rb2,r23,rb3
      double precision l13,la3,la1,lb2,l23,lb3
      double precision l13squared,la3squared,la1squared,lb2squared,
     C     l23squared,lb3squared      
      double precision ca,cf,oo4pi
      parameter (cf=4d0/3d0,ca=3d0,oo4pi=1d0/4d0/pi)

      double precision dotrr
      external dotrr
      double precision tr,nf,gammaq,gammag
      double precision rest12,rest34

      double precision ffunc
      external ffunc

      logical lcol
      logical lborn,lbox,lbox21,lbox43,lvirt,ltri
      logical ldebug 
      parameter(ldebug = .false.)

      logical mg_comp
      parameter (mg_comp=.true.)

      double precision p12, p1g, p2g
      integer ig 


c ====================================================
c
c initialize:
      ans     = 0d0
      ansc(:) = 0d0

C  Reset the coupling factors
      do i = 1,6
         do isig = 1,4
            fcpl(isig,i) = 0d0
         enddo
      enddo

c determine the Yukawa coupling Higgs bb from BR(H--->bb)
c division by 3 takes into account the color factor 3 for H--->qq
c
      do i = 1,6
         do isig = 1,4
            isig1 = sig1(isig,(i+3)/4)
            if ( isig1.ne.0 ) then
               isig3 = sig3(isig,(i+3)/4)
               v2 = v2ini(i)
               v3 = v3ini(i)
               fcpl(isig,i) = clr(fl1(i),v2,isig1)*
     &              clr(fl3(i),v3,isig3)*
     &              b(6,v2,v3)*xmw
            endif
         enddo
      enddo

c              
      imin = 1
      if(nlo.eq.0) then         !born
         lborn = .true.
         lvirt = .false.
         lbox = .false.
         ltri = .false. 
      elseif(nlo.eq.1) then     !box+tri (no Born!)
         lborn = .false.
         lvirt = .true.
         lbox  = .true.
         ltri = .true.
         imin = 0               !do gauge inv.chk.
      elseif(nlo.eq.4) then     !tri (no Born!)
         lborn = .false.
         lvirt = .true.
         lbox = .false.   
         ltri = .true. 
      else !nlo.eq.-4           !box only
         lborn = .false.
         lvirt = .true.
         lbox = .true.
         ltri = .false.
         imin = 0               !do gauge inv.chk.
      endif
c color-correlated Born computed only for LO call:
      if (nlo.eq.0) lcol = .true. 
c
C  Define the internal momenta
      do mu = 0,3
         do i = 1,4+nv
            p(mu,i) = pbar(mu,i)*sign(i)
         enddo
         q(mu) = gsign*qbar(mu)

         p21(mu) = p(mu,1)-p(mu,2)
         p43(mu) = p(mu,3)-p(mu,4)
         p21g(mu) = p21(mu)-q(mu)
         p43g(mu) = p43(mu)-q(mu)
         p65(mu) = p(mu,5)-p(mu,6)
         ph(mu) = p65(mu)
   
      enddo !mu
c
      q(4) = 0
      p21(4) = p21(0)**2-p21(1)**2-p21(2)**2-p21(3)**2
      p43(4) = p43(0)**2-p43(1)**2-p43(2)**2-p43(3)**2
      p21g(4) = p21g(0)**2-p21g(1)**2-p21g(2)**2-p21g(3)**2
      p43g(4) = p43g(0)**2-p43g(1)**2-p43g(2)**2-p43g(3)**2
      p65(4) = p65(0)**2-p65(1)**2-p65(2)**2-p65(3)**2
      ph(4) = ph(0)**2 - ph(1)**2 - ph(2)**2 - ph(3)**2 

c         
C  Get the vector boson propagator factors
c
      if ((p21(4).le.0d0).and.(.not.mg_comp)) then
         prop21(2) = 1/(p21(4)-xm2(2))
         prop21(3) = 1/(p21(4)-xm2(3))     
      else
         prop21(2) = 1/dcmplx(p21(4)-xm2(2),xmg(2))
         prop21(3) = 1/dcmplx(p21(4)-xm2(3),xmg(3))
      endif
      prop21(4) = prop21(3)
c
      if ((p43(4).le.0d0).and.(.not.mg_comp)) then
         prop43(2) = 1/(p43(4)-xm2(2))
         prop43(3) = 1/(p43(4)-xm2(3))
      else
         prop43(2) = 1/dcmplx(p43(4)-xm2(2),xmg(2))
         prop43(3) = 1/dcmplx(p43(4)-xm2(3),xmg(3))
      endif
      prop43(4) = prop43(3)

      if ((p21g(4).le.0d0).and.(.not.mg_comp)) then
         prop21g(2) = 1/(p21g(4)-xm2(2))
         prop21g(3) = 1/(p21g(4)-xm2(3))
      else
         prop21g(2) = 1/dcmplx(p21g(4)-xm2(2),xmg(2))
         prop21g(3) = 1/dcmplx(p21g(4)-xm2(3),xmg(3))
      endif
      prop21g(4) = prop21g(3)

      if ((p43g(4).le.0d0).and.(.not.mg_comp)) then
         prop43g(2) = 1/(p43g(4)-xm2(2))
         prop43g(3) = 1/(p43g(4)-xm2(3))
      else
         prop43g(2) = 1/dcmplx(p43g(4)-xm2(2),xmg(2))
         prop43g(3) = 1/dcmplx(p43g(4)-xm2(3),xmg(3))  
      endif
      prop43g(4) = prop43g(3)
C
C  Get the external spinors (including factor sqrt(2E) )
      call psi0m(4,pbar(0,1),sign(1),psi)

C  Get the gluon polarization vector and the gluon emission spinors
c  or qvec(mu) 
      do i = imin,2
            if(i.eq.0) then
               do mu =0,3
                  eps(mu,i) = q(mu)/q(0) !for gauge inv.check
               enddo
            else 
               call polvec(qbar,i,eps(0,i)) 
            endif 
c         endif
c
c     fill complex*16 eps1(0,i) needed for boxline
         do mu=0,3
            eps1(mu,i) = dcmplx(eps(mu,i),0d0)
         enddo
         do isig = -1,1,2
            call ket2r(psi(1,isig,1),.true.,p(0,1),isig,q,eps(0,i),
     1           braket(1,isig,1,i),pk(0,1))
            call bra2r(psi(1,isig,2),.true.,p(0,2),isig,q,eps(0,i),
     1           braket(1,isig,2,i),pk(0,2))
            call ket2r(psi(1,isig,3),.true.,p(0,3),isig,q,eps(0,i),
     1           braket(1,isig,3,i),pk(0,3))
            call bra2r(psi(1,isig,4),.true.,p(0,4),isig,q,eps(0,i),
     1           braket(1,isig,4,i),pk(0,4))
         enddo
      enddo

C  get the f-fbar currents J21^mu, J43^mu, E21^mu, E43^mu,
c  The calculation of one-loop qcd corrections will be done here.
      call curr(1,psi(1,-1,2),psi(1,-1,1),j21)
      call curr(1,psi(1,-1,4),psi(1,-1,3),j43)
      do i = imin,2
         call curr(1,psi(1,-1,2),braket(1,-1,1,i),jh1)
         call curr(1,braket(1,-1,2,i),psi(1,-1,1),jh2)
         do isig=-1,1,2
            do mu = 0,3
               e21(mu,isig,i) = jh1(mu,isig) + jh2(mu,isig)
            enddo
         enddo

         call curr(1,psi(1,-1,4),braket(1,-1,3,i),jh1)
         call curr(1,braket(1,-1,4,i),psi(1,-1,3),jh2)
         do isig = -1,1,2
            do mu = 0,3
               e43(mu,isig,i) = jh1(mu,isig) + jh2(mu,isig)
            enddo
         enddo
      enddo
c
c     initialize
c
      do i=0,2
         do isig1=-1,1,2
            do isig3=-1,1,2
               mbox43(isig1,isig3,i)=(0d0,0d0)
               mbox21(isig1,isig3,i)=(0d0,0d0)
            enddo
         enddo
      enddo

      
      if (nlo.ne.0) then
      if(lbox) then
         lbox21 = .true.
         lbox43 = .true.
C     GZ suppress coupling of Born gluon to "wrong" fermion line 
C     (particle 6 may be a gluon or a (anti)quark, but tag 1/2 will
C     always correctly tell us which side has the "extra" tag parton)
         ig = flst_borntags(6,flst_cur_iborn) 
C     in case of no tagging, the gluon can be in the initial state
C     so look for a 0 in position 1:2 or 4:6 
         if (any(flst_borntags(1:2,flst_cur_iborn) ==0)) ig = 0  
         if (any(flst_borntags(4:6,flst_cur_iborn) ==0)) ig = 0  
         if (ig .eq.1) then 
            lbox43 = .false.
         elseif (ig .eq.2) then 
            lbox21 = .false.
         elseif (ig .eq.0) then 
!     do nothing in this case (ie keep full Born) 
         else 
            write(*,*) 'ig=',ig 
            write(*,*) 'flst_cur_iborn',flst_cur_iborn
            write(*,*) 'flst_borntag',flst_borntags(:,flst_cur_iborn)
            stop 'qqhqqj-virt: ig should be 0, 1 or 2' 
         endif


         if (gsign.eq.-1) then
            if (sign(1).eq.-sign(2) .and. 
     1          sign(3).eq.sign(4) ) then
                lbox43 = .false. !initial gluon attached to 1-2 line only
            elseif (sign(1).eq.sign(2) .and. 
     1              sign(3).eq.-sign(4) ) then
                lbox21 = .false. !initial gluon attached to 3-4 line only
            endif
         endif
         
         p12=abs(dotrr(p(0,3),p(0,4)))   !no boxes if there are less than 3 jets
         p1g=abs(dotrr(p(0,3),q(0)))     !these terms do not contribute anyway
         p2g=abs(dotrr(p(0,4),q(0)))
         if((p12.lt.1d0).or.(p1g.lt.1d0).or.(p2g.lt.1d0)) then
            lbox43=.false.
         endif

         p12=abs(dotrr(p(0,1),p(0,2)))
         p1g=abs(dotrr(p(0,1),q(0)))
         p2g=abs(dotrr(p(0,2),q(0)))
         if((p12.lt.1d0).or.(p1g.lt.1d0).or.(p2g.lt.1d0)) then
            lbox21=.false.
         endif         
         
c     boxline contributions 
c    
c     1-loop corrections to 21 line with gluon emission
c     First, get the Dij's
         if(lbox21) then
            call bcd_fill_v(p(0,1),p(0,2),q,p21g)
c     
            do i=0,2 ! gluon polarization, 0 is for gauge inv.check
               call curr(1,psi(1,-1,2),braket(1,-1,1,i),jh1) !mborn1
               call curr(1,braket(1,-1,2,i),psi(1,-1,1),jh2) !mborn2
               do isig1 = -1,1,2 ! fermion 1 chirality
                  do isig3 = -1,1,2 ! fermion 3 chirality
c     born-like amplitudes
                     born(1) = jh1(0,isig1)*j43(0,isig3)-
     1                    jh1(1,isig1)*j43(1,isig3)-
     2                    jh1(2,isig1)*j43(2,isig3)-
     3                    jh1(3,isig1)*j43(3,isig3)
                     born(2) = jh2(0,isig1)*j43(0,isig3)-
     1                    jh2(1,isig1)*j43(1,isig3)-
     2                    jh2(2,isig1)*j43(2,isig3)-
     3                    jh2(3,isig1)*j43(3,isig3)
c     
                     call boxline_vg(psi(1,isig1,1),psi(1,isig1,2),
     1                    p(0,1),p(0,2),isig1,
     2                    eps1(0,i),j43(0,isig3),q,p21g,
     3                    born,virt)
c
                     mbox21(isig1,isig3,i) = virt(4)
                  enddo !isig3
               enddo !isig1
            enddo ! gluon pol (i)
         endif !lbox21

c     1-loop corrections to 43 line with gluon emission
c     First, get the Dij's
         if(lbox43) then
            call bcd_fill_v(p(0,3),p(0,4),q,p43g)

            do i=0,2               ! gluon polarizations
               call curr(1,psi(1,-1,4),braket(1,-1,3,i),jh1) !mborn1
               call curr(1,braket(1,-1,4,i),psi(1,-1,3),jh2) !mborn2
               do isig1 = -1,1,2    ! fermion 1 chirality
                  do isig3 = -1,1,2 ! fermion 3 chirality
c     born-like amplitudes
                     born(1) = jh1(0,isig3)*j21(0,isig1)-
     1                 jh1(1,isig3)*j21(1,isig1)-jh1(2,isig3)*j21(2,isig1)-
     2                 jh1(3,isig3)*j21(3,isig1)
                     born(2) = jh2(0,isig3)*j21(0,isig1)-
     1                 jh2(1,isig3)*j21(1,isig1)-jh2(2,isig3)*j21(2,isig1)-
     2                 jh2(3,isig3)*j21(3,isig1)

                     call boxline_vg(psi(1,isig3,3),psi(1,isig3,4),
     1                    p(0,3),p(0,4),isig3,
     2                    eps1(0,i),j21(0,isig1),q,p43g,
     3                    born,virt)
                     mbox43(isig1,isig3,i) = virt(4)

                  enddo !isig3
               enddo !isig1
            enddo !gluon pol (i)
         endif !lbox43

      endif !lbox             
      endif !nlo.ne.0

C  get the dot products of the currents for the 4 helicity combinations
      do isig1 = -1,1,2
         do isig3 = -1,1,2
            do i = imin,2
               e21j43(isig1,isig3,i) = e21(0,isig1,i)*j43(0,isig3)-
     1                                 e21(1,isig1,i)*j43(1,isig3)-
     2                                 e21(2,isig1,i)*j43(2,isig3)-
     3                                 e21(3,isig1,i)*j43(3,isig3)
               j21e43(isig1,isig3,i) = e43(0,isig3,i)*j21(0,isig1)-
     1                                 e43(1,isig3,i)*j21(1,isig1)-
     2                                 e43(2,isig3,i)*j21(2,isig1)-
     3                                 e43(3,isig3,i)*j21(3,isig1)
     
            enddo
         enddo
      enddo

C
C  now get the coupling*propagator factors for subprocess k, helicity
C  combination isig
c     Done for both born and virtual graphs
c
      colfac = cf               !color factors 
c      do k = 1,6
         do isig = 1,4
            isig1 = sig1(isig,(k+3)/4)
            if (isig1.ne.0) then
               isig3 = sig3(isig,(k+3)/4)
               v2 = v2ini(k)
               v3 = v3ini(k)
               do i = imin,2
c
                  mv21(k,isig,i)=0.0d0
                  mv43(k,isig,i)=0.0d0
c     born: 
                  mm21(k,isig,i) = fcpl(isig,k)*prop21g(v2)*
     $                 prop43(v3)*e21j43(isig1,isig3,i) !born term 21 line
c     
                  mm43(k,isig,i) = fcpl(isig,k)*prop21(v2)*
     $                 prop43g(v3)*j21e43(isig1,isig3,i) !born term 43 line
c     
c     boxes:
                  if (nlo.ne.0) then
c
                  if(lbox21) then !check that the correct factors are here
                     mv21(k,isig,i) = fcpl(isig,k)*prop21g(v2)*
     1                    prop43(v3)*mbox21(isig1,isig3,i) +
     2                    colfac*cvirtH3j*mm21(k,isig,i)     !NOT NEEDED HERE, ADD ALL TOGETHER AT THE END

                     if(i.gt.0) then !gauge invariance testing
                        ratio = abs(mv21(k,isig,0)/mm21(k,isig,i))
                        if(ratio.gt.0.1d0) then ! use constant factor
                           mv21(k,isig,i) =colfac*cvirtH3j*mm21(k,isig,i) !GETS ADDED BACK LATER
                        endif
                     endif
                  endif !lbox21
                  if(lbox43) then
                     mv43(k,isig,i) = fcpl(isig,k)*prop21(v2)*
     1                    prop43g(v3)*mbox43(isig1,isig3,i)
     2                   +colfac*cvirtH3j*mm43(k,isig,i)     
ccc      


                     if(i.gt.0) then ! gauge invariance testing
                        ratio = abs(mv43(k,isig,0)/mm43(k,isig,i))
                        if(ratio.gt.0.1d0) then ! use constant factor      
                           mv43(k,isig,i) =colfac*cvirtH3j*mm43(k,isig,i)
                        endif
                     endif
                  endif !lbox43
c     triangles:
                  if(ltri) then
                     mv21(k,isig,i) = mv21(k,isig,i)  
     1                   + cf*cvirtH3j*mm21(k,isig,i) !21 line 43 tri
                     mv43(k,isig,i) = mv43(k,isig,i)
     1                   + cf*cvirtH3j*mm43(k,isig,i) !41 line 21 tri
                  endif !ltri
                  endif !nlo.ne.0

               enddo !i
            else !isig1
               do i = imin,2
                  mm21(k,isig,i) = 0
                  mm43(k,isig,i) = 0
                  mv21(k,isig,i) = 0
                  mv43(k,isig,i) = 0                
               enddo
               
            endif !isig1
         enddo !isig
c      enddo !k
C
C  Now sum the contributions from the 4 helicity combinations for all 
C  subprocesses and multiply by the factor corresponding to
C
C
      fac = 1d0 
      fac = fac * 12         ! 12 = C_2*9 is the color factor
c     
c for initial gluon, i.e. gsign=-1, 
c eliminate diagrams with V-->q qbar decay
c
c     All alpha_s have the same renormalization scale!
c
      al21 = fpi*als(1,1)               
      al43 = fpi*als(2,1)               

      if (gsign.eq.-1) then
         if (sign(1).eq.-sign(2) .and. 
     1       sign(3).eq.sign(4) ) then
            al43 = 0            !  initial gluon attached to 1-2 line only
         elseif (sign(1).eq.sign(2) .and. 
     1           sign(3).eq.-sign(4) ) then
            al21 = 0            !  initial gluon attached to 3-4 line only
         endif
      endif

C     GZ suppress coupling of Born gluon to "wrong" fermion line 
C     (particle 6 may be a gluon or a (anti)quark, but tag 1/2 will
C     always correctly tell us which side has the "extra" tag parton)
      ig = flst_borntags(6,flst_cur_iborn) 
C     in case of no tagging, the gluon can be in the initial state
C     so look for a 0 in position 1:2 or 4:6 
      if (any(flst_borntags(1:2,flst_cur_iborn) ==0)) ig = 0  
      if (any(flst_borntags(4:6,flst_cur_iborn) ==0)) ig = 0  
      if (ig .eq.1) then 
         al43 = 0d0 
      elseif (ig .eq.2) then 
         al21 = 0d0 
      elseif (ig .eq.0) then 
!     do nothing in this case (ie keep full Born) 
      else 
         write(*,*) 'ig=',ig 
         write(*,*) 'flst_cur_iborn',flst_cur_iborn
         write(*,*) 'flst_borntag',flst_borntags(:,flst_cur_iborn)
         stop 'qqhqqj-virt: ig should be 0, 1 or 2' 
      endif


c
c     sum over helicity combinations
c     
c      do k = 1,6
         res(k,1) = 0
         res(k,2) = 0
         resv(k,1) = 0
         resv(k,2) = 0

ccccccccccccccccc

         do isig = 1,4
            do i = 1,2   !imin,2  imin=0 -> gauge check

               res(k,1) = res(k,1)  
     &                   + al21*(dreal(mm21(k,isig,i))**2
     &                   +       dimag(mm21(k,isig,i))**2)
               res(k,2) = res(k,2)  
     &                   + al43*(dreal(mm43(k,isig,i))**2
     &                   +       dimag(mm43(k,isig,i))**2)

               if (nlo.ne.0) then

c  add Born type term and multiply by F_q = alphas/4pi;
c  the factor for the born term is 
c  after adding the subtraction term
c  and the counter term for the renormalization of the pdfs

c
c     The factor  (4*Pi)^ep/Gamma(1-ep) IS NOT RETURNED by this subroutine
c     and it's thought as factorized in front of the real counterterms too.
c
c
cccccccc
c
c VBF_H3j:
c
C virtual contributions (2.29) in arXiv:0710.5621: 
c
c factoring out common scale factor gives rise 
c to extra terms of the form
c
c  (st_muren2/sij)^eps = (ratio)^eps 
c   ~  1+eps log(ratio)+1/2 eps^2 log^2(ratio)
c
c -> corrections to upper line:
c
c for q(pa) Q(pb) -> q(p1) Q(p2) g(p3) H
c

!          do i = 1,4+nv
!             p(mu,i) = pbar(mu,i)*sign(i)
!          enddo
!          q(mu) = gsign*qbar(mu)

c s13 -> 2 p2.q, sa3 -> 2 p1.q, sa1 -> 2 p1.p2, sb2 -> 2 p3.p4
c
                  s13 = 2d0*dotrr(pbar(0,2),qbar(0))
                  sa3 = 2d0*dotrr(pbar(0,1),qbar(0))
                  sa1 = 2d0*dotrr(pbar(0,1),pbar(0,2))
                  sb2 = 2d0*dotrr(pbar(0,3),pbar(0,4))

c corrections to lower line:
c
c s23 -> 2 p4.q, sb3 -> 2 p3.q, sb2 -> 2 p3.p4, sa1 -> 2 p1.p2
c
                  s23 = 2d0*dotrr(pbar(0,4),qbar(0))
                  sb3 = 2d0*dotrr(pbar(0,3),qbar(0))
                  
c
c build ratios of invariants to scale:    
                  r13 = abs(st_muren2/s13) 
                  l13 = log(r13)
                  if(s13.lt.0d0) then
                    l13squared=l13**2-pi**2
                    print*, 's13<0, should not happen: ', s13  
                  else
                    l13squared=l13**2
                  endif

                  ra3 = abs(st_muren2/sa3) 
                  la3 = log(ra3)
                  if(sa3.lt.0d0) then
                    la3squared=la3**2-pi**2
                    print*, 'sa3<0, should not happen: ', sa3  
                  else
                    la3squared=la3**2
                  endif

                  ra1 = abs(st_muren2/sa1) 
                  la1 = log(ra1)
                  if(sa1.lt.0d0) then
                    la1squared=la1**2-pi**2
                    print*, 'sa1<0, should not happen: ', sa1  
                  else
                    la1squared=la1**2
                  endif

                  rb2 = abs(st_muren2/sb2) 
                  lb2 = log(rb2)
                  if(sb2.lt.0d0) then
                    lb2squared=lb2**2-pi**2
                    print*, 'sb2<0, should not happen: ', sb2  
                  else
                    lb2squared=lb2**2
                  endif

                  r23 = abs(st_muren2/s23)
                  l23 = log(r23)
                  if(s23.lt.0d0) then
                    l23squared=l23**2-pi**2
                    print*, 's23<0, should not happen: ', s23  
                  else
                    l23squared=l23**2
                  endif

                  rb3 = abs(st_muren2/sb3)
                  lb3 = log(rb3)
                  if(sb3.lt.0d0) then
                    lb3squared=lb3**2-pi**2
                    print*, 'sb3<0, should not happen: ', sb3  
                  else
                    lb3squared=lb3**2
                  endif

          tr = 0.5d0
!          nf = 4d0
          nf = dble(st_nlight)
          gammag = 11d0/6d0*ca-2d0/3d0*tr*nf   
          gammaq = 1.5d0*cf

c BJ: need to adapt constant terms in restij!
          rest12 = -1d0/6d0*pi**2*(ca+4d0*cf) +   !-1d0/6d0*pi**2*ca-16d0*cf+
     &         ffunc(sa1,sa3,s13,st_muren2)  
     

          rest34 =  -1d0/6d0*pi**2*(ca+4d0*cf) +   !-1d0/6d0*pi**2*ca-16d0*cf+
     &         ffunc(sb2,sb3,s23,st_muren2)   

                  if (nlo.gt.0) then ! virtuals +log terms 

                     mv21(k,isig,i) = als(1,1)*oo4pi*     
     &                ( mv21(k,isig,i) + 
     &                  mm21(k,isig,i)*(
     &                  0.5d0*(-ca*0.5d0*l13squared-gammag*l13) 
     &                 +0.5d0*(-ca*0.5d0*la3squared-gammag*la3) 
     &                 +0.5d0*ca/cf*(-cf*0.5d0*l13squared-gammaq*l13)  
     &                 +0.5d0*ca/cf*(-cf*0.5d0*la3squared-gammaq*la3)  
     &                 -ca/cf*(-cf*0.5d0*la1squared-gammaq*la1) 
     &                 +2d0*(-cf*0.5d0*la1squared-gammaq*la1)  
     &                 +2d0*(-cf*0.5d0*lb2squared-gammaq*lb2) 
     &                  + rest12 
     &                    ) )

                     mv43(k,isig,i) = als(2,1)*oo4pi*    !al43*
     &                ( mv43(k,isig,i) + 
     &                  mm43(k,isig,i)*(
     &                  0.5d0*(-ca*0.5d0*l23squared-gammag*l23) 
     &                 +0.5d0*(-ca*0.5d0*lb3squared-gammag*lb3) 
     &                 +0.5d0*ca/cf*(-cf*0.5d0*l23squared-gammaq*l23)  
     &                 +0.5d0*ca/cf*(-cf*0.5d0*lb3squared-gammaq*lb3)  
     &                 -ca/cf*(-cf*0.5d0*lb2squared-gammaq*lb2) 
     &                 +2d0*(-cf*0.5d0*lb2squared-gammaq*lb2)  
     &                 +2d0*(-cf*0.5d0*la1squared-gammaq*la1) 
     &                  + rest34
     &                    ) )

                  else !virtuals only
                     mv21(k,isig,i) = 
     1                    als(1,1)*oo4pi*mv21(k,isig,i)
                     mv43(k,isig,i) = 
     1                    als(2,1)*oo4pi*mv43(k,isig,i)
                  endif
                  resv(k,1)= resv(k,1) + al21*2d0*dreal(
     1                 mm21(k,isig,i)   *
     1              conjg( mv21(k,isig,i) )  )
     
                  resv(k,2) = resv(k,2) + al43*2d0*dreal(
     1                 mm43(k,isig,i)   *
     1              conjg( mv43(k,isig,i) )  )
               endif !lnlo
            enddo !isig3
         enddo  !isig1

ccccccccccccc
         
c      enddo !k
      
c      do k =1,6
c     born**2
         if(lborn) then
            res(k,1) = res(k,1)*fac !21 line
            res(k,2) = res(k,2)*fac !43 line
         endif
c     2Re[mborn mvirt*]
         if (nlo.ne.0) then  
         if(lvirt) then
            resv(k,1) = resv(k,1)*fac !21 line
            resv(k,2) = resv(k,2)*fac !43 line 
         endif
         endif ! nlo

c      enddo !k      
c     
c    

      if(nlo.ne.0d0) then
          ansc(2) =  resv(k,1)
          ansc(3) =  resv(k,2)
          ansc(1) =  resv(k,1) + resv(k,2)
          ans = resv(k,1)+resv(k,2)
      else
          ansc(2) =  resv(k,1) + res(k,1)
          ansc(3) =  resv(k,2) + res(k,2)
          ansc(1) =  resv(k,1) + resv(k,2) + res(k,1) + res(k,2)
          ans = resv(k,1)+resv(k,2)+ res(k,1) + res(k,2)      
      endif
      
     
      
      
      return
      end
c------------------------------------------------------------------------

      real * 8 function ffunc(s,t,u,mu2)
      
      implicit none

      real * 8 pi,pi2
      parameter (pi=3.141592653589793238462643383279502884197D0,
     1           pi2=pi**2)         
      
      real*8 s,t,u,mu2
      real*8 ca,cf,tr,nf
      parameter (ca=3d0,cf=4d0/3d0,tr=0.5d0,nf=4d0)

      real*8 ratu,ratt,rats
      real*8 logu,logt,logs,logu2,logt2,logs2

      ratu = abs(u/mu2)
      ratt = abs(t/mu2)
      rats = abs(s/mu2)

      logu = log(ratu)
      logt = log(ratt)
      logs = log(rats)
      

      
      if(s.lt.0d0) then
        logs2=logs**2-pi**2
        print*, 's<0, should not happen: ', s
      else
        logs2=logs**2
      endif

 
      
      if(t.lt.0d0) then
        logt2=logt**2-pi**2
        print*, 't<0, should not happen: ', t        
      else
        logt2=logt**2
      endif      

      if(u.lt.0d0) then
        logu2=logu**2-pi**2
        print*, 'u<0, should not happen: ', u          
      else
        logu2=logu**2
      endif
      
      ffunc = 0.5d0*ca*(logu2+logt2)-0.5d0*(ca-2d0*cf)*logs2+
     &       1.5d0*(ca-2d0*cf)*logs+
     &       (tr*nf/3d0-5d0*ca/3d0)*(logu+logt)



      end
