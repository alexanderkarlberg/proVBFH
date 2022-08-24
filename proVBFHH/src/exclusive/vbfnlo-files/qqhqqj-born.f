c     taken from VBFNLO
c     last modified by Barbara Jaeger: Jan. 2013
C
      subroutine qqHqqj_born_channel(nlo,pbar,sign,qbar,gsign,k,ans,ansc,bmunu)
c
      implicit none 
      real * 8 pi,pi2
      parameter (pi=3.141592653589793238462643383279502884197D0,
     1           pi2=pi**2)   
      include 'pwhg_st.h' 
      include 'global_col.inc'
      include 'nlegborn.h' 
      include 'pwhg_flst.h' 
      include 'pwhg_flst_2.h' 
      include 'tags.h' 
      include 'leptens.h'
      double precision  fpi
      parameter ( fpi=4d0*pi)
c
C  calculate the LO matrix elements**2 for light Higgs production by
C  electroweak boson fusion in quark quark scattering
C
C        q1 q3    ---->   q2 q4 g H
c
      real*8 fcpl(4,6)
      real*8 clr,xm2,xmg,b,v,a
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),V(4,5),A(4,5)
      include 'koppln_ew.inc'  
C
C  Here fcpl(sig,i) contains the coupling constant factors. 
C     sig = 1,..,4    identifies the 4 different helicity combinations of the
C                     quarks
C     i   = 1,..,6    numbers the possible subprocesses (4 NC and 2 CC)
C  
C  The various processes are identified by the following
      integer v2ini(6), v3ini(6)
      integer fl1(6), fl3(6)
      integer sig1(4,2),sig3(4,2)

      data fl1 /3,3,4,4,3,4/
      data fl3 /3,4,3,4,4,3/

      save fl1,fl3
      common /cqqhqqj/ fcpl, sig1,sig3, v2ini,v3ini
      data v2ini /2,2,2,2,3,4/, v3ini /2,2,2,2,4,3/
      data sig1 /-1,-1, 1, 1,-1,3*0/
      data sig3 /-1, 1,-1, 1,-1,3*0/
c
c alfas, scales etc
      include 'scales.inc'
      include 'color.inc'
C  Other local variables for QQBQQI
      integer i,mu,nu,v2,v3, isig,isig1,isig3
      integer imin,nlo, pol1, pol2
C
C  Variables for the main part of the program
      real*8 ans,ansc(3)
      integer k
      real*8 pbar(0:3,4+nv),qbar(0:4),!uucc,uuss,ddcc,ddss,udsc,ducs,
     1       p(0:3,4+nv),q(0:4),
     2       eps(0:3,0:2),res(6,2), 
     3       l21, l43
      real*8 ph(0:4), bmunu(0:3,0:3)
      complex*16 complexbmunu(0:3,0:3)
      integer sign(4+nv), gsign
      complex*16 psi(2,-1:1,4), braket(2,-1:1,4,0:2), 
     1        j21(0:3,-1:1), j43(0:3,-1:1), j21e43(-1:1,-1:1,0:2,2),
     1        j21h(0:3,-1:1,2), j43h(0:3,-1:1,2),
     2        jh1(0:3,-1:1), jh2(0:3,-1:1), e21(0:3,-1:1,0:2),
     3        e43(0:3,-1:1,0:2), e21j43(-1:1,-1:1,0:2,2)
      complex*16 mm21(6,4,0:2), mm43(6,4,0:2)
      complex*16 born(2)
      complex*16 eps1(0:3,0:2)
      real*8 pk(0:4,4)
      real*8 fac,colfac

      logical lcol
      logical lborn,lbox,lvirt,ltri
      logical ldebug 
      parameter(ldebug = .false.)

      logical mg_comp
      parameter (mg_comp=.true.)
      integer ig,ig2 

      double complex contract_Tjj
      external contract_Tjj
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
     &              clr(fl3(i),v3,isig3)
            endif
         enddo
      enddo

      imin = 1
      if(nlo.eq.0) then         !born
         lborn = .true.
         lvirt = .false.
         lbox = .false.
         ltri = .false. 
      endif
c color-correlated Born computed only for LO call:
      if (nlo.eq.0) lcol = .true. 
c      
c
C  Define the internal momenta
      do mu = 0,3
         do i = 1,4+nv
            p(mu,i) = pbar(mu,i)*sign(i)
         enddo
         q(mu) = gsign*qbar(mu)
      enddo !mu

C  Get the external spinors (including factor sqrt(2E) )
      call psi0m(4,pbar(0,1),sign(1),psi)

C  Get the gluon polarization vector and the gluon emission spinors
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
c        fill complex*16 eps1(0,i) needed for boxline
         do mu=0,3
            eps1(mu,i) = dcmplx(eps(mu,i),0d0)
         enddo
c
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

      do isig = -1,1,2
        do mu=0,3
          j21h(mu,isig,(k+3)/4) =
     &      vvhh(0,mu,2,(k+3)/4)*j21(0,isig)
     &     -vvhh(1,mu,2,(k+3)/4)*j21(1,isig)
     &     -vvhh(2,mu,2,(k+3)/4)*j21(2,isig)
     &     -vvhh(3,mu,2,(k+3)/4)*j21(3,isig)
          j43h(mu,isig,(k+3)/4) =
     &      vvhh(mu,0,1,(k+3)/4)*j43(0,isig)
     &     -vvhh(mu,1,1,(k+3)/4)*j43(1,isig)
     &     -vvhh(mu,2,1,(k+3)/4)*j43(2,isig)
     &     -vvhh(mu,3,1,(k+3)/4)*j43(3,isig)
        enddo
      enddo
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
            enddo !mu
         enddo !isig
      enddo !i
C  get the dot products of the currents for the 4 helicity combinations
      do isig1 = -1,1,2
         do isig3 = -1,1,2
            do i = imin,2
              e21j43(isig1,isig3,i,(k+3)/4) = e21(0,isig1,i)*j43h(0,isig3,(k+3)/4)-
     1                                  e21(1,isig1,i)*j43h(1,isig3,(k+3)/4)-
     2                                  e21(2,isig1,i)*j43h(2,isig3,(k+3)/4)-
     3                                  e21(3,isig1,i)*j43h(3,isig3,(k+3)/4)
              j21e43(isig1,isig3,i,(k+3)/4) = e43(0,isig3,i)*j21h(0,isig1,(k+3)/4)-
     1                                  e43(1,isig3,i)*j21h(1,isig1,(k+3)/4)-
     2                                  e43(2,isig3,i)*j21h(2,isig1,(k+3)/4)-
     3                                  e43(3,isig3,i)*j21h(3,isig1,(k+3)/4)
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
c                born:
                  mm21(k,isig,i) = fcpl(isig,k)*e21j43(isig1,isig3,i,(k+3)/4) !born term 21 line
c     
                  mm43(k,isig,i) = fcpl(isig,k)*j21e43(isig1,isig3,i,(k+3)/4) !born term 43 line

               enddo !i
            else !isig1
               do i = 1,2
                  mm21(k,isig,i) = 0d0
                  mm43(k,isig,i) = 0d0
               enddo
            endif !isig1
         enddo !isig
c      enddo !k
C
C  Now sum the contributions from the 4 helicity combinations for all 
C  subprocesses and multiply by the factor corresponding to
C
      fac = 1d0   

      fac = fac * 12d0         ! 12 = C_2*9 is the color factor

c for initial gluon, i.e. gsign=-1, eliminate diagrams with V-->q qbar decay
c     All alpha_s have the same renormalization scale!
c
      l21 = fpi*als(1,1)                 !*als(gnlo) is in dipolesub.f
      l43 = fpi*als(2,1)                 !*als(gnlo) is in dipolesub.f
      
!       print*, als(1,2)
c     
      if (gsign.eq.-1) then
         if (sign(1).eq.-sign(2) .and. 
     1       sign(3).eq.sign(4) ) then
            l43 = 0d0            !  initial gluon attached to 1-2 line only
         elseif (sign(1).eq.sign(2) .and. 
     1           sign(3).eq.-sign(4) ) then
            l21 = 0d0            !  initial gluon attached to 3-4 line only
         endif
      endif
!     VBFHHMOD
!      ig = flst_borntags(6,flst_cur_iborn) 
      ig = flst_borntags(nlegborn,flst_cur_iborn) 
C     GZ sanity check: that performing physto_diag permutation of 
C     tags gives the same answer 
C     i.e. quarks with tag 1 (2) are always kept on line l21 (l43)
      ig2 = flst_borntags_todiag(nlegborn) 
      if (ig /= ig2 .and. ig*ig2 /=0) then 
         write(*,*) flst_borntags(:,flst_cur_iborn) 
         write(*,*) flst_borntags_todiag(:) 
         write(*,*) 'qqhqqj-born: ig and ig2 do not match', ig, ig2
         stop 
      endif

C     GZ suppress coupling of Born gluon to "wrong" fermion line 
C     (particle 6 may be a gluon or a (anti)quark, but tag 1/2 will
C     always correctly tell us which side has the "extra" parton)
      if (ig .eq.1) then 
         l43 = 0d0 
      elseif (ig .eq.2) then 
         l21 = 0d0 
      elseif (ig .eq.0) then 
!     do nothing in this case (ie keep full Born) 
      else 
         write(*,*) 'ig=',ig 
         write(*,*) 'flst_cur_iborn',flst_cur_iborn
         write(*,*) 'flst_borntag',flst_borntags(:,flst_cur_iborn)
         stop 'qqhqqj-born: ig should be 0, 1 or 2' 
      endif

c      do k = 1,6
         res(k,1) = 0d0
         res(k,2) = 0d0

         do isig = 1,4
            do i = imin,2

               res(k,1) = res(k,1)  
     &                   + l21*(dreal(mm21(k,isig,i))**2
     &                   +      dimag(mm21(k,isig,i))**2)
               res(k,2) = res(k,2)  
     &                   + l43*(dreal(mm43(k,isig,i))**2
     &                   +      dimag(mm43(k,isig,i))**2)
             enddo
          enddo     
c      enddo !k

         bmunu(:,:)=0d0
         complexbmunu(:,:)=(0d0,0d0)
         do mu=0,3
         do nu=0,3
         do pol1 = imin,2   !polarization 
         do pol2 = imin,2
         do isig = 1,4   ! 4 helicity combinations


            complexbmunu(mu,nu) =  complexbmunu(mu,nu) + 
     $                 fac*l21*( (mm21(k,isig,pol1))*
     $                 dconjg(mm21(k,isig,pol2)))*
     $                 dconjg(eps1(mu,pol1))*(eps1(nu,pol2) )                              !|mm21|^2     
     $                 + 
     $                 fac*l43*( (mm43(k,isig,pol1))*
     $                 dconjg(mm43(k,isig,pol2)))*
     $                 dconjg(eps1(mu,pol1))*(eps1(nu,pol2) )  
     
            enddo
            enddo !i
            
         enddo ! isig    

         enddo
         enddo      


         if(lborn) then
            res(k,1) = res(k,1)*fac !21 line
            res(k,2) = res(k,2)*fac !43 line
         endif
c      enddo !k
      bmunu(:,:) = dreal(complexbmunu(:,:))
      ansc(2) = res(k,1) 
      ansc(3) = res(k,2) 
      ansc(1) = res(k,1) + res(k,2) 

      ans = res(k,1)+res(k,2)

      
      return
      end
c------------------------------------------------------------------------

      subroutine calcleptens(p)
      implicit none
      
      include 'nlegborn.h' 
      include 'leptens.h'
      include "tensor.inc"

      double precision p(0:3,nlegborn)
      double precision q12(0:3), q12g(0:3), q34(0:3), q34g(0:3)

c momentum setup
      q12 = p(:,1)-p(:,5)
      q12g = q12-p(:,nlegborn)
      q34 = p(:,2)-p(:,6)
      q34g = q34-p(:,nlegborn)

c get leptonic tensors for Higgs pair production
      CALL SXXXXX(p(0,3),1,wp)
      CALL SXXXXX(p(0,4),1,wm)

      call ZZtoHH(q12g,q34,p(0,3),vvhh(0,0,1,1))
      call ZZtoHH(q12,q34g,p(0,3),vvhh(0,0,2,1))
      call WWtoHH(q12g,q34,p(0,3),vvhh(0,0,1,2))
      call WWtoHH(q12,q34g,p(0,3),vvhh(0,0,2,2))

      return
      end

