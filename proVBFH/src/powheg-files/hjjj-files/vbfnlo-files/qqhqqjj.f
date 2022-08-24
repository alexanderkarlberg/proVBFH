c     adapted from VBFNLO
c     last modified by Barbara Jaeger: Feb. 2013
C
c     All q(pbar1)+ q(pbar3) -> q(pbar2)+ q(pbar4) 
c     +g(qbar1)+ g(qbar2) + h
c                           |
c                           ---> pbar5 pbar6 
c     Decay products of the higgs are stored in 
c     pbar(0:3,5) and pbar(0:3,6)
c     

      subroutine qqHqqjj_channel(pbar,sign,qbar,gsign,k,ans)
c
c     ans gives color-sum, ansc gives individual color structures
c     (use switch lcol to choose between those)
c     probably don't need ansc at all.
c
      implicit none     
      real * 8 pi,pi2
      parameter (pi=3.141592653589793238462643383279502884197D0,
     1           pi2=pi**2)   
      include '../../include/pwhg_st.h'  
      include 'global_col.inc'          
      double precision  fpi
      parameter ( fpi=4d0*pi)
c     switch for interference  terms
      logical lintOff
      common/interference/ lintOff
c
c     Arguments
      real*8 pbar(0:3,4+nv),qbar(0:4,2)
      real*8 ans,ansc(6)
      integer sign(4+nv),gsign(2)
c and which are calculated from output of KOPPLN
c
      double precision  clr, xm2, xmg, b, v, a
      COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),
     1                  V(4,5),A(4,5)
      include 'koppln_ew.inc'
c     
c     Local variables
      real*8 ph(0:4),lpbar(0:3,4),lqbar(0:3,2),resnc(4,0:6),rescc(0:6)
      real*8 res(6,0:6)
      integer lfsign(4),lgsign(2)
      real*8 p65(0:4),p(0:3,4+nv),fac
c      real*8 p57,p68,den,betah,lamb
      integer mu,i,k,kl
      
      logical lcol
      parameter (lcol = .false.) ! compute color sum
c      parameter (lcol = .true.) ! compute individual color structures

c ====================================================
c
c initialize:
      ans     = 0d0
      ansc(:) = 0d0

c     switch off interference terms
c
      lintOff = .true.
c
c     initialize res(k)
      do i = 0,6
         do kl=1,6
            res(kl,i)=0d0
         enddo
         do kl=1,4
            resnc(kl,i)=0d0
         enddo
         rescc(i) = 0d0
      enddo
c     fill local momenta and sign factors
c     diagramatic momenta
      do mu = 0,3
         do i = 1,4+nv
            p(mu,i) = pbar(mu,i)*sign(i)
         enddo         
         p65(mu) = p(mu,5) - p(mu,6)
         ph(mu) = p65(mu)
      enddo
      p65(4) = p65(0)**2 - p65(1)**2 - p65(2)**2 - p65(3)**2
      ph(4) = ph(0)**2 - ph(1)**2 - ph(2)**2 - ph(3)**2
c
      do i=1,4
         lfsign(i)=sign(i)
         do mu=0,3
            lpbar(mu,i)=pbar(mu,i)
         enddo
      enddo
      do i=1,2
         lgsign(i)=gsign(i)
         do mu=0,3
            lqbar(mu,i)=qbar(mu,i)
         enddo
      enddo

c
c     call madgraph routines:
      if (k.le.4) then 
         call SQQ_QQGGH_NC(lpbar,lfsign,lqbar,lgsign,ph,k,resnc) !NC
      else   
         call SQQ_QQggH_CC(lpbar,lfsign,lqbar,lgsign,ph,rescc) !CC
      endif
c
c    do not include higgs decay:
      fac = 1d0
      
      do i=0,6                  !loop over color flows
         if (k.le.4) then !NC
c         do k=1,4               !loop over subprocesses
            res(k,i) = resnc(k,i)*fac
c         enddo
         else ! k=5 or 6 ... CC   
            res(k,i) = rescc(i)*fac
         endif
c         res(6,i) = rescc(i)*fac
      enddo

c     return matrix elements squared:
      if(lcol) then
         do i = 0,6
            ansc(i) = res(k,i)
         enddo
      else
         ans = res(k,0)
      endif
      return
      end
