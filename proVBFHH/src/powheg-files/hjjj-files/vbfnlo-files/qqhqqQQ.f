c
c     adapted from VBFNLO
c     last modified by Barbara Jaeger: Feb. 2013
C
c     All q(pbar1)+ q(pbar3) -> q(pbar2)+ q(pbar4) 
c     +q(qbar5)+ qbar(qbar6) + h
c                           |
c                           ---> pbar7 pbar8 
c     Decay products of the higgs are stored in 
c     pbar(0:3,7) and pbar(0:3,8)
c     
c     Neutral current weak boson fusion processes
c
      subroutine qqHqqQQ_nc_channel(pbar,sign,k,ans)

      implicit none     
      real * 8 pi,pi2
      parameter (pi=3.141592653589793238462643383279502884197D0,
     1           pi2=pi**2)   
      include '../../include/pwhg_st.h' 
      include 'global_col.inc'                         
      double precision  fpi
      parameter ( fpi=4d0*pi)
c
c     Arguments
      real*8 pbar(0:3,6+nv)
      real*8 ans,ansc(0:2)
      integer sign(6+nv)
c
c and which are calculated from output of KOPPLN
c
      double precision  clr, xm2, xmg, b, v, a
      COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),
     1                  V(4,5),A(4,5)
      include 'koppln_ew.inc'
c     
c     Local variables
      real*8 ph(0:4),lpbar(0:3,6),res(0:2,8),msq21,msq43
      integer lfsign(6)
      integer iflav(6,8)
      real*8 p87(0:4),p(0:3,6+nv),fac
      integer mu,i,k,kl
      
      logical lcol
      parameter (lcol = .false.) ! compute color sum
c      parameter (lcol = .true.) ! compute individual color structures

c     flavors of different subprocesses (1=up,2=down!)
c
      DATA(iflav(i,1),i=1,6)/1,1,1,1,1,1/ ! uuuutt (k=1)
      DATA(iflav(i,5),i=1,6)/1,1,1,1,2,2/ ! uuuubb (k=5)
      DATA(iflav(i,2),i=1,6)/1,1,2,2,1,1/ ! uusstt (k=2)
      DATA(iflav(i,6),i=1,6)/1,1,2,2,2,2/ ! uussbb (k=6)
      DATA(iflav(i,3),i=1,6)/2,2,1,1,1,1/ ! ddcctt (k=3)
      DATA(iflav(i,7),i=1,6)/2,2,1,1,2,2/ ! ddccbb (k=7)
      DATA(iflav(i,4),i=1,6)/2,2,2,2,1,1/ ! ddsstt (k=4)
      DATA(iflav(i,8),i=1,6)/2,2,2,2,2,2/ ! ddssbb (k=8)

c ====================================================
c
c initialize:
      ans     = 0d0
      ansc(:) = 0d0
c
c     initialize res(k)
      do i = 0,2
         do kl=1,8
            res(i,kl)=0d0
         enddo
      enddo
c     fill local momenta and sign factors
c     diagramatic momenta
      do mu=0,3
         do i=1,6+nv
            p(mu,i)=sign(i)*pbar(mu,i)
         enddo
         p87(mu) = p(mu,8) - p(mu,7)
         ph(mu) = -p87(mu)
      enddo
      p87(4)=p87(0)**2-p87(1)**2-p87(2)**2-p87(3)**2
      ph(4) = ph(0)**2 - ph(1)**2 - ph(2)**2 - ph(3)**2 
c 

      do i=1,6
         lfsign(i)=sign(i)
         do mu=0,3
            lpbar(mu,i)=pbar(mu,i)
         enddo
      enddo
c 
      call sqq_qqqqh_nc(lpbar,lfsign,iflav(1,k),ph,res(0,k))
c
c    do not include higgs decay:
      fac = 1d0
c
      do i=0,2                  !loop over color flow
c         do k=1,8               !loop over subprocesses
            res(i,k) = res(i,k)*fac
c         enddo
      enddo
c     return matrix element squareds (so far just neutral currents)
      if(lcol) then
         do i = 0,2
            ansc(i) = res(i,k)
         enddo
      else
         ans = res(0,k)
      endif
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     All q(pbar1)+ q(pbar3) -> q(pbar2)+ q(pbar4) 
c     +q(qbar5)+ qbar(qbar6) + h
c                           |
c                           ---> pbar7 pbar8 
c     Decay products of the higgs are stored in 
c     pbar(0:3,7) and pbar(0:3,8)
c     
c     charged current weak boson fusion processes
c
      subroutine qqHqqQQ_cc_channel(pbar,sign,k,ans)

c k = 1,2 for udsc,ducs
c
      implicit none    
      real * 8 pi,pi2
      parameter (pi=3.141592653589793238462643383279502884197D0,
     1           pi2=pi**2)   
      include '../../include/pwhg_st.h' 
      include 'global_col.inc'  
      double precision  fpi
      parameter ( fpi=4d0*pi)
c
c     Arguments
      real*8 pbar(0:3,6+nv)
      real*8 udsc,ducs
      real*8 udsc_c(0:2),ducs_c(0:2)
      integer sign(6+nv)
      real*8 ans,ansc(0:2)
c
c and which are calculated from output of KOPPLN
c
      double precision  clr, xm2, xmg, b, v, a
      COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),
     1                  V(4,5),A(4,5)
      include 'koppln_ew.inc'
c     
c     Local variables
      real*8 ph(0:4),lpbar(0:3,6),res(0:2,2),res2(0:2,2)
      integer lfsign(6)
      real*8 p87(0:4),p(0:3,6+nv),fac
      integer mu,i,k,kl
      logical lcol
      parameter (lcol = .false.) ! compute color sum
c      parameter (lcol = .true.) ! compute individual color structures

c ====================================================
c
c initialize:
      ans     = 0d0
      ansc(:) = 0d0
c
c     initialize res(k)
      do i = 0,2
         do kl=1,2
            res(i,kl)=0d0
         enddo
      enddo
c     fill local momenta and sign factors
c     diagramatic momenta
      do mu=0,3
         do i=1,6+nv
            p(mu,i)=sign(i)*pbar(mu,i)

         enddo
         p87(mu) = p(mu,8) - p(mu,7)
         ph(mu) = -p87(mu)
      enddo
      p87(4)=p87(0)**2-p87(1)**2-p87(2)**2-p87(3)**2
      ph(4) = ph(0)**2 - ph(1)**2 - ph(2)**2 - ph(3)**2 
c
      do i=1,6
         lfsign(i)=sign(i)
         do mu=0,3
            lpbar(mu,i)=pbar(mu,i)
         enddo
      enddo

      call sqq_qqqqh_cc(lpbar,lfsign,ph,res(0,k))
c
c    do not include higgs decay:
      fac = 1d0

      do i=0,2
c         do k=1,2               !loop over subprocesses
            res(i,k) = res(i,k)*fac
c         enddo
      enddo

      if(lcol) then
         do i=0,2
            ansc(i) = res(i,k)
         enddo
      else
         ans = res(0,k)
      endif
c      
      return
      end
ccc
