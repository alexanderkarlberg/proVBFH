      subroutine setlocalscales(iuborn,imode,rescfac)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_st.h'
      include 'pwhg_flg.h'
      include 'pwhg_math.h'
      include 'PhysPars.h'
      include 'pwhg_flst_2.h'
!     Input
      integer iuborn, imode
!     Output
      double precision rescfac
!     Internal
      double precision zero, one, two, half
      parameter (zero = 0d0, one = 1d0, two = 2d0, half = 0.5d0)
      logical ini
      data ini/.true./
      save ini
      integer i,j,k
      integer ltags(nlegreal), ltags_sub(nlegreal)
      double precision p(0:3,nlegreal)
      double precision as, asQ(2), renfact2, facfact2, mu2, muR2,
     $     KR2,b0,factsc2min,frensc2min
      save factsc2min,frensc2min
!     DIS variables
!     We define one set for each line. Upper: 1, Lower: 2.
      double precision EEC(1:2), Qsq(1:2), EECmax,EECmaxsq,cvirt
      integer npart(2)
      double precision FF(2), expsud
      parameter (cvirt = - 8d0*CF) ! From virtual.f in VBF_H
!     External
      double precision invmsq, dotp, sudakovFF, pwhg_alphas,
     $     sudakov_expansion,powheginput
      external invmsq, dotp, sudakovFF, pwhg_alphas, sudakov_expansion
     $     ,powheginput
!     Cache
      double precision orescfac(2),omufA2(2),omufB2(2)
      integer oltags(3,1:nlegreal)
      save orescfac, oltags,omufA2,omufB2
      logical cache_off
      common/cminlo/cache_off
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      renfact2=st_renfact**2
      facfact2=st_facfact**2
      b0 = (33d0-2d0*st_nlight)/(12*pi)
      
      rescfac = one
      Qsq = zero
      ltags = 0
      p = zero
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!     Uncomment to turn the cache off. Should only be done for
!     debug. 'cahce_off = .false.' will break the program, as the value
!     is stored in a common block and used elsewhere.
!
!      cache_off = .true.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      if(flg_minlo_real) then   ! Enter here only if doing real
         ltags(1:nlegreal)=flst_alrtags(:,flst_cur_alr)
         if(cache_off.or.(any(oltags(2,1:nlegreal).ne.ltags(1:nlegreal)))) then ! Recompute
!     This is the real momenta of the current real configuration
            p(:,1:nlegreal)=kn_cmpreal(:,1:nlegreal)
!     Compute the DIS variables EEC and Q for each line
            call VBF2DIS(nlegreal,p,ltags,Qsq)
!     Set the scales muF and muR and the scale mu entering into alphas
!     Compute form factors using maximum EEC and Q for the respective line
            asQ(1)=pwhg_alphas(Qsq(1)*renfact2,st_lambda5MSB,st_nlight)
            asQ(2)=pwhg_alphas(Qsq(2)*renfact2,st_lambda5MSB,st_nlight)

            st_mufactA2 = max(Qsq(1)*facfact2,1d0) 
            st_mufactB2 = max(Qsq(2)*facfact2,1d0) 
            
            if(sum(ltags).eq.8.or.sum(ltags).eq.28) then
               rescfac = rescfac * (asQ(1)/st_alpha)**2
            elseif(sum(ltags).eq.9) then
               rescfac = rescfac * asQ(1)/st_alpha * asQ(2)/st_alpha
            elseif(sum(ltags).eq.10.or.sum(ltags).eq.30) then
               rescfac = rescfac * (asQ(2)/st_alpha)**2
            else
               print*, 'Complain...'
            endif

!     Save variables for cache
            oltags(2,1:nlegreal) = ltags(1:nlegreal)
            orescfac(2) = rescfac
            omuFA2(2) = st_mufactA2
            omuFB2(2) = st_mufactB2
         else
            rescfac = orescfac(2)
            st_mufactA2 = max(omuFA2(2),1d0)
            st_mufactB2 = max(omuFB2(2),1d0)
         endif
      else
         ltags(1:nlegborn)=flst_borntags(:,iuborn)
         ltags_sub(1:nlegreal)=flst_alrtags(:,flst_cur_alr)
         if(cache_off.or.(any(oltags(1,1:nlegborn).ne.ltags(1:nlegborn)))
     $        .or.(any(oltags(3,1:nlegreal)
     $        .ne.ltags_sub(1:nlegreal)).and.imode.eq.3)) 
     $        then              ! Recompute
            
            p(:,1:nlegborn)=kn_cmpborn(:,1:nlegborn)
            
            call VBF2DIS(nlegborn,p,ltags,Qsq)
            
!     Alphas for the two lines
            asQ(1)=pwhg_alphas(Qsq(1)*renfact2,st_lambda5MSB ,st_nlight)
            asQ(2)=pwhg_alphas(Qsq(2)*renfact2,st_lambda5MSB ,st_nlight)
            
            st_mufactA2 = max(Qsq(1)*facfact2,1d0) 
            st_mufactB2 = max(Qsq(2)*facfact2,1d0) 
            
!     muR is computed according to where the Born emission is
            if(ltags(nlegborn).eq.1) then
               muR2 = max(Qsq(1)*renfact2,1d0) 
               as=asQ(1)
            elseif(ltags(nlegborn).eq.2) then
               muR2 = max(Qsq(2)*renfact2,1d0) 
               as=asQ(2)
            endif
            
            if(imode.eq.2) then ! virtual
               rescfac = rescfac * (as/st_alpha)**2 ! Virtual
            elseif(imode.eq.3) then ! subtraction
               if(sum(ltags_sub(:)).eq.8
     $              .or.sum(ltags_sub(:)).eq.28)
     $              then
                  rescfac = rescfac * (asQ(1)/st_alpha)**2
               elseif(sum(ltags_sub(:)).eq.9) then
                  rescfac = rescfac * asQ(1)/st_alpha * asQ(2)/st_alpha
               elseif(sum(ltags_sub(:)).eq.10
     $                 .or.sum(ltags_sub(:)).eq.30)
     $                 then
                  rescfac = rescfac * (asQ(2)/st_alpha)**2
               else
                  print*, 'Complain...'
               endif
            else                ! born
               rescfac = rescfac * as/st_alpha ! Born emission 
               if(.not.flg_bornonly) then
                  if(ltags(nlegborn).eq.1) then
                     rescfac = rescfac * (1d0+as*b0*log(muR2 /st_muren2)
     $                    -(asQ(1)-asQ(2))*cvirt/(2d0*pi))
                  elseif(ltags(nlegborn).eq.2) then
                     rescfac = rescfac * (1d0+as*b0*log(muR2 /st_muren2)
     $                    -(asQ(2)-asQ(1))*cvirt/(2d0*pi))
                  endif
               endif
            endif
         !     Save variables for cache
            oltags(1,1:nlegborn) = ltags(1:nlegborn)
            oltags(3,1:nlegreal) = ltags_sub(1:nlegreal)
            orescfac(1) = rescfac
            omuFA2(1) = st_mufactA2
            omuFB2(1) = st_mufactB2 
         else
            rescfac = orescfac(1)
            st_mufactA2 = max(omuFA2(1),1d0)
            st_mufactB2 = max(omuFB2(1),1d0)
         endif
      endif
      
      end

!     This routine computes the DIS variables EEC and Q from a set of
!     VBF like kinematics.
      subroutine VBF2DIS(nlegs,p,ltags,Qsq)
      implicit none
      double precision zero, one, two, half
      parameter (zero = 0d0, one = 1d0, two = 2d0, half = 0.5d0)
!     Input
      integer nlegs, ltags(1:nlegs)
      double precision p(0:3,nlegs)
!     Output
      double precision EEC(2), Qsq(2)
      integer npart(2)

      integer i,j,ileg
      double precision Pcurr(0:3,1:2), Premn(0:3,1:2), Q(0:3,2),
     $     Qval(1:2)
!     First index ith particle, second upper/lower
      double precision alpha(1:4,2), beta(1:4,2), Qb(0:3,1:2)
      double precision alphaprime(1:2)
      double precision pbreit(0:3,4,1:2),pin(0:3,1:4,1:2)
!     External
      double precision invmsq, dotp,heaviside,sinthetaij
      external invmsq, dotp,heaviside,sinthetaij
!     Counter
      integer total, large
      data total,large /0,0/
      save total, large

!     Initialise
      Q = zero
      npart = 0

      Q(:,1:2) = p(:,1:2)       ! Incoming partons
      npart = 1

      do i=4,nlegs ! Skip incoming and Higgs
         if(ltags(i).eq.1.or.ltags(i).eq.11) then ! Upper line
            Q(:,1) = Q(:,1) - p(:,i)
            npart(1) = npart(1) + 1
         elseif(ltags(i).eq.2.or.ltags(i).eq.12) then ! Lower line
            Q(:,2) = Q(:,2) - p(:,i)
            npart(2) = npart(2) + 1
         else                   ! No tag
            print*, 'ltags', ltags
            stop
         endif
      enddo

      if(any(npart.lt.2)) stop 'Not enough particles...'
      
      Qsq(1) = - invmsq(Q(:,1))
      Qsq(2) = - invmsq(Q(:,2))
      
      if(any(Qsq(:).lt.zero)) print*, 'Qsq', Qsq

!     We need to avoid extreme configurations. These are anyways
!     suppressed by the matrix elements.
      Qsq(1) = max(Qsq(1),1d0)
      Qsq(2) = max(Qsq(2),1d0)

      end

      double precision function pt(p)
      double precision p(0:3)

      pt = sqrt(p(1)**2 + p(2)**2)
      end

      double precision function invm(p)
      double precision p(0:3)

      invm = sqrt(p(0)**2 - p(1)**2 - p(2)**2 - p(3)**2)
      end

      double precision function transvm(p)
      double precision p(0:3)

      transvm = sqrt(p(0)**2d0-p(3)**2d0)
      end

      double precision function transvmsq(p)
      double precision p(0:3)
      double precision invmsq,pt

      transvmsq = invmsq(p) - pt(p)**2 
      end

      double precision function ptsq(p)
      double precision p(0:3)

      ptsq = p(1)**2 + p(2)**2
      end

      double precision function invmsq(p)
      double precision p(0:3)

      invmsq = p(0)**2 - p(1)**2 - p(2)**2 - p(3)**2
      end

!     Powheg needs this routine to compile....
      subroutine sudakov_exponent(q20,q2h,theExponent,isquark,theAccuracy)
      implicit none 
      logical isQuark 
      integer theAccuracy
      real * 8 q2h,q20,theExponent

      stop 'Dummy sudakov_exponent called!'
      end

      double precision function phidist(phi1,phi2)
      implicit none
      include 'pwhg_math.h'
      double precision phi1,phi2

      phidist = abs(phi2-phi1)

      if(phidist.gt.pi) then
         phidist = abs(2d0 * pi - phidist)
      endif
      
      end
