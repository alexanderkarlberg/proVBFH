****************subroutine qqh4q ****************************************
!     
!     Michael Rauch, <michael.rauch@kit.edu>
!     Initial version:  2017, February
!     Last modified: 2017, February
!     
!     qqh4q calculates the matrix elements**2 for
!     
!     /--<-- qbar5
!     8\-->-- q6
!     g 8
!     8
!     q1 -->---->------>-- q2
!     S
!     W/Z S
!     S
!     S- - -  h (p7+p8)
!     S
!     W/Z S
!     S
!     q3 -->---->------>-- q4
!     
!     and crossing related processes.
!     Interference terms between the two quark lines originating
!     from the proton are neglected, interference terms from quark lines
!     connected by a gluon are included.

      subroutine qqh4q(pbar,fsign,m2snc,m2scc,lsymmcontrib,cc_type,nc_type)

c     use globalvars, only: ldoblha

      implicit none

      
      real * 8 pi,pi2
      parameter (pi=3.141592653589793238462643383279502884197D0,
     1     pi2=pi**2)   
      include 'pwhg_st.h'
      
      include 'global.inc'
      include 'koppln_ew.inc'
      include 'scales.inc'
      include 'nlegborn.h'
      include 'pwhg_flst.h' 
      include 'pwhg_flst_2.h' 
      include 'tags.h' 
      include 'leptens.h' 
      integer tags(nlegreal)
      real * 8 tfac(2)
      integer tag_factor
      common/cdoubletag/tag_factor
! arguments
      double precision, intent(in) :: pbar(0:3,6+max_v)
      integer, intent(in) :: fsign(6+max_v)
! indices: flavor%2 of quarks 1, 3 and 5 (1-flavor%2 of quark 2/4 if 1/3 is an anti-quark),
!          color structure (1/2: radiation from upper/lower line,
!                           3/4: flipped diagram for different quark generations on upper/lower line)
!          interference
      double precision, intent(out) :: m2snc(0:1,0:1,0:1,1:2,0:1)
      double precision, intent(out) :: m2scc(0:1,0:1,0:1,1:4,0:1)
      logical, intent(in) :: lsymmcontrib, nc_type, cc_type

! options and checks
      logical lonlyinterference
      parameter (lonlyinterference = .false.)
!double complex qterm(1:4)

! local
      double precision fpials(2:3), fpi
      parameter (fpi=4d0*pi)
!double complex zm2i(2:3)
      integer i,isigg,isig1,isig3,isig5,flav1,flav3,flav5
      double precision p(0:3,6+max_v)
      double precision p21(0:4)
      double precision p43(0:4)
      double precision p65(0:4)
      double precision p87(0:4)
      double precision p109(0:4)
      double precision ph(0:4)
      double precision p2165(0:4)
      double precision p4365(0:4)
      double precision fac
      double precision facZlower,facZupper,facWlower,facWupper
      double complex psi(2,-1:1,6)
! last index are combinations (12,16,34,36,52,54,56)
      double complex gcurr(1:6,-1:1,7)
      double precision gmom(0:4,-1:1,7)
! last two index pairs: quark psi, gluon-forming pair combination
      double complex psig(2,-1:1,6,-1:1,7)
      double precision pmomg(0:4,6,7)
! last three indices are "upper/lower line", "normal or crossed
! fermions" and "gluon pair(1:56/52/54,2:12/34/16/36);
! vcurr(:,:,-1,:,1,3) is the other line
      double complex vcurr(1:6,-1:1,-1:1,2,2,3)
! indices: 3*isig,3*flav(3:u,4:d,5:g), upper/lower
      double complex mat(-1:1,-1:1,-1:1,3:4,3:4,3:4,2)
      double complex matflip(-1:1,-1:1,-1:1,3:4,3:4,3:4,2)
      double precision m2stmp(0:1)
      double precision colsq, colint
      integer issymm
      double precision multfact

! external
      double precision  clr, xm2, xmg, b, v, a
      COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),
     &     V(4,5),A(4,5)
! anomalous couplings (resulting from electroweak corrections)
      double complex ahvv(3,4,4), ahvvL(3,4,4)
      common/tensorhvv/ ahvv, ahvvL
! CPS scheme
      double precision qgammaq
      common /VBFNLO_HIGGS_CPS/ qgammaq

      double precision mass2
      double complex dotcc, dotrc, contract_Tjj
      external mass2, dotcc, dotrc, contract_Tjj

!     initialization
      tfac = 1d0
      if (realequiv_tag > 0) then 
         tags = flst_realtags(:,realequiv_tag)
      elseif (alr_tag > 0) then 
         tags = flst_alrtags(:,alr_tag)
      endif
      
C     kill amplitudes where gluons are connected to wrong fermion lines 
C     (based on comments above about colour structure of the amplitudes) 
      if (any(tags(1:2) == iqpairtag*tag_factor) .or. any(tags(5:8) == iqpairtag*tag_factor)) then 
C     running without fulltags, so do nothing here 

      elseif (sum(tags(1:nlegreal)) == 8+2*iqpairtag*tag_factor) then ! case 1 2 -> 1 2 1 1 or permutations of it 
         tfac(2) = 0d0
      elseif (sum(tags(1:nlegreal)) == 10+2*iqpairtag*tag_factor) then ! case 1 2 -> 1 2 2 2 or permutations of it 
         tfac(1) = 0d0
      else
         write(*,*) 'tags', tags 
         write(*,*) 'iqpairtag, tag_factor',iqpairtag, tag_factor
         stop 'compreal_hqqqq_new: invalid tags'
      endif

c     eliminate any higgsstrahlung graphs
      if (fsign(5).eq.fsign(6)) then !either q5 or q6 in i.s.
         if (fsign(1).eq.-fsign(2) .and. 
     1        fsign(3).eq.fsign(4) ) then
c     cf(2,2) = 0.0d0         !  initial gluon attached to 1-2 line only
            tfac(2) = 0d0
         elseif (fsign(1).eq.fsign(2) .and. 
     1           fsign(3).eq.-fsign(4) ) then
c     cf(1,1)= 0.0d0          !  initial gluon attached to 3-4 line only
            tfac(1) = 0d0
         else 
            stop
         endif
      endif

      
      m2snc = 0
      m2scc = 0

!     fix strong coupling gs**2 for the two quarks:
      fpials(2) = fpi*als(1,1)
      fpials(3) = fpi*als(2,1)
!     fpials(2) = fpi*alfas

!zm2i(2) = dcmplx(xm2(2),-xmg(2))
!zm2i(3) = dcmplx(xm2(3),-xmg(3))

!     color factors
c     if (ldoblha) then
c     colsq = blha_CF*blha_CA**2/2d0
c     colint = - blha_CF*blha_CA**2*(blha_CA/2d0-blha_CF)
c     else
      colsq = 6
      colint = -2
c     endif

!     fill local momenta and sign factors
!     diagramatic momenta
      do i=1,6+max(n_v,2)
         p(0:3,i)=fsign(i)*pbar(0:3,i)
      enddo

      p21(0:3) = p(0:3,2) - p(0:3,1)
      p43(0:3) = p(0:3,4) - p(0:3,3)
      p65(0:3) = p(0:3,6) - p(0:3,5)
      p87(0:3) = p(0:3,8) - p(0:3,7)
      if (n_v.eq.4) then
         stop 'nv=4 is not covered in this code version'
c     p109(0:3) = p(0:3,10) - p(0:3,9)
c     ph(0:3) = p87(0:3) + p109(0:3)
      elseif (n_v.eq.2 .or. n_v.eq.1) then
         ph(0:3) = p87(0:3)
c     BJ: maybe need to change sign of ph
         
      endif
      p2165(0:3) = p21(0:3)+p65(0:3)
      p4365(0:3) = p43(0:3)+p65(0:3)

      p21(4) = mass2(p21)
      p43(4) = mass2(p43)
      p65(4) = mass2(p65)
      p87(4) = mass2(p87)
      p109(4) = mass2(p109)
      ph(4) = mass2(ph)
      p2165(4) = mass2(p2165)
      p4365(4) = mass2(p4365)

      fac=1d0

!     Prepare momentum of V currents with gluon pair attached
      vcurr = 0
      vcurr(5,:,:,1,:,:) = dcmplx(
     &     p(0,1)-p(0,2)+p(0,5)-p(0,6),
     &     p(3,1)-p(3,2)+p(3,5)-p(3,6))
      vcurr(6,:,:,1,:,:) = dcmplx(
     &     p(1,1)-p(1,2)+p(1,5)-p(1,6),
     &     p(2,1)-p(2,2)+p(2,5)-p(2,6))
      vcurr(5,:,:,2,:,:) = dcmplx(
     &     p(0,3)-p(0,4)+p(0,5)-p(0,6),
     &     p(3,3)-p(3,4)+p(3,5)-p(3,6))
      vcurr(6,:,:,2,:,:) = dcmplx(
     &     p(1,3)-p(1,4)+p(1,5)-p(1,6),
     &     p(2,3)-p(2,4)+p(2,5)-p(2,6))

! VVH coupling, alphas and V propagators
      facZupper = fpials(2)*fpials(3)
      facZlower = fpials(2)*fpials(3)
      facWupper = fpials(2)*fpials(3)
      facWlower = fpials(2)*fpials(3)

!-----------------------------------------------------

!     Get fermion currents
      call psi0m(6,pbar(0,1),fsign(1),psi)

!     Construct gluon currents splitting into q-qbar
!     ! 1-2
      call curr6(+1,psi(1,-1,2),p(0,2),psi(1,-1,1),p(0,1),gcurr(1,-1,1))
!     ! 1-6
      call curr6(+1,psi(1,-1,6),p(0,6),psi(1,-1,1),p(0,1),gcurr(1,-1,2))
!     ! 3-4
      call curr6(+1,psi(1,-1,4),p(0,4),psi(1,-1,3),p(0,3),gcurr(1,-1,3))
!     ! 3-6
      call curr6(+1,psi(1,-1,6),p(0,6),psi(1,-1,3),p(0,3),gcurr(1,-1,4))
!     ! 5-2
      call curr6(+1,psi(1,-1,2),p(0,2),psi(1,-1,5),p(0,5),gcurr(1,-1,5))
!     ! 5-4
      call curr6(+1,psi(1,-1,4),p(0,4),psi(1,-1,5),p(0,5),gcurr(1,-1,6))
!     ! 5-6
      call curr6(+1,psi(1,-1,6),p(0,6),psi(1,-1,5),p(0,5),gcurr(1,-1,7))

!     ! Save non-gluon-pair line entries
! 12
      vcurr(:,:,-1,2,1,3) = gcurr(:,:,1)
! 34
      vcurr(:,:,-1,1,1,3) = gcurr(:,:,3)

!     ! Add gluon propagator
      do i=1,7
         call propagate(+1,gcurr(1,-1,i),gmom(0,-1,i))
      enddo

!     Attach to the fermion line
      call qqh4q_gtobraket(1,7,psi,p,gcurr,gmom,psig,pmomg) ! 1 - 56
      call qqh4q_gtobraket(2,7,psi,p,gcurr,gmom,psig,pmomg) ! 2 - 56
      call qqh4q_gtobraket(3,7,psi,p,gcurr,gmom,psig,pmomg) ! 3 - 56
      call qqh4q_gtobraket(4,7,psi,p,gcurr,gmom,psig,pmomg) ! 4 - 56
      call qqh4q_gtobraket(5,1,psi,p,gcurr,gmom,psig,pmomg) ! 5 - 12
      call qqh4q_gtobraket(6,1,psi,p,gcurr,gmom,psig,pmomg) ! 6 - 12
      call qqh4q_gtobraket(5,3,psi,p,gcurr,gmom,psig,pmomg) ! 5 - 34
      call qqh4q_gtobraket(6,3,psi,p,gcurr,gmom,psig,pmomg) ! 6 - 34
! crossed diagrams
      call qqh4q_gtobraket(1,5,psi,p,gcurr,gmom,psig,pmomg) ! 1 - 52
      call qqh4q_gtobraket(6,5,psi,p,gcurr,gmom,psig,pmomg) ! 6 - 52
      call qqh4q_gtobraket(3,6,psi,p,gcurr,gmom,psig,pmomg) ! 3 - 54
      call qqh4q_gtobraket(6,6,psi,p,gcurr,gmom,psig,pmomg) ! 6 - 54
      call qqh4q_gtobraket(5,2,psi,p,gcurr,gmom,psig,pmomg) ! 5 - 16
      call qqh4q_gtobraket(2,2,psi,p,gcurr,gmom,psig,pmomg) ! 2 - 16
      call qqh4q_gtobraket(5,4,psi,p,gcurr,gmom,psig,pmomg) ! 5 - 36
      call qqh4q_gtobraket(4,4,psi,p,gcurr,gmom,psig,pmomg) ! 4 - 36

!     Connect into current
      do isigg=-1,1,2
!! 2(56) - 1
         call curr6add(+1,psig(1,-1,2,isigg,7),pmomg(0,2,7),
     &        psi(1,-1,1),p(0,1),vcurr(1,-1,isigg,1,1,1))
!! 2 - 1(56)
         call curr6add(+1,psi(1,-1,2),p(0,2),
     &        psig(1,-1,1,isigg,7),pmomg(0,1,7),vcurr(1,-1,isigg,1,1,1))
!! 4(56) - 3
         call curr6add(+1,psig(1,-1,4,isigg,7),pmomg(0,4,7),
     &        psi(1,-1,3),p(0,3),vcurr(1,-1,isigg,2,1,1))
!! 4 - 3(56)
         call curr6add(+1,psi(1,-1,4),p(0,4),
     &        psig(1,-1,3,isigg,7),pmomg(0,3,7),vcurr(1,-1,isigg,2,1,1))
!! 6(12) - 5
         call curr6add(+1,psig(1,-1,6,isigg,1),pmomg(0,6,1),
     &        psi(1,-1,5),p(0,5),vcurr(1,-1,isigg,1,1,2))
!! 6 - 5(12)
         call curr6add(+1,psi(1,-1,6),p(0,6),
     &        psig(1,-1,5,isigg,1),pmomg(0,5,1),vcurr(1,-1,isigg,1,1,2))
!! 6(34) - 5
         call curr6add(+1,psig(1,-1,6,isigg,3),pmomg(0,6,3),
     &        psi(1,-1,5),p(0,5),vcurr(1,-1,isigg,2,1,2))
!! 6 - 5(34)
         call curr6add(+1,psi(1,-1,6),p(0,6),
     &        psig(1,-1,5,isigg,3),pmomg(0,5,3),vcurr(1,-1,isigg,2,1,2))
! crossed diagrams
!! 6(52) - 1
         call curr6add(+1,psig(1,-1,6,isigg,5),pmomg(0,6,5),
     &        psi(1,-1,1),p(0,1),vcurr(1,-1,isigg,1,2,1))
!! 6 - 1(52)
         call curr6add(+1,psi(1,-1,6),p(0,6),
     &        psig(1,-1,1,isigg,5),pmomg(0,1,5),vcurr(1,-1,isigg,1,2,1))
!! 6(54) - 3
         call curr6add(+1,psig(1,-1,6,isigg,6),pmomg(0,6,6),
     &        psi(1,-1,3),p(0,3),vcurr(1,-1,isigg,2,2,1))
!! 6 - 3(54)
         call curr6add(+1,psi(1,-1,6),p(0,6),
     &        psig(1,-1,3,isigg,6),pmomg(0,3,6),vcurr(1,-1,isigg,2,2,1))
!! 2(16) - 5
         call curr6add(+1,psig(1,-1,2,isigg,2),pmomg(0,2,2),
     &        psi(1,-1,5),p(0,5),vcurr(1,-1,isigg,1,2,2))
!! 2 - 5(16)
         call curr6add(+1,psi(1,-1,2),p(0,2),
     &        psig(1,-1,5,isigg,2),pmomg(0,5,2),vcurr(1,-1,isigg,1,2,2))
!! 4(36) - 5
         call curr6add(+1,psig(1,-1,4,isigg,4),pmomg(0,4,4),
     &        psi(1,-1,5),p(0,5),vcurr(1,-1,isigg,2,2,2))
!! 4 - 5(36)
         call curr6add(+1,psi(1,-1,4),p(0,4),
     &        psig(1,-1,5,isigg,4),pmomg(0,5,4),vcurr(1,-1,isigg,2,2,2))
      enddo
      if(nc_type) then
!     Connect with lower line -- NC
         do flav1=3,4
            do flav3=3,4
               do flav5=3,4
                  do isig1=-1,1,2
                     do isig3=-1,1,2
                        do isig5=-1,1,2
!     radiation from upper, normal
                           if(tfac(1).ne.0d0) then
                              mat(isig1,isig3,isig5,flav1,flav3,flav5,1) =
!     ! (12(56))-(34)
     &                        + contract_Tjj(vvhh_re(0,0,1,1),
     &                            vcurr(1,isig1,isig5,1,1,1),vcurr(1,isig3,-1,1,1,3))
     &                             *clr(flav1,2,isig1)*clr(flav3,2,isig3)
!     ! (56(12))-(34)
     &                        + contract_Tjj(vvhh_re(0,0,1,1),
     &                            vcurr(1,isig5,isig1,1,1,2),vcurr(1,isig3,-1,1,1,3))
     &                             *clr(flav5,2,isig5)*clr(flav3,2,isig3)
!     radiation from lower, normal
                           endif
                           if(tfac(2).ne.0d0) then
                              mat(isig1,isig3,isig5,flav1,flav3,flav5,2) =
!     ! (12)-(34(56))
     &                        + contract_Tjj(vvhh_re(0,0,2,1),
     &                            vcurr(1,isig3,isig5,2,1,1),vcurr(1,isig1,-1,2,1,3))
     &                             *clr(flav3,2,isig3)*clr(flav1,2,isig1)
!     ! (12)-(56(34))
     &                        + contract_Tjj(vvhh_re(0,0,2,1),
     &                            vcurr(1,isig5,isig3,2,1,2),vcurr(1,isig1,-1,2,1,3))
     &                             *clr(flav5,2,isig5)*clr(flav1,2,isig1)
                           endif
!     radiation from upper, flipped - minus signs due to interchange of fermions
                           if(tfac(1).ne.0d0) then
                              matflip(isig1,isig3,isig5,flav1,flav3,flav5,1) =
!     ! (16(52))-(34)
     &                        - contract_Tjj(vvhh_re(0,0,1,1),
     &                            vcurr(1,isig1,isig5,1,2,1),vcurr(1,isig3,-1,1,1,3))
     &                             *clr(flav1,2,isig1)*clr(flav3,2,isig3)
!     ! (52(16))-(34)
     &                        - contract_Tjj(vvhh_re(0,0,1,1),
     &                            vcurr(1,isig5,isig1,1,2,2),vcurr(1,isig3,-1,1,1,3))
     &                             *clr(flav5,2,isig5)*clr(flav3,2,isig3)
!     radiation from lower, flipped
                           endif
                           if(tfac(2).ne.0d0) then
                              matflip(isig1,isig3,isig5,flav1,flav3,flav5,2) =
!     ! (12)-(36(54))
     &                        - contract_Tjj(vvhh_re(0,0,2,1),
     &                            vcurr(1,isig3,isig5,2,2,1),vcurr(1,isig1,-1,2,1,3))
     &                             *clr(flav3,2,isig3)*clr(flav1,2,isig1)
!     ! (12)-(54(36))
     &                        - contract_Tjj(vvhh_re(0,0,2,1),
     &                            vcurr(1,isig5,isig3,2,2,2),vcurr(1,isig1,-1,2,1,3))
     &                             *clr(flav5,2,isig5)*clr(flav1,2,isig1)
                           endif
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
         
         
!     Sum everything up -- NC
         do flav1=0,1
            do flav3=0,1
               do flav5=0,1
                  do isig1=-1,1,2
                     do isig3=-1,1,2
                        do isig5=-1,1,2
!     standard diagram - upper line
                           if(tfac(1).ne.0d0) then
                              m2stmp = 0
                              if (.not.lonlyinterference) then
                                 m2stmp =
     &                                colsq* ! colour factor
     &                                abs(mat(isig1,isig3,isig5,3+flav1,3+flav3,3+flav5,1))**2
                              endif
                              if (lsymmcontrib) then
!     flipped
                                 if (.not.lonlyinterference) then
                                    m2stmp(1) = m2stmp(1)
     &                                   + colsq* ! colour factor
     &                                   abs(matflip(isig1,isig3,isig5,3+flav1,3+flav3,3+flav5,1))**2
                                 endif
                                 if (isig1.eq.isig5) then
!     interference
                                    m2stmp(1) = m2stmp(1)
     &                                   + 2*colint*dreal(
     &                                   mat(isig1,isig3,isig5,3+flav1,3+flav3,3+flav5,1)*
     &                                   dconjg(matflip(isig1,isig3,isig5,3+flav1,3+flav3,3+flav5,1))
     &                                   )
                                 endif
                              else
                                 if (isig1.eq.isig5) then
                                    m2stmp(1) = m2stmp(1)
     &                                   + colint*dreal(
     &                                   mat(isig1,isig3,isig5,3+flav1,3+flav3,3+flav5,1)*
     &                                   dconjg(matflip(isig1,isig3,isig5,3+flav1,3+flav3,3+flav5,1))
     &                                   )
                                 endif
                              endif
                              m2snc(flav1,flav3,flav5,1,:) = m2snc(flav1,flav3,flav5,1,:)
     &                             + m2stmp(:)*fac*facZupper
                           endif
!     standard diagram - lower line
                           if(tfac(2).ne.0d0) then
                              m2stmp = 0
                              if (.not.lonlyinterference) then
                                 m2stmp =
     &                                colsq* ! colour factor
     &                                abs(mat(isig1,isig3,isig5,3+flav1,3+flav3,3+flav5,2))**2
                              endif
                              if (lsymmcontrib) then
!     flipped
                                 if (.not.lonlyinterference) then
                                    m2stmp(1) = m2stmp(1)
     &                                   + colsq* ! colour factor
     &                                   abs(matflip(isig1,isig3,isig5,3+flav1,3+flav3,3+flav5,2))**2
                                 endif
                                 if (isig3.eq.isig5) then
!     interference
                                    m2stmp(1) = m2stmp(1)
     &                                   + 2*colint*dreal(
     &                                   mat(isig1,isig3,isig5,3+flav1,3+flav3,3+flav5,1)*
     &                                   dconjg(matflip(isig1,isig3,isig5,3+flav1,3+flav3,3+flav5,2))
     &                                   )
                                 endif
                              else
                                 if (isig3.eq.isig5) then
                                    m2stmp(1) = m2stmp(1)
     &                                   + colint*dreal(
     &                                   mat(isig1,isig3,isig5,3+flav1,3+flav3,3+flav5,1)*
     &                                   dconjg(matflip(isig1,isig3,isig5,3+flav1,3+flav3,3+flav5,2))
     &                                   )
                                 endif
                              endif
                              m2snc(flav1,flav3,flav5,2,:) = m2snc(flav1,flav3,flav5,2,:)
     &                             + m2stmp(:)*fac*facZlower
                           endif
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      endif
!     Connect with lower line -- CC
!     ! possible contributions for the upper line are
!     ! flav135: mat           |  matflip
!     ! 333/444: ---           | ---
!     ! 334/443: (56(12))-(34) | fsign(1)==+1 ? (52(16))-(34)*delta_15 : (16(52))-(34)*delta_15 
!     ! 343/434: (12(56))-(34) | (52(16))-(34)*(delta_15 || fsign(1)==+1)
!     ! 344/433: (12(56))-(34) | (16(52))-(34)*(delta_15 || fsign(1)==-1)
      mat = 0
      matflip = 0
      if(cc_type) then
         do flav1=3,4
            do flav3=3,4
               do flav5=3,4
                  do isig1=-1,1,2
                     do isig3=-1,1,2
                        do isig5=-1,1,2
                           if (flav1.ne.flav3 .and. isig1.eq.-1 .and.
     $                          isig3.eq.-1) then ! W only couples to left-handed fermions
!     radiation from upper, normal, W on 12
!     ! (12(56))-(34)
                              if(tfac(1).ne.0d0) then
                                 mat(isig1,isig3,isig5,flav1,flav3,flav5,1) =
     &                           + contract_Tjj(vvhh_re(0,0,1,2),
     &                               vcurr(1,isig1,isig5,1,1,1),vcurr(1,isig3,-1,1,1,3))
     &                                *clr(flav1,3,isig1)*clr(flav3,3,isig3)
!     radiation from lower, normal, W on 34
!     ! (12)-(34(56))
                              endif
                              if(tfac(2).ne.0d0) then
                                 mat(isig1,isig3,isig5,flav1,flav3,flav5,2) =
     &                           + contract_Tjj(vvhh_re(0,0,2,2),
     &                               vcurr(1,isig3,isig5,2,1,1),vcurr(1,isig1,-1,2,1,3))
     &                                *clr(flav3,3,isig3)*clr(flav1,3,isig1)
                              endif
                           elseif (flav1.eq.flav3 .and. flav1.ne.flav5
     $                             .and. isig5.eq.-1) then
                              if (isig3.eq.-1.and.tfac(1).ne.0d0) then
!     radiation from upper, normal, W on 56
!     ! (56(12))-(34)
                                 mat(isig1,isig3,isig5,flav1,flav3,flav5,1) =
     &                             + contract_Tjj(vvhh_re(0,0,1,2),
     &                                 vcurr(1,isig5,isig1,1,1,2),vcurr(1,isig3,-1,1,1,3))
     &                                *clr(flav5,3,isig5)*clr(flav3,3,isig3)
                              endif
                              if (isig1.eq.-1.and.tfac(2).ne.0d0) then
!     radiation from lower, normal, W on 56
!     ! (12)-(56(34))
                                 mat(isig1,isig3,isig5,flav1,flav3,flav5,2) =
     &                             + contract_Tjj(vvhh_re(0,0,2,2),
     &                                 vcurr(1,isig5,isig3,2,1,2),vcurr(1,isig1,-1,2,1,3))
     &                                *clr(flav5,3,isig5)*clr(flav1,3,isig1)
                              endif
                           endif
!     flipped lines -- minus signs due to interchange of fermions
                           if ((flav3.ne.flav5 .and. isig3.eq.-1 .and.
     $                          isig5.eq.-1 .and.( flav1.eq.flav5 .or.
     $                          fsign(1).eq.+1 )).and.tfac(1).ne.0d0 )
     $                          then
!     radiation from upper, flipped, W on 52
!     ! (52(16))-(34)
                              matflip(isig1,isig3,isig5,flav1,flav3,flav5,1) =
     &                           - contract_Tjj(vvhh_re(0,0,1,2),
     &                               vcurr(1,isig5,isig1,1,2,2),vcurr(1,isig3,-1,1,1,3))
     &                             *clr(flav5,3,isig5)*clr(flav3,3,isig3)
                           endif
                           if ((flav1.ne.flav5 .and. isig1.eq.-1 .and.
     $                          isig5.eq.-1 .and.( flav3.eq.flav5 .or.
     $                          fsign(3).eq.+1 )).and.tfac(2).ne.0d0 )
     $                          then
!     radiation from lower, flipped, W on 54
!     ! (12)-(54(36))
                              matflip(isig1,isig3,isig5,flav1,flav3,flav5,2) =
     &                           - contract_Tjj(vvhh_re(0,0,2,2),
     &                               vcurr(1,isig5,isig3,2,2,2),vcurr(1,isig1,-1,2,1,3))
     &                             *clr(flav5,3,isig5)*clr(flav1,3,isig1)
                           endif
                           if (isig1.eq.-1 .and. isig3.eq.-1) then
                              if ( (( flav1.ne.flav3 .and.
     $                             flav1.ne.flav5) .or.( flav1.eq.flav3
     $                             .and. flav1.ne.flav5 .and.
     $                             fsign(1).eq.-1 )).and.tfac(1).ne.0d0
     $                             ) then
!     radiation from upper, flipped, W on 16
!     ! (16(52))-(34)
                                 matflip(isig1,isig3,isig5,flav1,flav3,flav5,1) =
     &                             - contract_Tjj(vvhh_re(0,0,1,2),
     &                                 vcurr(1,isig1,isig5,1,2,1),vcurr(1,isig3,-1,1,1,3))
     &                                *clr(flav1,3,isig1)*clr(flav3,3,isig3)
                              endif
                              if ((flav1.ne.flav3 .and. flav1.eq.flav5
     $                             .or.flav1.eq.flav3 .and.
     $                             flav3.ne.flav5 .and. fsign(3).eq.
     $                             -1).and.tfac(2).ne.0d0 ) then
!     radiation from lower, flipped, W on 36
!     ! (12)-(36(54))
                                 matflip(isig1,isig3,isig5,flav1,flav3,flav5,2) =
     &                             - contract_Tjj(vvhh_re(0,0,2,2),
     &                                 vcurr(1,isig3,isig5,2,2,1),vcurr(1,isig1,-1,2,1,3))
     &                                *clr(flav3,3,isig3)*clr(flav1,3,isig1)
                              endif
                           endif
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
         
!     Sum everything up -- CC
         do flav1=0,1
            do flav3=0,1
               do flav5=0,1
                  do isig1=-1,1,2
                     do isig3=-1,1,2
                        do isig5=-1,1,2
!     upper line
                           if(tfac(1).ne.0d0) then
                              m2stmp = 0
                              if (.not.lonlyinterference) then
!     normal
                                 m2stmp =
     &                                colsq* ! colour factor
     &                                abs(mat(isig1,isig3,isig5,3+flav1,3+flav3,3+flav5,1))**2
!     flipped
                                 if (lsymmcontrib) then
                                    m2stmp(1) = m2stmp(1)
     &                                   + colsq* ! colour factor
     &                                   abs(matflip(isig1,isig3,isig5,3+flav1,3+flav3,3+flav5,1))**2
                                 endif
                              endif
                              if (lsymmcontrib) then
                                 if (isig1.eq.isig5) then
!     interference
                                    m2stmp(1) = m2stmp(1)
     &                                   + 2*colint*dreal(
     &                                   mat(isig1,isig3,isig5,3+flav1,3+flav3,3+flav5,1)*
     &                                   dconjg(matflip(isig1,isig3,isig5,3+flav1,3+flav3,3+flav5,1))
     &                                   )
                                 endif
                              else
                                 if (isig1.eq.isig5) then
                                    m2stmp(1) = m2stmp(1)
     &                                   + colint*dreal(
     &                                   mat(isig1,isig3,isig5,3+flav1,3+flav3,3+flav5,1)*
     &                                   dconjg(matflip(isig1,isig3,isig5,3+flav1,3+flav3,3+flav5,1))
     &                                   )
                                 endif
                              endif
                              m2scc(flav1,flav3,flav5,1,:) = m2scc(flav1,flav3,flav5,1,:)
     &                             + m2stmp(:)*fac*facWupper
                           endif
                           m2stmp = 0
                           if (.not. lonlyinterference) then
                              m2stmp(0) = m2stmp(0)
     &                             + colsq* ! colour factor
     &                             abs(matflip(isig1,isig3,isig5,3+flav1,3+flav3,3+flav5,1))**2
                           endif
                           m2scc(flav1,flav3,flav5,3,0) = m2scc(flav1,flav3,flav5,3,0)
     &                          + m2stmp(0)*fac*facWupper
                                                                                 
!     lower line
                           if(tfac(2).ne.0d0) then
                              m2stmp = 0
                              if (.not.lonlyinterference) then
!     normal
                                 m2stmp =
     &                                colsq* ! colour factor
     &                                abs(mat(isig1,isig3,isig5,3+flav1,3+flav3,3+flav5,2))**2
!     flipped
                                 if ( lsymmcontrib .and. (flav3-1)
     $                                /2.eq.(flav5-1)/2 ) then
                                    m2stmp(1) = m2stmp(1)
     &                                   + colsq* ! colour factor
     &                                   abs(matflip(isig1,isig3,isig5,3+flav1,3+flav3,3+flav5,2))**2
                                 endif
                              endif
                              if (lsymmcontrib) then
                                 if (isig3.eq.isig5) then
!     interference
                                    m2stmp(1) = m2stmp(1)
     &                                   + 2*colint*dreal(
     &                                   mat(isig1,isig3,isig5,3+flav1,3+flav3,3+flav5,2)*
     &                                   dconjg(matflip(isig1,isig3,isig5,3+flav1,3+flav3,3+flav5,2))
     &                                   )
                                 endif
                              else
                                 if (isig3.eq.isig5) then
                                    m2stmp(1) = m2stmp(1)
     &                                   + colint*dreal(
     &                                   mat(isig1,isig3,isig5,3+flav1,3+flav3,3+flav5,2)*
     &                                   dconjg(matflip(isig1,isig3,isig5,3+flav1,3+flav3,3+flav5,2))
     &                                   )
                                 endif
                              endif
                              m2scc(flav1,flav3,flav5,2,:) = m2scc(flav1,flav3,flav5,2,:)
     &                             + m2stmp(:)*fac*facWlower
                           endif
                           m2stmp = 0
                           if (.not. lonlyinterference ) then
                              m2stmp(0) = m2stmp(0)
     &                             + colsq* ! colour factor
     &                             abs(matflip(isig1,isig3,isig5,3+flav1,3+flav3,3+flav5,2))**2
                           endif
                           m2scc(flav1,flav3,flav5,4,0) = m2scc(flav1,flav3,flav5,4,0)
     &                          + m2stmp(0)*fac*facWlower
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      endif
c     eliminate any higgsstrahlung graphs
      if (fsign(5).eq.fsign(6)) then !either q5 or q6 in i.s.
         if (fsign(1).eq.-fsign(2) .and. 
     1        fsign(3).eq.fsign(4) ) then
c     cf(2,2) = 0.0d0         !  initial gluon attached to 1-2 line only
            m2scc(:,:,:,2,:)= 0d0
         elseif (fsign(1).eq.fsign(2) .and. 
     1           fsign(3).eq.-fsign(4) ) then
c     cf(1,1)= 0.0d0          !  initial gluon attached to 3-4 line only
            m2scc(:,:,:,1,:)= 0d0
         else 
            stop
         endif
      endif

      

      return
      end

!-----------------------------------------------------
      subroutine qqh4q_gtobraket(numf,numg,psi,pbar,gcurr,gmom,psig,pmomg)

      implicit none

c     include 'VBFNLO/utilities/global.inc'

      real * 8 pi,pi2
      parameter (pi=3.141592653589793238462643383279502884197D0,
     1     pi2=pi**2)
      include 'pwhg_st.h'

      include 'global.inc'

      integer, intent(in) :: numf, numg
      double precision, intent(in) :: pbar(0:3,6+max_v), gmom(0:4,-1:1,7)
      double complex, intent(in) :: psi(2,-1:1,6), gcurr(1:6,-1:1,7)
      double precision, intent(out) :: pmomg(0:4,6,7)
      double complex, intent(out) :: psig(2,-1:1,6,-1:1,7)
      integer isigf,isigg
      double precision minusgmom(0:3)

      minusgmom = -gmom(0:3,-1,numg) ! adjusting conventions
      do isigf=-1,1,2
         do isigg=-1,1,2
            if (mod(numf,2) .eq. 0) then
               call bra2c(psi(1,isigf,numf),.true.,pbar(0,numf),isigf,
     &              minusgmom,gcurr(1,isigg,numg),
     &              psig(1,isigf,numf,isigg,numg),pmomg(0,numf,numg))
            else
               call ket2c(psi(1,isigf,numf),.true.,pbar(0,numf),isigf,
     &              minusgmom,gcurr(1,isigg,numg),
     &              psig(1,isigf,numf,isigg,numg),pmomg(0,numf,numg))
            endif
         enddo
      enddo
      return
      end
