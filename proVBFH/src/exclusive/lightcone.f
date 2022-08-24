!     Fortran code for transformation to lightcone coordinates
!     and mapping of k1,k2,k3 to Born kinematics
      
!----------------------------------------------------------------------
!     Subroutine to go from Minkowski to lightcone coordinates
!     Given a 5-vector p_minkowski in format (px, py, pz, E , m)
!     transform it to  p_lightcone in format (px, py, p-, p+, m)
!     where p+- = 1/sqrt(2) (E +- pz)
      subroutine lightcone(p_minkowski,p_lightcone)
      implicit none
      real*8, intent(in)  :: p_minkowski(5)
      real*8, intent(out) :: p_lightcone(5)
      !real * 8 p_minkowski(5), p_lightcone(5)
!     p^+ = (1/sqrt(2)) (E + p_z)
      p_lightcone(4) = (1.d0/sqrt(2.d0))*(p_minkowski(4) 
     $     + p_minkowski(3))
!     p^- = (1/sqrt(2)) (E - p_z)
      p_lightcone(3) = (1.d0/sqrt(2.d0))*(p_minkowski(4)
     $     - p_minkowski(3))
      p_lightcone(1) = p_minkowski(1)
      p_lightcone(2) = p_minkowski(2)
      p_lightcone(5) = p_minkowski(5)
      end

!----------------------------------------------------------------------
!     Subroutine to go from Minkowski to lightcone coordinates
!     Given a 5-vector p_lightcone in format (px, py, p-, p+, m)
!     transform it to  p_minkowski in format (px, py, pz, E , m)
!     where p+- = 1/sqrt(2) (E +- pz)
      subroutine minkowski(p_lightcone, p_minkowski)
      implicit none
      real*8, intent(in)  :: p_lightcone(5)
      real*8, intent(out) :: p_minkowski(5)
! NB: watch out for confusion here between p_{0123} (Exyz) and
!     the fortran encoding which is (xyzE)
!     E = (1/sqrt(2)) (p^+ + p^-)
      p_minkowski(4) = (1.d0/sqrt(2.d0))*(p_lightcone(4)
     $     + p_lightcone(3))
!     p_z = (1/sqrt(2)) (p^+ - p^-)
      p_minkowski(3) = (1.d0/sqrt(2.d0))*(p_lightcone(4)
     $     - p_lightcone(3))
      p_minkowski(1) = p_lightcone(1)
      p_minkowski(2) = p_lightcone(2)
      p_minkowski(5) = p_lightcone(5)
      end

!----------------------------------------------------------------------
!     work out the squared mass from a lightcone-coordinates vector
!     (intended for tests, since mass is available as 5th component)
      function m2_lightcone(p_lc)
      implicit none
      real * 8 p_lc(5)
      real * 8 m2_lightcone
      m2_lightcone = -2.d0*p_lc(3)*p_lc(4) + p_lc(1)**2 + p_lc(2)**2
      return
      end function

!----------------------------------------------------------------------
!     work out the squared mass from a minkowski-coordinates vector
!     (intended for tests, since mass is available as 5th component)
      function m2_minkowski(p_m)
      implicit none
      real * 8 p_m(5)
      real * 8 m2_minkowski
      m2_minkowski = -p_m(4)**2 + p_m(1)**2 + p_m(2)**2 + p_m(3)**2
      return
      end function
      
!----------------------------------------------------------------------
!     Mapping from NNLO to Born for both gluons on upper line
!     k1 is the incoming quark, k2 is the outgoing quark and k3 and
!     k4 the radiated gluons
      subroutine mapping_1_1(k1, k2, k3, k4, p1, p2)
      real * 8 k1(5), k2(5), k3(5), k4(5), p1(5), p2(5)
      real * 8 k1_lc(5), k2_lc(5), k3_lc(5), p1_lc(5), p2_lc(5)
      call mapping_1(k1,k2,k3+k4,p1,p2)

      end

!----------------------------------------------------------------------
!     Mapping from NNLO to Born for both gluons on lower line
!     k1 is the incoming quark, k2 is the outgoing quark and k3 and
!     k4 the radiated gluons
      subroutine mapping_2_2(k1, k2, k3, k4, p1, p2)
      real * 8 k1(5), k2(5), k3(5), k4(5), p1(5), p2(5)
      real * 8 k1_lc(5), k2_lc(5), k3_lc(5), p1_lc(5), p2_lc(5)
      call mapping_2(k1,k2,k3+k4,p1,p2)

      end
!----------------------------------------------------------------------
!     Mapping from NNLO to Born for gluons on opposite lines
      subroutine mapping_1_2(ku1,ku2,ku3,kl1,kl2,kl3,pu1,pu2,pl1,pl2)
c     ku1,ku2,ku3 are the input momenta from the upper quark line and
C     pu1, pu2 the output momenta for the upper line
      real * 8 ku1(5), ku2(5), ku3(5), pu1(5), pu2(5)
c     kl1,kl2,kl3 are the input momenta from the lower quark line and
C     pl1, pl2 the output momenta for the lower line
      real * 8 kl1(5), kl2(5), kl3(5), pl1(5), pl2(5)
      call mapping_1(ku1, ku2, ku3, pu1, pu2)
      call mapping_2(kl1, kl2, kl3, pl1, pl2)
      end
!----------------------------------------------------------------------
!     Mapping from NLO to Born for gluon on upper line
!     k1 is the incoming quark, k2 is the outgoing quark and k3 the
!     radiated gluon
      subroutine mapping_1(k1, k2, k3, p1, p2)
      real * 8 m2_lightcone
      external m2_lightcone 
      real * 8 k1(5), k2(5), k3(5), p1(5), p2(5)
      real * 8 k1_lc(5), k2_lc(5), k3_lc(5), p1_lc(5), p2_lc(5)
      real * 8 x
!     k1 - k2 - k3 = p1 - p2
!     k1 = (0, 0, 0, k1^+)
!     k2 = (k2_x , k2_y , k2^- , k2^+ )
!     k3 = (k3_x , k3_y , k3^- , k3^+ )
!     p1 = (0    , 0    , 0    , p1^+ ) => p1^+ = k1^+ -k2^+ -k3^+ + p2^+
!     p2 = (p2_x , p2_y , p2^- , x    ) => p2^- = k2^- + k3^-
!     p2**2 = 0 => -2 x (k2^- + k3^-) + p2_x**2 + p2_y**2 = 0
!     x = ((k2_x + k3_x)**2 + (k2_y + k3_y)**2)/(2 (k2^- + k3^-))

      call lightcone(k1,k1_lc)
      call lightcone(k2,k2_lc)
      call lightcone(k3,k3_lc)

      x =  ((k2_lc(1) + k3_lc(1))**2 + (k3_lc(2) 
     $      + k2_lc(2))**2)/(2.d0*(k2_lc(3)+k3_lc(3)))

      p1_lc(1) = 0.d0
      p1_lc(2) = 0.d0
      p1_lc(3) = 0.d0
      p1_lc(4) = k1_lc(4) - k2_lc(4) - k3_lc(4) + x
      p1_lc(5) = m2_lightcone(p1_lc)

      p2_lc(1) = k2_lc(1) + k3_lc(1)
      p2_lc(2) = k2_lc(2) + k3_lc(2)
      p2_lc(3) = k2_lc(3) + k3_lc(3)
      p2_lc(4) = x
      p2_lc(5) = m2_lightcone(p2_lc)

      call minkowski(p1_lc,p1)
      call minkowski(p2_lc,p2)

      end

!----------------------------------------------------------------------
!     Mapping from NLO to Born for gluon on lower line
!     k1 is the incoming quark, k2 is the outgoing quark and k3 the
!     radiated gluon
      subroutine mapping_2(k1, k2, k3, p1, p2)
      real * 8 m2_lightcone
      external m2_lightcone
      real * 8 k1(5), k2(5), k3(5), p1(5), p2(5)
      real * 8 k1_lc(5), k2_lc(5), k3_lc(5), p1_lc(5), p2_lc(5)
      real * 8 x
!     k1 - k2 - k3 = p1 - p2
!     k1 = (0    ,0    ,k1^- ,0   )
!     k2 = (k2_x ,k2_y ,k2^- ,k2^+)
!     k3 = (k3_x ,k3_y ,k3^- ,k3^+)
!     p1 = (0    ,0    ,p1^- ,0   ) => p1^- = k1^- - k2^- - k3^- + x
!     p2 = (p2_x ,p2_y ,x    ,p2^+) => p2^+ = k2^+ + k3^+
!     p2**2 = 0 => -2x(k2^+ + k3^+) + p2_x**2 + p2_y**2 = 0
!     x = ( (k2_x + k3_x)**2 + (k2_y + k3_y)**2 ) / (2 (k2^+ + k3^+))
      call lightcone(k1,k1_lc)
      call lightcone(k2,k2_lc)
      call lightcone(k3,k3_lc)

      x =  ((k2_lc(1) + k3_lc(1))**2 + (k3_lc(2) 
     $      + k2_lc(2))**2)/(2.d0*(k2_lc(4)+k3_lc(4)))

      p1_lc(1) = 0.d0
      p1_lc(2) = 0.d0
      p1_lc(4) = 0.d0
      p1_lc(3) = k1_lc(3) - k2_lc(3) - k3_lc(3) + x
      p1_lc(5) = m2_lightcone(p1_lc)

      p2_lc(1) = k2_lc(1) + k3_lc(1)
      p2_lc(2) = k2_lc(2) + k3_lc(2)
      p2_lc(4) = k2_lc(4) + k3_lc(4)
      p2_lc(3) = x
      p2_lc(5) = m2_lightcone(p2_lc)

      call minkowski(p1_lc,p1)
      call minkowski(p2_lc,p2)

      end



