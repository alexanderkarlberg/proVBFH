!======================================================================
! Comments on coefficient functions
! ---------------------------------------------------------
!
! - F1 = (F2 - FL)/2x
!
! - hep-ph/0504042 (MMV): gives F2 and FL results (electromagnetic)
!
!   (1/x) F_a = C_{a,ns} \otimes q_{ns}
!                   +  <e^2> (C_{a,q} \otimes q_s  + C_{a,g} \otimes g)
!
!   with z = 2,L, where:
!
!   * <e^2> is the average electric charge
!   * q_{ns} = ???? (we're supposed to deduce it from Eq.(4.2))
!   * C_{a,q} = C_{a,ns} + C_{a,ps};
!   * C_{2,ps}^{(0)} = C_{2,ps}^{(1)} = 0 (and presumably for FL too)
!
! - http://www.sciencedirect.com/science/article/pii/055032139290087R#
!   (Zijlstra & van Neerven) has the second-order coefficient
!   functions. 
!
! - from 1109.3717 (Zaro & Co):
!
!   * q_{ns,i}^+ = (q_i + qbar_i) - q_s
!   * q_{ns,i}^- = (q_i - qbar_i) - q_{ns}^v
!   * q_{ns}^v   = \sum_{i=1}^{n_f} (q_i - qbar_i)
!   * q_s        = \sum_{i=1}^{n_f} (q_i + qbar_i)
!
!   That, together with [C_{a,q} = C_{a,ns} + C_{a,ps}] means that the combination
! 
!       q_{ns,j}^+ * C_{i,ns}^+ + q_s * C_{i,q}
!
!   reduces to 
!
!       (q_j+qbar_j) * C_{i,ns}^+ + q_s * C_{i,ps}
!
module matrix_element_dihiggs
  use hoppet
  use parameters
  use tensor
  implicit none

  private
  public :: eval_matrix_element
  public :: eval_matrix_element_tensor
  interface operator(.dot.)
     module procedure dot
  end interface operator(.dot.)

contains

  !----------------------------------------------------------------------
  function eval_matrix_element(order_start,order_stop, x1, x2, P1, P2, q1, q2, &
       pH1, pH2, ptH1H2) result(res)
    use ME_expressions
    integer , intent(in) :: order_start,order_stop
    real(dp), intent(in) :: x1, x2, P1(0:3), P2(0:3), q1(0:3), q2(0:3)
    real(dp), intent(in) :: pH1(0:3), pH2(0:3), ptH1H2
    real(dp)             :: res
    !----------------------------------------------------------------------
    real(dp) :: Q1sq, Q2sq, Q1val, Q2val
    real(dp) :: muR1val, muR2val, muF1val, muF2val
    real(dp) :: Fx1(-6:7,4), Fx2(-6:7,4)
    real(dp) :: F1F1_AA, F1F2_AA, F2F1_AA, F2F2_AA, F3F3_AA
    real(dp) :: F1F1_AB, F1F2_AB, F2F1_AB, F2F2_AB, F3F3_AB
    real(dp) :: F1F1_AC, F1F2_AC, F2F1_AC, F2F2_AC, F3F3_AC
    real(dp) :: F1F1_BB, F1F2_BB, F2F1_BB, F2F2_BB, F3F3_BB
    real(dp) :: F1F1_BC, F1F2_BC, F2F1_BC, F2F2_BC, F3F3_BC
    real(dp) :: F1F1_CC, F1F2_CC, F2F1_CC, F2F2_CC, F3F3_CC
    real(dp) :: WW_norm, ZZ_norm, overall_norm
    real(dp) :: q1q2, P1q1, P1q2, P2q1, P2q2, P1P2
    real(dp) :: q1pH1sq, q1pH2sq, pH1pH2sq
    complex(dp) :: WW_A, ZZ_A, WW_B, ZZ_B, WW_C, ZZ_C
    integer  :: i, j
    logical, parameter :: WpWm = .true., WmWp = .true., ZZ = .true.
    integer :: iorder
    real(dp) :: res2

    Q1sq = -(q1.dot.q1)
    Q2sq = -(q2.dot.q2)
    Q1val = sqrt(Q1sq)
    Q2val = sqrt(Q2sq)
    
    muR1val = muR1(Q1val, Q2val, ptH1H2)
    muR2val = muR2(Q1val, Q2val, ptH1H2)
    muF1val = muF1(Q1val, Q2val, ptH1H2)
    muF2val = muF2(Q1val, Q2val, ptH1H2)

    F1F1_AA = zero
    F2F1_AA = zero
    F1F2_AA = zero
    F2F2_AA = zero
    F3F3_AA = zero

    F1F1_AB = zero
    F2F1_AB = zero
    F1F2_AB = zero
    F2F2_AB = zero
    F3F3_AB = zero
    
    F1F1_AC = zero
    F2F1_AC = zero
    F1F2_AC = zero
    F2F2_AC = zero
    F3F3_AC = zero

    F1F1_BB = zero
    F2F1_BB = zero
    F1F2_BB = zero
    F2F2_BB = zero
    F3F3_BB = zero
    
    F1F1_BC = zero
    F2F1_BC = zero
    F1F2_BC = zero
    F2F2_BC = zero
    F3F3_BC = zero

    F1F1_CC = zero
    F2F1_CC = zero
    F1F2_CC = zero
    F2F2_CC = zero
    F3F3_CC = zero
 
    !VBFHHMOD: We start by implementing the first piece, which is
    !          proportional to g^\mu\nu, and thus equivalent up to
    !          normalisation to single-higgs
    !overall_norm = (GFermi**3)/(S) * four*sqrt(two) ! SINGLE HIGGS

    ! kinematics: pH1 and pH2 are four-vectors of final state Higgs
    ! dsigma = GF^2 Mv^4 / (S (Q1^2 + Mv^2)^2 (Q2^2 + Mv^2)^2)
    !          * W_uv M^ur M^ns W_rs
    ! M^uv = 2 sqrt(2) GF g^uv
    !        * [ 2 Mv^4/(Q1 + pH1)^2 + 2 Mv^4/(Q1 + pH2)^2
    !            + 6 v lambda_HHH Mv^2 / ((pH1 + pH2)^2 - Mh^2) + Mv^2 ]
    !      + (sqrt(2) GF Mv^2 / ((q1 + kh1)^2 - Mv^2))
    !        * (2 pH1^u + q1^u)(-k1^n + k2^n - q1^n)
                  
    overall_norm = four*two*(GFermi**4)/S

    q1pH1sq = ((q1+pH1).dot.(q1+pH1))
    ! I think there might be a typo in the paper, not sure if
    ! this should be q2 instead...
    q1pH2sq = ((q1+pH2).dot.(q1+pH2))
    pH1pH2sq = ((pH1+pH2).dot.(pH1+pH2))

    ! W or Z propagators
    WW_norm = (MW**4) / (((Q1sq + MW**2)**2 + W_WIDTH**2 * MW**2)& 
         &     * ((Q2sq + MW**2)**2 + W_WIDTH**2 * MW**2 ))
    ZZ_norm = (MZ**4) / (((Q1sq + MZ**2)**2 + Z_WIDTH**2 * MZ**2)& 
         &     * ((Q2sq + MZ**2)**2 + Z_WIDTH**2 * MZ**2 ))

    ! the first single-Higgs like piece, labelled A
    ! A = 2 Mv^4/((q1 + pH1)^2 - Mv^2) + 2 Mv^4/((q1 + pH2)^2 - Mv^2)
    !     + 6 v lambda Mv^2/((pH1 + pH2)^2 - Mh^2) + Mv^2
    WW_A = cVVHfact**2 * two*MW**4/complex(q1pH1sq - MW**2,MW&
         &*W_WIDTH) + cVVHfact**2 * two*MW**4/complex(q1pH2sq - MW**2&
         &,MW*W_WIDTH) + 6.0_dp * v_H * cVVHfact * lambda_HHH * (MW&
         &**2) /complex(pH1pH2sq - mh_sq,MH*HWIDTH) + cVVHHfact * MW&
         &**2
    ZZ_A = cVVHfact**2 * two*MZ**4/complex(q1pH1sq - MZ**2,MZ&
         &*Z_WIDTH) + cVVHfact**2 * two*MZ**4/complex(q1pH2sq - MZ**2&
         &,MZ*Z_WIDTH) + 6.0_dp * v_H * cVVHfact * lambda_HHH * (MZ&
         &**2) /complex(pH1pH2sq - mh_sq,MH*HWIDTH) + cVVHHfact * MZ&
         &**2

    
    ! the two new terms, B and C
    WW_B = cVVHfact**2*(MW**2/complex(q1pH1sq - MW**2,MW*W_WIDTH))/two &
         &           *MW**2/complex(MW**2,-MW*W_WIDTH)
    ZZ_B = cVVHfact**2*(MZ**2/complex(q1pH1sq - MZ**2,MZ*Z_WIDTH))/two &
         &           *MZ**2/complex(MZ**2,-MZ*Z_WIDTH)
    WW_C = cVVHfact**2*(MW**2/complex(q1pH2sq - MW**2,MW*W_WIDTH))/two &
         &           *MW**2/complex(MW**2,-MW*W_WIDTH)
    ZZ_C = cVVHfact**2*(MZ**2/complex(q1pH2sq - MZ**2,MZ*Z_WIDTH))/two &
         &           *MZ**2/complex(MZ**2,-MZ*Z_WIDTH)
    
    ! Below is a narrow-width propagator; This was used in earlier versions of the code
    ! ZZ_norm =  MZ**8 / ((Q1sq + MZ**2)**2 * (Q2sq + MZ**2)**2)
    ! WW_norm =  MW**8 / ((Q1sq + MW**2)**2 * (Q2sq + MW**2)**2)

    ! Compute the LO structure funtion by adding all the pieces
    ! from tables
    Fx1(:,1) = two*F_LO(x1, Q1val, muR1val, muF1val)
    Fx2(:,1) = two*F_LO(x2, Q2val, muR2val, muF2val)

    if (order_stop.ge.2) then
       ! Compute the NLO structure funtion by adding all the pieces
       ! from tables
       Fx1(:,2) = two*F_NLO(x1, Q1val, muR1val, muF1val)
       Fx2(:,2) = two*F_NLO(x2, Q2val, muR2val, muF2val)
    endif

    if (order_stop.ge.3) then
       ! Compute the NNLO structure funtion by adding all the pieces
       ! from tables
       Fx1(:,3) = two*F_NNLO(x1, Q1val, muR1val, muF1val)
       Fx2(:,3) = two*F_NNLO(x2, Q2val, muR2val, muF2val)
    endif

    if (order_stop.ge.4) then
       ! Compute the N3LO structure funtion by adding all the pieces
       ! from tables
       Fx1(:,4) = two*F_N3LO(x1, Q1val, muR1val, muF1val)
       Fx2(:,4) = two*F_N3LO(x2, Q2val, muR2val, muF2val)
    endif

    do iorder = order_start,order_stop
       do i = 1, iorder
          j = 1 + iorder - i
          if (WpWm) then
             F1F1_AA = F1F1_AA + WW_norm * real( WW_A * conjg(WW_A) ) * Fx1(iF1Wp,i)*Fx2(iF1Wm,j)
             F2F1_AA = F2F1_AA + WW_norm * real( WW_A * conjg(WW_A) ) * Fx1(iF2Wp,i)*Fx2(iF1Wm,j)
             F1F2_AA = F1F2_AA + WW_norm * real( WW_A * conjg(WW_A) ) * Fx1(iF1Wp,i)*Fx2(iF2Wm,j)
             F2F2_AA = F2F2_AA + WW_norm * real( WW_A * conjg(WW_A) ) * Fx1(iF2Wp,i)*Fx2(iF2Wm,j)
             F3F3_AA = F3F3_AA + WW_norm * real( WW_A * conjg(WW_A) ) * Fx1(iF3Wp,i)*Fx2(iF3Wm,j)

             F1F1_AB = F1F1_AB + WW_norm * real( WW_A * conjg(WW_B) ) * Fx1(iF1Wp,i)*Fx2(iF1Wm,j)
             F2F1_AB = F2F1_AB + WW_norm * real( WW_A * conjg(WW_B) ) * Fx1(iF2Wp,i)*Fx2(iF1Wm,j)
             F1F2_AB = F1F2_AB + WW_norm * real( WW_A * conjg(WW_B) ) * Fx1(iF1Wp,i)*Fx2(iF2Wm,j)
             F2F2_AB = F2F2_AB + WW_norm * real( WW_A * conjg(WW_B) ) * Fx1(iF2Wp,i)*Fx2(iF2Wm,j)
             F3F3_AB = F3F3_AB + WW_norm * real( WW_A * conjg(WW_B) ) * Fx1(iF3Wp,i)*Fx2(iF3Wm,j)

             F1F1_AC = F1F1_AC + WW_norm * real( WW_A * conjg(WW_C) ) * Fx1(iF1Wp,i)*Fx2(iF1Wm,j)
             F2F1_AC = F2F1_AC + WW_norm * real( WW_A * conjg(WW_C) ) * Fx1(iF2Wp,i)*Fx2(iF1Wm,j)
             F1F2_AC = F1F2_AC + WW_norm * real( WW_A * conjg(WW_C) ) * Fx1(iF1Wp,i)*Fx2(iF2Wm,j)
             F2F2_AC = F2F2_AC + WW_norm * real( WW_A * conjg(WW_C) ) * Fx1(iF2Wp,i)*Fx2(iF2Wm,j)
             F3F3_AC = F3F3_AC + WW_norm * real( WW_A * conjg(WW_C) ) * Fx1(iF3Wp,i)*Fx2(iF3Wm,j)

             F1F1_BB = F1F1_BB + WW_norm * real( WW_B * conjg(WW_B) ) * Fx1(iF1Wp,i)*Fx2(iF1Wm,j)
             F2F1_BB = F2F1_BB + WW_norm * real( WW_B * conjg(WW_B) ) * Fx1(iF2Wp,i)*Fx2(iF1Wm,j)
             F1F2_BB = F1F2_BB + WW_norm * real( WW_B * conjg(WW_B) ) * Fx1(iF1Wp,i)*Fx2(iF2Wm,j)
             F2F2_BB = F2F2_BB + WW_norm * real( WW_B * conjg(WW_B) ) * Fx1(iF2Wp,i)*Fx2(iF2Wm,j)
             F3F3_BB = F3F3_BB + WW_norm * real( WW_B * conjg(WW_B) ) * Fx1(iF3Wp,i)*Fx2(iF3Wm,j)

             F1F1_BC = F1F1_BC + WW_norm * real( WW_B * conjg(WW_C) ) * Fx1(iF1Wp,i)*Fx2(iF1Wm,j)
             F2F1_BC = F2F1_BC + WW_norm * real( WW_B * conjg(WW_C) ) * Fx1(iF2Wp,i)*Fx2(iF1Wm,j)
             F1F2_BC = F1F2_BC + WW_norm * real( WW_B * conjg(WW_C) ) * Fx1(iF1Wp,i)*Fx2(iF2Wm,j)
             F2F2_BC = F2F2_BC + WW_norm * real( WW_B * conjg(WW_C) ) * Fx1(iF2Wp,i)*Fx2(iF2Wm,j)
             F3F3_BC = F3F3_BC + WW_norm * real( WW_B * conjg(WW_C) ) * Fx1(iF3Wp,i)*Fx2(iF3Wm,j)
             
             F1F1_CC = F1F1_CC + WW_norm * real( WW_C * conjg(WW_C) ) * Fx1(iF1Wp,i)*Fx2(iF1Wm,j)
             F2F1_CC = F2F1_CC + WW_norm * real( WW_C * conjg(WW_C) ) * Fx1(iF2Wp,i)*Fx2(iF1Wm,j)
             F1F2_CC = F1F2_CC + WW_norm * real( WW_C * conjg(WW_C) ) * Fx1(iF1Wp,i)*Fx2(iF2Wm,j)
             F2F2_CC = F2F2_CC + WW_norm * real( WW_C * conjg(WW_C) ) * Fx1(iF2Wp,i)*Fx2(iF2Wm,j)
             F3F3_CC = F3F3_CC + WW_norm * real( WW_C * conjg(WW_C) ) * Fx1(iF3Wp,i)*Fx2(iF3Wm,j)
          end if

          if (WmWp) then
             F1F1_AA = F1F1_AA + WW_norm * real( WW_A * conjg(WW_A) ) * Fx1(iF1Wm,i)*Fx2(iF1Wp,j)
             F2F1_AA = F2F1_AA + WW_norm * real( WW_A * conjg(WW_A) ) * Fx1(iF2Wm,i)*Fx2(iF1Wp,j)
             F1F2_AA = F1F2_AA + WW_norm * real( WW_A * conjg(WW_A) ) * Fx1(iF1Wm,i)*Fx2(iF2Wp,j)
             F2F2_AA = F2F2_AA + WW_norm * real( WW_A * conjg(WW_A) ) * Fx1(iF2Wm,i)*Fx2(iF2Wp,j)
             F3F3_AA = F3F3_AA + WW_norm * real( WW_A * conjg(WW_A) ) * Fx1(iF3Wm,i)*Fx2(iF3Wp,j)

             F1F1_AB = F1F1_AB + WW_norm * real( WW_A * conjg(WW_B) ) * Fx1(iF1Wm,i)*Fx2(iF1Wp,j)
             F2F1_AB = F2F1_AB + WW_norm * real( WW_A * conjg(WW_B) ) * Fx1(iF2Wm,i)*Fx2(iF1Wp,j)
             F1F2_AB = F1F2_AB + WW_norm * real( WW_A * conjg(WW_B) ) * Fx1(iF1Wm,i)*Fx2(iF2Wp,j)
             F2F2_AB = F2F2_AB + WW_norm * real( WW_A * conjg(WW_B) ) * Fx1(iF2Wm,i)*Fx2(iF2Wp,j)
             F3F3_AB = F3F3_AB + WW_norm * real( WW_A * conjg(WW_B) ) * Fx1(iF3Wm,i)*Fx2(iF3Wp,j)

             F1F1_AC = F1F1_AC + WW_norm * real( WW_A * conjg(WW_C) ) * Fx1(iF1Wm,i)*Fx2(iF1Wp,j)
             F2F1_AC = F2F1_AC + WW_norm * real( WW_A * conjg(WW_C) ) * Fx1(iF2Wm,i)*Fx2(iF1Wp,j)
             F1F2_AC = F1F2_AC + WW_norm * real( WW_A * conjg(WW_C) ) * Fx1(iF1Wm,i)*Fx2(iF2Wp,j)
             F2F2_AC = F2F2_AC + WW_norm * real( WW_A * conjg(WW_C) ) * Fx1(iF2Wm,i)*Fx2(iF2Wp,j)
             F3F3_AC = F3F3_AC + WW_norm * real( WW_A * conjg(WW_C) ) * Fx1(iF3Wm,i)*Fx2(iF3Wp,j)

             F1F1_BB = F1F1_BB + WW_norm * real( WW_B * conjg(WW_B) ) * Fx1(iF1Wm,i)*Fx2(iF1Wp,j)
             F2F1_BB = F2F1_BB + WW_norm * real( WW_B * conjg(WW_B) ) * Fx1(iF2Wm,i)*Fx2(iF1Wp,j)
             F1F2_BB = F1F2_BB + WW_norm * real( WW_B * conjg(WW_B) ) * Fx1(iF1Wm,i)*Fx2(iF2Wp,j)
             F2F2_BB = F2F2_BB + WW_norm * real( WW_B * conjg(WW_B) ) * Fx1(iF2Wm,i)*Fx2(iF2Wp,j)
             F3F3_BB = F3F3_BB + WW_norm * real( WW_B * conjg(WW_B) ) * Fx1(iF3Wm,i)*Fx2(iF3Wp,j)

             F1F1_BC = F1F1_BC + WW_norm * real( WW_B * conjg(WW_C) ) * Fx1(iF1Wm,i)*Fx2(iF1Wp,j)
             F2F1_BC = F2F1_BC + WW_norm * real( WW_B * conjg(WW_C) ) * Fx1(iF2Wm,i)*Fx2(iF1Wp,j)
             F1F2_BC = F1F2_BC + WW_norm * real( WW_B * conjg(WW_C) ) * Fx1(iF1Wm,i)*Fx2(iF2Wp,j)
             F2F2_BC = F2F2_BC + WW_norm * real( WW_B * conjg(WW_C) ) * Fx1(iF2Wm,i)*Fx2(iF2Wp,j)
             F3F3_BC = F3F3_BC + WW_norm * real( WW_B * conjg(WW_C) ) * Fx1(iF3Wm,i)*Fx2(iF3Wp,j)

             F1F1_CC = F1F1_CC + WW_norm * real( WW_C * conjg(WW_C) ) * Fx1(iF1Wm,i)*Fx2(iF1Wp,j)
             F2F1_CC = F2F1_CC + WW_norm * real( WW_C * conjg(WW_C) ) * Fx1(iF2Wm,i)*Fx2(iF1Wp,j)
             F1F2_CC = F1F2_CC + WW_norm * real( WW_C * conjg(WW_C) ) * Fx1(iF1Wm,i)*Fx2(iF2Wp,j)
             F2F2_CC = F2F2_CC + WW_norm * real( WW_C * conjg(WW_C) ) * Fx1(iF2Wm,i)*Fx2(iF2Wp,j)
             F3F3_CC = F3F3_CC + WW_norm * real( WW_C * conjg(WW_C) ) * Fx1(iF3Wm,i)*Fx2(iF3Wp,j)
          end if

          if (ZZ) then
             F1F1_AA = F1F1_AA + ZZ_norm * real( ZZ_A * conjg(ZZ_A) ) * Fx1(iF1Z,i)*Fx2(iF1Z,j)
             F2F1_AA = F2F1_AA + ZZ_norm * real( ZZ_A * conjg(ZZ_A) ) * Fx1(iF2Z,i)*Fx2(iF1Z,j)
             F1F2_AA = F1F2_AA + ZZ_norm * real( ZZ_A * conjg(ZZ_A) ) * Fx1(iF1Z,i)*Fx2(iF2Z,j)
             F2F2_AA = F2F2_AA + ZZ_norm * real( ZZ_A * conjg(ZZ_A) ) * Fx1(iF2Z,i)*Fx2(iF2Z,j)
             F3F3_AA = F3F3_AA + ZZ_norm * real( ZZ_A * conjg(ZZ_A) ) * Fx1(iF3Z,i)*Fx2(iF3Z,j)

             F1F1_AB = F1F1_AB + ZZ_norm * real( ZZ_A * conjg(ZZ_B) ) * Fx1(iF1Z,i)*Fx2(iF1Z,j)
             F2F1_AB = F2F1_AB + ZZ_norm * real( ZZ_A * conjg(ZZ_B) ) * Fx1(iF2Z,i)*Fx2(iF1Z,j)
             F1F2_AB = F1F2_AB + ZZ_norm * real( ZZ_A * conjg(ZZ_B) ) * Fx1(iF1Z,i)*Fx2(iF2Z,j)
             F2F2_AB = F2F2_AB + ZZ_norm * real( ZZ_A * conjg(ZZ_B) ) * Fx1(iF2Z,i)*Fx2(iF2Z,j)
             F3F3_AB = F3F3_AB + ZZ_norm * real( ZZ_A * conjg(ZZ_B) ) * Fx1(iF3Z,i)*Fx2(iF3Z,j)

             F1F1_AC = F1F1_AC + ZZ_norm * real( ZZ_A * conjg(ZZ_C) ) * Fx1(iF1Z,i)*Fx2(iF1Z,j)
             F2F1_AC = F2F1_AC + ZZ_norm * real( ZZ_A * conjg(ZZ_C) ) * Fx1(iF2Z,i)*Fx2(iF1Z,j)
             F1F2_AC = F1F2_AC + ZZ_norm * real( ZZ_A * conjg(ZZ_C) ) * Fx1(iF1Z,i)*Fx2(iF2Z,j)
             F2F2_AC = F2F2_AC + ZZ_norm * real( ZZ_A * conjg(ZZ_C) ) * Fx1(iF2Z,i)*Fx2(iF2Z,j)
             F3F3_AC = F3F3_AC + ZZ_norm * real( ZZ_A * conjg(ZZ_C) ) * Fx1(iF3Z,i)*Fx2(iF3Z,j)

             F1F1_BB = F1F1_BB + ZZ_norm * real( ZZ_B * conjg(ZZ_B) ) * Fx1(iF1Z,i)*Fx2(iF1Z,j)
             F2F1_BB = F2F1_BB + ZZ_norm * real( ZZ_B * conjg(ZZ_B) ) * Fx1(iF2Z,i)*Fx2(iF1Z,j)
             F1F2_BB = F1F2_BB + ZZ_norm * real( ZZ_B * conjg(ZZ_B) ) * Fx1(iF1Z,i)*Fx2(iF2Z,j)
             F2F2_BB = F2F2_BB + ZZ_norm * real( ZZ_B * conjg(ZZ_B) ) * Fx1(iF2Z,i)*Fx2(iF2Z,j)
             F3F3_BB = F3F3_BB + ZZ_norm * real( ZZ_B * conjg(ZZ_B) ) * Fx1(iF3Z,i)*Fx2(iF3Z,j)

             F1F1_BC = F1F1_BC + ZZ_norm * real( ZZ_B * conjg(ZZ_C) ) * Fx1(iF1Z,i)*Fx2(iF1Z,j)
             F2F1_BC = F2F1_BC + ZZ_norm * real( ZZ_B * conjg(ZZ_C) ) * Fx1(iF2Z,i)*Fx2(iF1Z,j)
             F1F2_BC = F1F2_BC + ZZ_norm * real( ZZ_B * conjg(ZZ_C) ) * Fx1(iF1Z,i)*Fx2(iF2Z,j)
             F2F2_BC = F2F2_BC + ZZ_norm * real( ZZ_B * conjg(ZZ_C) ) * Fx1(iF2Z,i)*Fx2(iF2Z,j)
             F3F3_BC = F3F3_BC + ZZ_norm * real( ZZ_B * conjg(ZZ_C) ) * Fx1(iF3Z,i)*Fx2(iF3Z,j)

             F1F1_CC = F1F1_CC + ZZ_norm * real( ZZ_C * conjg(ZZ_C) ) * Fx1(iF1Z,i)*Fx2(iF1Z,j)
             F2F1_CC = F2F1_CC + ZZ_norm * real( ZZ_C * conjg(ZZ_C) ) * Fx1(iF2Z,i)*Fx2(iF1Z,j)
             F1F2_CC = F1F2_CC + ZZ_norm * real( ZZ_C * conjg(ZZ_C) ) * Fx1(iF1Z,i)*Fx2(iF2Z,j)
             F2F2_CC = F2F2_CC + ZZ_norm * real( ZZ_C * conjg(ZZ_C) ) * Fx1(iF2Z,i)*Fx2(iF2Z,j)
             F3F3_CC = F3F3_CC + ZZ_norm * real( ZZ_C * conjg(ZZ_C) ) * Fx1(iF3Z,i)*Fx2(iF3Z,j)
          end if
       end do
    end do

    res = zero

    ! ! the single-Higgs like component, the (g^mn)^2 piece of the ME
    ! q1q2 = q1 .dot. q2
    ! P1q1 = P1 .dot. q1
    ! P1q2 = P1 .dot. q2
    ! P2q1 = P2 .dot. q1
    ! P2q2 = P2 .dot. q2
    ! P1P2 = P1 .dot. P2
    ! res = res + F1F1_AA*(two + q1q2**2/(Q1sq*Q2sq))
    ! res = res + F1F2_AA/P2q2 * (P2q2**2/(-Q2sq) + (P2q1 - P2q2*q1q2/(-Q2sq))**2/(-Q1sq))
    ! res = res + F2F1_AA/P1q1 * (P1q1**2/(-Q1sq) + (P1q2 - P1q1*q1q2/(-Q1sq))**2/(-Q2sq))
    ! res = res + F2F2_AA/(P1q1*P2q2)*(P1P2 - P1q1*P2q1/(-Q1sq) - &
    !      &            P1q2*P2q2/(-Q2sq) + P1q1*P2q2*q1q2/(Q1sq*Q2sq))**2
    ! res = res + F3F3_AA/(two*P1q1*P2q2)*(P1P2*q1q2 - P1q2*P2q1)
    ! res2 = res*overall_norm
    ! res=zero
    
    ! the full ME (ie not only the piece proportional to g^mn g^mn)
    ! start with F1F1
    res = res + F1F1(F1F1_AA,F1F1_AB,F1F1_AC,F1F1_BB,F1F1_BC,F1F1_CC, &
         &           q1, q2, P1, P2, pH1, pH2)
    res = res + F1F2(F1F2_AA,F1F2_AB,F1F2_AC,F1F2_BB,F1F2_BC,F1F2_CC, &
         &           q1, q2, P1, P2, pH1, pH2)
    res = res + F2F1(F2F1_AA,F2F1_AB,F2F1_AC,F2F1_BB,F2F1_BC,F2F1_CC, &
         &           q1, q2, P1, P2, pH1, pH2)
    res = res + F2F2(F2F2_AA,F2F2_AB,F2F2_AC,F2F2_BB,F2F2_BC,F2F2_CC, &
         &           q1, q2, P1, P2, pH1, pH2)
    res = res + F3F3(F3F3_AA,F3F3_AB,F3F3_AC,F3F3_BB,F3F3_BC,F3F3_CC, &
         &           q1, q2, P1, P2, pH1, pH2)

    ! add the overall normalisation by (2 sqrt(2) Gf^2 / sqrt(s) )^2
    res = res * overall_norm

    !write(*,*) res2, res2/res
  end function eval_matrix_element

  !----------------------------------------------------------------------
  ! The same matrix element as eval_matrix_element, computed with
  ! numerical tensor contractions. This is the reference
  ! implementation: unlike eval_matrix_element it includes all F3
  ! interference terms, also those that arise from the imaginary parts
  ! of the propagators with finite widths.
  !
  ! The hadronic tensor of beam b at order i is
  !   W_b^(i) = sum_k Fk_b^(i) G_b(k),
  ! with the basis tensors G_b(1:3) of hadronic_basis, which depend
  ! neither on the order nor on the boson. The contraction
  ! Tr[(W_1 M)(M^* W_2)] is therefore computed once per phase-space
  ! point and boson for each pair of basis tensors (trace_matrix), and
  ! the orders are combined as in eval_matrix_element, i.e. summing
  ! the products of beam-1 order i and beam-2 order j with
  ! i + j = n + 1 for n = order_start..order_stop (order_sum).
  function eval_matrix_element_tensor(order_start,order_stop, x1, x2, P1, P2, q1, q2, &
       pH1, pH2, ptH1H2) result(res)
    integer , intent(in) :: order_start,order_stop
    real(dp), intent(in) :: x1, x2, P1(0:3), P2(0:3), q1(0:3), q2(0:3)
    real(dp), intent(in) :: ptH1H2,pH1(0:3),pH2(0:3)
    real(dp)             :: res
    !----------------------------------------------------------------------
    real(dp) :: Q1sq, Q2sq, Q1val, Q2val
    real(dp) :: muR1val, muR2val, muF1val, muF2val
    real(dp) :: Fx1(-6:7,4), Fx2(-6:7,4)
    real(dp) :: WW_norm, ZZ_norm, overall_norm
    real(dp) :: TW(3,3), TZ(3,3)
    logical, parameter :: WpWm = .true., WmWp = .true., ZZ = .true.
    integer, parameter :: Wp=1, Wm=-1, Z=0
    integer :: iWp(3), iWm(3), iZ(3)
    type(tensors) :: M(-1:1), Mstar(-1:1), G1(3), G2(3)
    type(tensors) :: g_mu_nu
    real(dp) :: sigma

    ! Structure-function indices of F1, F2, F3 for each boson
    iWp = (/ iF1Wp, iF2Wp, iF3Wp /)
    iWm = (/ iF1Wm, iF2Wm, iF3Wm /)
    iZ  = (/ iF1Z,  iF2Z,  iF3Z  /)

    if(.not.gmunu%initialised) then
       call SetMetric(1)
    endif
    g_mu_nu = gmunu
    g_mu_nu%up = .false. ! Lower all the indices

    ! Compute the VVHH currents
    M(Wp) = VVtoHH_tensor(Wp, q1, pH1, pH2, g_mu_nu)
    M(Z)  = VVtoHH_tensor(Z,  q1, pH1, pH2, g_mu_nu)

    ! We need to raise the indices of M. The matrix element is
    ! identical between W+ and W-, so M(Wm) is not needed.
    call raise(M(Wp),1) ! Raise first index
    call raise(M(Wp),2) ! Raise second index
    call raise(M(Z),1) ! Raise first index
    call raise(M(Z),2) ! Raise second index

    Mstar(Wp) = M(Wp)
    Mstar(Wp)%values = dconjg(M(Wp)%values) ! Complex conjugate
    Mstar(Z) = M(Z)
    Mstar(Z)%values = dconjg(M(Z)%values) ! Complex conjugate

    ! Contract the basis tensors of the two hadronic tensors with the
    ! matrix elements
    call hadronic_basis(P1, q1, G1)
    call hadronic_basis(P2, q2, G2)
    call trace_matrix(G1, G2, M(Wp), Mstar(Wp), TW)
    call trace_matrix(G1, G2, M(Z), Mstar(Z), TZ)

    Q1sq = -(q1 .dot. q1)
    Q2sq = -(q2 .dot. q2)
    Q1val = sqrt(Q1sq)
    Q2val = sqrt(Q2sq)
    
    muR1val = muR1(Q1val, Q2val, ptH1H2)
    muR2val = muR2(Q1val, Q2val, ptH1H2)
    muF1val = muF1(Q1val, Q2val, ptH1H2)
    muF2val = muF2(Q1val, Q2val, ptH1H2)

    ! Compute the overall numerical factors
    overall_norm = (GFermi**4)/(S) * two
    WW_norm =  MW**8 / (((Q1sq + MW**2)**2 + W_WIDTH**2 * MW**2)& 
         & * ((Q2sq + MW**2)**2 + W_WIDTH**2 * MW**2 ))
    ZZ_norm =  MZ**8 / (((Q1sq + MZ**2)**2 + Z_WIDTH**2 * MZ**2)& 
         & * ((Q2sq + MZ**2)**2 + Z_WIDTH**2 * MZ**2 ))

    ! Compute the LO structure funtion by adding all the pieces
    ! from tables
    Fx1(:,1) = two*F_LO(x1, Q1val, muR1val, muF1val)
    Fx2(:,1) = two*F_LO(x2, Q2val, muR2val, muF2val)

    if (order_stop.ge.2) then
       ! Compute the NLO structure funtion by adding all the pieces
       ! from tables
       Fx1(:,2) = two*F_NLO(x1, Q1val, muR1val, muF1val)
       Fx2(:,2) = two*F_NLO(x2, Q2val, muR2val, muF2val)
    endif

    if (order_stop.ge.3) then
       ! Compute the NNLO structure funtion by adding all the pieces
       ! from tables
       Fx1(:,3) = two*F_NNLO(x1, Q1val, muR1val, muF1val)
       Fx2(:,3) = two*F_NNLO(x2, Q2val, muR2val, muF2val)
    endif

    if (order_stop.ge.4) then
       ! Compute the N3LO structure funtion by adding all the pieces
       ! from tables
       Fx1(:,4) = two*F_N3LO(x1, Q1val, muR1val, muF1val)
       Fx2(:,4) = two*F_N3LO(x2, Q2val, muR2val, muF2val)
    endif

    ! Compute the 3 different contributions
    sigma = zero
    if (WpWm) sigma = sigma + WW_norm * order_sum(order_start, order_stop, Fx1, Fx2, iWp, iWm, TW)
    if (WmWp) sigma = sigma + WW_norm * order_sum(order_start, order_stop, Fx1, Fx2, iWm, iWp, TW)
    if (ZZ)   sigma = sigma + ZZ_norm * order_sum(order_start, order_stop, Fx1, Fx2, iZ,  iZ,  TZ)

    res = overall_norm * sigma
  end function eval_matrix_element_tensor

  !----------------------------------------------------------------------
  ! The VV -> HH current M_mu_nu (both indices down) for V = W (V=+-1)
  ! or Z (V=0), with the coupling modifiers cVVHfact, cVVHHfact and
  ! lambdafact. As VVtoHH_tensor in proVBFHH, which in addition has
  ! switches for the individual diagrams and returns them separately
  ! for the non-factorisable corrections.
  function VVtoHH_tensor(V, q1, pH1, pH2, g_mu_nu) result(res)
    integer, intent(in) :: V
    real(dp), intent(in) :: q1(0:3), pH1(0:3), pH2(0:3)
    type(tensors), intent(in) :: g_mu_nu ! The metric with all indices down
    type(tensors) :: res
    type(tensors) :: q1mu, k1mu, k2mu, VV_H_HH, VVHH, VHVHVt, VHVHVu
    real(dp) :: q1k1sq, q1k2sq, k1k2sq, MV, MVsq, V_WIDTH

    ! Copy the four-vectors into tensor types, with lower indices
    call InitFourVector(q1mu,q1,.false.)
    call InitFourVector(k1mu,pH1,.false.)
    call InitFourVector(k2mu,pH2,.false.)

    !Compute all combinations of dot-products
    k1k2sq = (pH1 + pH2) .dot. (pH1 + pH2)
    q1k1sq = (q1 + pH1) .dot. (q1 + pH1)
    q1k2sq = (q1 + pH2) .dot. (q1 + pH2)

    if(V.eq.0) then ! Z
       MV = MZ
       V_WIDTH = Z_WIDTH
    elseif(abs(V).eq.1) then ! W+/W-
       MV = MW
       V_WIDTH = W_WIDTH
    else
       stop 'Wrong boson in VVtoHH'
    endif
    MVsq = MV**2

    ! VV -> H -> HH part
    VV_H_HH = 6.0_dp * cVVHfact * lambdafact * MH**2/cmplx(k1k2sq - MH**2,MH*HWIDTH,kind=dp)*g_mu_nu
    ! Quartic vertex part
    VVHH = two * cVVHHfact * g_mu_nu 
    ! Double VBF part
    VHVHVt = two*((two*MVsq)/cmplx(q1k1sq-MVsq,MV*V_WIDTH,kind=dp))*g_mu_nu &
         & + (one/cmplx(q1k1sq -MVsq,MV*V_WIDTH,kind=dp)) &
         &         *MVsq/cmplx(MVsq,-MV*V_WIDTH,kind=dp) &
         &         *((two*k1mu+q1mu).otimes.(k2mu-k1mu-q1mu))
    VHVHVt = cVVHfact**2 * VHVHVt
    VHVHVu = two*((two*MVsq)/cmplx(q1k2sq-MVsq,MV*V_WIDTH,kind=dp))*g_mu_nu &
         & + (one/cmplx(q1k2sq -MVsq,MV*V_WIDTH,kind=dp)) &
         &         *MVsq/cmplx(MVsq,-MV*V_WIDTH,kind=dp) &
         &         *((two*k2mu+q1mu).otimes.(k1mu-k2mu-q1mu))
    VHVHVu = cVVHfact**2 * VHVHVu

    res = VV_H_HH + VVHH + VHVHVt + VHVHVu
  end function VVtoHH_tensor

  !----------------------------------------------------------------------
  ! Basis tensors of the hadronic tensor of a beam with momentum P and
  ! momentum transfer q (all indices down):
  !   W_mu_nu = F1 G(1) + F2 G(2) + F3 G(3), with
  !   G(1) = q_mu q_nu / q^2 - g_mu_nu
  !   G(2) = Phat_mu Phat_nu / (P.q),  Phat = P - (P.q)/q^2 q
  !   G(3) = i epsilon_mu_nu_rho_sigma P^rho q^sigma / (2 P.q)
  subroutine hadronic_basis(P, q, G)
    real(dp), intent(in) :: P(0:3), q(0:3)
    type(tensors), intent(inout) :: G(3)
    type(tensors) :: g_mu_nu, qmu, Pmu, Phatmu
    real(dp) :: qq, Pq
    integer :: i, j

    g_mu_nu = gmunu
    g_mu_nu%up = .false. ! Lower all the indices
    call InitFourVector(qmu,q,.false.)
    call InitFourVector(Pmu,P,.false.)
    qq = q .dot. q
    Pq = P .dot. q
    Phatmu = Pmu - Pq/qq * qmu

    G(1) = (one/qq)*(qmu.otimes.qmu) - g_mu_nu
    G(2) = (one/Pq)*(Phatmu.otimes.Phatmu)

    ! We perform the contraction between the levi-civita tensor and
    ! P and q explicitly, as implementing rank-4 tensors just for
    ! the levi-civita symbol is tedious. These are the lower-index
    ! components epsilon_mu_nu_rho_sigma P^rho q^sigma, with
    ! epsilon_0123 = +1 and P, q the contravariant components.
    call InitTensor(G(3),2,.false.)
    G(3)%values(0,1) =   P(2)*q(3) - P(3)*q(2)
    G(3)%values(0,2) = - P(1)*q(3) + P(3)*q(1)
    G(3)%values(0,3) =   P(1)*q(2) - P(2)*q(1)
    G(3)%values(1,2) = - P(3)*q(0) + P(0)*q(3)
    G(3)%values(1,3) =   P(2)*q(0) - P(0)*q(2)
    G(3)%values(2,3) = - P(1)*q(0) + P(0)*q(1)
    ! Use anti-symmetric property    
    do i = 0,3
       do j = i,3
          G(3)%values(j,i) = - G(3)%values(i,j)
       enddo
    enddo
    ! And now mulitply with overall factor i/(2 P.q)
    G(3)%values(:,:) = cmplx(zero,one/(two*Pq),kind=dp) * G(3)%values(:,:)
  end subroutine hadronic_basis

  !----------------------------------------------------------------------
  ! T(k,l) = Re[ G1(k)_{mu nu} G2(l)_{rho sigma} Ma^{mu rho} Mbstar^{nu sigma} ]
  ! (Ma and Mbstar with both indices up, G1 and G2 with both down),
  ! i.e. the coefficient of Fk(beam 1) Fl(beam 2) in the squared
  ! matrix element for Ma = M, Mbstar = M^*. Since the G are hermitian,
  ! T is real in that case; for Ma /= Mb only the real part is kept,
  ! as in the interference terms of the non-factorisable corrections.
  !
  ! This is the bulk of the work, so it is done directly on the
  ! component arrays: with X(k)^nu_rho = G1(k)_{mu nu} Ma^{mu rho}
  ! and Y(l)^nu_rho = Mbstar^{nu sigma} G2(l)_{rho sigma},
  ! T(k,l) = Re sum_{nu,rho} X(k)^nu_rho Y(l)^nu_rho. This is what
  ! ContractTensors(G1(k),1,Ma,1,X), ContractTensors(Mbstar,2,G2(l),2,Y),
  ! ContractTensors(X,1,Y,1,XY) and TensorTrace(XY) would give, at
  ! about a third of the cost.
  subroutine trace_matrix(G1, G2, Ma, Mbstar, T)
    type(tensors), intent(in) :: G1(3), G2(3), Ma, Mbstar
    real(dp), intent(out) :: T(3,3)
    complex(dp) :: X(0:3,0:3,3), Y(0:3,0:3,3)
    integer :: k, l

    ! The index positions assumed above
    if(.not.(all(Ma%up).and.all(Mbstar%up)).or.Ma%rank.ne.2.or.Mbstar%rank.ne.2) then
       stop 'trace_matrix: Ma and Mbstar must be rank 2 with both indices up'
    endif
    do k = 1, 3
       if(any(G1(k)%up).or.any(G2(k)%up).or.G1(k)%rank.ne.2.or.G2(k)%rank.ne.2) then
          stop 'trace_matrix: G1 and G2 must be rank 2 with both indices down'
       endif
    enddo

    do k = 1, 3
       X(:,:,k) = matmul(transpose(G1(k)%values), Ma%values)
       Y(:,:,k) = matmul(Mbstar%values, transpose(G2(k)%values))
    enddo
    do l = 1, 3
       do k = 1, 3
          T(k,l) = real(sum(X(:,:,k)*Y(:,:,l)), kind=dp)
       enddo
    enddo
  end subroutine trace_matrix

  !----------------------------------------------------------------------
  ! sum_{n=order_start}^{order_stop} sum_{i+j=n+1}
  !     sum_{k,l} Fk(beam 1, order i) Fl(beam 2, order j) T(k,l)
  ! where i1(k) and i2(l) are the indices of F1, F2, F3 for the bosons
  ! attached to beams 1 and 2.
  real(dp) function order_sum(order_start, order_stop, Fx1, Fx2, i1, i2, T) result(res)
    integer,  intent(in) :: order_start, order_stop, i1(3), i2(3)
    real(dp), intent(in) :: Fx1(-6:7,4), Fx2(-6:7,4), T(3,3)
    integer :: iorder, i, j, k, l

    res = zero
    do iorder = order_start, order_stop
       do i = 1, iorder
          j = 1 + iorder - i
          do l = 1, 3
             do k = 1, 3
                res = res + Fx1(i1(k),i) * Fx2(i2(l),j) * T(k,l)
             enddo
          enddo
       enddo
    enddo
  end function order_sum

  !----------------------------------------------------------------------
  ! dot product
  real(dp) function dot(p1,p2)
    real(dp), intent(in) :: p1(0:3), p2(0:3)
    dot = p1(0)*p2(0) - sum(p1(1:3)*p2(1:3))
  end function dot

  !----------------------------------------------------------------------
  ! mu_R1 as a function of Q1 and Q2
  real(dp) function muR1(Q1, Q2, ptH1H2)
    real(dp), intent(in) :: Q1, Q2, ptH1H2
    muR1 = zero
    if (scale_choice.le.1) then
       ! if scale_choice = 0,1 then muR1(Q1,Q2) = muR(Q1)
       muR1 = sf_muR(Q1)
    elseif (scale_choice.eq.2) then
       ! else if scale_choice=2, use sqrt(Q1*Q2)
       muR1 = xmur * sqrt(Q1 * Q2)
    elseif (scale_choice.eq.3) then
       ! else if scale_choice=3, use mixed scale
       muR1 = xmur * mixed_scale(Q1, Q2, ptH1H2)
    endif
  end function muR1
  
  !----------------------------------------------------------------------
  ! mu_R2 as a function of Q1 and Q2
  real(dp) function muR2(Q1, Q2, ptH1H2)
    real(dp), intent(in) :: Q1, Q2, ptH1H2
    muR2 = zero
    if (scale_choice.le.1) then
       ! if scale_choice = 0,1 then muR1(Q1,Q2) = muR(Q1)
       muR2 = sf_muR(Q2)
    elseif (scale_choice.eq.2) then
       ! else if scale_choice=2, use sqrt(Q1*Q2)
       muR2 = xmur * sqrt(Q1 * Q2)
    elseif (scale_choice.eq.3) then
       ! else if scale_choice=3, use mixed scale
       muR2 = xmur * mixed_scale(Q1, Q2, ptH1H2)
    endif
  end function muR2

  !----------------------------------------------------------------------
  ! mu_F1 as a function of Q1 and Q2
  real(dp) function muF1(Q1, Q2, ptH1H2)
    real(dp), intent(in) :: Q1, Q2, ptH1H2
    muF1 = zero
    if (scale_choice.le.1) then
       ! if scale_choice = 0,1 then muF1(Q1,Q2) = muF(Q1)
       muF1 = sf_muF(Q1)
    elseif (scale_choice.eq.2) then
       ! else if scale_choice=2, use sqrt(Q1*Q2)
       muF1 = xmuf * sqrt(Q1 * Q2)
    elseif (scale_choice.eq.3) then
       ! else if scale_choice=3, use mixed scale
       muF1 = xmuf * mixed_scale(Q1, Q2, ptH1H2)
    else
       call wae_error('muF1(Q)', 'illegal value for scale_choice', intval = scale_choice)
    endif
  end function muF1

  !----------------------------------------------------------------------
  ! mu_F2 as a function of Q1 and Q2
  real(dp) function muF2(Q1, Q2, ptH1H2)
    real(dp), intent(in) :: Q1, Q2, ptH1H2
    muF2 = zero
    if (scale_choice.le.1) then
       ! if scale_choice = 0,1 then muF1(Q1,Q2) = muF(Q1)
       muF2 = sf_muF(Q2)
    elseif (scale_choice.eq.2) then
       ! else if scale_choice=2, use sqrt(Q1*Q2)
       muF2 = xmuf * sqrt(Q1 * Q2)
    elseif (scale_choice.eq.3) then
       ! else if scale_choice=3, use mixed scale
       muF2 = xmuf * mixed_scale(Q1, Q2, ptH1H2)
    else
       call wae_error('muF2(Q)', 'illegal value for scale_choice', intval = scale_choice)
    endif
  end function muF2

  !----------------------------------------------------------------------
  ! Defines which scale to use for scale_choice = 3,
  ! which can be any function of Q1, Q2 and ptH1H2
  real(dp) function mixed_scale(Q1,Q2,ptH1H2)
    real(dp), intent(in) :: Q1, Q2, ptH1H2
    mixed_scale = ((mh*0.5d0)**4d0+(mh*ptH1H2*0.5d0)**2d0)**(0.25d0)
  end function mixed_scale

end module matrix_element_dihiggs
