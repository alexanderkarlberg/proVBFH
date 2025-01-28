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
  use hoppet_v1
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
  ! This function returns the same output as eval_matrix_element, but
  ! is done with numerical tensor manipulations. It is significantly
  ! slower than the analytically computed routine above.
  function eval_matrix_element_tensor(order_start,order_stop, x1, x2, P1, P2, q1, q2, &
       pH1, pH2, ptH1H2) result(res)
    integer , intent(in) :: order_start,order_stop
    real(dp), intent(in) :: x1, x2, P1(0:3), P2(0:3), q1(0:3), q2(0:3)
    real(dp), intent(in) :: ptH1H2,pH1(0:3),pH2(0:3)
    real(dp)             :: res
    !----------------------------------------------------------------------
    real(dp) :: Q1sq, Q2sq, Q1val, Q2val
    real(dp) :: muR1val, muR2val, muF1val, muF2val
    real(dp) :: Fx1(-6:7,4), Fx2(-6:7,4), F1(-6:7), F2(-6:7)
    real(dp) :: WW_norm, ZZ_norm, overall_norm
    real(dp) :: P1q1, P2q2, q1q1, q2q2, k1k2sq, q1k1sq, q1k2sq
    integer  :: i, j
    logical, parameter :: WpWm = .true., WmWp = .true., ZZ = .true.
    integer, parameter :: Wp=1, Wm=-1, Z=0
    integer :: iorder
    type(tensors) :: M(-1:1),Mstar(-1:1), Wx1(-1:1), Wx2(-1:1)
    type(tensors) :: g_mu_nu, q1mu, q2mu, P1mu, P2mu, k1mu, k2mu
    type(tensors) ::  P1hatmu, P2hatmu, dummy1, dummy2, T3(2), MWx1(-1:1), MstarWx2(-1:1)
    real * 8 :: sigma, sigmaWpWm, sigmaWmWp, sigmaZZ

    ! Start by initialising all the tensors needed
    if(.not.gmunu%initialised) then
       call SetMetric(1)
    endif
    
    do i = -1,1
       if(.not.Wx1(i)%initialised) then
          call InitTensor(Wx1(i),2,.false.)
       else
          call ResetTensor(Wx1(i))
       endif
       if(.not.Wx2(i)%initialised) then
          call InitTensor(Wx2(i),2,.false.)
       else
          call ResetTensor(Wx2(i))
       endif
    enddo
    
    if(.not.q1mu%initialised) then
       call InitTensor(q1mu,1,.false.)
    else
       call ResetTensor(q1mu)
    endif
    if(.not.q2mu%initialised) then
       call InitTensor(q2mu,1,.false.)
    else
       call ResetTensor(q2mu)
    endif
    
    if(.not.k1mu%initialised) then
       call InitTensor(k1mu,1,.false.)
    else
       call ResetTensor(k1mu)
    endif
    if(.not.k2mu%initialised) then
       call InitTensor(k2mu,1,.false.)
    else
       call ResetTensor(k2mu)
    endif
    
    if(.not.P1mu%initialised) then
       call InitTensor(P1mu,1,.false.)
    else
       call ResetTensor(P1mu)
    endif
    if(.not.P2mu%initialised) then
       call InitTensor(P2mu,1,.false.)
    else
       call ResetTensor(P2mu)
    endif

    if(.not.T3(1)%initialised) then
       call InitTensor(T3(1),2,.false.)
    else
       call ResetTensor(T3(1))
    endif
    if(.not.T3(2)%initialised) then
       call InitTensor(T3(2),2,.false.)
    else
       call ResetTensor(T3(2))
    endif

    
    g_mu_nu = gmunu
    g_mu_nu%up = .false. ! Lower all the indices

    ! Copy the four-vectors into tensor types
    q1mu%values(:,1) = q1(:)
    q2mu%values(:,1) = q2(:)
    k1mu%values(:,1) = pH1(:)
    k2mu%values(:,1) = pH2(:)
    P1mu%values(:,1) = P1(:)
    P2mu%values(:,1) = P2(:)

    !Compute all combinations of dot-products
    q1q1 = q1 .dot. q1
    q2q2 = q2 .dot. q2
    P1q1 = P1 .dot. q1
    P2q2 = P2 .dot. q2
    k1k2sq = (pH1 + pH2) .dot. (pH1 + pH2)
    q1k1sq = (q1 + pH1) .dot. (q1 + pH1)
    q1k2sq = (q1 + pH2) .dot. (q1 + pH2)


    P1hatmu = P1mu - P1q1/q1q1 * q1mu 
    P2hatmu = P2mu - P2q2/q2q2 * q2mu

    M(Wp) = (two*((two*MW**2)/complex(q1k1sq-MW**2,MW*W_WIDTH) &
         & + (two*MW**2)/complex(q1k2sq-MW**2,MW*W_WIDTH) & 
         & + three*MH**2/complex(k1k2sq - MH**2,MH*HWIDTH) + one))*g_mu_nu ! gmunu
                                                                               ! contribution
    M(Wp) = M(Wp) + (one/complex(q1k1sq -MW**2,MW*W_WIDTH)) &
         &           *MW**2/complex(MW**2,-MW*W_WIDTH) &
         &           *((two*k1mu+q1mu).otimes.(k2mu-k1mu-q1mu))
    M(Wp) = M(Wp) + (one/complex(q1k2sq -MW**2,MW*W_WIDTH)) &
         &           *MW**2/complex(MW**2,-MW*W_WIDTH) &
         &           *((two*k2mu+q1mu).otimes.(k1mu-k2mu-q1mu))
    M(Wm) = M(Wp) ! Matrix element is identical between W+/W-

    M(Z) = (two*((two*MZ**2)/complex(q1k1sq-MZ**2,MZ*Z_WIDTH) &
         & + (two*MZ**2)/complex(q1k2sq-MZ**2,MZ*Z_WIDTH) & 
         & + three*MH**2/complex(k1k2sq - MH**2,MH*HWIDTH) + one))*g_mu_nu ! gmunu
                                                                               ! contribution
    M(Z) = M(Z) + (one/complex(q1k1sq -MZ**2,MZ*Z_WIDTH)) &
         &         *MZ**2/complex(MZ**2,-MZ*Z_WIDTH) &
         &         *((two*k1mu+q1mu).otimes.(k2mu-k1mu-q1mu))
    M(Z) = M(Z) + (one/complex(q1k2sq -MZ**2,MZ*Z_WIDTH)) &
         &         *MZ**2/complex(MZ**2,-MZ*Z_WIDTH) &
         &         *((two*k2mu+q1mu).otimes.(k1mu-k2mu-q1mu))


    ! We need to raise the indeces of M and Mstar
    call raise(M(Wp),1) ! Raise first index
    call raise(M(Wp),2) ! Raise second index
    call raise(M(Wm),1) ! Raise first index
    call raise(M(Wm),2) ! Raise second index
    call raise(M(Z),1) ! Raise first index
    call raise(M(Z),2) ! Raise second index
    
    Mstar = M ! For initialisation
    Mstar(Wp)%values = dconjg(M(Wp)%values) ! Complex conjugate
    Mstar(Wm)%values = dconjg(M(Wm)%values) ! Complex conjugate
    Mstar(Z)%values = dconjg(M(Z)%values) ! Complex conjugate
    

    ! We perform the contraction between the levi-civita tensor and
    ! P_i and q_i explicitly, as implementing rank-4 tensors just for
    ! the levi-civita symbol is tedious. So T3(i) is given by
    ! epsilon_mu_nu_rho_sigma * P_i^rho * q_i^sigma

    T3(1)%values(0,1) =   P1mu%values(2,1)*q1mu%values(3,1) &
         &              - P1mu%values(3,1)*q1mu%values(2,1)
    T3(1)%values(0,2) = - P1mu%values(1,1)*q1mu%values(3,1) &
         &              + P1mu%values(3,1)*q1mu%values(1,1) 
    T3(1)%values(0,3) =   P1mu%values(1,1)*q1mu%values(2,1) &
         &              - P1mu%values(2,1)*q1mu%values(1,1) 
    T3(1)%values(1,2) = - P1mu%values(3,1)*q1mu%values(0,1) &
         &              + P1mu%values(0,1)*q1mu%values(3,1) 
    T3(1)%values(1,3) =   P1mu%values(2,1)*q1mu%values(0,1) &
         &              - P1mu%values(0,1)*q1mu%values(2,1) 
    T3(1)%values(2,3) = - P1mu%values(1,1)*q1mu%values(0,1) &
         &              + P1mu%values(0,1)*q1mu%values(1,1) 

    T3(2)%values(0,1) =   P2mu%values(2,1)*q2mu%values(3,1) &
         &              - P2mu%values(3,1)*q2mu%values(2,1)
    T3(2)%values(0,2) = - P2mu%values(1,1)*q2mu%values(3,1) &
         &              + P2mu%values(3,1)*q2mu%values(1,1) 
    T3(2)%values(0,3) =   P2mu%values(1,1)*q2mu%values(2,1) &
         &              - P2mu%values(2,1)*q2mu%values(1,1) 
    T3(2)%values(1,2) = - P2mu%values(3,1)*q2mu%values(0,1) &
         &              + P2mu%values(0,1)*q2mu%values(3,1) 
    T3(2)%values(1,3) =   P2mu%values(2,1)*q2mu%values(0,1) &
         &              - P2mu%values(0,1)*q2mu%values(2,1) 
    T3(2)%values(2,3) = - P2mu%values(1,1)*q2mu%values(0,1) &
         &              + P2mu%values(0,1)*q2mu%values(1,1) 

    ! Use anti-symmetric property    
    do i = 0,3
       do j = i,3
          T3(1)%values(j,i) = - T3(1)%values(i,j)
          T3(2)%values(j,i) = - T3(2)%values(i,j)
       enddo
    enddo

    ! And now mulitply with overall factor i
    T3(1)%values(:,:) = complex(zero,one) * T3(1)%values(:,:)
    T3(2)%values(:,:) = complex(zero,one) * T3(2)%values(:,:)
    
    Q1sq = -q1q1
    Q2sq = -q2q2
    Q1val = sqrt(Q1sq)
    Q2val = sqrt(Q2sq)
    
    muR1val = muR1(Q1val, Q2val, ptH1H2)
    muR2val = muR2(Q1val, Q2val, ptH1H2)
    muF1val = muF1(Q1val, Q2val, ptH1H2)
    muF2val = muF2(Q1val, Q2val, ptH1H2)

    !These will contain the full structure functions for the two protons
    F1 = zero
    F2 = zero
    
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
    F1(:) = Fx1(:,1)
    F2(:) = Fx2(:,1)

    if (order_stop.ge.2) then
       ! Compute the NLO structure funtion by adding all the pieces
       ! from tables
       Fx1(:,2) = two*F_NLO(x1, Q1val, muR1val, muF1val)
       Fx2(:,2) = two*F_NLO(x2, Q2val, muR2val, muF2val)
       F1(:) = F1(:) + Fx1(:,2)
       F2(:) = F2(:) + Fx2(:,2)
    endif

    if (order_stop.ge.3) then
       ! Compute the NNLO structure funtion by adding all the pieces
       ! from tables
       Fx1(:,3) = two*F_NNLO(x1, Q1val, muR1val, muF1val)
       Fx2(:,3) = two*F_NNLO(x2, Q2val, muR2val, muF2val)
       F1(:) = F1(:) + Fx1(:,3)
       F2(:) = F2(:) + Fx2(:,3)
    endif

    if (order_stop.ge.4) then
       ! Compute the N3LO structure funtion by adding all the pieces
       ! from tables
       Fx1(:,4) = two*F_N3LO(x1, Q1val, muR1val, muF1val)
       Fx2(:,4) = two*F_N3LO(x2, Q2val, muR2val, muF2val)
       F1(:) = F1(:) + Fx1(:,4)
       F2(:) = F2(:) + Fx2(:,4)
    endif

    ! Compute hadronic tensors
    Wx1(Wp) = F1(iF1Wp)*((one/q1q1)*(q1mu.otimes.q1mu)-g_mu_nu) &
         & + F1(iF2Wp)*(one/P1q1)*(P1hatmu.otimes.P1hatmu) &
         & + F1(iF3Wp)*(one/(two*P1q1))*T3(1)

    Wx1(Wm) = F1(iF1Wm)*((one/q1q1)*(q1mu.otimes.q1mu)-g_mu_nu) &
         & + F1(iF2Wm)*(one/P1q1)*(P1hatmu.otimes.P1hatmu) &
         & + F1(iF3Wm)*(one/(two*P1q1))*T3(1)

    Wx1(Z) = F1(iF1Z)*((one/q1q1)*(q1mu.otimes.q1mu)-g_mu_nu) &
         & + F1(iF2Z)*(one/P1q1)*(P1hatmu.otimes.P1hatmu) &
         & + F1(iF3Z)*(one/(two*P1q1))*T3(1)

    Wx2(Wp) = F2(iF1Wp)*((one/q2q2)*(q2mu.otimes.q2mu)-g_mu_nu) &
         & + F2(iF2Wp)*(one/P2q2)*(P2hatmu.otimes.P2hatmu) &
         & + F2(iF3Wp)*(one/(two*P2q2))*T3(2)

    Wx2(Wm) = F2(iF1Wm)*((one/q2q2)*(q2mu.otimes.q2mu)-g_mu_nu) &
         & + F2(iF2Wm)*(one/P2q2)*(P2hatmu.otimes.P2hatmu) &
         & + F2(iF3Wm)*(one/(two*P2q2))*T3(2)

    Wx2(Z) = F2(iF1Z)*((one/q2q2)*(q2mu.otimes.q2mu)-g_mu_nu) &
         & + F2(iF2Z)*(one/P2q2)*(P2hatmu.otimes.P2hatmu) &
         & + F2(iF3Z)*(one/(two*P2q2))*T3(2)

    ! Do first contraction between a hadronic tensor and a matrix element.
    call ContractTensors(Wx1(Wp),1,M(Wp),1,MWx1(Wp))
    call ContractTensors(Wx1(Wm),1,M(Wm),1,MWx1(Wm))
    call ContractTensors(Wx1(Z),1,M(Z),1,MWx1(Z))

    call ContractTensors(Mstar(Wp),2,Wx2(Wp),2,MstarWx2(Wp))
    call ContractTensors(Mstar(Wm),2,Wx2(Wm),2,MstarWx2(Wm))
    call ContractTensors(Mstar(Z),2,Wx2(Z),2,MstarWx2(Z))

    ! Compute the 3 different contributions
    if (WpWm) then
       call ContractTensors(MWx1(Wp),1,MstarWx2(Wm),1,Dummy1)
       sigmaWpWm = WW_norm * TensorTrace(Dummy1)
    end if
    
    if (WmWp) then
       call ContractTensors(MWx1(Wm),1,MstarWx2(Wp),1,Dummy1)
       sigmaWmWp = WW_norm * TensorTrace(Dummy1)
    end if
    
    if (ZZ) then
       call ContractTensors(MWx1(Z),1,MstarWx2(Z),1,Dummy1)
       sigmaZZ = ZZ_norm * TensorTrace(Dummy1)
    end if
 
    sigma = overall_norm * (sigmaWpWm + sigmaWmWp + sigmaZZ)
    res = sigma
  end function eval_matrix_element_tensor

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
