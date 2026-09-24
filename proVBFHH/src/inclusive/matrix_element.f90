!======================================================================
! List of issues relative to 1109.3717
! ------------------------------------
!
! - product of sums in Eq.(3.19) can't be right
!
! - getting W- from q_ns(-) + q_s would give 2u + 2dbar, which
!   interacts with W+, not W-; similar issue in F3
!
! - An additional factor of two seems to be needed in Eqs. (3.8),
!   (3.9), (3.17) and (3.18) in order get agreement with Marco
!   Zaro's code. (There he puts it in front of the vi's and ai's)
!
! - Just below Eq.(3.18) they say C_{3,ns}^V = C_{i,ns}^-: should
!   the "i" in the second term actually be a "3"?
!======================================================================
!
! Summary of our understanding of the coefficient functions
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
module matrix_element
  use hoppet
  use incl_parameters
  use streamlined_interface
  use tensor
  implicit none

  private
  public :: eval_matrix_element
  public :: eval_matrix_element_tensor
  public :: eval_matrix_element_vbfnlo
  interface operator(.dot.)
     module procedure dot
  end interface operator(.dot.)
  real(dp), public :: lambda = 1d-5 ! To test explicitly that
                                    ! expressions do not depend on
                                    ! lambda. This value is being
                                    ! overwritten in the below
                                    ! routines.

contains

  !----------------------------------------------------------------------
  function eval_matrix_element(order_start,order_stop, x1, x2, P1, P2, q1, q2, &
       pH1, pH2, ptH1H2) result(res)
    use ME_expressions
    use nonfact
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
    real(dp) :: res2, alphas_lcl,ptj1,ptj2,muR_nonfact

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
    WW_norm = (MW**4) / (((Q1sq + MWsq)**2 + W_WIDTH**2 * MWsq)& 
         &     * ((Q2sq + MWsq)**2 + W_WIDTH**2 * MWsq ))
    ZZ_norm = (MZ**4) / (((Q1sq + MZsq)**2 + Z_WIDTH**2 * MZsq)& 
         &     * ((Q2sq + MZsq)**2 + Z_WIDTH**2 * MZsq ))

    if(non_fact) then
        stop 'Entered eval_matrix_element with non_fact. Cannot proceed.'
    endif
    
    ! the first single-Higgs like piece, labelled A
    ! A = 2 Mv^4/((q1 + pH1)^2 - Mv^2) + 2 Mv^4/((q1 + pH2)^2 - Mv^2)
    !     + 6 v lambda Mv^2/((pH1 + pH2)^2 - Mh^2) + Mv^2
    WW_A = cVVHfact**2*two*MW**4/cmplx(q1pH1sq - MWsq,MW*W_WIDTH,kind=dp) &
         & + cVVHfact**2*two*MW**4/cmplx(q1pH2sq - MWsq,MW*W_WIDTH,kind=dp) &
         & + 6.0_dp * v_H * cVVHfact * lambda_HHH * (MWsq) /cmplx(pH1pH2sq - mh_sq,MH*HWIDTH,kind=dp) + cVVHHfact*MWsq
    ZZ_A = cVVHfact**2*two*MZ**4/cmplx(q1pH1sq - MZsq,MZ*Z_WIDTH,kind=dp) &
         & + cVVHfact**2*two*MZ**4/cmplx(q1pH2sq - MZsq,MZ*Z_WIDTH,kind=dp) &
         & + 6.0_dp * v_H * cVVHfact * lambda_HHH * (MZsq) /cmplx(pH1pH2sq - mh_sq,MH*HWIDTH,kind=dp) + cVVHHfact*MZsq

    
    ! the two new terms, B and C
    WW_B = cVVHfact**2*(MWsq/cmplx(q1pH1sq - MWsq,MW*W_WIDTH,kind=dp))/two &
         &           *MWsq/cmplx(MWsq,-MW*W_WIDTH,kind=dp)
    ZZ_B = cVVHfact**2*(MZsq/cmplx(q1pH1sq - MZsq,MZ*Z_WIDTH,kind=dp))/two &
         &           *MZsq/cmplx(MZsq,-MZ*Z_WIDTH,kind=dp)
    WW_C = cVVHfact**2*(MWsq/cmplx(q1pH2sq - MWsq,MW*W_WIDTH,kind=dp))/two &
         &           *MWsq/cmplx(MWsq,-MW*W_WIDTH,kind=dp)
    ZZ_C = cVVHfact**2*(MZsq/cmplx(q1pH2sq - MZsq,MZ*Z_WIDTH,kind=dp))/two &
         &           *MZsq/cmplx(MZsq,-MZ*Z_WIDTH,kind=dp)
    
    ! Below is a narrow-width propagator; This was used in earlier versions of the code
    ! ZZ_norm =  MZ**8 / ((Q1sq + MZsq)**2 * (Q2sq + MZsq)**2)
    ! WW_norm =  MW**8 / ((Q1sq + MWsq)**2 * (Q2sq + MWsq)**2)

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
  ! numerical tensor contractions. Unlike eval_matrix_element it also
  ! includes the F1 x F3 and F2 x F3 terms that arise from the
  ! imaginary parts of the propagators with finite widths. They are
  ! parity odd, so they integrate to zero for the cross section and
  ! any parity-even distribution; otherwise the two agree to machine
  ! precision. The diagram switches (tri_off, box_t_off, ...) and the
  ! non-factorisable corrections are only available here.
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
    use nonfact
    integer , intent(in) :: order_start,order_stop
    real(dp), intent(in) :: x1, x2, P1(0:3), P2(0:3), q1(0:3), q2(0:3)
    real(dp), intent(in) :: ptH1H2,pH1(0:3),pH2(0:3)
    real(dp)             :: res
    !----------------------------------------------------------------------
    real(dp) :: Q1sq, Q2sq, Q1val, Q2val
    real(dp) :: muR1val, muR2val, muF1val, muF2val
    real(dp) :: Fx1(-6:7,4), Fx2(-6:7,4)
    real(dp) :: WW_norm, ZZ_norm, overall_norm
    real(dp) :: TW(3,3), TZ(3,3), T1(3,3), T2(3,3)
    logical, parameter :: WpWm = .true., WmWp = .true., ZZ = .true.
    integer, parameter :: Wp=1, Wm=-1, Z=0
    integer :: iWp(3), iWm(3), iZ(3), i
    type(tensors) :: M(Z:Wp), Mstar(Z:Wp), G1(3), G2(3), g_mu_nu
    type(tensors) :: M1loop(Z:Wp), M1loopstar(Z:Wp), M2loop(Z:Wp), M2loopstar(Z:Wp)
    ! Tensors to be passed to nonfact loop routines
    type(tensors) :: VV_H_HH(Z:Wp), VVHH(Z:Wp), VHVHVt(Z:Wp), VHVHVu(Z:Wp)
    real(dp) :: sigma, alphas_lcl, muR_nonfact

    ! Structure-function indices of F1, F2, F3 for each boson
    iWp = (/ iF1Wp, iF2Wp, iF3Wp /)
    iWm = (/ iF1Wm, iF2Wm, iF3Wm /)
    iZ  = (/ iF1Z,  iF2Z,  iF3Z  /)

    if(.not.gmunu%initialised) then
       call SetMetric(1)
    endif
    g_mu_nu = gmunu
    g_mu_nu%up = .false. ! Lower all the indices

    do i = Z, Wp
       call InitTensor(VV_H_HH(i),2,.false.)
       call InitTensor(VVHH(i),2,.false.)
       call InitTensor(VHVHVt(i),2,.false.)
       call InitTensor(VHVHVu(i),2,.false.)
    enddo

    Q1sq = -(q1 .dot. q1)
    Q2sq = -(q2 .dot. q2)
    Q1val = sqrt(Q1sq)
    Q2val = sqrt(Q2sq)

    muR1val = muR1(Q1val, Q2val, ptH1H2)
    muR2val = muR2(Q1val, Q2val, ptH1H2)
    muF1val = muF1(Q1val, Q2val, ptH1H2)
    muF2val = muF2(Q1val, Q2val, ptH1H2)

    ! Compute VVHH currents. The matrix element is identical between
    ! W+ and W-, so only Wp and Z are needed.
    M(Wp) = VVtoHH_tensor(Wp,q1,q2,pH1,pH2,g_mu_nu,VV_H_HH(Wp), VVHH(Wp), VHVHVt(Wp), VHVHVu(Wp))
    M(Z) = VVtoHH_tensor(Z,q1,q2,pH1,pH2,g_mu_nu,VV_H_HH(Z), VVHH(Z), VHVHVt(Z), VHVHVu(Z))

    ! Compute the overall numerical factors
    overall_norm = (GFermi**4)/(S) * two
    WW_norm =  MW**8 / (((Q1sq + MWsq)**2 + W_WIDTH**2 * MWsq)& 
         & * ((Q2sq + MWsq)**2 + W_WIDTH**2 * MWsq ))
    ZZ_norm =  MZ**8 / (((Q1sq + MZsq)**2 + Z_WIDTH**2 * MZsq)& 
         & * ((Q2sq + MZsq)**2 + Z_WIDTH**2 * MZsq ))
    
    ! We need to raise the indices of M and Mstar
    do i = Z, Wp
       call raise(M(i),1) ! Raise first index
       call raise(M(i),2) ! Raise second index
       Mstar(i) = M(i)
       Mstar(i)%values = dconjg(M(i)%values) ! Complex conjugate
    enddo

    if(non_fact) then
       if(order_stop.gt.1) stop 'Cannot do non factorisable corrections'
       muR_nonfact = sqrt(muR1val * muR2val)
       alphas_lcl = Value(coupling,muR_nonfact)
       ! Eq. 6 in 1906.10899 for Nc = 3
       overall_norm = overall_norm * 2.0_dp/9.0_dp * alphas_lcl**2

       do i = Z, Wp
          ! Compute 1-loop and 2-loop nonfact VVHH currents
          M1loop(i) = nonfact_1loop_VVtoHH_tensor(i, q1, q2, pH1, pH2, g_mu_nu, VV_H_HH(i), &
               & VVHH(i), VHVHVt(i), VHVHVu(i))
          M2loop(i) = nonfact_2loop_VVtoHH_tensor(i, q1, q2, pH1, pH2, g_mu_nu, VV_H_HH(i), &
               & VVHH(i), VHVHVt(i), VHVHVu(i))
          ! We need to raise the indices of M and Mstar
          call raise(M1loop(i),1)
          call raise(M1loop(i),2)
          call raise(M2loop(i),1)
          call raise(M2loop(i),2)
          M1loopstar(i) = M1loop(i)
          M1loopstar(i)%values = dconjg(M1loop(i)%values)
          M2loopstar(i) = M2loop(i)
          M2loopstar(i)%values = dconjg(M2loop(i)%values)
       enddo
    endif

    ! Contract the basis tensors of the two hadronic tensors with the
    ! matrix elements
    call hadronic_basis(P1, q1, G1)
    call hadronic_basis(P2, q2, G2)
    if(.not.non_fact) then
       call trace_matrix(G1, G2, M(Wp), Mstar(Wp), TW)
       call trace_matrix(G1, G2, M(Z), Mstar(Z), TZ)
    else
       ! 1-loop squared, plus the interference of the 2-loop currents
       ! with the Born
       call trace_matrix(G1, G2, M1loop(Wp), M1loopstar(Wp), TW)
       call trace_matrix(G1, G2, M2loop(Wp), Mstar(Wp), T1)
       call trace_matrix(G1, G2, M(Wp), M2loopstar(Wp), T2)
       TW = TW + T1 + T2
       call trace_matrix(G1, G2, M1loop(Z), M1loopstar(Z), TZ)
       call trace_matrix(G1, G2, M2loop(Z), Mstar(Z), T1)
       call trace_matrix(G1, G2, M(Z), M2loopstar(Z), T2)
       TZ = TZ + T1 + T2
    endif

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
       error stop 'trace_matrix: Ma and Mbstar must be rank 2 with both indices up'
    endif
    do k = 1, 3
       if(any(G1(k)%up).or.any(G2(k)%up).or.G1(k)%rank.ne.2.or.G2(k)%rank.ne.2) then
          error stop 'trace_matrix: G1 and G2 must be rank 2 with both indices down'
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
  
  function VVtoHH_tensor(V, q1, q2, pH1, pH2, g_mu_nu, VV_H_HH, VVHH, VHVHVt,VHVHVu) result(res)
    implicit none
    double precision, intent(in) :: q1(0:3), q2(0:3), pH1(0:3), pH2(0:3)
    integer, intent(in) :: V
    type(tensors), intent(in) :: g_mu_nu ! The metric with all indices down
    type(tensors) :: VV_H_HH, VVHH, VHVHVt, VHVHVu
    type(tensors) :: res
    type(tensors) :: q1mu, q2mu, k1mu, k2mu
    double precision q1k1sq, q1k2sq, k1k2sq, MV, MVsq, V_WIDTH

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

    ! Copy the four-vectors into tensor types. The indices are down,
    ! so the components are lowered with the metric (storing the
    ! contravariant components directly would describe the
    ! parity-flipped vectors, which gave wrong F3 interference terms
    ! between the g^mu_nu and the t/u-channel parts of M)
    call InitFourVector(q1mu,q1,.false.)
    call InitFourVector(q2mu,q2,.false.)
    call InitFourVector(k1mu,pH1,.false.)
    call InitFourVector(k2mu,pH2,.false.)

    !Compute all combinations of dot-products
    k1k2sq = (pH1 + pH2) .dot. (pH1 + pH2)
    q1k1sq = (q1 + pH1) .dot. (q1 + pH1)
    q1k2sq = (q1 + pH2) .dot. (q1 + pH2)

    if(V.eq.0) then ! Z
       MV = MZ
       MVsq = MZsq
       V_WIDTH = Z_WIDTH
    elseif(abs(V).eq.1) then ! W=/W-
       MV = MW
       MVsq = MWsq
       V_WIDTH = W_WIDTH
    else
       stop 'Wrong boson in VVtoHH'
    endif

    if(tri_on) then
       ! VV -> H -> HH part
       VV_H_HH = 6.0_dp * cVVHfact * lambdafact * MH**2/cmplx(k1k2sq - MH**2,MH*HWIDTH,kind=dp)*g_mu_nu
       
       ! Quartic vertex part
       VVHH = two * cVVHHfact * g_mu_nu 
    else
       VV_H_HH%values = (zero,zero)
       VVHH%values = (zero,zero)
    endif
    ! Double VBF part
    if(box_t_on) then
       VHVHVt = two*((two*MVsq)/cmplx(q1k1sq-MVsq,MV*V_WIDTH,kind=dp))*g_mu_nu &
            & + (one/cmplx(q1k1sq -MVsq,MV*V_WIDTH,kind=dp)) &
            &         *MVsq/cmplx(MVsq,-MV*V_WIDTH,kind=dp) &
            &         *((two*k1mu+q1mu).otimes.(k2mu-k1mu-q1mu))
       VHVHVt = cVVHfact**2 * VHVHVt
    else
       VHVHVt%values = (zero,zero)
    endif
    if(box_u_on) then
       VHVHVu = two*((two*MVsq)/cmplx(q1k2sq-MVsq,MV*V_WIDTH,kind=dp))*g_mu_nu &
            & + (one/cmplx(q1k2sq -MVsq,MV*V_WIDTH,kind=dp)) &
            &         *MVsq/cmplx(MVsq,-MV*V_WIDTH,kind=dp) &
            &         *((two*k2mu+q1mu).otimes.(k1mu-k2mu-q1mu))
       VHVHVu = cVVHfact**2 * VHVHVu
    else
       VHVHVu%values = (zero,zero)
    endif
    
    res = VV_H_HH + VVHH + VHVHVt + VHVHVu
  end function VVtoHH_tensor

  function nonfact_1loop_VVtoHH_tensor(V,q1,q2,pH1,pH2,g_mu_nu,VV_H_HH, VVHH, VHVHVt, VHVHVu) &
       & result(res)
    use nonfact
    use nonfact_expressions
    implicit none
    double precision, intent(in) :: q1(0:3), q2(0:3), pH1(0:3), pH2(0:3)
    double precision :: q1rot(0:3), q2rot(0:3), pH1rot(0:3), pH2rot(0:3)
    type(tensors), intent(in) :: g_mu_nu ! The metric with all indices down
    integer, intent(in) :: V
    type(tensors) :: res, VV_H_HH, VVHH, VHVHVt, VHVHVu, VHVHV
    double precision sval, tval, uval
    ! For non factorisable contributions
    double precision pt1,pt2,ptH1,ptH2,qT1(1:2),qT2(1:2),qTH(1:2), MV, MVH2
    double precision pH1m, q1p, pH2m
    double precision, parameter :: z(1:3) = (/zero, zero, one/)
    double precision :: sinphi, cosphi, norm

    if(V.eq.0) then ! Z
       MV = MZ
       MVsq = MZsq
    elseif(abs(V).eq.1) then ! W=/W-
       MV = MW
       MVsq = MWsq
    else
       stop 'Wrong boson in VVtoHH'
    endif

    if(.not.res%initialised) then
       call InitTensor(res,2,.false.)
    else
       call ResetTensor(res)
    endif

    q1rot = q1
    q2rot = q2
    pH1rot = pH1
    pH2rot = pH2
    if(abs(q1(1)).gt.0d0) then
       cosphi = q1(1)/sqrt(q1(1)**2 + q1(2)**2)
       sinphi = sqrt(one - cosphi**2)
       if(q1(2).gt.zero) then
          sinphi = -sinphi
       endif
    else
       cosphi = one
       sinphi = zero
    endif
    call mrotate(z, sinphi, cosphi, q1rot(1:3))
    call mrotate(z, sinphi, cosphi, q2rot(1:3))
    call mrotate(z, sinphi, cosphi, pH1rot(1:3))
    call mrotate(z, sinphi, cosphi, pH2rot(1:3))
    
    qTH(1:2) = - q1(1:2) - q2(1:2)
    pt1 = q1(1)**2 + q1(2)**2
    q1p = (q1(0) + q1(3))/sqrt(2d0)
    pt2 = q2(1)**2 + q2(2)**2
    ptH1 = pH1(1)**2 + pH1(2)**2
    pH1m = (pH1(0) - pH1(3))/sqrt(two)
    pH2m = (pH2(0) - pH2(3))/sqrt(two)
    ptH2 = pH2(1)**2 + pH2(2)**2
    sval = (q1(1)+q2(1))**2+(q1(2)+q2(2))**2
    tval = (q1(1)+pH1(1))**2+(q1(2)+pH1(2))**2
    uval = (q1(1)+pH2(1))**2+(q1(2)+pH2(2))**2

    lambda = nf_regfact * MVsq ! nf_regfact = 1 (default) sets log(lambda/MV^2) = 0

    res%values = (zero,zero)
    if(tri1_on) then
       ! Triangle
       res = VV_H_HH + VVHH 
       res%values = res%values * chi_tri1(pt1,pt2,sval,MV,lambda) * cmplx(zero,one,kind=dp)
    endif
    if(box1_on) then
       ! Box
       ! Fix. Minus in front of 2*ph1m*q1p as we define
       ! q1 = p5 - p1
       MVH2 = MVsq - MH**2 - ptH1 - two * pH1m * q1p
       res%values = res%values + VHVHVt%values * cmplx(zero,one,kind=dp) * &
            box_1loop_new(MV, MVH2, q1rot(1), q2rot(1), q2rot(2), pH1rot(1), pH1rot(2), lambda)
       ! Fix. Minus in front of 2*ph1m*q1p as we define
       ! q1 = p5 - p1
       MVH2 = MVsq - MH**2 - ptH2 - two * pH2m * q1p
       res%values = res%values + VHVHVu%values * cmplx(zero,one,kind=dp) * &
            box_1loop_new(MV, MVH2, q1rot(1), q2rot(1), q2rot(2), pH2rot(1), pH2rot(2), lambda)
    endif
  end function nonfact_1loop_VVtoHH_tensor

  function nonfact_2loop_VVtoHH_tensor(V, q1, q2, pH1, pH2, g_mu_nu, VV_H_HH, VVHH, VHVHVt, VHVHVu) &
       & result(res)
    use nonfact
    implicit none
    double precision, intent(in) :: q1(0:3), q2(0:3), pH1(0:3), pH2(0:3)
    double precision :: q1rot(0:3), q2rot(0:3), pH1rot(0:3), pH2rot(0:3)
    type(tensors), intent(in) :: g_mu_nu ! The metric with all indices down
    integer, intent(in) :: V
    type(tensors) :: res, VV_H_HH, VVHH, VHVHVt, VHVHVu
    ! For non factorisable contributions
    double precision pt1,pt2,sval,ptH1,ptH2,qT1(1:2),qT2(1:2),qTH(1:2), MV, MVH2
    double precision pH1m, q1p, pH2m
    double precision tval, uval
    double precision, parameter :: z(1:3) = (/zero, zero, one/)
    double precision :: sinphi, cosphi, norm
        
    if(V.eq.0) then ! Z
       MV = MZ
       MVsq = MZsq
    elseif(abs(V).eq.1) then ! W=/W-
       MV = MW
       MVsq = MWsq
    else
       stop 'Wrong boson in VVtoHH'
    endif

    if(.not.res%initialised) then
       call InitTensor(res,2,.false.)
    else
       call ResetTensor(res)
    endif
    q1rot = q1
    q2rot = q2
    pH1rot = pH1
    pH2rot = pH2
    if(abs(q1(1)).gt.0d0) then
       cosphi = q1(1)/sqrt(q1(1)**2 + q1(2)**2)
       sinphi = sqrt(one - cosphi**2)
       if(q1(2).gt.zero) then
          sinphi = -sinphi
       endif
    else
       cosphi = one
       sinphi = zero
    endif
    call mrotate(z, sinphi, cosphi, q1rot(1:3))
    call mrotate(z, sinphi, cosphi, q2rot(1:3))
    call mrotate(z, sinphi, cosphi, pH1rot(1:3))
    call mrotate(z, sinphi, cosphi, pH2rot(1:3))

    ! Kinematical quantities
    qTH(1:2) = - q1(1:2) - q2(1:2)
    pt1 = q1(1)**2 + q1(2)**2
    q1p = (q1(0) + q1(3))/sqrt(2d0)
    pt2 = q2(1)**2 + q2(2)**2
    ptH1 = pH1(1)**2 + pH1(2)**2
    pH1m = (pH1(0) - pH1(3))/sqrt(two)
    pH2m = (pH2(0) - pH2(3))/sqrt(two)
    ptH2 = pH2(1)**2 + pH2(2)**2
    sval = (q1(1)+q2(1))**2+(q1(2)+q2(2))**2
    tval = (q1(1)+pH1(1))**2+(q1(2)+pH1(2))**2
    uval = (q1(1)+pH2(1))**2+(q1(2)+pH2(2))**2

    lambda = nf_regfact * MVsq ! nf_regfact = 1 (default) sets log(lambda/MV^2) = 0

    res%values = (zero,zero)
    if(tri2_on) then
       ! Triangle
       res = VV_H_HH + VVHH 
       res%values = - half * res%values * tri_2loop(MV, q1rot(1), q2rot(1), q2rot(2), lambda)
    endif
    if(box2_on) then
       ! Box
       ! Fix
       MVH2 = MVsq - MH**2 - ptH1 - two * pH1m * q1p
       res%values = res%values - half * VHVHVt%values * &
            box_2loop(MV, MVH2, q1rot(1), q2rot(1), q2rot(2), pH1rot(1), pH1rot(2), lambda)
       ! Fix
       MVH2 = MVsq - MH**2 - ptH2 - two * pH2m * q1p
       res%values = res%values - half * VHVHVu%values * &
            box_2loop(MV, MVH2, q1rot(1), q2rot(1), q2rot(2), pH2rot(1), pH2rot(2), lambda)
   endif
  end function nonfact_2loop_VVtoHH_tensor

!----------------------------------------------------------------------
  ! dot product
  real(dp) function dot(p1,p2)
    real(dp), intent(in) :: p1(0:3), p2(0:3)
    dot = p1(0)*p2(0) - sum(p1(1:3)*p2(1:3))
  end function dot
  
  ! !----------------------------------------------------------------------
  ! ! mu_R with scale_choice < 2
  ! real(dp) function muR(Q)
  !   real(dp), intent(in) :: Q
  !   muR = zero
  !   if (scale_choice.eq.0) then
  !      ! scale_choice = 0 : mu_R = xmur * mh
  !      muR = xmur * mh
  !   elseif (scale_choice.eq.-1) then
  !      ! scale_choice = -1 : mu_R = xmur * mw
  !      muR = xmur * mw
  !   elseif (scale_choice.eq.1) then
  !      ! scale_choice = 1 : mu_R = xmur * Q
  !      muR = xmur * Q
  !   endif
  ! end function muR
  

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
       ! else if scale_choice=3, use mixedscale
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

  ! !----------------------------------------------------------------------
  ! ! mu_F with scale_choice < 2
  ! real(dp) function muF(Q)
  !   real(dp), intent(in) :: Q
  !   muF = zero
  !   if (scale_choice.eq.0) then
  !      ! scale_choice = 0 : muF = xmuf * mh
  !      muF = xmuf * mh
  !   elseif (scale_choice.eq.-1) then
  !      ! scale_choice = -1 : mu_F = xmuf * mw
  !      muF = xmuf * mw
  !   elseif (scale_choice.eq.1) then
  !      ! scale_choice = 1 : muF = xmuf * Q
  !      muF = xmuf * Q
  !   else
  !      call wae_error('muF(Q)', 'illegal value for scale_choice', intval = scale_choice)
  !   endif
  ! end function muF
  
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

  ! NOT FUNCTIONAL! VBFNLO TENSORS USE DIFFERENT SPINOR CONVENTION
  !----------------------------------------------------------------------
  ! This function returns the same output as eval_matrix_element, but
  ! is done with numerical tensor manipulations. It is significantly
  ! slower than the analytically computed routine above.
  function eval_matrix_element_vbfnlo(order_start,order_stop, x1, x2, P1, P2, q1, q2, &
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
    real(dp) :: sigma, sigmaWpWm, sigmaWmWp, sigmaZZ, pH1H2(0:3,3:4)
    complex(dp) :: AA(0:3,0:3)

    stop 'eval_matrix_element_vbfnlo does not work'
    ! Start by initialising all the tensors needed
    if(.not.gmunu%initialised) then
       call init_couplings
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

    pH1H2(:,3) = pH1(:) 
    pH1H2(:,4) = pH2(:)

!    print*, 'Q1', q1
!    print*, 'Q2', q2
!    print*, 'pH1', pH1H2(:,3)
    !    print*, 'pH2', pH1H2(:,4)
    
    call WWtoHH(-Q1,-Q2,pH1H2,AA)
!    print*, AA
!    stop
    M(Wp) = g_mu_nu ! Initialise
    M(Wp)%values(0:3,0:3) = AA(0:3,0:3) 

!    print*, 'VBFNLO W tensor'
!    call PrintTensor(M(Wp))
    
    call ZZtoHH(-Q1,-Q2,pH1H2,AA)
    M(Z) = g_mu_nu ! Initialise
    M(Z)%values(0:3,0:3) = AA(0:3,0:3) 

!    print*, 'VBFNLO Z tensor'
!    call PrintTensor(M(Z))
    
!    M(Wp) = (two*((two*MWsq)/cmplx(q1k1sq-MWsq,MW*W_WIDTH) &
!         & + (two*MWsq)/cmplx(q1k2sq-MWsq,MW*W_WIDTH) & 
!         & + three*MH**2/cmplx(k1k2sq - MH**2,MH*HWIDTH) + one))*g_mu_nu ! gmunu
!                                                                               ! contribution
!    M(Wp) = M(Wp) + (one/cmplx(q1k1sq -MWsq,MW*W_WIDTH)) &
!         &           *MWsq/cmplx(MWsq,-MW*W_WIDTH) &
!         &           *((two*k1mu+q1mu).otimes.(k2mu-k1mu-q1mu))
!    M(Wp) = M(Wp) + (one/cmplx(q1k2sq -MWsq,MW*W_WIDTH)) &
!         &           *MWsq/cmplx(MWsq,-MW*W_WIDTH) &
!         &           *((two*k2mu+q1mu).otimes.(k1mu-k2mu-q1mu))
!    M(Wm) = M(Wp) ! Matrix element is identical between W+/W-
!
!    M(Z) = (two*((two*MZsq)/cmplx(q1k1sq-MZsq,MZ*Z_WIDTH) &
!         & + (two*MZsq)/cmplx(q1k2sq-MZsq,MZ*Z_WIDTH) & 
!         & + three*MH**2/cmplx(k1k2sq - MH**2,MH*HWIDTH) + one))*g_mu_nu ! gmunu
!                                                                               ! contribution
!    M(Z) = M(Z) + (one/cmplx(q1k1sq -MZsq,MZ*Z_WIDTH)) &
!         &         *MZsq/cmplx(MZsq,-MZ*Z_WIDTH) &
!         &         *((two*k1mu+q1mu).otimes.(k2mu-k1mu-q1mu))
!    M(Z) = M(Z) + (one/cmplx(q1k2sq -MZsq,MZ*Z_WIDTH)) &
!         &         *MZsq/cmplx(MZsq,-MZ*Z_WIDTH) &
!         &         *((two*k2mu+q1mu).otimes.(k1mu-k2mu-q1mu))

! OLD - no width
!    M(Wp) = (two*((two*MWsq)/(q1k1sq-MWsq) &
!         & + (two*MWsq)/(q1k2sq-MWsq) & 
!         & + three*MH**2/(k1k2sq - MH**2) + one))*g_mu_nu ! gmunu
!                                                                               ! contribution
!    M(Wp) = M(Wp) + (one/(q1k1sq -MWsq)) &
!         &           *((two*k1mu+q1mu).otimes.(k2mu-k1mu-q1mu))
!    M(Wp) = M(Wp) + (one/(q1k2sq -MWsq)) &
!         &           *((two*k2mu+q1mu).otimes.(k1mu-k2mu-q1mu))
!    M(Wm) = M(Wp) ! Matrix element is identical between W+/W-
!!
!    M(Z) = (two*((two*MZsq)/(q1k1sq-MZsq) &
!         & + (two*MZsq)/(q1k2sq-MZsq) & 
!         & + three*MH**2/(k1k2sq - MH**2) + one))*g_mu_nu ! gmunu
!                                                                               ! contribution
!    M(Z) = M(Z) + (one/(q1k1sq -MZsq)) &
!         &         *((two*k1mu+q1mu).otimes.(k2mu-k1mu-q1mu))
!    M(Z) = M(Z) + (one/(q1k2sq -MZsq)) &
!         &         *((two*k2mu+q1mu).otimes.(k1mu-k2mu-q1mu))
!
!    overall_norm = (GFermi**4)/(S) * two
!    WW_norm =  MW**8 / (((Q1sq + MWsq)**2 + W_WIDTH**2 * MWsq)& 
!         & * ((Q2sq + MWsq)**2 + W_WIDTH**2 * MWsq ))
!    ZZ_norm =  MZ**8 / (((Q1sq + MZsq)**2 + Z_WIDTH**2 * MZsq)& 
!         & * ((Q2sq + MZsq)**2 + Z_WIDTH**2 * MZsq ))
!!
!    M(Wp)%values = M(Wp)%values * sqrt(overall_norm * WW_norm)
!    M(Z)%values = M(Z)%values * sqrt(overall_norm * ZZ_norm)
!    print*, 'proVBFHH W tensor'
!    call PrintTensor(M(Wp))
!    print*, 'proVBFHH Z tensor'
!    call PrintTensor(M(Z))
!    stop
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
    T3(1)%values(:,:) = cmplx(zero,one,kind=dp) * T3(1)%values(:,:)
    T3(2)%values(:,:) = cmplx(zero,one,kind=dp) * T3(2)%values(:,:)
    
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
!    overall_norm = one!(GFermi**4)/(S) * two
    overall_norm = (GFermi**2)/(S) !* two
    WW_norm =  MW**4! / (((Q1sq + MWsq)**2 + W_WIDTH**2 * MWsq)& 
         !& * ((Q2sq + MWsq)**2 + W_WIDTH**2 * MWsq ))
    ZZ_norm =  MZ**4! / (((Q1sq + MZsq)**2 + Z_WIDTH**2 * MZsq)& 
         !& * ((Q2sq + MZsq)**2 + Z_WIDTH**2 * MZsq ))

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

  end function eval_matrix_element_vbfnlo
end module matrix_element
