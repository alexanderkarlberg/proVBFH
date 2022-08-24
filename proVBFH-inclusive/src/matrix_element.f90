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
module matrix_element
  use hoppet_v1
  use parameters
  implicit none

  private
  public :: eval_matrix_element
  interface operator(.dot.)
     module procedure dot
  end interface operator(.dot.)

contains

  !----------------------------------------------------------------------
  function eval_matrix_element(order_start,order_stop, x1, x2, P1, P2, q1, q2, ptH) result(res)
    integer , intent(in) :: order_start,order_stop
    real(dp), intent(in) :: x1, x2, P1(0:3), P2(0:3), q1(0:3), q2(0:3), ptH
    real(dp)             :: res
    !----------------------------------------------------------------------
    real(dp) :: y1, y2, Q1sq, Q2sq, Q1val, Q2val
    real(dp) :: muR1val, muR2val, muF1val, muF2val
    real(dp) :: Fx1(-6:7,4), Fx2(-6:7,4)
    real(dp) :: F1F1, F1F2, F2F1, F2F2, F3F3
    real(dp) :: WW_norm, ZZ_norm, overall_norm
    real(dp) :: q1q2, P1q1, P1q2, P2q1, P2q2, P1P2
    integer  :: i, j
    logical, parameter :: WpWm = .true., WmWp = .true., ZZ = .true.
    integer :: iorder

    y1 = -log(x1)
    y2 = -log(x2)

    Q1sq = -(q1.dot.q1)
    Q2sq = -(q2.dot.q2)
    Q1val = sqrt(Q1sq)
    Q2val = sqrt(Q2sq)
    
    muR1val = muR1(Q1val, Q2val, ptH)
    muR2val = muR2(Q1val, Q2val, ptH)
    muF1val = muF1(Q1val, Q2val, ptH)
    muF2val = muF2(Q1val, Q2val, ptH)

    F1F1 = zero
    F2F1 = zero
    F1F2 = zero
    F2F2 = zero
    F3F3 = zero

    overall_norm = (GFermi**3)/(S) * four*sqrt(two)
    WW_norm =  MW**8 / (((Q1sq + MW**2)**2 + W_WIDTH**2 * MW**2)& 
         & * ((Q2sq + MW**2)**2 + W_WIDTH**2 * MW**2 ))
    ZZ_norm =  MZ**8 / (((Q1sq + MZ**2)**2 + Z_WIDTH**2 * MZ**2)& 
         & * ((Q2sq + MZ**2)**2 + Z_WIDTH**2 * MZ**2 ))
    ! Below is a narrow-width propagator; This was used in earlier versions of the code
    ! ZZ_norm =  MZ**8 / ((Q1sq + MZ**2)**2 * (Q2sq + MZ**2)**2)
    ! WW_norm =  MW**8 / ((Q1sq + MW**2)**2 * (Q2sq + MW**2)**2)

    ! Compute the LO structure funtion by adding all the pieces
    ! from tables
    Fx1(:,1) = F_LO(y1, Q1val, muR1val, muF1val)
    Fx2(:,1) = F_LO(y2, Q2val, muR2val, muF2val)

    if (order_stop.ge.2) then
       ! Compute the NLO structure funtion by adding all the pieces
       ! from tables
       Fx1(:,2) = F_NLO(y1, Q1val, muR1val, muF1val)
       Fx2(:,2) = F_NLO(y2, Q2val, muR2val, muF2val)
    endif

    if (order_stop.ge.3) then
       ! Compute the NNLO structure funtion by adding all the pieces
       ! from tables
       Fx1(:,3) = F_NNLO(y1, Q1val, muR1val, muF1val)
       Fx2(:,3) = F_NNLO(y2, Q2val, muR2val, muF2val)
    endif

    if (order_stop.ge.4) then
       ! Compute the N3LO structure funtion by adding all the pieces
       ! from tables
       Fx1(:,4) = F_N3LO(y1, Q1val, muR1val, muF1val)
       Fx2(:,4) = F_N3LO(y2, Q2val, muR2val, muF2val)
    endif

    do iorder = order_start,order_stop
       do i = 1, iorder
          j = 1 + iorder - i
          if (WpWm) then
             F1F1 = F1F1 + WW_norm * Fx1(F1Wp,i)*Fx2(F1Wm,j)
             F2F1 = F2F1 + WW_norm * Fx1(F2Wp,i)*Fx2(F1Wm,j)
             F1F2 = F1F2 + WW_norm * Fx1(F1Wp,i)*Fx2(F2Wm,j)
             F2F2 = F2F2 + WW_norm * Fx1(F2Wp,i)*Fx2(F2Wm,j)
             F3F3 = F3F3 + WW_norm * Fx1(F3Wp,i)*Fx2(F3Wm,j)
          end if
          
          if (WmWp) then
             F1F1 = F1F1 + WW_norm * Fx1(F1Wm,i)*Fx2(F1Wp,j)
             F2F1 = F2F1 + WW_norm * Fx1(F2Wm,i)*Fx2(F1Wp,j)
             F1F2 = F1F2 + WW_norm * Fx1(F1Wm,i)*Fx2(F2Wp,j)
             F2F2 = F2F2 + WW_norm * Fx1(F2Wm,i)*Fx2(F2Wp,j)
             F3F3 = F3F3 + WW_norm * Fx1(F3Wm,i)*Fx2(F3Wp,j)
          end if
          
          if (ZZ) then
             F1F1 = F1F1 + ZZ_norm * Fx1(F1Z,i)*Fx2(F1Z,j)
             F2F1 = F2F1 + ZZ_norm * Fx1(F2Z,i)*Fx2(F1Z,j)
             F1F2 = F1F2 + ZZ_norm * Fx1(F1Z,i)*Fx2(F2Z,j)
             F2F2 = F2F2 + ZZ_norm * Fx1(F2Z,i)*Fx2(F2Z,j)
             F3F3 = F3F3 + ZZ_norm * Fx1(F3Z,i)*Fx2(F3Z,j)
          end if
       end do
    end do

    q1q2 = q1 .dot. q2
    P1q1 = P1 .dot. q1
    P1q2 = P1 .dot. q2
    P2q1 = P2 .dot. q1
    P2q2 = P2 .dot. q2
    P1P2 = P1 .dot. P2

    res = zero
    res = res + F1F1*(two + q1q2**2/(Q1sq*Q2sq))
    res = res + F1F2/P2q2 * (P2q2**2/(-Q2sq) + (P2q1 - P2q2*q1q2/(-Q2sq))**2/(-Q1sq))
    res = res + F2F1/P1q1 * (P1q1**2/(-Q1sq) + (P1q2 - P1q1*q1q2/(-Q1sq))**2/(-Q2sq))
    res = res + F2F2/(P1q1*P2q2)*(P1P2 - P1q1*P2q1/(-Q1sq) - P1q2*P2q2/(-Q2sq) + P1q1*P2q2*q1q2/(Q1sq*Q2sq))**2
    res = res + F3F3/(two*P1q1*P2q2)*(P1P2*q1q2 - P1q2*P2q1)
    
    res = res * overall_norm

  end function eval_matrix_element

  !----------------------------------------------------------------------
  ! dot product
  real(dp) function dot(p1,p2)
    real(dp), intent(in) :: p1(0:3), p2(0:3)
    dot = p1(0)*p2(0) - sum(p1(1:3)*p2(1:3))
  end function dot

  !----------------------------------------------------------------------
  ! mu_R1 as a function of Q1 and Q2
  real(dp) function muR1(Q1, Q2, ptH)
    real(dp), intent(in) :: Q1, Q2, ptH
    muR1 = zero
    if (scale_choice.le.1) then
       ! if scale_choice = 0,1 then muR1(Q1,Q2) = muR(Q1)
       muR1 = muR(Q1)
    elseif (scale_choice.eq.2) then
       ! else if scale_choice=2, use sqrt(Q1*Q2)
       muR1 = xmur * sqrt(Q1 * Q2)
    elseif (scale_choice.eq.3) then
       ! else if scale_choice=3, use mixed scale
       muR1 = xmur * mixed_scale(Q1, Q2, ptH)
    endif
  end function muR1
  
  !----------------------------------------------------------------------
  ! mu_R2 as a function of Q1 and Q2
  real(dp) function muR2(Q1, Q2, ptH)
    real(dp), intent(in) :: Q1, Q2, ptH
    muR2 = zero
    if (scale_choice.le.1) then
       ! if scale_choice = 0,1 then muR1(Q1,Q2) = muR(Q1)
       muR2 = muR(Q2)
    elseif (scale_choice.eq.2) then
       ! else if scale_choice=2, use sqrt(Q1*Q2)
       muR2 = xmur * sqrt(Q1 * Q2)
    elseif (scale_choice.eq.3) then
       ! else if scale_choice=3, use mixed scale
       muR2 = xmur * mixed_scale(Q1, Q2, ptH)
    endif
  end function muR2

  !----------------------------------------------------------------------
  ! mu_F1 as a function of Q1 and Q2
  real(dp) function muF1(Q1, Q2, ptH)
    real(dp), intent(in) :: Q1, Q2, ptH
    muF1 = zero
    if (scale_choice.le.1) then
       ! if scale_choice = 0,1 then muF1(Q1,Q2) = muF(Q1)
       muF1 = muF(Q1)
    elseif (scale_choice.eq.2) then
       ! else if scale_choice=2, use sqrt(Q1*Q2)
       muF1 = xmuf * sqrt(Q1 * Q2)
    elseif (scale_choice.eq.3) then
       ! else if scale_choice=3, use mixed scale
       muF1 = xmuf * mixed_scale(Q1, Q2, ptH)
    else
       call wae_error('muF1(Q)', 'illegal value for scale_choice', intval = scale_choice)
    endif
  end function muF1

  !----------------------------------------------------------------------
  ! mu_F2 as a function of Q1 and Q2
  real(dp) function muF2(Q1, Q2, ptH)
    real(dp), intent(in) :: Q1, Q2, ptH
    muF2 = zero
    if (scale_choice.le.1) then
       ! if scale_choice = 0,1 then muF1(Q1,Q2) = muF(Q1)
       muF2 = muF(Q2)
    elseif (scale_choice.eq.2) then
       ! else if scale_choice=2, use sqrt(Q1*Q2)
       muF2 = xmuf * sqrt(Q1 * Q2)
    elseif (scale_choice.eq.3) then
       ! else if scale_choice=3, use mixed scale
       muF2 = xmuf * mixed_scale(Q1, Q2, ptH)
    else
       call wae_error('muF2(Q)', 'illegal value for scale_choice', intval = scale_choice)
    endif
  end function muF2

  !----------------------------------------------------------------------
  ! Defines which scale to use for scale_choice = 3,
  ! which can be any function of Q1, Q2 and ptH
  real(dp) function mixed_scale(Q1,Q2,ptH)
    real(dp), intent(in) :: Q1, Q2, ptH
    mixed_scale = ((mh*0.5d0)**4d0+(mh*ptH*0.5d0)**2d0)**(0.25d0)
  end function mixed_scale

end module matrix_element
