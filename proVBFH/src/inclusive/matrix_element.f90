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
  use hoppet_v1
  use incl_parameters
  use streamlined_interface
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
  function eval_matrix_element(order_start,order_stop, x1, x2, P1, P2, q1, q2, ptH) result(res)
    use nonfact
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
    real(dp) :: alphas_lcl,ptj1,ptj2,muR_nonfact
    
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
    if(non_fact) then
       if(order_stop.gt.1) stop 'Cannot do non factorisable corrections'
       !       ptj1 = q1(1)**2 + q1(2)**2
       !       ptj2 = q2(1)**2 + q2(2)**2
       !       muR_nonfact = max((ptj1*ptj2)**0.25_dp,Qmin)
       !       muR_nonfact = muR_nonfact * xmur
       muR_nonfact = sqrt(muR1val * muR2val)
       alphas_lcl = Value(coupling,muR_nonfact)
       ! Eq. 6 in 1906.10899 for Nc = 3
       WW_norm = WW_norm * (2.0_dp/9.0_dp * alphas_lcl**2 * chinf(q1(1:2),q2(1:2),MW))
       ZZ_norm = ZZ_norm * (2.0_dp/9.0_dp * alphas_lcl**2 * chinf(q1(1:2),q2(1:2),MZ))
    endif
    ! Below is a narrow-width propagator; This was used in earlier versions of the code
    ! ZZ_norm =  MZ**8 / ((Q1sq + MZ**2)**2 * (Q2sq + MZ**2)**2)
    ! WW_norm =  MW**8 / ((Q1sq + MW**2)**2 * (Q2sq + MW**2)**2)

    ! Compute the LO structure funtion by adding all the pieces
    ! from tables
    Fx1(:,1) = two*F_LO(y1, Q1val, muR1val, muF1val)
    Fx2(:,1) = two*F_LO(y2, Q2val, muR2val, muF2val)

    if (order_stop.ge.2) then
       ! Compute the NLO structure funtion by adding all the pieces
       ! from tables
       Fx1(:,2) = two*F_NLO(y1, Q1val, muR1val, muF1val)
       Fx2(:,2) = two*F_NLO(y2, Q2val, muR2val, muF2val)
    endif

    if (order_stop.ge.3) then
       ! Compute the NNLO structure funtion by adding all the pieces
       ! from tables
       Fx1(:,3) = two*F_NNLO(y1, Q1val, muR1val, muF1val)
       Fx2(:,3) = two*F_NNLO(y2, Q2val, muR2val, muF2val)
    endif

    if (order_stop.ge.4) then
       ! Compute the N3LO structure funtion by adding all the pieces
       ! from tables
       Fx1(:,4) = two*F_N3LO(y1, Q1val, muR1val, muF1val)
       Fx2(:,4) = two*F_N3LO(y2, Q2val, muR2val, muF2val)
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
  ! This function returns the same output as eval_matrix_element, but
  ! is done with numerical tensor manipulations. It is significantly
  ! slower than the analytically computed routine above.
  function eval_matrix_element_tensor(order_start,order_stop, x1, x2, P1, P2, q1, q2, ptH) result(res)
    integer , intent(in) :: order_start,order_stop
    real(dp), intent(in) :: x1, x2, P1(0:3), P2(0:3), q1(0:3), q2(0:3), ptH
    real(dp)             :: res
    !----------------------------------------------------------------------
    real(dp) :: y1, y2, Q1sq, Q2sq, Q1val, Q2val
    real(dp) :: muR1val, muR2val, muF1val, muF2val
    real(dp) :: Fx1(-6:7,4), Fx2(-6:7,4), F1(-6:7), F2(-6:7)
    real(dp) :: WW_norm, ZZ_norm, overall_norm
    real(dp) :: q1q2, P1q1, P1q2, P2q1, P2q2, P1P2,q1q1,q2q2
    integer  :: i, j
    logical, parameter :: WpWm = .true., WmWp = .true., ZZ = .true.
    integer, parameter :: Wp=1, Wm=-1, Z=0
    integer :: iorder
    type(tensors) :: M,Mstar, Wx1(-1:1), Wx2(-1:1), g_mu_nu, q1mu, q2mu, P1mu, P2mu
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
    P1mu%values(:,1) = P1(:)
    P2mu%values(:,1) = P2(:)

    !Compute all combinations of dot-products
    q1q1 = q1 .dot. q1
    q2q2 = q2 .dot. q2
    q1q2 = q1 .dot. q2
    P1q1 = P1 .dot. q1
    P1q2 = P1 .dot. q2
    P2q1 = P2 .dot. q1
    P2q2 = P2 .dot. q2
    P1P2 = P1 .dot. P2

    P1hatmu = P1mu - P1q1/q1q1 * q1mu 
    P2hatmu = P2mu - P2q2/q2q2 * q2mu

    M = gmunu ! Trivial for now
    Mstar = M ! For initialisation
    Mstar%values = dconjg(M%values) ! Complex conjugate
    

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
    T3(1)%values(1,2) =   P1mu%values(3,1)*q1mu%values(0,1) &
         &              - P1mu%values(0,1)*q1mu%values(3,1) 
    T3(1)%values(1,3) = - P1mu%values(2,1)*q1mu%values(0,1) &
         &              + P1mu%values(0,1)*q1mu%values(2,1) 
    T3(1)%values(2,3) =   P1mu%values(1,1)*q1mu%values(0,1) &
         &              - P1mu%values(0,1)*q1mu%values(1,1) 

    T3(2)%values(0,1) =   P2mu%values(2,1)*q2mu%values(3,1) &
         &              - P2mu%values(3,1)*q2mu%values(2,1)
    T3(2)%values(0,2) = - P2mu%values(1,1)*q2mu%values(3,1) &
         &              + P2mu%values(3,1)*q2mu%values(1,1) 
    T3(2)%values(0,3) =   P2mu%values(1,1)*q2mu%values(2,1) &
         &              - P2mu%values(2,1)*q2mu%values(1,1) 
    T3(2)%values(1,2) =   P2mu%values(3,1)*q2mu%values(0,1) &
         &              - P2mu%values(0,1)*q2mu%values(3,1) 
    T3(2)%values(1,3) = - P2mu%values(2,1)*q2mu%values(0,1) &
         &              + P2mu%values(0,1)*q2mu%values(2,1) 
    T3(2)%values(2,3) =   P2mu%values(1,1)*q2mu%values(0,1) &
         &              - P2mu%values(0,1)*q2mu%values(1,1) 

    ! Use anti-symmetric property    
    do i = 0,3
       do j = i,3
          T3(1)%values(j,i) = - T3(1)%values(i,j)
          T3(2)%values(j,i) = - T3(2)%values(i,j)
       enddo
    enddo

    ! And now mulitply with overall factor i
    T3(1)%values(:,:) = DCMPLX(zero,one) * T3(1)%values(:,:)
    T3(2)%values(:,:) = DCMPLX(zero,one) * T3(2)%values(:,:)
    
    y1 = -log(x1)
    y2 = -log(x2)

    Q1sq = -q1q1
    Q2sq = -q2q2
    Q1val = sqrt(Q1sq)
    Q2val = sqrt(Q2sq)
    
    muR1val = muR1(Q1val, Q2val, ptH)
    muR2val = muR2(Q1val, Q2val, ptH)
    muF1val = muF1(Q1val, Q2val, ptH)
    muF2val = muF2(Q1val, Q2val, ptH)

    !These will contain the full structure functions for the two protons
    F1 = zero
    F2 = zero
    
    ! Compute the overall numerical factors
    overall_norm = (GFermi**3)/(S) * four*sqrt(two)
    WW_norm =  MW**8 / (((Q1sq + MW**2)**2 + W_WIDTH**2 * MW**2)& 
         & * ((Q2sq + MW**2)**2 + W_WIDTH**2 * MW**2 ))
    ZZ_norm =  MZ**8 / (((Q1sq + MZ**2)**2 + Z_WIDTH**2 * MZ**2)& 
         & * ((Q2sq + MZ**2)**2 + Z_WIDTH**2 * MZ**2 ))

    ! Compute the LO structure funtion by adding all the pieces
    ! from tables
    Fx1(:,1) = two*F_LO(y1, Q1val, muR1val, muF1val)
    Fx2(:,1) = two*F_LO(y2, Q2val, muR2val, muF2val)
    F1(:) = Fx1(:,1)
    F2(:) = Fx2(:,1)

    if (order_stop.ge.2) then
       ! Compute the NLO structure funtion by adding all the pieces
       ! from tables
       Fx1(:,2) = two*F_NLO(y1, Q1val, muR1val, muF1val)
       Fx2(:,2) = two*F_NLO(y2, Q2val, muR2val, muF2val)
       F1(:) = F1(:) + Fx1(:,2)
       F2(:) = F2(:) + Fx2(:,2)
    endif

    if (order_stop.ge.3) then
       ! Compute the NNLO structure funtion by adding all the pieces
       ! from tables
       Fx1(:,3) = two*F_NNLO(y1, Q1val, muR1val, muF1val)
       Fx2(:,3) = two*F_NNLO(y2, Q2val, muR2val, muF2val)
       F1(:) = F1(:) + Fx1(:,3)
       F2(:) = F2(:) + Fx2(:,3)
    endif

    if (order_stop.ge.4) then
       ! Compute the N3LO structure funtion by adding all the pieces
       ! from tables
       Fx1(:,4) = two*F_N3LO(y1, Q1val, muR1val, muF1val)
       Fx2(:,4) = two*F_N3LO(y2, Q2val, muR2val, muF2val)
       F1(:) = F1(:) + Fx1(:,4)
       F2(:) = F2(:) + Fx2(:,4)
    endif

    ! Compute hadronic tensors
    Wx1(Wp) = F1(F1Wp)*((one/q1q1)*(q1mu.otimes.q1mu)-g_mu_nu) &
         & + F1(F2Wp)*(one/P1q1)*(P1hatmu.otimes.P1hatmu) &
         & + F1(F3Wp)*(one/(two*P1q1))*T3(1)

    Wx1(Wm) = F1(F1Wm)*((one/q1q1)*(q1mu.otimes.q1mu)-g_mu_nu) &
         & + F1(F2Wm)*(one/P1q1)*(P1hatmu.otimes.P1hatmu) &
         & + F1(F3Wm)*(one/(two*P1q1))*T3(1)

    Wx1(Z) = F1(F1Z)*((one/q1q1)*(q1mu.otimes.q1mu)-g_mu_nu) &
         & + F1(F2Z)*(one/P1q1)*(P1hatmu.otimes.P1hatmu) &
         & + F1(F3Z)*(one/(two*P1q1))*T3(1)

    Wx2(Wp) = F2(F1Wp)*((one/q2q2)*(q2mu.otimes.q2mu)-g_mu_nu) &
         & + F2(F2Wp)*(one/P2q2)*(P2hatmu.otimes.P2hatmu) &
         & + F2(F3Wp)*(one/(two*P2q2))*T3(2)

    Wx2(Wm) = F2(F1Wm)*((one/q2q2)*(q2mu.otimes.q2mu)-g_mu_nu) &
         & + F2(F2Wm)*(one/P2q2)*(P2hatmu.otimes.P2hatmu) &
         & + F2(F3Wm)*(one/(two*P2q2))*T3(2)

    Wx2(Z) = F2(F1Z)*((one/q2q2)*(q2mu.otimes.q2mu)-g_mu_nu) &
         & + F2(F2Z)*(one/P2q2)*(P2hatmu.otimes.P2hatmu) &
         & + F2(F3Z)*(one/(two*P2q2))*T3(2)

    ! Do first contraction between a hadronic tensor and a matrix element.
    call ContractTensors(Wx1(Wp),1,M,1,MWx1(Wp))
    call ContractTensors(Wx1(Wm),1,M,1,MWx1(Wm))
    call ContractTensors(Wx1(Z),1,M,1,MWx1(Z))

    call ContractTensors(Mstar,2,Wx2(Wp),2,MstarWx2(Wp))
    call ContractTensors(Mstar,2,Wx2(Wm),2,MstarWx2(Wm))
    call ContractTensors(Mstar,2,Wx2(Z),2,MstarWx2(Z))

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
