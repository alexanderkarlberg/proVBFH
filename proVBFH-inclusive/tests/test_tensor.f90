! Unit test for src/tensor.f90.
!
! Compares every operation of the tensor module with the previous,
! allocatable implementation (tests/tensor_legacy.f90) and with
! explicit component formulae, on random inputs. It then evaluates
! the chain of operations used by the VBF HH matrix elements
! (current M, hadronic tensors W1 and W2, Tr[(W1 M) (M* W2)]) with both
! implementations, and times it.
!
! Build and run from proVBFH-inclusive/ with "make check".
!
! The two implementations perform the same floating-point operations
! in the same order, so they are expected to agree exactly; the test
! fails if any relative difference exceeds tol. The one intended
! difference is the rank-2 x rank-1 contraction over the second index
! of the rank-2 tensor, which the legacy code got wrong (it is not used
! by any of the matrix elements); there the new code is checked
! against the explicit formula only.
program test_tensor
  use hoppet, only: dp, zero, one, two
  use tensor
  use tensor_legacy, only: ltensors => tensors, lInitTensor => InitTensor, &
       & lSetMetric => SetMetric, lContractTensors => ContractTensors, &
       & lraise => raise, llower => lower, lTensorTrace => TensorTrace, &
       & lgmunu => gmunu, operator(+), operator(-), operator(*), operator(.otimes.)
  implicit none

  real(dp), parameter :: tol = 1e-15_dp
  integer,  parameter :: ntrials = 10000, nchain = 1000, ntime = 200000
  real(dp) :: maxdev_legacy, maxdev_formula, maxdev_chain
  integer  :: nfail, itrial

  nfail = 0
  maxdev_legacy = zero
  maxdev_formula = zero
  maxdev_chain = zero

  call init_random()
  call SetMetric(1)
  call lSetMetric(1)
  if (reldev(gmunu%values, lgmunu%values) > zero) call fail('SetMetric')

  do itrial = 1, ntrials
     call test_contractions()
     call test_raise_lower()
     call test_arithmetic()
     call test_products()
     call test_four_vector()
  end do

  do itrial = 1, nchain
     call test_chain()
  end do

  write(*,'(a,es10.3)') ' max rel. deviation from legacy code:      ', maxdev_legacy
  write(*,'(a,es10.3)') ' max rel. deviation from explicit formulae:', maxdev_formula
  write(*,'(a,es10.3)') ' max rel. deviation, matrix-element chain: ', maxdev_chain

  call time_chain()

  if (nfail > 0) then
     write(*,'(a,i0,a)') ' FAILED: ', nfail, ' check(s) above tolerance'
     error stop 1
  end if
  write(*,'(a)') ' All tensor tests passed.'

contains

  !----------------------------------------------------------------------
  ! Contractions of all rank combinations and index choices
  subroutine test_contractions()
    type(tensors)  :: a, b, c
    type(ltensors) :: la, lb, lc
    complex(dp)    :: ref(0:3,0:3)
    integer :: r1, r2, i1, i2, i, j, k
    logical :: up1(2), up2(2)

    do r1 = 1, 2
       do r2 = 1, 2
          do i1 = 1, r1
             do i2 = 1, r2
                call random_up(up1)
                call random_up(up2)
                up2(i2) = .not. up1(i1)
                call random_pair(a, la, r1, up1)
                call random_pair(b, lb, r2, up2)
                call ContractTensors(a, i1, b, i2, c)
                ! The legacy code does not resize an output tensor that
                ! is already initialised (and then writes out of
                ! bounds), so give it one of the right rank.
                call lInitTensor(lc, r1+r2-2, .true.)
                call lContractTensors(la, i1, lb, i2, lc)

                ! Explicit formula
                ref = zero
                if (r1 == 1 .and. r2 == 1) then
                   do i = 0, 3
                      ref(1,1) = ref(1,1) + a%values(i,1)*b%values(i,1)
                   end do
                else if (r1 == 1 .and. r2 == 2) then
                   do j = 0, 3
                      do i = 0, 3
                         if (i2 == 1) ref(j,1) = ref(j,1) + a%values(i,1)*b%values(i,j)
                         if (i2 == 2) ref(j,1) = ref(j,1) + a%values(i,1)*b%values(j,i)
                      end do
                   end do
                else if (r1 == 2 .and. r2 == 1) then
                   do j = 0, 3
                      do i = 0, 3
                         if (i1 == 1) ref(j,1) = ref(j,1) + a%values(i,j)*b%values(i,1)
                         if (i1 == 2) ref(j,1) = ref(j,1) + a%values(j,i)*b%values(i,1)
                      end do
                   end do
                else
                   do k = 0, 3
                      do j = 0, 3
                         do i = 0, 3
                            ref(j,k) = ref(j,k) + elem(a, i1, i, j) * elem(b, i2, i, k)
                         end do
                      end do
                   end do
                end if
                call check_formula(c, ref, r1+r2-2, 'ContractTensors')
                call check_up(c, contracted_up(up1, r1, i1, up2, r2, i2), 'ContractTensors')

                if (.not. (r1 == 2 .and. r2 == 1 .and. i1 == 2)) then
                   call check_legacy(c, lc, 'ContractTensors')
                end if
             end do
          end do
       end do
    end do
  end subroutine test_contractions

  !----------------------------------------------------------------------
  subroutine test_raise_lower()
    type(tensors)  :: a
    type(ltensors) :: la
    complex(dp)    :: ref(0:3,0:3)
    real(dp)       :: g(0:3)
    integer :: r, idx
    logical :: upv(2)

    g = (/ one, -one, -one, -one /)
    do r = 1, 2
       do idx = 1, r
          call random_up(upv)
          call random_pair(a, la, r, upv)
          ref = a%values
          if (upv(idx) .eqv. .false.) then
             call raise(a, idx)
             call lraise(la, idx)
          else
             call lower(a, idx)
             call llower(la, idx)
          end if
          if (idx == 1) ref = spread(g, 2, 4) * ref
          if (idx == 2) ref = spread(g, 1, 4) * ref
          upv(idx) = .not. upv(idx)
          call check_formula(a, ref, r, 'raise/lower')
          call check_up(a, upv, 'raise/lower')
          call check_legacy(a, la, 'raise/lower')
          ! A second call must not change anything
          if (upv(idx)) then
             call raise(a, idx)
          else
             call lower(a, idx)
          end if
          call check_formula(a, ref, r, 'raise/lower (no-op)')
       end do
    end do
  end subroutine test_raise_lower

  !----------------------------------------------------------------------
  subroutine test_arithmetic()
    type(tensors)  :: a, b, c
    type(ltensors) :: la, lb, lc
    real(dp)       :: x, y
    complex(dp)    :: z
    integer :: r
    logical :: upv(2)

    do r = 0, 2
       call random_up(upv)
       call random_pair(a, la, r, upv)
       call random_pair(b, lb, r, upv)
       call random_number(x); call random_number(y)
       x = x - 0.5_dp; z = cmplx(x, y - 0.5_dp, kind=dp)

       c = a + b; lc = la + lb
       call check_formula(c, a%values + b%values, r, 'AddTensors')
       call check_legacy(c, lc, 'AddTensors')

       c = a - b; lc = la - lb
       call check_formula(c, a%values - b%values, r, 'SubtractTensors')
       call check_legacy(c, lc, 'SubtractTensors')

       c = x * a; lc = x * la
       call check_formula(c, a%values * x, r, 'MultiplyWithScalar')
       call check_legacy(c, lc, 'MultiplyWithScalar')

       c = z * a; lc = z * la
       call check_formula(c, a%values * z, r, 'MultiplyWithComplexScalar')
       call check_legacy(c, lc, 'MultiplyWithComplexScalar')

       if (r == 2) then
          if (reldev_scalar(TensorTrace(a), lTensorTrace(la)) > tol) call fail('TensorTrace')
          maxdev_legacy = max(maxdev_legacy, reldev_scalar(TensorTrace(a), lTensorTrace(la)))
       end if
    end do
  end subroutine test_arithmetic

  !----------------------------------------------------------------------
  subroutine test_products()
    type(tensors)  :: a, b, c
    type(ltensors) :: la, lb, lc
    complex(dp)    :: ref(0:3,0:3)
    integer :: i, j
    logical :: up1(2), up2(2)

    ! rank 1 x rank 1
    call random_up(up1); call random_up(up2)
    call random_pair(a, la, 1, up1)
    call random_pair(b, lb, 1, up2)
    c = a .otimes. b; lc = la .otimes. lb
    do j = 0, 3
       do i = 0, 3
          ref(i,j) = a%values(i,1) * b%values(j,1)
       end do
    end do
    call check_formula(c, ref, 2, 'TensorProduct 1x1')
    call check_up(c, (/ up1(1), up2(1) /), 'TensorProduct 1x1')
    call check_legacy(c, lc, 'TensorProduct 1x1')

    ! rank 0 x rank 1 and rank 1 x rank 0. The legacy code sets all
    ! indices of the result up here, so only the values are compared.
    call random_pair(a, la, 0, up1)
    call random_pair(b, lb, 1, up2)
    c = a .otimes. b; lc = la .otimes. lb
    ref = zero
    ref(:,1) = a%values(1,1) * b%values(:,1)
    call check_formula(c, ref, 1, 'TensorProduct 0x1')
    call check_up(c, up2, 'TensorProduct 0x1')
    maxdev_legacy = max(maxdev_legacy, reldev(c%values(:,1:1), lc%values))
    if (reldev(c%values(:,1:1), lc%values) > tol) call fail('TensorProduct 0x1 (legacy)')

    c = b .otimes. a; lc = lb .otimes. la
    call check_formula(c, ref, 1, 'TensorProduct 1x0')
    call check_up(c, up2, 'TensorProduct 1x0')
    maxdev_legacy = max(maxdev_legacy, reldev(c%values(:,1:1), lc%values))
    if (reldev(c%values(:,1:1), lc%values) > tol) call fail('TensorProduct 1x0 (legacy)')
  end subroutine test_products

  !----------------------------------------------------------------------
  ! InitFourVector: the index positions, the lowered components, and
  ! that contracting p_mu with q^mu gives the Minkowski product
  subroutine test_four_vector()
    type(tensors) :: pu, pd, qu, c
    real(dp) :: p(0:3), q(0:3)
    complex(dp) :: ref(0:3,0:3)

    call random_number(p); call random_number(q)
    p = p - 0.5_dp; q = q - 0.5_dp
    call InitFourVector(pu, p, .true.)
    call InitFourVector(pd, p, .false.)
    call InitFourVector(qu, q, .true.)
    ref = zero
    ref(:,1) = p
    call check_formula(pu, ref, 1, 'InitFourVector (up)')
    call check_up(pu, (/ .true., .true. /), 'InitFourVector (up)')
    ref(1:3,1) = -p(1:3)
    call check_formula(pd, ref, 1, 'InitFourVector (down)')
    call check_up(pd, (/ .false., .true. /), 'InitFourVector (down)')
    ! p.q is summed in a different order than in dot(), and may
    ! involve cancellations, so compare to the scale sum |p_i q_i|
    call ContractTensors(pd, 1, qu, 1, c)
    if (c%rank /= 0) call fail('InitFourVector (p.q) (wrong rank)')
    maxdev_formula = max(maxdev_formula, abs(c%values(1,1) - dot(p, q)) / sum(abs(p*q)))
    if (abs(c%values(1,1) - dot(p, q)) > tol * sum(abs(p*q))) &
         & call fail('InitFourVector (p.q) (differs from the Minkowski product)')
  end subroutine test_four_vector

  !----------------------------------------------------------------------
  ! The chain of operations of eval_matrix_element_tensor, with random
  ! momenta, propagators and structure functions.
  subroutine test_chain()
    real(dp)    :: q1(0:3), q2(0:3), P1(0:3), P2(0:3), k1(0:3), k2(0:3), F1(3), F2(3)
    complex(dp) :: cA, cB, cC
    real(dp)    :: s_new, s_old

    call random_inputs(q1, q2, P1, P2, k1, k2, F1, F2, cA, cB, cC)
    s_new = chain_new(q1, q2, P1, P2, k1, k2, F1, F2, cA, cB, cC)
    s_old = chain_legacy(q1, q2, P1, P2, k1, k2, F1, F2, cA, cB, cC)
    maxdev_chain = max(maxdev_chain, abs(s_new - s_old) / abs(s_old))
    if (abs(s_new - s_old) > tol * abs(s_old)) call fail('matrix-element chain')
  end subroutine test_chain

  subroutine time_chain()
    real(dp)    :: q1(0:3), q2(0:3), P1(0:3), P2(0:3), k1(0:3), k2(0:3), F1(3), F2(3)
    complex(dp) :: cA, cB, cC
    real(dp)    :: acc_new, acc_old, t0, t1, t2
    integer :: i

    call random_inputs(q1, q2, P1, P2, k1, k2, F1, F2, cA, cB, cC)
    acc_new = zero
    acc_old = zero
    call cpu_time(t0)
    do i = 1, ntime
       q1(0) = q1(0) + 1e-9_dp
       acc_new = acc_new + chain_new(q1, q2, P1, P2, k1, k2, F1, F2, cA, cB, cC)
    end do
    call cpu_time(t1)
    do i = 1, ntime
       q1(0) = q1(0) - 1e-9_dp
       acc_old = acc_old + chain_legacy(q1, q2, P1, P2, k1, k2, F1, F2, cA, cB, cC)
    end do
    call cpu_time(t2)
    write(*,'(a,i0,a,f8.3,a,f8.3,a,f6.2)') ' timing of ', ntime, ' chains: new ', t1-t0, &
         & ' s, legacy ', t2-t1, ' s, speed-up ', (t2-t1)/max(t1-t0, 1e-9_dp)
    ! Use the results, so that the loops are not optimised away
    if (acc_new /= acc_new .or. acc_old /= acc_old) call fail('timing loop gave NaN')
  end subroutine time_chain

  real(dp) function chain_new(q1, q2, P1, P2, k1, k2, F1, F2, cA, cB, cC) result(res)
    real(dp), intent(in)    :: q1(0:3), q2(0:3), P1(0:3), P2(0:3), k1(0:3), k2(0:3), F1(3), F2(3)
    complex(dp), intent(in) :: cA, cB, cC
    type(tensors) :: g, q1mu, q2mu, P1mu, P2mu, k1mu, k2mu, P1hat, P2hat
    type(tensors) :: M, Mstar, W1, W2, T3(2), MW1, MstarW2, dummy

    g = gmunu
    g%up = .false.
    call vec(q1mu, q1); call vec(q2mu, q2); call vec(P1mu, P1); call vec(P2mu, P2)
    call vec(k1mu, k1); call vec(k2mu, k2)
    P1hat = P1mu - dot(P1,q1)/dot(q1,q1) * q1mu
    P2hat = P2mu - dot(P2,q2)/dot(q2,q2) * q2mu

    M = cA * g
    M = M + cB * ((two*k1mu+q1mu).otimes.(k2mu-k1mu-q1mu))
    M = M + cC * ((two*k2mu+q1mu).otimes.(k1mu-k2mu-q1mu))
    call raise(M, 1)
    call raise(M, 2)
    Mstar = M
    Mstar%values = conjg(M%values)

    call InitTensor(T3(1), 2, .false.)
    call InitTensor(T3(2), 2, .false.)
    call eps_contract(P1, q1, T3(1)%values)
    call eps_contract(P2, q2, T3(2)%values)

    W1 = F1(1)*((one/dot(q1,q1))*(q1mu.otimes.q1mu)-g) &
         & + F1(2)*(one/dot(P1,q1))*(P1hat.otimes.P1hat) &
         & + F1(3)*(one/(two*dot(P1,q1)))*T3(1)
    W2 = F2(1)*((one/dot(q2,q2))*(q2mu.otimes.q2mu)-g) &
         & + F2(2)*(one/dot(P2,q2))*(P2hat.otimes.P2hat) &
         & + F2(3)*(one/(two*dot(P2,q2)))*T3(2)

    call ContractTensors(W1, 1, M, 1, MW1)
    call ContractTensors(Mstar, 2, W2, 2, MstarW2)
    call ContractTensors(MW1, 1, MstarW2, 1, dummy)
    res = real(TensorTrace(dummy), kind=dp)
  end function chain_new

  real(dp) function chain_legacy(q1, q2, P1, P2, k1, k2, F1, F2, cA, cB, cC) result(res)
    real(dp), intent(in)    :: q1(0:3), q2(0:3), P1(0:3), P2(0:3), k1(0:3), k2(0:3), F1(3), F2(3)
    complex(dp), intent(in) :: cA, cB, cC
    type(ltensors) :: g, q1mu, q2mu, P1mu, P2mu, k1mu, k2mu, P1hat, P2hat
    type(ltensors) :: M, Mstar, W1, W2, T3(2), MW1, MstarW2, dummy

    g = lgmunu
    g%up = .false.
    call lvec(q1mu, q1); call lvec(q2mu, q2); call lvec(P1mu, P1); call lvec(P2mu, P2)
    call lvec(k1mu, k1); call lvec(k2mu, k2)
    P1hat = P1mu - dot(P1,q1)/dot(q1,q1) * q1mu
    P2hat = P2mu - dot(P2,q2)/dot(q2,q2) * q2mu

    M = cA * g
    M = M + cB * ((two*k1mu+q1mu).otimes.(k2mu-k1mu-q1mu))
    M = M + cC * ((two*k2mu+q1mu).otimes.(k1mu-k2mu-q1mu))
    call lraise(M, 1)
    call lraise(M, 2)
    Mstar = M
    Mstar%values = conjg(M%values)

    call lInitTensor(T3(1), 2, .false.)
    call lInitTensor(T3(2), 2, .false.)
    call eps_contract(P1, q1, T3(1)%values)
    call eps_contract(P2, q2, T3(2)%values)

    W1 = F1(1)*((one/dot(q1,q1))*(q1mu.otimes.q1mu)-g) &
         & + F1(2)*(one/dot(P1,q1))*(P1hat.otimes.P1hat) &
         & + F1(3)*(one/(two*dot(P1,q1)))*T3(1)
    W2 = F2(1)*((one/dot(q2,q2))*(q2mu.otimes.q2mu)-g) &
         & + F2(2)*(one/dot(P2,q2))*(P2hat.otimes.P2hat) &
         & + F2(3)*(one/(two*dot(P2,q2)))*T3(2)

    call lContractTensors(W1, 1, M, 1, MW1)
    call lContractTensors(Mstar, 2, W2, 2, MstarW2)
    call lContractTensors(MW1, 1, MstarW2, 1, dummy)
    res = real(lTensorTrace(dummy), kind=dp)
  end function chain_legacy

  !----------------------------------------------------------------------
  ! i epsilon_{mu nu rho sigma} P^rho q^sigma, as in the matrix elements
  subroutine eps_contract(P, q, T)
    real(dp),    intent(in)  :: P(0:3), q(0:3)
    complex(dp), intent(out) :: T(0:3,0:3)
    integer :: i, j
    T = zero
    T(0,1) =   P(2)*q(3) - P(3)*q(2)
    T(0,2) = - P(1)*q(3) + P(3)*q(1)
    T(0,3) =   P(1)*q(2) - P(2)*q(1)
    T(1,2) = - P(3)*q(0) + P(0)*q(3)
    T(1,3) =   P(2)*q(0) - P(0)*q(2)
    T(2,3) = - P(1)*q(0) + P(0)*q(1)
    do i = 0, 3
       do j = i, 3
          T(j,i) = - T(i,j)
       end do
    end do
    T = cmplx(zero, one, kind=dp) * T
  end subroutine eps_contract

  subroutine random_inputs(q1, q2, P1, P2, k1, k2, F1, F2, cA, cB, cC)
    real(dp),    intent(out) :: q1(0:3), q2(0:3), P1(0:3), P2(0:3), k1(0:3), k2(0:3), F1(3), F2(3)
    complex(dp), intent(out) :: cA, cB, cC
    real(dp) :: r(6)
    call random_number(q1); call random_number(q2)
    call random_number(P1); call random_number(P2)
    call random_number(k1); call random_number(k2)
    q1 = 200.0_dp*(q1 - 0.5_dp); q2 = 200.0_dp*(q2 - 0.5_dp)
    P1 = 1000.0_dp*P1; P2 = 1000.0_dp*P2
    k1 = 300.0_dp*(k1 - 0.5_dp); k2 = 300.0_dp*(k2 - 0.5_dp)
    call random_number(F1); call random_number(F2)
    F1 = F1 - 0.3_dp; F2 = F2 - 0.3_dp
    call random_number(r)
    cA = cmplx(r(1)-0.5_dp, r(2)-0.5_dp, kind=dp)
    cB = cmplx(r(3)-0.5_dp, r(4)-0.5_dp, kind=dp) * 1e-3_dp
    cC = cmplx(r(5)-0.5_dp, r(6)-0.5_dp, kind=dp) * 1e-3_dp
  end subroutine random_inputs

  !----------------------------------------------------------------------
  ! Helpers
  subroutine init_random()
    integer :: n
    integer, allocatable :: seed(:)
    call random_seed(size=n)
    allocate(seed(n))
    seed = 20260924
    call random_seed(put=seed)
  end subroutine init_random

  subroutine random_up(upv)
    logical, intent(out) :: upv(2)
    real(dp) :: r(2)
    call random_number(r)
    upv = r > 0.5_dp
  end subroutine random_up

  ! A new-style and a legacy tensor with the same random components
  subroutine random_pair(t, lt, rank, upv)
    type(tensors),  intent(inout) :: t
    type(ltensors), intent(inout) :: lt
    integer, intent(in) :: rank
    logical, intent(in) :: upv(2)
    real(dp) :: re(0:3,0:3), im(0:3,0:3)
    call random_number(re); call random_number(im)
    call InitTensor(t, rank, .true.)
    call lInitTensor(lt, rank, .true.)
    if (rank >= 1) then
       t%up(1:rank) = upv(1:rank)
       lt%up(1:rank) = upv(1:rank)
    end if
    select case (rank)
    case (0)
       t%values(1,1) = cmplx(re(0,0)-0.5_dp, im(0,0)-0.5_dp, kind=dp)
       lt%values(1,1) = t%values(1,1)
    case (1)
       t%values(:,1) = cmplx(re(:,0)-0.5_dp, im(:,0)-0.5_dp, kind=dp)
       lt%values(:,1) = t%values(:,1)
    case (2)
       t%values = cmplx(re-0.5_dp, im-0.5_dp, kind=dp)
       lt%values = t%values
    end select
  end subroutine random_pair

  subroutine vec(t, p)
    type(tensors), intent(inout) :: t
    real(dp), intent(in) :: p(0:3)
    call InitTensor(t, 1, .false.)
    t%values(:,1) = p
  end subroutine vec

  subroutine lvec(t, p)
    type(ltensors), intent(inout) :: t
    real(dp), intent(in) :: p(0:3)
    call lInitTensor(t, 1, .false.)
    t%values(:,1) = p
  end subroutine lvec

  ! Component of a rank-2 tensor with the contracted index (position
  ! idx) equal to i and the other one equal to j
  complex(dp) function elem(t, idx, i, j)
    type(tensors), intent(in) :: t
    integer, intent(in) :: idx, i, j
    if (idx == 1) then
       elem = t%values(i,j)
    else
       elem = t%values(j,i)
    end if
  end function elem

  function contracted_up(up1, r1, i1, up2, r2, i2) result(upv)
    logical, intent(in) :: up1(2), up2(2)
    integer, intent(in) :: r1, i1, r2, i2
    logical :: upv(2)
    integer :: i, n
    upv = .true.
    n = 0
    do i = 1, r1
       if (i /= i1) then
          n = n + 1; upv(n) = up1(i)
       end if
    end do
    do i = 1, r2
       if (i /= i2) then
          n = n + 1; upv(n) = up2(i)
       end if
    end do
  end function contracted_up

  real(dp) function dot(p, q)
    real(dp), intent(in) :: p(0:3), q(0:3)
    dot = p(0)*q(0) - sum(p(1:3)*q(1:3))
  end function dot

  real(dp) function reldev(a, b)
    complex(dp), intent(in) :: a(:,:), b(:,:)
    real(dp) :: scale
    scale = max(maxval(abs(b)), tiny(one))
    reldev = maxval(abs(a - b)) / scale
  end function reldev

  real(dp) function reldev_scalar(a, b)
    complex(dp), intent(in) :: a, b
    reldev_scalar = abs(a - b) / max(abs(b), tiny(one))
  end function reldev_scalar

  ! The components of t that belong to its rank
  function used(t) result(v)
    type(tensors), intent(in) :: t
    complex(dp), allocatable :: v(:,:)
    select case (t%rank)
    case (0); v = t%values(1:1,1:1)
    case (1); v = t%values(:,1:1)
    case default; v = t%values
    end select
  end function used

  function lused(t) result(v)
    type(ltensors), intent(in) :: t
    complex(dp), allocatable :: v(:,:)
    v = t%values
  end function lused

  subroutine check_legacy(t, lt, what)
    type(tensors),  intent(in) :: t
    type(ltensors), intent(in) :: lt
    character(len=*), intent(in) :: what
    real(dp) :: d
    if (t%rank /= lt%rank) then
       call fail(what//' (rank differs from legacy)'); return
    end if
    if (t%rank > 0) then
       if (any(t%up(1:t%rank) .neqv. lt%up(1:t%rank))) call fail(what//' (index positions differ from legacy)')
    end if
    d = reldev(used(t), lused(lt))
    maxdev_legacy = max(maxdev_legacy, d)
    if (d > tol) call fail(what//' (values differ from legacy)')
  end subroutine check_legacy

  subroutine check_formula(t, ref, rank, what)
    type(tensors), intent(in) :: t
    complex(dp),   intent(in) :: ref(0:3,0:3)
    integer,       intent(in) :: rank
    character(len=*), intent(in) :: what
    real(dp) :: d
    if (t%rank /= rank) then
       call fail(what//' (wrong rank)'); return
    end if
    select case (rank)
    case (0); d = reldev(t%values(1:1,1:1), ref(1:1,1:1))
    case (1); d = reldev(t%values(:,1:1), ref(:,1:1))
    case default; d = reldev(t%values, ref)
    end select
    maxdev_formula = max(maxdev_formula, d)
    if (d > tol) call fail(what//' (values differ from explicit formula)')
  end subroutine check_formula

  subroutine check_up(t, upv, what)
    type(tensors), intent(in) :: t
    logical,       intent(in) :: upv(2)
    character(len=*), intent(in) :: what
    if (t%rank > 0) then
       if (any(t%up(1:t%rank) .neqv. upv(1:t%rank))) call fail(what//' (wrong index positions)')
    end if
  end subroutine check_up

  subroutine fail(what)
    character(len=*), intent(in) :: what
    nfail = nfail + 1
    if (nfail <= 20) write(*,'(a,a)') ' FAIL: ', what
  end subroutine fail

end program test_tensor
