! This is a module to facilitate computations with tensors. Currently
! we are working in flat metric but one can easily generalise by
! modifying gmunu. The maximum rank of a tensor is currently 2.
!
! The components live in fixed-size arrays, so creating, copying and
! combining tensors never allocates memory. A rank-0 tensor is stored
! in values(1,1), a rank-1 tensor in values(0:3,1) and a rank-2
! tensor in values(0:3,0:3); the remaining entries are not used.
!
! The same file is used in proVBFH, proVBFHH and proVBFH-inclusive.
! tests/test_tensor.f90 in proVBFH-inclusive checks it against the
! previous, allocatable implementation (tests/tensor_legacy.f90).

module tensor
  use hoppet
  implicit none

  integer, parameter :: rankmax = 2 ! Currently this is the maximum rank implemented

  private :: CheckInitialised
  private :: CheckSameStructure
  private :: CheckIndex

  public :: InitTensor
  public :: ResetTensor
  public :: SetMetric
  public :: PrintTensor
  public :: ContractTensors
  public :: raise
  public :: lower
  public :: AddTensors
  public :: SubtractTensors
  public :: MultiplyWithScalar
  public :: MultiplyWithComplexScalar
  public :: TensorProduct
  public :: TensorTrace

  type, public :: tensors
     integer :: rank = 0
     logical :: up(rankmax) = .true. ! True if index at position i is up
                                     ! and false if it is down
     complex(dp) :: values(0:3,0:3) = (0.0_dp,0.0_dp) ! The actual values of the tensor
     logical :: initialised=.false.
  end type tensors

  INTERFACE OPERATOR (+)
     PROCEDURE AddTensors
  END INTERFACE OPERATOR (+)

  INTERFACE OPERATOR (*)
     PROCEDURE MultiplyWithScalar
     PROCEDURE MultiplyWithComplexScalar
  END INTERFACE OPERATOR (*)

  INTERFACE OPERATOR (-)
     PROCEDURE SubtractTensors
  END INTERFACE OPERATOR (-)

  INTERFACE OPERATOR (.otimes.)
     PROCEDURE TensorProduct
  END INTERFACE OPERATOR (.otimes.)

  type(tensors) :: gmunu, emunurhosigma
contains
  subroutine InitTensor(tensor,rank,up)
    integer rank
    logical up
    type(tensors) :: tensor

    if(rank.lt.0) then
       print*, 'Negative rank. Exiting'
       stop
    elseif(rank.gt.rankmax) then
       print*, 'rank not implemented', rank
       stop
    endif

    tensor%rank = rank
    tensor%up(:) = up
    tensor%values = zero
    tensor%initialised = .true.
  end subroutine InitTensor

  subroutine ResetTensor(tensor)
    type(tensors) :: tensor
    if(.not.tensor%initialised) then
       print*, 'Tensor not initialised. Exiting.'
       stop
    endif

    tensor%values = zero
  end subroutine ResetTensor

  ! This routine initialises the metric with the specified signature
  ! of the 00 component.
  subroutine SetMetric(signature)
    integer signature
    integer, parameter :: rank = 2
    logical up

    up = .true.
    if(.not.gmunu%initialised) then
       call initTensor(gmunu,rank,up)
    endif
    call ResetTensor(gmunu)

    gmunu%values(0,0) = one * signature
    gmunu%values(1,1) = -one * signature
    gmunu%values(2,2) = -one * signature
    gmunu%values(3,3) = -one * signature
  end subroutine SetMetric

  ! Routine to print the tensor in a readable format
  subroutine PrintTensor(tensor)
    type(tensors) :: tensor

    if(tensor%rank.eq.0) then
       write(*,'(4("[",d11.4,", ",d11.4,"I","] "))') tensor%values(1,1)
    elseif(tensor%rank.eq.1) then
       write(*,'(4("[",d11.4,", ",d11.4,"I","] "))') tensor%values(:,1)
    else
       write(*,'(4("[",d11.4,", ",d11.4,"I","] "))') tensor%values
    endif
    write(*,*) ''
    write(*,*) 'Index structure is', tensor%up(1:tensor%rank)
    write(*,*) ''

  end subroutine PrintTensor

  ! This routine is the bulk of the work. It takes two tensors and
  ! contract them into one tensor. The result is accumulated locally,
  ! so tout may be one of the inputs.
  subroutine ContractTensors(tin1,index1,tin2,index2,tout)
    type(tensors), intent(in) :: tin1,tin2
    integer, intent(in) :: index1,index2
    type(tensors), intent(inout) :: tout
    complex(dp) :: res(0:3,0:3)
    logical :: up(rankmax)
    integer :: i,j,k,rank

    ! First some sanity checks
    call CheckInitialised(tin1)
    call CheckInitialised(tin2)
    if(index1.lt.1.or.index1.gt.tin1%rank) then
       print*, 'index1 out of bounds', index1, tin1%rank
       stop
    endif
    if(index2.lt.1.or.index2.gt.tin2%rank) then
       print*, 'index2 out of bounds', index2, tin2%rank
       stop
    endif

    rank = tin1%rank+tin2%rank-2
    if(rank.gt.rankmax) then
       print*, 'Trying to contract two tensors into a tensor of too high rank', &
            & tin1%rank, tin2%rank, rank
       stop
    endif

    ! The routine assumes that the two indices are not in the same
    ! position.
    if(tin1%up(index1).eqv.tin2%up(index2)) then
       print*, 'ERROR: Trying to contract two indices in the same position.', &
            & tin1%up(index1), tin2%up(index2)
       stop
    endif

    ! Here we contract the tensors. Right now I don't see a smart way
    ! of doing it, so we are brute-forcing the problem as rank-2
    ! tensors are manageable. If we want to extend to higher ranks it
    ! may be useful to return here...
    res = zero
    if(tin1%rank.eq.1.and.tin2%rank.eq.1) then ! rank-1 with rank-1
       do i=0,3
          res(1,1) = res(1,1) + tin1%values(i,1)*tin2%values(i,1)
       enddo
    elseif(tin1%rank.eq.1.and.tin2%rank.eq.2) then ! rank-1 with rank-2
       if(index2.eq.1) then
          do i=0,3
             do j=0,3
                res(j,1) = res(j,1) + tin1%values(i,1)*tin2%values(i,j)
             enddo
          enddo
       else
          do i=0,3
             do j=0,3
                res(j,1) = res(j,1) + tin1%values(i,1)*tin2%values(j,i)
             enddo
          enddo
       endif
    elseif(tin1%rank.eq.2.and.tin2%rank.eq.1) then ! rank-2 with rank-1
       if(index1.eq.1) then
          do i=0,3
             do j=0,3
                res(j,1) = res(j,1) + tin1%values(i,j)*tin2%values(i,1)
             enddo
          enddo
       else
          do i=0,3
             do j=0,3
                res(j,1) = res(j,1) + tin1%values(j,i)*tin2%values(i,1)
             enddo
          enddo
       endif
    elseif(tin1%rank.eq.2.and.tin2%rank.eq.2) then ! rank-2 with rank-2
       if(index1.eq.1.and.index2.eq.1) then
          do i=0,3
             do k=0,3
                do j=0,3
                   res(j,k) = res(j,k) + tin1%values(i,j)*tin2%values(i,k)
                enddo
             enddo
          enddo
       elseif(index1.eq.1.and.index2.eq.2) then
          do i=0,3
             do k=0,3
                do j=0,3
                   res(j,k) = res(j,k) + tin1%values(i,j)*tin2%values(k,i)
                enddo
             enddo
          enddo
       elseif(index1.eq.2.and.index2.eq.1) then
          do i=0,3
             do k=0,3
                do j=0,3
                   res(j,k) = res(j,k) + tin1%values(j,i)*tin2%values(i,k)
                enddo
             enddo
          enddo
       else
          do i=0,3
             do k=0,3
                do j=0,3
                   res(j,k) = res(j,k) + tin1%values(j,i)*tin2%values(k,i)
                enddo
             enddo
          enddo
       endif
    endif

    ! Set the index position of the new tensor according to the old tensor.
    up = .true.
    do i=1,tin1%rank
       if(i.lt.index1) then
          up(i) = tin1%up(i)
       elseif(i.gt.index1) then
          up(i-1) = tin1%up(i)
       endif
    enddo

    do i=1,tin2%rank
       if(i.lt.index2) then
          up(tin1%rank-1+i) = tin2%up(i)
       elseif(i.gt.index2) then
          up(tin1%rank-1+i-1) = tin2%up(i)
       endif
    enddo

    tout%rank = rank
    tout%up = up
    tout%values = res
    tout%initialised = .true.
  end subroutine ContractTensors

  ! This routine raises the index of a tensor
  subroutine raise(tensor,index)
    type(tensors) :: tensor
    integer :: index, i

    call CheckIndex(tensor,index,'raise')
    if(.not.tensor%up(index)) then
       ! For now the metric is always diagonal. This saves time...
       if(index.eq.1) then
          do i=0,3
             tensor%values(i,:) = gmunu%values(i,i)*tensor%values(i,:)
          enddo
       else
          do i=0,3
             tensor%values(:,i) = gmunu%values(i,i)*tensor%values(:,i)
          enddo
       endif
       tensor%up(index) = .true.
    endif
  end subroutine raise

  ! This routine lowers the index of a tensor
  subroutine lower(tensor,index)
    type(tensors) :: tensor
    integer :: index, i

    call CheckIndex(tensor,index,'lower')
    if(tensor%up(index)) then
       ! For now the metric is always diagonal. This saves time...
       if(index.eq.1) then
          do i=0,3
             tensor%values(i,:) = gmunu%values(i,i)*tensor%values(i,:)
          enddo
       else
          do i=0,3
             tensor%values(:,i) = gmunu%values(i,i)*tensor%values(:,i)
          enddo
       endif
       tensor%up(index) = .false.
    endif
  end subroutine lower

  ! This routine checks if a tensor has been initialised.
  subroutine CheckInitialised(tensor)
    type(tensors), intent(in) :: tensor

    if(.not.tensor%initialised) then
       print*, 'tensor not initialised. Exiting.'
       stop
    endif
  end subroutine CheckInitialised

  ! Checks before raising or lowering an index
  subroutine CheckIndex(tensor,index,operation)
    type(tensors), intent(in) :: tensor
    integer, intent(in) :: index
    character(len=*), intent(in) :: operation

    call CheckInitialised(tensor)
    if(.not.gmunu%initialised) then
       print*, 'Cannot ', operation, ' index. Metric not set (call SetMetric)'
       stop
    endif
    if(index.lt.1.or.index.gt.tensor%rank) then
       print*, 'Cannot ', operation, ' index. Index out of bounds', index, tensor%rank
       stop
    endif
  end subroutine CheckIndex

  ! Checks that two tensors can be added or subtracted
  subroutine CheckSameStructure(t1,t2,operation)
    type(tensors), intent(in) :: t1,t2
    character(len=*), intent(in) :: operation

    call CheckInitialised(t1)
    call CheckInitialised(t2)

    if(t1%rank.ne.t2%rank) then
       print*, 'Trying to ', operation, ' two tensors with different ranks', t1%rank, t2%rank
       stop
    endif
    if(any(t1%up(1:t1%rank).neqv.t2%up(1:t1%rank))) then
       print*, 'Trying to ', operation, ' two tensors with different index structure', &
            & t1%up(1:t1%rank), t2%up(1:t1%rank)
       stop
    endif
  end subroutine CheckSameStructure

  function AddTensors(t1,t2) result(res)
    type(tensors), intent(in) :: t1,t2
    type(tensors) :: res

    call CheckSameStructure(t1,t2,'add')

    res%rank = t1%rank
    res%up = t1%up
    res%initialised = .true.
    res%values = t1%values + t2%values
  end function AddTensors

  function SubtractTensors(t1,t2) result(res)
    type(tensors), intent(in) :: t1,t2
    type(tensors) :: res

    call CheckSameStructure(t1,t2,'subtract')

    res%rank = t1%rank
    res%up = t1%up
    res%initialised = .true.
    res%values = t1%values - t2%values
  end function SubtractTensors

  function MultiplyWithScalar(r1,t1) result(res)
    type(tensors), intent(in) :: t1
    real(dp),intent(in) :: r1
    type(tensors) :: res

    call CheckInitialised(t1)

    res%rank = t1%rank
    res%up = t1%up
    res%initialised = .true.
    res%values = t1%values * r1
  end function MultiplyWithScalar

  function MultiplyWithComplexScalar(c1,t1) result(res)
    type(tensors), intent(in) :: t1
    complex(dp),intent(in) :: c1
    type(tensors) :: res

    call CheckInitialised(t1)

    res%rank = t1%rank
    res%up = t1%up
    res%initialised = .true.
    res%values = t1%values * c1
  end function MultiplyWithComplexScalar

  function TensorProduct(t1,t2) result(res)
    type(tensors), intent(in) :: t1,t2
    type(tensors) :: res
    integer :: i,j

    call CheckInitialised(t1)
    call CheckInitialised(t2)

    if((t1%rank+t2%rank).gt.rankmax) then
       print*, 'Rank not supported!', t1%rank+t2%rank
       stop
    endif

    res%rank = t1%rank+t2%rank
    res%initialised = .true.

    ! Special case if one of them is rank 0
    if(t1%rank.eq.0) then
       res%up = t2%up
       res%values = t1%values(1,1) * t2%values
       return
    elseif(t2%rank.eq.0) then
       res%up = t1%up
       res%values = t2%values(1,1) * t1%values
       return
    endif

    ! Otherwise both are rank 1
    res%up(1) = t1%up(1)
    res%up(2) = t2%up(1)
    do j = 0,3
       do i = 0,3
          res%values(i,j) = t1%values(i,1) * t2%values(j,1)
       enddo
    enddo

  end function TensorProduct

  function TensorTrace(tensor) result(res)
    type(tensors), intent(in) :: tensor
    complex(dp) :: res

    if(tensor%rank.ne.2) then
       print*, 'Trace only implemented for rank 2 tensor'
       stop
    endif
    res = tensor%values(0,0) + tensor%values(1,1) + tensor%values(2,2) + tensor%values(3,3)
  end function TensorTrace

end module tensor
