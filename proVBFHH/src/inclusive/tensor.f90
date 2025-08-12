! This is a module to facilitate computations with tensors. Currently
! we are working in flat metric but one can easily generalise by
! modifying gmunu. The maximum rank of a tensor is currently 2.

module tensor
  use hoppet
  implicit none

  private :: CheckInitialised

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
  public :: TensorProduct
  public :: TensorTrace
  
  type, public :: tensors
     ! Inputs
     integer :: rank
     logical, allocatable :: up(:) ! True if index at position i is up
                                   ! and false if it is down
     complex(dp), allocatable :: values(:,:) ! The actual values of the tensor
     logical :: initialised=.false.
   contains
!     procedure :: trace => TensorTrace
!     procedure :: reset => ResetTensor
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
  integer, parameter :: rankmax = 2 ! Currently this is the maximum rank implemented
contains
  subroutine InitTensor(tensor,rank,up)
    integer rank
    logical up
    type(tensors) :: tensor

    if(tensor%initialised) then
       deallocate(tensor%up)
       deallocate(tensor%values)
       !       print*, 'Tensor already initialised. Exiting'
       !       stop
    endif

    if(rank.lt.0) then
       print*, 'Zero or negative rank. Exiting'
       stop
    endif
    tensor%rank = rank
    allocate(tensor%up(1:rank))

    if(rank.eq.0) then
       allocate(tensor%values(1,1))
    elseif(rank.eq.1) then
       allocate(tensor%values(0:3,1))
    elseif(rank.eq.2) then
       allocate(tensor%values(0:3,0:3))
    else
       print*, 'rank not implemented', rank
       stop
    endif
    
    tensor%up(:) = up
    tensor%values = zero
    tensor%initialised = .true.

!    print*, 'InitTensor called with ', rank,up
  end subroutine InitTensor

  subroutine ResetTensor(tensor)
    type(tensors) :: tensor
    if(.not.tensor%initialised) then
       print*, 'Tensor not initialised. Exiting.'
       stop
    endif
    
    tensor%values = zero
!    tensor%up = .true.
  end subroutine ResetTensor
  
  ! This routine initialises the metric with the specified signature
  ! of the 00 component. We also initialise the levi-civita symbols
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
    
!    call PrintTensor(gmunu)

  end subroutine SetMetric

  ! Routine to print the tensor in a readable format
  subroutine PrintTensor(tensor)
    type(tensors) :: tensor

!    write(*,'(4(d14.5))') tensor%values
    write(*,'(4("[",d11.4,", ",d11.4,"I","] "))') tensor%values
    write(*,*) ''
    write(*,*) 'Index structure is', tensor%up
    write(*,*) ''
    
  end subroutine PrintTensor
  
  ! This routine is the bulk of the work. It takes two tensors and
  ! contract them into one tensor.
  subroutine ContractTensors(tin1,index1,tin2,index2,tout)
    type(tensors), intent(in) :: tin1,tin2
    integer, intent(in) :: index1,index2
    type(tensors) :: tout
    type(tensors) :: t1,t2
    integer :: i,j,k,l,m

    ! First some sanity checks    
    call CheckInitialised(tin1)
    call CheckInitialised(tin2)
!    call CheckInitialised(tout)
!    call ResetTensor(tout)
    if(index1.gt.tin1%rank) then
       print*, 'index1 out of bounds', index1, tin1%rank
       stop
    endif
    if(index2.gt.tin2%rank) then
       print*, 'index2 out of bounds', index2, tin2%rank
       stop
    endif

    if(.not.tout%initialised) then
       call InitTensor(tout,tin1%rank+tin2%rank-2,.true.)
    else
       call ResetTensor(tout)
    endif
    
    if((tin1%rank+tin2%rank-2).gt.rankmax) then
       print*, 'Trying to contract two tensors into a tensor of too high rank', &
            & tin1%rank, tin2%rank, tout%rank
       stop
    endif
    
    ! Local copies of inputs
    t1 = tin1
    t2 = tin2

    ! At this point we need to check the status of the various
    ! indeces. The routine will assume that the two indeces are not in
    ! the same position, but if they are, the routine will
    ! automatically raise/lower an index.

    if(t1%up(index1).eqv.t2%up(index2)) then
       print*, 'ERROR: Trying to contract two indices in the same position.', t1%up(index1), t2%up(index2)
       stop
!       if(t1%up(index1)) then
!          call lower(t1,index1)
!       else
!          call raise(t1,index1)
!       endif
    endif

    ! Here we contract the tensors. Right now I don't see a smart way
    ! of doing it, so we are brute-forcing the problem as rank-2
    ! tensors are manageable. If we want to extend to higher ranks it
    ! may be useful to return here...

!    print*, t1%rank, t2%rank
    if(t1%rank.eq.1.and.t2%rank.eq.1) then ! rank-1 with rank-1
       do i=0,3
          tout%values(1,1) = tout%values(1,1) + t1%values(i,1)*t2%values(i,1)
       enddo
    elseif(t1%rank.eq.1.and.t2%rank.eq.2) then ! rank-1 with rank-2
       if(index2.eq.1) then
          do i=0,3
             do j=0,3
                tout%values(j,1) = tout%values(j,1) +  t1%values(i,1)*t2%values(i,j)
             enddo
          enddo
       elseif(index2.eq.2) then
          do i=0,3
             do j=0,3
                tout%values(j,1) = tout%values(j,1) +  t1%values(i,1)*t2%values(j,i)
             enddo
          enddo
       endif
    elseif(t1%rank.eq.2.and.t2%rank.eq.1) then ! rank-2 with rank-1
       if(index1.eq.1) then
          do i=0,3
             do j=0,3
                tout%values(j,1) = tout%values(j,1) +  t1%values(i,j)*t2%values(i,1)
             enddo
          enddo
       elseif(index1.eq.2) then
          do i=0,3
             do j=0,3
                tout%values(j,1) = tout%values(j,1) +  t1%values(j,1)*t2%values(i,1)
             enddo
          enddo
       endif
    elseif(t1%rank.eq.2.and.t2%rank.eq.2) then ! rank-2 with rank-2
       if(index1.eq.1.and.index2.eq.1) then
          do i=0,3
             do k=0,3
                do j=0,3
                   tout%values(j,k) = tout%values(j,k) + t1%values(i,j)*t2%values(i,k)
                enddo
             enddo
          enddo
       elseif(index1.eq.1.and.index2.eq.2) then
          do i=0,3
             do k=0,3
                do j=0,3
                   tout%values(j,k) = tout%values(j,k) + t1%values(i,j)*t2%values(k,i)
                enddo
             enddo
          enddo
       elseif(index1.eq.2.and.index2.eq.1) then
          do i=0,3
             do k=0,3
                do j=0,3
                   tout%values(j,k) = tout%values(j,k) + t1%values(j,i)*t2%values(i,k)
                enddo
             enddo
          enddo
       elseif(index1.eq.2.and.index2.eq.2) then
          do i=0,3
             do k=0,3
                do j=0,3
                   tout%values(j,k) = tout%values(j,k) + t1%values(j,i)*t2%values(k,i)
                enddo
             enddo
          enddo
       endif
    endif
    ! Set the index position of the new tensor according to the old tensor.
    do i=1,tin1%rank
       if(i.lt.index1) then
          tout%up(i) = tin1%up(i)
       elseif(i.gt.index1) then
          tout%up(i-1) = tin1%up(i)
       endif
    enddo
    
    do i=1,tin2%rank
       if(i.lt.index2) then
          tout%up(tin1%rank-1+i) = tin2%up(i)
       elseif(i.gt.index2) then
          tout%up(tin1%rank-1+i-1) = tin2%up(i)
       endif
    enddo
    
  end subroutine ContractTensors

  ! This routine raises the index of a tensor
  subroutine raise(tensor,index)
    type(tensors) :: tensor,tout
    integer :: index, i

    call CheckInitialised(tensor)
    if(index.gt.tensor%rank) then
       print*, 'Cannot raise index. Index out of bounds', index, tensor%rank
       stop
    endif
    tout = tensor ! Copy original tensor
    if(.not.tensor%up(index)) then
!       if(index.eq.1) then
!          do i=0,3
!             tout%values(i,:) = gmunu%values(i,0)*tensor%values(0,:) &
!                  &           + gmunu%values(i,1)*tensor%values(1,:) &
!                  &           + gmunu%values(i,2)*tensor%values(2,:) &
!                  &           + gmunu%values(i,3)*tensor%values(3,:)
!          enddo
!       elseif(index.eq.2) then
!          do i=0,3
!             tout%values(:,i) = gmunu%values(0,i)*tensor%values(:,0) &
!                  &           + gmunu%values(1,i)*tensor%values(:,1) &
!                  &           + gmunu%values(2,i)*tensor%values(:,2) &
!                  &           + gmunu%values(3,i)*tensor%values(:,3)
!          enddo
!       endif
       ! For now the metric is always diagonal. This saves time...
       if(index.eq.1) then
          do i=0,3
             tout%values(i,:) = gmunu%values(i,i)*tensor%values(i,:) 
          enddo
       elseif(index.eq.2) then
          do i=0,3
             tout%values(:,i) = gmunu%values(i,i)*tensor%values(:,i) 
          enddo
       endif
       tout%up(index) = .true.
    endif
    tensor = tout
  end subroutine raise
  
  ! This routine lowers the index of a tensor
  subroutine lower(tensor,index)
    type(tensors) :: tensor,tout
    integer :: index, i

    call CheckInitialised(tensor)
    if(index.gt.tensor%rank) then
       print*, 'Cannot raise index. Index out of bounds', index, tensor%rank
       stop
    endif
    tout = tensor ! Copy original tensor
    if(tensor%up(index)) then
!       if(index.eq.1) then
!          do i=0,3
!             tout%values(i,:) = gmunu%values(i,0)*tensor%values(0,:) &
!                  &           + gmunu%values(i,1)*tensor%values(1,:) &
!                  &           + gmunu%values(i,2)*tensor%values(2,:) &
!                  &           + gmunu%values(i,3)*tensor%values(3,:)
!          enddo
!       elseif(index.eq.2) then
!          do i=0,3
!             tout%values(:,i) = gmunu%values(0,i)*tensor%values(:,0) &
!                  &           + gmunu%values(1,i)*tensor%values(:,1) &
!                  &           + gmunu%values(2,i)*tensor%values(:,2) &
!                  &           + gmunu%values(3,i)*tensor%values(:,3)
!          enddo
!       endif
          ! For now the metric is always diagonal. This saves time...
       if(index.eq.1) then
          do i=0,3
             tout%values(i,:) = gmunu%values(i,i)*tensor%values(i,:) 
          enddo
       elseif(index.eq.2) then
          do i=0,3
             tout%values(:,i) = gmunu%values(i,i)*tensor%values(:,i) 
          enddo
       endif
       tout%up(index) = .false.
    endif
    tensor = tout
  end subroutine lower
  
  ! This routine checks if a tensor has been initialised.
  subroutine CheckInitialised(tensor)
    type(tensors) :: tensor

    if(.not.tensor%initialised) then
       print*, 'tensor not initialised. Exiting.'
       stop
    endif
  end subroutine CheckInitialised

  function AddTensors(t1,t2) result(res) 
    type(tensors), intent(in) :: t1,t2
    type(tensors) :: res

    call CheckInitialised(t1)
    call CheckInitialised(t2)

    if(t1%rank.ne.t2%rank) then
       print*, 'Trying to add two tensors with different ranks', t1%rank, t2%rank
       stop
    endif
    if(any(t1%up.neqv.t2%up)) then
       print*, 'Trying to add two tensors with different index structure', t1%up, t2%up
       stop
    endif

    res = t1 ! Copy tensor

    res%values = t1%values + t2%values
    
  end function AddTensors

  function SubtractTensors(t1,t2) result(res)
    type(tensors), intent(in) :: t1,t2
    type(tensors) :: res

    call CheckInitialised(t1)
    call CheckInitialised(t2)

    if(t1%rank.ne.t2%rank) then
       print*, 'Trying to add two tensors with different ranks', t1%rank, t2%rank
       stop
    endif
    if(any(t1%up.neqv.t2%up)) then
       print*, 'Trying to subtract two tensors with different index structure', t1%up, t2%up
       stop
    endif

    res = t1 ! Copy tensor

    res%values = t1%values - t2%values
    
  end function SubtractTensors

  function MultiplyWithScalar(r1,t1) result(res) 
    type(tensors), intent(in) :: t1
    real(dp),intent(in) :: r1
    type(tensors) :: res

    call CheckInitialised(t1)

    res = t1 ! Copy tensor

    res%values = t1%values * r1
  end function MultiplyWithScalar

  function MultiplyWithComplexScalar(c1,t1) result(res) 
    type(tensors), intent(in) :: t1
    complex(dp),intent(in) :: c1
    type(tensors) :: res

    call CheckInitialised(t1)

    res = t1 ! Copy tensor

    res%values = t1%values * c1
  end function MultiplyWithComplexScalar

  function TensorProduct(t1,t2) result(res)
    type(tensors), intent(in) :: t1,t2
    type(tensors) :: res
    type(tensors), save :: rank1, rank2
    integer :: i,j,k
    
    call CheckInitialised(t1)
    call CheckInitialised(t2)
    if(.not.rank1%initialised) call InitTensor(rank1,1,.true.)
    if(.not.rank2%initialised) call InitTensor(rank2,2,.true.)
    
    if((t1%rank+t2%rank).gt.rankmax) then
       print*, 'Rank not supported!', t1%rank+t2%rank
       stop
    elseif((t1%rank+t2%rank).eq.1) then
       res = rank1
    elseif((t1%rank+t2%rank).eq.2) then
       res = rank2
    endif

    ! Initialise tensor of appropriate rank
    !call InitTensor(res,t1%rank+t2%rank,.true.)

    ! Special case if one is a rank 0 and the other rank 2.
    if(t1%rank.eq.0) then
       res%values(:,:) = t1%values(1,1) * t2%values(:,:)
       return
    elseif(t2%rank.eq.0) then
       res%values(:,:) = t2%values(1,1) * t1%values(:,:)
       return
    endif

    ! Set index structure
    do k=1,t1%rank
          res%up(k) = t1%up(k)
    enddo
    do k=1,t2%rank
          res%up(t1%rank+k) = t2%up(k)
    enddo

    ! Assume rank 1 for each of them (only supported case for the moment)
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

