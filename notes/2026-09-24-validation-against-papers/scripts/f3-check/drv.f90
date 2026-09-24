program drv
  use ME_expressions
  implicit none
  real*8 :: p(0:3,6), c(6)
  character(len=2) :: lab(6) = (/'AA','AB','AC','BB','BC','CC'/)
  integer :: i
  open(10,file='point.txt'); do i=1,6; read(10,*) p(:,i); enddo; close(10)
  do i=1,6
     c=0d0; c(i)=1d0
     write(*,'(a,a,a,es16.8,a,es16.8)') ' analytic ',lab(i),': F1F1 ', &
          F1F1(c(1),c(2),c(3),c(4),c(5),c(6),p(:,1),p(:,2),p(:,3),p(:,4),p(:,5),p(:,6)), &
          '  F3F3 ', F3F3(c(1),c(2),c(3),c(4),c(5),c(6),p(:,1),p(:,2),p(:,3),p(:,4),p(:,5),p(:,6))
  enddo
end program drv
