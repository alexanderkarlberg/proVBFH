program csp
  use incl_vbfh, only: inclusive_init
  implicit none
  double precision :: p(0:3,6), m1, m2
  integer :: f(6,6), j
  data f / 1,2,25,2,1,0,  1,4,25,2,3,0,  2,1,25,1,2,0,  1,1,25,1,1,0,  1,2,25,1,2,0,  1,3,25,1,3,0 /
  call inclusive_init()
  call cs_init_vbfnlo()
  call testpoint(p)
  do j=1,6
     call cs_hjjj_lines(p, f(:,j), m1, m2)
     write(*,'(a,6i4,a,es16.8,a,2es14.6)') ' flav', f(:,j), '  m1+m2 ', m1+m2, '  m1,m2 ', m1, m2
  enddo
end program
