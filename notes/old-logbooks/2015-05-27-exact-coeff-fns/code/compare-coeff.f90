program compare_coeff
  use types; use consts_dp
  use xc2ns2e
  use xc3ns2e
  use xclns2e
  use xcdiff2e
  
  use xc2ns2p
  use xc3ns2p
  use xclns2p
  implicit none
  !----------------------------------------------------------------------
  integer :: ix
  integer, parameter :: nx = 25, nf = 5
  real(dp), parameter :: xmin = 1e-3_dp, xmax = 0.99_dp
  real(dp) :: x, res_exact, res_apprx

  call set_c2soft(nf)
  call set_c3soft(nf)

  !---------------------------------------------------------------
  write(6,'(a)') "# F2NSplus-A"
  do ix = 0, nx
     x = xmin*(xmax/xmin)**(ix*one/nx)
     res_exact = x2np2a(x,nf)
     res_apprx = c2nn2a(x,nf)
     write(6,*) x, res_exact, res_apprx
  end do
  write(6,*); write(6,*)

  write(6,'(a)') "# F2NSplus-B"
  do ix = 0, nx
     x = xmin*(xmax/xmin)**(ix*one/nx)
     res_exact = x2ns2b(x,nf)
     res_apprx = c2ns2b(x,nf)
     write(6,*) x, res_exact, res_apprx
  end do
  write(6,*); write(6,*)

  !---------------------------------------------------------------
  write(6,'(a)') "# F2NSminus-A"
  do ix = 0, nx
     x = xmin*(xmax/xmin)**(ix*one/nx)
     res_exact = x2np2a(x,nf) - XC2DFF2(x)
     res_apprx = c2nc2a(x,nf)
     write(6,*) x, res_exact, res_apprx
  end do
  write(6,*); write(6,*)

  write(6,'(a)') "# F2NSminus-B"
  do ix = 0, nx
     x = xmin*(xmax/xmin)**(ix*one/nx)
     res_exact = x2ns2b(x,nf)
     res_apprx = c2ns2b(x,nf)
     write(6,*) x, res_exact, res_apprx
  end do
  write(6,*); write(6,*)

  !---------------------------------------------------------------
  write(6,'(a)') "# FLNSplus-A"
  do ix = 0, nx
     x = xmin*(xmax/xmin)**(ix*one/nx)
     res_exact = xlnp2a(x,nf)
     res_apprx = clnn2a(x,nf)
     write(6,*) x, res_exact, res_apprx
  end do
  write(6,*); write(6,*)
  !---------------------------------------------------------------
  write(6,'(a)') "# FLNSminus-A"
  do ix = 0, nx
     x = xmin*(xmax/xmin)**(ix*one/nx)
     res_exact = xlnp2a(x,nf) - xcldff2(x)
     res_apprx = clnc2a(x,nf)
     write(6,*) x, res_exact, res_apprx
  end do
  write(6,*); write(6,*)
  

  !---------------------------------------------------------------
  write(6,'(a)') "# F3NSplus-A"
  do ix = 0, nx
     x = xmin*(xmax/xmin)**(ix*one/nx)
     res_exact = x3nm2a(x,nf) + XC3DFF2(x)
     res_apprx = c3np2a(x,nf)
     write(6,*) x, res_exact, res_apprx
  end do
  write(6,*); write(6,*)

  write(6,'(a)') "# F3NSplus-B"
  do ix = 0, nx
     x = xmin*(xmax/xmin)**(ix*one/nx)
     res_exact = x3ns2b(x,nf)
     res_apprx = c3ns2b(x,nf)
     write(6,*) x, res_exact, res_apprx
  end do
  write(6,*); write(6,*)

  !---------------------------------------------------------------
  write(6,'(a)') "# F3NSminus-A"
  do ix = 0, nx
     x = xmin*(xmax/xmin)**(ix*one/nx)
     res_exact = x3nm2a(x,nf)
     res_apprx = c3nm2a(x,nf)
     write(6,*) x, res_exact, res_apprx
  end do
  write(6,*); write(6,*)

  write(6,'(a)') "# F3NSminus-B"
  do ix = 0, nx
     x = xmin*(xmax/xmin)**(ix*one/nx)
     res_exact = x3ns2b(x,nf)
     res_apprx = c3ns2b(x,nf)
     write(6,*) x, res_exact, res_apprx
  end do
  write(6,*); write(6,*)

  !---------------------------------------------------------------
  write(6,'(a)') "# F3NSWrongSignDiffplus-A"
  do ix = 0, nx
     x = xmin*(xmax/xmin)**(ix*one/nx)
     res_exact = x3nm2a(x,nf) - XC3DFF2(x)
     res_apprx = c3np2a(x,nf)
     write(6,*) x, res_exact, res_apprx
  end do
  write(6,*); write(6,*)


end program compare_coeff
