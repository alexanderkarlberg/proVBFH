! Writes the Les Houches benchmark PDF (as provided by hoppet,
! pdfs_for_benchmarks), evolved with hoppet, as an LHAPDF6 grid:
!   <basename>.info, <basename>_0000.dat
! Initial condition at Q0 = sqrt(2) GeV with alpha_s(Q0) = 0.35, three-loop
! (NNLO) evolution, variable flavour number scheme with
! mc = sqrt(2), mb = 4.5, mt = 175 GeV (the Les Houches benchmark values).
! Used by ci/benchmark-pdf.sh, so that the CI tests do not depend on
! downloading a PDF set from the LHAPDF server.
!
! Usage: make_benchmark_pdf <basename>
program make_benchmark_pdf
  use hoppet
  use pdfs_for_benchmarks
  implicit none
  integer,  parameter :: nloop = 3
  real(dp), parameter :: Q0 = sqrt(2.0_dp), asQ0 = 0.35_dp
  real(dp), parameter :: mc = sqrt(2.0_dp), mb = 4.5_dp, mt = 175.0_dp
  character(len=512) :: basename
  if (command_argument_count() /= 1) error stop "usage: make_benchmark_pdf <basename>"
  call get_command_argument(1, basename)
  ! same grid parameters as the hoppet_N2LO_dy_0.05 LHAPDF sets of hoppet's benchmarking
  call hoppetStartExtended(12.0_dp, 0.05_dp, 1.0_dp, 28000.0_dp, 0.0125_dp, nloop, -6, factscheme_MSbar)
  call hoppetSetVFN(mc, mb, mt)
  call hoppetEvolve(asQ0, Q0, nloop, 1.0_dp, benchmark_pdf_unpolarized_lha, Q0)
  call hoppetWriteLHAPDFGrid(trim(basename), 0)
end program make_benchmark_pdf
