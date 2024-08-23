!----------------------------------------------------------------------
! A module to define all parameters of the run, which can be read from
! command line arguments
module parameters
  use sub_defs_io
  use integration
  use types
  implicit none

  private
  real(dp), parameter, public :: gfermi = 1.16638e-5_dp
  real(dp), parameter, public :: gev2pb = 389379660.0_dp
  real(dp), parameter, public :: eps    = 1.0e-14_dp
  real(dp), parameter, public :: scales_mur(1:7) = &
       & (/1.0_dp, 2.0_dp, 0.5_dp, 1.0_dp, 1.0_dp, 2.0_dp, 0.5_dp/)
  real(dp), parameter, public :: scales_muf(1:7) = &
       & (/1.0_dp, 2.0_dp, 0.5_dp, 2.0_dp, 0.5_dp, 1.0_dp, 1.0_dp/)
  real(dp), public :: xmuf, xmur, Qmin
  real(dp), public :: mh, mh_sq, hwidth
  real(dp), public :: sin_thw, mw, mz, w_width, z_width
  real(dp), public :: v_H, lambda_HHH, lambdafact, cVVHHfact, cVVHfact
  real(dp), public :: sqrts, S, pi, Q0_cut_sq
  real(dp), public :: hmasswindow, toyas
  integer,  public :: order_min, order_max, scale_choice, scale_choice_hoppet, iseed
  integer,  public :: nflav, ipdf, itmx1, itmx2, ncall1, ncall2
  logical,  public :: higgs_use_BW, higgsfixwdth, toypdf
  character * 4, public :: seedstr
  character(len=50), public :: pdfname
  integer, public :: nmempdf
  logical, public :: pdfuncert, scaleuncert3, scaleuncert7, alphasuncert
  real(dp), public :: toy_Q0, test_Q0, muR_PDF
  real(dp), public :: dy, dlnlnQ, minQval, maxQval, ymax
  integer, public :: nloop, order
  public :: set_parameters

contains

  ! set all parameters from input card or command line arguments
  subroutine set_parameters ()
    real(dp) :: ebeam1, ebeam2
    integer  :: idum
    common/ranno/idum
    character * 6 WHCPRG
    common/cWHCPRG/WHCPRG
    mw           = dble_val_opt("-mw",80.379_dp)
    mz           = dble_val_opt("-mz",91.1876_dp)
    nflav        = int_val_opt ("-nf",5)
    mh           = dble_val_opt("-mh",125.0_dp)
    w_width      = dble_val_opt("-wwidth",2.085_dp)
    z_width      = dble_val_opt("-zwidth",2.4952_dp)
    hwidth       = dble_val_opt("-hwidth",0.004029643852284941_dp)
    order_min    = int_val_opt ('-order-min',1)
    order_max    = int_val_opt ('-order-max',4)
    ! if "-lo/-nlo/-nnlo/-n3lo" command is given, overwrite order_min and order_max accordingly
    if (log_val_opt("-lo")) then
       order_min = 1
       order_max = 1
    elseif (log_val_opt("-nlo")) then
       order_min = 1
       order_max = 2
    elseif (log_val_opt("-nnlo")) then
       order_min = 1
       order_max = 3
    elseif (log_val_opt("-n3lo")) then
       order_min = 1
       order_max = 4
    endif    
    sqrts        = dble_val_opt("-sqrts",13600.0_dp)
    scale_choice = int_val_opt ('-scale-choice',3)
    readin       = log_val_opt ("-readingrid")
    higgs_use_BW = log_val_opt ("-higgsbreitwigner")
    hmasswindow  = dble_val_opt("-higgsmasswindow",30.0_dp)
    xmuf         = dble_val_opt("-xmuf",1.0_dp)
    xmur         = dble_val_opt("-xmur",1.0_dp)
    pdfname      = string_val_opt("-pdf", "PDF4LHC21_40")
    nmempdf      = int_val_opt ("-nmempdf",0)
    pdfuncert    = log_val_opt ("-pdfuncert")
    alphasuncert    = log_val_opt ("-alphasuncert")
    scaleuncert3 = log_val_opt ("-3scaleuncert")
    scaleuncert7 = log_val_opt ("-7scaleuncert")
    if(scaleuncert3.and.scaleuncert7) then
       write(*,*) 'WARNING: Have to do either 3 or 7 scales. Cannot do both. Doing none.'
       scaleuncert3 = .false.
       scaleuncert7 = .false.
    endif
    ncall1       = int_val_opt ("-ncall1",100000)
    ncall2       = int_val_opt ("-ncall2",100000)
    itmx1        = int_val_opt ("-itmx1",3)
    itmx2        = int_val_opt ("-itmx2",3)
    iseed        = int_val_opt ("-iseed",10)
    ! for debugging only:
    toypdf       = log_val_opt("-toy")  ! use a toy PDF
    toy_Q0       = dble_val_opt("-toy-Q0",-1.0_dp) ! Q0 where the initial condition of toy PDF is set 
    toyas        = dble_val_opt("-toy-alphas",-1.0_dp) ! toy alphas value at specified scale
    test_Q0      = dble_val_opt("-dglap-Q0",-1.0_dp) ! Q0 for hoppet DGLAP evolution
    muR_PDF      = dble_val_opt("-muR_PDF",1.0_dp) ! for debugging only

    ! various setup
    write(seedstr,"(I4.4)") iseed
    idum = -iseed
    Qmin = 1.0_dp
    pi   = 4.0_dp*atan(1.0_dp)
    Q0_cut_sq = 4.0_dp
    mh_sq = mh**2
    S = sqrts**2
    ! compute couplings needed for dihiggs
    lambdafact = dble_val_opt("-lambdafact",1.0_dp)
    cVVHHfact  = dble_val_opt("-cVVHHfact",1.0_dp)
    cVVHfact   = dble_val_opt("-cVVHfact",1.0_dp)
    v_H = 1.0_dp / (sqrt(sqrt(2.0_dp)*gfermi)) ! v = 246 GeV
    lambda_HHH = lambdafact * mh_sq / (2.0_dp * v_H) ! SM trilinear Higgs self-coupling
    ! compute sin(\theta_w) from W/Z mass
    sin_thw = 1.0_dp - (mw/mz)**2

    ! For hoppetStartExtended. Could think of putting on commandline...
    ! Streamlined initialization
    ! including  parameters for x-grid
    order = -6 
    ymax  = 16.0_dp
    dy    = 0.05_dp  ! dble_val_opt("-dy",0.1_dp)
    dlnlnQ = dy/4.0_dp
    nloop = 3 
    minQval = min(xmuF*Qmin, Qmin)
    maxQval = max(xmuF*sqrts, sqrts)
    scale_choice_hoppet = min(scale_choice,2)
    if (.not.CheckAllArgsUsed(0)) call exit()
  end subroutine set_parameters
  
end module parameters
