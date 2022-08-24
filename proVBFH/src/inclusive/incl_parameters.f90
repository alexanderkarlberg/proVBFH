!----------------------------------------------------------------------
! A module to define all parameters of the run, which can either be read from
! an input card (powheg.input, vbfnlo.input) or from command line arguments
module incl_parameters
  use integration
  use types
  implicit none

  private
  logical, save, public :: param_initialised = .false.
  real(dp), parameter, public :: gfermi = 1.16637e-5_dp
  real(dp), parameter, public :: gev2pb = 389379660.0_dp
  real(dp), parameter, public :: eps    = 1.0e-14_dp
  real(dp), public :: xmuf, xmur, Qmin
  real(dp), public :: mh, mh_sq, hwidth
  real(dp), public :: sin_thw, mw, mz, w_width, z_width
  real(dp), public :: sqrts, S, pi, Q0_cut_sq
  real(dp), public :: hmasswindow, toyas
  integer,  public :: order_min, order_max, scale_choice, iwhichseed
  integer,  public :: nflav, ipdf, itmx1, itmx2, ncall1, ncall2
  logical,  public :: do_analysis_call, manyseeds,non_fact,oneloop_on,twoloop_on
  logical,  public :: higgs_use_BW, higgsfixwdth, complexpole, toypdf
  character * 5, public :: seedstr
  real(dp), public :: test_Q0, muR_PDF
  integer,  public :: niter
  public :: set_parameters

contains

  ! set all parameters from input card or command line arguments
  subroutine set_parameters ()
    real(dp) :: powheginput, vbfnloinput
    real(dp) :: ebeam1, ebeam2
    integer  :: idum, iseed, iun, numseeds, ios, j
    common/ranno/idum
    character * 6 WHCPRG
    common/cWHCPRG/WHCPRG
    character * 20 pwgprefix
    integer lprefix
    common/cpwgprefix/pwgprefix,lprefix

    ! read all parameters from input cards
    mw=vbfnloinput('#WMASS')
    mz=vbfnloinput('#ZMASS')
    w_width=vbfnloinput('#WWIDTH')
    z_width=vbfnloinput('#ZWIDTH')
    nflav=vbfnloinput('#NFLAVOUR')
    mh=vbfnloinput('#HMASS')
    hwidth = vbfnloinput('#HWIDTH')
    order_min= powheginput('#order_min')
    if(order_min.lt.1.) order_min = 1
    order_max=powheginput('#qcd_order')
    if(order_min.gt.order_max) stop 'order_min > order_max'
    ebeam1=powheginput('#ebeam1')
    ebeam2=powheginput("#ebeam2")
    higgs_use_BW=.false.
    if(powheginput('#higgsbreitwigner').eq.1) then
       higgs_use_BW=.true.
    endif
    hmasswindow = 30.0_dp
    if(powheginput('#higgsmasswindow').gt.0d0) then
       hmasswindow = powheginput('#higgsmasswindow')
    endif
    sqrts=ebeam1+ebeam2
    scale_choice = 0 ! MH
    if(powheginput('#runningscales').eq.1d0) then
       scale_choice=3 ! Scale of 1506.02660
    elseif(powheginput('#runningscales').eq.2d0) then
       scale_choice=1 ! Q1, Q2 scale
    endif
    readin=.false.
    writeout=.false.
    if(powheginput('#readingrid').eq.1) then
       readin=.true.
    endif
    if(powheginput('#writeoutgrid').eq.1) then
       writeout=.true.
    endif
    if(powheginput('#nonfact').eq.1) then
       non_fact = .true.
       niter = 100
       if(powheginput('#niter').gt.0d0) then
          niter = int(powheginput('#niter'))
       endif
       if(order_max.ne.1) then
          print*, 'In order to compute non factorisable corrections qcd_order &
               & has to be set to 1 in the input card'
          stop
       endif
       print*, 'Non factorisable contributions are being computed on their own. &
            & All factorisable contributions have been turned off.'
       oneloop_on = .true.
       if(powheginput("#oneloop_off").eq.1d0) oneloop_on = .false.
       twoloop_on = .true.
       if(powheginput("#twoloop_off").eq.1d0) twoloop_on = .false.
    else
       non_fact = .false.
    endif
    xmuf=powheginput('#facscfact')
    xmur=powheginput('#renscfact')
    ipdf=powheginput('#lhans1')
    ncall1=powheginput('#ncall1')
    ncall2=powheginput('#ncall2') 
    itmx1=powheginput('#itmx1')
    itmx2=powheginput('#itmx2')
    iseed=powheginput('#iseed')
    seedstr=''
    !   if manyseeds flag is on, ignore iseed and read from pwgseeds.dat instead
    manyseeds=.false.
    iwhichseed=-1
    if(powheginput('#manyseeds').eq.1) then
       manyseeds=.true.
       iun=6
       call newunit(iun)
       open(unit=iun,status='old',iostat=ios, file=pwgprefix(1:lprefix)//'seeds.dat')
       if(ios.ne.0) then
          write(*,*) 'option manyseeds required but '
          write(*,*) 'file ',pwgprefix(1:lprefix)//'seeds.dat not found'
          call exit(-1)
       endif
       do j=1,1000000
          read(iun,*,iostat=ios) iseed
          if(ios.ne.0) exit
       enddo
       numseeds=j-1
       write(*,*) 'enter which seed'
       read(*,*) iwhichseed
       if(iwhichseed.gt.numseeds) then
          write(*,*) 'no more than',numseeds,'seeds in',pwgprefix(1:lprefix)//'seeds.dat'
          call exit(-1)
       endif
       rewind(iun)
       do j=1,iwhichseed
          read(iun,*) iseed
       enddo
       close(iun)
       write(seedstr,"(I4.4)") iseed
       seedstr='-'//seedstr
    endif

    ! decide whether we need to call the analysis based on inputs
    ! call analysis only at LO, when inclusive_only is set to 0, i.e.:
    !  inclusive_only 0
    !  qcd_order      1
    if (powheginput('#qcd_order').eq.1d0 &
         & .and. powheginput('#inclusive_only').ne.1d0) then
       do_analysis_call = .true.
    else
       do_analysis_call = .false.
    endif

    do_analysis_call = .true.
    
    ! debugging and advanced options
    toypdf=.false.
    if(powheginput('#toypdf').eq.1) then
       toypdf=.true.
    endif
    higgsfixwdth=powheginput("#higgsfixedwidth").gt.0
    complexpole=powheginput("#complexpolescheme").gt.0
    toyas=powheginput('#toyalphas')
    test_Q0=powheginput('#test_Q0')
    muR_PDF=powheginput('#muR_PDF')
    
    ! now set up a few other things
    idum = -iseed
    Qmin = 1.0_dp
    pi   = 4.0_dp*atan(1.0_dp)
    Q0_cut_sq = 4.0_dp
    mh_sq = mh**2
    S = sqrts**2
    ! compute sin(\theta_w) from W/Z mass
    sin_thw = 1.0_dp - (mw/mz)**2

    param_initialised = .true.
  end subroutine set_parameters

end module incl_parameters
