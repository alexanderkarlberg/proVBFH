!------------------------------------------------------------
!------------------------------------------------------------
!------------------------------------------------------------
!     
!     Main module for calculation of inclusive VBFH 
!     
!------------------------------------------------------------
!------------------------------------------------------------

module incl_vbfh
  use hoppet_v1, EvolvePDF_hoppet => EvolvePDF, InitPDF_hoppet => InitPDF
  use dummy_pdfs
  use matrix_element
  use incl_parameters
  use phase_space
  use integration
  use types
  implicit none

contains
  !------------------------------------------------------------
  ! inclusive_init sets up the structure functions and PDFs
  subroutine inclusive_init()
    ! set up all constants and parameters from input card
    call set_parameters()

    ! initialize PDF set
    call PDFSET('DEFAULT', dble(ipdf))
    call getQ2min(0,Qmin)
    Qmin = sqrt(Qmin)

    ! !! === DEBUGGING ONLY
    ! ! the following is for the case where
    ! ! we want to use a toy PDF and a toy alphas
    ! if (toypdf) then
    !    write(6,*) 'WARNING: Setting toy pdf values to non-default'
    !    toy_Q0 = 100.0_dp
    !    toy_alphas_Q0 = toyas
    !    Qmin = 1d0
    ! end if
    ! !! ==== END DEBUGGING

    ! initialise hoppet
    call hoppetStartExtended(ymax,dy,minQval,maxQval,dlnlnQ,nloop,&
         &         order,factscheme_MSbar)
    call StartStrFct(sqrts, order_max, nflav, xmur, xmuf,&
         & scale_choice, mh, .true., Qmin)
    call read_PDF()
    call InitStrFct(order_max, separate_orders = .true.)

    ! !! === DEBUGGING ONLY
    ! ! write the structure functions to file, for debugging purposes
    ! call debugging(100.0_dp, 1,1)
    ! stop
    ! !! === END DEBUGGING
    
  end subroutine inclusive_init
  
  subroutine run_inclusive()
    use sub_defs_io
    integer, parameter :: ndim = 7 ! dimension for vegas integration
    real(dp) :: integ, error_int, proba
    real(dp) :: sigma_tot, error_tot, region(1:14)
    logical  :: writeout_tmp
    character(len=30) :: analysis_name, histname
    character * 6 WHCPRG
    common/cWHCPRG/WHCPRG
    integer ilast
    common/last_integ/ilast
    character * 20 pwgprefix
    integer lprefix
    common/cpwgprefix/pwgprefix,lprefix

    call inclusive_init()
    
    ! set beams for phase space generation
    call set_beams(sqrts)

    region(1:ndim)        = 0.0_dp
    region(ndim+1:2*ndim) = 1.0_dp
    sigma_tot = 0.0_dp
    error_tot = 0.0_dp

    outgridfile='grids'//trim(seedstr)//'.dat'
    outgridtopfile='grids'//trim(seedstr)//'.top'
    ilast=0
    
    ! vegas warmup call
    if(.not.readin) then
       ! Skip grid generation if grid is being read from file
       writeout_tmp=writeout
       writeout   = .true.
       fill_plots = .false.
       call vegas(region,ndim,dsigma,0,ncall1,itmx1,0,integ,error_int,proba)
       write(6,*) integ
       writeout=writeout_tmp
    endif
    
    fill_plots = do_analysis_call
    ! intialize the powheg analysis if necessary
    if (fill_plots) call init_hist 

    ! vegas main call
    ingridfile = 'grids'//trim(seedstr)//'.dat'
    call vegas(region,ndim,dsigma,1,ncall2,itmx2,0,integ,error_int,proba)

    ! add integral to the total cross section
    sigma_tot = sigma_tot + integ
    error_tot = error_tot + error_int**2

    !------------------------------------------------------------
    ! output results
    
    ! write histograms to file if they have been requested
    if (fill_plots) then
       call pwhgsetout
       histname=pwgprefix(1:lprefix)//'-LO'//trim(seedstr)
       call pwhgtopout(histname)
    endif

    ! print total cross section and error into file  
    ! construct name of output file
    analysis_name='xsct'
    if (order_max.eq.1) then
       analysis_name='xsct-lo'//trim(seedstr)//'.dat'
    else if (order_max.eq.2) then
       analysis_name='xsct-nlo'//trim(seedstr)//'.dat'
    else if (order_max.eq.3) then
       analysis_name='xsct-nnlo'//trim(seedstr)//'.dat'
    else if (order_max.eq.4) then
       analysis_name='xsct-n3lo'//trim(seedstr)//'.dat'
    endif

    OPEN(UNIT=11, FILE=trim(analysis_name), ACTION="write")
    write(11,'(a)',advance='no') '# '
    call time_stamp(11)
    write(11,'(a)') '# '//trim(command_line())
    if (order_max.eq.1) then 
       write(11,'(a,es13.6,a,es13.6,a)') '# Total LO cross-section'
       write(11,'(es13.6,a,es13.6,a)') sigma_tot, ' +/- ', sqrt(error_tot),' (pb)'
    else if (order_max.eq.2) then 
       write(11,'(a,es13.6,a,es13.6,a)') '# Total NLO cross-section'
       write(11,'(es13.6,a,es13.6,a)') sigma_tot, ' +/- ', sqrt(error_tot),' (pb)'
    else if (order_max.eq.3) then 
       write(11,'(a,es13.6,a,es13.6,a)') '# Total NNLO cross-section'
       write(11,'(es13.6,a,es13.6,a)') sigma_tot, ' +/- ', sqrt(error_tot),' (pb)'
    else if (order_max.eq.4) then 
       write(11,'(a)') '# Total N3LO cross-section'
       write(11,'(es13.6,a,es13.6,a)') sigma_tot, ' +/- ', sqrt(error_tot),' (pb)'
    endif
    write(6,'(a)')

  end subroutine run_inclusive


  !----------------------------------------------------------------------
  ! Debugging routine to output structure functions at
  ! (xmur,xmuf) = (1,1),(1,4),(4,1),(4,3)
  subroutine debugging (Qtest, sc_min, sc_max)
    use incl_parameters
    use hoppet_v1
    real(dp), intent(in) :: Qtest
    integer,  intent(in) :: sc_min, sc_max
    integer  :: ihopf1,ihopf2,ihopf3, ix, ny, sc
    real(dp) :: xmufvals(1:4), xmurvals(1:4), ytest_max
    character*20 :: dir
    character*3  :: sxmur,sxmuf
    character*6  :: str_scale_ch
    
    dir = "str_fct/"
    ytest_max = log(1e5)
    ny = 100
    
    !scale choices
    xmufvals(1) = 1.0d0
    xmurvals(1) = 1.0d0

    xmufvals(2) = 4.0d0
    xmurvals(2) = 1.0d0

    xmufvals(3) = 1.0d0
    xmurvals(3) = 4.0d0

    xmufvals(4) = 3.0d0
    xmurvals(4) = 4.0d0

    do sc = sc_min, sc_max
      scale_choice = sc
      if (scale_choice.eq.0) then
         str_scale_ch = '_scMH'
      elseif (scale_choice.eq.1) then
         str_scale_ch = '_scQ'
      elseif (scale_choice.eq.2) then
         str_scale_ch = '_scMIX'
      endif
      ! loop over scale variations in xmuf, xmur
      do ix = 1, 4
         xmur = xmurvals(ix)
         xmuf = xmufvals(ix)
         write(sxmur,"(F3.1)") xmur
         write(sxmuf,"(F3.1)") xmuf
         ihopf1 = 21
         OPEN(UNIT=ihopf1, &
              &        FILE=trim(dir)//"hoppet-F1WpmZ_mur"//sxmur//"_muf"// &
              &             sxmuf//trim(str_scale_ch)//".dat", &
              &        ACTION="write")
         ihopf2 = 22
         OPEN(UNIT=ihopf2, & 
              &        FILE=trim(dir)//"hoppet-F2WpmZ_mur"//sxmur//"_muf"// &
              &             sxmuf//trim(str_scale_ch)//".dat", &
              &        ACTION="write")
         ihopf3 = 23
         OPEN(UNIT=ihopf3, &
              &        FILE=trim(dir)//"hoppet-F3WpmZ_mur"//sxmur//"_muf"// &
              &             sxmuf//trim(str_scale_ch)//".dat", &
              &        ACTION="write")

         !call set_structure_functions_up_to(order_max)
         call write_f1(ihopf1, Qtest, ytest_max, ny)
         call write_f2(ihopf2, Qtest, ytest_max, ny)
         call write_f3(ihopf3, Qtest, ytest_max, ny)

         close(ihopf1)
         close(ihopf2)
         close(ihopf3)
      enddo
   enddo
  end subroutine debugging

  
  !----------------------------------------------------------------------
  ! fill the streamlined interface PDF table (possibly using hoppet's
  ! evolution)
  subroutine read_PDF(toyQ0, dglapQ0, xR_PDF)
    use streamlined_interface
    real(dp), optional :: toyQ0, dglapQ0, xR_PDF
    real(dp), external :: alphasPDF
    interface
       subroutine EvolvePDF(x,Q,res)
         use types; implicit none
         real(dp), intent(in)  :: x,Q
         real(dp), intent(out) :: res(*)
       end subroutine EvolvePDF
    end interface
    !----------------
    real(dp) :: toy_Q0, Q0pdf, xmuR_PDF
    real(dp) :: res_lhapdf(-6:6), x, Q
    real(dp) :: res_hoppet(-6:6)
    real(dp) :: toy_pdf_at_Q0(0:grid%ny,ncompmin:ncompmax)
    real(dp), parameter :: mz = 91.2_dp
    real(dp) :: pdf_at_Q0(0:grid%ny,ncompmin:ncompmax)

    toy_Q0       = -one
    Q0pdf        = -one
    xmuR_PDF     = one  
    if(present(toyQ0)) toy_Q0=toyQ0
    if(present(dglapQ0)) Q0pdf=dglapQ0
    if(present(xR_PDF)) xmuR_PDF=xR_PDF

    if (toy_Q0 > zero) then
       write(6,*) "WARNING: Using toy PDF"
       toy_pdf_at_Q0 = unpolarized_dummy_pdf(xValues(grid))
       call InitRunningCoupling(toy_coupling, alfas=toy_alphas_Q0, &
            &                   nloop = 3, Q = toy_Q0, fixnf=nf_int)
       call EvolvePdfTable(tables(0), toy_Q0, toy_pdf_at_Q0, dh, toy_coupling, nloop=3)
    elseif (Q0pdf > zero) then

       write(6,*) "WARNING: Using internal HOPPET DGLAP evolution"
       call InitPDF_LHAPDF(grid, pdf_at_Q0, EvolvePDF, Q0pdf)
       call InitRunningCoupling(coupling, alphasPDF(MZ) , MZ , 4,&
            & -1000000045, sf_quark_masses(4:6), .true.)
       call EvolvePdfTable(tables(0), Q0pdf, pdf_at_Q0, dh, coupling, &
            &  muR_Q=xmuR_PDF, nloop=3)

    else
       ! InitRunningCoupling has to be called for the HOPPET coupling to be initialised 
       ! Default is to ask for 4 loop running and threshold corrections at quark masses.  
       call InitRunningCoupling(coupling, alphasPDF(MZ) , MZ , 4,&
            & -1000000045, sf_quark_masses(4:6), .true.)
       ! fixnf can be set to a positive number for
       ! fixed nf. -1000000045 gives variable nf
       ! and threshold corrections at quarkmasses.
       call hoppetAssign(EvolvePDF)
    end if

    ! quickly test that we have read in the PDFs correctly
    write(6,*) "Quick test that PDFs have been read in correctly"
    x = 0.08_dp
    Q = 17.0_dp
    call EvolvePDF(x, Q, res_lhapdf)
    call EvalPdfTable_xQ(tables(0), x, Q, res_hoppet)
    write(6,*) 'lhapdf: ', res_lhapdf
    write(6,*) 'hoppet: ', res_hoppet
  end subroutine read_PDF

end module incl_vbfh
