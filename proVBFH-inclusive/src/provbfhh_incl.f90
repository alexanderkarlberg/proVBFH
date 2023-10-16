!------------------------------------------------------------
!------------------------------------------------------------
!---- proVBFH-inclusive, Version 2.0.2 ----------------------
!     
!  Program for the computation of inclusive VBFH pair production
!  up to N3LO.
!
!  Written by
!  Matteo Cacciari cacciari@lpthe.jussieu.fr 
!  Frederic Dreyer, frederic.dreyer@physics.ox.ac.uk
!  Alexander Karlberg, alexander.karlberg@cern.ch
!  Gavin Salam, gavin.salam@cern.ch
!  Giulia Zanderighi, g.zanderighi1@physics.ox.ac.uk
!
!  based on  
!  [arXiv:1811.07918], submitted to Phys.Rev.Lett.
!  [arXiv:1811.07906], submitted to Phys.Rev.D
!  Phys.Rev.Lett. 115 (2015) no.8, 082002 [arXiv:1506.02660]
!  Phys.Rev.Lett. 117 (2016) no.7, 072001 [arXiv:1606.00840]
!
!  as well as the phase space from VBFNLO
!
!  First version: 10/02/2015
!  Last edited: 10/10/2022
!     
!------------------------------------------------------------
!------------------------------------------------------------

program provbfh_incl
  use hoppet_v1, EvolvePDF_hoppet => EvolvePDF, InitPDF_hoppet => InitPDF
  use sub_defs_io
  use dummy_pdfs
  use matrix_element_dihiggs
  use parameters
  use phase_space_dihiggs
  use integration
  use types
  implicit none
  integer, parameter :: ndim = 11 ! dimension for vegas integration
  integer, parameter :: nmempdfMAX = 200 ! max number of pdfs
  integer  :: nmempdf_start, nmempdf_end, imempdf
  real(dp) :: integ, error_int, proba, tini, tfin
  real(dp) :: sigma_tot, error_tot, region(1:2*ndim)
  real(dp) :: res(0:nmempdfMAX), central, errminus, errplus, errsymm, respdf(0:nmempdfMAX),resas(0:nmempdfMAX), central_dummy
  real(dp) :: res_scales(1:7),maxscale,minscale
  integer :: iscales,Nscales
  character * 30 :: analysis_name
  character * 6 WHCPRG
  common/cWHCPRG/WHCPRG
  integer ilast
  common/last_integ/ilast
  integer saveseed,idum
  COMMON /ranno/ idum
  call cpu_time(tini)

  ! set up all constants and parameters from command line arguments
  call set_parameters()

  ! initialize PDF set
  ! call PDFSET('DEFAULT', dble(ipdf))
  call initPDFSetByName(pdfname)
  
  if (pdfuncert) then
     nmempdf_start = 0
     call numberPDF(nmempdf_end)
  else
     nmempdf_start = nmempdf
     nmempdf_end   = nmempdf
  endif

  if (nmempdf_end .gt. nmempdfMAX) stop "ERROR: increase nmempdfMAX"

  if(.not.scaleuncert3.and..not.scaleuncert7) then
     Nscales = 1
  elseif(scaleuncert3) then
     Nscales = 3
  elseif(scaleuncert7) then
     Nscales = 7
  endif

  do imempdf = nmempdf_start,nmempdf_end
     write(6,*) "PDF member:",imempdf
     call InitPDF(imempdf)
     call getQ2min(0,Qmin)
     Qmin = sqrt(Qmin)

     ! !! === DEBUGGING ONLY
     ! ! the following is for the case where
     ! ! we want to use a toy PDF and a toy alphas
     ! if (toypdf.gt.0) then
     !    write(6,*) 'WARNING: Setting toy pdf values to non-default'
     !    toy_Q0 = 100.0_dp
     !    toy_alphas_Q0 = toyas
     !    Qmin = 1d0
     ! end if
     ! !! ==== END DEBUGGING
     do iscales = 1,Nscales
        ! initialise hoppet
        xmur = scales_mur(iscales)
        xmuf = scales_muf(iscales)
        call hoppetStartExtended(ymax,dy,minQval,maxQval,dlnlnQ,nloop,&
             &         order,factscheme_MSbar)
        call StartStrFct(sqrts, order_max, nflav, xmur, &
             & xmuf, scale_choice, mh, .true., Qmin, mw, mz)
        call read_PDF(toy_Q0, test_Q0, mur_PDF)
        call InitStrFct(order_max, separate_orders = .true.)
        
        ! !! === DEBUGGING ONLY
        ! ! write the structure functions to file, for debugging purposes
        ! call debugging(100.0_dp, 1,1)
        ! stop
        ! !! === END DEBUGGING
        
        ! set beams for phase space generation
        call set_beams(sqrts)
        
        region(1:ndim)        = 0.0_dp
        region(ndim+1:2*ndim) = 1.0_dp
        sigma_tot = 0.0_dp
        error_tot = 0.0_dp
        
        outgridfile='grids_'//seedstr//'.dat'
        outgridtopfile='grids_'//seedstr//'.top'
        ilast=0
        
        ! vegas warmup call
        if(.not.readin) then 
           ! Skip grid generation if grid is being read from file
           writeout=.true.
           call vegas(region,ndim,dsigma,0,ncall1,itmx1,0,integ,error_int,proba)
           writeout=.false.
           ! set random seed to current idum value
           saveseed = idum
        elseif (imempdf.eq.nmempdf_start) then
           ! if reading in grids from first loop iteration, make sure
           ! saveseed is initialized to correct value
           saveseed = iseed
        endif
        ! vegas main call
        ingridfile ='grids_'//seedstr//'.dat'
        ! set random seed to saved value 
        idum     = -saveseed
        call vegas(region,ndim,dsigma,1,ncall2,itmx2,0,integ,error_int,proba)
        readin = .true.
        ! add integral to the total cross section
        sigma_tot = sigma_tot + integ
        error_tot = error_tot + error_int**2

        ! convert pb => fb
        sigma_tot = sigma_tot*1000.0_dp
        error_tot = error_tot*1000000.0_dp        
        
        res(imempdf) = sigma_tot
        res_scales(iscales) = sigma_tot
        ! end loop over pdfs
     enddo
     if(imempdf.eq.nmempdf_start) then ! First PDF, this is where we compute scale uncertainties
       maxscale = maxval(res_scales(1:Nscales)) 
       minscale = minval(res_scales(1:Nscales))
       Nscales = 1
       res(imempdf) = res_scales(1) ! Copy central scale
    endif
  enddo

  ! print total cross section and error into file  
  ! construct name of output file
  analysis_name='xsct'
  if (order_max.eq.1) then
     analysis_name="xsct_lo_seed"//seedstr//".dat"
  else if (order_max.eq.2) then
     analysis_name="xsct_nlo_seed"//seedstr//".dat"
  else if (order_max.eq.3) then
     analysis_name="xsct_nnlo_seed"//seedstr//".dat"
  else if (order_max.eq.4) then
     analysis_name="xsct_n3lo_seed"//seedstr//".dat"
  endif
  
  !------------------------------------------------------------
  ! output results
  OPEN(UNIT=11, FILE=analysis_name, ACTION="write")
  write(11,'(a)',advance='no') '# '
  call time_stamp(11)
  write(11,'(a)') '# '//trim(command_line())
  if (order_max.eq.1) then 
     write(11,'(a,es13.6,a,es13.6,a)') '# Total LO cross-section (fb)'
  else if (order_max.eq.2) then 
     write(11,'(a,es13.6,a,es13.6,a)') '# Total NLO cross-section (fb)'
  else if (order_max.eq.3) then 
     write(11,'(a,es13.6,a,es13.6,a)') '# Total NNLO cross-section (fb)'
  else if (order_max.eq.4) then 
     write(11,'(a)') '# Total N3LO cross-section (fb)'
  endif

  if (nmempdf_start.eq.nmempdf_end.and..not.scaleuncert3.and..not.scaleuncert7) then
     write(11,'(a)') '# central     MC_error'
     write(11,'(es13.6,es13.6)') sigma_tot, sqrt(error_tot)
  elseif(nmempdf_start.eq.nmempdf_end) then
     central = res(0)
     write(11,'(a)') '# central     max          min          MC_error'
     write(11,'(7(es13.6))') central, maxscale, minscale,sqrt(error_tot)
     write(11,'(a)') ''
     write(11,'(a)') ''
     write(11,'(a)') '# Summary:'
     write(11,'(a,f10.5,a)') '# sigma =', central,' fb'
     write(11,'(a,f9.5,a,a,f9.3,a)') '# QCD scale uncertainty (+) =', maxscale-central, ' fb', &
          & ' (', ((maxscale-central)/central)*100.0_dp, ' %)'
     write(11,'(a,f9.5,a,a,f9.3,a)') '# QCD scale uncertainty (-) =', minscale-central, ' fb', &
          & ' (', ((minscale-central)/central)*100.0_dp, ' %)'
     write(11,'(a,f10.3,a)') '# MC integration uncertainty =', sqrt(error_tot)/central*100.0_dp, ' %'
  elseif(.not.scaleuncert3.and..not.scaleuncert7) then
     call getpdfuncertainty(res(nmempdf_start:nmempdf_end),central,errplus,errminus,errsymm)
     write(11,'(a)') '# central     max          min          MC_error'
     write(11,'(7(es13.6))') central, maxscale, minscale,sqrt(error_tot)
     write(11,'(a)') ''
     write(11,'(a)') ''
     write(11,'(a)') '# Summary:'
     write(11,'(a,f10.5,a)') '# sigma =', central,' fb'
     write(11,'(a,f10.3,a)') '# MC integration uncertainty =', sqrt(error_tot)/central*100.0_dp, ' %'
     write(11,'(a,f10.3,a)') '# PDF symmetric uncertainty* =', errsymm/central*100.0_dp, ' %'
     if(alphasuncert) then
        central_dummy = central
        respdf = central
        resas = central
        respdf(0:nmempdf_end-2) = res(0:nmempdf_end-2)
        resas(nmempdf_end-1:nmempdf_end) = res(nmempdf_end-1:nmempdf_end)
        call getpdfuncertainty(respdf(nmempdf_start:nmempdf_end),central_dummy,errplus,errminus,errsymm)
        write(11,'(a,f10.3,a)') '# Pure PDF symmetric uncertainty =', errsymm/central*100.0_dp, ' %'
        call getpdfuncertainty(resas(nmempdf_start:nmempdf_end),central_dummy,errplus,errminus,errsymm)
        write(11,'(a,f10.3,a)') '# Pure αS symmetric uncertainty =', errsymm/central*100.0_dp, ' %'
     endif
     write(11,'(a)') '# (*PDF uncertainty contains alphas uncertainty if using a'
     write(11,'(a)') '#   PDF set that supports it (eg PDF4LHC15_nnlo_100_pdfas))'
  else
     call getpdfuncertainty(res(nmempdf_start:nmempdf_end),central,errplus,errminus,errsymm)
     write(11,'(a)') '# central     max          min          MC_error     PDF_err_plus PDF_err_min  PDF_err_symm'
     write(11,'(7(es13.6))') central, maxscale, minscale,sqrt(error_tot), errplus, errminus, errsymm
     write(11,'(a)') ''
     write(11,'(a)') ''
     write(11,'(a)') '# Summary:'
     write(11,'(a,f10.5,a)') '# sigma =', central,' fb'
     write(11,'(a,f9.5,a,a,f9.3,a)') '# QCD scale uncertainty (+) =', maxscale-central, ' fb', &
          & ' (', ((maxscale-central)/central)*100.0_dp, ' %)'
     write(11,'(a,f9.5,a,a,f9.3,a)') '# QCD scale uncertainty (-) =', minscale-central, ' fb', &
          & ' (', ((minscale-central)/central)*100.0_dp, ' %)'
     write(11,'(a,f10.3,a)') '# MC integration uncertainty =', sqrt(error_tot)/central*100.0_dp, ' %'
     write(11,'(a,f10.3,a)') '# PDF symmetric uncertainty* =', errsymm/central*100.0_dp, ' %'
     if(alphasuncert) then
        central_dummy = central
        respdf = central
        resas = central
        respdf(0:nmempdf_end-2) = res(0:nmempdf_end-2)
        resas(nmempdf_end-1:nmempdf_end) = res(nmempdf_end-1:nmempdf_end)
        call getpdfuncertainty(respdf(nmempdf_start:nmempdf_end),central_dummy,errplus,errminus,errsymm)
        write(11,'(a,f10.3,a)') '# Pure PDF symmetric uncertainty =', errsymm/central*100.0_dp, ' %'
        call getpdfuncertainty(resas(nmempdf_start:nmempdf_end),central_dummy,errplus,errminus,errsymm)
        write(11,'(a,f10.3,a)') '# Pure αS symmetric uncertainty =', errsymm/central*100.0_dp, ' %'
     endif
     write(11,'(a)') '# (*PDF uncertainty contains alphas uncertainty if using a'
     write(11,'(a)') '#   PDF set that supports it (eg PDF4LHC15_nnlo_100_pdfas))'
  endif

  call cpu_time(tfin)
  write(6,'(a)')
  write(6,'(a,es9.2,a)') '==================== TOTAL TIME : ', &
       &     tfin-tini, ' s.'
  write(6,'(a)')
contains

  !----------------------------------------------------------------------
  ! Debugging routine to output structure functions at
  ! (xmur,xmuf) = (1,1),(1,4),(4,1),(4,3)
  subroutine debugging (Qtest, sc_min, sc_max)
    use parameters
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
            &  muR_Q=muR_PDF, nloop=3)

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

  end program provbfh_incl
