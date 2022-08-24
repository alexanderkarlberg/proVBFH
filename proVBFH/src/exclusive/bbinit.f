      subroutine bbinit
      implicit none
      integer iret
      real * 8 powheginput
      integer parallelstages
      external powheginput
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_flg.h'
      include 'cgengrids.h'
      include 'pwhg_rad.h'
      include 'pwhg_rnd.h'
      include 'phspcuts.h'
      real * 8 xx(ndiminteg)
      integer mcalls,icalls,j
      integer btilde,sigremnant
      external btilde,sigremnant
c parallelstages:
c 1   prepare the importance sampling grids
c 2   prepare the upper bounding envelopes for the
c     generation of the b_tilde function
c 3   prepare the upper bound for the generation of radiation
c     4   generate events
      phspcuts = .false.
      if(powheginput('#phspcuts').eq.1d0) phspcuts = .true.
      parallelstages =  powheginput('#parallelstage')
      if(flg_newweight .and. parallelstages .gt. 0) then
         write(*,*) ' Since we are running in reweighting mode '
         write(*,*) ' we set parallelstages to 4 '
         parallelstages = 4
      endif
      if(parallelstages.gt.0.and.rnd_cwhichseed.eq.'none') then
         write(*,*) ' with parallelstage also manyseeds '
         write(*,*) ' must be set'
         call pwhg_exit(-1)
      endif
      if(parallelstages.le.0) then
c Do all stages in one go if needed;
c first look for <prefix>xgrid.dat file      
         call loadxgrids(iret)
         if(iret.ne.0) then
c if not there look for pwggridinfo-* files from parallel runs
            call loadlatestxgridinfo(iret)
            if(iret.ne.0) then
c if not there generate the xgrid; the argument 0 means:
c generate *xgrid.dat file
               call bbinitxgrids(0)
            endif
         endif
c make sure the ifold arrays are read in
         do j=1,ndiminteg
            ifold(j)=1
            ifoldrm(j)=1
         enddo
c Override real integration parameters with powheg.input values         
         ifold(ndiminteg-2) = powheginput("foldcsi")
         ifold(ndiminteg-1) = powheginput("foldy")
         ifold(ndiminteg)   = powheginput("foldphi")

         if(flg_storemintupb) then
            call loadgrids(iret,xgrid,ymax,ymaxrat,xgridrm,ymaxrm,
     1           ymaxratrm,ifold,ifoldrm,'fullgrid')
         else
            call loadgrids(iret,xgrid,ymax,ymaxrat,xgridrm,ymaxrm,
     2           ymaxratrm,ifold,ifoldrm,'grid')
         endif
         if(iret.eq.0) then
            write(*,*) 'upper bound grids successfully loaded'
            write(*,*) 'btilde pos.   weights:', rad_totposbtl,' +-',
     1           rad_etotposbtl
            write(*,*) 'btilde |neg.| weights:', rad_totnegbtl,' +-',
     1           rad_etotnegbtl
            write(*,*)
     1           'btilde total (pos.-|neg.|):', rad_totbtl,' +-',
     2           rad_etotbtl
         else
            call bbinitgrids
            if(flg_storemintupb) then
               call loadmintupb(ndiminteg,'btildeupb',ymax,ymaxrat)
               write(*,*) ' Upper bounding envelope for btilde computed'
               write(*,*)
     1              ' Efficiency for btilde generation is printed above'
               if((flg_withreg.or.flg_withdamp)
     1              .and..not.flg_bornonly) then
                  call loadmintupb(ndiminteg,'remnupb',ymaxrm,ymaxratrm)
                  write(*,*)
     1                 ' Upper bounding envelope for remnants computed'
                  write(*,*)
     1             ' Efficiency for remnant generation is printed above'
               endif
               call storegrids(xgrid,ymax,ymaxrat,xgridrm,ymaxrm,
     1              ymaxratrm,ifold,ifoldrm,-1,-1,'fullgrid')
            endif
         endif


c initialize gen; the array xmmm is set up at this stage.
         call gen(btilde,ndiminteg,xgrid,ymax,ymaxrat,xmmm,ifold,0,
     1    mcalls,icalls,xx)
         
         if (.not.flg_LOevents) then
c load or compute normalization of upper bounding function for radiation
c (iret ignored)
            call do_maxrat(mcalls,icalls,-1,iret)
         endif
c print statistics
         call gen(btilde,ndiminteg,xgrid,ymax,ymaxrat,xmmm,ifold,3,
     1        mcalls,icalls,xx)
         if(xx(1).gt.0)
     1        write(*,*) 'POWHEG: efficiency in the generation'
     2        //' of the Born variables =',xx(1)
         if((flg_withreg.or.flg_withdamp).and..not.flg_bornonly) then
c initialize gen for remnants
            call gen(sigremnant,ndiminteg,xgridrm,ymaxrm,ymaxratrm,
     1           xmmmrm,ifoldrm,0,mcalls,icalls,xx)
c     save random number seeds
            if(rad_totrm/rad_totgen.gt.1d-4.and.
     1           powheginput('#skipextratests').lt.0) then
               call randomsave
c     generate few events from remnants, just to determine the generation efficiency
               do j=1,int(min(powheginput('nubound'),10d0))
                  call gen(sigremnant,ndiminteg,xgridrm,ymaxrm,
     1                ymaxratrm,xmmmrm,ifoldrm,1,mcalls,icalls,xx)
               enddo
c     restore  random number seeds
               call randomrestore
c     print statistics
               call gen(sigremnant,ndiminteg,xgridrm,ymaxrm,ymaxratrm,
     1              xmmmrm,ifoldrm,3,mcalls,icalls,xx)
               if(xx(1).gt.0) then
                  write(*,*) 'POWHEG: efficiency in the generation'
     1                 //' of the remnant variables =',xx(1)
               endif
            endif
         endif
         return
      endif
      if(parallelstages.eq.1) then
         call bbinitxgrids(1)
         call pwhg_exit(0)
      endif
! AK DEBUG. Here we load a uniform grid to see how the phasespace
! performs on its own
!      call initxgrid(xgrid,ndiminteg)
!      xacc = 1d0 
!      nhits = 1d0
!      xint = 1d0
!      call storegridinfo('btl-'//trim('xg1'),xgrid,xint,xacc,nhits
!     $     ,ndiminteg)
! AK DEBUG END
      call loadxgrids(iret)
      if(iret.ne.0) then
c if not there look for pwggridinfo-* files from parallel runs
         call loadlatestxgridinfo(iret)
         if(iret.ne.0) then
            write(*,*) ' cannot load xgrid or gridinfo files;'
            write(*,*) ' cannot perform stage 2'
            call pwhg_exit(-1)
         endif
      endif
c make sure the ifold arrays are read in
      do j=1,ndiminteg
         ifold(j)=1
         ifoldrm(j)=1
      enddo
      ifold(ndiminteg-2) = powheginput("foldcsi")
      ifold(ndiminteg-1) = powheginput("foldy")
      ifold(ndiminteg)   = powheginput("foldphi")
      if(parallelstages.eq.2) then
         call bbinitgrids
         call writestat('st2')
         call pwhg_exit(0)
      endif
c in all cases pwggrids files must be loaded, to include information
c on total cross sections components (rad_tot* and rad_etot* variables)
      call loadgrids(iret,xgrid,ymax,ymaxrat,xgridrm,ymaxrm,
     1     ymaxratrm,ifold,ifoldrm,'grid')
      if(iret.ne.0) then
         write(*,*) ' cannot load grid files'
         call pwhg_exit(-1)
      endif

      if(flg_storemintupb) then
c     Set up better upper bounding envelopes for btilde
         write(*,*) ' Loading the upper bounding envelope for btilde:'
         call loadgrids(iret,xgrid,ymax,ymaxrat,xgridrm,ymaxrm,
     1        ymaxratrm,ifold,ifoldrm,'fullgrid')
         if(iret.ne.0) then
c this aborts if no files are found
            call loadmintupb(ndiminteg,'btildeupb',ymax,ymaxrat)
            write(*,*) ' Upper bounding envelope for btilde computed'
            write(*,*)
     1           ' Efficiency for btilde generation is printed above'
            if((flg_withreg.or.flg_withdamp).and..not.flg_bornonly) then
               call loadmintupb(ndiminteg,'remnupb',ymaxrm,ymaxratrm)
               write(*,*)
     1             ' Upper bounding envelope for remnants computed'
               write(*,*)
     1             ' Efficiency for remnant generation is printed above'
            endif
            call storegrids(xgrid,ymax,ymaxrat,xgridrm,ymaxrm,
     1        ymaxratrm,ifold,ifoldrm,-1,-1,'fullgrid')
  
         endif
      endif

c initialize gen; the array xmmm is set up at this stage.
      call gen(btilde,ndiminteg,xgrid,ymax,ymaxrat,xmmm,ifold,0,
     1     mcalls,icalls,xx)
      
      if((flg_withreg.or.flg_withdamp).and..not.flg_bornonly) then
c initialize gen for remnants
         call gen(sigremnant,ndiminteg,xgridrm,ymaxrm,ymaxratrm,
     1        xmmmrm,ifoldrm,0,mcalls,icalls,xx)
      endif

      if(parallelstages.eq.3) then
         call writestat('st3')
c force compute upper bound for radiation (iret ignored)
         if(rnd_cwhichseed.ne.'none')
     1     call setrandom(rnd_initialseed,rnd_i1,rnd_i2)
         if (.not.flg_LOevents) then
            call do_maxrat(mcalls,icalls,0,iret)
         endif
         call pwhg_exit(0)
      endif

c force loading bound for radiation (iret ignored)
      iret=0
      if (.not.flg_LOevents) then
         call do_maxrat(mcalls,icalls,1,iret)
      endif
      if(iret.ne.0) then
         call pwhg_exit(0)
      endif
      if(rnd_cwhichseed.ne.'none')
     1     call setrandom(rnd_initialseed,rnd_i1,rnd_i2)
      end

      subroutine bbinitgrids
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_flg.h'
      include 'pwhg_rnd.h'
      include 'pwhg_rad.h'
      include 'multigrid.h'
      integer iret1,iret2,iun
      real * 8 sigbtl,errbtl,sigrm,errrm
      real * 8 btilde,sigremnant
      integer ncall1,ncall1rm,ncall2,ncall2rm,itmx1,itmx1rm,
     1     itmx2,itmx2rm
      real * 8 xx(ndiminteg)
      include 'cgengrids.h'
      character * 40 mergelabels
      character * 40 filename
      character * 20 pwgprefix
      integer lprefix
      common/cpwgprefix/pwgprefix,lprefix
      integer j,k,mcalls,icalls,imode,iunstat
      logical savewithnegweights,multigrid
      real * 8 powheginput
      external btilde,sigremnant,powheginput
      call newunit(iunstat)
      open(unit=iunstat,file=mergelabels(pwgprefix,rnd_cwhichseed,
     1     'stat.dat',' '),status='unknown')
      if(rnd_cwhichseed.ne.'none')
     1     call setrandom(rnd_initialseed,rnd_i1,rnd_i2)
      if (powheginput('#testplots').eq.1d0) call init_hist
      if(flg_withnegweights) then
         write(*,*)' POWHEG: Computing pos.+|neg.| '
     1        //' weight contribution to inclusive cross section' 
      else
         write(*,*)' POWHEG: Computing positive weight'
     1        //' contribution to inclusive cross section' 
      endif
      ncall2=powheginput('ncall2')
      itmx2=powheginput('itmx2')
      ncall2rm=powheginput('#ncall2rm')
      if(ncall2rm.lt.0) ncall2rm=ncall2
      itmx2rm=powheginput('#itmx2rm')
      if(itmx2rm.lt.0) itmx2rm=itmx2
      if(ncall2*itmx2.le.0) then
         write(*,*) 'ncall2=',ncall2,'itmx2=',itmx2
         write(*,*) 'Cannot continue'
         call pwhg_exit(0)
      endif
      flg_nlotest=.true.
      imode=1
c Totals will also be made available in the rad_tot???btl variables.
c The output in sigbtl is: positive weight only (flg_withnegweights=.false.)
c                          pos-|neg|            (flg_withnegweights=.true.)
c Results in rad_tot???btl do not depend upon flg_withnegweights
      call resettotals
      if(flg_storemintupb) call startstoremintupb('btildeupb')
      call setstage2init
      call mint(btilde,ndiminteg,ncall2,itmx2,ifold,imode,iun,
     1        xgrid,xint,xacc,nhits,ymax,ymaxrat,sigbtl,errbtl)
      if(flg_storemintupb) call stopstoremintupb
      call finaltotals
c finalize btilde output in histograms
      call pwhgaddout
      flg_nlotest=.false.
      write(*,*) 'btilde pos.   weights:', rad_totposbtl,' +-',
     1     rad_etotposbtl
      write(*,*) 'btilde |neg.| weights:', rad_totnegbtl,' +-',
     1     rad_etotnegbtl
      write(*,*) 'btilde total (pos.-|neg.|):', rad_totbtl,' +-',
     1     rad_etotbtl
      write(iunstat,*) 'btilde pos.   weights:', rad_totposbtl,' +-',
     1     rad_etotposbtl
      write(iunstat,*) 'btilde |neg.| weights:', rad_totnegbtl,' +-',
     1     rad_etotnegbtl
      write(iunstat,*) 'btilde Total (pos.-|neg.|):', rad_totbtl,
     1     ' +-',rad_etotbtl
c Now compute the remnant contributions
      if((flg_withreg.or.flg_withdamp).and..not.flg_bornonly) then
         write(*,*)' POWHEG: Computing remnant'//
     1        ' and/or regular remnants'
         flg_nlotest=.true.
         imode=1
         if(flg_storemintupb) call startstoremintupb('remnupb')
         call samegridasbtilde(ndiminteg,
     1        ncall2,itmx2,ifold,xgrid,
     1        ncall2rm,itmx2rm,ifoldrm,xgridrm)
         call mint(sigremnant,ndiminteg,ncall2rm,itmx2rm,ifoldrm,imode,
     1        iun,xgridrm,xintrm,xaccrm,nhitsrm,
     2        ymaxrm,ymaxratrm,sigrm,errrm)
         if(flg_storemintupb) call stopstoremintupb
c     add finalized remnant contributions in histograms
         call pwhgaddout
         flg_nlotest=.false.
      else
         sigrm=0
         errrm=0
      endif
      rad_totrm=sigrm
      rad_etotrm=errrm
c rad_totgen is used for the generation of the events.
c btilde and remnant event are chosen in proportion to
c rad_totbtlgen and rad_totrm.
      if(flg_withnegweights) then
         rad_totbtlgen=rad_totabsbtl
         rad_etotbtlgen=rad_etotabsbtl
      else
c notice: this is correct only if the negative fraction is
c negligible
         rad_totbtlgen=rad_totbtl
         rad_etotbtlgen=rad_etotbtl
      endif
      rad_totgen=rad_totrm+rad_totbtlgen
      rad_etotgen=sqrt(rad_etotbtlgen**2+rad_etotrm**2)
      rad_tot=rad_totrm+rad_totbtl
      rad_etot=sqrt(rad_etotbtl**2+rad_etotrm**2)
      
c Grids are stored in all cases; they contain informations in
c rad_tot* and rad_etot* variables. The upper bound informations
c in pwggrid files is used only if mintupb files have not been saved.
      call storegrids(xgrid,ymax,ymaxrat,
     1 xgridrm,ymaxrm,ymaxratrm,ifold,ifoldrm,ncall2,itmx2,'grid')
c Output NLO histograms
      if (powheginput('#testplots').eq.1d0) then
         if (powheginput('#qcd_order').eq.3d0) then
            filename=mergelabels(pwgprefix,rnd_cwhichseed,'NNLO',' ')
         else
            filename=mergelabels(pwgprefix,rnd_cwhichseed,'NLO',' ')
         endif
         call pwhgtopout(filename)
      endif
         
      if(flg_withreg.or.flg_withdamp) then
         write(iunstat,*) ' Remnant cross section in pb',
     1        rad_totrm,'+-',rad_etotrm
      endif
      
      if (flg_weightedev) then
         write(iunstat,*) 
     1        ' total (btilde+remnants) cross section times '
         write(iunstat,*) ' suppression factor in pb',
     1        rad_tot,'+-',rad_etot         
      else
         write(iunstat,*) 
     1        ' total (btilde+remnants) cross section in pb',
     2        rad_tot,'+-',rad_etot
      endif
      write(iunstat,*) ' negative weight fraction:',
     1     rad_totnegbtl/(2*rad_totnegbtl+rad_tot)
      
      if(flg_withreg.or.flg_withdamp) then
         write(*,*) ' Remnant cross section in pb',
     1        rad_totrm,'+-',rad_etotrm
      endif
      
      if (flg_weightedev) then
         write(*,*) ' total (btilde+remnants) cross section times '
         write(*,*) ' suppression factor in pb',
     1        rad_tot,'+-',rad_etot         
      else
         write(*,*) ' total (btilde+remnants) cross section in pb',
     1        rad_tot,'+-',rad_etot
      endif
      write(*,*) ' negative weight fraction:',
     1     rad_totnegbtl/(2*rad_totnegbtl+rad_tot)
      close(iunstat)
      end

      subroutine samegridasbtilde(ndiminteg,
     1        ncall2,itmx2,ifold,xgrid,
     1        ncall2rm,itmx2rm,ifoldrm,xgridrm)
c if the flag stage2init is set to 1 in powheg.input, the remnant
c cross section is computed with the same grid as the btilde cross
c section. In this way, the total cross section (times the suppression
c factor) computed with powheg (when using the same importance sampling
c grids) is identical whether withdamp is set or
c not, irrespective of the statistics. Useful for debugging.
      implicit none
      integer nintervals
      parameter (nintervals=50)
      integer ndiminteg,
     1        ncall2,itmx2,ifold(ndiminteg),
     1        ncall2rm,itmx2rm,ifoldrm(ndiminteg)
      real * 8 xgrid(0:nintervals,ndiminteg),
     1       xgridrm(0:nintervals,ndiminteg)
      real * 8 powheginput
      if(powheginput("#stage2init").eq.1d0) then
         ncall2rm = ncall2
         itmx2rm = itmx2
         ifoldrm = ifold
         xgridrm = xgrid
         call resetrandom
      endif
      end

      subroutine setstage2init
      implicit none
      real * 8 powheginput
      if(powheginput("#stage2init").eq.1d0) then
         call resetrandom
      endif
      end

      subroutine bbinitxgrids(iparallel)
c iparallel = 0: generate xgrid file
c iparallel = 1: generate gridinfo files for parallel grid      
      implicit none
      character * 6 tag,prevtag
      integer iparallel,itmp
      integer iteration,imode,j,itmx1,ncall1,ncall1rm,itmx1rm,iun,iret
      integer btilde,sigremnant
      include 'nlegborn.h'
      include 'pwhg_flg.h'
      include 'pwhg_rnd.h'
      include 'pwhg_flst.h'
      include 'pwhg_rad.h'
      include 'pwhg_par.h'
      include 'cgengrids.h'
      include 'phspcuts.h'
      character * 20 pwgprefix
      character * 40 mergelabels
      integer lprefix
      common/cpwgprefix/pwgprefix,lprefix
      real * 8 sigbtl,errbtl,sigrm,errrm
      real * 8 random,powheginput
      logical savewithnegweights
      external btilde,sigremnant,random,powheginput
c No folding while generating importance sampling grids
      do j=1,ndiminteg
         ifold(j)=1
         ifoldrm(j)=1
      enddo
c         
      if(iparallel.eq.1) then
         iteration = powheginput('#xgriditeration')
         if(iteration.gt.par_maxxgriditerations) then
            write(*,*) ' POWHEG is compiled with a maximum'
            write(*,*) ' number of x-grid iterations=',
     1           par_maxxgriditerations
            write(*,*) ' increase the par_maxxgriditerations parameter'
            write(*,*) ' in the pwhg_par.h file to use more.'
            write(*,*) ' (BUT YOU SHOULD NOT NEED MORE THAN 3 or 4!'
            write(*,*) ' increase ncall1 instead)'
            call pwhg_exit(-1)
         endif
         if(iteration.le.0) iteration = 1
         write(tag,'(i3)') iteration
         tag=adjustl(tag)
         tag='xg'//tag
         if(iteration.gt.1) then
            write(prevtag,'(i3)') iteration-1
            prevtag=adjustl(prevtag)
            prevtag='xg'//prevtag
         endif
      else
         iteration = 0
         tag=' '
         prevtag=' '
      endif
c In this block we compute the importance sampling grid
      write(*,*)
      write(*,*)' POWHEG: initialization'
      write(*,*)' Computing the integral of the absolute value'
      if(flg_weightedev) then
         write(*,*)' of the cross section times the suppression' 
         write(*,*)' factor to set up the adaptive grid'
      else
         write(*,*)' of the cross section to set up the adaptive grid'
      endif
c If parallel operations, set a different random number also for
c different grid iterations
      if(iparallel.eq.1) then
         itmp=0
         call resetrandom
         do j=1,iteration
            itmp = random()*1d9
         enddo
         call setrandom(rnd_initialseed+itmp,rnd_i1,rnd_i2)
      endif
      if(iteration.le.1) then
         call initxgrid(xgrid,ndiminteg)
      else
         call loadgridinfo(mergelabels('btl',prevtag,' ',' '),.false.,
     1        xgrid,xint,iret)
         if(iret.ne.0) then
            write(*,*) 'cannot find level '//prevtag//'gridinfo files'
            call pwhg_exit()
         endif
      endif
c The actual value of the grid is the one to be saved in the
c gridinfo file, since loadgridinfo updates the grid by itself
      xgrid0 = xgrid
      ncall1 = powheginput("ncall1")
      ncall1rm = powheginput("#ncall1rm")
c ncall1rm is additional, if it is not present use the standard one:
      if (ncall1rm.lt.0d0) ncall1rm = ncall1
      itmx1 = powheginput("itmx1")
      itmx1rm = powheginput("#itmx1rm")
      if(itmx1rm.lt.0) itmx1rm=itmx1
c with parallel grids only one iteration is allowed
      if(iparallel.eq.1) then
         itmx1 = 1
         itmx1rm = 1
      endif
      call newunit(iun)
      call regridplotopen(mergelabels(pwgprefix(1:lprefix),tag,
     1     rnd_cwhichseed,'btlgrid.top'))
      write(*,*)' result +- errtot (picobarn) for each iteration'
      flg_nlotest=.false.
      imode=0
c     AK for phspcuts
      if(phspcuts) call init_hist
      savewithnegweights=flg_withnegweights
      flg_withnegweights=.true.
c ********** CALL to mint for btilde

      call mint(btilde,ndiminteg,ncall1,itmx1,ifold,imode,iun,
     1     xgrid,xint,xacc,nhits,ymax,ymaxrat,sigbtl,errbtl)

c **********
      call regridplotclose
      if(iteration.ge.1) then
         call storegridinfo('btl-'//trim(tag),xgrid0,xint,
     1        xacc,nhits,ndiminteg)
      endif
      flg_withnegweights=savewithnegweights
      close(iun)
      if((flg_withreg.or.flg_withdamp).and..not.flg_bornonly) then
         write(*,*) ' Computing the integral of the'//
     1        ' remnant cross section'
         if(flg_weightedev) then
            write(*,*) 'times the suppression factor'
         endif
         write(*,*) ' to set up the adaptive grid'
         flg_nlotest=.false.
         imode=0
         if(iteration.le.1) then
            call initxgrid(xgridrm,ndiminteg)
         else
            call loadgridinfo(mergelabels('rmn',prevtag,' ',' '),
     1           .false.,xgridrm,xintrm,iret)
         endif
         call newunit(iun)
         call regridplotopen(mergelabels(pwgprefix(1:lprefix),tag,
     1     rnd_cwhichseed,'rmngrid.top'))
         xgrid0rm=xgridrm
c ********** CALL to mint for remnants

         call mint(sigremnant,ndiminteg,ncall1rm,itmx1rm,ifoldrm,imode,
     1        iun,xgridrm,xintrm,xaccrm,nhitsrm,
     1        ymaxrm,ymaxratrm,sigrm,errrm)

c **********
         call regridplotclose
         if(iteration.ge.1) then
            call storegridinfo('rmn-'//trim(tag),xgrid0rm,xintrm,
     1           xaccrm,nhitsrm,ndiminteg)
         endif
      endif
      if(iparallel.eq.0) then
         call storexgrid(xgrid,xint,xgridrm,xintrm)
c      else
c         call newunit(iun)
c         open(unit=iun,
c     1        file=mergelabels('stage-xiterations',tag,'.dat',' '),
c     1        status='unknown')
c         write(iun,*)'done'
c         close(iun)
      endif
      write(*,*)' Importance sampling x grids generated and stored'
      end

      subroutine pwhg_openoutput(iun,string,suffix)
      implicit none
      integer iun
      character * (*) string
      character * (*) suffix
      include 'pwhg_rnd.h'
      character * 20 pwgprefix
      integer lprefix
      common/cpwgprefix/pwgprefix,lprefix
      character * 40 mergelabels
      open(unit=iun,file=mergelabels(pwgprefix(1:lprefix)
     1     //trim(string),rnd_cwhichseed,suffix,' '),status='unknown')
      end

      subroutine gen_btilde(mcalls,icalls)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_rad.h'
      include 'pwhg_flg.h'
      integer mcalls,icalls
      include 'cgengrids.h'
      real * 8 xx(ndiminteg)      
      real * 8 btilde
      external btilde
c use these to provide an estimate of the cross section while generating an event
      real * 8 sigma, sigma2
      integer isigma
      common/gencommon/sigma,sigma2,isigma
      if(mcalls == 0) then
         gen_sigma  = 0
         gen_sigma2 = 0
         gen_isigma = 0
         gen_totev  = 0
      endif
      call gen(btilde,ndiminteg,xgrid,ymax,ymaxrat,xmmm,ifold,1,
     1     mcalls,icalls,xx)
      gen_sigma  = gen_sigma  + sigma
      gen_sigma2 = gen_sigma2 + sigma2
      gen_isigma = gen_isigma + isigma
      gen_totev  = gen_totev  + rad_genubexceeded
      gen_mcalls    = mcalls
      call setcnt("btilde cross section used:", rad_totgen-rad_totrm)
      call setcnt("btilde cross section estimate:",gen_sigma/gen_isigma)
      call setcnt("btilde cross section estimate num. points:",
     $     dble(gen_isigma))
      call setcnt("btilde cross section error estimate:",
     1     sqrt(((gen_sigma2/gen_isigma)-(gen_sigma/gen_isigma)**2)/
     $     gen_isigma))
      if(flg_ubexcess_correct) then
         call setcnt("btilde bound violation correction factor:",
     1        gen_totev/gen_mcalls)
      endif
      end
      
      subroutine gen_sigremnant
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_flg.h'
      include 'pwhg_rad.h'
      include 'cgengrids.h'
      real * 8 xx(ndiminteg)
      integer mcalls,icalls
      logical savelogical
      real * 8 sigremnant
      external sigremnant
c use these to provide an estimate of the cross section while generating an event
      real * 8 sigma, sigma2
      integer isigma
      common/gencommon/sigma,sigma2,isigma
      if(mcalls == 0) then
         gen_sigmarm  = 0
         gen_sigma2rm = 0
         gen_isigmarm = 0
         gen_totevrm  = 0
      endif
c communicate file to load upper bound data
      savelogical=flg_fastbtlbound
      flg_fastbtlbound=.false.
      call gen(sigremnant,ndiminteg,xgridrm,ymaxrm,ymaxratrm,
     1    xmmmrm,ifoldrm,1,mcalls,icalls,xx)
      flg_fastbtlbound=savelogical
      gen_sigmarm  = gen_sigmarm  + sigma
      gen_sigma2rm = gen_sigma2rm + sigma2
      gen_isigmarm = gen_isigmarm + isigma
      gen_totevrm  = gen_totevrm  + rad_genubexceeded
      gen_mcallsrm    = mcalls
      call setcnt("remnant cross section used:", rad_totrm)
      call setcnt("remnant cross section estimate:",
     $     gen_sigmarm/gen_isigmarm)
      call setcnt("remnant cross section error estimate:",
     1     sqrt(((gen_sigma2rm/gen_isigmarm)
     2     -(gen_sigmarm/gen_isigmarm)**2)/gen_isigmarm))
      call setcnt("remnant cross section estimate num. points:",
     1     dble(gen_isigmarm))
      if(flg_ubexcess_correct) then
         call setcnt("remnant bound violation correction factor:",
     1        gen_totevrm/gen_mcallsrm)
      endif
      end


      subroutine storegrids(xgrid,ymax,ymaxrat,xgridrm,ymaxrm,
     1                ymaxratrm,ifold,ifoldrm,ncall2,itmx2,gridtag)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_pdf.h'
      include 'pwhg_rad.h'
      include 'pwhg_rnd.h'
      real * 8 xgrid(0:50,ndiminteg),ymax(50,ndiminteg),
     1     ymaxrat(50,ndiminteg),xgridrm(0:50,ndiminteg),
     2     ymaxrm(50,ndiminteg),ymaxratrm(50,ndiminteg)
      integer nbins
      parameter (nbins=50)
      integer ifold(ndiminteg),ifoldrm(ndiminteg),ncall2,itmx2
      character *(*) gridtag
      character * 20 pwgprefix
      integer lprefix
      common/cpwgprefix/pwgprefix,lprefix
      integer j,k,iun
      call newunit(iun)
      if(rnd_cwhichseed.eq.'none') then
         open(unit=iun,file=pwgprefix(1:lprefix)//gridtag//'.dat',
     #     form='unformatted',status='unknown')
      else
         open(unit=iun,
     1        file=pwgprefix(1:lprefix)//gridtag//'-'//rnd_cwhichseed//
     2        '.dat',form='unformatted',status='unknown')
      endif
      write(iun) ((xgrid(j,k),k=1,ndiminteg),j=0,nbins)
      write(iun) ((ymax(j,k),k=1,ndiminteg),j=1,nbins)
      write(iun) ((ymaxrat(j,k),k=1,ndiminteg),j=1,nbins)
      write(iun) ((xgridrm(j,k),k=1,ndiminteg),j=0,nbins)
      write(iun) ((ymaxrm(j,k),k=1,ndiminteg),j=1,nbins)
      write(iun) ((ymaxratrm(j,k),k=1,ndiminteg),j=1,nbins)
      write(iun) (ifold(k),k=1,ndiminteg)
      write(iun) (ifoldrm(k),k=1,ndiminteg)
      write(iun) ncall2*itmx2
      write(iun) kn_sbeams, pdf_ih1, pdf_ih2, pdf_ndns1, pdf_ndns2
      write(iun)
     1     rad_totbtl,rad_etotbtl,
     2     rad_totabsbtl,rad_etotabsbtl,
     3     rad_totposbtl,rad_etotposbtl,
     4     rad_totnegbtl,rad_etotnegbtl,
     5     rad_totrm,rad_etotrm,
     6     rad_totbtlgen,rad_etotbtlgen,
     7     rad_totgen,rad_etotgen,
     8     rad_tot,rad_etot
      close(iun)
      end

      subroutine loadgrids(iret,xgrid,ymax,ymaxrat,xgridrm,ymaxrm,
     #           ymaxratrm,ifold,ifoldrm,gridtag)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_pdf.h'
      include 'pwhg_rad.h'
      include 'pwhg_rnd.h'
      include 'pwhg_par.h'
      real * 8 xgrid(0:50,ndiminteg),ymax(50,ndiminteg),
     1     ymaxrat(50,ndiminteg),xgridrm(0:50,ndiminteg),
     2     ymaxrm(50,ndiminteg),ymaxratrm(50,ndiminteg)
      real * 8 xxgrid(0:50,ndiminteg),xymax(50,ndiminteg),
     1     xymaxrat(50,ndiminteg),xxgridrm(0:50,ndiminteg),
     2     xymaxrm(50,ndiminteg),xymaxratrm(50,ndiminteg)
      real * 8 tot(2,8),rtot(2,8)
      integer ifold(ndiminteg),ifoldrm(ndiminteg)
      integer iifold(ndiminteg),iifoldrm(ndiminteg)
      integer iret,iretcheck,jfound
      character *(*) gridtag
c
      integer ios
      integer nbins
      parameter (nbins=50)
      character * 20 pwgprefix
      integer lprefix
      common/cpwgprefix/pwgprefix,lprefix
      character * 4 chseed, firstfound
      real * 8 shx
      integer ih1x, ih2x, ndns1x, ndns2x
      integer j,k,iun,jfile,nfiles,ncall2,itmx2
      logical lpresent,manyfiles,filefound
      double precision, allocatable :: totarr(:,:)
      logical, allocatable :: goodentries(:)
      logical firsttime
      logical, save :: check_bad_st2
      real * 8 rjfound, rncall2
      real * 8 powheginput
      external powheginput
      if(powheginput('use-old-grid').eq.0) then
         iret=1
         return
      endif
      iret=0
      call newunit(iun)
      open(unit=iun,file=pwgprefix(1:lprefix)//gridtag//'.dat',
     #     form='unformatted',status='old',iostat=ios)
      if(ios.eq.0) then
         nfiles=1
         check_bad_st2 = .false.
      else
         nfiles=par_maxseeds
         manyfiles=.true.
         check_bad_st2 =
     1        powheginput("#check_bad_st2") == 1
         if(check_bad_st2) then
            allocate(totarr(2,nfiles),goodentries(nfiles))
            goodentries = .true.
            totarr = 0
         endif
      endif
c Try to open and merge a set of grid files, generated with different
c random seeds
      firsttime = .true.
 11   continue
      filefound=.false.
      jfound=0
      do jfile=1,nfiles
         if(manyfiles) then
            if(check_bad_st2) then
               if(.not.goodentries(jfile)) cycle
            endif
            write(chseed,'(i4)') jfile
            do k=1,4
               if(chseed(k:k).eq.' ') chseed(k:k)='0'
            enddo
            inquire(file=pwgprefix(1:lprefix)//gridtag//'-'//
     1           chseed//'.dat',exist=lpresent)
            if(.not.lpresent) goto 111
            open(unit=iun,file=pwgprefix(1:lprefix)//gridtag//'-'//
     1           chseed//'.dat',
     2           form='unformatted',status='old',iostat=ios)
            if(ios.ne.0) then
               iret=-1
               return
            else
               write(*,*)
     1              ' Opened ',pwgprefix(1:lprefix)//gridtag//'-'//
     2              chseed//'.dat'
            endif
         endif
         filefound=.true.
         read(iun,iostat=ios) ((xxgrid(j,k),k=1,ndiminteg),j=0,nbins)
         if(ios.ne.0) goto 998
         read(iun,iostat=ios) ((xymax(j,k),k=1,ndiminteg),j=1,nbins)
         if(ios.ne.0) goto 998
         read(iun,iostat=ios) ((xymaxrat(j,k),k=1,ndiminteg),j=1,nbins)
         if(ios.ne.0) goto 998
         read(iun,iostat=ios) ((xxgridrm(j,k),k=1,ndiminteg),j=0,nbins)
         if(ios.ne.0) goto 998
         read(iun,iostat=ios) ((xymaxrm(j,k),k=1,ndiminteg),j=1,nbins)
         if(ios.ne.0) goto 998
         read(iun,iostat=ios)((xymaxratrm(j,k),k=1,ndiminteg),j=1,nbins)
         if(ios.ne.0) goto 998
         read(iun,iostat=ios) (iifold(k),k=1,ndiminteg)
         if(ios.ne.0) goto 998
         read(iun,iostat=ios) (iifoldrm(k),k=1,ndiminteg)
         if(ios.ne.0) goto 998
         read(iun,iostat=ios) ncall2
         if(powheginput("#ncallfrominput").eq.1) then
            ncall2=powheginput("ncall2")
            itmx2=powheginput("itmx2")
            ncall2=ncall2*itmx2
         endif
         if(ios.ne.0) goto 998
         read(iun,iostat=ios) shx, ih1x, ih2x, ndns1x, ndns2x
         if(ios.ne.0) goto 998
         if(shx.ne.kn_sbeams.or.ih1x.ne.pdf_ih1.or.ih2x.ne.pdf_ih2
     1      .or.ios.ne.0)
     2        goto 998
         read(iun,iostat=ios) ((tot(k,j),k=1,2),j=1,8)
         if(check_bad_st2) then
            totarr(:,jfile)=tot(:,1)
         endif
         if(ios.ne.0) goto 998
         jfound=jfound+1
         if(jfound.lt.2) then
            firstfound = chseed
            do k=1,ndiminteg
               do j=0,nbins
                  xgrid(j,k)=xxgrid(j,k)
                  xgridrm(j,k)=xxgridrm(j,k)
               enddo
               ifold(k)=iifold(k)
               ifoldrm(k)=iifoldrm(k)
            enddo
            do k=1,ndiminteg
               do j=1,nbins
                  ymax(j,k)=xymax(j,k)
                  ymaxrm(j,k)=xymaxrm(j,k)
                  ymaxrat(j,k)=xymaxrat(j,k)
                  ymaxratrm(j,k)=xymaxratrm(j,k)
               enddo
            enddo
            do k=1,2
               do j=1,8
                  rtot(k,j)=tot(k,j)
               enddo
            enddo
         else
            do k=1,ndiminteg
               do j=0,nbins
                  if(xgrid(j,k).ne.xxgrid(j,k).or.
     1                 xgridrm(j,k).ne.xxgridrm(j,k)) then
                     write(*,*) ' error loading grids: '
                     write(*,*)  pwgprefix(1:lprefix)//gridtag//'-'//
     1          rnd_cwhichseed//'.dat does not have the same importance'
                    write(*,*) 'sampling grid as ',pwgprefix(1:lprefix)
     1              //gridtag//'-'//firstfound//'.dat'
                     call pwhg_exit(-1)
                  endif
               enddo
               if(ifold(k).ne.iifold(k)
     1              .or.ifoldrm(k).ne.iifoldrm(k)) then
                  write(*,*) ' error loading grids: '
                  write(*,*)  pwgprefix(1:lprefix)//gridtag//'-'//
     1                 rnd_cwhichseed//
     2                 '.dat does not have the same folding as '
                  write(*,*) pwgprefix(1:lprefix)//gridtag//'-'
     1                 //firstfound//'.dat'
                  call pwhg_exit(-1)
               endif
            enddo
            do k=1,ndiminteg
               do j=1,nbins
                  ymax(j,k)=max(ymax(j,k),xymax(j,k))
                  ymaxrm(j,k)=max(ymaxrm(j,k),xymaxrm(j,k))
                  ymaxrat(j,k)=max(ymaxrat(j,k),xymaxrat(j,k))
                  ymaxratrm(j,k)=max(ymaxratrm(j,k),xymaxratrm(j,k))
               enddo
            enddo
            do j=1,8
c turn these to reals; very large integer can overflow in fortran
               rjfound = jfound
               rncall2 = ncall2
               rtot(2,j)=sqrt((rtot(2,j)**2*(rjfound-1)**2+tot(2,j)**2)
     1              /rjfound**2+(rjfound-1)*(rtot(1,j)-tot(1,j))**2/
     2              (rjfound**3*rncall2))
               rtot(1,j)=(rtot(1,j)*(rjfound-1)+tot(1,j))/rjfound
            enddo
         endif
         rad_totbtl     =rtot(1,1)
         rad_etotbtl    =rtot(2,1)
         rad_totabsbtl  =rtot(1,2)
         rad_etotabsbtl =rtot(2,2)
         rad_totposbtl  =rtot(1,3)
         rad_etotposbtl =rtot(2,3)
         rad_totnegbtl  =rtot(1,4)
         rad_etotnegbtl =rtot(2,4)
         rad_totrm      =rtot(1,5)
         rad_etotrm     =rtot(2,5)
         rad_totbtlgen  =rtot(1,6)
         rad_etotbtlgen =rtot(2,6)
         rad_totgen     =rtot(1,7)
         rad_etotgen    =rtot(2,7)
         rad_tot        =rtot(1,8)
         rad_etot       =rtot(2,8)
         close(iun)
 111     continue
      enddo
      if(filefound) then
         if(manyfiles .and. check_bad_st2) then
c check that the different runs are more or less consistent
            if(firsttime) then
               call check_stat_consistency(nfiles,totarr,
     1              goodentries,iretcheck)
               if(iretcheck.eq.-1) then
                  firsttime = .false.
                  goto 11
               endif
            endif
            deallocate(goodentries,totarr)
         endif
         return
      endif
 998  continue
      iret=-1
      end

      subroutine  check_stat_consistency(nentries,res,goodentries,iret)
      implicit none
      integer nentries,iret
      integer indices(nentries)
      double precision res(2,nentries)
      logical goodentries(nentries)
      double precision average,weight,tmpav,tmpweight
      double precision tmp2(2),ow,oav
      integer nonz,j,k,itmp
      logical :: ex=.true.
c DEBUG
c      res(1,7) = 3.7
c      res(2,7) = 1
c      res(1,5) = 3.7
c      res(2,5) = 1
c      res(1,3) = 3.7
c      res(2,3) = 1
c end DEBUG
      do j=1,nentries
         indices(j)=j
      enddo
c bubblesort
      do while(ex)
         ex = .false.
         do j=1,nentries-1
c swap in growing order, but put zeros at the end
            if((res(2,j+1).ne.0.and.res(2,j+1).lt.res(2,j)).or.
     1        (res(2,j+1).ne.0.and.res(2,j).eq.0)) then
               tmp2 = res(:,j)
               res(:,j) = res(:,j+1)
               res(:,j+1) = tmp2
               itmp = indices(j)
               indices(j) = indices(j+1)
               indices(j+1) = itmp
               ex = .true.
            endif
         enddo
      enddo
      do j=1,nentries
         if(res(2,j).eq.0) then
            nonz = j - 1
            exit
         endif
      enddo
c Compute the average. Neglect zero entries.
      average = 0
      weight = 0
      do j=1,nonz
         if(res(1,j).ne.0) then
            oav = average
            ow = weight
            average = (average*(j-1) + res(1,j))/j
            weight = sqrt((weight*(j-1))**2 + res(2,j)**2)/j
            if(j.gt.1) then
c               write(*,*) ' old,new av.',oav,average
c               write(*,*) ' old,new err.',ow,weight
c               write(*,*) ' deviation:',abs(oav-average)/ow
c after half of the runs
               if(abs(oav-average)/ow.gt.10) then
                  write(*,*) ' check_stat_consistency:'
                  write(*,*)
     1                 ' The program has detected inconsistent results'
                  write(*,*) ' among different runs. The runs:'
                  write(*,*) (indices(k),k=j,nonz)
                  write(*,*) ' look suspicious. '
                  write(*,*) ' Inspect your runs at stage 2 '
                  if(nonz-j+1.gt. nonz/10) then
                     write(*,*)
     1               ' The fraction of inconsistent runs is too large'
                     write(*,*) ' exiting ...'
                     call exit(-1)
                  else
                     write(*,*)
     1               ' The fraction of inconsistent file is < 10%'
                     write(*,*) ' We discard the following files:'
                     write(*,*) (indices(k),k=j,nonz)
                     write(*,*) ' and reload the others'                     
                     do k=j,nonz
                        goodentries(indices(k)) = .false.
                     enddo
                     iret = -1
                     return
                  endif
               endif
            endif
         endif
      enddo
      iret = 0
      end


      subroutine storexgrid(xgrid,xint,xgridrm,xintrm)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_pdf.h'
      include 'pwhg_rad.h'
      real * 8 xgrid(0:50,ndiminteg),xgridrm(0:50,ndiminteg),
     1     xint,xintrm
      integer nbins
      parameter (nbins=50)
      character * 20 pwgprefix
      integer lprefix
      common/cpwgprefix/pwgprefix,lprefix
      integer j,k,iun
      call newunit(iun)
      open(unit=iun,file=pwgprefix(1:lprefix)//'xgrid.dat',
     #     form='unformatted',status='unknown')
      write(iun) ((xgrid(j,k),k=1,ndiminteg),j=0,nbins),xint
      write(iun) ((xgridrm(j,k),k=1,ndiminteg),j=0,nbins),xintrm
      write(iun) kn_sbeams, pdf_ih1, pdf_ih2, pdf_ndns1, pdf_ndns2
      close(iun)
      end

      subroutine loadxgrids(iret)
c loads an integration grid from file <prefix>xgrid.dat, if found
c (returns iret=0). If not found iret=-1. It checks that some key
c parameters are the same as for the current run
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_pdf.h'
      include 'pwhg_rad.h'
      include 'cgengrids.h'
      integer iret
c
      integer ios
      integer nbins
      parameter (nbins=50)
      character * 20 pwgprefix
      integer lprefix
      common/cpwgprefix/pwgprefix,lprefix
      real * 8 shx
      integer ih1x, ih2x, ndns1x, ndns2x
      integer j,k,iun
      real * 8 powheginput
      external powheginput
      if(powheginput('use-old-grid').eq.0) then
         iret=1
         return
      endif
      call newunit(iun)
      open(unit=iun,file=pwgprefix(1:lprefix)//'xgrid.dat',
     #     form='unformatted',status='old',iostat=ios)
      if(ios.ne.0) then
         iret=-1
         return
      endif
      read(iun,iostat=ios) ((xgrid(j,k),k=1,ndiminteg),j=0,nbins),xint
      read(iun,iostat=ios) ((xgridrm(j,k),k=1,ndiminteg),j=0,nbins),
     1     xintrm
      read(iun,iostat=ios) shx, ih1x, ih2x, ndns1x, ndns2x
      if(shx.ne.kn_sbeams.or.ih1x.ne.pdf_ih1.or.ih2x.ne.pdf_ih2
     #  .or.ios.ne.0) then
         iret=-1
         close(iun)
         return
      endif
      close(iun)
      iret=0
      end


      subroutine loadlatestxgridinfo(iret)
c loads x information from gridinfo files, and computes the xgrid.
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flg.h'
      include 'pwhg_par.h'
      integer iret
      character * 40 mergelabels
      character * 6 chnum
      integer j,iteration
      logical probe
      real * 8 ans
      include 'cgengrids.h'
      probe=.true.
      do j=par_maxxgriditerations,1,-1
         write(chnum,'(i4)') j
         chnum='xg'//adjustl(chnum)
c remember, tags look like '-1', '-2', etc.
         call loadgridinfo
     1        (mergelabels('btl',chnum,' ',' '),probe,xgrid,ans,iret)
         if(iret.eq.0) then
            iteration = j
            exit
         endif
      enddo
      if(j.eq.0) then
         iret = -1
         return
      endif
      probe = .false.
      write(*,*) ' loading '//chnum//' iteration file'
      call loadgridinfo
     1     (mergelabels('btl',chnum,' ',' '),probe,xgrid,xint,iret)
      if(iret.ne.0) return
      if((flg_withreg.or.flg_withdamp).and..not.flg_bornonly) then
         call loadgridinfo
     1     (mergelabels('rmn',chnum,' ',' '),probe,xgridrm,xintrm,iret)
      endif
      end
      
      subroutine loadgridinfo(storelabel,probe,xgrid,ans,iret)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_rnd.h'
      include 'pwhg_par.h'
      character *(*) storelabel
      logical probe
      integer nintervals,ndimmax
      parameter (nintervals=50,ndimmax=ndiminteg)
      real * 8 xgrid(0:nintervals,*),ans
      integer iret
      real * 8 xacc(0:nintervals,ndimmax)
      integer nhits(1:nintervals,ndimmax)
      real * 8 xacc0(0:nintervals,ndimmax),tmp
      integer nhits0(1:nintervals,ndimmax),ndim0
      character * 4 chseed
      character * 20 pwgprefix
      character * 100 file
      character * 40 mergelabels
      integer lprefix
      common/cpwgprefix/pwgprefix,lprefix
      integer j,k,iun,ios,kdim
      logical filefound,lpresent
      integer nfiles,jfile,jfound

c first probe if there are files to load
      nfiles=par_maxseeds
      do jfile=1,nfiles
         write(chseed,'(i4)') jfile
         do k=1,4
            if(chseed(k:k).eq.' ') chseed(k:k)='0'
         enddo
         file = mergelabels(pwgprefix(1:lprefix)//'gridinfo',
     1        storelabel,chseed,'.dat')
         inquire(file= file,exist=lpresent)
         if(lpresent) then
            if(probe) then
               iret = 0
               return
            endif
            exit
         endif
      enddo
      if(probe) then
         iret = -1
         return
      endif

      nfiles=par_maxseeds
      filefound=.false.
      jfound=0
      xacc(:,1:ndiminteg)=0
      nhits(:,1:ndiminteg)=0
      ans=0
      do jfile=1,nfiles
         write(chseed,'(i4)') jfile
         do k=1,4
            if(chseed(k:k).eq.' ') chseed(k:k)='0'
         enddo
         file = mergelabels(pwgprefix(1:lprefix)//'gridinfo',
     1        storelabel,chseed,'.dat')
         inquire(file= file,exist=lpresent)
         if(.not.lpresent) cycle
         call newunit(iun)
         open(unit=iun,file= file,
     2        form='unformatted',status='old',iostat=ios)
         if(ios.ne.0) then
            iret=-1
            return
         else
            write(*,*) ' Opened ', file
         endif
         filefound=.true.
         jfound=jfound+1
         read(iun,iostat=ios) ndim0
         if(ios.ne.0.or.ndim0.ne.ndiminteg) goto 111
         read(iun,iostat=ios) tmp
         ans = ans + tmp
         if(ios.ne.0) goto 111
         read(iun,iostat=ios) xgrid(0:nintervals,1:ndiminteg)
         if(ios.ne.0) goto 111
         read(iun,iostat=ios) xacc0(0:nintervals,1:ndiminteg)
         if(ios.ne.0) goto 111
         xacc(:,1:ndiminteg) =
     1        xacc(:,1:ndiminteg) + xacc0(:,1:ndiminteg)
         read(iun,iostat=ios) nhits0(1:nintervals,1:ndiminteg)
         if(ios.ne.0) goto 111
         nhits(:,1:ndiminteg) =
     1        nhits(:,1:ndiminteg) + nhits0(:,1:ndiminteg)
         close(iun)
      enddo
c      write(*,*) ' loadgridinfo: found ',jfound,
c     1     ' files with grid information'
      if(filefound) then
         if(rnd_cwhichseed.eq.'none') then
            call regridplotopen(pwgprefix(1:lprefix)//'-'//
     1          storelabel(1:3) //'grid.top')
         else
            call regridplotopen(pwgprefix(1:lprefix)//'-'//
     1           rnd_cwhichseed//'-'//storelabel(1:3)//'grid.top')
         endif
         ans = ans/jfound
         do kdim=1,ndiminteg
            call regrid(xacc(0,kdim),xgrid(0,kdim),
     1           nhits(1,kdim),kdim,nintervals)
         enddo
         iret = 0
         call regridplotclose
         return
      endif

      iret = -1
      return
 111  continue
      write(*,*) ' loadgridinfo: problems loading files'
      call pwhg_exit(-1)
      end


      subroutine storegridinfo(storelabel,xgrid,xint,xacc,nhits,ndim)
c stores accumulated results and number of hits in gridinfo files
      implicit none
      include 'nlegborn.h'
      include 'pwhg_rnd.h'
      character *(*) storelabel
      integer nintervals,ndimmax
      parameter (nintervals=50,ndimmax=ndiminteg)
      real * 8 xgrid(0:nintervals,*),xint
      real * 8 xacc(0:nintervals,*)
      integer nhits(1:nintervals,*),ndim
      character * 20 pwgprefix
      character * 40 mergelabels
      integer lprefix
      common/cpwgprefix/pwgprefix,lprefix
      integer j,k,iun
      call newunit(iun)
      open(unit=iun,
     1     file=mergelabels(pwgprefix(1:lprefix)//'gridinfo',
     2     storelabel,rnd_cwhichseed,'.dat'),
     3     form='unformatted',status='unknown')
      write(iun) ndim
      write(iun) xint
      write(iun) xgrid(0:nintervals,1:ndim)
      write(iun) xacc(0:nintervals,1:ndim)
      write(iun) nhits(1:nintervals,1:ndim)
      close(iun)
      end


      function mergelabels(lab1,lab2,lab3,lab4)
c puts together up to 4 labels, separating them with '-'.
c Empty labels are ignored.
      character * 40 mergelabels
      character *(*) lab1,lab2,lab3,lab4
      integer where
      mergelabels=' '
      if(lab1.ne.' '.and.lab1.ne.'none') then
         mergelabels=adjustl(lab1)
         where=index(mergelabels,' ')
         if(where.eq.0) goto 999
      endif
      if(lab2.ne.' '.and.lab2.ne.'none') then
         mergelabels(where:)='-'//adjustl(lab2)
         where=index(mergelabels,' ')
         if(where.eq.0) goto 999
      endif
      if(lab3.ne.' '.and.lab3.ne.'none') then
         mergelabels(where:)='-'//adjustl(lab3)
         where=index(mergelabels,' ')
         if(where.eq.0) goto 999
      endif
      if(lab4.ne.' '.and.lab4.ne.'none') then
         mergelabels(where:)='-'//adjustl(lab4)
         where=index(mergelabels,' ')
         if(where.eq.0) goto 999
      endif
c get rid of hiphen before extension
      where=index(mergelabels,'-.',.true.)
      if(where.ne.0) then
         mergelabels(where:)=mergelabels(where+1:)
      endif
      return
 999  write(*,*) ' mergelabels: strings too long'
      call pwhg_exit(-1)
      end

      subroutine writestat(stage)
      implicit none
      character *(*) stage
      character * 40 mergelabels
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_rad.h'
      include 'pwhg_rnd.h'
      include 'pwhg_flg.h'
      integer iunstat
      character * 20 pwgprefix
      integer lprefix
      common/cpwgprefix/pwgprefix,lprefix
      call newunit(iunstat)
      open(unit=iunstat,file=mergelabels(pwgprefix,stage,rnd_cwhichseed,
     1     'stat.dat'),status='unknown')

      write(iunstat,*) 'btilde pos.   weights:', rad_totposbtl,' +-',
     1     rad_etotposbtl
      write(iunstat,*) 'btilde |neg.| weights:', rad_totnegbtl,' +-',
     1     rad_etotnegbtl
      write(iunstat,*) 'btilde Total (pos.-|neg.|):', rad_totbtl,
     1     ' +-',rad_etotbtl
      if(flg_withreg.or.flg_withdamp) then
         write(iunstat,*) ' Remnant cross section in pb',
     1        rad_totrm,'+-',rad_etotrm
      endif
      if (flg_weightedev) then
         write(iunstat,*) 
     1        ' total (btilde+remnants) cross section times '
         write(iunstat,*) ' suppression factor in pb',
     1        rad_tot,'+-',rad_etot         
      else
         write(iunstat,*) 
     1        ' total (btilde+remnants) cross section in pb',
     2        rad_tot,'+-',rad_etot
      endif
      write(iunstat,*) ' negative weight fraction:',
     1     rad_totnegbtl/(2*rad_totnegbtl+rad_tot)
      close(iunstat)
      end
