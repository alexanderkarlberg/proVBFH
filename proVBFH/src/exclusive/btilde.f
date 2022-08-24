      function btilde(xx,www0,ifirst,imode,
     1     retval,retval0)
c retval is the function return value retval0 is an 'avatar' function
c that has similar value, but is much easier to compute (i.e. the Born
c term in this case) imode = 0 compute retval0 only.  imode = 1 compute
c retval, retval0 return value: 0: success; 1: retval0 was not computed
c (this function does not support an avatar function)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_rad.h'
      include 'pwhg_flg.h'
      include 'pwhg_math.h'
      include 'phspcuts.h'
      include 'incl2pwhg.h'
c     independent variables for real graph: number of final state legs
c     times 3, take away 4 for 4-momentum conservation, add 2 for x_1
c     and x_2, and take away an overall azimuth
      real * 8 xx(ndiminteg),www0,retval,retval0
      real * 8 xrad(3)
      real * 8 xborn(ndiminteg-3)
      integer btilde,ifirst,imode,iret
      real * 8 resborn(maxprocborn),resvirt(maxprocborn),
     #     resreal(maxprocborn),rescoll(maxprocborn)
      real * 8 results(maxprocborn)
      real * 8 tmp,suppfact,www,wwwtot
      integer j
      save results,resborn,resvirt,wwwtot,suppfact
      real * 8 seconds
      real *8 totborn,totvirt,ptotborn
      logical pwhg_isfinite 
      external pwhg_isfinite
      real * 8 powheginput
      logical ini,btildebornon,btildevirton,btildecollon,btilderealon
      data ini/.true./
      save ini,btildebornon,btildevirton,btildecollon,btilderealon
      
      if(ini) then
         if(powheginput('#qcd_order').eq.1d0) then
            ! In this case we are actually only running the inclusive part of the code
            btildebornon = .false.
            btildevirton = .false.
            btildecollon = .false.
            btilderealon = .false.
            phspcuts = .false.
            order = 1
         else
            if(powheginput('#qcd_order').eq.2d0) then
               order = 2
            elseif(powheginput('#qcd_order').eq.3d0) then
               order = 3
            endif
            btildebornon = .not.(powheginput("#btildeborn").eq.0)
            btildevirton = .not.(powheginput("#btildevirt").eq.0)
            btildecollon = .not.(powheginput("#btildecoll").eq.0)
            btilderealon = .not.(powheginput("#btildereal").eq.0)
            btildennloon = .not.(powheginput("#btildennlo").eq.0)
         endif
         ini = .false.
!         print*, btildebornon, btildevirton, btilderealon, btildecollon
!         stop
      endif
      btilde=0
      www=www0*hc2
      do j=1,ndiminteg-3
         xborn(j)=xx(j)
      enddo
      do j=1,3
         xrad(j)=xx(ndiminteg-3 + j)
      enddo
      if(ifirst.eq.0) then
         vbftot = 0d0
         dummy_analysis=.false.
         wwwtot=www
c     sets born momenta in kin. common block
         call reset_timer
         call gen_born_phsp(xborn)
         call born_suppression(suppfact)
c set scales
         call setscalesbtilde
         call allborn
         call get_timer(seconds)
         call addtocnt('born time (sec)',seconds)
c     Here we check whether or not events will pass the analysis. If
c     not, we don't bother computing the matrix element. Here we only
c     check for events with bornlike kinematics.
         if(phspcuts) then
            dummy_analysis = .true.
            passed_cuts_ig = .false.
            call analysis_extrainfo('born',flst_nborn,1d0,wwwtot)
            call analysis_driver(1d0,0)
            if(all(.not.passed_cuts_ig(0:2))) then ! In this case nothing is going to pass the cuts
               btildebornon = .false.
               btildevirton = .false.
               btildecollon = .false.
            else ! At least one passes the cut
               btildebornon = .true.
               btildevirton = .true..and.(.not.flg_bornonly)
               btildecollon = .true..and.(.not.flg_bornonly)
            endif
!     For NNLO part
!     Doesn't speed up things as the NNLO computation is super fast
!     compared to the subtraction terms, but we keep it for transparancy.
            passed_cuts_ig = .false.
            call analysis_extrainfo('nnlo',1,0d0,0d0)
            call analysis_driver(1d0,0)
            if(passed_cuts_ig(1)) then
               btildennloon = .true.
            else
               btildennloon = .false.
            endif
            
            dummy_analysis = .false. ! Next time analysis is called, it is for real
            vbftot = 0d0 ! Initialise
         endif
         call btildeborn(resborn)
         if(.not.btildebornon) resborn = 0
         if (.not.flg_bornonly.and..not.imode.eq.0) then
            call reset_timer
            if(btildevirton) then
               call btildevirt(resvirt)
            else
               resvirt = 0
            endif
            call get_timer(seconds)
            call addtocnt('virt time (sec)',seconds)
            call reset_timer
            if(btildecollon) then
               call btildecoll(xrad,rescoll,www)
            else
               rescoll = 0
            endif
            if(btilderealon) then
               call btildereal(xrad,resreal,www)
            else
               resreal = 0
            endif
            call get_timer(seconds)
            call addtocnt('real time (sec)',seconds)
         endif
c     accumulate values
         retval=0
         do j=1,flst_nborn
c     jacobians are already included in rescoll and resreal
            tmp=resborn(j)
            if (.not.flg_bornonly.and..not.imode.eq.0) then
               tmp = tmp   
     $              + resvirt(j)
     $              + rescoll(j) 
     $              + resreal(j)
            endif
c     initial value in results
            results(j)=tmp*www*suppfact
            retval=retval+tmp*www*suppfact
         enddo
      elseif(ifirst.eq.1) then
c     subsequent calls: In case of folding the call to btildeborn and
c     btildevirt can be avoided, since results are the same.  If the NLO
c     calculation is performed also (flg_nlotest is set) we need to
c     accumulate all weight within a single folding sequence in order to
c     later output the correct Born and Virtual contribution to the NLO
c     analysis routine.
         wwwtot=wwwtot+www
         if (.not.flg_bornonly.and..not.imode.eq.0) then
c btildecoll and btildereal take care themselves to invoke the NLO
c analysis if required.
            call reset_timer
c     in case btlscalereal is set, we need to reset the scales to the
c underlying Born value for the computation of the collinear remnants.
            call setscalesbtilde
            if(btildecollon) then
               call btildecoll(xrad,rescoll,www)
            else
               rescoll = 0
            endif
            if(btilderealon) then
               call btildereal(xrad,resreal,www)
            else
               resreal = 0
            endif
            call get_timer(seconds)
            call addtocnt('real time (sec)',seconds)
         endif
         retval=0
         do j=1,flst_nborn
            tmp=resborn(j)
            if (.not.flg_bornonly.and..not.imode.eq.0) then
               tmp = tmp   
     #              + resvirt(j)
     #              + rescoll(j) 
     #              + resreal(j)
            endif
c     accumulate values in results
            results(j)=results(j)+tmp*www*suppfact
            retval=retval+tmp*www*suppfact
         enddo
      elseif(ifirst.eq.2) then
         totborn=0d0
c compute Born, to return in retval0
         ptotborn=0
         do j=1,flst_nborn
            totborn=totborn+resborn(j)
            ptotborn=ptotborn+abs(resborn(j))
         enddo
         totborn=totborn*wwwtot
         ptotborn=ptotborn*wwwtot
         if (.not.pwhg_isfinite(totborn)) then 
            totborn = 0d0
            ptotborn = 0d0
            resborn = 0d0 
         endif
c         print*, 'phspcuts, imode', phspcuts, imode
c         stop
         if(flg_nlotest.or.phspcuts) then
c     output Born
            if(.not.imode.eq.0) then
               dummy_analysis = .false.
               call incl2pwhg(wwwtot*suppfact)
               call analysis_extrainfo('nnlo',1,0d0,0d0)
               call analysis_driver(1d0,0)
               call analysis_extrainfo('born',flst_nborn,resborn,wwwtot)
               call analysis_driver(totborn,0)
            endif
            if(.not.flg_bornonly.and..not.imode.eq.0) then
c     output virtual
               totvirt=0d0
               do j=1,flst_nborn
                  totvirt=totvirt+resvirt(j)
               enddo
               totvirt=totvirt*wwwtot
               if (.not.pwhg_isfinite(totvirt)) then 
                  totvirt = 0d0 
                  resvirt = 0d0 
               endif
               call analysis_extrainfo('virt',flst_nborn,resvirt,wwwtot)
               call analysis_driver(totvirt,0)
            endif
            call pwhgaccumup
         endif
c Make the born part of the result available; (to test, if bornonly is
c set, should equal the output btilde when ifirst=2)
         retval0=ptotborn*suppfact
         btilde=0
c closing call to end a sequence of correlated events in the analysis
c routines.  closing call: accumulate values with correct signs
         retval=0
         call adduptotals(results,flst_nborn)
         do j=1,flst_nborn
c this is only useful if withnegweights on (i.e. =1 in powheg.input,
c logical true here). However, better set a default (Les Houches
c interface will simply output this sign for the event.)
            rad_btilde_sign(j)=1
            if(flg_withnegweights) then
               if(results(j).lt.0) then
                  results(j)=-results(j)
                  rad_btilde_sign(j)=-1
               endif                  
            else
               if(results(j).lt.0) then
                  results(j)=0
               endif
            endif
            retval=retval+results(j)
c     Transfer all flavour components of btilde to the array in common
c     block; will be used to decide the underlying flavour of the event
            rad_btilde_arr(j)=results(j)
         enddo
c     AK replace retval0 and retval with corresponding values from analysis
         if(phspcuts) then
            retval0 = retval0
            retval = abs(vbftot)*suppfact
         endif
      else
         write(*,*) 'wrong value of ifirst in btilde => ',ifirst
         call exit(-1)
      endif
      end
      

      subroutine adduptotals(results,n)
      implicit none
      include 'nlegborn.h'
      integer n
      real * 8 results(n)
      real * 8 tot,totabs,totpos,totneg,etot,etotabs,etotpos,etotneg
      real * 8 totj(maxprocborn),totabsj(maxprocborn),
     1     totposj(maxprocborn),totnegj(maxprocborn),
     2     etotj(maxprocborn),etotabsj(maxprocborn),
     3     etotposj(maxprocborn),etotnegj(maxprocborn)
      integer nentries
      common/cadduptotals/tot,totabs,totpos,totneg,etot,etotabs,
     1     etotpos,etotneg,totj,totabsj,totposj,totnegj,
     1         etotj,etotabsj,etotposj,etotnegj,nentries
      real * 8 dtot,dtotabs,dtotpos,dtotneg
      integer j
      nentries=nentries+1
      dtot=0
      dtotabs=0
      dtotpos=0
      dtotneg=0
      do j=1,n
         dtot=dtot+results(j)
         dtotabs=dtotabs+abs(results(j))
         if(results(j).gt.0) then
            dtotpos=dtotpos+results(j)
         else
            dtotneg=dtotneg-results(j)
         endif
      enddo
      tot=tot+dtot
      totabs=totabs+dtotabs
      totpos=totpos+dtotpos
      totneg=totneg+dtotneg      
      etot=etot+dtot**2
      etotabs=etotabs+dtotabs**2
      etotpos=etotpos+dtotpos**2
      etotneg=etotneg+dtotneg**2
c j contributions
      do j=1,n
         dtot=results(j)
         dtotabs=abs(results(j))
         if(results(j).gt.0) then
            dtotpos=results(j)
         else
            dtotneg=-results(j)
         endif
         totj(j)=totj(j)+dtot
         totabsj(j)=totabsj(j)+dtotabs
         totposj(j)=totposj(j)+dtotpos
         totnegj(j)=totnegj(j)+dtotneg     
         etotj(j)=etotj(j)+dtot**2
         etotabsj(j)=etotabsj(j)+dtotabs**2
         etotposj(j)=etotposj(j)+dtotpos**2
         etotnegj(j)=etotnegj(j)+dtotneg**2
      enddo
      end

      subroutine resettotals
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      real * 8 tot,totabs,totpos,totneg,etot,etotabs,etotpos,etotneg
      real * 8 totj(maxprocborn),totabsj(maxprocborn),
     1     totposj(maxprocborn),totnegj(maxprocborn),
     2     etotj(maxprocborn),etotabsj(maxprocborn),
     3     etotposj(maxprocborn),etotnegj(maxprocborn)
      integer nentries
      common/cadduptotals/tot,totabs,totpos,totneg,etot,etotabs,
     1     etotpos,etotneg,totj,totabsj,totposj,totnegj,
     1         etotj,etotabsj,etotposj,etotnegj,nentries
      integer j
      nentries=0
      tot=0
      etot=0
      totabs=0
      etotabs=0
      totpos=0
      etotpos=0
      totneg=0
      etotneg=0
      do j=1,flst_nborn
         totj(j)=0
         etotj(j)=0
         totabsj(j)=0
         etotabsj(j)=0
         totposj(j)=0
         etotposj(j)=0
         totnegj(j)=0
         etotnegj(j)=0
      enddo
      end

      subroutine finaltotals
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_rad.h'
      real * 8 tot,totabs,totpos,totneg,etot,etotabs,etotpos,etotneg
      real * 8 totj(maxprocborn),totabsj(maxprocborn),
     1     totposj(maxprocborn),totnegj(maxprocborn),
     2     etotj(maxprocborn),etotabsj(maxprocborn),
     3     etotposj(maxprocborn),etotnegj(maxprocborn)
      integer nentries
      common/cadduptotals/tot,totabs,totpos,totneg,etot,etotabs,
     1     etotpos,etotneg,totj,totabsj,totposj,totnegj,
     1         etotj,etotabsj,etotposj,etotnegj,nentries
      integer n,j,k
      character * 80 format
      real * 8 tmp_totbtlj,tmp_etotbtlj,tmp_totabsbtlj,
     1     tmp_etotabsbtlj,tmp_totposbtlj,tmp_etotposbtlj,
     2     tmp_totnegbtlj,tmp_etotnegbtlj
      real * 8 tmp
      integer iun
      character * 20 pwgprefix
      integer lprefix
      common/cpwgprefix/pwgprefix,lprefix
      real * 8 powheginput
      external powheginput
      n=nentries
      rad_totbtl=tot/n
      rad_etotbtl=sqrt((etot/n-(tot/n)**2)/n)
      rad_totabsbtl=totabs/n
      rad_etotabsbtl=sqrt((etotabs/n-(totabs/n)**2)/n)
      rad_totposbtl=totpos/n
      rad_etotposbtl=sqrt((etotpos/n-(totpos/n)**2)/n)
      rad_totnegbtl=totneg/n
      rad_etotnegbtl=sqrt((etotneg/n-(totneg/n)**2)/n)
      write(*,*) 'tot:',rad_totbtl,'+-',rad_etotbtl
      write(*,*) 'abs:',rad_totabsbtl,'+-',rad_etotabsbtl
      write(*,*) 'pos:',rad_totposbtl,'+-',rad_etotposbtl
      write(*,*) 'neg:',rad_totnegbtl,'+-',rad_etotnegbtl
c
      if(powheginput('#ubsigmadetails').eq.1) then
         call newunit(iun)
         open(iun,file=pwgprefix(1:lprefix)//'ubsigma.dat')
         format='(      (i8,1x),4(a,1x,d10.4,a,d7.1))'
         write(format(2:4),'(i3)') nlegborn
         tmp=0
         do j=1,flst_nborn
            tmp_totbtlj=totj(j)/n
            tmp_etotbtlj=sqrt((etotj(j)/n-(totj(j)/n)**2)/n)
            tmp_totabsbtlj=totabsj(j)/n
            tmp_etotabsbtlj=sqrt((etotabsj(j)/n-(totabsj(j)/n)**2)/n)
            tmp_totposbtlj=totposj(j)/n
            tmp_etotposbtlj=sqrt((etotposj(j)/n-(totposj(j)/n)**2)/n)
            tmp_totnegbtlj=totnegj(j)/n
            tmp_etotnegbtlj=sqrt((etotnegj(j)/n-(totnegj(j)/n)**2)/n)
            write(iun,format) (flst_born(k,j),k=1,nlegborn),
     1           'tot:',tmp_totbtlj,' +- ',tmp_etotbtlj,
     2           '; abs:',tmp_totabsbtlj,' +- ',tmp_etotabsbtlj,
     3           '; pos:',tmp_totposbtlj,' +- ',tmp_etotposbtlj,
     4           '; neg:',tmp_totnegbtlj,' +- ',tmp_etotnegbtlj
            tmp=tmp+tmp_totbtlj
         enddo
         write(iun,*) tmp
      endif
      end

