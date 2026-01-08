c
c Example analysis for "p p > h h j j [QCD]" process.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine analysis_begin(nwgt,weights_info)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer nwgt
      character*(*) weights_info(*)
      integer i,kk,l,nwgt_analysis
      common/c_analysis/nwgt_analysis
      character*8 cc(2)
      data cc/'        ','vbfcuts '/
      double precision pi
      PARAMETER (PI=3.14159265358979312D0)
      include 'dbook.inc'
      real*8 vetomin, vetomax
      integer nbinveto
      common /to_veto_hist/vetomin,vetomax,nbinveto
c      
      call inihist
      vetomin = 0d0
      vetomax = 100d0
      nbinveto = 50

      nwgt_analysis=nwgt
      if (nwgt_analysis*50.gt.nplots/4) then
         write (*,*) 'error in analysis_begin: '/
     &        /'too many histograms, increase NPLOTS to',
     &        nwgt_analysis*50*4
         stop 1
      endif
      do kk=1,nwgt_analysis
      do i=1,2
        l=(kk-1)*50+(i-1)*25 !replace 25, 50, if total of hists changes
         call bookup(l+  1,'total '
     &       //cc(i)//weights_info(kk),1.0d0,0d0,1d0)

         call bookup(l+  2,'ptj1 '
     &       //cc(i)//weights_info(kk),20.0d0,0.d0,300.d0)
         call bookup(l+  3,'ptj2 '
     &       //cc(i)//weights_info(kk),10.0d0,0.d0,160.d0)
         call bookup(l+  4,'ptH1 '
     &       //cc(i)//weights_info(kk),30.0d0,0.d0,300.d0)
         call bookup(l+  5,'ptH2 '
     &       //cc(i)//weights_info(kk),30.0d0,0.d0,300.d0)
         call bookup(l+  6,'ptH_h '
     &       //cc(i)//weights_info(kk),30.0d0,0.d0,300.d0)
         call bookup(l+  7,'ptH_s '
     &        //cc(i)//weights_info(kk),30.0d0,0.d0,300.d0)
         
         call bookup(l+  8,'yj1 '
     &       //cc(i)//weights_info(kk),5d-1,-6.d0,6.d0)
         call bookup(l+  9,'yj2 '
     &       //cc(i)//weights_info(kk),5d-1,-6.d0,6.d0)
         call bookup(l+ 10,'yH1 '
     &       //cc(i)//weights_info(kk),1d-1,-5.d0,5.d0)
         call bookup(l+ 11,'yH2 '
     &       //cc(i)//weights_info(kk),1d-1,-5.d0,5.d0)
         call bookup(l+ 12,'yH_h '
     &       //cc(i)//weights_info(kk),1d-1,-5.d0,5.d0)
         call bookup(l+ 13,'yH_s '
     &        //cc(i)//weights_info(kk),1d-1,-5.d0,5.d0)
         
         call bookup(l+ 14,'mH1 '
     &       //cc(i)//weights_info(kk),200d0,0d0,2000d0)
         call bookup(l+ 15,'mH2 '
     &        //cc(i)//weights_info(kk),200d0,0d0,2000d0)
         
         call bookup(l+ 16,'DphiH1H2 '
     &       //cc(i)//weights_info(kk),pi/10d0,0d0,pi)
         call bookup(l+ 17,'yHH '
     &       //cc(i)//weights_info(kk),1d-1,-5d0,5d0)
         call bookup(l+ 18,'mHH '
     &       //cc(i)//weights_info(kk),200d0,0d0,2000d0)
         call bookup(l+ 19,'ptHH '
     &        //cc(i)//weights_info(kk),30.0d0,0.d0,300.d0)
         
         call bookup(l+ 20,'M_jj '
     &       //cc(i)//weights_info(kk),200d0,0d0,6000d0)
         call bookup(l+ 21,'Drapj1j2 '
     &       //cc(i)//weights_info(kk),5d-1,0d0,10d0)
         call bookup(l+ 22,'Dphij1j2 '
     &       //cc(i)//weights_info(kk),pi/10d0,0d0,pi)
         
         call bookup(l+ 23,'ptj3 '
     &        //cc(i)//weights_info(kk),5d0,0.d0,100.d0)
         call bookup(l+ 24,'yj3 '
     &       //cc(i)//weights_info(kk),5d-1,-6d0,6d0)
         call bookup(l+ 25,'y*j3 '
     &       //cc(i)//weights_info(kk),5d-1,-6d0,6d0)
         

      enddo
      enddo
      return
      end




cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine analysis_end(xnorm)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      character*14 ytit
      double precision xnorm
      integer i
      integer kk,l,nwgt_analysis
      common/c_analysis/nwgt_analysis
      include 'dbook.inc'
      call open_topdrawer_file
      call mclear
      do i=1,NPLOTS
         call mopera(i,'+',i,i,xnorm,0.d0)
         call mfinal(i)
      enddo
      ytit='sigma per bin '
      do kk=1,nwgt_analysis
      do i=1,2
        l=(kk-1)*50+(i-1)*25 ! replace 25 if total of hists changes
         call multitop(l+  1,3,2,'total rate       ',' ','LIN')

         call multitop(l+  2,3,2,'ptj1             ',' ','LIN')
         call multitop(l+  3,3,2,'ptj2             ',' ','LIN')
         call multitop(l+  4,3,2,'ptH1             ',' ','LIN')
         call multitop(l+  5,3,2,'ptH2             ',' ','LIN')
         call multitop(l+  6,3,2,'ptH_h            ',' ','LIN')
         call multitop(l+  7,3,2,'ptH_s            ',' ','LIN')
         
         call multitop(l+  8,3,2,'yj1              ',' ','LIN')
         call multitop(l+  9,3,2,'yj2              ',' ','LIN')
         call multitop(l+ 10,3,2,'yH1              ',' ','LIN')
         call multitop(l+ 11,3,2,'yH2              ',' ','LIN')
         call multitop(l+ 12,3,2,'yH_h             ',' ','LIN')
         call multitop(l+ 13,3,2,'yH_s             ',' ','LIN')

         call multitop(l+ 14,3,2,'mH1              ',' ','LIN')
         call multitop(l+ 15,3,2,'mH2              ',' ','LIN')

         call multitop(l+ 16,3,2,'DphiH1H2         ',' ','LIN')
         call multitop(l+ 17,3,2,'yHH              ',' ','LIN')
         call multitop(l+ 18,3,2,'mHH              ',' ','LIN')
         call multitop(l+ 19,3,2,'ptHH             ',' ','LIN')

         call multitop(l+ 20,3,2,'M_jj             ',' ','LIN')
         call multitop(l+ 21,3,2,'Drapj1j2         ',' ','LIN')
         call multitop(l+ 22,3,2,'Dphij1j2         ',' ','LIN')

         call multitop(l+ 23,3,2,'ptj3             ',' ','LIN')
         call multitop(l+ 24,3,2,'yj3              ',' ','LIN')
         call multitop(l+ 25,3,2,'y*j3             ',' ','LIN')
      enddo
      enddo
      call close_topdrawer_file
      return                
      end




cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine analysis_fill(p,istatus,ipdg,wgts,ibody)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      include 'nexternal.inc'
      integer istatus(nexternal)
      integer iPDG(nexternal)
      double precision p(0:4,nexternal)
      double precision wgts(*)
      integer ibody
      double precision wgt,var
      integer i,j,k,kk,l,nwgt_analysis
      common/c_analysis/nwgt_analysis
      double precision www,pQCD(0:3,nexternal),palg,rfj,sycut,yjmax
     $     ,pjet(0:3,nexternal),ptjet(nexternal),yjet(nexternal)
     $     ,ptj_tag,deltay12,mj1j2min
     $     ,ph1(0:3),ph2(0:3),mH1,mH2,mHH
     $     ,pj1(0:3),pj2(0:3),pj3(0:3),pjj(0:3),pHH(0:3)
     $     ,ptH1,yH1,ptH2,yH2,ptHd,yHd,ptHs,yHs,ptHH,yHH,phiHH
     $     ,ptj1,yj1,ptj2,yj2
     $     ,Dphij1j2,mj1j2,Dyj1j2,ptj3,yj3,yj3rel
      logical pass_tag_cuts,flag
      integer nQCD,jet(nexternal),ij1y,ij2y,ij3y,njet,njety
     $     ,ij1,ij2,ij3
      double precision getptv4,getrapidityv4,getpseudorapv4,getdelphiv4
     $     ,getdrv4,getinvmv4,getmod
      external getptv4,getrapidityv4,getpseudorapv4,getdelphiv4,getdrv4
     $     ,getinvmv4,getmod
      real*8 vetomin, vetomax
      integer nbinveto
      common /to_veto_hist/vetomin,vetomax,nbinveto
      real*8 xd(1:3)
      data (xd(i),i=1,3)/0d0,0d0,1d0/
      if (nexternal.ne.7) then
         write (*,*) 'error #1 in analysis_fill: '/
     &        /'only for process "p p > h h j j [QED]"'
         stop 1
      endif
      if (.not.(abs(ipdg(1)).le.5 .or. ipdg(1).eq.21.
     &     or. ipdg(1).eq.22)) then
         write (*,*) 'error #2 in analysis_fill: '/
     &        /'only for process "p p > h h j j [QED]"'
         stop 1
      endif
      if (.not.(abs(ipdg(2)).le.5 .or. ipdg(2).eq.21.
     &     or. ipdg(2).eq.22)) then
         write (*,*) 'error #3 in analysis_fill: '/
     &        /'only for process "p p > h h j j [QED]"'
         stop 1
      endif
      if (.not.(abs(ipdg(5)).le.5 .or. ipdg(5).eq.21.
     &     or. ipdg(5).eq.22)) then
         write (*,*) 'error #4 in analysis_fill: '/
     &        /'only for process "p p > h h j j [QED]"'
         stop 1
      endif
      if (.not.(abs(ipdg(6)).le.5 .or. ipdg(6).eq.21.
     &     or. ipdg(6).eq.22)) then
         write (*,*) 'error #5 in analysis_fill: '/
     &        /'only for process "p p > h h j j [QED]"'
         stop 1
      endif
      if (.not.(abs(ipdg(7)).le.5 .or. ipdg(7).eq.21.
     &     or. ipdg(7).eq.22)) then
         write (*,*) 'error #6 in analysis_fill: '/
     &        /'only for process "p p > h h j j [QED]"'
         stop 1
      endif
      if (ipdg(3).ne.25) then
         write (*,*) 'error #7 in analysis_fill: '/
     &        /'only for process "p p > h h j j [QED]"'
         stop 1
      endif
      if (ipdg(4).ne.25) then
         write (*,*) 'error #8 in analysis_fill: '/
     &        /'only for process "p p > h h j j [QED]"'
         stop 1
      endif
c

c Put all (light) QCD partons in momentum array for jet clustering.
      nQCD=0
      do j=nincoming+2,nexternal
         if (abs(ipdg(j)).le.5 .or. ipdg(j).eq.21.or.ipdg(j).eq.22) then
            nQCD=nQCD+1
            do i=0,3
               pQCD(i,nQCD)=p(i,j) 
            enddo
         endif
      enddo
      
C---CLUSTER THE EVENT
      palg  = -1.d0
      rfj   = 0.4d0
      sycut = 25d0
      yjmax = 4.5d0
      do i=1,nexternal
         do j=0,3
            pjet(j,i)=0d0
         enddo
         ptjet(i)=0d0
         yjet(i)=0d0
         jet(i)=0
      enddo
      ij1y=0
      ij2y=0
      ij3y=0
      njet=0
      njety=0



c******************************************************************************
c     call FASTJET to get all the jets
c
c     INPUT:
c     input momenta:               pQCD(0:3,nexternal), energy is 0th component
c     number of input momenta:     nQCD
c     radius parameter:            rfj
c     minumum jet pt:              sycut
c     jet algorithm:               palg, 1.0=kt, 0.0=C/A, -1.0 = anti-kt
c
c     OUTPUT:
c     jet momenta:                             pjet(0:3,nexternal), E is 0th cmpnt
c     the number of jets (with pt > SYCUT):    njet
c     the jet for a given particle 'i':        jet(i),   note that this is
c     the particle in pQCD, which doesn't necessarily correspond to the particle
c     label in the process
c
      call amcatnlo_fastjetppgenkt(pQCD,nQCD,rfj,sycut,palg,pjet,njet
     $     ,jet)
c
c******************************************************************************
      do i=1,njet
         ptjet(i)=getptv4(pjet(0,i))
         if(i.gt.1)then
            if (ptjet(i).gt.ptjet(i-1)) then
               write (*,*) "Error 1: jets should be ordered in pt"
               stop
            endif
         endif
         yjet(i)=getrapidityv4(pjet(0,i))
C now look for jets within the rapidity cuts
         if (dabs(yjet(i)).lt.yjmax) then
            njety=njety+1
            if (ij1y.eq.0) then
               ij1y = i
            else if (ij2y.eq.0) then
               ij2y = i
            else if (ij3y.eq.0) then
               ij3y = i
            endif
         endif
      enddo

c Nason-Oleari cuts (hep-ph/0911.5299)
      ptj_tag  = 25d0
      deltay12 = 4.5d0
      mj1j2min = 600d0

c this is the loop for w-o / w vbf cuts
      do i=1,2
      if(i.eq.1) then 
         ij1 = 1
         ij2 = 2
         ij3 = 3
      endif
      if(i.eq.2) then
         njet = njety
         ij1 = ij1y
         ij2 = ij2y
         ij3 = ij3y
      endif

c Load momenta
      do k=0,3
         pH1(k)  =p(k,3)
         pH2(k)  =p(k,4)
         pHH(k)  =p(k,3) + p(k,4)
         pj1(k)  =pjet(k,ij1)
         pj2(k)  =pjet(k,ij2)
         pj3(k)  =pjet(k,ij3)
         pjj(k)  =pjet(k,ij1)+pjet(k,ij2)
      enddo

c Define observables
c Higgs
      ptH1     = getptv4(pH1)
      yH1      = getrapidityv4(pH1)
      ptH2     = getptv4(pH2)
      yH2      = getrapidityv4(pH2)
      ptHH     = getptv4(pHH)
      yHH      = getrapidityv4(pHH)
      phiHH    = getdelphiv4(pH1,pH2)
      if (ptH1.gt.ptH2) then
         ptHd = ptH1
         yHd  = yH1
         ptHs = ptH2
         yHs  = yH2
      else
         ptHd = ptH2
         yHd  = yH2
         ptHs = ptH1
         yHs  = yH1
      endif
c At least one jet
      if(njet.ge.1)then
        ptj1    = getptv4(pj1)
        yj1     = getrapidityv4(pj1)
      endif
c At least two jets
      if(njet.ge.2)then
        ptj2    = getptv4(pj2)
        yj2     = getrapidityv4(pj2)
        Dphij1j2= getdelphiv4(pj1,pj2)
        mj1j2   = getinvmv4(pjj)
        Dyj1j2  = abs(yj1-yj2)
      endif
c At least three jets
      if(njet.ge.3)then
        ptj3    = getptv4(pj3)
        yj3     = getrapidityv4(pj3)
        yj3rel  = yj3-(yj1+yj2)/2d0
      endif
c

      pass_tag_cuts = njety.ge.2 .and.
     &                ptj1.ge.ptj_tag .and.
     &                ptj2.ge.ptj_tag .and.
     &                abs(yj1-yj2).ge.deltay12 .and.
     &                yj1*yj2.le.0d0 .and.
     &                mj1j2.ge.mj1j2min 

      if(i.eq.1) then 
         flag=.true.
      endif

      if(i.eq.2) then
         flag=pass_tag_cuts
      endif



      do kk=1,nwgt_analysis
      www=wgts(kk)
      l=(kk-1)*50+(i-1)*25 ! replace 50 and 25 if number of hists changes
      if(flag)then
         call mfill(l+  1,0.5d0, www)
         call mfill(l+  4,ptH1,  www)
         call mfill(l+  5,ptH2,  www)
         call mfill(l+  6,ptHd, www)
         call mfill(l+  7,ptHs, www)
         call mfill(l+ 10,yH1,   www)
         call mfill(l+ 11,yH2,   www)
         call mfill(l+ 12,yHd,  www)
         call mfill(l+ 13,yHs,  www)
         call mfill(l+ 14,mH1,   www)
         call mfill(l+ 15,mH2,   www)
         call mfill(l+ 16,phiHH, www)
         call mfill(l+ 17,yHH,   www)
         call mfill(l+ 18,mHH,   www)
         call mfill(l+ 19,ptHH,  www)

         if(njet.ge.1)then
            call mfill(l+  2,ptj1,  www)
            call mfill(l+  8,yj1,   www)
         endif 
         
         if(njet.ge.2)then
            call mfill(l+  3,ptj2,  www)
            call mfill(l+  9,yj2,   www)
            call mfill(l+ 20,mj1j2, www)
            call mfill(l+ 21,Dyj1j2,www)
            call mfill(l+ 22,Dphij1j2,www)
         endif
         
         if(njet.ge.3)then
            call mfill(l+ 23,ptj3,  www)
            call mfill(l+ 24,yj3,   www)
            call mfill(l+ 25,yj3rel,www)
         endif
      endif
      enddo
      enddo

 999  return      
      end


      function getrapidity(en,pl)
      implicit none
      real*8 getrapidity,en,pl,tiny,xplus,xminus,y
      parameter (tiny=1.d-8)
c
      xplus=en+pl
      xminus=en-pl
      if(xplus.gt.tiny.and.xminus.gt.tiny)then
        if( (xplus/xminus).gt.tiny.and.(xminus/xplus).gt.tiny )then
          y=0.5d0*log( xplus/xminus )
        else
          y=sign(1.d0,pl)*1.d8
        endif
      else
        y=sign(1.d0,pl)*1.d8
      endif
      getrapidity=y
      return
      end


      function getpseudorap(en,ptx,pty,pl)
      implicit none
      real*8 getpseudorap,en,ptx,pty,pl,tiny,pt,eta,th
      parameter (tiny=1.d-5)
c
      pt=sqrt(ptx**2+pty**2)
      if(pt.lt.tiny.and.abs(pl).lt.tiny)then
        eta=sign(1.d0,pl)*1.d8
      else
        th=atan2(pt,pl)
        eta=-log(tan(th/2.d0))
      endif
      getpseudorap=eta
      return
      end


      function getinvm(en,ptx,pty,pl)
      implicit none
      real*8 getinvm,en,ptx,pty,pl,tiny,tmp
      parameter (tiny=1.d-5)
c
      tmp=en**2-ptx**2-pty**2-pl**2
      if(tmp.gt.0.d0)then
        tmp=sqrt(tmp)
      elseif(tmp.gt.-tiny)then
        tmp=0.d0
      else
        write(*,*)'Attempt to compute a negative mass'
        stop
      endif
      getinvm=tmp
      return
      end


      function getdelphi(ptx1,pty1,ptx2,pty2)
      implicit none
      real*8 getdelphi,ptx1,pty1,ptx2,pty2,tiny,pt1,pt2,tmp
      parameter (tiny=1.d-5)
c
      pt1=sqrt(ptx1**2+pty1**2)
      pt2=sqrt(ptx2**2+pty2**2)
      if(pt1.ne.0.d0.and.pt2.ne.0.d0)then
        tmp=ptx1*ptx2+pty1*pty2
        tmp=tmp/(pt1*pt2)
        if(abs(tmp).gt.1.d0+tiny)then
          write(*,*)'Cosine larger than 1'
          stop
        elseif(abs(tmp).ge.1.d0)then
          tmp=sign(1.d0,tmp)
        endif
        tmp=acos(tmp)
      else
        tmp=1.d8
      endif
      getdelphi=tmp
      return
      end


      function getdr(en1,ptx1,pty1,pl1,en2,ptx2,pty2,pl2)
      implicit none
      real*8 getdr,en1,ptx1,pty1,pl1,en2,ptx2,pty2,pl2,deta,dphi,
     # getpseudorap,getdelphi
c
      deta=getpseudorap(en1,ptx1,pty1,pl1)-
     #     getpseudorap(en2,ptx2,pty2,pl2)
      dphi=getdelphi(ptx1,pty1,ptx2,pty2)
      getdr=sqrt(dphi**2+deta**2)
      return
      end


      function getdry(en1,ptx1,pty1,pl1,en2,ptx2,pty2,pl2)
      implicit none
      real*8 getdry,en1,ptx1,pty1,pl1,en2,ptx2,pty2,pl2,deta,dphi,
     # getrapidity,getdelphi
c
      deta=getrapidity(en1,pl1)-
     #     getrapidity(en2,pl2)
      dphi=getdelphi(ptx1,pty1,ptx2,pty2)
      getdry=sqrt(dphi**2+deta**2)
      return
      end


      function getptv(p)
      implicit none
      real*8 getptv,p(5)
c
      getptv=sqrt(p(1)**2+p(2)**2)
      return
      end


      function getpseudorapv(p)
      implicit none
      real*8 getpseudorapv,p(5)
      real*8 getpseudorap
c
      getpseudorapv=getpseudorap(p(4),p(1),p(2),p(3))
      return
      end


      function getrapidityv(p)
      implicit none
      real*8 getrapidityv,p(5)
      real*8 getrapidity
c
      getrapidityv=getrapidity(p(4),p(3))
      return
      end


      function getdrv(p1,p2)
      implicit none
      real*8 getdrv,p1(5),p2(5)
      real*8 getdr
c
      getdrv=getdr(p1(4),p1(1),p1(2),p1(3),
     #             p2(4),p2(1),p2(2),p2(3))
      return
      end


      function getinvmv(p)
      implicit none
      real*8 getinvmv,p(5)
      real*8 getinvm
c
      getinvmv=getinvm(p(4),p(1),p(2),p(3))
      return
      end


      function getdelphiv(p1,p2)
      implicit none
      real*8 getdelphiv,p1(5),p2(5)
      real*8 getdelphi
c
      getdelphiv=getdelphi(p1(1),p1(2),
     #                     p2(1),p2(2))
      return
      end


      function getptv4(p)
      implicit none
      real*8 getptv4,p(0:3)
c
      getptv4=sqrt(p(1)**2+p(2)**2)
      return
      end


      function getpseudorapv4(p)
      implicit none
      real*8 getpseudorapv4,p(0:3)
      real*8 getpseudorap
c
      getpseudorapv4=getpseudorap(p(0),p(1),p(2),p(3))
      return
      end


      function getrapidityv4(p)
      implicit none
      real*8 getrapidityv4,p(0:3)
      real*8 getrapidity
c
      getrapidityv4=getrapidity(p(0),p(3))
      return
      end


      function getdrv4(p1,p2)
      implicit none
      real*8 getdrv4,p1(0:3),p2(0:3)
      real*8 getdr
c
      getdrv4=getdr(p1(0),p1(1),p1(2),p1(3),
     #              p2(0),p2(1),p2(2),p2(3))
      return
      end


      function getinvmv4(p)
      implicit none
      real*8 getinvmv4,p(0:3)
      real*8 getinvm
c
      getinvmv4=getinvm(p(0),p(1),p(2),p(3))
      return
      end


      function getdelphiv4(p1,p2)
      implicit none
      real*8 getdelphiv4,p1(0:3),p2(0:3)
      real*8 getdelphi
c
      getdelphiv4=getdelphi(p1(1),p1(2),
     #                      p2(1),p2(2))
      return
      end


      function getcosv4(q1,q2)
      implicit none
      real*8 getcosv4,q1(0:3),q2(0:3)
      real*8 xnorm1,xnorm2,tmp
c
      if(q1(0).lt.0.d0.or.q2(0).lt.0.d0)then
        getcosv4=-1.d10
        return
      endif
      xnorm1=sqrt(q1(1)**2+q1(2)**2+q1(3)**2)
      xnorm2=sqrt(q2(1)**2+q2(2)**2+q2(3)**2)
      if(xnorm1.lt.1.d-6.or.xnorm2.lt.1.d-6)then
        tmp=-1.d10
      else
        tmp=q1(1)*q2(1)+q1(2)*q2(2)+q1(3)*q2(3)
        tmp=tmp/(xnorm1*xnorm2)
        if(abs(tmp).gt.1.d0.and.abs(tmp).le.1.001d0)then
          tmp=sign(1.d0,tmp)
        elseif(abs(tmp).gt.1.001d0)then
          write(*,*)'Error in getcosv4',tmp
          stop
        endif
      endif
      getcosv4=tmp
      return
      end



      function getmod(p)
      implicit none
      double precision p(0:3),getmod

      getmod=sqrt(p(1)**2+p(2)**2+p(3)**2)

      return
      end



      subroutine getperpenv4(q1,q2,qperp)
c Normal to the plane defined by \vec{q1},\vec{q2}
      implicit none
      real*8 q1(0:3),q2(0:3),qperp(0:3)
      real*8 xnorm1,xnorm2
      integer i
c
      xnorm1=sqrt(q1(1)**2+q1(2)**2+q1(3)**2)
      xnorm2=sqrt(q2(1)**2+q2(2)**2+q2(3)**2)
      if(xnorm1.lt.1.d-6.or.xnorm2.lt.1.d-6)then
        do i=1,4
          qperp(i)=-1.d10
        enddo
      else
        qperp(1)=q1(2)*q2(3)-q1(3)*q2(2)
        qperp(2)=q1(3)*q2(1)-q1(1)*q2(3)
        qperp(3)=q1(1)*q2(2)-q1(2)*q2(1)
        do i=1,3
          qperp(i)=qperp(i)/(xnorm1*xnorm2)
        enddo
        qperp(0)=1.d0
      endif
      return
      end




      subroutine getwedge(p1,p2,pout)
      implicit none
      real*8 p1(0:3),p2(0:3),pout(0:3)

      pout(1)=p1(2)*p2(3)-p1(3)*p2(2)
      pout(2)=p1(3)*p2(1)-p1(1)*p2(3)
      pout(3)=p1(1)*p2(2)-p1(2)*p2(1)
      pout(0)=0d0

      return
      end

