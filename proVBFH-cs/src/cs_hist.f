c proVBFH-cs: adding the histograms of the inclusive and the exclusive part
c (the arrays of pwhg_bookhist-multi). cs_hist_save keeps the snapshot of
c the finished part (yhistarr2, errhistarr2: call pwhgsetout first),
c cs_hist_add adds it to the snapshot of the second part (errors in
c quadrature) once the second part's histograms are booked and finished;
c cs_hist_total gives the integral of histogram 1 (the cross section
c with the inclusive cuts of the analysis) for weight set k.
      subroutine cs_hist_save
      implicit none
      include 'pwhg_bookhist-multi.h'
      real * 8 ysave(maxmulti,0:maxbins+1,nhist),
     1     esave(maxmulti,0:maxbins+1,nhist)
      integer jsave
      common/cs_histsave/ysave,esave,jsave
      save /cs_histsave/
      jsave=jhist
      ysave=yhistarr2
      esave=errhistarr2
      end

      subroutine cs_hist_add
      implicit none
      include 'pwhg_bookhist-multi.h'
      real * 8 ysave(maxmulti,0:maxbins+1,nhist),
     1     esave(maxmulti,0:maxbins+1,nhist)
      integer jsave
      common/cs_histsave/ysave,esave,jsave
      save /cs_histsave/
      integer j,k
      if(jsave.ne.jhist) then
         write(*,*) 'cs_hist_add: different numbers of histograms'
         call exit(-1)
      endif
      do j=1,jhist
         do k=0,nbins(j)+1
            yhistarr2(1:nmulti,k,j)=yhistarr2(1:nmulti,k,j)
     1           +ysave(1:nmulti,k,j)
            errhistarr2(1:nmulti,k,j)=sqrt(errhistarr2(1:nmulti,k,j)**2
     1           +esave(1:nmulti,k,j)**2)
         enddo
      enddo
      end

      subroutine cs_hist_total(k,sigma,err)
      implicit none
      include 'pwhg_bookhist-multi.h'
      integer k
      real * 8 sigma,err
      integer i
      real * 8 w
      sigma=0
      err=0
      do i=1,nbins(1)
         w=xhistarr(i+1,1)-xhistarr(i,1)
         sigma=sigma+yhistarr2(k,i,1)*w
         err=err+(errhistarr2(k,i,1)*w)**2
      enddo
      err=sqrt(err)
      end
