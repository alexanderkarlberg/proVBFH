      subroutine setborn(p,bflav,born,bornjk,bmunu,bornsub)

      implicit none
      include 'nlegborn.h'
      integer nlegs
      parameter (nlegs=nlegborn)
      real * 8 p(0:3,nlegs),bornjk(nlegs,nlegs), borncol(2)
      real * 8 bmunu(0:3,0:3,nlegs),bbmunu(0:3,0:3),born,colcf
      integer bflav(nlegs), bflav_mg(nlegs)
      integer j,k,mu,nu
      real * 8 bornsub(2)
c
cccccccccccccccccccccccccccccccccccccccccccccccc
c
c ordering of momenta is q(1) q(2) -> H(3) q(4)q(5)"g(6)"
c
      call calc_als
      call ctrans(1)

      call compborn_hhjjj(p,bflav,born,bmunu,bornjk,borncol) 
      
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine compborn_hhjjj(pin,bflav,born,bmunu,bornjk, born_col)
      implicit none
c
      include 'nlegborn.h'
      include 'pwhg_math.h'
      include 'pwhg_flst.h'
      include 'pwhg_flst_2.h'
      include 'tags.h' 
c
      integer nlegs
      parameter (nlegs=nlegborn)
      real*8 pin(0:3,nlegs)  
      integer bflav(nlegs)
      real*8 born, borncol(2)
      real*8 born_col(1:2), bbmunu(0:3,0:3)
      real*8 bmunu(0:3,0:3,nlegs),bornjk(nlegs,nlegs)
      integer j,jj
c
c vbfnlo stuff:
      include 'global_col.inc'
      real*8 p(0:3,np), v(0:3,nv)
      real*8 pbar(0:3,nlegborn),qbar(0:4)
      real*8 polcol,polcolq,polcolg
      real*8 res,resc(3),ress(3)
      real*8 qvec(0:3,2), cqvec(0:3,2)

      real*8 N ! color factors
      parameter(N=3d0)

      complex*16 zero
      parameter (zero=(0d0,0d0))
c
c declare local variables
c
      integer i,mu,nu
      integer FSIGN(nlegborn),gsign,physToDiag(5)
      
      real*8 nans(2,2,3),cans(2,3)
      logical nc_type
      integer k,icc

      integer ftype(1:7)
      integer j1,j2,j3,j4 !identify fermion color structures
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      j1 = 0
      j2 = 0
      j3 = 0
      j4 = 0

      polcol = 0d0
      polcolq = 1d0/(4d0*N**2)
      polcolg = 1d0/(4d0*N*(N**2-1))

      ftype(1:7) = 0
      do mu = 0,3
c fermions:
         p(mu,1) = pin(mu,1)
         p(mu,2) = pin(mu,2)
         
      ! VBFHHMOD: 4->5, 5->6, 6->7
         if (bflav(1)*bflav(2).ne.0) then ! final-state gluon 
            if (bflav(7).eq.0) then
               p(mu,3) = pin(mu,5)
               p(mu,4) = pin(mu,6) 
               p(mu,5) = pin(mu,7)  !gluon
            elseif (bflav(5).eq.0) then
               p(mu,3) = pin(mu,5)
               p(mu,4) = pin(mu,6) 
               p(mu,5) = pin(mu,7) 
               print*, 'This should not happen'
            elseif (bflav(6).eq.0) then
               p(mu,3) = pin(mu,5)
               p(mu,4) = pin(mu,6) 
               p(mu,5) = pin(mu,7) 
               print*, 'This should not happen'
            endif
         else   ! initial-state gluon 
            p(mu,3) = pin(mu,5)
            p(mu,4) = pin(mu,6) 
            p(mu,5) = pin(mu,7) 
         endif ! fin/in state gluon
      enddo ! mu

      if (bflav(1).gt.0.and.bflav(2).gt.0) then

C*******************  q1 q3 ---> q2 q4 g H H  **********************

      physToDiag(1)=1    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=3
      physToDiag(3)=2    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=4
      physToDiag(5)=5   ! gluon

      fsign(1) = 1
      fsign(2) = 1
      fsign(3) = 1
      fsign(4) = 1
      gsign    = 1

      polcol = polcolq

c up- or down-type quark:      
      ! VBFHHMOD: bflav(4)->bflav(5),ftype(4)->ftype(5), ... 
      ftype(1) = 2-mod(abs(bflav(1)),2)         
      ftype(2) = 2-mod(abs(bflav(2)),2)      
      ftype(5) = 2-mod(abs(bflav(5)),2)   
      ftype(6) = 2-mod(abs(bflav(6)),2)

      if ((ftype(1).eq.ftype(5))) then
         k = -2*ftype(1)-ftype(2)+7
      else !cc 
         icc = (-4*(bflav(1)/abs(bflav(1)))+6)/2
         k = 7-ftype(icc)
      endif   

      ! VBFHHMOD: 4->5, 5->6, 6->7
      j1 = 1
      j2 = 5
      j3 = 2
      j4 = 6


      elseif (bflav(1).gt.0.and.bflav(2).lt.0) then
c            
C******************* q1 qb4 ---> q2 qb3 g H H  **********************
      
      physToDiag(1)=1    
      physToDiag(2)=4
      physToDiag(3)=2    
      physToDiag(4)=3
      physToDiag(5)=5   ! gluon
c
      fsign(1) =  1
      fsign(2) =  1
      fsign(3) = -1
      fsign(4) = -1
      gsign    =  1
      
      polcol = polcolq
      
c up- or down-type quark:      
      ! VBFHHMOD: bflav(4)->bflav(5), ... 
      ftype(1) = 2-mod(abs(bflav(1)),2)         
      ftype(2) = 2-mod(abs(bflav(2)),2)      
      ftype(5) = 2-mod(abs(bflav(5)),2)   
      ftype(6) = 2-mod(abs(bflav(6)),2)

      if ((ftype(1).eq.ftype(5))) then !nc
         k = -2*ftype(1)-ftype(2)+7
      else !cc 
         icc = (-4*(bflav(1)/abs(bflav(1)))+6)/2
         k = 7-ftype(icc)
      endif !nc/cc  

      ! VBFHHMOD: 4->5, 5->6, 6->7
      j1 = 1
      j2 = 5
      j3 = 2!5
      j4 = 6!2

      elseif (bflav(1).lt.0.and.bflav(2).gt.0) then
      
C******************* qbar2 q3 ---> qbar1 q4 g H H  **********************
      
      physToDiag(1)=2    
      physToDiag(2)=3
      physToDiag(3)=1    
      physToDiag(4)=4
      physToDiag(5)=5   ! gluon
c
      fsign(1) = -1
      fsign(2) = -1
      fsign(3) =  1
      fsign(4) =  1
      gsign    =  1
      
      polcol = polcolq

c up- or down-type quark:       
      ! VBFHHMOD: bflav(4)->bflav(5), ...
      ftype(1) = 2-mod(abs(bflav(1)),2)         
      ftype(2) = 2-mod(abs(bflav(2)),2)      
      ftype(5) = 2-mod(abs(bflav(5)),2)   
      ftype(6) = 2-mod(abs(bflav(6)),2)

      if ((ftype(1).eq.ftype(5))) then !nc
         k = -2*ftype(1)-ftype(2)+7
      else !cc 
         icc = (-4*(bflav(1)/abs(bflav(1)))+6)/2
         k = 7-ftype(icc)
      endif !nc/cc  

      ! VBFHHMOD: 4->5, 5->6, 6->7
      j1 = 1!4
      j2 = 5!1 
      j3 = 2 
      j4 = 6

      elseif (bflav(1).lt.0.and.bflav(2).lt.0) then

C*******************  qbar2 qb4 ---> qbar1 qb3 g H H  **********************

      physToDiag(1)=2    
      physToDiag(2)=4
      physToDiag(3)=1    
      physToDiag(4)=3
      physToDiag(5)=5   ! gluon
c
      fsign(1) = -1
      fsign(2) = -1
      fsign(3) = -1
      fsign(4) = -1
      gsign    =  1
      
      polcol = polcolq

c up- or down-type quark:  
      ! VBFHHMOD: bflav(4)->bflav(5), ...
      ftype(1) = 2-mod(abs(bflav(1)),2)         
      ftype(2) = 2-mod(abs(bflav(2)),2)      
      ftype(5) = 2-mod(abs(bflav(5)),2)   
      ftype(6) = 2-mod(abs(bflav(6)),2)

      if ((ftype(1).eq.ftype(5))) then !nc
         k = -2*ftype(1)-ftype(2)+7
      else !cc 
         icc = (-4*(bflav(1)/abs(bflav(1)))+6)/2
         k = 7-ftype(icc)
      endif !nc/cc  

      ! VBFHHMOD: 4->5, 5->6, 6->7
      j1 = 1!4
      j2 = 5!1
      j3 = 2!5
      j4 = 6!2

      elseif (bflav(1).eq.0.and.bflav(2).gt.0.and.bflav(5).gt.0) then
         
c*******************  g q ---> q q qb H H  **********************

      physToDiag(1)=5          
      physToDiag(2)=3           
      physToDiag(3)=2           
      physToDiag(4)=4
      physToDiag(5)=1
c
      fsign(1) = -1
      fsign(2) =  1
      fsign(3) =  1
      fsign(4) =  1
      gsign    = -1
      
      polcol = polcolg  

c up- or down-type quark:   
      ! VBFHHMOD: bflav(4)->bflav(5), ...   
      ftype(1) = 2-mod(abs(bflav(7)),2)         
      ftype(2) = 2-mod(abs(bflav(2)),2)      
      ftype(5) = 2-mod(abs(bflav(5)),2)   
      ftype(6) = 2-mod(abs(bflav(6)),2)

      if ((ftype(1).eq.ftype(5))) then !nc
         k = -2*ftype(1)-ftype(2)+7
      else !cc
         icc = (4*(bflav(2)/abs(bflav(2)))+8)/2
         k = 7-ftype(icc)
      endif !nc/cc  

      ! VBFHHMOD: 4->5, 5->6, 6->7
      j1 = 5!6
      j2 = 7!4
      j3 = 2
      j4 = 6

      elseif (bflav(1).eq.0.and.bflav(2).gt.0.and.bflav(5).lt.0) then
         
c*******************  g q ---> qb q q H H  **********************
c
      physToDiag(1)=5          
      physToDiag(2)=3           
      physToDiag(3)=1           
      physToDiag(4)=4
      physToDiag(5)=2
c
      fsign(1) = -1
      fsign(2) =  1
      fsign(3) =  1
      fsign(4) =  1
      gsign    = -1
      
      polcol = polcolg

c up- or down-type quark:   
      ! VBFHHMOD: bflav(4)->bflav(5), ...   
      ftype(1) = 2-mod(abs(bflav(7)),2)         
      ftype(2) = 2-mod(abs(bflav(2)),2)      
      ftype(5) = 2-mod(abs(bflav(5)),2)   
      ftype(6) = 2-mod(abs(bflav(6)),2)     

      if ((ftype(1).eq.ftype(5))) then !nc
         k = -2*ftype(1)-ftype(2)+7
      else !cc 
         icc = (4*(bflav(2)/abs(bflav(2)))+8)/2
         k = 7-ftype(icc)
      endif !nc/cc  

      ! VBFHHMOD: 4->5, 5->6, 6->7
      j1 = 5
      j2 = 7
      j3 = 2
      j4 = 6

      elseif (bflav(1).gt.0.and.bflav(2).eq.0.
     &    and.bflav(5).gt.0.and.bflav(6).gt.0) then
      
C*******************  q g ---> q q qb H H  **********************
      
      physToDiag(1)=1             
      physToDiag(2)=5
      physToDiag(3)=2               
      physToDiag(4)=4
      physToDiag(5)=3 
c
      fsign(1) =  1
      fsign(2) =  1
      fsign(3) = -1
      fsign(4) =  1
      gsign    = -1

      polcol = polcolg
c
c up- or down-type quark:     
      ! VBFHHMOD: bflav(4)->bflav(5), ... 
      ftype(1) = 2-mod(abs(bflav(1)),2)         
      ftype(2) = 2-mod(abs(bflav(7)),2)      
      ftype(5) = 2-mod(abs(bflav(5)),2)   
      ftype(6) = 2-mod(abs(bflav(6)),2)

      if ((ftype(1).eq.ftype(5))) then !nc
         k = -2*ftype(1)-ftype(2)+7
      else !cc 
         icc = (-4*(bflav(1)/abs(bflav(1)))+6)/2
         k = 7-ftype(icc)
      endif !nc/cc  

      ! VBFHHMOD: 4->5, 5->6, 6->7
      j1 = 1
      j2 = 5
      j3 = 6!6
      j4 = 7!5

      elseif (bflav(1).gt.0.and.bflav(2).eq.0.and.
     &        bflav(5).gt.0.and.bflav(6).lt.0) then
      
C*******************  q g ---> q qb q H H  **********************

      physToDiag(1)=1             
      physToDiag(2)=5            
      physToDiag(3)=2                  
      physToDiag(4)=3
      physToDiag(5)=4
c
      fsign(1) =  1
      fsign(2) =  1
      fsign(3) = -1
      fsign(4) =  1
      gsign    = -1

      polcol = polcolg

c up- or down-type quark:      
      ! VBFHHMOD: bflav(4)->bflav(5), ...
      ftype(1) = 2-mod(abs(bflav(1)),2)         
      ftype(2) = 2-mod(abs(bflav(7)),2)      
      ftype(5) = 2-mod(abs(bflav(5)),2)   
      ftype(6) = 2-mod(abs(bflav(6)),2)

      if ((ftype(1).eq.ftype(5))) then !nc
         k = -2*ftype(1)-ftype(2)+7
      else !cc 
         icc = (-4*(bflav(1)/abs(bflav(1)))+6)/2
         k = 7-ftype(icc)
      endif !nc/cc  

      ! VBFHHMOD: 4->5, 5->6, 6->7
      j1 = 1
      j2 = 5
      j3 = 6!6
      j4 = 7!5

      elseif (bflav(1).eq.0.and.bflav(2).lt.0.and.bflav(5).gt.0) then
        
C*******************  g qbar ---> q qb qb H H **********************

      physToDiag(1)=5
      physToDiag(2)=4
      physToDiag(3)=2
      physToDiag(4)=3
      physToDiag(5)=1
c
      fsign(1) = -1
      fsign(2) =  1
      fsign(3) = -1
      fsign(4) = -1
      gsign    = -1

      polcol = polcolg  

c up- or down-type quark:  
      ! VBFHHMOD: bflav(4)->bflav(5), ...    
      ftype(1) = 2-mod(abs(bflav(7)),2)         
      ftype(2) = 2-mod(abs(bflav(2)),2)      
      ftype(5) = 2-mod(abs(bflav(5)),2)   
      ftype(6) = 2-mod(abs(bflav(6)),2)

      if ((ftype(1).eq.ftype(5))) then !nc
         k = -2*ftype(1)-ftype(2)+7
      else !cc 
         icc = (4*(bflav(2)/abs(bflav(2)))+8)/2
         k = 7-ftype(icc)
      endif !nc/cc  

      ! VBFHHMOD: 4->5, 5->6, 6->7
      j1 = 5!6
      j2 = 7!4
      j3 = 2!5
      j4 = 6!2

      elseif (bflav(1).eq.0.and.bflav(2).lt.0.and.bflav(5).lt.0) then
        
C*******************  g qbar ---> qb qb q H H **********************

      physToDiag(1)=5
      physToDiag(2)=4
      physToDiag(3)=1
      physToDiag(4)=3
      physToDiag(5)=2
c
      fsign(1) = -1
      fsign(2) =  1
      fsign(3) = -1
      fsign(4) = -1
      gsign    = -1

      polcol = polcolg

c up- or down-type quark:      
      ! VBFHHMOD: bflav(4)->bflav(5), ...
      ftype(1) = 2-mod(abs(bflav(7)),2)         
      ftype(2) = 2-mod(abs(bflav(2)),2)      
      ftype(5) = 2-mod(abs(bflav(5)),2)   
      ftype(6) = 2-mod(abs(bflav(6)),2)

      if ((ftype(1).eq.ftype(5))) then !nc
         k = -2*ftype(1)-ftype(2)+7
      else !cc 
         icc = (4*(bflav(2)/abs(bflav(2)))+8)/2
         k = 7-ftype(icc)
      endif !nc/cc  

      ! VBFHHMOD: 4->5, 5->6, 6->7
      j1 = 5
      j2 = 7
      j3 = 2!5
      j4 = 6!2

      ! VBFHHMOD: bflav(4)->bflav(5), ...
      elseif (bflav(1).lt.0.and.bflav(2).eq.0.and.bflav(5).lt.0.and.bflav(6).lt.0) then
 
C*******************  qbar2 g ---> qbar1 qb3 q4 H H  **********************
c
      physToDiag(1)=2             
      physToDiag(2)=5           
      physToDiag(3)=1               
      physToDiag(4)=3!4
      physToDiag(5)=4!3 

      fsign(1) = -1
      fsign(2) = -1
      fsign(3) = -1
      fsign(4) =  1
      gsign    = -1
      
      polcol = polcolg

c up- or down-type quark:  
      ! VBFHHMOD: bflav(4)->bflav(5), ...    
      ftype(1) = 2-mod(abs(bflav(1)),2)         
      ftype(2) = 2-mod(abs(bflav(7)),2)      
      ftype(5) = 2-mod(abs(bflav(5)),2)   
      ftype(6) = 2-mod(abs(bflav(6)),2)

      if ((ftype(1).eq.ftype(5))) then !nc
         k = -2*ftype(1)-ftype(2)+7
      else !cc 
         icc = (-4*(bflav(1)/abs(bflav(1)))+6)/2
         k = 7-ftype(icc)
      endif !nc/cc  

      ! VBFHHMOD: 4->5, 5->6, 6->7
      j1 = 1!4
      j2 = 5!1
      j3 = 6!6
      j4 = 7!5

      elseif (bflav(1).lt.0.and.bflav(2).eq.0.and.bflav(5).lt.0.and.bflav(6).gt.0) then
 
C*******************  qbar2 g ---> qbar1 q4 q3bar H H  **********************
c
      physToDiag(1)=2 
      physToDiag(2)=5            
      physToDiag(3)=1             
      physToDiag(4)=4
      physToDiag(5)=3

      fsign(1) = -1
      fsign(2) = -1
      fsign(3) = -1
      fsign(4) =  1
      gsign    = -1
      
      polcol = polcolg

c up- or down-type quark:     
      ! VBFHHMOD: bflav(4)->bflav(5), ... 
      ftype(1) = 2-mod(abs(bflav(1)),2)         
      ftype(2) = 2-mod(abs(bflav(7)),2)      
      ftype(5) = 2-mod(abs(bflav(5)),2)   
      ftype(6) = 2-mod(abs(bflav(6)),2)

      if ((ftype(1).eq.ftype(5))) then !nc
         k = -2*ftype(1)-ftype(2)+7
      else !cc 
         icc = (-4*(bflav(1)/abs(bflav(1)))+6)/2
         k = 7-ftype(icc)
      endif !nc/cc  

      ! VBFHHMOD: 4->5, 5->6, 6->7
      j1 = 1!4
      j2 = 5!1
      j3 = 6!6
      j4 = 7!5

      else
         
         print*,'this flav combination is not implemented'
         print*,'bflav=',bflav

      endif
         
C*****************  end of process evaluation  **********************

      do mu = 0,3
         do i = 1,5
            pbar(mu,physToDiag(i))=p(mu,i)
         enddo
         qbar(mu) = pbar(mu,5)
      enddo 
      qbar(4) = 0d0


      do i=1,2
         if (phystodiag(i)  .le. 2) then 
         flst_borntags_todiag(phystodiag(i)) = flst_borntags(i,flst_cur_iborn) 
      else
         flst_borntags_todiag(phystodiag(i)+1) = flst_borntags(i,flst_cur_iborn) 
      endif
      enddo
      flst_borntags_todiag(3) = flst_borntags(3,flst_cur_iborn) 
      do i=5,nlegborn
         if (phystodiag(i-2) .le. 2) then 
            flst_borntags_todiag(phystodiag(i-2)) = flst_borntags(i,flst_cur_iborn) 
         else
            flst_borntags_todiag(phystodiag(i-2)+1) = flst_borntags(i,flst_cur_iborn) 
         endif
      enddo

         
c Higgs:
      do mu = 0,3
         pbar(mu,5) = pin(mu,3)
         ! VBFHHMOD: added second higgs
         pbar(mu,6) = pin(mu,4)
c dummy:
         ! VBFHHMOD: 6->7
         pbar(mu,7:nlegborn) = 0d0
      enddo
c
      fsign(5:nlegborn) = 1

      if (k.lt.1.or.k.gt.6) then 
         print*,'in Born you obtained k=',k
         stop 'something went wrong'
      endif   

c k-ordering is: uucc,uuss,ddcc,ddss,udsc,ducs
      call qqhqqj_born_channel(0,pbar,fsign,qbar,gsign,k,res,resc,bbmunu)

      born = res*polcol
c symmetry factor for identical decay particles:
c      born = born*wsymfact

      if (j1*j2*j3*j4.eq.0) then
         print*,'identification of color structures '
         print*,'in Born does not work'
         stop
      endif   


c color correlated Born amplitudes:
      born_col(1:2) = resc(2:3)*polcol  
      bornjk(:,:)=0d0    

      
      ! VBFHHMOD: bflav(4)->bflav(5), ... 
      if(bflav(7).eq.0) then
      do j=1,nlegs   
         if(abs(bflav(j)).le.6) then
            do jj=j+1,nlegs
               if ((j.eq.j1).and.(jj.eq.j2)) then 
                  bornjk(j,jj)=(cf-ca/2d0)*born_col(1)+cf*born_col(2)
                  bornjk(jj,j)=bornjk(j,jj)
               elseif ((j.eq.j3).and.(jj.eq.j4)) then  
                  bornjk(j,jj)=cf*born_col(1)+(cf-ca/2d0)*born_col(2)
                  bornjk(jj,j)=bornjk(j,jj)
               elseif((bflav(j).eq.0).and.(abs(bflav(jj)).le.5)) then
                  if(jj.eq.j1) bornjk(j,jj)=ca/2d0*born_col(1)!cf*born - bornjk(jj,j2)
                  if(jj.eq.j2) bornjk(j,jj)=ca/2d0*born_col(1)!cf*born - bornjk(jj,j1)       
                  if(jj.eq.j3) bornjk(j,jj)=ca/2d0*born_col(2)!cf*born - bornjk(jj,j4)
                  if(jj.eq.j4) bornjk(j,jj)=ca/2d0*born_col(2)!cf*born - bornjk(jj,j3) 
                  bornjk(jj,j)=bornjk(j,jj)
               elseif((bflav(jj).eq.0).and.(abs(bflav(j)).le.5)) then
                  if(j.eq.j1) bornjk(j,jj)=ca/2d0*born_col(1)!cf*born - bornjk(j,j2)
                  if(j.eq.j2) bornjk(j,jj)=ca/2d0*born_col(1)!cf*born - bornjk(j,j1)      
                  if(j.eq.j3) bornjk(j,jj)=ca/2d0*born_col(2)!cf*born - bornjk(j,j4)
                  if(j.eq.j4) bornjk(j,jj)=ca/2d0*born_col(2)!cf*born - bornjk(j,j3) 
                  bornjk(jj,j)=bornjk(j,jj)
               else
                  bornjk(j,jj)=0d0
                  bornjk(jj,j)=bornjk(j,jj)
               endif
               
            enddo !k
         endif !abs
      enddo
      elseif(bflav(1).eq.0) then !Gluon on upper line
      do j=1,nlegs  
         if(abs(bflav(j)).le.6) then
            do jj=j+1,nlegs
               if ((j.eq.j1).and.(jj.eq.j2)) then ! Eventuell Minuszeichen wegen 1 antiquark
                  bornjk(j,jj)=(cf-ca/2d0)*born_col(1)
                  bornjk(jj,j)=bornjk(j,jj)
               elseif ((j.eq.j3).and.(jj.eq.j4)) then  
                  bornjk(j,jj)=cf*born_col(1)
                  bornjk(jj,j)=bornjk(j,jj)
               elseif((bflav(j).eq.0).and.(abs(bflav(jj)).le.5)) then
                  if(jj.eq.j1) bornjk(j,jj)=ca/2d0*born_col(1)!- bornjk(jj,j2)
                  if(jj.eq.j2) bornjk(j,jj)=ca/2d0*born_col(1)! - bornjk(jj,j1)       
                  if(jj.eq.j3) bornjk(j,jj)=0d0!cf*born_col(1) - bornjk(jj,j4)
                  if(jj.eq.j4) bornjk(j,jj)=0d0!cf*born_col(1) - bornjk(jj,j3) 
                  bornjk(jj,j)=bornjk(j,jj)
               elseif((bflav(jj).eq.0).and.(abs(bflav(j)).le.5)) then
                  if(jj.eq.j1) bornjk(j,jj)=ca/2d0*born_col(1)!- bornjk(jj,j2)
                  if(jj.eq.j2) bornjk(j,jj)=ca/2d0*born_col(1)! - bornjk(jj,j1)       
                  if(jj.eq.j3) bornjk(j,jj)=0d0!cf*born_col(1) - bornjk(jj,j4)
                  if(jj.eq.j4) bornjk(j,jj)=0d0!cf*born_col(1) - bornjk(jj,j3) 
                  bornjk(jj,j)=bornjk(j,jj)                  
               else
                  bornjk(j,jj)=0d0
                  bornjk(jj,j)=bornjk(j,jj)
               endif
               
            enddo !k
         endif !abs
      enddo
      elseif(bflav(2).eq.0) then !Gluon on lower line
      do j=1,nlegs  
         if(abs(bflav(j)).le.6) then
            do jj=j+1,nlegs
               if ((j.eq.j1).and.(jj.eq.j2)) then 
                  bornjk(j,jj)=cf*born_col(2)
                  bornjk(jj,j)=bornjk(j,jj)
               elseif ((j.eq.j3).and.(jj.eq.j4)) then  ! gluon on lower line
                  bornjk(j,jj)=(cf-ca/2d0)*born_col(2)
                  bornjk(jj,j)=bornjk(j,jj)
               elseif((bflav(j).eq.0).and.(abs(bflav(jj)).le.5)) then
                  if(jj.eq.j1) bornjk(j,jj)=0d0!cf*born_col(2) - bornjk(jj,j2)
                  if(jj.eq.j2) bornjk(j,jj)=0d0!cf*born_col(2) - bornjk(jj,j1)       
                  if(jj.eq.j3) bornjk(j,jj)=ca/2d0*born_col(2) !- bornjk(jj,j4)
                  if(jj.eq.j4) bornjk(j,jj)=ca/2d0*born_col(2) !- bornjk(jj,j3) 
                  bornjk(jj,j)=bornjk(j,jj)
               elseif((bflav(jj).eq.0).and.(abs(bflav(j)).le.5)) then
                  if(j.eq.j1) bornjk(j,jj)=0d0!cf*born_col(2) - bornjk(j,j2)
                  if(j.eq.j2) bornjk(j,jj)=0d0!cf*born_col(2) - bornjk(j,j1)      
                  if(j.eq.j3) bornjk(j,jj)=ca/2d0*born_col(2)! - bornjk(j,j4)
                  if(j.eq.j4) bornjk(j,jj)=ca/2d0*born_col(2)! - bornjk(j,j3) 
                  bornjk(jj,j)=bornjk(j,jj)
               else
                  bornjk(j,jj)=0d0
                  bornjk(jj,j)=bornjk(j,jj)
               endif
               
            enddo !k
         endif !abs
      enddo
       endif
c spin-correlated Born amplitudes:
c
      bmunu(:,:,:)  = 0d0
c in each relevant case: determine qvec, then call qqHqqj_spc
      do j=1,nlegs
         if(bflav(j).eq.0) then
               do mu=0,3
                  do nu=0,3
                       bmunu(mu,nu,j)=bbmunu(mu,nu)*polcol
                  enddo
               enddo
         endif                 
      enddo 

      return
      end      

c================================================


      subroutine borncolour_lh
c     Sets up the colour for the given flavour configuration
c     already filled in the Les Houches interface.
c     In case there are several colour structure, one
c     should pick one with a probability proportional to
c     the value of the corresponding cross section, for the
c     kinematics defined in the Les Houches interface
      implicit none 
      include 'LesHouches.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_math.h'

c
      integer nlegs
      parameter (nlegs=nlegborn)
      real*8 pin(0:3,nlegs)  
      integer bflav(nlegs)
      real*8 born
      real*8 bbmunu(0:3,0:3)
      real*8 bmunu(0:3,0:3,nlegs),bornjk(nlegs,nlegs)      
      real * 8 borncol(2), ransum
      
      real * 8 random
      integer i
      integer HWW,HZZ


c     fix CKM for VBF WW -> H 
      call particle_identif(HWW,HZZ)
      if (idup(3).eq.HWW) then   
         if ((idup(1).ne.21).and.(idup(2).ne.21)) then
            ! VBFHHMOD: 4->5, 5->6, 6->7
            call CKM_reshuffling(idup(1),idup(5),0)
            call CKM_reshuffling(idup(2),idup(6),0)
         elseif (idup(1).eq.21) then 
            ! VBFHHMOD: 4->5, 5->6, 6->7
            call CKM_reshuffling(idup(7),idup(5),1)
            call CKM_reshuffling(idup(2),idup(6),0)
         elseif (idup(2).eq.21) then 
            ! VBFHHMOD: 4->5, 5->6, 6->7
            call CKM_reshuffling(idup(1),idup(5),0)
            call CKM_reshuffling(idup(7),idup(6),1)
         else
            print*, 'wrong idup assignment:', idup(1:7)
            stop
         endif
      endif
            
c colours of incoming quarks, antiquarks

c     Higgs Bosons             
      icolup(1,3)=0
      icolup(2,3)=0
      !VBFHHMOD: second higgs
      icolup(1,4)=0
      icolup(2,4)=0
c     change Higgs boson codification back to PDG
      idup(3) = 25
      !VBFHHMOD: second higgs
      idup(4) = 25

c     -- colored particles
      !VBFHHMOD: 4->5, 5->6, 6->7
      icolup(1,1)=0
      icolup(2,1)=0
      icolup(1,2)=0
      icolup(2,2)=0
      icolup(1,5)=0
      icolup(2,5)=0
      icolup(1,6)=0
      icolup(2,6)=0
      icolup(1,7)=0
      icolup(2,7)=0
      
c Gluon initiated processes -> unique color combination  

      if(idup(1).eq.21) then
        icolup(1,1)=501
        icolup(2,1)=502
c-----lower quarkline: qq or qbqb?
        if(idup(2).gt.0) then
          icolup(1,2)=503
!     VBFHHMOD: 4->5, 5->6, 6->7
          icolup(1,6)=503        
        elseif(idup(2).lt.0) then
          icolup(2,2)=503
      !VBFHHMOD: 4->5, 5->6, 6->7
          icolup(2,6)=503        
        else
          print*, 'something wrong with gluon initiated process, lower line:', idup(2)
        endif
c-----upper quarkline: qqb or qbq?
      !VBFHHMOD: 4->5, 5->6, 6->7
        if(idup(5).gt.0) then
          icolup(1,5)=501 !outgoing quark 
          icolup(2,7)=502 !-> antiquark from initial -> final state      
        elseif(idup(5).lt.0) then
          icolup(2,5)=502
          icolup(1,7)=501        
        else
          print*, 'something wrong with gluon initiated process, upper line:', idup(5)
        endif     
     
      elseif(idup(2).eq.21) then
        icolup(1,2)=501
        icolup(2,2)=502
c-----upper quarkline: qq or qbqb?
        if(idup(1).gt.0) then
          icolup(1,1)=503
      !VBFHHMOD: 4->5, 5->6, 6->7
          icolup(1,5)=503        
        elseif(idup(1).lt.0) then
          icolup(2,1)=503
      !VBFHHMOD: 4->5, 5->6, 6->7
          icolup(2,5)=503        
        else
          print*, 'something wrong with gluon initiated process, lower line:', idup(2)
        endif
c-----lower quarkline: qqb or qbq?
      !VBFHHMOD: 4->5, 5->6, 6->7
        if(idup(6).gt.0) then
          icolup(1,6)=501 !outgoing quark 
          icolup(2,7)=502 !-> antiquark from initial -> final state      
        elseif(idup(6).lt.0) then
          icolup(2,6)=502
          icolup(1,7)=501        
        else
          print*, 'something wrong with gluon initiated process, upper line:', idup(5)
        endif          
     
      else
         do i=1, nlegborn
           bflav(i)=idup(i)
         enddo      
         call compborn_hhjjj(kn_cmpborn,bflav,born,bmunu,bornjk, borncol)
         
         ransum=random()*born
         
         if(borncol(1).ge.ransum) then !gluon emitted on upper line
            if(idup(1).gt.0) then !quark
              icolup(1,1)=501
              icolup(2,1)=0 
      !VBFHHMOD: 4->5, 5->6, 6->7 
              icolup(1,7)=501
              icolup(2,7)=503
              icolup(1,5)=503
              icolup(2,5)=0
            else          !antiquark
              icolup(2,1)=501
              icolup(1,1)=0  
      !VBFHHMOD: 4->5, 5->6, 6->7
              icolup(2,7)=501
              icolup(1,7)=503
              icolup(2,5)=503
              icolup(1,5)=0
            endif
            
            if(idup(2).gt.0) then !quark
              icolup(1,2)=511
              icolup(2,2)=0
              icolup(1,6)=511
              icolup(2,6)=0
            else
              icolup(2,2)=511
              icolup(1,2)=0
              icolup(2,6)=511
              icolup(1,6)=0            
            endif
         else  !lower line
            if(idup(2).gt.0) then !quark
              icolup(1,2)=501
              icolup(2,2)=0  
              icolup(1,7)=501
              icolup(2,7)=503
              icolup(1,6)=503
              icolup(2,6)=0
            else          !antiquark
              icolup(2,2)=501
              icolup(1,2)=0  
              icolup(2,7)=501
              icolup(1,7)=503
              icolup(2,6)=503
              icolup(1,6)=0
            endif
            
            if(idup(1).gt.0) then !quark
              icolup(1,1)=511
              icolup(2,1)=0
              icolup(1,5)=511
              icolup(2,5)=0
            else
              icolup(2,1)=511
              icolup(1,1)=0
              icolup(2,5)=511
              icolup(1,5)=0            
            endif         
         endif
         
      
      endif
     

      end    
      
c
      subroutine finalize_lh
c     Set up the resonances whose mass must be preserved
c     on the Les Houches interface. 

      implicit none 

      include 'pwhg_physpar.h'   
      logical ini
      data ini/.true./
      save ini
c     lepton masses
      real *8 lepmass(3),decmass
      common/clepmass/lepmass,decmass

      end

      function signn(j)
      implicit none
      integer signn, j
      if (j.eq.0) then
         signn =0
      else
         signn =j/abs(j)
      endif
      end


      subroutine CKM_reshuffling(flin,flout,gtype)
      implicit none
      integer flin,flout
      integer i,j
      logical ini
      real * 8 CKM_sq(6,6)
      include 'PhysPars.h' 
      integer signn
      external signn
      save CKM_sq,ini
      data ini/.true./
      integer nf_max
      integer gtype 
c     no initial or final state t quark!!
      nf_max = 4
      if (ini) then
         do j=1,nf_max
            do i=1,nf_max
               CKM_sq(i,j) = ph_CKM_matrix(i,j)**2
            enddo
         enddo
         ini = .false.
      endif      
      call pick_random(nf_max,CKM_sq(1,abs(flin)),flout)

      if (gtype.ne.1) flout = flout * signn(flin)
      end
