      subroutine setvirtual(p,vflav,virtual)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_math.h'
      include 'pwhg_st.h'
      integer nleg
      parameter (nleg=nlegborn)
      real * 8 p(0:3,nleg)
      integer vflav(nleg)
      real * 8 virtual, borncol(2)
      real * 8 bmunu_dum(0:3,0:3,nleg),bornjk_dum(nleg,nleg),born

      real *8 powheginput
      external powheginput 
      logical, save :: firsttime = .true. 

      integer fakevirt
      save fakevirt 

c================================================

c numbering of momenta is q(1) q(2) -> H(3) q(4)q(5)"g(6)"
c
      virtual = 0d0
      if (firsttime) then
         fakevirt=powheginput("#fakevirt")
         if (fakevirt == 1) write(*,*) 'WARNING: Using fakevirt !'
         firsttime = .false.
      endif

      call calc_als

      if(fakevirt.eq.1) then  
         call compborn_hjjj(p,vflav,born,bmunu_dum,bornjk_dum, borncol) 
         virtual = 0.2d0*born

      else    
         call compvirt_hjjj(p,vflav,virtual) 

c        cancel as/(2pi) associated with amp2. 
         virtual = virtual/(st_alpha/(2d0*pi))

      endif
         
      return
      end

c
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine compvirt_hjjj(pin,bflav,virtual)
      implicit none
c
      include 'nlegborn.h'
      include 'pwhg_math.h'
      include 'pwhg_flst.h'
c
      integer nlegs
      parameter (nlegs=nlegborn)
      real*8 pin(0:3,nlegs)  
      integer bflav(nlegs)
      real*8 virtual
c
c vbfnlo stuff:
      include 'global_col.inc'
      integer nlo
      real*8 p(0:3,np), v(0:3,nv)
      real*8 pbar(0:3,4+nv),qbar(0:4)
      real*8 polcol,polcolq,polcolg
      real*8 res,resc_dum(3)
      real*8 tri,box

      real*8 N ! color factors
      parameter(N=3d0)

      complex*16 zero
      parameter (zero=(0d0,0d0))
c
c declare local variables
c
      integer i,mu
      integer FSIGN(4+nv),gsign,physToDiag(5)
      
      real*8 nans(2,2,3),cans(2,3)
      logical nc_type
      integer k,icc

      integer ftype(1:6)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      polcol = 0d0
      polcolq = 1d0/(4d0*N**2)
      polcolg = 1d0/(4d0*N*(N**2-1))

      virtual = 0d0

      ftype(1:6) = 0
      do mu = 0,3
c fermions:
         p(mu,1) = pin(mu,1)
         p(mu,2) = pin(mu,2)
         
         if (bflav(1)*bflav(2).ne.0) then ! final-state gluon 
            if (bflav(6).eq.0) then
               p(mu,3) = pin(mu,4)
               p(mu,4) = pin(mu,5) 
               p(mu,5) = pin(mu,6)  !gluon
            elseif (bflav(4).eq.0) then
               p(mu,3) = pin(mu,5)
               p(mu,4) = pin(mu,6) 
               p(mu,5) = pin(mu,4)  
            elseif (bflav(5).eq.0) then
               p(mu,3) = pin(mu,4)
               p(mu,4) = pin(mu,6) 
               p(mu,5) = pin(mu,5)  
            endif
         else   ! initial-state gluon 
            p(mu,3) = pin(mu,4)
            p(mu,4) = pin(mu,5) 
            p(mu,5) = pin(mu,6) 
         endif ! fin/in state gluon
      enddo ! mu

      if (bflav(1).gt.0.and.bflav(2).gt.0) then

C*******************  q1 q3 ---> q2 q4 g H   **********************

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
      ftype(1) = 2-mod(abs(bflav(1)),2)         
      ftype(2) = 2-mod(abs(bflav(2)),2)      
      ftype(4) = 2-mod(abs(bflav(4)),2)   
      ftype(5) = 2-mod(abs(bflav(5)),2)

      if ((ftype(1).eq.ftype(4))) then
         k = -2*ftype(1)-ftype(2)+7
      else !cc 
         icc = (-3*(bflav(1)/abs(bflav(1)))+5)/2
         k = 7-ftype(icc)
      endif   

      elseif (bflav(1).gt.0.and.bflav(2).lt.0) then
c            
C******************* q1 qb4 ---> q2 qb3 g H   **********************
      
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
      ftype(1) = 2-mod(abs(bflav(1)),2)         
      ftype(2) = 2-mod(abs(bflav(2)),2)      
      ftype(4) = 2-mod(abs(bflav(4)),2)   
      ftype(5) = 2-mod(abs(bflav(5)),2)

      if ((ftype(1).eq.ftype(4))) then !nc
         k = -2*ftype(1)-ftype(2)+7
      else !cc 
         icc = (-3*(bflav(1)/abs(bflav(1)))+5)/2
         k = 7-ftype(icc)
      endif !nc/cc  

      elseif (bflav(1).lt.0.and.bflav(2).gt.0) then
      
C******************* qbar2 q3 ---> qbar1 q4 g H   **********************
      
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
      ftype(1) = 2-mod(abs(bflav(1)),2)         
      ftype(2) = 2-mod(abs(bflav(2)),2)      
      ftype(4) = 2-mod(abs(bflav(4)),2)   
      ftype(5) = 2-mod(abs(bflav(5)),2)

      if ((ftype(1).eq.ftype(4))) then !nc
         k = -2*ftype(1)-ftype(2)+7
      else !cc 
         icc = (-3*(bflav(1)/abs(bflav(1)))+5)/2
         k = 7-ftype(icc)
      endif !nc/cc  

      elseif (bflav(1).lt.0.and.bflav(2).lt.0) then

C*******************  qbar2 qb4 ---> qbar1 qb3 g H   **********************

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
      ftype(1) = 2-mod(abs(bflav(1)),2)         
      ftype(2) = 2-mod(abs(bflav(2)),2)      
      ftype(4) = 2-mod(abs(bflav(4)),2)   
      ftype(5) = 2-mod(abs(bflav(5)),2)

      if ((ftype(1).eq.ftype(4))) then !nc
         k = -2*ftype(1)-ftype(2)+7
      else !cc 
         icc = (-3*(bflav(1)/abs(bflav(1)))+5)/2
         k = 7-ftype(icc)
      endif !nc/cc  

      elseif (bflav(1).eq.0.and.bflav(2).gt.0.and.bflav(4).gt.0) then
         
c*******************  g q ---> q q qb H   **********************

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
      ftype(1) = 2-mod(abs(bflav(6)),2)         
      ftype(2) = 2-mod(abs(bflav(2)),2)      
      ftype(4) = 2-mod(abs(bflav(4)),2)   
      ftype(5) = 2-mod(abs(bflav(5)),2)

      if ((ftype(1).eq.ftype(4))) then !nc
         k = -2*ftype(1)-ftype(2)+7
      else !cc
         icc = (3*(bflav(2)/abs(bflav(2)))+7)/2
         k = 7-ftype(icc)
      endif !nc/cc  

      elseif (bflav(1).eq.0.and.bflav(2).gt.0.and.bflav(4).lt.0) then
         
c*******************  g q ---> qb q q H   **********************
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
      ftype(1) = 2-mod(abs(bflav(6)),2)         
      ftype(2) = 2-mod(abs(bflav(2)),2)      
      ftype(4) = 2-mod(abs(bflav(4)),2)   
      ftype(5) = 2-mod(abs(bflav(5)),2)     

      if ((ftype(1).eq.ftype(4))) then !nc
         k = -2*ftype(1)-ftype(2)+7
      else !cc 
         icc = (3*(bflav(2)/abs(bflav(2)))+7)/2
         k = 7-ftype(icc)
      endif !nc/cc  

      elseif (bflav(1).gt.0.and.bflav(2).eq.0.
     &    and.bflav(4).gt.0.and.bflav(5).gt.0) then
      
C*******************  q g ---> q q qb H   **********************
      
      physToDiag(1)=1             
      physToDiag(5)=3            
      physToDiag(3)=2             
      physToDiag(2)=5
      physToDiag(4)=4
c
      fsign(1) =  1
      fsign(2) =  1
      fsign(3) = -1
      fsign(4) =  1
      gsign    = -1

      polcol = polcolg
c
c up- or down-type quark:      
      ftype(1) = 2-mod(abs(bflav(1)),2)         
      ftype(2) = 2-mod(abs(bflav(6)),2)      
      ftype(4) = 2-mod(abs(bflav(4)),2)   
      ftype(5) = 2-mod(abs(bflav(5)),2)

      if ((ftype(1).eq.ftype(4))) then !nc
         k = -2*ftype(1)-ftype(2)+7
      else !cc 
         icc = (-3*(bflav(1)/abs(bflav(1)))+5)/2
         k = 7-ftype(icc)
      endif !nc/cc  

      elseif (bflav(1).gt.0.and.bflav(2).eq.0.and.
     &        bflav(4).gt.0.and.bflav(5).lt.0) then
      
C*******************  q g ---> q qb q H   **********************

      physToDiag(1)=1             
      physToDiag(4)=3            
      physToDiag(3)=2             
      physToDiag(2)=5
      physToDiag(5)=4
c
      fsign(1) =  1
      fsign(2) =  1
      fsign(3) = -1
      fsign(4) =  1
      gsign    = -1

      polcol = polcolg

c up- or down-type quark:      
      ftype(1) = 2-mod(abs(bflav(1)),2)         
      ftype(2) = 2-mod(abs(bflav(6)),2)      
      ftype(4) = 2-mod(abs(bflav(4)),2)   
      ftype(5) = 2-mod(abs(bflav(5)),2)

      if ((ftype(1).eq.ftype(4))) then !nc
         k = -2*ftype(1)-ftype(2)+7
      else !cc 
         icc = (-3*(bflav(1)/abs(bflav(1)))+5)/2
         k = 7-ftype(icc)
      endif !nc/cc  

      elseif (bflav(1).eq.0.and.bflav(2).lt.0.and.bflav(4).gt.0) then
        
C*******************  g qbar ---> q qb qb H  **********************

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
      ftype(1) = 2-mod(abs(bflav(6)),2)         
      ftype(2) = 2-mod(abs(bflav(2)),2)      
      ftype(4) = 2-mod(abs(bflav(4)),2)   
      ftype(5) = 2-mod(abs(bflav(5)),2)

      if ((ftype(1).eq.ftype(4))) then !nc
         k = -2*ftype(1)-ftype(2)+7
      else !cc 
         icc = (3*(bflav(2)/abs(bflav(2)))+7)/2
         k = 7-ftype(icc)
      endif !nc/cc  

      elseif (bflav(1).eq.0.and.bflav(2).lt.0.and.bflav(4).lt.0) then
        
C*******************  g qbar ---> qb qb q H  **********************

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
      ftype(1) = 2-mod(abs(bflav(6)),2)         
      ftype(2) = 2-mod(abs(bflav(2)),2)      
      ftype(4) = 2-mod(abs(bflav(4)),2)   
      ftype(5) = 2-mod(abs(bflav(5)),2)

      if ((ftype(1).eq.ftype(4))) then !nc
         k = -2*ftype(1)-ftype(2)+7
      else !cc 
         icc = (3*(bflav(2)/abs(bflav(2)))+7)/2
         k = 7-ftype(icc)
      endif !nc/cc  

      elseif (bflav(1).lt.0.and.bflav(2).eq.0.and.bflav(5).lt.0) then
 
C*******************  qbar2 g ---> qbar1 qb3 q4 H   **********************
c
      physToDiag(1)=2             
      physToDiag(5)=4!3            
      physToDiag(3)=1             
      physToDiag(2)=5
      physToDiag(4)=3!4

      fsign(1) = -1
      fsign(2) = -1
      fsign(3) = -1
      fsign(4) =  1
      gsign    = -1
      
      polcol = polcolg

c up- or down-type quark:      
      ftype(1) = 2-mod(abs(bflav(1)),2)         
      ftype(2) = 2-mod(abs(bflav(6)),2)      
      ftype(4) = 2-mod(abs(bflav(4)),2)   
      ftype(5) = 2-mod(abs(bflav(5)),2)

      if ((ftype(1).eq.ftype(4))) then !nc
         k = -2*ftype(1)-ftype(2)+7
      else !cc 
         icc = (-3*(bflav(1)/abs(bflav(1)))+5)/2
         k = 7-ftype(icc)
      endif !nc/cc  

      elseif (bflav(1).lt.0.and.bflav(2).eq.0.and.bflav(5).gt.0) then
 
C*******************  qbar2 g ---> qbar1 q4 q3bar H   **********************
c
      physToDiag(1)=2             
      physToDiag(5)=3            
      physToDiag(3)=1             
      physToDiag(2)=5
      physToDiag(4)=4

      fsign(1) = -1
      fsign(2) = -1
      fsign(3) = -1
      fsign(4) =  1
      gsign    = -1
      
      polcol = polcolg

c up- or down-type quark:      
      ftype(1) = 2-mod(abs(bflav(1)),2)         
      ftype(2) = 2-mod(abs(bflav(6)),2)      
      ftype(4) = 2-mod(abs(bflav(4)),2)   
      ftype(5) = 2-mod(abs(bflav(5)),2)

      if ((ftype(1).eq.ftype(4))) then !nc
         k = -2*ftype(1)-ftype(2)+7
      else !cc 
         icc = (-3*(bflav(1)/abs(bflav(1)))+5)/2
         k = 7-ftype(icc)
      endif !nc/cc  

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

c Higgs:
      do mu = 0,3
         pbar(mu,5) = pin(mu,3)
c dummy:
         pbar(mu,6:4+nv) = 0d0
      enddo
c
      fsign(5:4+nv) = 1

c     k-ordering is: uucc,uuss,ddcc,ddss,udsc,ducs
      if (k.lt.1.or.k.gt.6) then 
         print*,'in virtuals you obtained k=',k
         stop 'something went wrong'
      endif   
	 

      call qqhqqj_vonly_channel(1,pbar,fsign,qbar,gsign,k,res,resc_dum) 
      

      virtual = res*polcol
c symmetry factor for identical decay particles:
c      virtual = virtual*wsymfact

      return
      end      
 
