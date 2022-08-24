      function theta(x)
      implicit none
      double precision theta,x
      if (x.gt.0d0) then
         theta = 1d0
      else
         theta = 0d0
      endif
      end


c     imaginary part of a/b
c     IF a/b > 0 then returns 0
c     IF a/b < 0 then
c        if a < 0 then return +1
c        else return -1
      function im_part(a,b)
      implicit none
      double precision a,b,im_part
      if ((a/b).gt.0d0) then
         im_part = 0d0
      else
         if (a.lt.0d0) then
            im_part = +1d0
         else
            im_part = -1d0
         endif
      endif

      end


c    ***********   D01m_fin(s,t,m4sq,musq)   *************
c
c Scalar box with MASSLESS PROPAGATORS and with external kinematics:
c
c                     k
c     p1 ->-----------<-------------<-- p4
c             |                 |
c             |                 |
c             |                 |
c     p2 ->-------------------------<-- p3
c
c
c    s = (p1+p2)^2
c    t = (p2+p3)^2
c    p1^2 = 0, p2^2 = 0, p3^2 = 0, p4^2 = m4sq <>0
c    musq = mu^2 = reference dimensional scale
c
c  int d^dk/(2 pi)^d 1/(k^2)/(k+p1)^2/(k+p1+p2)^2/(k+p1+p2+p3)^2 =
c          N_ep * D01m_fin(s,t,m3sq,m4sq,musq);
c          N_ep = i/(4 pi)^2 (4 pi)^ep Gamma(1+ep) (musq)^(-ep)

      function D01m_fin(s,t,m4sq,musq)
      implicit none
      double complex D01m_fin
      double precision s,t,m4sq,musq
      double precision ms,mt,mm4sq
      double complex lnms,lnmt,lnmm4sq,lnsot
      double complex li2arg1,li2arg2
      double precision arg1,arg2
      double complex ipi
      parameter (ipi=(0d0,3.14159265358979323846264338328d0))
      double precision pi,pi2, pi2o3t2
      parameter (pi=3.14159265358979323846264338328d0,
     -     pi2 = 9.86960440108935861883449099988d0,
     -   pi2o3t2= 6.57973626739290574588966066658410076d0 )

      double complex ris
      double precision prefactor, theta, dilog, im_part
      external theta,im_part,dilog


      if (musq.lt.0d0) then
         write(*,*)
     -        'POSSIBLE ERROR IN D01m_fin: SCALE MUSQ LESS THAN ZERO!!'
      endif

      prefactor = 1d0/(s*t)
      ms = -s/musq
      mt = -t/musq
      mm4sq = - m4sq/musq

      lnms = log(abs(ms)) - ipi*theta(-ms)
      lnmt = log(abs(mt)) - ipi*theta(-mt)
      lnmm4sq = log(abs(mm4sq)) - ipi*theta(-mm4sq)
      lnsot = log(abs(s/t)) + ipi * im_part(s,t)

      arg1 = 1-m4sq/s
      arg2 = 1-m4sq/t

c     (m4sq/s.lt.0d0)
      if (arg1.gt.1d0) then
         li2arg1 = -dilog(1d0/arg1) - log(arg1)**2/2+pi2/3
     -        - ipi*log(arg1)*im_part(m4sq,s)
      else
         li2arg1 = dilog(arg1)
      endif

c     (m4sq/t.lt.0d0)
      if (arg2.gt.1d0)  then
         li2arg2 = -dilog(1d0/arg2) - log(arg2)**2/2+pi2/3
     -        - ipi*log(arg2)*im_part(m4sq,t)
      else
         li2arg2 = dilog(arg2)
      endif

      ris = -lnmm4sq**2+lnms**2+lnmt**2-lnsot**2-2d0*(li2arg1+li2arg2)
     -   -pi2o3t2

      D01m_fin = prefactor * ris
      end

c
c---------------- BCD_fill_v(k1,k2,q1,q2) -------------------------------
c
      subroutine BCD_fill_v(k1,k2,q1,q2)
      implicit none
      double precision k1(0:3),k2(0:3),q1(0:4),q2(0:4)
c
c     q1s = 0 for  the massless gluon
c
c This subroutine does the same as BCD_fill for a single momentum
c assignment as in the diagram
c
c     ID = 1  D0(k2,q2,q1),Dij(k2,q2,q1)
c     k1  -->---------->------------- k2
c                 $         $
c                 $         $
c                 $ q1      $  q2
c                 $         $ 
c               gluon     vector       
c 
c     ID = 2  D0(k2,q1,q2),Dij(k2,q1,q2)
c     k1  -->---------->------------- k2
c                 $         $
c                 $         $
c                 $ q2      $  q1
c                 $         $ 
c               vector     gluon      
c
c     ID = 3 D0(q1,k2,q2),Dij(q1,k2,q2)
c
c
c results are stored in Dijv(*,*,ID) etc. i.e.
c
c     renormalization scale
c
c     #include "VBFNLO/utilities/scales.inc"
c
Common block for the output
      complex*16 D0v(3), Dijv(3,13,3), Teps1(2), Teps2(2), Tborn(2,2),
     $     Tg1(2),Tg2(2),Ceps(2)
      common /bcd_qqvg/ D0v,Dijv,Teps1,Teps2,Tborn,Tg1,Tg2,Ceps
      double precision s, t, q1sq, q2sq,u
c
c local variables, external functions etc.
      double precision p1p2, p1p3, p2p3
      double precision musq     !renormalization scale
      complex*16 C0123(3),  C0124(3),  C0134(3),  C0234(3)
      complex*16 Cij123(2,4), Cij124(2,4), 
     1           Cij134(2,4), Cij234(2,4)
      complex*16 B0a(3), B0aq1, B0aq2, B0at, B0as,B0au
      complex*16 B0t, C0t, D01m_fin
      external B0t, C0t, D01m_fin
      double precision tmq2i,umq2i
      double precision dot0p, psumsq
      external dot0p, psumsq
      double precision pi,pi2o3
      parameter (pi=3.14159 26535 89793d0,pi2o3=pi**2/3d0)

      include 'scales.inc'

c
      s = -2*dot0p(k1,k2)
      t = psumsq(k2,q2)
      u = psumsq(k2,q1)
      q1sq = dot0p(q1,q1)        !q1s = 0 
c      print*,'q1sq=',q1sq
      q2sq = dot0p(q2,q2)
c      print*,'q2sq=',q2sq
      musq = mursq(1,1)        !set in scales.inc
c
c now determine the required B0, C0 and D0 functions
c D0v(1) = D0t(s,t,0,q2sq,musq) = D0t(s,t,q2sq,0,musq)
c D0v(2) = D0t(s,u,0,q2sq,musq) = D0t(s,u,q2sq,0,musq)
      D0v(1) = D01m_fin(s,t,q2sq,musq)
      D0v(2) = D01m_fin(s,u,q2sq,musq)
      D0v(3) = D01m_fin(t,u,q2sq,musq)
c ID=1
c C0123 = C0(0,q2sq,t), C0124 = C0(0,s,0), C0134 = C0(t,0,0),
c C0234 = C0(0,q2sq,s)
c ID=2
c C0123 = C0(0,0,u), C0124 =C0(0,s,0) , C0134 = C0(u,q2s,0),
c C0234 = C0(0,q2s,s)
c ID=3 
c  C0123 = C0(0,0,u), C0124 =C0(0,t,0) , C0134 = C0(u,q2s,0),
c C0234 = C0(0,q2s,t)
      
      C0124(1) = C0t(0d0,s,0d0,musq)
      C0234(1) = C0t(q2sq,0d0,s,musq)
      C0123(1) = C0t(0d0,q2sq,t,musq)
      C0134(1) = C0t(t,0d0,0d0,musq) !
c
      C0124(2) = C0t(0d0,s,0d0,musq) !=C0124(1)
      C0234(2) = C0t(0d0,q2sq,s,musq) !=C0234(1)
      C0123(2) = C0t(0d0,0d0,u,musq)
      C0134(2) = C0t(u,q2sq,0d0,musq)
c     
      C0124(3) = C0t(0d0,t,0d0,musq) !=C0134(1) 
      C0234(3) = C0t(0d0,q2sq,t,musq) !=C0123(1) 
      C0123(3) = C0t(0d0,0d0,u,musq) !=C0123(2)
      C0134(3) = C0t(u,q2sq,0d0,musq) !=C0134(2)    
c
c     and now the B0 functions: B0t(q2sq), B0t(t)
c     Note: B0t(s,musq=s) = 2
      B0aq1 = 0.0d0
      B0aq2 = B0t(q2sq,musq)
      B0at  = B0t(t,musq)
      B0as  = B0t(s,musq)
      B0au  = B0t(u,musq)
c
c Now everything is ready to calculate the C_ij and D_ij functions
c     ID = 1
c     p1 = k2,p2 = q2,p3 = q1
      B0a(1) = 0
      B0a(2) = B0aq2
      B0a(3) = B0at
      call tens3(0d0,0d0,q2sq,t,B0a,C0123(1),Cij123(1,1))
      B0a(1) = 0d0
      B0a(2) = B0as
      B0a(3) = 0d0
      call tens3(0d0,0d0,s,0d0,B0a,C0124(1),Cij124(1,1))
      B0a(1) = B0at
      B0a(2) = 0d0
      B0a(3) = 0d0
      call tens3(0d0,t,0d0,0d0,B0a,C0134(1),Cij134(1,1))
      B0a(1) = b0aq2
      B0a(2) = 0d0
      B0a(3) = B0as
      call tens3(0d0,q2sq,0d0,s,B0a,C0234(1),Cij234(1,1))
      p1p2 = dot0p(k2,q2)
      p1p3 = dot0p(k2,q1)
      p2p3 = dot0p(q2,q1)
      call tens4( 0d0, 0d0, q2sq, 0d0, p1p2, p1p3, p2p3, 
     &            C0123(1),  C0124(1),  C0134(1),  C0234(1),
     &            Cij123(1,1), Cij124(1,1), 
     &            Cij134(1,1), Cij234(1,1), 
     &            D0v(1), Dijv(1,1,1) )
c
c     ID = 2
c     p1 = k2,p2=q1,p3=q2 
c
      B0a(1) = 0
      B0a(2) = 0
      B0a(3) = B0au
      call tens3(0d0,0d0,0d0,u,B0a,C0123(2),Cij123(1,1))
      B0a(1) = 0d0
      B0a(2) = B0as
      B0a(3) = 0d0
      call tens3(0d0,0d0,s,0d0,B0a,C0124(2),Cij124(1,1))
      B0a(1) = B0au
      B0a(2) = B0aq2
      B0a(3) = 0d0
      call tens3(0d0,u,q2sq,0d0,B0a,C0134(2),Cij134(1,1))
      B0a(1) = 0d0
      B0a(2) = B0aq2
      B0a(3) = B0as
      call tens3(0d0,0d0,q2sq,s,B0a,C0234(2),Cij234(1,1))
      p1p2 = dot0p(k2,q1)
      p1p3 = dot0p(k2,q2)
      p2p3 = dot0p(q2,q1)
      call tens4( 0d0, 0d0, 0d0, q2sq, p1p2, p1p3, p2p3, 
     &            C0123(2),  C0124(2),  C0134(2),  C0234(2),
     &            Cij123(1,1), Cij124(1,1), 
     &            Cij134(1,1), Cij234(1,1), 
     &            D0v(2), Dijv(1,1,2) )
c
c      ID = 3
c     p1 = q1,p2=k2,p3=q2
c
      B0a(1) = 0
      B0a(2) = 0
      B0a(3) = B0au
      call tens3(0d0,0d0,0d0,u,B0a,C0123(3),Cij123(1,1))
      B0a(1) = 0d0
      B0a(2) = B0at
      B0a(3) = 0d0
      call tens3(0d0,0d0,t,0d0,B0a,C0124(3),Cij124(1,1))
      B0a(1) = B0au
      B0a(2) = B0aq2
      B0a(3) = 0d0
      call tens3(0d0,u,q2sq,0d0,B0a,C0134(3),Cij134(1,1))
      B0a(1) = 0d0
      B0a(2) = B0aq2
      B0a(3) = B0at
      call tens3(0d0,0d0,q2sq,t,B0a,C0234(3),Cij234(1,1))
      p1p2 = dot0p(k2,q1)
      p1p3 = dot0p(q2,q1)
      p2p3 = dot0p(k2,q2)
      call tens4( 0d0, 0d0, 0d0, q2sq, p1p2, p1p3, p2p3, 
     &            C0123(3),  C0124(3),  C0134(3),  C0234(3),
     &            Cij123(1,1), Cij124(1,1), 
     &            Cij134(1,1), Cij234(1,1), 
     &            D0v(3), Dijv(1,1,3) )
c     
c     
c and the coefficients of the tensors in M_6:
c     
c     check these formulas
      tmq2i = 1d0/(t-q2sq)
      umq2i = 1d0/(u-q2sq)
      
      Teps1(1) = (2*B0at+1)/t
      Teps2(1) = ( (B0at-b0aq2)*(2*t+3*q2sq)*tmq2i
     &     +2*b0aq2+1-2*q2sq*C0123(1) ) * tmq2i
      Tg1(1) = B0at/t
      Tg2(1) = (B0at-B0aq2)*tmq2i

      Ceps(1) = 2*(t*C0134(1)+1d0)/t
      Ceps(2) = 2*(u*C0123(2)+1d0)/u
c
      Tborn(1,1) =  ( 2*q2sq*(B0at-b0aq2) + t*B0at
     &     -q2sq*b0aq2 )*tmq2i - 2*q2sq*C0123(1)
     &     +2*B0at -5d0 +pi2o3 !-5d0 + pi2o3  
c
c      print*,'C0123-c0123=',C0123(1)-C0123(1)
      Teps1(2) = (2*B0au+1)/u
      Teps2(2) = ( (B0au-b0aq2)*(2*u+3*q2sq)*umq2i
     &     +2*b0aq2+1-2*q2sq*C0134(2) ) * umq2i
      Tg1(2) = B0au/u
      Tg2(2) = (B0au-B0aq2)*umq2i
      Tborn(2,1) =  ( 2*q2sq*(B0au-b0aq2) + u*B0au
     &     -q2sq*b0aq2 )*umq2i - 2*q2sq*C0134(2)
     &     +2*B0au -5d0 + pi2o3 !-5d0 + pi2o3 
c 
c     5-pi2o3 + cvirt = -3 see expression for div piece  
c     
      Tborn(1,2) = -2*(t*C0134(1)+1) +( 2*q2sq*(B0at-b0aq2) + t*B0at
     &     -q2sq*b0aq2 )*tmq2i - 2*q2sq*C0123(1)
     &     + B0at -5d0 + pi2o3 !-5d0 + pi2o3 
c
      Tborn(2,2) = -2*(u*C0123(2)+1) +( 2*q2sq*(B0au-b0aq2) + u*B0au
     &     -q2sq*b0aq2 )*umq2i - 2*q2sq*C0134(2)
     &     + B0au -5d0 + pi2o3 !-5d0 + pi2o3
c      print*,'bcd_fill'
c      print*,'Teps1=',Teps1
c      print*,'Teps2=',Teps2
c      print*,'Tg1=',Tg1
c      print*,'Tg2=',Tg2
c      print*,'Tborn(1,1)=',Tborn(1,1)
c      print*,'Tborn(2,1)=',Tborn(2,1)
c      print*,'Tborn(1,2)=',Tborn(1,2)
c      print*,'Tborn(2,2)=',Tborn(2,2)
      return
      end

      
c$$$c-------------  spence function  vli2(z) ------------------------
c
         complex*16 function vli2(zin)
         implicit none
         complex*16 zin, z, u, u2, unpo, ans, zext
         double precision r, r2, r2n, fac
c
c determine the value of the dilogarithm 
c
c    vli2(z) = - int_0^1  log(1-zt)/t dt  with cut along the positive 
c                                        real axis, z>1
c
c	Dieter Zeppenfeld, <dieter@pheno.physics.wisc.edu>
c	Initial version:  2000 November 6
c	Last modified:    2000 November 12
c  
         integer i
         double precision c0,c1,c2,c4,c6,c8,c10,c12,c14,c16,c18,c20,c22
         double precision b0,b1,b2,b4,b6,b8,b10,b12,b14,b16,b18,b20,b22
         double precision d0,d1,d2,d4,d6,d8,d10,d12,d14,d16,d18,d20,d22
         parameter (b0=1d0,            d0 =1d0,      c0= b0/d0)
         parameter (b1=-1d0/2d0,       d1 =d0*2d0,   c1= b1/d1)
         parameter (b2= 1d0/6d0,       d2 =d1*3d0,   c2= b2/d2)
         parameter (b4=-1d0/30d0,      d4 =d2*20d0,  c4= b4/d4)
         parameter (b6=1d0/42d0,       d6 =d4*42d0,  c6= b6/d6)
         parameter (b8=-1d0/30d0,      d8 =d6*72d0,  c8= b8/d8)
         parameter (b10=5d0/66d0,      d10=d8*110d0, c10=b10/d10)
         parameter (b12=-691d0/2730d0, d12=d10*156d0,c12=b12/d12)
         parameter (b14=7d0/6d0,       d14=d12*210d0,c14=b14/d14)
         parameter (b16=-3617d0/510d0, d16=d14*272d0,c16=b16/d16)
         parameter (b18=43867d0/798d0, d18=d16*342d0,c18=b18/d18)
         parameter (b20=-174611d0/330d0,d20=d18*420d0,c20=b20/d20)
         parameter (b22=854513d0/138d0,d22=d20*506d0,c22=b22/d22)
         double precision eps, epst, pi, pi2o6
         parameter (eps=1d-16, epst=1d-3)
         parameter (pi=3.14159 26535 89793d0, pi2o6=pi**2/6)
c
c debug information
         logical ldebug
         parameter (ldebug=.false.)
   
         z = zin
c      print*,' vli2 call with z = ',z
         u = z**2
         r2 = dreal(z)**2+dimag(z)**2 
         if (r2.lt.eps) then
            vli2 = z + u/4d0
            return
         elseif (r2.lt.epst) then
            ans = z + u/4
            do i = 3,11
               u = u*z
               ans = ans + u/i**2
            enddo
            vli2 = ans
            return
         endif
         if (dreal(z).ge.1d0 .and. dimag(z).eq.0 ) then
            z = z + (0d0,1d0)*eps
         endif
c
c use z-->1/z and z--> 1-z mappings of the spence function to restrict 
c agument to unit circle in the complex plane with Re(z) <= 0.5
c
         zext = (0d0,0d0)
         fac = 1
         if (r2.gt.1d0) then     ! map z ---> 1/z
            fac = -fac
            zext = -pi2o6 - 0.5d0*(log(-z))**2
            z = 1/z
         endif
         if (dreal(z).gt.0.5d0) then     ! map new z ---> 1-z
            zext = zext + fac*(pi2o6-log(z)*log(1-z))
            fac = -fac
            z = 1-z
         endif
c
c now use t = 1 - exp(-u) mapping to write Li(z) in terms of Bernoulli 
c numbers
c
         u = - log(1-z)
         r2 = abs(u)**2
         u2 = u*u
         ans = u*(c0 + u*(c1+c2*u))
         r2n = r2*r2       !r^4
   
         unpo = u2*u2*u
         ans = ans + c4*unpo
   
         unpo = unpo*u2
         ans = ans + c6*unpo
   
         r = r2n*r2n       !r^8
         unpo = unpo*u2
         ans = ans + c8*unpo
   
         r2n = r*r2        !r^10 
         if ((r2n*c10).gt.eps) then
            unpo = unpo*u2
            ans = ans + c10*unpo
         else
            vli2 = fac * ans + zext
            if (ldebug) print*,' exit vli2s at n=8 '
            return
         endif
   
         unpo = unpo*u2
         ans = ans + c12*unpo
   
         unpo = unpo*u2
         ans = ans + c14*unpo
   
         unpo = unpo*u2
         ans = ans + c16*unpo
   
         r2n = r2n*r
         if ((r2n*c18).gt.eps) then
            unpo = unpo*u2
            ans = ans + c18*unpo
         else
            vli2 = fac * ans + zext
            if (ldebug) print*,' exit vli2s at n=16 '
            return
         endif
   
         unpo = unpo*u2
         ans = ans + c20*unpo
   
         unpo = unpo*u2
         ans = ans + c22*unpo
   
         vli2 = fac * ans + zext
         end
