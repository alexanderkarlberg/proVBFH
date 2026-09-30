c-----------------------------------------------------------------------
c proVBFH-cs stage 3: the VBF H + 2 parton Born |M|^2 (averaged over
c initial spins and colours) in VBFNLO's conventions, for the line-level
c Catani-Seymour dipoles of the (1,1) contribution. VBFNLO has no
c H + 2 parton routine in proVBFH; the amplitude is the one of its
c H + 3 parton Born (qqhqqj-born.f) without the gluon:
c   M = clr(t1,V1,h1) clr(t2,V2,h2) b(6,V1,V2) xmw D_V1(q1) D_V2(q2)
c       [ubar gamma^mu u]_1 [ubar gamma_mu u]_2 ,
c colour factor 9/9. The squared currents are 16 (in1.in2)(out1.out2)
c for equal chiralities and 16 (in1.out2)(out1.in2) otherwise, with
c (in, out) the incoming and outgoing momenta of a quark line and
c exchanged for an antiquark line.
c
c p: proVBFH order (1, 2 incoming; 3 Higgs; 4, 5 outgoing quarks of
c lines 1 and 2); bflav: PDG codes in the same order (bflav(3) ignored).
c Checked against the H + 3 parton Born of VBFNLO in the soft limit
c (cs_born2_check).
c-----------------------------------------------------------------------
      double precision function cs_born2(p,bflav)
      implicit none
      include 'koppln_ew.inc'
      double precision clr,xm2,xmg,b,v,a
      common /bkopou/ clr(4,5,-1:1),xm2(6),xmg(6),b(6,6,6),v(4,5),a(4,5)
      double precision p(0:3,5)
      integer bflav(5)
      integer t1,t2,v1,v2,h1,h2,ich1,ich2
      double precision in1(0:3),out1(0:3),in2(0:3),out2(0:3),s,b2dot
      double precision q1(0:3),q2(0:3),q1s,q2s,c
      double complex d1,d2
      external b2dot
      cs_born2 = 0
      if (bflav(1).eq.0 .or. bflav(2).eq.0) return
c     fermion-line types (3 up, 4 down) of the incoming (anti)quarks and
c     the charge they give to the boson (quark: flavour charge change
c     of the line)
      call line_boson(bflav(1),bflav(4),t1,ich1)
      call line_boson(bflav(2),bflav(5),t2,ich2)
      if (ich1+ich2.ne.0) return
      if (ich1.eq.0) then
         v1 = 2
         v2 = 2
      elseif (ich1.eq.1) then
         v1 = 3
         v2 = 4
      else
         v1 = 4
         v2 = 3
      endif
      if (bflav(1).gt.0) then
         in1 = p(:,1)
         out1 = p(:,4)
      else
         in1 = p(:,4)
         out1 = p(:,1)
      endif
      if (bflav(2).gt.0) then
         in2 = p(:,2)
         out2 = p(:,5)
      else
         in2 = p(:,5)
         out2 = p(:,2)
      endif
      q1 = p(:,1) - p(:,4)
      q2 = p(:,2) - p(:,5)
      q1s = b2dot(q1,q1)
      q2s = b2dot(q2,q2)
      d1 = 1/dcmplx(q1s-xm2(v1),xmg(v1))
      d2 = 1/dcmplx(q2s-xm2(v2),xmg(v2))
      do h1 = -1,1,2
         do h2 = -1,1,2
            c = clr(t1,v1,h1)*clr(t2,v2,h2)*b(6,v1,v2)*xmw
            if (c.eq.0) cycle
            if (h1.eq.h2) then
               s = b2dot(in1,in2)*b2dot(out1,out2)
            else
               s = b2dot(in1,out2)*b2dot(out1,in2)
            endif
            cs_born2 = cs_born2 + c**2*16*s
         enddo
      enddo
      cs_born2 = cs_born2*abs(d1)**2*abs(d2)**2/4
      end

c line of incoming a and outgoing b: VBFNLO fermion type t (3 up-type,
c 4 down-type) of the line's quark flavour (for W: of the incoming one)
c and the charge ich (0, +1, -1) of the boson it emits into the vertex
      subroutine line_boson(a,b,t,ich)
      implicit none
      integer a,b,t,ich,ca,cb
      ca = charge3(a)
      cb = charge3(b)
      ich = (ca - cb)/3
      if (mod(abs(a),2).eq.0) then
         t = 3
      else
         t = 4
      endif
      contains
      integer function charge3(f)
      integer f
      if (mod(abs(f),2).eq.0) then
         charge3 = 2
      else
         charge3 = -1
      endif
      if (f.lt.0) charge3 = -charge3
      end function charge3
      end

      double precision function b2dot(a,b)
      implicit none
      double precision a(0:3),b(0:3)
      b2dot = a(0)*b(0)-a(1)*b(1)-a(2)*b(2)-a(3)*b(3)
      end
