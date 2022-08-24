      module nonfact_expressions
      implicit none

      public
      contains

c.....................................................................
! Below are the roots r1 - r6 as defined in Lorenzo's logbook
      function r1(MV2, p1x, xi)
      real*8 MV2, p1x, xi
      complex*8 r1
      r1 = p1x*Cos(xi) - ((0,1)*
     -     Sqrt(2*MV2 + p1x**2 - p1x**2*Cos(2*xi)))/
     -   Sqrt(2.d0)
      end function

      function r2(MV2, p1x, xi)
      real*8 MV2, p1x, xi
      complex*8 r2
      r2 = p1x*Cos(xi) + ((0,1)*
     -     Sqrt(2*MV2 + p1x**2 - p1x**2*Cos(2*xi)))/
     -   Sqrt(2.d0)
      end function

      function r3(MV2, p2x, p2y, xi)
      real*8 MV2, p2x, p2y, xi
      complex*8 r3
      r3 = -(p2x*Cos(xi)) - p2y*Sin(xi) - 
     -  (0,0.5)*Sqrt(4*(MV2 + p2x**2 + p2y**2) - 
     -     4*(p2x*Cos(xi) + p2y*Sin(xi))**2)
      end function

      function r4(MV2, p2x, p2y, xi)
      real*8 MV2, p2x, p2y, xi
      complex*8 r4
      r4 = -(p2x*Cos(xi)) - p2y*Sin(xi) + 
     -  (0,0.5)*Sqrt(4*(MV2 + p2x**2 + p2y**2) - 
     -     4*(p2x*Cos(xi) + p2y*Sin(xi))**2)
      end function

      function r5(MVH2, p1x, p3x, p3y, xi)
      real*8 MVH2, p1x, p3x, p3y, xi
      complex*8 r5
      r5 = p1x*Cos(xi) + p3x*Cos(xi) + p3y*Sin(xi) - 
     -  (0,0.5)*Sqrt(4*
     -      (MVH2 + p1x**2 + 2*p1x*p3x + p3x**2 + p3y**2)
     -       - 4*((p1x + p3x)*Cos(xi) + p3y*Sin(xi))**2)
      end function

      function r6(MVH2, p1x, p3x, p3y, xi)
      real*8 MVH2, p1x, p3x, p3y, xi
      complex*8 r6
      r6 = p1x*Cos(xi) + p3x*Cos(xi) + p3y*Sin(xi) + 
     -  (0,0.5)*Sqrt(4*
     -      (MVH2 + p1x**2 + 2*p1x*p3x + p3x**2 + p3y**2)
     -       - 4*((p1x + p3x)*Cos(xi) + p3y*Sin(xi))**2)
      end function

! 1-loop Triangle
      function t01(MV2, pi, p1x, p2x, p2y, xi)
      real*8 MV, MV2, pi, p1x, p2x, p2y, xi
      complex*8 r1v, r2v, r3v, r4v
      real*8 t01
      MV = sqrt(MV2)
      r1v=r1(MV2, p1x, xi)
      r2v=r2(MV2, p1x, xi)
      r3v=r3(MV2, p2x, p2y, xi)
      r4v=r4(MV2, p2x, p2y, xi)

      t01 = -(Log(-(r1v/Mv))/(Pi*r1v*(r1v - r2v)*(r1v - r3v)*(r1v -
     $     r4v))) + Log(-(r3v/Mv))/(Pi*(r2v - r3v)*r3v*(-r1v + r3v)*(r3v
     $     - r4v))
      end function

      function t11(MV2, pi, p1x, p2x, p2y)
      real*8 MV2, pi, p1x, p2x, p2y
      real*8 t11

      t11 = -(1/((MV2 + p1x**2)*(MV2 + p2x**2 + p2y**2)))
      end
      
!     1-loop Box
      function b01(MV2, MVH2, pi, p1x, p2x, p2y, p3x, p3y,xi)
      real*8 MV, MV2, MVH2, pi, p1x, p2x, p2y, p3x, p3y, xi
      complex*8 r1v, r2v, r3v, r4v, r5v, r6v
      real*8 b01
      MV = sqrt(MV2)
      r1v=r1(MV2, p1x, xi)
      r2v=r2(MV2, p1x, xi)
      r3v=r3(MV2, p2x, p2y, xi)
      r4v=r4(MV2, p2x, p2y, xi)
      r5v=r5(MVH2, p1x, p3x, p3y, xi)
      r6v=r6(MVH2, p1x, p3x, p3y, xi)

      b01 = -(Log(-(r1v/Mv))/(Pi*r1v*(r1v - r2v)*(r1v - r3v)*(r1v - r4v)
     $     *(r1v - r5v)*(r1v - r6v))) +Log(-(r3v/Mv))/(Pi*(r2v - r3v)
     $     *r3v*(-r1v + r3v)*(r3v - r4v)*(r3v - r5v)*(r3v - r6v)) +Log(
     $     -(r5v/Mv))/(Pi*(r2v - r5v)*r5v*(-r1v + r5v)*(-r3v + r5v)*(
     $     -r4v + r5v)*(r5v - r6v))
      end function

      function b11(MV2, MVH2, pi, p1x, p2x, p2y, p3x, p3y)
      real*8 MV2, MVH2, pi, p1x, p2x, p2y, p3x, p3y
      real*8 b11
      b11 = -(1/((MV2 + p1x**2)*(MV2 + p2x**2 + p2y**2)*(MVH2 + p1x**2 +
     $     2*p1x*p3x + p3x**2 + p3y**2)))
      end function

!     2-loop Triangle expressions
      
      function t012(MV2, pi, p1x, p2x, p2y)
      real*8 MV2, pi, p1x, p2x, p2y
      real*8 t012
      t012 = (4*Pi**2)/(3D0*(MV2 + p1x**2)*(MV2 + p2x**2 + p2y**2))
      end function
      
      function t022(MV2, pi, p1x, p2x, p2y, xi)
      real*8 MV, MV2, pi, p1x, p2x, p2y, xi
      complex*8 r1v, r2v, r3v, r4v
      real*8 t022
      MV = sqrt(MV2)
      r1v=r1(MV2, p1x, xi)
      r2v=r2(MV2, p1x, xi)
      r3v=r3(MV2, p2x, p2y, xi)
      r4v=r4(MV2, p2x, p2y, xi)
      t022 = (-2*Log(-(r1v/Mv))**2)/(Pi*r1v*(r1v - r2v)*(r1v - r3v)*(r1v
     $     - r4v)) -(2*Log(-(r3v/Mv))**2)/(Pi*r3v*(-r1v + r3v)*(-r2v +
     $     r3v)*(r3v - r4v))
      end function
      
      function t12(MV2, pi, p1x, p2x, p2y, xi)
      real*8 MV, MV2, pi, p1x, p2x, p2y, xi
      complex*8 r1v, r2v, r3v, r4v
      real*8 t12
      MV = sqrt(MV2)
      r1v=r1(MV2, p1x, xi)
      r2v=r2(MV2, p1x, xi)
      r3v=r3(MV2, p2x, p2y, xi)
      r4v=r4(MV2, p2x, p2y, xi)
      t12 = (2*Log(-(r1v/Mv)))/(Pi*r1v*(r1v - r2v)*(r1v - r3v)*(r1v -
     $     r4v)) +(2*Log(-(r3v/Mv)))/(Pi*r3v*(-r1v + r3v)*(-r2v + r3v)
     $     *(r3v - r4v))
      end function
      
      function t22(MV2, pi, p1x, p2x, p2y)
      real*8 MV2, pi, p1x, p2x, p2y
      real*8 t22
      t22 = 1D0/((MV2 + p1x**2)*(MV2 + p2x**2 + p2y**2))
      end function

      ! 2-loop Box expressions

      function b012(MV2, MVH2, pi, p1x, p2x, p2y, p3x, p3y)
      real*8 MV2, MVH2, pi, p1x, p2x, p2y, p3x, p3y, xi
      real*8 b012
      b012 = (4*Pi**2)/(3.*(MV2 + p1x**2)*(MV2 + p2x**2 + p2y**2)*(MVH2 +
     $     (p1x + p3x)**2 + p3y**2))
      end function

      function b022(MV2, MVH2, pi, p1x, p2x, p2y, p3x, p3y, xi)
      real*8 MV, MV2, MVH2, pi, p1x, p2x, p2y, p3x, p3y, xi
      complex*8 r1v, r2v, r3v, r4v, r5v, r6v
      complex*8 b022
      MV = sqrt(MV2)
      r1v=r1(MV2, p1x, xi)
      r2v=r2(MV2, p1x, xi)
      r3v=r3(MV2, p2x, p2y, xi)
      r4v=r4(MV2, p2x, p2y, xi)
      r5v=r5(MVH2, p1x, p3x, p3y, xi)
      r6v=r6(MVH2, p1x, p3x, p3y, xi)
      
      b022 =  (-2d0*Log(-(r1v/Mv))**2)/ (Pi*r1v*(r1v - r2v)*(r1v - r3v)
     $     *(r1v - r4v)*(r1v - r5v)*(r1v - r6v)) - (2d0*Log(-(r3v/Mv))
     $     **2) / (Pi*r3v*(-r1v + r3v)*(-r2v + r3v)*(r3v - r4v)*(r3v -
     $     r5v) *(r3v - r6v)) - (2d0*Log(-(r5v/Mv))**2)/ (Pi*r5v*(-r1v +
     $     r5v)*( -r2v + r5v)*(-r3v + r5v)*(-r4v + r5v)*(r5v - r6v))
      end function

      function b12(MV2, MVH2, pi, p1x, p2x, p2y, p3x, p3y, xi)
      real*8 MV, MV2, MVH2, pi, p1x, p2x, p2y, p3x, p3y, xi
      complex*8 r1v, r2v, r3v, r4v, r5v, r6v
      complex*8 b12
      MV = sqrt(MV2)
      r1v=r1(MV2, p1x, xi)
      r2v=r2(MV2, p1x, xi)
      r3v=r3(MV2, p2x, p2y, xi)
      r4v=r4(MV2, p2x, p2y, xi)
      r5v=r5(MVH2, p1x, p3x, p3y, xi)
      r6v=r6(MVH2, p1x, p3x, p3y, xi)
      b12 = (2d0*Log(-(r1v/Mv)))/(Pi*r1v*(r1v - r2v)*(r1v - r3v)*(r1v -
     $     r4v)*(r1v - r5v)*(r1v - r6v)) +(2d0*Log(-(r3v/Mv)))/(Pi*r3v*(
     $     -r1v + r3v)*(-r2v + r3v)*(r3v - r4v)*(r3v - r5v)*(r3v - r6v))
     $     +(2d0*Log(-(r5v/Mv)))/(Pi*r5v*(-r1v + r5v)*(-r2v + r5v)*(-r3v
     $     + r5v)*(-r4v + r5v)*(r5v - r6v))
      end function
      
      function b22(MV2, MVH2, pi, p1x, p2x, p2y, p3x, p3y)
      real*8 MV2, MVH2, pi, p1x, p2x, p2y, p3x, p3y
      real*8 b22
      b22 = 1d0/((MV2 + p1x**2)*(MV2 + p2x**2 + p2y**2)*(MVH2 + p1x**2 +
     $     2d0*p1x*p3x + p3x**2 + p3y**2))
      end function

c......................................................................
      end module nonfact_expressions
