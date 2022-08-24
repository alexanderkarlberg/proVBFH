      module ME_expressions
      implicit none

      public
      private::dot
      contains

!----------------------------------------------------------------------
!     dot product
      function dot(p1,p2)
      real*8 p1(0:3), p2(0:3)
      real*8 dot
      dot = p1(0)*p2(0) - sum(p1(1:3)*p2(1:3))
      end function dot

!----------------------------------------------------------------------
!     F1F1 
      function F1F1(F1F1AA, F1F1AB, F1F1AC, F1F1BB, F1F1BC, F1F1CC,
     &     q1, q2, P1, P2, k1, k2)
      real*8 F1F1AA, F1F1AB, F1F1AC, F1F1BB, F1F1BC, F1F1CC
      real*8 q1(0:3),q2(0:3),P1(0:3),P2(0:3),k1(0:3),k2(0:3)
      real*8 F1F1
      ! copied from mathematica output, formatted with FortranForm


      F1F1 = (4*(F1F1BB*dot(k1,k1) + F1F1CC*dot(k2,k2))*dot(q1,q1)**3*dot(q2,q2)**2 + 
     -    dot(q1,q1)*dot(q2,q2)*(4*F1F1CC*dot(k1,q2)**2*dot(k2,q1)**2 + 
     -       4*F1F1CC*dot(k2,q1)**2*dot(k2,q2)**2 + 8*F1F1CC*dot(k2,q1)**2*dot(k2,q2)*dot(q1,q2) + 
     -       F1F1AA*dot(q1,q2)**2 + 4*F1F1CC*dot(k2,q1)**2*dot(q1,q2)**2 - 
     -       8*F1F1CC*dot(k1,q2)*dot(k2,q1)**2*(dot(k2,q2) + dot(q1,q2)) - 
     -       8*F1F1BB*dot(k1,q1)**3*dot(q2,q2) + 8*F1F1CC*dot(k1,q1)*dot(k2,q1)**2*dot(q2,q2) - 
     -       4*F1F1CC*dot(k2,q1)**2*(dot(k1,k1) - 2*dot(k1,k2) + dot(k2,k2) + 2*dot(k2,q1))*dot(q2,q2) + 
     -       4*F1F1BB*dot(k1,q1)**2*((dot(k1,q2) - dot(k2,q2) + dot(q1,q2))**2 - 
     -          (dot(k1,k1) - 2*dot(k1,k2) + dot(k2,k2) - 2*dot(k2,q1))*dot(q2,q2))) - 
     -    2*dot(q1,q1)*(2*dot(q1,q1)*(F1F1BB*dot(k1,k1)*(dot(k1,q2) - dot(k2,q2) + dot(q1,q2))**2 + 
     -          F1F1CC*dot(k2,k2)*(-dot(k1,q2) + dot(k2,q2) + dot(q1,q2))**2)*dot(q2,q2) + 
     -       dot(q2,q2)*(4*F1F1BC*dot(k1,q1)*dot(k2,q1)*(dot(k1,q2) - dot(k2,q2))**2 - 
     -          2*(dot(k1,q2) - dot(k2,q2))*
     -           ((F1F1AB + 2*F1F1BC*dot(k1,k2))*dot(k1,q2) - (F1F1AC + 2*F1F1BC*dot(k1,k2))*dot(k2,q2))*
     -           dot(q1,q1) + 2*((F1F1AB*dot(k1,q1) - F1F1AC*dot(k2,q1))*(dot(k1,q2) - dot(k2,q2)) - 
     -             (F1F1AB*dot(k1,q2) + F1F1AC*dot(k2,q2))*dot(q1,q1))*dot(q1,q2) + 
     -          2*(F1F1AC*dot(k2,q1) + dot(k1,q1)*(F1F1AB - 2*F1F1BC*dot(k2,q1)) + 
     -             2*F1F1BC*dot(k1,k2)*dot(q1,q1))*dot(q1,q2)**2 - 
     -          (2*F1F1AB*dot(k1,q1)**2 + 2*F1F1AC*dot(k2,q1)**2 - 
     -             2*dot(k1,q1)*dot(k2,q1)*
     -              (F1F1AB + F1F1AC - 2*F1F1BC*(dot(k1,k1) - 2*dot(k1,k2) + dot(k2,k2) - dot(q1,q1))) + 
     -             (F1F1AA - 2*dot(k1,k1)*(F1F1AB + 2*F1F1BC*dot(k1,k2)) + 
     -                2*dot(k1,k2)*(F1F1AB + F1F1AC + 4*F1F1BC*dot(k1,k2)) - 
     -                2*(F1F1AC + 2*F1F1BC*dot(k1,k2))*dot(k2,k2) + 
     -                2*(F1F1BB*dot(k1,k1)**2 - F1F1BB*dot(k1,q1)**2 - 2*F1F1CC*dot(k1,q1)*dot(k2,k2) + 
     -                   dot(k1,k1)*(-2*F1F1BB*dot(k1,k2) + 2*F1F1BB*dot(k1,q1) + 
     -                      (F1F1BB + F1F1CC)*dot(k2,k2) - 2*F1F1BB*dot(k2,q1)) + 
     -                   F1F1CC*(-2*dot(k1,k2)*dot(k2,k2) + dot(k2,k2)**2 + 2*dot(k2,k2)*dot(k2,q1) - 
     -                      dot(k2,q1)**2)))*dot(q1,q1) + 4*F1F1BC*dot(k1,k2)*dot(q1,q1)**2)*dot(q2,q2)))
     -    )/(dot(q1,q1)**2*dot(q2,q2)**2)
      
      end function

!----------------------------------------------------------------------
!     F1F2
      function F1F2(F1F2AA, F1F2AB, F1F2AC, F1F2BB, F1F2BC, F1F2CC,
     &     q1, q2, P1, P2, k1, k2)
      real*8 F1F2AA, F1F2AB, F1F2AC, F1F2BB, F1F2BC, F1F2CC
      real*8 q1(0:3),q2(0:3),P1(0:3),P2(0:3),k1(0:3),k2(0:3)
      real*8 F1F2

      ! copied from mathematica output, formatted with FortranForm
      F1F2 = (dot(P2,q2)*dot(q1,q1)*(4*F1F2BB*dot(k1,P2)**2*dot(k1,q1)**2 - 
     -       8*F1F2BB*dot(k1,P2)*dot(k1,q1)**2*dot(k2,P2) + 
     -       4*F1F2BB*dot(k1,q1)**2*dot(k2,P2)**2 + 4*F1F2CC*dot(k1,P2)**2*dot(k2,q1)**2 - 
     -       8*F1F2CC*dot(k1,P2)*dot(k2,P2)*dot(k2,q1)**2 + 
     -       4*F1F2CC*dot(k2,P2)**2*dot(k2,q1)**2 + 
     -       8*F1F2BB*dot(k1,P2)*dot(k1,q1)**2*dot(P2,q1) - 
     -       8*F1F2BB*dot(k1,q1)**2*dot(k2,P2)*dot(P2,q1) - 
     -       8*F1F2CC*dot(k1,P2)*dot(k2,q1)**2*dot(P2,q1) + 
     -       8*F1F2CC*dot(k2,P2)*dot(k2,q1)**2*dot(P2,q1) + F1F2AA*dot(P2,q1)**2 + 
     -       4*F1F2BB*dot(k1,q1)**2*dot(P2,q1)**2 + 4*F1F2CC*dot(k2,q1)**2*dot(P2,q1)**2 - 
     -       (F1F2AA*dot(P2,P2) + 4*F1F2BB*dot(k1,k1)*
     -           (dot(k1,P2) - dot(k2,P2) + dot(P2,q1))**2 + 
     -          4*F1F2CC*dot(k2,k2)*(-dot(k1,P2) + dot(k2,P2) + dot(P2,q1))**2)*dot(q1,q1))
     -      *dot(q2,q2)**4 + dot(P2,q2)**3*dot(q1,q1)*dot(q2,q2)**2*
     -     (4*F1F2CC*dot(k2,q1)**2*dot(k2,q2)**2 - 
     -       4*F1F2BB*dot(k1,k1)*dot(k2,q2)**2*dot(q1,q1) - 
     -       4*F1F2CC*dot(k2,k2)*dot(k2,q2)**2*dot(q1,q1) + 
     -       4*dot(k1,q2)**2*(F1F2CC*dot(k2,q1)**2 - 
     -          (F1F2BB*dot(k1,k1) + F1F2CC*dot(k2,k2))*dot(q1,q1)) + 
     -       8*F1F2CC*dot(k2,q1)**2*dot(k2,q2)*dot(q1,q2) + 
     -       8*F1F2BB*dot(k1,k1)*dot(k2,q2)*dot(q1,q1)*dot(q1,q2) - 
     -       8*F1F2CC*dot(k2,k2)*dot(k2,q2)*dot(q1,q1)*dot(q1,q2) + F1F2AA*dot(q1,q2)**2 + 
     -       4*F1F2CC*dot(k2,q1)**2*dot(q1,q2)**2 - 
     -       4*F1F2BB*dot(k1,k1)*dot(q1,q1)*dot(q1,q2)**2 - 
     -       4*F1F2CC*dot(k2,k2)*dot(q1,q1)*dot(q1,q2)**2 + 
     -       4*F1F2BB*dot(k1,q1)**2*(dot(k1,q2) - dot(k2,q2) + dot(q1,q2))**2 + 
     -       8*dot(k1,q2)*(-(F1F2CC*dot(k2,q1)**2*(dot(k2,q2) + dot(q1,q2))) + 
     -          dot(q1,q1)*(F1F2BB*dot(k1,k1)*(dot(k2,q2) - dot(q1,q2)) + 
     -             F1F2CC*dot(k2,k2)*(dot(k2,q2) + dot(q1,q2)))) + 
     -       F1F2AA*dot(q1,q1)*dot(q2,q2)) - 
     -    2*dot(P2,q2)*dot(q2,q2)*(dot(P2,q2)*dot(q1,q1)*
     -        (4*(dot(k1,q2) - dot(k2,q2))*
     -           ((dot(k1,P2) - dot(k2,P2))*
     -              (F1F2BB*dot(k1,q1)**2 + F1F2CC*dot(k2,q1)**2) + 
     -             (F1F2BB*dot(k1,q1)**2 - F1F2CC*dot(k2,q1)**2)*dot(P2,q1)) + 
     -          (4*(dot(k1,P2) - dot(k2,P2))*
     -              (F1F2BB*dot(k1,q1)**2 - F1F2CC*dot(k2,q1)**2) + 
     -             (F1F2AA + 4*F1F2BB*dot(k1,q1)**2 + 4*F1F2CC*dot(k2,q1)**2)*dot(P2,q1))*
     -           dot(q1,q2))*dot(q2,q2)**2 + 
     -       2*dot(q1,q1)*(-2*dot(P2,q2)*dot(q1,q1)*
     -           (F1F2CC*dot(k2,k2)*(dot(k1,P2) - dot(k2,P2) - dot(P2,q1))*
     -              (dot(k1,q2) - dot(k2,q2) - dot(q1,q2)) + 
     -             F1F2BB*dot(k1,k1)*(dot(k1,P2) - dot(k2,P2) + dot(P2,q1))*
     -              (dot(k1,q2) - dot(k2,q2) + dot(q1,q2)))*dot(q2,q2)**2 + 
     -          dot(q2,q2)*(dot(P2,q2)**2*
     -              (-((dot(k1,q2) - dot(k2,q2))*
     -                   ((F1F2AB + 2*F1F2BC*dot(k1,k2))*dot(k1,q2) - 
     -                     (F1F2AC + 2*F1F2BC*dot(k1,k2))*dot(k2,q2))*dot(q1,q1)) - 
     -                (F1F2AC*dot(k2,q2)*(-dot(k2,q1) + dot(q1,q1)) + 
     -                   dot(k1,q2)*(F1F2AC*dot(k2,q1) + F1F2AB*dot(q1,q1)))*dot(q1,q2) + 
     -                (F1F2AC*dot(k2,q1) + 2*F1F2BC*dot(k1,k2)*dot(q1,q1))*dot(q1,q2)**2)
     -              + dot(P2,q2)*(dot(k1,q2)*
     -                 (F1F2AC*dot(k2,q1)*dot(P2,q1) + 
     -                   (2*(F1F2AB + 2*F1F2BC*dot(k1,k2))*dot(k1,P2) - 
     -                      (F1F2AB + F1F2AC + 4*F1F2BC*dot(k1,k2))*dot(k2,P2) + 
     -                      F1F2AB*dot(P2,q1))*dot(q1,q1)) + 
     -                dot(k2,q2)*(-(F1F2AC*dot(k2,q1)*dot(P2,q1)) + 
     -                   (-((F1F2AB + F1F2AC + 4*F1F2BC*dot(k1,k2))*dot(k1,P2)) + 
     -                      2*(F1F2AC + 2*F1F2BC*dot(k1,k2))*dot(k2,P2) + F1F2AC*dot(P2,q1)
     -                      )*dot(q1,q1)) + 
     -                (F1F2AC*dot(k2,q1)*(dot(k1,P2) - dot(k2,P2) - 2*dot(P2,q1)) + 
     -                   (F1F2AB*dot(k1,P2) + F1F2AC*dot(k2,P2) - 
     -                      4*F1F2BC*dot(k1,k2)*dot(P2,q1))*dot(q1,q1))*dot(q1,q2))*
     -              dot(q2,q2) - (F1F2AC*dot(k2,q1)*(dot(k1,P2) - dot(k2,P2) - dot(P2,q1))*
     -                 dot(P2,q1) + ((dot(k1,P2) - dot(k2,P2))*
     -                    ((F1F2AB + 2*F1F2BC*dot(k1,k2))*dot(k1,P2) - 
     -                      (F1F2AC + 2*F1F2BC*dot(k1,k2))*dot(k2,P2)) + 
     -                   (F1F2AB*dot(k1,P2) + F1F2AC*dot(k2,P2))*dot(P2,q1) - 
     -                   2*F1F2BC*dot(k1,k2)*dot(P2,q1)**2)*dot(q1,q1))*dot(q2,q2)**2 + 
     -             dot(k1,q1)*(dot(P2,q2)*(dot(k1,q2) - dot(k2,q2) + dot(q1,q2)) - 
     -                (dot(k1,P2) - dot(k2,P2) + dot(P2,q1))*dot(q2,q2))*
     -              (dot(P2,q2)*(2*F1F2BC*dot(k2,q1)*
     -                    (dot(k1,q2) - dot(k2,q2) - dot(q1,q2)) + F1F2AB*dot(q1,q2)) + 
     -                (-(F1F2AB*dot(P2,q1)) + 
     -                   2*F1F2BC*dot(k2,q1)*(-dot(k1,P2) + dot(k2,P2) + dot(P2,q1)))*
     -                 dot(q2,q2))))))/(dot(P2,q2)**2*dot(q1,q1)**2*dot(q2,q2)**4)


      end function
      
!----------------------------------------------------------------------
!     F2F1
      function F2F1(F2F1AA, F2F1AB, F2F1AC, F2F1BB, F2F1BC, F2F1CC,
     &     q1, q2, P1, P2, k1, k2)
      real*8 F2F1AA, F2F1AB, F2F1AC, F2F1BB, F2F1BC, F2F1CC
      real*8 q1(0:3),q2(0:3),P1(0:3),P2(0:3),k1(0:3),k2(0:3)
      real*8 F2F1

      ! copied from mathematica output, formatted with FortranForm
      F2F1 = (dot(P1,q1)**3*dot(q1,q1)**2*dot(q2,q2)*
     -     (4*F2F1CC*dot(k1,q2)**2*dot(k2,q1)**2 + 4*F2F1CC*dot(k2,q1)**2*dot(k2,q2)**2 + 
     -       8*F2F1CC*dot(k2,q1)**2*dot(k2,q2)*dot(q1,q2) + F2F1AA*dot(q1,q2)**2 + 
     -       4*F2F1CC*dot(k2,q1)**2*dot(q1,q2)**2 - 
     -       8*F2F1CC*dot(k1,q2)*dot(k2,q1)**2*(dot(k2,q2) + dot(q1,q2)) - 
     -       8*F2F1BB*dot(k1,q1)**3*dot(q2,q2) + 
     -       8*F2F1CC*dot(k1,q1)*dot(k2,q1)**2*dot(q2,q2) + 
     -       (F2F1AA*dot(q1,q1) - 4*F2F1CC*dot(k2,q1)**2*
     -           (dot(k1,k1) - 2*dot(k1,k2) + dot(k2,k2) + 2*dot(k2,q1) + dot(q1,q1)))*
     -        dot(q2,q2) + 4*F2F1BB*dot(k1,q1)**2*
     -        ((dot(k1,q2) - dot(k2,q2) + dot(q1,q2))**2 - 
     -          (dot(k1,k1) - 2*dot(k1,k2) + dot(k2,k2) - 2*dot(k2,q1) + dot(q1,q1))*
     -           dot(q2,q2))) + dot(P1,q1)*dot(q1,q1)**4*dot(q2,q2)*
     -     (F2F1AA*dot(P1,q2)**2 + 4*F2F1CC*dot(k2,P1)**2*
     -        (-dot(k1,q2) + dot(k2,q2) + dot(q1,q2))**2 - 
     -       (F2F1AA*dot(P1,P1) + 4*F2F1CC*dot(k2,P1)**2*
     -           (dot(k1,k1) - 2*dot(k1,k2) - 2*dot(k1,q1) + dot(k2,k2) + 2*dot(k2,q1) + 
     -             dot(q1,q1)))*dot(q2,q2) + 
     -       4*F2F1BB*dot(k1,P1)**2*((dot(k1,q2) - dot(k2,q2) + dot(q1,q2))**2 - 
     -          (dot(k1,k1) - 2*dot(k1,k2) + 2*dot(k1,q1) + dot(k2,k2) - 2*dot(k2,q1) + 
     -             dot(q1,q1))*dot(q2,q2))) - 
     -    2*dot(P1,q1)*dot(q1,q1)*(dot(P1,q1)*dot(q1,q1)**2*dot(q2,q2)*
     -        (4*F2F1CC*dot(k1,q2)**2*dot(k2,P1)*dot(k2,q1) + 
     -          4*F2F1CC*dot(k2,P1)*dot(k2,q1)*dot(k2,q2)**2 + 
     -          8*F2F1CC*dot(k2,P1)*dot(k2,q1)*dot(k2,q2)*dot(q1,q2) + 
     -          F2F1AA*dot(P1,q2)*dot(q1,q2) + 
     -          4*F2F1CC*dot(k2,P1)*dot(k2,q1)*dot(q1,q2)**2 - 
     -          8*F2F1CC*dot(k1,q2)*dot(k2,P1)*dot(k2,q1)*(dot(k2,q2) + dot(q1,q2)) - 
     -          4*F2F1CC*dot(k2,P1)*dot(k2,q1)*
     -           (dot(k1,k1) - 2*dot(k1,k2) - 2*dot(k1,q1) + dot(k2,k2) + 2*dot(k2,q1))*
     -           dot(q2,q2) + 4*F2F1BB*dot(k1,P1)*dot(k1,q1)*
     -           ((dot(k1,q2) - dot(k2,q2) + dot(q1,q2))**2 - 
     -             (dot(k1,k1) - 2*dot(k1,k2) + 2*dot(k1,q1) + dot(k2,k2) - 2*dot(k2,q1))*
     -              dot(q2,q2))) - 2*dot(q1,q1)*dot(q2,q2)*
     -        (F2F1AC*(dot(k2,q1)*dot(P1,q1) - dot(k2,P1)*dot(q1,q1))*
     -           (dot(k1,q2) - dot(k2,q2) - dot(q1,q2))*
     -           (-(dot(P1,q2)*dot(q1,q1)) + dot(P1,q1)*dot(q1,q2)) + 
     -          F2F1AB*dot(k1,q1)**2*dot(P1,q1)**2*dot(q2,q2) + 
     -          F2F1AB*dot(k1,P1)**2*dot(q1,q1)**2*dot(q2,q2) + 
     -          (F2F1AC*dot(k2,q1)**2*dot(P1,q1)**2 - 
     -             2*F2F1AC*dot(k2,P1)*dot(k2,q1)*dot(P1,q1)*dot(q1,q1) + 
     -             (F2F1AC*dot(k2,P1)**2 + 
     -                2*(F2F1BB*dot(k1,P1)*dot(k1,q1) + F2F1CC*dot(k2,P1)*dot(k2,q1))*
     -                 dot(P1,q1))*dot(q1,q1)**2)*dot(q2,q2) + 
     -          dot(k1,P1)*dot(q1,q1)*
     -           ((dot(k1,q2) - dot(k2,q2) + dot(q1,q2))*
     -              ((2*F2F1BC*dot(k2,P1)*dot(k2,q2) - F2F1AB*dot(P1,q2))*dot(q1,q1) + 
     -                2*F2F1BC*dot(k1,q2)*
     -                 (dot(k2,q1)*dot(P1,q1) - dot(k2,P1)*dot(q1,q1)) + 
     -                (F2F1AB*dot(P1,q1) + 2*F2F1BC*dot(k2,P1)*dot(q1,q1))*dot(q1,q2) - 
     -                2*F2F1BC*dot(k2,q1)*dot(P1,q1)*(dot(k2,q2) + dot(q1,q2))) + 
     -             (F2F1AB + F2F1AC - 
     -                2*F2F1BC*(dot(k1,k1) - 2*dot(k1,k2) + dot(k2,k2) - dot(q1,q1)))*
     -              (dot(k2,q1)*dot(P1,q1) - dot(k2,P1)*dot(q1,q1))*dot(q2,q2)) + 
     -          dot(k1,q1)*dot(P1,q1)*
     -           (-((dot(k1,q2) - dot(k2,q2) + dot(q1,q2))*
     -                ((2*F2F1BC*dot(k2,P1)*dot(k2,q2) - F2F1AB*dot(P1,q2))*dot(q1,q1) + 
     -                  2*F2F1BC*dot(k1,q2)*
     -                   (dot(k2,q1)*dot(P1,q1) - dot(k2,P1)*dot(q1,q1)) + 
     -                  (F2F1AB*dot(P1,q1) + 2*F2F1BC*dot(k2,P1)*dot(q1,q1))*dot(q1,q2) - 
     -                  2*F2F1BC*dot(k2,q1)*dot(P1,q1)*(dot(k2,q2) + dot(q1,q2)))) + 
     -             (-(dot(k2,q1)*dot(P1,q1)*
     -                   (F2F1AB + F2F1AC - 
     -                     2*F2F1BC*(dot(k1,k1) - 2*dot(k1,k2) + dot(k2,k2) - dot(q1,q1))))
     -                  + (-2*F2F1AB*dot(k1,P1) + 
     -                   dot(k2,P1)*(F2F1AB + F2F1AC - 
     -                      2*F2F1BC*(dot(k1,k1) - 2*dot(k1,k2) + dot(k2,k2) - dot(q1,q1)))
     -                   )*dot(q1,q1))*dot(q2,q2)))))/
     -  (dot(P1,q1)**2*dot(q1,q1)**4*dot(q2,q2)**2)
      
      end function

!----------------------------------------------------------------------
!     F2F2
      function F2F2(F2F2AA, F2F2AB, F2F2AC, F2F2BB, F2F2BC, F2F2CC,
     &     q1, q2, P1, P2, k1, k2)
      real*8 F2F2AA, F2F2AB, F2F2AC, F2F2BB, F2F2BC, F2F2CC
      real*8 q1(0:3),q2(0:3),P1(0:3),P2(0:3),k1(0:3),k2(0:3)
      real*8 F2F2

      ! copied from mathematica output, formatted with FortranForm
      F2F2 =        (dot(P1,q1)*dot(P2,q2)*dot(q1,q1)**4*dot(q2,q2)**2*
     -     (4*F2F2CC*dot(k1,q2)**2*dot(k2,P1)**2*dot(P2,q2)**2 + 
     -       4*F2F2CC*dot(k2,P1)**2*dot(k2,q2)**2*dot(P2,q2)**2 + 
     -       F2F2AA*dot(P1,q2)**2*dot(P2,q2)**2 + 
     -       8*F2F2CC*dot(k2,P1)**2*dot(k2,q2)*dot(P2,q2)**2*dot(q1,q2) + 
     -       4*F2F2CC*dot(k2,P1)**2*dot(P2,q2)**2*dot(q1,q2)**2 + 
     -       8*F2F2CC*dot(k1,P2)*dot(k2,P1)**2*dot(k2,q2)*dot(P2,q2)*dot(q2,q2) - 
     -       8*F2F2CC*dot(k2,P1)**2*dot(k2,P2)*dot(k2,q2)*dot(P2,q2)*dot(q2,q2) - 
     -       2*F2F2AA*dot(P1,P2)*dot(P1,q2)*dot(P2,q2)*dot(q2,q2) - 
     -       8*F2F2CC*dot(k2,P1)**2*dot(k2,q2)*dot(P2,q1)*dot(P2,q2)*dot(q2,q2) + 
     -       8*F2F2CC*dot(k1,P2)*dot(k2,P1)**2*dot(P2,q2)*dot(q1,q2)*dot(q2,q2) - 
     -       8*F2F2CC*dot(k2,P1)**2*dot(k2,P2)*dot(P2,q2)*dot(q1,q2)*dot(q2,q2) - 
     -       8*F2F2CC*dot(k2,P1)**2*dot(P2,q1)*dot(P2,q2)*dot(q1,q2)*dot(q2,q2) + 
     -       4*F2F2CC*dot(k1,P2)**2*dot(k2,P1)**2*dot(q2,q2)**2 - 
     -       8*F2F2CC*dot(k1,P2)*dot(k2,P1)**2*dot(k2,P2)*dot(q2,q2)**2 + 
     -       4*F2F2CC*dot(k2,P1)**2*dot(k2,P2)**2*dot(q2,q2)**2 + 
     -       F2F2AA*dot(P1,P2)**2*dot(q2,q2)**2 - 
     -       8*F2F2CC*dot(k1,P2)*dot(k2,P1)**2*dot(P2,q1)*dot(q2,q2)**2 + 
     -       8*F2F2CC*dot(k2,P1)**2*dot(k2,P2)*dot(P2,q1)*dot(q2,q2)**2 + 
     -       4*F2F2CC*dot(k2,P1)**2*dot(P2,q1)**2*dot(q2,q2)**2 - 
     -       8*F2F2CC*dot(k1,q2)*dot(k2,P1)**2*dot(P2,q2)*
     -        (dot(k2,q2)*dot(P2,q2) + dot(P2,q2)*dot(q1,q2) + 
     -          (dot(k1,P2) - dot(k2,P2) - dot(P2,q1))*dot(q2,q2)) + 
     -       4*F2F2BB*dot(k1,P1)**2*(dot(k1,q2)*dot(P2,q2) - dot(k2,q2)*dot(P2,q2) + 
     -           dot(P2,q2)*dot(q1,q2) - dot(k1,P2)*dot(q2,q2) + dot(k2,P2)*dot(q2,q2) - 
     -           dot(P2,q1)*dot(q2,q2))**2) + 
     -    dot(P1,q1)**3*dot(P2,q2)*dot(q1,q1)**2*dot(q2,q2)**2*
     -     (4*F2F2CC*dot(k1,q2)**2*dot(k2,q1)**2*dot(P2,q2)**2 + 
     -       4*F2F2CC*dot(k2,q1)**2*dot(k2,q2)**2*dot(P2,q2)**2 + 
     -       8*F2F2CC*dot(k2,q1)**2*dot(k2,q2)*dot(P2,q2)**2*dot(q1,q2) + 
     -       F2F2AA*dot(P2,q2)**2*dot(q1,q2)**2 + 
     -       4*F2F2CC*dot(k2,q1)**2*dot(P2,q2)**2*dot(q1,q2)**2 + 
     -       8*F2F2CC*dot(k1,P2)*dot(k2,q1)**2*dot(k2,q2)*dot(P2,q2)*dot(q2,q2) - 
     -       8*F2F2CC*dot(k2,P2)*dot(k2,q1)**2*dot(k2,q2)*dot(P2,q2)*dot(q2,q2) - 
     -       8*F2F2CC*dot(k2,q1)**2*dot(k2,q2)*dot(P2,q1)*dot(P2,q2)*dot(q2,q2) + 
     -       8*F2F2CC*dot(k1,P2)*dot(k2,q1)**2*dot(P2,q2)*dot(q1,q2)*dot(q2,q2) - 
     -       8*F2F2CC*dot(k2,P2)*dot(k2,q1)**2*dot(P2,q2)*dot(q1,q2)*dot(q2,q2) - 
     -       2*F2F2AA*dot(P2,q1)*dot(P2,q2)*dot(q1,q2)*dot(q2,q2) - 
     -       8*F2F2CC*dot(k2,q1)**2*dot(P2,q1)*dot(P2,q2)*dot(q1,q2)*dot(q2,q2) + 
     -       4*F2F2CC*dot(k1,P2)**2*dot(k2,q1)**2*dot(q2,q2)**2 - 
     -       8*F2F2CC*dot(k1,P2)*dot(k2,P2)*dot(k2,q1)**2*dot(q2,q2)**2 + 
     -       4*F2F2CC*dot(k2,P2)**2*dot(k2,q1)**2*dot(q2,q2)**2 - 
     -       8*F2F2CC*dot(k1,P2)*dot(k2,q1)**2*dot(P2,q1)*dot(q2,q2)**2 + 
     -       8*F2F2CC*dot(k2,P2)*dot(k2,q1)**2*dot(P2,q1)*dot(q2,q2)**2 + 
     -       F2F2AA*dot(P2,q1)**2*dot(q2,q2)**2 + 
     -       4*F2F2CC*dot(k2,q1)**2*dot(P2,q1)**2*dot(q2,q2)**2 - 
     -       8*F2F2CC*dot(k1,q2)*dot(k2,q1)**2*dot(P2,q2)*
     -        (dot(k2,q2)*dot(P2,q2) + dot(P2,q2)*dot(q1,q2) + 
     -          (dot(k1,P2) - dot(k2,P2) - dot(P2,q1))*dot(q2,q2)) + 
     -       4*F2F2BB*dot(k1,q1)**2*(dot(k1,q2)*dot(P2,q2) - dot(k2,q2)*dot(P2,q2) + 
     -           dot(P2,q2)*dot(q1,q2) - dot(k1,P2)*dot(q2,q2) + dot(k2,P2)*dot(q2,q2) - 
     -           dot(P2,q1)*dot(q2,q2))**2) - 
     -    2*dot(P1,q1)*dot(q1,q1)*(dot(P1,q1)*dot(P2,q2)**3*dot(q1,q1)**2*
     -        (4*F2F2CC*dot(k1,q2)**2*dot(k2,P1)*dot(k2,q1) + 
     -          4*F2F2CC*dot(k2,P1)*dot(k2,q1)*dot(k2,q2)**2 + 
     -          8*F2F2CC*dot(k2,P1)*dot(k2,q1)*dot(k2,q2)*dot(q1,q2) + 
     -          F2F2AA*dot(P1,q2)*dot(q1,q2) + 
     -          4*F2F2CC*dot(k2,P1)*dot(k2,q1)*dot(q1,q2)**2 + 
     -          4*F2F2BB*dot(k1,P1)*dot(k1,q1)*(dot(k1,q2) - dot(k2,q2) + dot(q1,q2))**2 - 
     -          8*F2F2CC*dot(k1,q2)*dot(k2,P1)*dot(k2,q1)*(dot(k2,q2) + dot(q1,q2)))*
     -        dot(q2,q2)**2 + dot(P1,q1)*
     -        (4*F2F2CC*dot(k1,P2)**2*dot(k2,P1)*dot(k2,q1) + 
     -          4*F2F2CC*dot(k2,P1)*dot(k2,P2)**2*dot(k2,q1) + 
     -          8*F2F2CC*dot(k2,P1)*dot(k2,P2)*dot(k2,q1)*dot(P2,q1) + 
     -          F2F2AA*dot(P1,P2)*dot(P2,q1) + 
     -          4*F2F2CC*dot(k2,P1)*dot(k2,q1)*dot(P2,q1)**2 + 
     -          4*F2F2BB*dot(k1,P1)*dot(k1,q1)*(dot(k1,P2) - dot(k2,P2) + dot(P2,q1))**2 - 
     -          8*F2F2CC*dot(k1,P2)*dot(k2,P1)*dot(k2,q1)*(dot(k2,P2) + dot(P2,q1)))*
     -        dot(P2,q2)*dot(q1,q1)**2*dot(q2,q2)**4 + 
     -       dot(P2,q2)*dot(q2,q2)*(-(dot(P1,q1)*dot(P2,q2)*dot(q1,q1)**2*
     -             (-8*F2F2CC*dot(k1,q2)*dot(k2,P1)*dot(k2,P2)*dot(k2,q1) + 
     -               8*F2F2CC*dot(k2,P1)*dot(k2,P2)*dot(k2,q1)*dot(k2,q2) - 
     -               8*F2F2CC*dot(k1,q2)*dot(k2,P1)*dot(k2,q1)*dot(P2,q1) + 
     -               8*F2F2CC*dot(k2,P1)*dot(k2,q1)*dot(k2,q2)*dot(P2,q1) + 
     -               F2F2AA*dot(P1,q2)*dot(P2,q1) + 
     -               8*F2F2CC*dot(k1,P2)*dot(k2,P1)*dot(k2,q1)*
     -                (dot(k1,q2) - dot(k2,q2) - dot(q1,q2)) + 
     -               8*F2F2CC*dot(k2,P1)*dot(k2,P2)*dot(k2,q1)*dot(q1,q2) + 
     -               F2F2AA*dot(P1,P2)*dot(q1,q2) + 
     -               8*F2F2CC*dot(k2,P1)*dot(k2,q1)*dot(P2,q1)*dot(q1,q2) + 
     -               8*F2F2BB*dot(k1,P1)*dot(k1,q1)*(dot(k1,P2) - dot(k2,P2) + dot(P2,q1))*
     -                (dot(k1,q2) - dot(k2,q2) + dot(q1,q2)))*dot(q2,q2)**2) + 
     -          2*dot(q1,q1)*dot(q2,q2)*
     -           (dot(k1,q1)*dot(P1,q1)*
     -              (dot(k1,q2)*dot(P2,q2) - dot(k2,q2)*dot(P2,q2) + 
     -                dot(P2,q2)*dot(q1,q2) - dot(k1,P2)*dot(q2,q2) + 
     -                dot(k2,P2)*dot(q2,q2) - dot(P2,q1)*dot(q2,q2))*
     -              (2*F2F2BC*dot(k2,P1)*dot(k2,q2)*dot(P2,q2)*dot(q1,q1) - 
     -                F2F2AB*dot(P1,q2)*dot(P2,q2)*dot(q1,q1) + 
     -                2*F2F2BC*dot(k1,q2)*dot(P2,q2)*
     -                 (dot(k2,q1)*dot(P1,q1) - dot(k2,P1)*dot(q1,q1)) + 
     -                F2F2AB*dot(P1,q1)*dot(P2,q2)*dot(q1,q2) + 
     -                2*F2F2BC*dot(k2,P1)*dot(P2,q2)*dot(q1,q1)*dot(q1,q2) - 
     -                F2F2AB*dot(P1,q1)*dot(P2,q1)*dot(q2,q2) + 
     -                2*F2F2BC*dot(k1,P2)*dot(k2,P1)*dot(q1,q1)*dot(q2,q2) - 
     -                2*F2F2BC*dot(k2,P1)*dot(k2,P2)*dot(q1,q1)*dot(q2,q2) + 
     -                F2F2AB*dot(P1,P2)*dot(q1,q1)*dot(q2,q2) - 
     -                2*F2F2BC*dot(k2,P1)*dot(P2,q1)*dot(q1,q1)*dot(q2,q2) - 
     -                2*F2F2BC*dot(k2,q1)*dot(P1,q1)*
     -                 (dot(k2,q2)*dot(P2,q2) + dot(P2,q2)*dot(q1,q2) + 
     -                   (dot(k1,P2) - dot(k2,P2) - dot(P2,q1))*dot(q2,q2))) + 
     -             dot(k1,P1)*dot(q1,q1)*
     -              (dot(k1,q2)*dot(P2,q2) - dot(k2,q2)*dot(P2,q2) + 
     -                dot(P2,q2)*dot(q1,q2) - dot(k1,P2)*dot(q2,q2) + 
     -                dot(k2,P2)*dot(q2,q2) - dot(P2,q1)*dot(q2,q2))*
     -              (-2*F2F2BC*dot(k2,P1)*dot(k2,q2)*dot(P2,q2)*dot(q1,q1) + 
     -                F2F2AB*dot(P1,q2)*dot(P2,q2)*dot(q1,q1) + 
     -                2*F2F2BC*dot(k1,q2)*dot(P2,q2)*
     -                 (-(dot(k2,q1)*dot(P1,q1)) + dot(k2,P1)*dot(q1,q1)) - 
     -                F2F2AB*dot(P1,q1)*dot(P2,q2)*dot(q1,q2) - 
     -                2*F2F2BC*dot(k2,P1)*dot(P2,q2)*dot(q1,q1)*dot(q1,q2) + 
     -                F2F2AB*dot(P1,q1)*dot(P2,q1)*dot(q2,q2) - 
     -                2*F2F2BC*dot(k1,P2)*dot(k2,P1)*dot(q1,q1)*dot(q2,q2) + 
     -                2*F2F2BC*dot(k2,P1)*dot(k2,P2)*dot(q1,q1)*dot(q2,q2) - 
     -                F2F2AB*dot(P1,P2)*dot(q1,q1)*dot(q2,q2) + 
     -                2*F2F2BC*dot(k2,P1)*dot(P2,q1)*dot(q1,q1)*dot(q2,q2) + 
     -                2*F2F2BC*dot(k2,q1)*dot(P1,q1)*
     -                 (dot(k2,q2)*dot(P2,q2) + dot(P2,q2)*dot(q1,q2) + 
     -                   (dot(k1,P2) - dot(k2,P2) - dot(P2,q1))*dot(q2,q2))) - 
     -             F2F2AC*(dot(k2,q1)*dot(P1,q1) - dot(k2,P1)*dot(q1,q1))*
     -              (dot(k1,q2)*dot(P2,q2) - dot(k2,q2)*dot(P2,q2) - 
     -                dot(P2,q2)*dot(q1,q2) - dot(k1,P2)*dot(q2,q2) + 
     -                dot(k2,P2)*dot(q2,q2) + dot(P2,q1)*dot(q2,q2))*
     -              (-(dot(P1,q2)*dot(P2,q2)*dot(q1,q1)) + 
     -                dot(P1,P2)*dot(q1,q1)*dot(q2,q2) + 
     -                dot(P1,q1)*(dot(P2,q2)*dot(q1,q2) - dot(P2,q1)*dot(q2,q2)))))))/
     -     (dot(P1,q1)**2*dot(P2,q2)**2*dot(q1,q1)**4*dot(q2,q2)**4)

      end function

      
!----------------------------------------------------------------------
!     F3F3: implemented using 1401.7754
      function F3F3(F3F3AA, F3F3AB, F3F3AC, F3F3BB, F3F3BC, F3F3CC,
     &     q1in, q2in, P1in, P2in, k1in, k2in)
      real*8 F3F3AA, F3F3AB, F3F3AC, F3F3BB, F3F3BC, F3F3CC
      real*8 q1in(0:3),q2in(0:3),P1in(0:3),P2in(0:3),k1in(0:3),k2in(0:3)
      real*8 q1(1:4),q2(1:4),P1(1:4),P2(1:4),k1(1:4),k2(1:4)
      real*8 F3F3
      real*8 S, MH2,q1q1,P1q1,P1q2,P2q1,P2q2,k1k2,q1k1,q1k2,q2k1,q2k2
      real*8 q1q2, Q1sq,Q2sq,P1k1,P1k2,P2k1,P2k2
      q1(1:4) = q1in(0:3)
      q2(1:4) = q2in(0:3)
      P1(1:4) = P1in(0:3)
      P2(1:4) = P2in(0:3)
      k1(1:4) = k1in(0:3)
      k2(1:4) = k2in(0:3)
      S    = dot(P1+P2,P1+P2)
      MH2  = dot(k1,k1)
      q1q2 = dot(q1,q2)
      P1q1 = dot(P1,q1)
      P1q2 = dot(P1,q2)
      P2q1 = dot(P2,q1)
      P2q2 = dot(P2,q2)
      P1k1 = dot(P1,k1)
      P1k2 = dot(P1,k2)
      P2k1 = dot(P2,k1)
      P2k2 = dot(P2,k2)
      k1k2 = dot(k1,k2)
      q1k1 = dot(q1,k1)
      q1k2 = dot(q1,k2)
      q2k1 = dot(q2,k1)
      q2k2 = dot(q2,k2)
      Q1sq = -dot(q1,q1)
      Q2sq = -dot(q2,q2)
      
      F3F3 =(F3F3AA*((q1q2)*S-2*(P1q2)*
     -      (P2q1))-2*F3F3AB*(S*(-(k1k2)*(q1q2)               
     -      +MH2* (q1q2)+(q1k1)* (q1q2)-(q1k1)* (q2k1)+(q1k2)*
     -      (q2k1)+(q2k1)* Q1sq)                                          
     -      -2 *(MH2 *(P1q2)*(P2q1)-(k1k2)* (P1q2)* (P2q1)
     -      +(P1k1)* (P2k1)* (q1q2)                                     
     -      -(P1k1)* (P2q1)* (q2k1)-(P1k2)* (P2k1)*(q1q2)
     -      +(P1k2)* (P2q1)* (q2k1)                                     
     -      +(P1q1)*(P2k1)* (q1q2)-(P1q1)*(P2q1)*(q2k1)
     -      -(P1q2)*(P2k1)*(q1k1)                                     
     -      +(P1q2)*(P2k1)*(q1k2)+(P1q2)*(P2k1)* Q1sq
     -      +(P1q2)*(P2q1)*(q1k1)))                           
     -      -2 *F3F3AC* (S* (MH2* (q1q2)-(k1k2)*(q1q2)+(q1k1)*(q2k2)
     -      +(q1k2)*(q1q2)                                                    
     -      -(q1k2)*(q2k2)+(q2k2) *Q1sq)
     -      -2 *(MH2* (P1q2)*(P2q1)                                       
     -      -(k1k2)*(P1q2)*(P2q1)-(P1k1)*(P2k2)*(q1q2)
     -      +(P1k1)*(P2q1)*(q2k2)                                     
     -      +(P1k2)*(P2k2)*(q1q2) -(P1k2)*(P2q1)*(q2k2)
     -      +(P1q1)*(P2k2)*(q1q2)                                     
     -      -(P1q1)*(P2q1)*(q2k2)+(P1q2)*(P2k2)*(q1k1)
     -      -(P1q2)*(P2k2)*(q1k2)                                     
     -      +(P1q2)*(P2k2) *Q1sq+(P1q2)*(P2q1)*(q1k2)))
     -      +8 *F3F3BC *(2 *(k1k2)*                                                     
     -      ((P1q1)*(q1q2)*(P2k1+P2k2)
     -      -(P1q1)*(P2q1)*(q2k1+q2k2)                        
     -      +Q1sq *(P1q2)*(P2k1+P2k2)+(P1q2)*(P2q1)*
     -      (q1k1+q1k2))                                                  
     -      -(k1k2) *S* ((q1q2)*(q1k1+q1k2)+Q1sq*
     -      (q2k1+q2k2))                                                  
     -      -2* MH2* ((P1q1)*(q1q2)*(P2k1+P2k2)
     -      -(P1q1)*(P2q1)*(q2k1+q2k2)                        
     -      +(P1q2) *Q1sq* (P2k1+P2k2)+(P1q2) * (P2q1)*
     -      (q1k1+q1k2))                                                  
     -      +MH2 * S *((q1q2)*(q1k1+q1k2)+Q1sq*(q2k1+q2k2))    
     -      +2 *(Q1sq *(P1k1-P1k2)*(P2k2)*(q2k1)
     -      -Q1sq *(P1k1-P1k2)*(P2k1)*(q2k2)                  
     -      -(P1k1)*(P2k1)*(q1k2)*(q1q2)
     -      +(P1k1)*(P2k2)*(q1k1)*(q1q2)                      
     -      -(P1k1)*(P2q1)*(q1k1)*(q2k2)
     -      +(P1k1)*(P2q1)*(q1k2)*(q2k1)                      
     -      +(P1k2)*(P2k1)*(q1k2)*(q1q2)
     -      -(P1k2)*(P2k2)*(q1k1)*(q1q2)                      
     -      +(P1k2)*(P2q1)*(q1k1)*(q2k2)
     -      -(P1k2)*(P2q1)*(q1k2)*(q2k1)                      
     -      -(P1q1)*(P2k1)*(q1k1)*(q2k2)
     -      +(P1q1)*(P2k1)*(q1k2)*(q2k2)                      
     -      +(P1q1)*(P2k2)*(q1k1)*(q2k1)
     -      -(P1q1)*(P2k2)*(q1k2)*(q2k1)                      
     -      -(P1q2)*(P2k2)*(q1k1)**2
     -      +(P1q2)*(q1k2)*(P2k2)*(q1k1)                        
     -      +(P1q2)*(q1k1)*(P2k1)*(q1k2)
     -      -(P1q2)*(P2k1)*(q1k2)**2)                   
     -      +S *((q1k1)-(q1k2))*((q1k1)*(q2k2)
     -      -(q1k2) * (q2k1))))
     -      /(4*(P1q1)*(P2q2))
      
      end function
c......................................................................
      end module ME_expressions
