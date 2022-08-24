      module nonfact_expressions
      implicit none

      public
      contains

      function tri1110 (MV2, q1sq, q2sq, s, l2)
      real*8 MV2, q1sq, q2sq, s,l2
      real*8 tri1110
      ! extracted from logbook-dihiggs/2020-02-06-non-fact/fortran_format.nb
      !lambda^2 piece: -(Log(l2/MV_sq)/((MV_sq + q1sq)*(MV_sq + q2sq))) 
      tri1110 = - (((MV2*(q1sq - q2sq + s) + q1sq*(-q1sq + q2sq + s))
     $     *Log(MV2/(MV2 + q1sq)))/(MV2 + q1sq) +((q2sq*(q1sq - q2sq
     $     + s) + MV2*(-q1sq + q2sq + s))*Log(MV2/(MV2 + q2sq)))
     $     /(MV2 + q2sq) -(s*(2*MV2 - q1sq - q2sq + s)*Log((1 - Sqrt(s
     $     /(4*MV2 + s)))/(1 + Sqrt(s/(4*MV2 + s)))))/Sqrt(s*(4*MV2 +
     $     s)))/(MV2**2*s + q1sq*q2sq*s + MV2*(-q1sq**2 + q2sq*(-q2sq
     $     + s) + q1sq*(2*q2sq + s)))
!     To conform with Kirill's triangle and lambda piece
      tri1110 = tri1110 - log(l2/MV2)/((MV2 + q1sq)*(MV2 + q2sq))
!      print*, 'tri1110', tri1110, MV2, q1sq, q2sq, s, l2
      end function

      function tri0111 (MV2, qH1sq, qH2sq, s)
      real*8 MV2, qH1sq, qH2sq, s
      real*8 tri0111
      ! extracted from logbook-dihiggs/2020-02-06-non-fact/fortran_format.nb
      tri0111 = -((qH1sq*(qH1sq - qH2sq - s)*Log((1 - Sqrt(qH1sq/(4*MV2
     $     + qH1sq)))/(1 + Sqrt(qH1sq/(4*MV2 + qH1sq)))))/Sqrt(qH1sq*(4
     $     *MV2 + qH1sq)) + (qH2sq*(-qH1sq + qH2sq - s)*Log((1 -
     $     Sqrt(qH2sq/(4*MV2 + qH2sq)))/(1 + Sqrt(qH2sq/(4*MV2 +
     $     qH2sq)))))/Sqrt(qH2sq*(4*MV2 + qH2sq)) +(s*(-qH1sq - qH2sq +
     $     s)*Log((1 - Sqrt(s/(4*MV2 + s)))/(1 + Sqrt(s/(4*MV2 +
     $     s)))))/Sqrt(s*(4*MV2 + s)))/(-(qH1sq*qH2sq*s) + MV2*(qH1sq
     $     **2 + (qH2sq - s)**2 - 2*qH1sq*(qH2sq + s)))
!     To conform with Kirill's triangle
      tri0111 = tri0111
!      print*, 'tri0111', tri0111, MV2, qH1sq, qH2sq, s
      end function


      function tri1011 (MV2, q2sq, qH2sq, t, l2)
      real*8 MV2, q2sq, qH2sq, t, l2
      real*8 tri1011
      ! extracted from logbook-dihiggs/2020-02-06-non-fact/fortran_format.nb
      ! lambda^2 piece: -(Log(l2/MV_sq)/((MV_sq + q2sq)*(MV_sq + t)))
      tri1011 =   -(((MV2*(q2sq + qH2sq - t) + q2sq*(-q2sq + qH2sq + t))
     $     *Log(MV2/(MV2 + q2sq)))/(MV2 + q2sq) -(qH2sq*(2*MV2 -
     $     q2sq + qH2sq - t)*Log((1 - Sqrt(qH2sq/(4*MV2 + qH2sq)))/(1 +
     $     Sqrt(qH2sq/(4*MV2 + qH2sq)))))/Sqrt(qH2sq*(4*MV2 + qH2sq))
     $     + (((q2sq + qH2sq - t)*t + MV2*(-q2sq + qH2sq + t))*Log(MV2
     $     /(MV2 + t)))/(MV2 + t))/(MV2**2*qH2sq + q2sq*qH2sq*t +
     $     MV2*(-q2sq**2 + (qH2sq - t)*t + q2sq*(qH2sq + 2*t)))
!     To conform with Kirill's triangle and lambda piece
      tri1011 = tri1011 - log(l2/MV2)/((MV2 + q2sq)*(MV2 + t))
!      print*, 'tri1011', tri1011, MV2, q2sq, qH2sq, t, l2
      end function
      
      function tri1101 (MV2, q1sq, qH1sq, t, l2)
      real*8 MV2, q1sq, qH1sq, t, l2
      real*8 tri1101
      ! extracted from logbook-dihiggs/2020-02-06-non-fact/fortran_format.nb
      ! lambda^2 piece: -(Log(l2/MV2)/((MV2 + q1sq)*(MV2 + t))) 
      tri1101 = - (((MV2*(q1sq + qH1sq - t) + q1sq*(-q1sq + qH1sq + t))
     $     *Log(MV2/(MV2 + q1sq)))/(MV2 + q1sq) -(qH1sq*(2*MV2 -
     $     q1sq + qH1sq - t)*Log((1 - Sqrt(qH1sq/(4*MV2 + qH1sq)))/(1 +
     $     Sqrt(qH1sq/(4*MV2 + qH1sq)))))/Sqrt(qH1sq*(4*MV2 + qH1sq))
     $     + (((q1sq + qH1sq - t)*t + MV2*(-q1sq + qH1sq + t))*Log(MV2
     $     /(MV2 + t)))/(MV2 + t))/(MV2**2*qH1sq + q1sq*qH1sq*t +
     $     MV2*(-q1sq**2 + (qH1sq - t)*t + q1sq*(qH1sq + 2*t)))     
!     To conform with Kirill's triangle and lambda piece
      tri1101 = tri1101 - log(l2/MV2)/((MV2 + q1sq)*(MV2 + t))
!      print*, 'tri1101', tri1101, MV2, q1sq, qH1sq, t, l2
      end function

!!! New implementation using explicit coefficients for the box

      function box_1loop (MV, q1sq, q2sq, qH1sq, qH2sq, s, t, l2)
      real*8 MV, q1sq, q2sq, qH1sq, qH2sq, s, t, l2
      real*8 box_1loop
      
      box_1loop =
     $       C1110(MV**2, q1sq, q2sq, qH1sq, qH2sq, s, t) *
     $                         tri1110(MV**2, q1sq, q2sq, s, l2)
     $     + C1101(MV**2, q1sq, q2sq, qH1sq, qH2sq, s, t) *
     $                         tri1101(MV**2, q1sq, qH1sq, t, l2)
     $     + C0111(MV**2, q1sq, q2sq, qH1sq, qH2sq, s, t) *
     $                         tri0111(MV**2, qH1sq, qH2sq, s)
     $     + C1011(MV**2, q1sq, q2sq, qH1sq, qH2sq, s, t) *
     $                         tri1011(MV**2, q2sq, qH2sq, t, l2)
      end function

      function C1011(MV2, q1sq, q2sq, qH1sq, qH2sq, s, t)
      real*8 MV2, q1sq, q2sq, qH1sq, qH2sq, s, t
      real*8 C1011

      C1011 = (q2sq**2*qH1sq + (qH2sq - t)*(q1sq*qH2sq - s*t) - MV2*(-2
     $     *q1sq*qH2sq + qH1sq*qH2sq - qH2sq**2 + q2sq*(qH1sq + qH2sq -
     $     s) + qH2sq*s - qH1sq*t + qH2sq*t + s*t) - q2sq*(q1sq*qH2sq +
     $     (-2*qH2sq + s)*t + qH1sq*(qH2sq + t)))/ (q2sq**2*qH1sq**2 +
     $     MV2**2*(qH1sq**2 + (qH2sq - s)**2 - 2*qH1sq*(qH2sq + s)) +
     $     (q1sq*qH2sq - s*t)**2 - 2*q2sq*qH1sq*(q1sq*qH2sq + s*t) + 2
     $     *MV2*(2*q2sq**2*qH1sq + 2*q1sq**2*qH2sq + s*t*(-qH1sq -
     $     qH2sq + s + 2*t) - q1sq*(qH1sq*qH2sq - qH2sq**2 + 2*q2sq
     $     *(qH1sq + qH2sq - s) + qH2sq*s - 2*qH1sq*t + 2*qH2sq*t + 2*s
     $     *t) + q2sq*(qH1sq**2 + 2*(qH2sq - s)*t - qH1sq*(qH2sq + s + 2
     $     *t))))
      end function

      function C0111(MV2, q1sq, q2sq, qH1sq, qH2sq, s, t)
      real*8 MV2, q1sq, q2sq, qH1sq, qH2sq, s, t
      real*8 C0111

      C0111 = -((q1sq*qH1sq*qH2sq - q1sq*qH2sq**2 + q1sq*qH2sq*s - 2
     $     *qH1sq*qH2sq*s + q2sq*qH1sq*(-qH1sq + qH2sq + s) +MV2*(qH1sq
     $     **2 + (qH2sq - s)**2 - 2*qH1sq*(qH2sq + s)) + qH1sq*s*t +
     $     qH2sq*s*t - s**2*t)/(q2sq**2*qH1sq**2 + MV2**2*(qH1sq**2 +
     $     (qH2sq - s)**2 - 2*qH1sq*(qH2sq + s)) + (q1sq*qH2sq - s*t)**2
     $     -2*q2sq*qH1sq*(q1sq*qH2sq + s*t) + 2*MV2*(2*q2sq**2*qH1sq +
     $     2*q1sq**2*qH2sq + s*t*(-qH1sq - qH2sq + s + 2*t) -q1sq*(qH1sq
     $     *qH2sq - qH2sq**2 + 2*q2sq*(qH1sq + qH2sq - s) + qH2sq*s - 2
     $     *qH1sq*t + 2*qH2sq*t + 2*s*t) +q2sq*(qH1sq**2 + 2*(qH2sq - s)
     $     *t - qH1sq*(qH2sq + s + 2*t)))))
      end function


      function C1101(MV2, q1sq, q2sq, qH1sq, qH2sq, s, t)
      real*8 MV2, q1sq, q2sq, qH1sq, qH2sq, s, t
      real*8 C1101

      C1101 = (q1sq**2*qH2sq + (qH1sq - t)*(q2sq*qH1sq - s*t) - MV2*(-2
     $     *q2sq*qH1sq - qH1sq**2 + qH1sq*qH2sq + q1sq*(qH1sq + qH2sq -
     $     s) + qH1sq*s + qH1sq*t - qH2sq*t + s*t) - q1sq*(q2sq*qH1sq +
     $     qH1sq*(qH2sq - 2*t) + (qH2sq + s)*t))/ (q2sq**2*qH1sq**2 +
     $     MV2**2*(qH1sq**2 + (qH2sq - s)**2 - 2*qH1sq*(qH2sq + s)) +
     $     (q1sq*qH2sq - s*t)**2 - 2*q2sq*qH1sq*(q1sq*qH2sq + s*t) + 2
     $     *MV2*(2*q2sq**2*qH1sq + 2*q1sq**2*qH2sq + s*t*(-qH1sq -
     $     qH2sq + s + 2*t) - q1sq*(qH1sq*qH2sq - qH2sq**2 + 2*q2sq
     $     *(qH1sq + qH2sq - s) + qH2sq*s - 2*qH1sq*t + 2*qH2sq*t + 2*s
     $     *t) + q2sq*(qH1sq**2 + 2*(qH2sq - s)*t - qH1sq*(qH2sq + s + 2
     $     *t))))
      end function


      function C1110(MV2, q1sq, q2sq, qH1sq, qH2sq, s, t)
      real*8 MV2, q1sq, q2sq, qH1sq, qH2sq, s, t
      real*8 C1110

      C1110 = (q1sq**2*qH2sq - MV2*(q2sq*(qH1sq - qH2sq + s) + q1sq*(
     $     -qH1sq + qH2sq + s) + s*(qH1sq + qH2sq - s - 2*t)) +(q2sq -
     $     s)*(q2sq*qH1sq - s*t) - q1sq*(q2sq*(qH1sq + qH2sq - 2*s) + s
     $     *(qH2sq + t)))/(q2sq**2*qH1sq**2 + MV2**2*(qH1sq**2 + (qH2sq
     $     - s)**2 - 2*qH1sq*(qH2sq + s)) + (q1sq*qH2sq - s*t)**2 -2
     $     *q2sq*qH1sq*(q1sq*qH2sq + s*t) + 2*MV2*(2*q2sq**2*qH1sq + 2
     $     *q1sq**2*qH2sq + s*t*(-qH1sq - qH2sq + s + 2*t) -q1sq*(qH1sq
     $     *qH2sq - qH2sq**2 + 2*q2sq*(qH1sq + qH2sq - s) + qH2sq*s - 2
     $     *qH1sq*t + 2*qH2sq*t + 2*s*t) +q2sq*(qH1sq**2 + 2*(qH2sq - s)
     $     *t - qH1sq*(qH2sq + s + 2*t))))
      end function

      function box_1loop_detailed (MV2, q1sq, q2sq, qH1sq, qH2sq, s, t, l2)
      real*8 MV2, q1sq, q2sq, qH1sq, qH2sq, s, t, l2
      real*8 box_1loop_detailed
      ! extracted from logbook-dihiggs/2020-02-06-non-fact/fortran_format.nb
      box_1loop_detailed =  -(((q1sq*qH1sq*qH2sq - q1sq*qH2sq**2 + q1sq*qH2sq*s - 
     -         2*qH1sq*qH2sq*s + q2sq*qH1sq*(-qH1sq + qH2sq + s) + 
     -         MV2*(qH1sq**2 + (qH2sq - s)**2 - 
     -            2*qH1sq*(qH2sq + s)) + qH1sq*s*t + qH2sq*s*t - 
     -         s**2*t)*tri0111 (MV2, qH1sq, qH2sq, s))/
     -     (q2sq**2*qH1sq**2 + 
     -       MV2**2*(qH1sq**2 + (qH2sq - s)**2 - 
     -          2*qH1sq*(qH2sq + s)) + (q1sq*qH2sq - s*t)**2 - 
     -       2*q2sq*qH1sq*(q1sq*qH2sq + s*t) + 
     -       2*MV2*(2*q2sq**2*qH1sq + 2*q1sq**2*qH2sq + 
     -          s*t*(-qH1sq - qH2sq + s + 2*t) - 
     -          q1sq*(qH1sq*qH2sq - qH2sq**2 + 
     -             2*q2sq*(qH1sq + qH2sq - s) + qH2sq*s - 
     -             2*qH1sq*t + 2*qH2sq*t + 2*s*t) + 
     -          q2sq*(qH1sq**2 + 2*(qH2sq - s)*t - 
     -             qH1sq*(qH2sq + s + 2*t))))) + 
     -  ((q2sq**2*qH1sq + (qH2sq - t)*(q1sq*qH2sq - s*t) - 
     -       MV2*(-2*q1sq*qH2sq + qH1sq*qH2sq - qH2sq**2 + 
     -          q2sq*(qH1sq + qH2sq - s) + qH2sq*s - qH1sq*t + 
     -          qH2sq*t + s*t) - 
     -       q2sq*(q1sq*qH2sq + (-2*qH2sq + s)*t + 
     -          qH1sq*(qH2sq + t)))*tri1011 (MV2, q2sq, qH2sq, t, l2))/
     -   (q2sq**2*qH1sq**2 + 
     -     MV2**2*(qH1sq**2 + (qH2sq - s)**2 - 
     -        2*qH1sq*(qH2sq + s)) + (q1sq*qH2sq - s*t)**2 - 
     -     2*q2sq*qH1sq*(q1sq*qH2sq + s*t) + 
     -     2*MV2*(2*q2sq**2*qH1sq + 2*q1sq**2*qH2sq + 
     -        s*t*(-qH1sq - qH2sq + s + 2*t) - 
     -        q1sq*(qH1sq*qH2sq - qH2sq**2 + 
     -           2*q2sq*(qH1sq + qH2sq - s) + qH2sq*s - 2*qH1sq*t + 
     -           2*qH2sq*t + 2*s*t) + 
     -        q2sq*(qH1sq**2 + 2*(qH2sq - s)*t - 
     -           qH1sq*(qH2sq + s + 2*t)))) + 
     -  ((q1sq**2*qH2sq + (qH1sq - t)*(q2sq*qH1sq - s*t) - 
     -       MV2*(-2*q2sq*qH1sq - qH1sq**2 + qH1sq*qH2sq + 
     -          q1sq*(qH1sq + qH2sq - s) + qH1sq*s + qH1sq*t - 
     -          qH2sq*t + s*t) - 
     -       q1sq*(q2sq*qH1sq + qH1sq*(qH2sq - 2*t) + (qH2sq + s)*t)
     -       )*tri1101 (MV2, q1sq, qH1sq, t, l2))/
     -   (q2sq**2*qH1sq**2 + 
     -     MV2**2*(qH1sq**2 + (qH2sq - s)**2 - 
     -        2*qH1sq*(qH2sq + s)) + (q1sq*qH2sq - s*t)**2 - 
     -     2*q2sq*qH1sq*(q1sq*qH2sq + s*t) + 
     -     2*MV2*(2*q2sq**2*qH1sq + 2*q1sq**2*qH2sq + 
     -        s*t*(-qH1sq - qH2sq + s + 2*t) - 
     -        q1sq*(qH1sq*qH2sq - qH2sq**2 + 
     -           2*q2sq*(qH1sq + qH2sq - s) + qH2sq*s - 2*qH1sq*t + 
     -           2*qH2sq*t + 2*s*t) + 
     -        q2sq*(qH1sq**2 + 2*(qH2sq - s)*t - 
     -           qH1sq*(qH2sq + s + 2*t)))) + 
     -  ((q1sq**2*qH2sq - MV2*
     -        (q2sq*(qH1sq - qH2sq + s) + 
     -          q1sq*(-qH1sq + qH2sq + s) + 
     -          s*(qH1sq + qH2sq - s - 2*t)) + 
     -       (q2sq - s)*(q2sq*qH1sq - s*t) - 
     -       q1sq*(q2sq*(qH1sq + qH2sq - 2*s) + s*(qH2sq + t)))*
     -     tri1110 (MV2, q1sq, q2sq, s, l2))/
     -   (q2sq**2*qH1sq**2 + 
     -     MV2**2*(qH1sq**2 + (qH2sq - s)**2 - 
     -        2*qH1sq*(qH2sq + s)) + (q1sq*qH2sq - s*t)**2 - 
     -     2*q2sq*qH1sq*(q1sq*qH2sq + s*t) + 
     -     2*MV2*(2*q2sq**2*qH1sq + 2*q1sq**2*qH2sq + 
     -        s*t*(-qH1sq - qH2sq + s + 2*t) - 
     -        q1sq*(qH1sq*qH2sq - qH2sq**2 + 
     -           2*q2sq*(qH1sq + qH2sq - s) + qH2sq*s - 2*qH1sq*t + 
     -           2*qH2sq*t + 2*s*t) + 
     -        q2sq*(qH1sq**2 + 2*(qH2sq - s)*t - 
     -     qH1sq*(qH2sq + s + 2*t))))

      end function

      ! No lambda dependence!!!
      function box_1loop_expanded (MV2, q1sq, q2sq, qH1sq, qH2sq, s, t)
      real*8 MV2, q1sq, q2sq, qH1sq, qH2sq, s, t
      real*8 box_1loop_expanded

      box_1loop_expanded = -((((-(q1sq**2*qH2sq) + q1sq*q2sq*(qH1sq +
     $     qH2sq - 2*s) + MV2*(q2sq*(qH1sq - qH2sq + s) + q1sq*(-qH1sq
     $     + qH2sq + s) + s*(qH1sq + qH2sq - s - 2*t)) + q1sq*s*(qH2sq +
     $     t) - (q2sq - s)*(q2sq*qH1sq - s*t))*(((MV2*(q1sq - q2sq + s)
     $     + q1sq*(-q1sq + q2sq + s))*Log(MV2/(MV2 + q1sq)))/(MV2 +
     $     q1sq) + ((q2sq*(q1sq - q2sq + s) + MV2*(-q1sq + q2sq + s))
     $     *Log(MV2/(MV2 + q2sq)))/(MV2 + q2sq) -(s*(2*MV2 - q1sq -
     $     q2sq + s)*Log((1 - Sqrt(s/(4*MV2 + s)))/(1 + Sqrt(s/(4*MV2
     $     + s)))))/Sqrt(s*(4*MV2 + s))))/(MV2**2*s + q1sq*q2sq*s +
     $     MV2*(-q1sq**2 + q2sq*(-q2sq + s) + q1sq*(2*q2sq + s)))+
     $     ((q1sq*qH1sq*qH2sq - q1sq*qH2sq**2 + q1sq*qH2sq*s - 2*qH1sq
     $     *qH2sq*s + q2sq*qH1sq*(-qH1sq + qH2sq + s) + MV2*(qH1sq**2 +
     $     (qH2sq - s)**2 - 2*qH1sq*(qH2sq + s)) + qH1sq*s*t + qH2sq*s*t
     $     - s**2*t)*(-(qH1sq*Sqrt(qH2sq*(4*MV2 + qH2sq))*Sqrt(s*(4
     $     *MV2 + s))*(-qH1sq + qH2sq + s)*Log((1 - Sqrt(qH1sq/(4*MV2
     $     + qH1sq)))/(1 + Sqrt(qH1sq/(4*MV2 + qH1sq))))) -Sqrt(qH1sq
     $     *(4*MV2 + qH1sq))*(qH2sq*Sqrt(s*(4*MV2 + s))*(qH1sq - qH2sq
     $     + s)*Log((1 - Sqrt(qH2sq/(4*MV2 + qH2sq)))/(1 + Sqrt(qH2sq
     $     /(4*MV2 + qH2sq)))) +Sqrt(qH2sq*(4*MV2 + qH2sq))*(qH1sq +
     $     qH2sq - s)*s*Log((1 - Sqrt(s/(4*MV2 + s)))/(1 + Sqrt(s/(4
     $     *MV2 + s)))))))/(Sqrt(qH1sq*(4*MV2 + qH1sq))*Sqrt(qH2sq*(4
     $     *MV2 + qH2sq))*Sqrt(s*(4*MV2 + s))*(-(qH1sq*qH2sq*s) + MV2
     $     *(qH1sq**2 + (qH2sq - s)**2 - 2*qH1sq*(qH2sq + s)))) +((
     $     -(q1sq**2*qH2sq) - (qH1sq - t)*(q2sq*qH1sq - s*t) + MV2*(-2
     $     *q2sq*qH1sq - qH1sq**2 + qH1sq*qH2sq + q1sq*(qH1sq + qH2sq -
     $     s) + qH1sq*s + qH1sq*t - qH2sq*t + s*t) +q1sq*(q2sq*qH1sq +
     $     qH1sq*(qH2sq - 2*t) + (qH2sq + s)*t))*(((MV2*(q1sq + qH1sq -
     $     t) + q1sq*(-q1sq + qH1sq + t))*Log(MV2/(MV2 + q1sq)))/(MV2
     $     + q1sq) -(qH1sq*(2*MV2 - q1sq + qH1sq - t)*Log((1 -
     $     Sqrt(qH1sq/(4*MV2 + qH1sq)))/(1 + Sqrt(qH1sq/(4*MV2 +
     $     qH1sq)))))/Sqrt(qH1sq*(4*MV2 + qH1sq)) +(((q1sq + qH1sq - t)
     $     *t + MV2*(-q1sq + qH1sq + t))*Log(MV2/(MV2 + t)))/(MV2 +
     $     t)))/(MV2**2*qH1sq + q1sq*qH1sq*t + MV2*(-q1sq**2 + (qH1sq
     $     - t)*t + q1sq*(qH1sq + 2*t))) +((-(q2sq**2*qH1sq) - (qH2sq -
     $     t)*(q1sq*qH2sq - s*t) + MV2*(-2*q1sq*qH2sq + qH1sq*qH2sq -
     $     qH2sq**2 + q2sq*(qH1sq + qH2sq - s) + qH2sq*s - qH1sq*t +
     $     qH2sq*t + s*t) +q2sq*(q1sq*qH2sq + (-2*qH2sq + s)*t + qH1sq
     $     *(qH2sq + t)))*(((MV2*(q2sq + qH2sq - t) + q2sq*(-q2sq +
     $     qH2sq + t))*Log(MV2/(MV2 + q2sq)))/(MV2 + q2sq) -(qH2sq*(2
     $     *MV2 - q2sq + qH2sq - t)*Log((1 - Sqrt(qH2sq/(4*MV2 +
     $     qH2sq)))/(1 + Sqrt(qH2sq/(4*MV2 + qH2sq)))))/Sqrt(qH2sq*(4
     $     *MV2 + qH2sq)) +(((q2sq + qH2sq - t)*t + MV2*(-q2sq + qH2sq
     $     + t))*Log(MV2/(MV2 + t)))/(MV2 + t)))/(MV2**2*qH2sq +
     $     q2sq*qH2sq*t + MV2*(-q2sq**2 + (qH2sq - t)*t + q2sq*(qH2sq +
     $     2*t))))/(q2sq**2*qH1sq**2 + MV2**2*(qH1sq**2 + (qH2sq - s)
     $     **2 - 2*qH1sq*(qH2sq + s)) + (q1sq*qH2sq - s*t)**2 - 2*q2sq
     $     *qH1sq*(q1sq*qH2sq + s*t) +2*MV2*(2*q2sq**2*qH1sq + 2*q1sq
     $     **2*qH2sq + s*t*(-qH1sq - qH2sq + s + 2*t) - q1sq*(qH1sq
     $     *qH2sq - qH2sq**2 + 2*q2sq*(qH1sq + qH2sq - s) + qH2sq*s - 2
     $     *qH1sq*t + 2*qH2sq*t + 2*s*t) +q2sq*(qH1sq**2 + 2*(qH2sq - s)
     $     *t - qH1sq*(qH2sq + s + 2*t)))))
      box_1loop_expanded = -box_1loop_expanded
      
      end function
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

      function r5(MV2, p1x, p3x, p3y, xi)
      real*8 MV2, p1x, p3x, p3y, xi
      complex*8 r5
      r5 = p1x*Cos(xi) + p3x*Cos(xi) + p3y*Sin(xi) - 
     -  (0,0.5)*Sqrt(4*
     -      (MV2 + p1x**2 + 2*p1x*p3x + p3x**2 + p3y**2)
     -       - 4*((p1x + p3x)*Cos(xi) + p3y*Sin(xi))**2)
      end function

      function r6(MV2, p1x, p3x, p3y, xi)
      real*8 MV2, p1x, p3x, p3y, xi
      complex*8 r6
      r6 = p1x*Cos(xi) + p3x*Cos(xi) + p3y*Sin(xi) + 
     -  (0,0.5)*Sqrt(4*
     -      (MV2 + p1x**2 + 2*p1x*p3x + p3x**2 + p3y**2)
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
      function b01(MV2, pi, p1x, p2x, p2y, p3x, p3y,xi)
      real*8 MV, MV2, pi, p1x, p2x, p2y, p3x, p3y, xi
      complex*8 r1v, r2v, r3v, r4v, r5v, r6v
      real*8 b01
      MV = sqrt(MV2)
      r1v=r1(MV2, p1x, xi)
      r2v=r2(MV2, p1x, xi)
      r3v=r3(MV2, p2x, p2y, xi)
      r4v=r4(MV2, p2x, p2y, xi)
      r5v=r5(MV2, p1x, p3x, p3y, xi)
      r6v=r6(MV2, p1x, p3x, p3y, xi)

      b01 = -(Log(-(r1v/Mv))/(Pi*r1v*(r1v - r2v)*(r1v - r3v)*(r1v - r4v)
     $     *(r1v - r5v)*(r1v - r6v))) +Log(-(r3v/Mv))/(Pi*(r2v - r3v)
     $     *r3v*(-r1v + r3v)*(r3v - r4v)*(r3v - r5v)*(r3v - r6v)) +Log(
     $     -(r5v/Mv))/(Pi*(r2v - r5v)*r5v*(-r1v + r5v)*(-r3v + r5v)*(
     $     -r4v + r5v)*(r5v - r6v))
      end function

      function b11(MV2, pi, p1x, p2x, p2y, p3x, p3y)
      real*8 MV2, pi, p1x, p2x, p2y, p3x, p3y
      real*8 b11
      b11 = -(1/((MV2 + p1x**2)*(MV2 + p2x**2 + p2y**2)*(MV2 + p1x**2 +
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

      function b012(MV2, pi, p1x, p2x, p2y, p3x, p3y)
      real*8 MV2, pi, p1x, p2x, p2y, p3x, p3y, xi
      real*8 b012
      b012 = (4*Pi**2)/(3.*(MV2 + p1x**2)*(MV2 + p2x**2 + p2y**2)*(MV2 +
     $     (p1x + p3x)**2 + p3y**2))
      end function

      function b022(MV2, pi, p1x, p2x, p2y, p3x, p3y, xi)
      real*8 MV, MV2, pi, p1x, p2x, p2y, p3x, p3y, xi
      complex*8 r1v, r2v, r3v, r4v, r5v, r6v
      complex*8 b022
      MV = sqrt(MV2)
      r1v=r1(MV2, p1x, xi)
      r2v=r2(MV2, p1x, xi)
      r3v=r3(MV2, p2x, p2y, xi)
      r4v=r4(MV2, p2x, p2y, xi)
      r5v=r5(MV2, p1x, p3x, p3y, xi)
      r6v=r6(MV2, p1x, p3x, p3y, xi)
      
      b022 =  (-2d0*Log(-(r1v/Mv))**2)/ (Pi*r1v*(r1v - r2v)*(r1v - r3v)
     $     *(r1v - r4v)*(r1v - r5v)*(r1v - r6v)) - (2d0*Log(-(r3v/Mv))
     $     **2) / (Pi*r3v*(-r1v + r3v)*(-r2v + r3v)*(r3v - r4v)*(r3v -
     $     r5v) *(r3v - r6v)) - (2d0*Log(-(r5v/Mv))**2)/ (Pi*r5v*(-r1v +
     $     r5v)*( -r2v + r5v)*(-r3v + r5v)*(-r4v + r5v)*(r5v - r6v))
      end function

      function b12(MV2, pi, p1x, p2x, p2y, p3x, p3y, xi)
      real*8 MV, MV2, pi, p1x, p2x, p2y, p3x, p3y, xi
      complex*8 r1v, r2v, r3v, r4v, r5v, r6v
      complex*8 b12
      MV = sqrt(MV2)
      r1v=r1(MV2, p1x, xi)
      r2v=r2(MV2, p1x, xi)
      r3v=r3(MV2, p2x, p2y, xi)
      r4v=r4(MV2, p2x, p2y, xi)
      r5v=r5(MV2, p1x, p3x, p3y, xi)
      r6v=r6(MV2, p1x, p3x, p3y, xi)
      b12 = (2d0*Log(-(r1v/Mv)))/(Pi*r1v*(r1v - r2v)*(r1v - r3v)*(r1v -
     $     r4v)*(r1v - r5v)*(r1v - r6v)) +(2d0*Log(-(r3v/Mv)))/(Pi*r3v*(
     $     -r1v + r3v)*(-r2v + r3v)*(r3v - r4v)*(r3v - r5v)*(r3v - r6v))
     $     +(2d0*Log(-(r5v/Mv)))/(Pi*r5v*(-r1v + r5v)*(-r2v + r5v)*(-r3v
     $     + r5v)*(-r4v + r5v)*(r5v - r6v))
      end function
      
      function b22(MV2, pi, p1x, p2x, p2y, p3x, p3y)
      real*8 MV2, pi, p1x, p2x, p2y, p3x, p3y
      real*8 b22
      b22 = 1d0/((MV2 + p1x**2)*(MV2 + p2x**2 + p2y**2)*(MV2 + p1x**2 +
     $     2d0*p1x*p3x + p3x**2 + p3y**2))
      end function

c......................................................................
      end module nonfact_expressions
