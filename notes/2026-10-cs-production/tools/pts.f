c     common test point: p(0:3,6), 1,2 incoming, 3 Higgs, 4,5,6 outgoing
      subroutine testpoint(p)
      implicit none
      real*8 p(0:3,6), mh, ft(0:3,4), kb(0:3), a(0:3), b(0:3)
      integer mu
      mh = 125d0
      call mom(60d0, 2.2d0, 0.3d0, ft(0,2))
      call mom(55d0, -2.6d0, 2.9d0, ft(0,3))
      call mom(40d0, 0.4d0, 4.6d0, ft(0,4))
      ft(1,1) = -(ft(1,2)+ft(1,3)+ft(1,4))
      ft(2,1) = -(ft(2,2)+ft(2,3)+ft(2,4))
      ft(3,1) = 30d0
      ft(0,1) = sqrt(mh**2+ft(1,1)**2+ft(2,1)**2+ft(3,1)**2)
      do mu=0,3
         kb(mu) = ft(mu,1)+ft(mu,2)+ft(mu,3)+ft(mu,4)
      enddo
      a = 0; b = 0
      a(0) = (kb(0)+kb(3))/2; a(3) = a(0)
      b(0) = (kb(0)-kb(3))/2; b(3) = -b(0)
      p(:,1) = a; p(:,2) = b; p(:,3) = ft(:,1)
      p(:,4) = ft(:,2); p(:,5) = ft(:,3); p(:,6) = ft(:,4)
      end
      subroutine mom(pt, y, phi, p)
      implicit none
      real*8 pt, y, phi, p(0:3)
      p(0) = pt*cosh(y); p(1) = pt*cos(phi); p(2) = pt*sin(phi)
      p(3) = pt*sinh(y)
      end
