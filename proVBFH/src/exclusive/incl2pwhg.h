c -*- Fortran -*-

      double precision current_www,incl_tot,x1,x2,vq1(0:3),vq2(0:3)
      double precision ptH, Q1_sq, Q2_sq, Qmin
      integer order
      logical btildennloon
      common/Cincl2pwhg/current_www,incl_tot,x1,x2,vq1,vq2,ptH,
     $     Q1_sq, Q2_sq, Qmin,order, btildennloon
      save /Cincl2pwhg/
