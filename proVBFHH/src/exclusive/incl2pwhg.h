c -*- Fortran -*-

      double precision current_www,incl_tot,x1,x2,vq1(0:3),vq2(0:3)
      double precision vpH1(0:3), vpH2(0:3), ptH1H2, Q1_sq, Q2_sq, Qmin
      integer order
      logical btildennloon
      common/Cincl2pwhg/current_www,incl_tot,x1,x2,vq1,vq2,vpH1,vpH2,
     $     ptH1H2, Q1_sq, Q2_sq, Qmin,order, btildennloon
      save /Cincl2pwhg/
