c -*- Fortran -*-

      real * 8 br_born(maxprocborn),
     #  br_bmunu(0:3,0:3,nlegborn,maxprocborn),
     #  br_bornjk(nlegborn,nlegborn,maxprocborn)
      common/pwhg_br/br_born,br_bmunu,br_bornjk

      save /pwhg_br/
