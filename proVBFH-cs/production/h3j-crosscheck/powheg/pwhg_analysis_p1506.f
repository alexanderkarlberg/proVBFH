c     Hook-up of proVBFH-cs/analysis/p1506_analysis.f (copied unchanged
c     into this directory, with proVBFH's phspcuts.h) as the POWHEG-BOX-V2
c     analysis: pwhg_analysis_driver calls analysis(dsig0) once per point
c     with the full weight; p1506's user_analysis does the rest.
      subroutine analysis(dsig0)
      implicit none
      real * 8 dsig0
      call user_analysis(dsig0)
      end
