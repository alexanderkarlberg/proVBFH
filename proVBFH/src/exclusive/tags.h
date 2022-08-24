c -*- Fortran -*-

      integer, parameter  :: iqpairtag = 1
C      common/ciqtag/iqpairtag 

C     flst_borntags_todiag used only for sanity checks (could be removed) 
      integer flst_borntags_todiag(nlegborn)
C     tags used to set part of the real amplitudes to zero 
      integer flst_realquarktags(4),flst_realgluontags(2)
      common/ctags/flst_borntags_todiag,flst_realquarktags,
     C     flst_realgluontags
