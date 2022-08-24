c -*- Fortran -*-

      logical phspcuts, passed_cuts,passed_cuts_ig(0:4),dummy_analysis
      logical jet_opphem
      double precision  vbftot, vbftot_save
      double precision ptalljetmin, ptjetmin,yjetmax,mjjmin,deltay_jjmin
      double precision Rsep_jjmin
      common/cphspcuts/vbftot, vbftot_save,ptalljetmin, ptjetmin,yjetmax
     $     ,mjjmin,deltay_jjmin, Rsep_jjmin,phspcuts,passed_cuts
     $     ,passed_cuts_ig, dummy_analysis,jet_opphem
      save /cphspcuts/
