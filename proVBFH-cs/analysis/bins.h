c     -*- Fortran -*-
      integer stxsmjjnbins, stxspthnbins, stxspthjjnbins,
     $     stxsdelphijjnbins
      integer mjjnbins, pthnbins, delphijjnbins, delyjjnbins,
     $     njetnbins, ptjnbins
      parameter (stxsmjjnbins = 4, stxspthnbins = 2, stxspthjjnbins = 2,
     $     stxsdelphijjnbins =2)
      parameter (mjjnbins = 5, pthnbins = 6, delphijjnbins = 4)
      parameter (delyjjnbins = 5, njetnbins = 4, ptjnbins = 2)
      character * 100 stxsmjjstr(stxsmjjnbins), stxspthstr(stxspthnbins)
      character * 100 stxspthjjstr(stxspthjjnbins),
     $     stxsdelphijjstr(stxsdelphijjnbins)
      character * 100 mjjstr(0:mjjnbins), pthstr(0:pthnbins)
      character * 100 delphijjstr(0:delphijjnbins)
      character * 100 delyjjstr(0:delyjjnbins), njetstr(0:njetnbins),
     $     ptjstr(0:ptjnbins)

      common/binning/stxsmjjstr, stxspthstr, stxspthjjstr,
     $     stxsdelphijjstr, mjjstr,pthstr, delphijjstr, delyjjstr,
     $     njetstr, ptjstr

      save /binning/


