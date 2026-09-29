c Copies of POWHEG-BOX routines (from ../proVBFH/src/powheg-files/find_regions.f,
c unchanged) needed by proVBFH's init_processes.f to encode the flavour-list
c tags; the rest of find_regions.f (FKS regions) is not used by proVBFH-cs.
      subroutine doubletag_entry(born_or_real,tag_f,tag,leg,graph)
c     This soubroutine is invoked if "double tags" are needed, i.e.
c     lines are tagged with a flavour conserved number (tag_f) and with
c     a second tag that also applies to gluon lines (tag).
c     Instead of assigning tags in the user init_processes one calls
c     call doubletag_entry('born',tag_f,tag,leg,graph)
c     or
c     call doubletag_entry('born',tag_f,tag,leg,graph)
c     The POWHEG BOX takes care of encoding the tags in the tag arrays 
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_flg.h'
      character * 4 born_or_real
      integer tag_f,tag,leg,graph
      integer, allocatable, save :: tagborn(:,:,:),tagreal(:,:,:)
      integer maxtag_f,ileg,igraph
      logical :: ini=.true.
      integer tag_factor
      common/cdoubletag/tag_factor
      if(ini) then
         flg_doubletags = .true.
         allocate(tagborn(2,nlegborn,maxprocborn))
         allocate(tagreal(2,nlegreal,maxprocreal))
         tagborn = 0
         tagreal = 0
         ini = .false.
      endif
      if(born_or_real.eq.'born') then
         tagborn(1,leg,graph) = tag_f
         tagborn(2,leg,graph) = tag
      elseif(born_or_real.eq.'real') then
         tagreal(1,leg,graph) = tag_f
         tagreal(2,leg,graph) = tag
      elseif(born_or_real.eq.'fin ') then
c finalize tag array
         maxtag_f = 0
         do igraph=1,flst_nborn
            do ileg=1,nlegborn
               maxtag_f = max(maxtag_f,tagborn(2,ileg,igraph))
            enddo
         enddo
         do igraph=1,flst_nreal
            do ileg=1,nlegreal
               maxtag_f = max(maxtag_f,tagreal(2,ileg,igraph))
            enddo
         enddo
         if(maxtag_f.eq.0) then
            tag_factor = 1
         else
c     find nearest power of 10 larger than maxtag_f
            tag_factor = 10**(1+int(log10(dble(maxtag_f))))
         endif
         do igraph=1,flst_nborn
            do ileg=1,nlegborn
               flst_borntags(ileg,igraph) = tagborn(2,ileg,igraph)
     1              + tag_factor * tagborn(1,ileg,igraph)
            enddo
         enddo
         do igraph=1,flst_nreal
            do ileg=1,nlegreal
               flst_realtags(ileg,igraph) = tagreal(2,ileg,igraph)
     1              + tag_factor * tagreal(1,ileg,igraph)
            enddo
         enddo
         deallocate(tagborn,tagreal)
      endif
      end

      subroutine finalize_tags
      implicit none
      integer dum1,dum2,dum3,dum4
      call doubletag_entry('fin ',dum1,dum2,dum3,dum4)
      end
