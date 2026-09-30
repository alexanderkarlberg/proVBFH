c Glue for stage 2 of proVBFH-cs (the (2,0) and (0,2) contributions,
c docs/DESIGN.md): proVBFH's real flavour list with its line tags (built
c by its init_processes), and the VBFNLO matrix elements for one line:
c H+4j (both extra partons on the line), H+3j with the spin correlations
c of the gluon, and the one-loop H+3j.
c All matrix elements are at alpha_s = 1 (st_alpha = 1), averaged over
c initial spins and colours, with proVBFH's symmetry factors.
c H+3j momenta and flavours in the order of cs_glue.f (1, 2 incoming
c partons of lines 1 and 2, 3 Higgs, 4, 5 outgoing quarks of lines 1 and
c 2, 6 the extra parton); real entries in their own leg order.
c---------------------------------------------------------------------
      subroutine cs_init_reals(nreal)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      integer nreal
      call init_processes
      call finalize_tags
      nreal = flst_nreal
      end
c---------------------------------------------------------------------
      subroutine cs_get_real(j,flav,tags)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      integer j,flav(nlegreal),tags(nlegreal),k
      do k=1,nlegreal
         flav(k) = flst_real(k,j)
         tags(k) = flst_realtags(k,j)
      enddo
      end
c---------------------------------------------------------------------
c |M|^2 of real entry j (setreal, with its tags), momenta p in the
c entry's leg order
      subroutine cs_real_me(j,p,amp2)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_flst_2.h'
      include 'pwhg_math.h'
      include 'phspcuts.h'
      double precision p(0:3,nlegreal),amp2
      integer j,rflav(nlegreal),k
      logical save
      do k=1,nlegreal
         rflav(k) = flst_real(k,j)
      enddo
      realequiv_tag = j
      alr_tag = 0
c     setreal skips entries by the cut flags of proVBFH's phspcuts
c     (which would need alr_tag); proVBFH-cs applies its own cuts
      save = phspcuts
      phspcuts = .false.
      call setreal(p,rflav,amp2)
      phspcuts = save
c     setreal returns |M|^2/(alpha_s/(2 pi)); alpha_s = 1
      amp2 = amp2/(2*pi)
      end
c---------------------------------------------------------------------
c tags of the H+3j Born that select the extra parton on line (1 or 2):
c the VBFNLO Born and virtual then keep only that line's contribution
      subroutine cs_set_borntags(line)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      integer line
      flst_cur_iborn = 1
      flst_borntags(1,1) = 1
      flst_borntags(2,1) = 2
      flst_borntags(3,1) = 0
      flst_borntags(4,1) = 1
      flst_borntags(5,1) = 2
      flst_borntags(6,1) = line
      end

      subroutine cs_clear_borntags
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      integer k
      do k=1,nlegborn
         flst_borntags(k,1) = 0
      enddo
      end
c---------------------------------------------------------------------
c H+3j Born with the extra parton emitted from line: |M|^2 (born) and,
c if the line has a gluon, B^{mu nu} of that gluon (bmunu, with
c -g_{mu nu} B^{mu nu} = born; zero otherwise)
      subroutine cs_hjjj_born_line(p,bflav,line,born,bmunu)
      implicit none
      include 'nlegborn.h'
      double precision p(0:3,nlegborn),born,bmunu(0:3,0:3)
      integer bflav(nlegborn),line
      double precision b,bmn(0:3,0:3,nlegborn),bjk(nlegborn,nlegborn),
     $     bcol(2)
      integer mu,nu,j
      call calc_als
      call ctrans(1)
      call cs_set_borntags(line)
      call compborn_hjjj(p,bflav,b,bmn,bjk,bcol)
      call cs_clear_borntags
      born = bcol(line)
      do mu=0,3
         do nu=0,3
            bmunu(mu,nu) = 0
         enddo
      enddo
      do j=1,nlegborn
         if (bflav(j).eq.0) then
            do mu=0,3
               do nu=0,3
                  bmunu(mu,nu) = bmn(mu,nu,j)
               enddo
            enddo
         endif
      enddo
      end
c---------------------------------------------------------------------
c one-loop H+3j with the extra parton emitted from line: the finite part
c 2 Re(M0 M1*) of proVBFH's VBFNLO virtual (CDR, the factor
c (4 pi)^eps/Gamma(1-eps) (mu^2)^eps removed, mu^2 = mur2), both loops:
c on the radiating line and the vertex correction of the other line
      subroutine cs_hjjj_virt_line(p,bflav,line,mur2,virt)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_st.h'
      double precision p(0:3,nlegborn),mur2,virt
      integer bflav(nlegborn),line
      double precision save
      save = st_muren2
      st_muren2 = mur2
      call calc_als
      call ctrans(1)
      call cs_set_borntags(line)
      call compvirt_hjjj(p,bflav,virt)
      call cs_clear_borntags
      st_muren2 = save
      end
c---------------------------------------------------------------------
c number of real entries (after cs_init_reals)
      integer function cs_nreal()
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      cs_nreal = flst_nreal
      end
c---------------------------------------------------------------------
c cs_estimate 5: switch the original (buggy) NC gg pair type of
c real_vbfnlo.f on or off
      subroutine cs_set_ggbug(flag)
      implicit none
      logical flag
      logical cs_ggbug
      common/csggbug/cs_ggbug
      cs_ggbug = flag
      end
