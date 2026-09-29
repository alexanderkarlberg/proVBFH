c One-loop H+3j of the PUBLIC POWHEG-BOX-V2 VBF_HJJJ (svn r4135), as
c distributed, at the points of the NNLOJET comparison with isolated
c helicities: s(1) c(2) -> H s(4) c(5) g(6), Z fusion, mu_R = 100 GeV.
c Helicity ih1 (1 = L, 2 = R) is kept on line 1 (d-type, VBFNLO clr type
c 4) and ih2 on line 2 (u-type, type 3), as in harness_virt.f90.
c Reads heli_points.dat (ipt, p1, p2, k1 = s out, k2 = c out, k3 = g) and
c fort80_nf44.dat (ipt, ih1, ih2, ..., vN = NNLOJET V/B with nfC1 = the
c flavour of the non-radiating line, vV = VBFNLO V/B of the tested copy)
c and prints V/B of the public code (units alpha_s/2pi, alpha_s = 1 as
c in the harness), NNLOJET - public - (pi^2/6)(2 CF + CA/2), and public -
c tested copy.
      program virt_heli_public
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      include 'PhysPars.h'
      real * 8 clr, xm2, xmg, bb, vv, aa
      common /BKOPOU/ clr(4,5,-1:1), xm2(6), xmg(6), bb(6,6,6),
     $     vv(4,5), aa(4,5)
      real * 8 clr0(4,5,-1:1), pts(20,3), row(9), p(0:3,nlegborn)
      real * 8 virt, born, bmunu(0:3,0:3,nlegborn),
     $     bornjk(nlegborn,nlegborn), borncol(2), vpub, cst, dmax1, dmax2
      integer HWW, HZZ, ipt, i, j, ih1, ih2, bflav(nlegborn), iu, ir
      character(len=16) wenv
      cst = pi**2/6*(2*4d0/3d0 + 1.5d0)

      call init_flsttag
      call init_phys
      call particle_identif(HWW, HZZ)
      st_alpha = 1d0
      st_muren2 = 100d0**2
      flst_cur_iborn = 1
      clr0 = clr
      write(*,'(a,6f12.4)') ' BKOPOU masses:   ', (sqrt(xm2(i)), i=1,6)
      write(*,'(a,6f12.4)') ' BKOPOU widths:   ',
     $     (xmg(i)/sqrt(max(xm2(i),1d-30)), i=1,6)
c     WIDTHS=harness: the Z and W widths of the NNLOJET harness (proVBFH
c     reads them from vbfnlo.input; this code computes its own)
      call get_environment_variable('WIDTHS', wenv)
      if (trim(wenv).eq.'harness') then
         xmg(2) = sqrt(xm2(2))*2.4952d0
         xmg(3) = sqrt(xm2(3))*2.141d0
         xmg(4) = sqrt(xm2(4))*2.141d0
         write(*,'(a,6f12.4)') ' widths now:      ',
     $        (xmg(i)/sqrt(max(xm2(i),1d-30)), i=1,6)
      endif
      bflav = (/ 3, 4, HZZ, 3, 4, 0 /)

      open(newunit=iu, file='heli_points.dat', status='old')
      do i = 1, 3
         read(iu,*) ipt, (pts(j,i), j=1,20)
      enddo
      close(iu)
      dmax1 = 0; dmax2 = 0
      open(newunit=ir, file='fort80_nf44.dat', status='old')
      do i = 1, 12
         read(ir,*) ipt, ih1, ih2, (row(j), j=4,9)
         do j = 0, 3
            p(j,1) = pts(1+j,ipt)
            p(j,2) = pts(5+j,ipt)
            p(j,4) = pts(9+j,ipt)
            p(j,5) = pts(13+j,ipt)
            p(j,6) = pts(17+j,ipt)
            p(j,3) = p(j,1) + p(j,2) - p(j,4) - p(j,5) - p(j,6)
         enddo
         clr = clr0
         if (ih1.eq.1) then
            clr(4,2,1) = 0
         else
            clr(4,2,-1) = 0
         endif
         if (ih2.eq.1) then
            clr(3,2,1) = 0
         else
            clr(3,2,-1) = 0
         endif
         call calc_als
         call ctrans(1)
         call compvirt_hjjj(p, bflav, virt)
         call compborn_hjjj(p, bflav, born, bmunu, bornjk, borncol)
         vpub = 2*pi*virt/born
         write(*,'(a,i2,a,2a2,a,f11.5,a,es11.3,a,es11.3)') ' point', ipt,
     $        ' hel ', merge('L','R',ih1.eq.1), merge('L','R',ih2.eq.1),
     $        '  V/B public:', vpub, '  NNLOJET - public - const:',
     $        row(8) - vpub - cst, '  public - tested copy:',
     $        vpub - row(9)
         dmax1 = max(dmax1, abs(row(8) - vpub - cst))
         dmax2 = max(dmax2, abs(vpub - row(9)))
      enddo
      close(ir)
      clr = clr0
      write(*,'(a,es10.2,a,es10.2)') ' max |NNLOJET - public - const|:',
     $     dmax1, '   max |public - tested copy|:', dmax2
      end
