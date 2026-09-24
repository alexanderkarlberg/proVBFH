# Debug harness for proVBFH-inclusive (scratch builds only, not for
# commit). Patches a copy of proVBFH-inclusive/src (at edc7a24 or later)
# so that, with -tensorME, every phase-space point is also evaluated
# with the analytic matrix element and with the tensor matrix element
# as it was before the fixes, and prints (every 200k points) the
# ratios of the integrals, their MC errors, pointwise deviations and
# the fraction of points / of |sigma| above given deviations.
#
#   git show 8d00b5d:proVBFH-inclusive/src/matrix_element_dihiggs.f90 > $SP/me_old.f90
#   python3 dbg_patch.py <copy>/src $SP
#   make provbfhh_incl, run with -tensorME
#
# Environment: NOF3=1 sets F3 to zero in all three matrix elements;
# MIRROR=1 averages the (new) tensor matrix element over the point and
# its mirror image y -> -y, which removes its parity-odd part.
import sys
src=sys.argv[1]; sp=sys.argv[2]
old=open(sp+'/me_old.f90').read()
a=old.index("  function eval_matrix_element_tensor(")
b=old.index("  end function eval_matrix_element_tensor\n")+len("  end function eval_matrix_element_tensor\n")
f=old[a:b].replace("eval_matrix_element_tensor","eval_matrix_element_tensor_old")
f=f.replace("    ! Compute hadronic tensors\n","""    if (dbg_noF3) then
       F1(iF3Wp)=zero; F1(iF3Wm)=zero; F1(iF3Z)=zero
       F2(iF3Wp)=zero; F2(iF3Wm)=zero; F2(iF3Z)=zero
    endif
    ! Compute hadronic tensors\n""")
p=src+'/matrix_element_dihiggs.f90'
s=open(p).read()
s=s.replace("  public :: eval_matrix_element_tensor\n","  public :: eval_matrix_element_tensor\n  public :: eval_matrix_element_tensor_old\n  real(dp), public :: dbg_scale\n  logical, public :: dbg_noF3 = .false.\n",1)
zero3="""    if (dbg_noF3) then
       Fx1(iF3Wp,:)=zero; Fx1(iF3Wm,:)=zero; Fx1(iF3Z,:)=zero
       Fx2(iF3Wp,:)=zero; Fx2(iF3Wm,:)=zero; Fx2(iF3Z,:)=zero
    endif
"""
a="    do iorder = order_start,order_stop\n       do i = 1, iorder\n          j = 1 + iorder - i\n          if (WpWm) then"
assert s.count(a)==1; s=s.replace(a, zero3+a)
b="    ! Compute the 3 different contributions\n    sigma = zero"
assert s.count(b)==1; s=s.replace(b, zero3+b)
s=s.replace("  !----------------------------------------------------------------------\n  ! Basis tensors of the hadronic tensor", f+"\n  !----------------------------------------------------------------------\n  ! Basis tensors of the hadronic tensor",1)
s=s.replace("    res = overall_norm * sigma\n  end function eval_matrix_element_tensor\n","""    res = overall_norm * sigma
    dbg_scale = overall_norm * (WW_norm * order_sum(order_start, order_stop, abs(Fx1), abs(Fx2), iWp, iWm, abs(TW)) &
         & + WW_norm * order_sum(order_start, order_stop, abs(Fx1), abs(Fx2), iWm, iWp, abs(TW)) &
         & + ZZ_norm * order_sum(order_start, order_stop, abs(Fx1), abs(Fx2), iZ, iZ, abs(TZ)))
  end function eval_matrix_element_tensor
""",1)
open(p,'w').write(s)
p=src+'/phase_space_dihiggs.f'
s=open(p).read()
a="""C     convert to [pb] and add in jacobian
         dsigma = dsigma * gev2pb * jacobian"""
assert s.count(a)==1
s=s.replace(a, """         if(tensorME) then
            ds2 = dbg_scale
            call dbg_compare(dsigma, eval_matrix_element(
     $        order_min,order_max, x1, x2, kn_beams(:,1), kn_beams(:,2),
     $        vq1, vq2, pH1, pH2, ptH1H2),
     $        eval_matrix_element_tensor_old(
     $        order_min,order_max, x1, x2, kn_beams(:,1), kn_beams(:,2),
     $        vq1, vq2, pH1, pH2, ptH1H2), ds2, vegas_weight*jacobian)
         endif
"""+a)
s=s.replace("""      end function dsigma""","""      end function dsigma

      subroutine dbg_compare(dt, dp, do, sc, w)
      use matrix_element_dihiggs
      implicit none
      double precision dt, dp, do, sc, w
      double precision st, sp, so, sd, sd2, sdo, sdo2, dmax, dsmax, domax, wt(3)
      integer*8 n
      character*8 val
      integer ist
      data st/0d0/, sp/0d0/, so/0d0/, sd/0d0/, sd2/0d0/, dmax/0d0/
      data sdo/0d0/, sdo2/0d0/
      data dsmax/0d0/, domax/0d0/, n/0/
      save st, sp, so, sd, sd2, sdo, sdo2, dmax, dsmax, domax, n, wt
      double precision cnt(5), wsum(5), thr(5), sabs
      integer k
      data cnt/5*0d0/, wsum/5*0d0/, sabs/0d0/
      data thr/1d-14, 1d-12, 1d-10, 1d-8, 1d-4/
      save cnt, wsum, sabs
      if (n.eq.0) then
         call get_environment_variable('NOF3', val, status=ist)
         dbg_noF3 = (ist.eq.0)
         write(*,*) 'DBG noF3 =', dbg_noF3
      endif
      n = n + 1
      st = st + w*dt
      sp = sp + w*dp
      so = so + w*do
      sd = sd + w*(dp-dt)
      sd2 = sd2 + (w*(dp-dt))**2
      sdo = sdo + w*(do-dt)
      sdo2 = sdo2 + (w*(do-dt))**2
      if (dt.ne.0d0.and.abs(dp/dt-1d0).gt.dmax) then
         dmax = abs(dp/dt-1d0)
         wt = (/dt, dp, sc/)
      endif
      if (sc.gt.0d0) dsmax = max(dsmax, abs(dp-dt)/sc)
      if (sc.gt.0d0) domax = max(domax, abs(do-dt)/sc)
      sabs = sabs + abs(w*dt)
      do k = 1, 5
         if (dt.ne.0d0.and.abs(dp/dt-1d0).gt.thr(k)) then
            cnt(k) = cnt(k) + 1d0
            wsum(k) = wsum(k) + abs(w*dt)
         endif
      enddo
      if (mod(n,200000_8).eq.0) then
         write(*,'(a,i10,4es13.5)') ' DBG n, plain/tensor-1, err, '//
     $        'old/new tensor-1, err:', n, sp/st-1d0, sqrt(sd2)/abs(st),
     $        so/st-1d0, sqrt(sdo2)/abs(st)
         write(*,'(a,2es13.5)') ' DBG max |plain-tensor|/scale, '//
     $        '|oldtensor-tensor|/scale:', dsmax, domax
         write(*,'(a,5es10.2)') ' DBG frac of points |p/t-1|>'//
     $        '1e-14,-12,-10,-8,-4:', cnt/dble(n)
         write(*,'(a,5es10.2)') ' DBG frac of |sigma| from them:     ',
     $        wsum/sabs
         write(*,'(a,4es13.5)') ' DBG worst rel point: dev, t, p, '//
     $        'scale:', dmax, wt
      endif
      end subroutine""")
open(p,'w').write(s)

# mirror option and its declarations
s=open(p).read()
assert s.count('         if(tensorME) then\n            ds2 = dbg_scale\n            call dbg_compare(dsigma, eval_matrix_element(')==1 and s.count('      double precision ptH1H2, dsig_temp,ds2')==1
s=s.replace('         if(tensorME) then\n            ds2 = dbg_scale\n            call dbg_compare(dsigma, eval_matrix_element(',"         if(tensorME) then\n            ds2 = dbg_scale\n            call get_environment_variable('MIRROR', mval, status=mst)\n            if (mst.eq.0) then\n               mq1 = vq1; mq2 = vq2; mh1 = pH1; mh2 = pH2\n               mq1(2) = -mq1(2); mq2(2) = -mq2(2)\n               mh1(2) = -mh1(2); mh2(2) = -mh2(2)\n               mb = kn_beams(0:3,1:2); mb(2,:) = -mb(2,:)\n               dsigma = 0.5d0*(dsigma + eval_matrix_element_tensor(\n     $              order_min,order_max, x1, x2, mb(:,1), mb(:,2),\n     $              mq1, mq2, mh1, mh2, ptH1H2))\n            endif\n            call dbg_compare(dsigma, eval_matrix_element(").replace('      double precision ptH1H2, dsig_temp,ds2','      double precision ptH1H2, dsig_temp,ds2\n      double precision mq1(0:3), mq2(0:3), mh1(0:3), mh2(0:3), mb(0:3,2)\n      character*8 mval\n      integer mst')
open(p,'w').write(s)
