!----------------------------------------------------------------------
! Exclusive part of the line-by-line projection-to-Born at NLO
! (docs/DESIGN.md, stage 1): the (1,0) + (0,1) contributions, i.e. VBF
! H + 3 partons at tree level with the extra parton on line 1 or line 2,
! each with its counterevent at the VBF Born kinematics of the same point.
!
! One point (13 random numbers): a VBF Born point (1:7, the POWHEG VBF_H
! phase space of the inclusive code), radiation on line 1 (8:10) and on
! line 2 (11:13) with q1, q2, the Higgs and the other line kept fixed
! (cs_kinematics). The analysis receives
!   the line-1 event with weight w1, the line-2 event with weight w2,
!   and the Born event with weight -(w1 + w2),
! accumulated as one entry (pwhgaccumup), so that they are correlated in
! the error estimate. The function returns |w1| + |w2| for VEGAS's
! adaptation; the signed total of a point is zero.
!
! Flavours: the matrix element depends only on the class of each line
! (NC: u- or d-type quark or antiquark; CC: W+ or W- emitting quark or
! antiquark; for an incoming gluon the class of the q qbar pair), for
! massless quarks and a diagonal CKM matrix. It is evaluated once per pair
! of classes with a representative flavour, and multiplied by the PDF sum
! over the class members.
!----------------------------------------------------------------------
module cs_exclusive
  use types, only: dp
  use hoppet
  use incl_parameters
  use phase_space, only: gen_phsp, set_phsp, x1, x2, vq1, vq2, Q1_sq, Q2_sq
  use cs_kinematics
  implicit none
  private
  public :: cs_excl_dsigma, cs_excl_setup, excl_fill, excl_npow, excl_cutoff, excl_stats, excl_flavcheck

  logical, save :: excl_fill = .false.
  integer, save :: excl_npow = 2
  real(dp), save :: excl_cutoff = 1d-8
  ! counters: points, rejected by cutoff, NaN
  integer(8), save :: excl_stats(3) = 0
  ! debug: compare the class sums with an explicit flavour loop on the
  ! first excl_flavcheck points
  integer, save :: excl_flavcheck = 0

  ! classes of a line transition a -> b: representative flavours and
  ! the member list (up to 3)
  integer, parameter :: maxmem = 3
  type line_class
     integer :: a, b            ! representative incoming, outgoing flavour
     integer :: n               ! number of members
     integer :: ma(maxmem), mb(maxmem)
     integer :: w               ! 0 NC, +1 emits W+, -1 emits W-
  end type line_class
  type(line_class), save :: cls(8)
  integer, save :: ncls = 8
  ! gluon-initiated classes (g -> q qbar'): representative (quark, antiquark)
  type gluon_class
     integer :: q, qb, n, w
  end type gluon_class
  type(gluon_class), save :: gcls(4)
  integer, save :: ngcls = 4

contains

  subroutine cs_excl_setup()
    ! NC classes
    cls(1) = line_class( 2, 2, 2, [ 2, 4, 0], [ 2, 4, 0],  0)   ! u, c
    cls(2) = line_class( 1, 1, 3, [ 1, 3, 5], [ 1, 3, 5],  0)   ! d, s, b
    cls(3) = line_class(-2,-2, 2, [-2,-4, 0], [-2,-4, 0],  0)   ! ubar, cbar
    cls(4) = line_class(-1,-1, 3, [-1,-3,-5], [-1,-3,-5],  0)   ! dbar, sbar, bbar
    ! CC classes (diagonal CKM, no top): emits W+ (charge of the line drops by 1)
    cls(5) = line_class( 2, 1, 2, [ 2, 4, 0], [ 1, 3, 0], +1)   ! u -> d, c -> s
    cls(6) = line_class(-1,-2, 2, [-1,-3, 0], [-2,-4, 0], +1)   ! dbar -> ubar, sbar -> cbar
    ! emits W-
    cls(7) = line_class( 1, 2, 2, [ 1, 3, 0], [ 2, 4, 0], -1)   ! d -> u, s -> c
    cls(8) = line_class(-2,-1, 2, [-2,-4, 0], [-1,-3, 0], -1)   ! ubar -> dbar, cbar -> sbar
    ! gluon-initiated: g -> q (outgoing quark) + qbar (outgoing antiquark);
    ! the line runs from the antiquark (a crossed incoming quark) to the quark
    gcls(1) = gluon_class( 2, -2, 2,  0)    ! g -> u ubar, c cbar
    gcls(2) = gluon_class( 1, -1, 3,  0)    ! g -> d dbar, s sbar, b bbar
    gcls(3) = gluon_class( 1, -2, 2, +1)    ! g -> d ubar (as u -> d, W+), s cbar
    gcls(4) = gluon_class( 2, -1, 2, -1)    ! g -> u dbar (as d -> u, W-), c sbar
  end subroutine cs_excl_setup

  ! can a line of W charge w1 be combined with one of w2 in VBF H?
  logical function compatible(w1, w2)
    integer, intent(in) :: w1, w2
    compatible = (w1 == 0 .and. w2 == 0) .or. (w1 /= 0 .and. w1 == -w2)
  end function compatible

  real(dp) function pdfsum(f, c)
    real(dp), intent(in) :: f(-6:6)
    type(line_class), intent(in) :: c
    integer :: i
    pdfsum = 0
    do i = 1, c%n
       pdfsum = pdfsum + f(c%ma(i))
    enddo
  end function pdfsum

  double precision function cs_excl_dsigma(xrand, vegas_weight)
    real(dp), intent(in) :: xrand(13), vegas_weight
    real(dp) :: pb(0:3,5), xb1, xb2, jacb, sbeams, common, ptH
    real(dp) :: Q1, Q2, mur(2), muf(2), as(2), fB(-6:6,2), fE(-6:6,2)
    real(dp) :: pin(0:3), a(0:3), b(0:3), xp, z, wrad(2), p6(0:3,6,2), w(2)
    real(dp) :: m1, m2, me
    logical :: ok(2)
    integer :: line, i1, i2, bflav(6)
    real(dp), external :: hoppetAlphaS
    integer vegas_ncall
    common/vegas_ncall/vegas_ncall

    cs_excl_dsigma = 0
    excl_stats(1) = excl_stats(1) + 1
    call gen_phsp(xrand(1:7))
    call set_phsp()
    call cs_get_born(pb, xb1, xb2, jacb, sbeams)
    if (jacb == 0 .or. min(Q1_sq, Q2_sq) <= Qmin**2) return   ! as the inclusive dsigma
    Q1 = sqrt(Q1_sq); Q2 = sqrt(Q2_sq)
    ptH = sqrt(pb(1,3)**2 + pb(2,3)**2)
    mur = [cs_mu(xmur, 1, Q1, Q2, ptH, .true.), cs_mu(xmur, 2, Q1, Q2, ptH, .true.)]
    muf = [cs_mu(xmuf, 1, Q1, Q2, ptH, .false.), cs_mu(xmuf, 2, Q1, Q2, ptH, .false.)]
    ! Born-level flux, phase space and conversion to pb
    common = jacb/(2*x1*x2*S)*gev2pb
    ! Born-level PDFs of both lines at their factorisation scales
    call hoppetEval(x1, muf(1), fB(:,1))
    call hoppetEval(x2, muf(2), fB(:,2))
    fB = fB/spread([x1, x2], 1, 13)       ! hoppet returns x f(x)

    ! radiation on each line
    do line = 1, 2
       if (line == 1) then
          call line_radiation(pb(:,1), pb(:,4), x1, xrand(8:10), excl_npow, excl_cutoff, &
               & pin, a, b, xp, z, wrad(1), ok(1))
       else
          call line_radiation(pb(:,2), pb(:,5), x2, xrand(11:13), excl_npow, excl_cutoff, &
               & pin, a, b, xp, z, wrad(2), ok(2))
       endif
       if (.not. ok(line)) then
          excl_stats(2) = excl_stats(2) + 1
          cycle
       endif
       ! H+3j momenta in proVBFH's order: 1, 2 incoming, 3 H, 4, 5 outgoing
       ! quarks of lines 1, 2, 6 the extra parton (on this line)
       p6(:,1:5,line) = pb(:,1:5)
       if (line == 1) then
          p6(:,1,line) = pin; p6(:,4,line) = a
       else
          p6(:,2,line) = pin; p6(:,5,line) = a
       endif
       p6(:,6,line) = b
       as(line) = hoppetAlphaS(mur(line))
       ! PDFs of the radiating line at its new momentum fraction
       call hoppetEval(merge(x1, x2, line == 1)/xp, muf(line), fE(:,line))
       fE(:,line) = fE(:,line)/(merge(x1, x2, line == 1)/xp)
    enddo

    w = 0
    do line = 1, 2
       if (.not. ok(line)) cycle
       ! quark-initiated on the radiating line
       do i1 = 1, ncls
          do i2 = 1, ncls
             if (.not. compatible(cls(i1)%w, cls(i2)%w)) cycle
             bflav = [cls(i1)%a, cls(i2)%a, 25, cls(i1)%b, cls(i2)%b, 0]
             call cs_hjjj_lines(p6(:,:,line), bflav, m1, m2)
             me = merge(m1, m2, line == 1)
             if (line == 1) then
                w(1) = w(1) + me*pdfsum(fE(:,1), cls(i1))*pdfsum(fB(:,2), cls(i2))
             else
                w(2) = w(2) + me*pdfsum(fB(:,1), cls(i1))*pdfsum(fE(:,2), cls(i2))
             endif
          enddo
       enddo
       ! gluon-initiated on the radiating line: g -> q (at 4 or 5) + qbar (at 6)
       do i1 = 1, ngcls
          do i2 = 1, ncls
             if (.not. compatible(gcls(i1)%w, cls(i2)%w)) cycle
             if (line == 1) then
                bflav = [0, cls(i2)%a, 25, gcls(i1)%q, cls(i2)%b, gcls(i1)%qb]
                call cs_hjjj_lines(p6(:,:,line), bflav, m1, m2)
                w(1) = w(1) + m1*gcls(i1)%n*fE(0,1)*pdfsum(fB(:,2), cls(i2))
             else
                bflav = [cls(i2)%a, 0, 25, cls(i2)%b, gcls(i1)%q, gcls(i1)%qb]
                call cs_hjjj_lines(p6(:,:,line), bflav, m1, m2)
                w(2) = w(2) + m2*gcls(i1)%n*pdfsum(fB(:,1), cls(i2))*fE(0,2)
             endif
          enddo
       enddo
       w(line) = w(line)*common*wrad(line)*as(line)
    enddo
    if (excl_flavcheck > 0) then
       excl_flavcheck = excl_flavcheck - 1
       call flavour_check(p6, ok, fB, fE, common*wrad*as, w)
    endif
    if (any(w /= w)) then
       excl_stats(3) = excl_stats(3) + 1
       return
    endif

    if (excl_fill) then
       if (ok(1)) call cs_analysis(6, p6(:,:,1), w(1)*vegas_ncall*vegas_weight)
       if (ok(2)) call cs_analysis(6, p6(:,:,2), w(2)*vegas_ncall*vegas_weight)
       call cs_analysis(5, pb, -(w(1) + w(2))*vegas_ncall*vegas_weight)
       call pwhgaccumup
    endif
    cs_excl_dsigma = abs(w(1)) + abs(w(2))
  end function cs_excl_dsigma

  ! explicit loop over all flavour combinations (diagonal CKM, no top),
  ! independent of the classes, compared with the class sums wcls
  subroutine flavour_check(p6, ok, fB, fE, norm, wcls)
    real(dp), intent(in) :: p6(0:3,6,2), fB(-6:6,2), fE(-6:6,2), norm(2), wcls(2)
    logical, intent(in) :: ok(2)
    real(dp) :: wb(2), m1, m2, f1, f2
    integer :: line, a1, a2, b1, b2, w1, w2, bflav(6), c
    wb = 0
    do line = 1, 2
       if (.not. ok(line)) cycle
       do a1 = -5, 5
          do a2 = -5, 5
             if (a1 == 0 .or. a2 == 0) cycle
             ! NC on both lines, or CC on both with opposite W charges
             do c = 1, 2
                if (c == 1) then
                   b1 = a1; b2 = a2; w1 = 0; w2 = 0
                else
                   call cc_partner(a1, b1, w1); call cc_partner(a2, b2, w2)
                   if (w1 == 0 .or. w2 == 0 .or. w1 /= -w2) cycle
                endif
                bflav = [a1, a2, 25, b1, b2, 0]
                call cs_hjjj_lines(p6(:,:,line), bflav, m1, m2)
                if (line == 1) then
                   wb(1) = wb(1) + m1*fE(a1,1)*fB(a2,2)
                else
                   wb(2) = wb(2) + m2*fB(a1,1)*fE(a2,2)
                endif
             enddo
          enddo
       enddo
       ! incoming gluon on the radiating line: g -> q (b) qbar (-a) with
       ! a -> b a quark transition (NC: b = a; CC: b = partner of a)
       do a1 = 1, 5
          do c = 1, 2
             if (c == 1) then
                b1 = a1; w1 = 0
             else
                call cc_partner(a1, b1, w1)
                if (w1 == 0) cycle
             endif
             do a2 = -5, 5
                if (a2 == 0) cycle
                if (c == 1) then
                   b2 = a2; w2 = 0
                else
                   call cc_partner(a2, b2, w2)
                   if (w2 == 0 .or. w1 /= -w2) cycle
                endif
                if (line == 1) then
                   bflav = [0, a2, 25, b1, b2, -a1]
                   call cs_hjjj_lines(p6(:,:,line), bflav, m1, m2)
                   wb(1) = wb(1) + m1*fE(0,1)*fB(a2,2)
                else
                   bflav = [a2, 0, 25, b2, b1, -a1]
                   call cs_hjjj_lines(p6(:,:,line), bflav, m1, m2)
                   wb(2) = wb(2) + m2*fB(a2,1)*fE(0,2)
                endif
             enddo
          enddo
       enddo
    enddo
    wb = wb*norm
    write(6,'(a,2es16.8,a,2es10.2)') ' flavcheck: w1, w2 = ', wcls, '  explicit/classes - 1: ', &
         & merge(wb/wcls - 1, 0.0_dp*wb, wcls /= 0)
  end subroutine flavour_check

  ! CC transition a -> b of a quark line and the charge w of the emitted W
  ! (w = 0 if none: b quark, no top)
  subroutine cc_partner(a, b, w)
    integer, intent(in) :: a
    integer, intent(out) :: b, w
    integer :: aa
    aa = abs(a); b = 0; w = 0
    if (aa == 5) return
    if (mod(aa,2) == 0) then
       b = sign(aa-1, a)          ! u -> d, c -> s (ubar -> dbar)
    else
       b = sign(aa+1, a)          ! d -> u, s -> c
    endif
    ! charge of the line drops by w: u -> d emits W+, ubar -> dbar W-
    if (a > 0) then
       w = merge(+1, -1, mod(aa,2) == 0)
    else
       w = merge(-1, +1, mod(aa,2) == 0)
    endif
  end subroutine cc_partner

  ! the scales of the inclusive code (muR1, muR2, muF1, muF2 in
  ! ../proVBFH/src/inclusive/matrix_element.f90, private there; keep in step)
  real(dp) function cs_mu(xfac, line, Q1, Q2, ptH, ren)
    real(dp), intent(in) :: xfac, Q1, Q2, ptH
    integer, intent(in) :: line
    logical, intent(in) :: ren
    real(dp) :: Q
    Q = merge(Q1, Q2, line == 1)
    if (scale_choice <= 1) then
       if (ren) then
          cs_mu = sf_muR(Q)
       else
          cs_mu = sf_muF(Q)
       endif
    elseif (scale_choice == 2) then
       cs_mu = xfac*sqrt(Q1*Q2)
    elseif (scale_choice == 3) then
       cs_mu = xfac*((mh*0.5d0)**4 + (mh*ptH*0.5d0)**2)**0.25d0
    else
       stop 'cs_mu: illegal scale_choice'
    endif
  end function cs_mu

  subroutine cs_analysis(n, p, wgt)
    integer, intent(in) :: n
    real(dp), intent(in) :: p(0:3,n), wgt
    call cs_fill_phep(n, p)
    call user_analysis(wgt)
  end subroutine cs_analysis

end module cs_exclusive
