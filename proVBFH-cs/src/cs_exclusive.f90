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
! Stage 2 (excl_order >= 2): in addition the (2,0) and (0,2)
! contributions (cs_nlo2), i.e. the NLO correction of the line that
! radiated: per line, at the same three-parton event, the one-loop H+3j
! with the I operator and the K + P terms, and a four-parton event (seven
! more random numbers, 14:20, not adapted by VEGAS) with its six dipole
! counterevents. Stage 3 (excl_order = 3): the (1,1) contribution, i.e.
! each line's O(alpha_s) pieces times the other line's, and the H+4j event
! with one extra parton on each line. All have their counterevent at the
! VBF Born point, which receives minus the sum.
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
  use cs_nlo2
  use cs_memo
  use cs_smallx
  use integration, only: vegas_gdens
  implicit none
  private
  public :: cs_excl_dsigma, cs_excl_setup, excl_fill, excl_npow, excl_cutoff, excl_stats
  public :: excl_hardfrac, excl_order, cs_excl_setup2

  ! 1: (1,0) + (0,1) (NLO); 2: also (2,0) + (0,2); 3: also (1,1) (NNLO)
  integer, save :: excl_order = 1
  ! stage 3: for each pair of line states (1..8: q -> q g of class cls(s);
  ! 9..12: g -> q qbar of gcls(s-8)) the class-12 real entry, the leg of
  ! line 1's extra parton in it (6 or 7), and a normalisation factor
  integer, save :: e3_entry(12,12) = 0, e3_leg1(12,12) = 0
  real(dp), save :: e3_fac(12,12) = 1
  ! g -> q qbar on a line with the antiquark in the outgoing slot (4 or 5)
  ! and the quark as the extra parton (the NC gg entries have this on line 1)
  logical, save :: e3_swap(2,12,12) = .false.
  ! on-the-fly scale variations (cs_scales 3 or 7; 1: off): the event
  ! weights at mu_R = sc_r(k) mu_R0, mu_F = sc_f(k) mu_F0, k = 1..excl_nscale,
  ! written as weights 1..n (files -W1.top ... -Wn.top); 3: symmetric
  ! (1,1), (1/2,1/2), (2,2); 7: the same, then (1/2,1), (1,1/2), (2,1), (1,2)
  integer, parameter :: max_scales = 7
  integer, public, save :: excl_nscale = 1
  real(dp), public, save :: sc_r(max_scales) = [1.0_dp, 0.5_dp, 2.0_dp, 0.5_dp, 1.0_dp, 2.0_dp, 1.0_dp]
  real(dp), public, save :: sc_f(max_scales) = [1.0_dp, 0.5_dp, 2.0_dp, 1.0_dp, 0.5_dp, 1.0_dp, 2.0_dp]
  public :: cs_scales_setup
  real(dp), save :: tb_pb(0:3,5,2), tb_xb(2,2)
  integer, save :: tb_n = 0

  logical, save :: excl_fill = .false.
  integer, save :: excl_npow = 0     ! 0: logarithmic sampling
  ! fraction of the line radiations sampled in the hard channel (small
  ! xp, z ~ 1/2; line_radiation), for the high-pT tails; 0: off;
  ! default 0.3 (production)
  real(dp), save :: excl_hardfrac = 0.3_dp
  real(dp), save :: excl_cutoff = 1d-5
  ! counters: points, line radiations rejected by the cutoff, NaN, points
  ! where neither event nor Born passes the cuts
  integer(8), save :: excl_stats(4) = 0
  ! Points where neither the events nor the Born pass the VBF cuts are
  ! skipped (as proVBFH's phspcuts; exact for all VBF-cut histograms, and
  ! the exclusive part does not contribute to histograms without cuts: the
  ! Higgs momentum is unchanged by the projection); loose cuts (tagging
  ! jet pt and mjj only) while building the grid.

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

  ! VEGAS integrand. The histograms (pwhg_bookhist-multi) are normalised
  ! by the number of pwhgaccumup calls, so every point must be
  ! accumulated, also points that return early with zero weight.
  double precision function cs_excl_dsigma(xrand, vegas_weight)
    real(dp), intent(in) :: xrand(*), vegas_weight
    cs_excl_dsigma = excl_point(xrand, vegas_weight)
    if (excl_fill) call pwhgaccumup
  end function cs_excl_dsigma

  double precision function excl_point(xrand, vegas_weight_in) result(cs_excl_dsigma)
    real(dp), intent(in) :: xrand(*), vegas_weight_in
    real(dp) :: vegas_weight
    ! small-x channel (cs_smallx): the Born variables of channel B, and the factor g_D/rho
    real(dp) :: sx_fac, sx_xr(7), sx_b(5), sx_gA, sx_gD, sx_lnt, sx_u(5), sx_rho_tot
    integer :: sx_ch, sx_j
    real(dp) :: pb(0:3,5), xb1, xb2, jacb, sbeams, common, ptH
    real(dp) :: Q1, Q2, mur(2), muf(2), as(2), fB(-6:6,2), fE(-6:6,2)
    real(dp) :: pin(0:3), a(0:3), b(0:3), xp, z, wrad(2), p6(0:3,6,2), w(2)
    real(dp) :: m1, m2, me
    real(dp) :: xi3(2), wv(2), wr(2), wd(6,2), pr(0:3,7,2), pd(0:3,6,6,2)
    logical :: ok(2), need(2), pass(2), passB
    logical :: okr(2), dok(6,2), needr(2), passr(2), passd(6,2)
    real(dp) :: prt(0:3,7), pdt(0:3,6,6)
    real(dp) :: p7(0:3,7), e12(2), e3
    logical :: ok3, need3, pass3
    logical :: dokt(6), okt
    integer :: m
    logical, external :: cs_passes
    integer :: line, i1, i2, bflav(6)
    real(dp), external :: hoppetAlphaS
    integer vegas_ncall
    common/vegas_ncall/vegas_ncall
    ! weights per scale point (cs_scales; point 1 is the central scale)
    integer :: ks, ns
    real(dp) :: sw(2,max_scales), swv(2,max_scales), swr(2,max_scales), swd(6,2,max_scales)
    real(dp) :: se12(2,max_scales), se3(max_scales), fB0(-6:6,2), fE0(-6:6,2), as0(2)

    cs_excl_dsigma = 0
    vegas_weight = vegas_weight_in
    sx_fac = 1
    sx_xr = xrand(1:7)
    excl_stats(1) = excl_stats(1) + 1
    if (sx_h > 0) then
       ! the Higgs mass, hence taumin, does not depend on b = (r2..r6)
       call gen_phsp(sx_xr)
       call set_phsp()
       call cs_get_born(pb, xb1, xb2, jacb, sbeams)
       sx_lnt = log(mdot(pb(:,3), pb(:,3))/sbeams)
       call sx_choose(xrand(sx_dim0 + 1), sx_h, sx_ch, sx_u(1))
       sx_u(2:5) = xrand(sx_dim0 + 2:sx_dim0 + 5)
       sx_b = xrand(2:6)
       sx_gD = 1
       do sx_j = 1, 5
          sx_gD = sx_gD*vegas_gdens(sx_j + 1, xrand(sx_j + 1))
       enddo
       if (sx_ch > 0) call sx_gen(sx_ch, sx_u, sx_lnt, sx_x0, sx_b)
       sx_gA = 1
       do sx_j = 1, 5
          sx_gA = sx_gA*vegas_gdens(sx_j + 1, sx_b(sx_j))
       enddo
       sx_rho_tot = (1 - sx_h)*sx_gA + sx_h*sx_rho(sx_b, sx_lnt, sx_x0)
       sx_fac = sx_gD/sx_rho_tot
       sx_xr(2:6) = sx_b
       vegas_weight = vegas_weight_in*sx_fac
    endif
    call gen_phsp(sx_xr)
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
               & pin, a, b, xp, z, wrad(1), ok(1), excl_hardfrac)
       else
          call line_radiation(pb(:,2), pb(:,5), x2, xrand(11:13), excl_npow, excl_cutoff, &
               & pin, a, b, xp, z, wrad(2), ok(2), excl_hardfrac)
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
       xi3(line) = merge(x1, x2, line == 1)/xp
       call hoppetEval(xi3(line), muf(line), fE(:,line))
       fE(:,line) = fE(:,line)/xi3(line)
    enddo

    ! four-parton events and their dipole counterevents
    okr = .false.; dok = .false.
    if (excl_order >= 2) then
       do line = 1, 2
          call nlo2_real_kin(line, pb, [x1, x2], xrand(14:20), excl_cutoff, pr(:,:,line), &
               & pd(:,:,:,line), dok(:,line), okr(line))
          if (.not. ok(line)) okr(line) = .false.
       enddo
    endif
    ! stage 3: the (1,1) four-parton event, both lines radiated (their
    ! three-parton points; proVBFH order: 4, 5 outgoing quarks, 6, 7 the
    ! extra partons of lines 1, 2)
    ok3 = excl_order >= 3 .and. ok(1) .and. ok(2)
    if (ok3) then
       p7(:,1) = p6(:,1,1); p7(:,2) = p6(:,2,2); p7(:,3) = pb(:,3)
       p7(:,4) = p6(:,4,1); p7(:,5) = p6(:,5,2); p7(:,6) = p6(:,6,1); p7(:,7) = p6(:,6,2)
    endif

    ! cut decisions before the matrix elements: a line's weight is needed
    ! only if its event or the Born event passes
    need = ok
    needr = okr
    need3 = ok3
    pass3 = .false.
    passB = cs_passes(5, pb, .not. excl_fill)
    if (ok3) then
       pass3 = cs_passes(7, p7, .not. excl_fill)
       need3 = pass3 .or. passB
    endif
    do line = 1, 2
       if (ok(line)) pass(line) = cs_passes(6, p6(:,:,line), .not. excl_fill)
       need(line) = ok(line) .and. (pass(line) .or. passB)
       if (okr(line)) then
          passr(line) = cs_passes(7, pr(:,:,line), .not. excl_fill)
          do m = 1, 6
             passd(m,line) = .false.
             if (dok(m,line)) passd(m,line) = cs_passes(6, pd(:,:,m,line), .not. excl_fill)
          enddo
          needr(line) = passr(line) .or. any(passd(:,line)) .or. passB
       endif
    enddo
    if (.not. any(need) .and. .not. any(needr) .and. .not. need3) then
       excl_stats(4) = excl_stats(4) + 1
       return
    endif

    ! the event weights at each scale point; the matrix elements are
    ! cached (memo) for the points after the first
    ns = excl_nscale
    if (ns > 1) then
       call memo_new_point()
       fB0 = fB; fE0 = fE; as0 = as
    endif
    do ks = 1, ns
    if (ks > 1) then
       mur = [cs_mu(xmur*sc_r(ks), 1, Q1, Q2, ptH, .true.), cs_mu(xmur*sc_r(ks), 2, Q1, Q2, ptH, .true.)]
       muf = [cs_mu(xmuf*sc_f(ks), 1, Q1, Q2, ptH, .false.), cs_mu(xmuf*sc_f(ks), 2, Q1, Q2, ptH, .false.)]
       call hoppetEval(x1, muf(1), fB(:,1))
       call hoppetEval(x2, muf(2), fB(:,2))
       fB = fB/spread([x1, x2], 1, 13)
       do line = 1, 2
          if (.not. ok(line)) cycle
          as(line) = hoppetAlphaS(mur(line))
          call hoppetEval(xi3(line), muf(line), fE(:,line))
          fE(:,line) = fE(:,line)/xi3(line)
       enddo
    endif
    w = 0
    do line = 1, 2
       if (.not. need(line)) cycle
       ! quark-initiated on the radiating line
       do i1 = 1, ncls
          do i2 = 1, ncls
             if (.not. compatible(cls(i1)%w, cls(i2)%w)) cycle
             bflav = [cls(i1)%a, cls(i2)%a, 25, cls(i1)%b, cls(i2)%b, 0]
             call memo_lines(p6(:,:,line), bflav, m1, m2)
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
                call memo_lines(p6(:,:,line), bflav, m1, m2)
                w(1) = w(1) + m1*gcls(i1)%n*fE(0,1)*pdfsum(fB(:,2), cls(i2))
             else
                bflav = [cls(i2)%a, 0, 25, cls(i2)%b, gcls(i1)%q, gcls(i1)%qb]
                call memo_lines(p6(:,:,line), bflav, m1, m2)
                w(2) = w(2) + m2*gcls(i1)%n*pdfsum(fB(:,1), cls(i2))*fE(0,2)
             endif
          enddo
       enddo
       w(line) = w(line)*common*wrad(line)*as(line)
    enddo

    ! (2,0) and (0,2)
    wv = 0; wr = 0; wd = 0
    if (excl_order >= 2) then
       do line = 1, 2
          if (need(line)) wv(line) = virt_point(line, p6(:,:,line), fB, fE(:,line), xi3(line), &
               & mur(line), muf(line), merge(Q2, Q1, line == 1))*common*wrad(line)*as(line)**2
          if (needr(line)) then
             ! nlo2_real_me works on the state of the last nlo2_real_kin call
             ! (line 2 after the loop above): set up this line's point again
             call nlo2_real_kin(line, pb, [x1, x2], xrand(14:20), excl_cutoff, prt, pdt, dokt, okt)
             if (.not. okt .or. any(dokt .neqv. dok(:,line))) stop 'cs_exclusive: nlo2_real_kin not reproducible'
             call nlo2_real_me(muf(line), fB(:,3-line), wr(line), wd(:,line))
             wr(line) = wr(line)*common*as(line)**2
             wd(:,line) = wd(:,line)*common*as(line)**2
          endif
       enddo
    endif
    ! (1,1): E1, E2 (a line radiated, the other's O(alpha_s) pieces) and E3
    e12 = 0; e3 = 0
    if (excl_order >= 3) then
       if (need(1)) e12(1) = nlo11_e12(1, pb, [x1, x2], fB, p6(:,:,1), fE(:,1), wrad(1), ok(2), &
            & p6(:,:,2), fE(:,2), wrad(2), mur, muf, [Q1, Q2], common)
       if (need(2)) e12(2) = nlo11_e12(2, pb, [x1, x2], fB, p6(:,:,2), fE(:,2), wrad(2), ok(1), &
            & p6(:,:,1), fE(:,1), wrad(1), mur, muf, [Q1, Q2], common)
       if (need3) e3 = nlo11_e3(pb, p6, fE, wrad, mur, common, p7)
    endif
    sw(:,ks) = w; swv(:,ks) = wv; swr(:,ks) = wr; swd(:,:,ks) = wd; se12(:,ks) = e12; se3(ks) = e3
    enddo
    ! the central-scale weights for everything below except the analysis
    if (ns > 1) then
       mur = [cs_mu(xmur, 1, Q1, Q2, ptH, .true.), cs_mu(xmur, 2, Q1, Q2, ptH, .true.)]
       muf = [cs_mu(xmuf, 1, Q1, Q2, ptH, .false.), cs_mu(xmuf, 2, Q1, Q2, ptH, .false.)]
       w = sw(:,1); wv = swv(:,1); wr = swr(:,1); wd = swd(:,:,1); e12 = se12(:,1); e3 = se3(1)
       fB = fB0; fE = fE0; as = as0
    endif
    if (any(sw(:,1:ns) /= sw(:,1:ns)) .or. any(swv(:,1:ns) /= swv(:,1:ns)) .or. any(swr(:,1:ns) /= swr(:,1:ns)) &
         & .or. any(swd(:,:,1:ns) /= swd(:,:,1:ns)) .or. any(se12(:,1:ns) /= se12(:,1:ns)) .or. any(se3(1:ns) /= se3(1:ns))) then
       excl_stats(3) = excl_stats(3) + 1
       return
    endif

    if (excl_fill .and. ns == 1) then
       do line = 1, 2
          if (need(line)) call cs_analysis(6, p6(:,:,line), (w(line) + wv(line) + e12(line))*vegas_ncall*vegas_weight)
          if (needr(line)) then
             call cs_analysis(7, pr(:,:,line), wr(line)*vegas_ncall*vegas_weight)
             do m = 1, 6
                if (dok(m,line)) call cs_analysis(6, pd(:,:,m,line), wd(m,line)*vegas_ncall*vegas_weight)
             enddo
          endif
       enddo
       if (need3) call cs_analysis(7, p7, e3*vegas_ncall*vegas_weight)
       call cs_analysis(5, pb, -(sum(w) + sum(wv) + sum(wr) + sum(wd) + sum(e12) + e3)*vegas_ncall*vegas_weight)
    elseif (excl_fill) then
       do line = 1, 2
          if (need(line)) call cs_analysis_m(6, p6(:,:,line), &
               & (sw(line,1:ns) + swv(line,1:ns) + se12(line,1:ns))*vegas_ncall*vegas_weight)
          if (needr(line)) then
             call cs_analysis_m(7, pr(:,:,line), swr(line,1:ns)*vegas_ncall*vegas_weight)
             do m = 1, 6
                if (dok(m,line)) call cs_analysis_m(6, pd(:,:,m,line), swd(m,line,1:ns)*vegas_ncall*vegas_weight)
             enddo
          endif
       enddo
       if (need3) call cs_analysis_m(7, p7, se3(1:ns)*vegas_ncall*vegas_weight)
       call cs_analysis_m(5, pb, -(sum(sw(:,1:ns), 1) + sum(swv(:,1:ns), 1) + sum(swr(:,1:ns), 1) &
            & + sum(sum(swd(:,:,1:ns), 1), 1) + sum(se12(:,1:ns), 1) + se3(1:ns))*vegas_ncall*vegas_weight)
    endif
    cs_excl_dsigma = abs(w(1) + wv(1) + e12(1)) + abs(w(2) + wv(2) + e12(2)) + abs(e3) + abs(wr(1) + sum(wd(:,1))) &
         & + abs(wr(2) + sum(wd(:,2)))
    cs_excl_dsigma = cs_excl_dsigma*sx_fac
  end function excl_point

  ! (2,0) (line 1) or (0,2) (line 2) at the three-parton event p6 of the
  ! line: one-loop H+3j with the loop on the line, plus the I operator,
  ! plus K + P, summed over the flavour classes with the PDFs (fB: Born
  ! PDFs of both lines, fE: the line's at its momentum fraction xi3), at
  ! alpha_s = 1. The VBFNLO virtual also contains the vertex correction of
  ! the other line (the (1,1) virtual, B alpha_s/(2 pi) CF (-8 - L^2 - 3 L),
  ! L = ln(mur^2/Qo^2)), which is removed here.
  real(dp) function virt_point(line, p6, fB, fE, xi3, mur, muf, Qo) result(v)
    integer, intent(in) :: line
    real(dp), intent(in) :: p6(0:3,6), fB(-6:6,2), fE(-6:6), xi3, mur, muf, Qo
    real(dp) :: fq(-6:6), fg, iq, ig, lo, born, bmn(0:3,0:3), virt, v20, c11, fgq, vi
    integer :: i1, i2, bflav(6), o, il, io
    real(dp), parameter :: CF = 4.0_dp/3.0_dp, twopi = 6.283185307179586476925286766559005768394_dp
    o = 3 - line
    call nlo2_kp(xi3, p6(:,line), p6(:,3+line), p6(:,6), muf, fq, fg, fgq)
    iq = nlo2_ifin(1, p6(:,line), p6(:,3+line), p6(:,6), mur**2)
    ig = nlo2_ifin(2, p6(:,line), p6(:,3+line), p6(:,6), mur**2)
    lo = log(mur**2/Qo**2)
    c11 = CF*(-8 - lo**2 - 3*lo)/twopi
    v = 0
    ! quark-initiated on the line
    do i1 = 1, ncls
       do i2 = 1, ncls
          if (.not. compatible(cls(i1)%w, cls(i2)%w)) cycle
          bflav = [cls(i1)%a, cls(i2)%a, 25, cls(i1)%b, cls(i2)%b, 0]
          il = merge(i1, i2, line == 1); io = merge(i2, i1, line == 1)
          call memo_born_line(p6, bflav, line, born, bmn)
          if (memo_on) then
             vi = memo_vi(p6, bflav, line, mur**2, born, c11, iq)
          else
             call cs_hjjj_virt_line(p6, bflav, line, mur**2, virt)
             v20 = virt - born*c11
             vi = v20 + born*iq/twopi
          endif
          v = v + (vi*pdfsum(fE, cls(il)) + born/twopi*pdfsum(fq, cls(il))) &
               & *pdfsum(fB(:,o), cls(io))
       enddo
    enddo
    ! gluon-initiated on the line
    do i1 = 1, ngcls
       do i2 = 1, ncls
          if (.not. compatible(gcls(i1)%w, cls(i2)%w)) cycle
          if (line == 1) then
             bflav = [0, cls(i2)%a, 25, gcls(i1)%q, cls(i2)%b, gcls(i1)%qb]
          else
             bflav = [cls(i2)%a, 0, 25, cls(i2)%b, gcls(i1)%q, gcls(i1)%qb]
          endif
          call memo_born_line(p6, bflav, line, born, bmn)
          if (memo_on) then
             vi = memo_vi(p6, bflav, line, mur**2, born, c11, ig)
          else
             call cs_hjjj_virt_line(p6, bflav, line, mur**2, virt)
             v20 = virt - born*c11
             vi = v20 + born*ig/twopi
          endif
          v = v + gcls(i1)%n*(vi*fE(0) + born/twopi*fg)*pdfsum(fB(:,o), cls(i2))
       enddo
    enddo
  end function virt_point

  ! stage-2 set-up: two VBF Born test points (for the grouping of the
  ! real flavour list) and cs_nlo2
  subroutine cs_excl_setup2(verbose)
    logical, intent(in) :: verbose
    real(dp) :: xr(7), pb(0:3,5), xb1, xb2, jacb, sbeams
    integer :: it
    tb_n = 0
    do it = 1, 1000
       xr = modulo(0.1234567_dp*it*[1.0_dp, 1.7_dp, 2.3_dp, 3.1_dp, 3.7_dp, 4.3_dp, 5.9_dp], 1.0_dp)
       call gen_phsp(xr)
       call set_phsp()
       call cs_get_born(pb, xb1, xb2, jacb, sbeams)
       if (jacb == 0 .or. min(Q1_sq, Q2_sq) <= max(Qmin**2, 100.0_dp)) cycle
       if (max(xb1, xb2) > 0.5_dp) cycle
       tb_n = tb_n + 1
       tb_pb(:,:,tb_n) = pb
       tb_xb(:,tb_n) = [xb1, xb2]
       if (tb_n == 2) exit
    enddo
    if (tb_n < 2) stop 'cs_excl_setup2: no test points'
    call nlo2_init(tb_pb, tb_xb, excl_cutoff, verbose)
    call nlo11_init()
  end subroutine cs_excl_setup2

  ! stage 3: map pairs of line states to the class-12 real entries
  ! (tags 1, 2 on the lines' legs 1/4 and 2/5; the two extra partons at
  ! legs 6, 7 with tags 1 and 2 in either order)
  subroutine nlo11_init()
    integer, external :: cs_nreal
    integer :: j, flav(7), tags(7), l1, s1, s2, nf, nmiss
    logical :: sw1, sw2
    e3_entry = 0; e3_leg1 = 0; e3_fac = 1
    do j = 1, cs_nreal()
       call cs_get_real(j, flav, tags)
       if (tags(1) /= 1 .or. tags(2) /= 2 .or. tags(4) /= 1 .or. tags(5) /= 2) cycle
       if (tags(6) == 1 .and. tags(7) == 2) then
          l1 = 6
       elseif (tags(6) == 2 .and. tags(7) == 1) then
          l1 = 7
       else
          cycle
       endif
       s1 = line_state(flav(1), flav(4), flav(l1), sw1)
       s2 = line_state(flav(2), flav(5), flav(13 - l1), sw2)
       if (s1 == 0 .or. s2 == 0) cycle
       if (e3_entry(s1,s2) /= 0) cycle
       e3_entry(s1,s2) = j; e3_leg1(s1,s2) = l1
       e3_swap(:,s1,s2) = [sw1, sw2]
    enddo
    ! normalisation: the tagged entries give the full matrix element for
    ! the two lines' extra partons (no symmetry factor for two gluons on
    ! different lines; checked against the H+3j times the other line's
    ! dipoles in its singular limits)
    e3_fac = 1
    nf = 0; nmiss = 0
    do s1 = 1, ncls + ngcls
       do s2 = 1, ncls + ngcls
          if (.not. compatible(state_w(s1), state_w(s2))) cycle
          if (e3_entry(s1,s2) > 0) then
             nf = nf + 1
          else
             nmiss = nmiss + 1
             write(6,'(a,2i3)') ' nlo11_init: no class-12 entry for line states', s1, s2
          endif
       enddo
    enddo
    write(6,'(a,i4,a,i4)') ' nlo11_init: class-12 entries for pairs of line states:', nf, ', missing:', nmiss
  end subroutine nlo11_init

  ! state of a line from its incoming, outgoing and extra flavour
  integer function line_state(a, b, e, swap) result(st)
    integer, intent(in) :: a, b, e
    logical, intent(out) :: swap
    integer :: i
    st = 0; swap = .false.
    if (e == 0 .and. a /= 0) then
       do i = 1, ncls
          if (cls(i)%a == a .and. cls(i)%b == b) st = i
       enddo
    elseif (a == 0) then
       do i = 1, ngcls
          if (gcls(i)%q == b .and. gcls(i)%qb == e) st = ncls + i
       enddo
       if (st == 0) then
          do i = 1, ngcls
             if (gcls(i)%q == e .and. gcls(i)%qb == b) then
                st = ncls + i; swap = .true.
             endif
          enddo
       endif
    endif
  end function line_state

  integer function state_w(st)
    integer, intent(in) :: st
    if (st <= ncls) then
       state_w = cls(st)%w
    else
       state_w = gcls(st - ncls)%w
    endif
  end function state_w

  ! PDF factor of a radiating line in state st at its momentum fraction
  real(dp) function state_pdf(st, f)
    integer, intent(in) :: st
    real(dp), intent(in) :: f(-6:6)
    if (st <= ncls) then
       state_pdf = pdfsum(f, cls(st))
    else
       state_pdf = gcls(st - ncls)%n*f(0)
    endif
  end function state_pdf

  ! H+3j flavours (cs_glue order) for line l in state st and the other
  ! line at Born level in class ib
  function h3_flav(l, st, ib) result(bflav)
    integer, intent(in) :: l, st, ib
    integer :: bflav(6)
    bflav = [cls(ib)%a, cls(ib)%a, 25, cls(ib)%b, cls(ib)%b, 0]
    if (st <= ncls) then
       bflav(l) = cls(st)%a; bflav(3+l) = cls(st)%b
    else
       bflav(l) = 0; bflav(3+l) = gcls(st - ncls)%q; bflav(6) = gcls(st - ncls)%qb
    endif
  end function h3_flav

  ! stage 3, events E1 (l = 1) and E2 (l = 2): line l radiated (its
  ! three-parton point p6l, PDFs fEl at its momentum fraction), the other
  ! line o at its Born with its O(alpha_s) pieces: (V + I) f + (K + P) (x) f
  ! at the Born, minus its dipoles at its own three-parton point (if
  ! oko; p6o, fEo, wrado). Weight in pb, without the VEGAS weight.
  real(dp) function nlo11_e12(l, pb, xb, fB, p6l, fEl, wradl, oko, p6o, fEo, wrado, mur, muf, Qs, common) result(w)
    use cs_dipoles, only: h_if_qg, h_fi_qg
    integer, intent(in) :: l
    real(dp), intent(in) :: pb(0:3,5), xb(2), fB(-6:6,2), p6l(0:3,6), fEl(-6:6), wradl, p6o(0:3,6), fEo(-6:6), wrado
    real(dp), intent(in) :: mur(2), muf(2), Qs(2), common
    logical, intent(in) :: oko
    real(dp), external :: hoppetAlphaS
    real(dp), parameter :: CF = 4.0_dp/3.0_dp, TR = 0.5_dp
    real(dp) :: twopi, lo, vi, fq(-6:6), me(8), m1, m2, pdfl, wb, wd, x, u, pa(0:3), pk(0:3), pt(0:3), kq, kg
    integer :: o, st, ib, ig, b5(5), jq, jqb
    twopi = 8*atan(1.0_dp)
    o = 3 - l
    lo = log(mur(o)**2/Qs(o)**2)
    vi = CF*(-8 - lo**2 - 3*lo)/twopi + nlo2_ifin_dis(pb(:,o), pb(:,3+o), mur(o)**2)/twopi
    call nlo2_kp_dis(xb(o), pb(:,o), pb(:,3+o), muf(o), fq)
    if (oko) then
       pa = p6o(:,o); pk = p6o(:,3+o); pt = p6o(:,6)
       x = 1 - mdot(pt, pk)/mdot(pa, pt + pk)
       u = mdot(pa, pt)/mdot(pa, pt + pk)
       kq = CF*4*twopi*(h_if_qg(x, u, 1.0_dp)/(2*mdot(pa, pt)*x) + h_fi_qg(1 - u, x, 1.0_dp)/(2*mdot(pk, pt)*x))
       kg = TR*(1 - 2*x*(1 - x))*4*twopi/x
    endif
    wb = 0; wd = 0
    do st = 1, ncls + ngcls
       pdfl = state_pdf(st, fEl)
       if (pdfl == 0) cycle
       me = 0
       do ib = 1, ncls
          if (.not. compatible(state_w(st), cls(ib)%w)) cycle
          call memo_lines(p6l, h3_flav(l, st, ib), m1, m2)
          me(ib) = merge(m1, m2, l == 1)
          wb = wb + me(ib)*pdfl*(vi*pdfsum(fB(:,o), cls(ib)) + pdfsum(fq, cls(ib))/twopi)
          if (oko) wd = wd - kq*me(ib)*pdfl*pdfsum(fEo, cls(ib))*wrado
       enddo
       if (.not. oko) cycle
       ! the other line gluon-initiated: its two IF dipoles, Born -qbar -> q
       ! (antiquark emitted, pt) and -q -> qbar (quark emitted, pk)
       do ig = 1, ngcls
          if (.not. compatible(state_w(st), gcls(ig)%w)) cycle
          jq = 0; jqb = 0
          do ib = 1, ncls
             if (cls(ib)%a == -gcls(ig)%qb .and. cls(ib)%b == gcls(ig)%q) jq = ib
             if (cls(ib)%a == -gcls(ig)%q .and. cls(ib)%b == gcls(ig)%qb) jqb = ib
          enddo
          if (jq == 0 .or. jqb == 0) stop 'nlo11_e12: Born classes of a gluon class not found'
          wd = wd - (kg/(2*mdot(pa, pt))*me(jq) + kg/(2*mdot(pa, pk))*me(jqb))*pdfl*gcls(ig)%n*fEo(0)*wrado
       enddo
    enddo
    w = (wb + wd)*common*wradl*hoppetAlphaS(mur(l))*hoppetAlphaS(mur(o))
  end function nlo11_e12

  ! stage 3, event E3: both lines radiated (the class-12 H+4j at the
  ! two three-parton points); p7: proVBFH order (1, 2 incoming, 3 Higgs,
  ! 4, 5 outgoing quarks of lines 1, 2, 6, 7 the extra partons of lines
  ! 1, 2). Weight in pb, without the VEGAS weight.
  real(dp) function nlo11_e3(pb, p6, fE, wrad, mur, common, p7) result(w)
    real(dp), intent(in) :: pb(0:3,5), p6(0:3,6,2), fE(-6:6,2), wrad(2), mur(2), common
    real(dp), intent(out) :: p7(0:3,7)
    real(dp), external :: hoppetAlphaS
    real(dp) :: pe(0:3,7), me, p1, p2
    integer :: s1, s2, j, l1
    p7(:,1) = p6(:,1,1); p7(:,2) = p6(:,2,2); p7(:,3) = pb(:,3)
    p7(:,4) = p6(:,4,1); p7(:,5) = p6(:,5,2); p7(:,6) = p6(:,6,1); p7(:,7) = p6(:,6,2)
    w = 0
    do s1 = 1, ncls + ngcls
       p1 = state_pdf(s1, fE(:,1))
       if (p1 == 0) cycle
       do s2 = 1, ncls + ngcls
          j = e3_entry(s1,s2)
          if (j == 0) cycle
          p2 = state_pdf(s2, fE(:,2))
          if (p2 == 0) cycle
          call e3_momenta(s1, s2, p7, pe)
          call memo_real_me(j, pe, me)
          w = w + e3_fac(s1,s2)*me*p1*p2
       enddo
    enddo
    w = w*common*wrad(1)*wrad(2)*hoppetAlphaS(mur(1))*hoppetAlphaS(mur(2))
  end function nlo11_e3

  ! momenta of the class-12 entry of the line states (s1, s2) from p7
  ! (proVBFH order: 4, 5 outgoing quarks, 6, 7 extra partons of lines 1, 2)
  subroutine e3_momenta(s1, s2, p7, pe)
    integer, intent(in) :: s1, s2
    real(dp), intent(in) :: p7(0:3,7)
    real(dp), intent(out) :: pe(0:3,7)
    integer :: l1
    l1 = e3_leg1(s1,s2)
    pe(:,1:3) = p7(:,1:3)
    if (e3_swap(1,s1,s2)) then
       pe(:,4) = p7(:,6); pe(:,l1) = p7(:,4)
    else
       pe(:,4) = p7(:,4); pe(:,l1) = p7(:,6)
    endif
    if (e3_swap(2,s1,s2)) then
       pe(:,5) = p7(:,7); pe(:,13-l1) = p7(:,5)
    else
       pe(:,5) = p7(:,5); pe(:,13-l1) = p7(:,7)
    endif
  end subroutine e3_momenta

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

  ! cs_scales: number of scale points (1, 3 or 7)
  subroutine cs_scales_setup(n)
    integer, intent(in) :: n
    if (n /= 1 .and. n /= 3 .and. n /= 7) stop 'cs_scales must be 1, 3 or 7'
    excl_nscale = n
    memo_on = n > 1
    if (n > 1) then
       if (scale_choice /= 2 .and. scale_choice /= 3) &
            & stop 'cs_scales: on-the-fly variations need scale_choice 2 or 3 (runningscales 1)'
    endif
  end subroutine cs_scales_setup

  ! the matrix-element routines of the event weights through the per-point
  ! cache (memo_on, scale variations), else directly
  subroutine memo_lines(p6, bflav, m1, m2)
    real(dp), intent(in) :: p6(0:3,6)
    integer, intent(in) :: bflav(6)
    real(dp), intent(out) :: m1, m2
    integer :: key(memo_nkey), slot
    real(dp) :: val(2)
    if (.not. memo_on) then
       call cs_hjjj_lines(p6, bflav, m1, m2)
       return
    endif
    key = 0; key(1) = 1; key(2:7) = bflav
    if (memo_get(key, reshape(p6, [24]), 24, val, 2, slot)) then
       m1 = val(1); m2 = val(2)
    else
       call cs_hjjj_lines(p6, bflav, m1, m2)
       call memo_put(slot, key, reshape(p6, [24]), 24, [m1, m2], 2)
    endif
  end subroutine memo_lines

  subroutine memo_born_line(p6, bflav, line, born, bmn)
    real(dp), intent(in) :: p6(0:3,6)
    integer, intent(in) :: bflav(6), line
    real(dp), intent(out) :: born, bmn(0:3,0:3)
    integer :: key(memo_nkey), slot
    real(dp) :: val(17)
    if (.not. memo_on) then
       call cs_hjjj_born_line(p6, bflav, line, born, bmn)
       return
    endif
    key = 0; key(1) = 2; key(2:7) = bflav; key(8) = line
    if (memo_get(key, reshape(p6, [24]), 24, val, 17, slot)) then
       born = val(1); bmn = reshape(val(2:17), [4, 4])
    else
       call cs_hjjj_born_line(p6, bflav, line, born, bmn)
       call memo_put(slot, key, reshape(p6, [24]), 24, [born, reshape(bmn, [16])], 17)
    endif
  end subroutine memo_born_line

  ! the virtual at mu_R^2 = mur2 (cached per value of mur2)
  subroutine memo_virt_line(p6, bflav, line, mur2, virt)
    real(dp), intent(in) :: p6(0:3,6), mur2
    integer, intent(in) :: bflav(6), line
    real(dp), intent(out) :: virt
    integer :: key(memo_nkey), slot
    real(dp) :: val(1), mom(25)
    if (.not. memo_on) then
       call cs_hjjj_virt_line(p6, bflav, line, mur2, virt)
       return
    endif
    key = 0; key(1) = 3; key(2:7) = bflav; key(8) = line
    mom(1:24) = reshape(p6, [24]); mom(25) = mur2
    if (memo_get(key, mom, 25, val, 1, slot)) then
       virt = val(1)
    else
       call cs_hjjj_virt_line(p6, bflav, line, mur2, virt)
       call memo_put(slot, key, mom, 25, [virt], 1)
    endif
  end subroutine memo_virt_line

  ! V + I of the radiating line, (virtual - born c11) + born i/(2 pi), at
  ! mu_R^2 = mur2 (c11: the other line's vertex term, i: the line's I
  ! operator, both at mur2), through the per-point cache: evaluated at the
  ! first mu_R of the point, then shifted by its renormalisation-group
  ! dependence, born b0/(4 pi) ln(mur2/mur2_ref), b0 = 11/3 CA - 4/3 TR nf
  ! (the infrared logarithms of V and I cancel; the H+3j tree is O(g_s))
  real(dp) function memo_vi(p6, bflav, line, mur2, born, c11, ii) result(vi)
    real(dp), intent(in) :: p6(0:3,6), mur2, born, c11, ii
    integer, intent(in) :: bflav(6), line
    real(dp), parameter :: CA = 3, TR = 0.5_dp, twopi = 6.283185307179586476925286766559005768394_dp
    real(dp), external :: cs_nlight
    integer :: key(memo_nkey), slot
    real(dp) :: val(2), virt, b0_4pi
    key = 0; key(1) = 6; key(2:7) = bflav; key(8) = line
    if (memo_get(key, reshape(p6, [24]), 24, val, 2, slot)) then
       b0_4pi = (11.0_dp/3*CA - 4.0_dp/3*TR*cs_nlight())/(2*twopi)
       vi = val(1) + born*b0_4pi*log(mur2/val(2))
    else
       call cs_hjjj_virt_line(p6, bflav, line, mur2, virt)
       vi = (virt - born*c11) + born*ii/twopi
       call memo_put(slot, key, reshape(p6, [24]), 24, [vi, mur2], 2)
    endif
  end function memo_vi

  subroutine memo_real_me(j, pe, me)
    integer, intent(in) :: j
    real(dp), intent(in) :: pe(0:3,7)
    real(dp), intent(out) :: me
    integer :: key(memo_nkey), slot
    real(dp) :: val(1)
    if (.not. memo_on) then
       call cs_real_me(j, pe, me)
       return
    endif
    key = 0; key(1) = 4; key(2) = j
    if (memo_get(key, reshape(pe, [28]), 28, val, 1, slot)) then
       me = val(1)
    else
       call cs_real_me(j, pe, me)
       call memo_put(slot, key, reshape(pe, [28]), 28, [me], 1)
    endif
  end subroutine memo_real_me

  ! cs_analysis with one weight per scale point (wgt(1): the central one);
  ! a single scale point goes through cs_analysis unchanged
  subroutine cs_analysis_m(n, p, wgt)
    integer, intent(in) :: n
    real(dp), intent(in) :: p(0:3,n), wgt(:)
    if (size(wgt) == 1) then
       call cs_analysis(n, p, wgt(1))
       return
    endif
    call cs_set_weights(size(wgt), wgt)
    call cs_fill_phep(n, p)
    call user_analysis(wgt(1))
  end subroutine cs_analysis_m

  subroutine cs_analysis(n, p, wgt)
    integer, intent(in) :: n
    real(dp), intent(in) :: p(0:3,n), wgt
    call cs_fill_phep(n, p)
    call user_analysis(wgt)
  end subroutine cs_analysis

end module cs_exclusive
