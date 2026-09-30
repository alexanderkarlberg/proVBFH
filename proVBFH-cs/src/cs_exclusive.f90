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
! Stage 2 (excl_order = 2, docs/DESIGN.md): in addition the (2,0) and
! (0,2) contributions (cs_nlo2), i.e. the NLO correction of the line
! that radiated. Per line, at the same three-parton event, the one-loop
! H+3j (the loop on the radiating line) with the I operator and the K + P
! terms; and a four-parton event (seven more random numbers, 14:20, not
! adapted by VEGAS) with its six dipole counterevents. All have their
! counterevent at the VBF Born point, which receives minus the sum.
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
  use cs_nlo2, dbg_mvar => mvar, dbg_mz => mz, dbg_w4 => cur_w4
  implicit none
  private
  public :: cs_excl_dsigma, cs_excl_setup, excl_fill, excl_npow, excl_cutoff, excl_stats, excl_flavcheck
  public :: excl_order, excl_only2, excl_dump2, dump2_min, cs_excl_replay, cs_excl_setup2, cs_excl_testlimits, cs_excl_testvirt

  ! 1: (1,0) + (0,1) (NLO); 2: also (2,0) + (0,2)
  integer, save :: excl_order = 1
  ! with excl_order = 2: keep only the (2,0) + (0,2) weights (tests of the
  ! O(alpha_s^2) part alone: cutoff independence, variance)
  logical, save :: excl_only2 = .false.
  ! with excl_order = 2: write points whose signed contribution to the
  ! first histogram (sigma with the inclusive cuts) exceeds dump2_min pb
  ! to unit 79 (file cs_dump2.dat)
  logical, save :: excl_dump2 = .false.
  ! replay (cs_replay): print the stage-2 real and dipoles of each point
  logical, public, save :: excl_verbose = .false.
  real(dp), save :: dump2_min = 1d-3
  ! with excl_order = 2: points whose integrand times VEGAS weight exceeds
  ! spike_min (default 10 pb, about 1e4 times a typical point) are written
  ! to unit 80 (cs_spikes.dat) with their random numbers, for cs_replay
  real(dp), public, save :: spike_min = 10
  ! estimate (debug, notes/2026-09-29-cs-p2b-stage2): weight the (1,0) and
  ! (0,1) events by alpha_s/(2 pi) Delta, Delta = F(mu_e) - F(Q) -
  ! b0 ln(mu_e^2/Q^2), the change of the one-loop H+3j when VBFNLO's F
  ! function (ffunc) is evaluated at mu_e instead of the line's Q;
  ! excl_estimu 1: mu_e = mu0(pt,H), 2: mu_e = M_H.
  ! excl_estimate 2: the coefficient of ln(1/cutoff) of the NC graphs in
  ! which an incoming quark emits a gluon that fuses with the Z into a
  ! pair, if their initial-state collinear singularity is not subtracted
  ! (no FKS region in proVBFH's POWHEG part) and with the pair type swapped
  ! (compreal_hqqqq_new.f): the NC gluon-initiated (1,0) events with the
  ! gluon density replaced by alpha_s/(2 pi) [P_gq (x) sum_q f_q] and the
  ! Born of the other pair type; excl_estimate 3: the same with the
  ! correct pair type.
  ! excl_estimate 4: the effect of nf = 4 in VBFNLO's ffunc (nf = 5 elsewhere)
  ! on the one-loop H+3j: alpha_s/(2 pi) B (1/6) [ln(mu^2/|t|) + ln(mu^2/|u|)],
  ! t, u the gluon's invariants with the two quarks of its line, at mu_e
  integer, public, save :: excl_estimate = 0, excl_estimu = 1
  real(dp), save :: tb_pb(0:3,5,2), tb_xb(2,2)
  integer, save :: tb_n = 0

  logical, save :: excl_fill = .false.
  integer, save :: excl_npow = 0     ! 0: logarithmic sampling
  real(dp), save :: excl_cutoff = 1d-6
  ! counters: points, line radiations rejected by the cutoff, NaN, points
  ! where neither event nor Born passes the cuts
  integer(8), save :: excl_stats(4) = 0
  ! skip points where neither the events nor the Born pass the VBF cuts
  ! (as proVBFH's phspcuts; exact for all VBF-cut histograms, and the
  ! exclusive part does not contribute to histograms without cuts: the
  ! Higgs momentum is unchanged by the projection); loose cuts (tagging
  ! jet pt and mjj only) while building the grid
  logical, public, save :: excl_phspcuts = .true.
  ! debug: write per-point contributions to the VBF-cut cross section
  ! (unit 77, file cs_dump.dat) when excl_dump is set
  logical, public, save :: excl_dump = .false.
  real(dp), save :: dxp(2), dz(2)
  integer, save :: ndbg = 0
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

  ! VEGAS integrand. The histograms (pwhg_bookhist-multi) are normalised
  ! by the number of pwhgaccumup calls, so every point must be
  ! accumulated, also points that return early with zero weight.
  double precision function cs_excl_dsigma(xrand, vegas_weight)
    real(dp), intent(in) :: xrand(*), vegas_weight
    cs_excl_dsigma = excl_point(xrand, vegas_weight)
    if (excl_fill) call pwhgaccumup
  end function cs_excl_dsigma

  ! replay: the points listed in cs_replay.dat (20 random numbers each)
  subroutine cs_excl_replay()
    real(dp) :: xr(20), r
    integer :: ios, n
    open(81, file='cs_replay.dat', status='old')
    excl_verbose = .true.
    n = 0
    do
       read(81, *, iostat=ios) xr
       if (ios /= 0) exit
       n = n + 1
       write(6,'(a,i4)') ' ===== replay point', n
       r = excl_point(xr, 1.0_dp)
       write(6,'(a,es12.4)') '   integrand (vegas weight 1):', r
    enddo
    close(81)
    excl_verbose = .false.
  end subroutine cs_excl_replay

  double precision function excl_point(xrand, vegas_weight) result(cs_excl_dsigma)
    real(dp), intent(in) :: xrand(*), vegas_weight
    real(dp) :: pb(0:3,5), xb1, xb2, jacb, sbeams, common, ptH
    real(dp) :: Q1, Q2, mur(2), muf(2), as(2), fB(-6:6,2), fE(-6:6,2)
    real(dp) :: pin(0:3), a(0:3), b(0:3), xp, z, wrad(2), p6(0:3,6,2), w(2)
    real(dp) :: m1, m2, me
    real(dp) :: xi3(2), wv(2), wr(2), wd(6,2), pr(0:3,7,2), pd(0:3,6,6,2), wg(2), mue, dlq, dlg, fsing
    integer :: ib
    logical :: ok(2), need(2), pass(2), passB
    logical :: okr(2), dok(6,2), needr(2), passr(2), passd(6,2)
    real(dp) :: prt(0:3,7), pdt(0:3,6,6), csig, dmv(6,2), dmz(6,2), dw4(2)
    logical :: dokt(6), okt
    integer :: m
    logical, external :: cs_passes
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
       dxp(line) = xp; dz(line) = z
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

    ! cut decisions before the matrix elements: a line's weight is needed
    ! only if its event or the Born event passes
    need = ok
    needr = okr
    if (excl_phspcuts) then
       passB = cs_passes(5, pb, .not. excl_fill)
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
       if (.not. any(need) .and. .not. any(needr)) then
          excl_stats(4) = excl_stats(4) + 1
          return
       endif
    endif

    w = 0; wg = 0
    do line = 1, 2
       if (.not. need(line)) cycle
       ! quark-initiated on the radiating line
       do i1 = 1, merge(0, ncls, excl_estimate == 2 .or. excl_estimate == 3)
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
       if (excl_estimate == 2 .or. excl_estimate == 3) then
          ! only the NC gluon-initiated classes, with the collinear density
          w(line) = 0
          fsing = pgq_conv(xi3(line), muf(line))*as(line)/(8*atan(1.0_dp))
          do i1 = 1, 2
             do i2 = 1, ncls
                if (.not. compatible(gcls(i1)%w, cls(i2)%w)) cycle
                ! the Born of the pair type the swapped index selects
                ib = i1
                if (excl_estimate == 2) ib = 3 - i1
                if (line == 1) then
                   bflav = [0, cls(i2)%a, 25, gcls(ib)%q, cls(i2)%b, gcls(ib)%qb]
                   call cs_hjjj_lines(p6(:,:,line), bflav, m1, m2)
                   w(1) = w(1) + m1*gcls(i1)%n*fsing*pdfsum(fB(:,2), cls(i2))
                else
                   bflav = [cls(i2)%a, 0, 25, cls(i2)%b, gcls(ib)%q, gcls(ib)%qb]
                   call cs_hjjj_lines(p6(:,:,line), bflav, m1, m2)
                   w(2) = w(2) + m2*gcls(i1)%n*pdfsum(fB(:,1), cls(i2))*fsing
                endif
             enddo
          enddo
          w(line) = w(line)*common*wrad(line)*as(line)
          cycle
       endif
       do i1 = 1, ngcls
          do i2 = 1, ncls
             if (.not. compatible(gcls(i1)%w, cls(i2)%w)) cycle
             if (line == 1) then
                bflav = [0, cls(i2)%a, 25, gcls(i1)%q, cls(i2)%b, gcls(i1)%qb]
                call cs_hjjj_lines(p6(:,:,line), bflav, m1, m2)
                w(1) = w(1) + m1*gcls(i1)%n*fE(0,1)*pdfsum(fB(:,2), cls(i2))
                wg(1) = wg(1) + m1*gcls(i1)%n*fE(0,1)*pdfsum(fB(:,2), cls(i2))
             else
                bflav = [cls(i2)%a, 0, 25, cls(i2)%b, gcls(i1)%q, gcls(i1)%qb]
                call cs_hjjj_lines(p6(:,:,line), bflav, m1, m2)
                w(2) = w(2) + m2*gcls(i1)%n*pdfsum(fB(:,1), cls(i2))*fE(0,2)
                wg(2) = wg(2) + m2*gcls(i1)%n*pdfsum(fB(:,1), cls(i2))*fE(0,2)
             endif
          enddo
       enddo
       w(line) = w(line)*common*wrad(line)*as(line)
       wg(line) = wg(line)*common*wrad(line)*as(line)
       if (excl_estimate == 1 .or. excl_estimate == 4) then
          ! the line's partons: in p6(:,line), out p6(:,3+line), extra p6(:,6);
          ! quark-initiated: s = in.out, t, u = gluon with in, out;
          ! gluon-initiated: s = the two outgoing quarks, t, u = with the gluon
          if (excl_estimu == 1) then
             mue = cs_mu(1.0_dp, line, Q1, Q2, ptH, .true.)
          else
             mue = mh
          endif
          if (excl_estimate == 1) then
             dlq = ffunc4(2*mdot(p6(:,line,line), p6(:,3+line,line)), 2*mdot(p6(:,line,line), p6(:,6,line)), &
                  & 2*mdot(p6(:,3+line,line), p6(:,6,line)), mue, merge(Q1, Q2, line == 1))
             dlg = ffunc4(2*mdot(p6(:,3+line,line), p6(:,6,line)), 2*mdot(p6(:,line,line), p6(:,6,line)), &
                  & 2*mdot(p6(:,line,line), p6(:,3+line,line)), mue, merge(Q1, Q2, line == 1))
          else
             ! quark-initiated: gluon at 6 with in (line) and out (3+line);
             ! gluon-initiated: gluon incoming (line) with Q (3+line), Qbar (6)
             dlq = (log(mue**2/abs(2*mdot(p6(:,line,line), p6(:,6,line)))) &
                  & + log(mue**2/abs(2*mdot(p6(:,3+line,line), p6(:,6,line)))))/6
             dlg = (log(mue**2/abs(2*mdot(p6(:,line,line), p6(:,3+line,line)))) &
                  & + log(mue**2/abs(2*mdot(p6(:,line,line), p6(:,6,line)))))/6
          endif
          w(line) = (w(line) - wg(line))*dlq*as(line)/(8*atan(1.0_dp)) + wg(line)*dlg*as(line)/(8*atan(1.0_dp))
       endif
    enddo

    ! (2,0) and (0,2)
    wv = 0; wr = 0; wd = 0; dmv = 0; dmz = 0; dw4 = 0
    if (excl_order >= 2) then
       do line = 1, 2
          if (need(line)) wv(line) = virt_point(line, p6(:,:,line), fB, fE(:,line), xi3(line), &
               & mur(line), muf(line), merge(Q2, Q1, line == 1))*common*wrad(line)*as(line)**2
          if (needr(line)) then
             ! nlo2_real_me works on the state of the last nlo2_real_kin call
             ! (line 2 after the loop above): set up this line's point again
             call nlo2_real_kin(line, pb, [x1, x2], xrand(14:20), excl_cutoff, prt, pdt, dokt, okt)
             if (.not. okt .or. any(dokt .neqv. dok(:,line))) stop 'cs_exclusive: nlo2_real_kin not reproducible'
             dmv(:,line) = dbg_mvar; dmz(:,line) = dbg_mz; dw4(line) = dbg_w4
             call nlo2_real_me(muf(line), fB(:,3-line), wr(line), wd(:,line))
             if (excl_verbose) call nlo2_debug_point(muf(line), fB(:,3-line), 6)
             wr(line) = wr(line)*common*as(line)**2
             wd(:,line) = wd(:,line)*common*as(line)**2
          endif
       enddo
    endif
    if (excl_only2) w = 0
    if (excl_flavcheck > 0) then
       excl_flavcheck = excl_flavcheck - 1
       call flavour_check(p6, need, fB, fE, common*wrad*as, w)
    endif
    if (any(w /= w) .or. any(wv /= wv) .or. any(wr /= wr) .or. any(wd /= wd)) then
       excl_stats(3) = excl_stats(3) + 1
       return
    endif

    if (excl_dump .and. excl_phspcuts) then
       do line = 1, 2
          if (ok(line) .and. (pass(line) .neqv. passB) .and. 1 - dxp(line) < 1d-3 .and. ndbg < 5) then
             ndbg = ndbg + 1
             write(6,'(a,i2,a,2es12.4,a,2l2)') ' DBG line', line, ' 1-xp, z = ', 1-dxp(line), dz(line), &
                  & '  pass evt, Born:', pass(line), passB
             call cs_debug_jets(5, pb)
             call cs_debug_jets(6, p6(:,:,line))
          endif
       enddo
    endif
    if (excl_dump .and. excl_phspcuts) then
       write(77,'(2l2,l3,4es14.6,2es14.6,2es12.4)') merge(pass, [.false.,.false.], ok), passB, &
            & w, w(1)*(merge(1,0,pass(1).and.ok(1)) - merge(1,0,passB)), &
            & w(2)*(merge(1,0,pass(2).and.ok(2)) - merge(1,0,passB)), dxp, dz, Q1, Q2
    endif
    if (excl_dump2 .and. excl_phspcuts .and. excl_fill) then
       csig = 0
       if (passB) csig = -(sum(w) + sum(wv) + sum(wr) + sum(wd))
       do line = 1, 2
          if (need(line) .and. pass(line)) csig = csig + w(line) + wv(line)
          if (needr(line)) then
             if (passr(line)) csig = csig + wr(line)
             do m = 1, 6
                if (dok(m,line) .and. passd(m,line)) csig = csig + wd(m,line)
             enddo
          endif
       enddo
       csig = csig*vegas_weight
       if (abs(csig) > dump2_min) write(79,'(es12.4,l3,2l2,2l2,4es11.3,2(2x,6l1),6es11.3,2(2x,12es10.2),2es10.2,2x,20es24.16)') &
            & csig, passB, merge(pass, [.false.,.false.], need), merge(passr, [.false.,.false.], needr), &
            & 1-dxp, dz, dok, (w + wv)*vegas_weight, wr*vegas_weight, sum(wd,1)*vegas_weight, &
            & (dmv(:,line), dmz(:,line), line=1,2), dw4*vegas_weight, xrand(1:20)
    endif
    if (excl_fill) then
       do line = 1, 2
          if (need(line)) call cs_analysis(6, p6(:,:,line), (w(line) + wv(line))*vegas_ncall*vegas_weight)
          if (needr(line)) then
             call cs_analysis(7, pr(:,:,line), wr(line)*vegas_ncall*vegas_weight)
             do m = 1, 6
                if (dok(m,line)) call cs_analysis(6, pd(:,:,m,line), wd(m,line)*vegas_ncall*vegas_weight)
             enddo
          endif
       enddo
       call cs_analysis(5, pb, -(sum(w) + sum(wv) + sum(wr) + sum(wd))*vegas_ncall*vegas_weight)
    endif
    cs_excl_dsigma = abs(w(1) + wv(1)) + abs(w(2) + wv(2)) + abs(wr(1) + sum(wd(:,1))) &
         & + abs(wr(2) + sum(wd(:,2)))
    if (excl_order >= 2 .and. cs_excl_dsigma*vegas_weight > spike_min .and. .not. excl_verbose) then
       write(80,'(es12.4,6es11.3,2x,20es24.16)') cs_excl_dsigma*vegas_weight, (w + wv)*vegas_weight, &
            & wr*vegas_weight, sum(wd,1)*vegas_weight, xrand(1:20)
       flush(80)
    endif
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
    real(dp) :: fq(-6:6), fg, iq, ig, lo, born, bmn(0:3,0:3), virt, v20, c11
    integer :: i1, i2, bflav(6), o, il, io
    real(dp), parameter :: CF = 4.0_dp/3.0_dp, twopi = 6.283185307179586476925286766559005768394_dp
    o = 3 - line
    call nlo2_kp(xi3, p6(:,line), p6(:,3+line), p6(:,6), muf, fq, fg)
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
          call cs_hjjj_born_line(p6, bflav, line, born, bmn)
          call cs_hjjj_virt_line(p6, bflav, line, mur**2, virt)
          v20 = virt - born*c11
          v = v + ((v20 + born*iq/twopi)*pdfsum(fE, cls(il)) + born/twopi*pdfsum(fq, cls(il))) &
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
          call cs_hjjj_born_line(p6, bflav, line, born, bmn)
          call cs_hjjj_virt_line(p6, bflav, line, mur**2, virt)
          v20 = virt - born*c11
          v = v + gcls(i1)%n*((v20 + born*ig/twopi)*fE(0) + born/twopi*fg)*pdfsum(fB(:,o), cls(i2))
       enddo
    enddo
  end function virt_point

  ! [P_gq (x) sum_q f_q](xi) = int_xi^1 dz/z CF (1 + (1-z)^2)/z sum_q f_q(xi/z),
  ! number densities, 5 flavours of quarks and antiquarks (Gauss-Legendre in
  ! ln z)
  real(dp) function pgq_conv(xi, muf) result(c)
    real(dp), intent(in) :: xi, muf
    integer, parameter :: n = 48
    real(dp), save :: gx(n), gw(n)
    logical, save :: ini = .true.
    real(dp) :: z, f(-6:6), wz
    integer :: i, j, it
    real(dp) :: t, t1, p1, p2, p3, pp
    if (ini) then
       do i = 1, (n + 1)/2
          t = cos(4*atan(1.0_dp)*(i - 0.25_dp)/(n + 0.5_dp))
          do it = 1, 100
             p1 = 1; p2 = 0
             do j = 1, n
                p3 = p2; p2 = p1
                p1 = ((2*j - 1)*t*p2 - (j - 1)*p3)/j
             enddo
             pp = n*(t*p1 - p2)/(t*t - 1)
             t1 = t; t = t1 - p1/pp
             if (abs(t - t1) < 1d-15) exit
          enddo
          gx(i) = (1 - t)/2; gx(n+1-i) = (1 + t)/2
          gw(i) = 1/((1 - t*t)*pp*pp); gw(n+1-i) = gw(i)
       enddo
       ini = .false.
    endif
    c = 0
    do i = 1, n
       z = xi**gx(i)                  ! ln z uniform in [ln xi, 0]
       wz = gw(i)*z*log(1/xi)
       call hoppetEval(xi/z, muf, f)
       f = f/(xi/z)                   ! number densities at xi/z
       c = c + wz/z*4.0_dp/3*(1 + (1 - z)**2)/z*(sum(f(-5:-1)) + sum(f(1:5)))
    enddo
  end function pgq_conv

  ! F(mu_e) - F(Q) - b0 ln(mu_e^2/Q^2) with VBFNLO's ffunc (qqhqqj-virt.f,
  ! nf = 4 there) and b0 = 11/6 CA - 2/3 TR nf with nf = 5 (as its log terms)
  real(dp) function ffunc4(s, t, u, mue, Q) result(d)
    real(dp), intent(in) :: s, t, u, mue, Q
    real(dp), parameter :: ca = 3, cf = 4.0_dp/3, tr = 0.5_dp, b0 = 11.0_dp/6*3 - 2.0_dp/3*0.5_dp*5
    d = f(mue**2) - f(Q**2) - b0*log(mue**2/Q**2)
  contains
    real(dp) function f(mu2)
      real(dp), intent(in) :: mu2
      real(dp) :: ls, lt, lu
      ls = log(abs(s/mu2)); lt = log(abs(t/mu2)); lu = log(abs(u/mu2))
      f = 0.5_dp*ca*(lu**2 + lt**2) - 0.5_dp*(ca - 2*cf)*ls**2 + 1.5_dp*(ca - 2*cf)*ls &
           & + (tr*4/3.0_dp - 5*ca/3.0_dp)*(lu + lt)
    end function f
  end function ffunc4

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
  end subroutine cs_excl_setup2

  ! limit tests of the real matrix elements against their dipoles
  subroutine cs_excl_testlimits(mode)
    integer, intent(in) :: mode
    if (mode == 3) then
       call debug_is3()
       call nlo2_debug_is4(tb_pb(:,:,1), tb_xb(:,1), 1)
       call nlo2_debug_is4(tb_pb(:,:,1), tb_xb(:,1), 9)
       call nlo2_debug_is4(tb_pb(:,:,1), tb_xb(:,1), 13)
       call nlo2_debug_isg(tb_pb(:,:,1), tb_xb(:,1), 37)
       call nlo2_debug_isg(tb_pb(:,:,1), tb_xb(:,1), 38)
       call nlo2_debug_isg(tb_pb(:,:,1), tb_xb(:,1), 41)
       call nlo2_debug_isg(tb_pb(:,:,1), tb_xb(:,1), 42)
       call nlo2_debug_isg(tb_pb(:,:,1), tb_xb(:,1), 45)
       call nlo2_debug_isg(tb_pb(:,:,1), tb_xb(:,1), 57)
       return
    endif
    if (mode == 2) then
       ! debug: initial-state collinear limits of groups 1 (S2) and 3 (S1)
       call nlo2_debug_is(tb_pb(:,:,1), tb_xb(:,1), 1, 2)
       call nlo2_debug_is(tb_pb(:,:,1), tb_xb(:,1), 3, 3)
       call nlo2_debug_is(tb_pb(:,:,1), tb_xb(:,1), 3, 1)
       return
    endif
    call nlo2_test_limits(tb_pb(:,:,1), tb_xb(:,1))
  end subroutine cs_excl_testlimits

  ! debug: initial-state collinear limit of the H+3j matrix element
  ! (u d -> d u g H, W fusion, gluon from the incoming u) against the VBF
  ! Born for W fusion, |M|^2 ~ (pa.pb)(p1.p2)/((q1^2-MW^2)^2 (q2^2-MW^2)^2):
  ! R3 (2 pa.pg x)/(P_qq(x) B2) must not depend on x
  subroutine debug_is3()
    use cs_dipoles, only: split_fi
    real(dp) :: pb(0:3,5), p6(0:3,6), pa(0:3), pg(0:3), pj(0:3), x, m1, m2, b2, q1s, q2s, r(5)
    real(dp) :: xs(5)
    integer :: ix, iphi
    xs = [0.2_dp, 0.4_dp, 0.6_dp, 0.8_dp, 0.95_dp]
    pb = tb_pb(:,:,1)
    q1s = mdot(pb(:,4) - pb(:,1), pb(:,4) - pb(:,1))
    q2s = mdot(pb(:,5) - pb(:,2), pb(:,5) - pb(:,2))
    b2 = mdot(pb(:,1), pb(:,2))*mdot(pb(:,4), pb(:,5))/((q1s - 80.398_dp**2)**2*(q2s - 80.398_dp**2)**2)
    do iphi = 1, 2
       do ix = 1, 5
          x = xs(ix)
          call split_fi(pb(:,4), pb(:,1), x, 1d-7, 1.0_dp*iphi, pg, pj, pa)
          p6(:,1:5) = pb
          p6(:,1) = pa; p6(:,4) = pj; p6(:,6) = pg
          call cs_hjjj_lines(p6, [2, 1, 25, 1, 2, 0], m1, m2)
          r(ix) = m1*2*mdot(pa, pg)*x/((1 + x**2)/(1 - x)*b2)
       enddo
       write(6,'(a,i2,a,5es13.5)') ' debug_is3 phi', iphi, '  R3 (2 pa.pg x)/(P(x) B2) at x = .2,.4,.6,.8,.95:', r
    enddo
  end subroutine debug_is3

  ! renormalisation-scale dependence of V + I at three-parton points:
  ! d(V + I)/d ln mu^2 = (11/6 CA - 2/3 TR nf) B alpha_s/(2 pi) for a Born
  ! with one power of alpha_s
  subroutine cs_excl_testvirt()
    real(dp) :: pb(0:3,5), pin(0:3), a(0:3), b(0:3), xp, z, wrad, p6(0:3,6), r(3)
    real(dp) :: mu(3), vi(3), born, bmn(0:3,0:3), virt, lo, Qo, pred, worst, iop
    real(dp), parameter :: CF = 4.0_dp/3.0_dp, twopi = 6.283185307179586476925286766559005768394_dp
    real(dp), parameter :: b0 = 11.0_dp/6.0_dp*3 - 2.0_dp/3.0_dp*0.5_dp*5
    integer :: line, it, i1, i2, im, bflav(6), btype, o
    logical :: ok
    worst = 0
    do it = 1, 4
       pb = tb_pb(:,:,1 + mod(it,2))
       do line = 1, 2
          o = 3 - line
          r = [0.2_dp + 0.15_dp*it, 0.3_dp + 0.1_dp*it, 0.7_dp]
          call line_radiation(pb(:,line), pb(:,3+line), tb_xb(line, 1 + mod(it,2)), r, 1, 1d-9, &
               & pin, a, b, xp, z, wrad, ok)
          p6(:,1:5) = pb
          p6(:,line) = pin; p6(:,3+line) = a; p6(:,6) = b
          Qo = sqrt(2*mdot(pb(:,o), pb(:,3+o)))
          do i1 = 1, ncls + ngcls
             i2 = 1 + mod(it, 4)
             if (i1 <= ncls) then
                if (.not. compatible(cls(i1)%w, cls(i2)%w)) cycle
                if (line == 1) bflav = [cls(i1)%a, cls(i2)%a, 25, cls(i1)%b, cls(i2)%b, 0]
                if (line == 2) bflav = [cls(i2)%a, cls(i1)%a, 25, cls(i2)%b, cls(i1)%b, 0]
                btype = 1
             else
                if (.not. compatible(gcls(i1-ncls)%w, cls(i2)%w)) cycle
                if (line == 1) bflav = [0, cls(i2)%a, 25, gcls(i1-ncls)%q, cls(i2)%b, gcls(i1-ncls)%qb]
                if (line == 2) bflav = [cls(i2)%a, 0, 25, cls(i2)%b, gcls(i1-ncls)%q, gcls(i1-ncls)%qb]
                btype = 2
             endif
             mu = [50.0_dp, 200.0_dp, 800.0_dp]
             if (it == 1 .and. line == 1 .and. (i1 == 1 .or. i1 == ncls + 1)) then
                do im = 1, 7
                   call cs_hjjj_born_line(p6, bflav, line, born, bmn)
                   call cs_hjjj_virt_line(p6, bflav, line, (25.0_dp*2**(im-1))**2, virt)
                   lo = log((25.0_dp*2**(im-1))**2/Qo**2)
                   iop = nlo2_ifin(btype, p6(:,line), p6(:,3+line), p6(:,6), (25.0_dp*2**(im-1))**2)
                   write(6,'(a,i2,a,f8.1,a,4es14.6)') ' pieces btype', btype, ' mu', 25.0_dp*2**(im-1), &
                        & '  virt/B, c11, I/(2pi), sum:', virt/born, CF*(-8 - lo**2 - 3*lo)/twopi, iop/twopi, &
                        & (virt - born*CF*(-8 - lo**2 - 3*lo)/twopi + born*iop/twopi)/born
                enddo
             endif
             do im = 1, 3
                call cs_hjjj_born_line(p6, bflav, line, born, bmn)
                call cs_hjjj_virt_line(p6, bflav, line, mu(im)**2, virt)
                lo = log(mu(im)**2/Qo**2)
                iop = nlo2_ifin(btype, p6(:,line), p6(:,3+line), p6(:,6), mu(im)**2)
                vi(im) = (virt - born*CF*(-8 - lo**2 - 3*lo)/twopi + born*iop/twopi)/born
             enddo
             pred = b0/twopi*log(mu(3)**2/mu(1)**2)
             worst = max(worst, abs((vi(3) - vi(1))/pred - 1))
             write(6,'(a,i2,a,i2,a,6i4,a,3es13.5,a,f10.6)') ' testvirt line', line, ' btype', btype, &
                  & ' flav', bflav, '  (V+I)/B at mu = 50, 200, 800:', vi, &
                  & '  [VI(800)-VI(50)]/(b0 ln/(2pi)):', (vi(3) - vi(1))/pred
          enddo
       enddo
    enddo
    write(6,'(a,es10.2)') ' testvirt: max deviation of the mu dependence from b0: ', worst
  end subroutine cs_excl_testvirt

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
