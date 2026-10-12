!----------------------------------------------------------------------
! Stage 2 of proVBFH-cs (docs/DESIGN.md): the (2,0) and (0,2)
! contributions, the NLO correction to VBF H + 3 partons from the line
! that radiated, with Catani-Seymour subtraction for that line (a DIS
! process: one incoming parton, its partons a colour singlet with the
! other line, the Higgs and the other line fixed by all maps).
!
! Real: VBF H + 4 partons with both extra partons on the line, proVBFH's
! real flavour list with its line tags (init_processes; classes 11 and
! 22), evaluated with its VBFNLO matrix elements (setreal). The line's
! part of an entry has one of four structures (roles 1, 2, 3 of its final
! partons; T = line tag, P = the pair tag 10 + T):
!   S1  g(T)  -> Q(T) Qbar(T) g(T)       incoming gluon
!   S2  q(T)  -> Q(T) g(T) g(T)
!   S3  q(T)  -> Q(T) q'(P) qbar'(P)     pair from a gluon
!   S4  q'(P) -> Q(T) Qbar(T) q'(P)      incoming q' emits the gluon that
!                                        fuses with the boson into Q Qbar
! (Q, Qbar: the quarks the boson couples to). Its dipoles are those of
! POWHEG's region finder for the same tags (mergetags): pairs that merge
! into a valid Born with the extra parton on the line. Entries are
! grouped by the value of their matrix element and dipoles, checked at
! initialisation; the matrix elements are evaluated once per group and
! multiplied by the sum of the members' PDFs.
!
! The dipoles need the H+3j Born with the gluon on the line and its spin
! correlations (cs_hjjj_born_line). Colour: -T_k.T_ij is a number for a
! line (Born partons: two quarks, one gluon).
!
! Integrated dipoles: I (nlo2_ifin) and K + P (nlo2_kp, a convolution
! done by Gauss quadrature, with the functions of DISENT's KPFUNS in the
! MSbar scheme), at the three-parton events of stage 1.
!
! All matrix elements at alpha_s = 1; weights here are sums over the
! flavours of |M|^2 x PDFs (number densities) times the phase-space
! weight; the caller multiplies by the Born flux and Jacobian and alpha_s.
! Momenta (E, px, py, pz), index 0:3.
!----------------------------------------------------------------------
module cs_nlo2
  use types, only: dp
  use cs_kinematics, only: mdot, line_radiation
  use cs_dipoles
  use cs_memo, only: memo_on, memo_get, memo_put, memo_nkey
  implicit none
  private
  public :: nlo2_init, nlo2_real_kin, nlo2_real_me, nlo2_ifin, nlo2_kp
  public :: nlo2_kp_dis, nlo2_ifin_dis
  public :: nlo2_ngroups, nlo2_ncount
  ! four-parton points dropped by the technical cut
  integer, public, save :: nlo2_ncut = 0

  real(dp), parameter :: pi = 3.141592653589793238462643383279502884197_dp
  real(dp), parameter :: CF = 4.0_dp/3.0_dp, CA = 3.0_dp, TR = 0.5_dp
  integer, parameter :: nf = 5
  real(dp), parameter :: gam_q = 1.5_dp*CF, gam_g = 11.0_dp/6.0_dp*CA - 2.0_dp/3.0_dp*TR*nf
  real(dp), parameter :: K_q = (3.5_dp - pi**2/6)*CF
  real(dp), parameter :: K_g = (67.0_dp/18.0_dp - pi**2/6)*CA - 10.0_dp/9.0_dp*TR*nf

  ! real groups
  integer, parameter :: maxgrp = 400, maxmem = 400
  type real_group
     integer :: line = 0, s = 0, rep = 0, n = 0
     ! legs of the representative: 1 incoming of the line, 2 incoming of
     ! the other line, 3 Higgs, 4:6 final partons of the line in role
     ! order, 7 outgoing of the other line
     integer :: leg(7) = 0
     integer :: fl(7) = 0          ! flavours of the representative, same order
     integer :: fa(maxmem) = 0, fo(maxmem) = 0   ! members: incoming flavours
     real(dp) :: fp(8) = 0         ! fingerprint (matrix elements at test points)
     logical :: nc = .true.        ! NC: the Born-level line keeps its flavour
  end type real_group
  type(real_group), allocatable, save :: grp(:)
  integer, save :: ngrp = 0
  ! counters: real ME calls, Born ME calls (for timing studies)
  integer(8), save :: nlo2_ncount(2) = 0
  external :: hoppetEval

  ! dipoles of each structure (see the table in nlo2_init)
  integer, parameter :: maxdip = 10
  type dipole
     integer :: typ = 0      ! 1 FF, 2 FI, 3 IF
     integer :: i = 0, j = 0 ! FF, FI: the pair (i the quark of qg kernels); IF: i emitted, j spectator
     integer :: kern = 0     ! 1 qg, 2 gg, 3 q qbar, 4 IF q->g, 5 IF g->q, 6 IF q->q, 7 IF g->g
     real(dp) :: col = 0     ! -T_k.T_ij
     integer :: bin = 0      ! Born incoming: 1 the real's, 0 gluon, -r minus the flavour of role r
     integer :: bout = 0     ! role whose flavour the Born's outgoing quark has
     integer :: bext = 0     ! Born extra parton: 0 gluon, r > 0 the flavour of role r
     integer :: mout = 0, mext = 0 ! momenta: -1 merged pair (for IF: the spectator), -2 FF spectator, r role r
  end type dipole
  ! structure 5: S3 with the line quark of the incoming flavour (NC), whose
  ! VBFNLO matrix element also has the graphs of S4 (the incoming quark
  ! emits a t-channel gluon that fuses with the boson into the pair):
  ! the S3 dipoles plus the initial-state q -> q ones of those graphs
  type(dipole), save :: dip(maxdip,5)
  integer, save :: ndip(5) = 0
  real(dp), parameter :: symf(5) = [1.0_dp, 0.5_dp, 1.0_dp, 1.0_dp, 1.0_dp]
  integer, parameter :: pairs(3,3) = reshape([1,2,3, 1,3,2, 2,3,1], [3,3])

  ! kinematics of the current point (nlo2_real_kin -> nlo2_real_me)
  integer, save :: cur_line = 0
  real(dp), save :: cur_pb(0:3,5), cur_xb(2), cur_pa(0:3), cur_k(0:3,3), cur_w4, cur_cutoff = 0
  real(dp), save :: mt(0:3,6), mk(0:3,6), ma(0:3,6), mvar(6), mz(6)
  logical, save :: mok(6)

  ! Born cache of the current point
  integer, parameter :: maxcache = 200
  integer, save :: ncache = 0
  integer, save :: ckey(10,maxcache)
  real(dp), save :: cb(maxcache), cbmn(0:3,0:3,maxcache)

  ! Gauss-Legendre nodes on [0,1] for the K + P convolutions
  integer, parameter :: nq = 32
  real(dp), save :: gx(nq), gw(nq)

contains

  integer function nlo2_ngroups(line)
    integer, intent(in) :: line
    nlo2_ngroups = count(grp(1:ngrp)%line == line)
  end function nlo2_ngroups

  !-------------------------------------------------------------------
  ! set up: the dipole tables, the real groups from proVBFH's flavour list
  ! (with the check that all members of a group give the same matrix
  ! elements and dipoles at test points), the quadrature
  subroutine nlo2_init(pbtest, xbtest, cutoff, verbose)
    real(dp), intent(in) :: pbtest(0:3,5,2), xbtest(2,2), cutoff
    logical, intent(in) :: verbose
    integer :: nreal, j, line, s, leg(7), fl(7), flav(7), tags(7), ig, ng0, nmem, ntot(2)
    real(dp) :: fp(8), fprep(8)
    logical :: found
    call dipole_tables()
    call gauleg(nq, gx, gw)
    call cs_init_reals(nreal)
    allocate(grp(maxgrp))
    ngrp = 0; ntot = 0
    do j = 1, nreal
       call cs_get_real(j, flav, tags)
       do line = 1, 2
          if (.not. line_entry(flav, tags, line, s, leg, fl)) cycle
          ntot(line) = ntot(line) + 1
          call fingerprint(j, line, s, leg, fl, pbtest, xbtest, cutoff, fp)
          found = .false.
          do ig = 1, ngrp
             if (grp(ig)%line /= line .or. grp(ig)%s /= s) cycle
             if (maxval(abs(fp - grp(ig)%fp)/max(abs(grp(ig)%fp), 1d-300)) < 1d-9) then
                nmem = grp(ig)%n + 1
                if (nmem > maxmem) stop 'cs_nlo2: too many members'
                grp(ig)%n = nmem
                grp(ig)%fa(nmem) = fl(1)
                grp(ig)%fo(nmem) = fl(2)
                if (grp(ig)%nc .neqv. (fl(2) == fl(7))) stop 'cs_nlo2: NC and CC entries in one group'
                found = .true.
                exit
             endif
          enddo
          if (.not. found) then
             ngrp = ngrp + 1
             if (ngrp > maxgrp) stop 'cs_nlo2: too many groups'
             grp(ngrp)%line = line; grp(ngrp)%s = s; grp(ngrp)%rep = j
             grp(ngrp)%leg = leg; grp(ngrp)%fl = fl; grp(ngrp)%fp = fp
             grp(ngrp)%n = 1; grp(ngrp)%fa(1) = fl(1); grp(ngrp)%fo(1) = fl(2)
             grp(ngrp)%nc = fl(2) == fl(7)
          endif
       enddo
    enddo
    if (verbose) then
       do line = 1, 2
          write(6,'(a,i2,a,i6,a,i5,a,4i5)') ' cs_nlo2: line', line, ': real entries', ntot(line), &
               & ', groups', count(grp(1:ngrp)%line == line), ', per structure', &
               & (count(grp(1:ngrp)%line == line .and. grp(1:ngrp)%s == s), s = 1, 4)
       enddo
       ng0 = 0
       do ig = 1, ngrp
          if (ng0 < 60) write(6,'(a,i4,a,i2,a,i2,a,7i4,a,i4)') '   group', ig, ' line', grp(ig)%line, &
               & ' S', grp(ig)%s, ' flav (in, other in, H, roles, other out)', grp(ig)%fl, &
               & '  members', grp(ig)%n
          ng0 = ng0 + 1
       enddo
    endif
  end subroutine nlo2_init

  ! is real entry (flav, tags) one with both extra partons on line? Its
  ! structure s, the legs and flavours in the group order
  logical function line_entry(flav, tags, line, s, leg, fl)
    integer, intent(in) :: flav(7), tags(7), line
    integer, intent(out) :: s, leg(7), fl(7)
    integer :: o, i, nfin, fin(4), ng, np, nq_, lpair
    line_entry = .false.
    s = 0; leg = 0; fl = 0
    o = 3 - line
    lpair = 10 + line
    ! both extra partons on the line: tags sum 8 or 28 (line 1), 10 or 30 (line 2)
    if (line == 1 .and. sum(tags) /= 8 .and. sum(tags) /= 28) return
    if (line == 2 .and. sum(tags) /= 10 .and. sum(tags) /= 30) return
    if (tags(o) /= o) return
    nfin = 0
    do i = 4, 7
       if (mod(tags(i), 10) == line) then
          nfin = nfin + 1
          if (nfin > 3) return
          fin(nfin) = i
       elseif (tags(i) == o) then
          leg(7) = i
       endif
    enddo
    if (nfin /= 3 .or. leg(7) == 0) return
    leg(1) = line; leg(2) = o; leg(3) = 3
    ng = count(flav(fin(1:3)) == 0)
    np = count(tags(fin(1:3)) == lpair)
    if (flav(line) == 0) then
       s = 1
       if (ng /= 1 .or. np /= 0) stop 'cs_nlo2: unexpected S1 entry'
       do i = 1, 3
          if (flav(fin(i)) > 0) leg(4) = fin(i)
          if (flav(fin(i)) < 0) leg(5) = fin(i)
          if (flav(fin(i)) == 0) leg(6) = fin(i)
       enddo
    elseif (tags(line) == line .and. ng == 2) then
       s = 2
       nq_ = 0
       do i = 1, 3
          if (flav(fin(i)) /= 0) then
             leg(4) = fin(i)
          else
             nq_ = nq_ + 1
             leg(4+nq_) = fin(i)
          endif
       enddo
    elseif (tags(line) == line .and. np == 2) then
       s = 3
       do i = 1, 3
          if (tags(fin(i)) == line) leg(4) = fin(i)
          if (tags(fin(i)) == lpair .and. flav(fin(i)) > 0) leg(5) = fin(i)
          if (tags(fin(i)) == lpair .and. flav(fin(i)) < 0) leg(6) = fin(i)
       enddo
    elseif (tags(line) == lpair .and. np == 1) then
       s = 4
       do i = 1, 3
          if (tags(fin(i)) == line .and. flav(fin(i)) > 0) leg(4) = fin(i)
          if (tags(fin(i)) == line .and. flav(fin(i)) < 0) leg(5) = fin(i)
          if (tags(fin(i)) == lpair) leg(6) = fin(i)
       enddo
       if (flav(leg(6)) /= flav(line)) stop 'cs_nlo2: unexpected S4 entry'
    else
       write(6,*) 'flav', flav, ' tags', tags
       stop 'cs_nlo2: unknown line structure'
    endif
    if (any(leg == 0)) stop 'cs_nlo2: roles not found'
    fl = flav(leg)
    line_entry = .true.
  end function line_entry

  ! the table of dipoles per structure (roles 1, 2, 3; see the header)
  subroutine dipole_tables()
    real(dp) :: c1, c2
    c1 = CF - CA/2      ! -T.T of the two quarks of a Born line
    c2 = CA/2           ! -T.T of a quark and the gluon
    ! S1: g -> Q(1) Qbar(2) g(3); Born g -> Q Qbar, or qbar -> Qbar g (IF of Q),
    ! q -> Q g (IF of Qbar)
    ndip(1) = 10
    dip(1,1) = dipole(1, 1,3, 1, c1,  0, 1, 2, -1, -2)
    dip(2,1) = dipole(2, 1,3, 1, c2,  0, 1, 2, -1,  2)
    dip(3,1) = dipole(1, 2,3, 1, c1,  0, 1, 2, -2, -1)
    dip(4,1) = dipole(2, 2,3, 1, c2,  0, 1, 2,  1, -1)
    dip(5,1) = dipole(3, 3,1, 7, c2,  0, 1, 2, -1,  2)
    dip(6,1) = dipole(3, 3,2, 7, c2,  0, 1, 2,  1, -1)
    dip(7,1) = dipole(3, 1,2, 5, c1, -1, 2, 0, -1,  3)
    dip(8,1) = dipole(3, 1,3, 5, c2, -1, 2, 0,  2, -1)
    dip(9,1) = dipole(3, 2,1, 5, c1, -2, 1, 0, -1,  3)
    dip(10,1)= dipole(3, 2,3, 5, c2, -2, 1, 0,  1, -1)
    ! S2: q -> Q(1) g(2) g(3); Born q -> Q g
    ndip(2) = 10
    dip(1,2) = dipole(1, 1,2, 1, c2,  1, 1, 0, -1, -2)
    dip(2,2) = dipole(2, 1,2, 1, c1,  1, 1, 0, -1,  3)
    dip(3,2) = dipole(1, 1,3, 1, c2,  1, 1, 0, -1, -2)
    dip(4,2) = dipole(2, 1,3, 1, c1,  1, 1, 0, -1,  2)
    dip(5,2) = dipole(1, 2,3, 2, c2,  1, 1, 0, -2, -1)
    dip(6,2) = dipole(2, 2,3, 2, c2,  1, 1, 0,  1, -1)
    dip(7,2) = dipole(3, 2,1, 4, c1,  1, 1, 0, -1,  3)
    dip(8,2) = dipole(3, 2,3, 4, c2,  1, 1, 0,  1, -1)
    dip(9,2) = dipole(3, 3,1, 4, c1,  1, 1, 0, -1,  2)
    dip(10,2)= dipole(3, 3,2, 4, c2,  1, 1, 0,  1, -1)
    ! S3: q -> Q(1) q'(2) qbar'(3); Born q -> Q g
    ndip(3) = 2
    dip(1,3) = dipole(1, 2,3, 3, c2,  1, 1, 0, -2, -1)
    dip(2,3) = dipole(2, 2,3, 3, c2,  1, 1, 0,  1, -1)
    ! S3 (NC, Q of the incoming flavour): also q -> Q(1) || q, Born g -> q'(2) qbar'(3)
    ndip(5) = 4
    dip(1:2,5) = dip(1:2,3)
    dip(3,5) = dipole(3, 1,2, 6, c2,  0, 2, 3, -1,  3)
    dip(4,5) = dipole(3, 1,3, 6, c2,  0, 2, 3,  2, -1)
    ! S4: q' -> Q(1) Qbar(2) q'(3); Born g -> Q Qbar
    ndip(4) = 2
    dip(1,4) = dipole(3, 3,1, 6, c2,  0, 1, 2, -1,  2)
    dip(2,4) = dipole(3, 3,2, 6, c2,  0, 1, 2,  1, -1)
  end subroutine dipole_tables

  ! index of the map of the pair (i, j) of roles: FF 1..3, FI-type 4..6
  integer function map_index(i, j, typ)
    integer, intent(in) :: i, j, typ
    integer :: lo, hi, ip
    lo = min(i, j); hi = max(i, j)
    ip = 3
    if (lo == 1 .and. hi == 2) ip = 1
    if (lo == 1 .and. hi == 3) ip = 2
    map_index = ip + merge(0, 3, typ == 1)
  end function map_index

  !-------------------------------------------------------------------
  ! kinematics of the (2,0) (line 1) or (0,2) (line 2) real contribution
  ! at a point: the four-parton line from gen_four (seven uniform random
  ! numbers r) and the six mapped three-parton configurations (FF of the
  ! pairs (1,2), (1,3), (2,3) with the third as spectator; FI-type of the
  ! same pairs, which is also the IF map of either parton with the other as
  ! spectator). Events for the analysis (1, 2 incoming, 3 Higgs, partons):
  ! pr (7 momenta) and pd(:,:,m) (6 momenta); dok(m): the mapped
  ! configuration m is inside the three-parton cutoff (dipole used).
  subroutine nlo2_real_kin(line, pb, xb, r, cutoff, pr, pd, dok, ok, nocut)
    integer, intent(in) :: line
    logical, intent(in), optional :: nocut
    real(dp), intent(in) :: pb(0:3,5), xb(2), r(7), cutoff
    real(dp), intent(out) :: pr(0:3,7), pd(0:3,6,6)
    logical, intent(out) :: dok(6), ok
    integer :: m, i, j, l, o
    real(dp) :: xp3, z3, y, z, x
    cur_line = line; cur_pb = pb; cur_xb = xb; cur_cutoff = cutoff
    o = 3 - line
    pr = 0; pd = 0; dok = .false.
    call gen_four(pb(:,line), pb(:,3+line), xb(line), r, cutoff, cur_pa, cur_k, cur_w4, ok)
    if (.not. ok) return
    pr(:,line) = cur_pa; pr(:,o) = pb(:,o); pr(:,3) = pb(:,3)
    pr(:,4:6) = cur_k; pr(:,7) = pb(:,3+o)
    do m = 1, 6
       i = pairs(1, mod(m-1,3)+1); j = pairs(2, mod(m-1,3)+1); l = pairs(3, mod(m-1,3)+1)
       if (m <= 3) then
          call map_ff(cur_k(:,i), cur_k(:,j), cur_k(:,l), mt(:,m), mk(:,m), y, z)
          ma(:,m) = cur_pa
          mvar(m) = y; mz(m) = z
       else
          call map_fi(cur_k(:,i), cur_k(:,j), cur_pa, mt(:,m), ma(:,m), x, z)
          mk(:,m) = cur_k(:,l)
          mvar(m) = x; mz(m) = z
       endif
       ! the three-parton cutoff of line_radiation
       xp3 = 1 - mdot(mt(:,m), mk(:,m))/mdot(ma(:,m), mt(:,m) + mk(:,m))
       z3 = mdot(ma(:,m), mt(:,m))/mdot(ma(:,m), mt(:,m) + mk(:,m))
       mok(m) = 1 - xp3 >= cutoff .and. min(z3, 1 - z3) >= cutoff
       dok(m) = mok(m)
       pd(:,line,m) = ma(:,m); pd(:,o,m) = pb(:,o); pd(:,3,m) = pb(:,3)
       pd(:,4,m) = mt(:,m); pd(:,5,m) = mk(:,m); pd(:,6,m) = pb(:,3+o)
    enddo
    ! technical cut, consistent with the support of gen_four: the whole
    ! four-parton point (real and counterevents) is dropped if any map is
    ! closer to its singular limit than the cutoff (FF: y, z, 1-z; FI: 1-x,
    ! z, 1-z). Without it the real's singular limits beyond the cutoff are
    ! reached by other paths with a tiny density (weights up to 1e13).
    if (present(nocut)) then
       if (nocut) return
    endif
    do m = 1, 6
       if (m <= 3) then
          if (mvar(m) < cutoff) ok = .false.
       else
          if (1 - mvar(m) < cutoff) ok = .false.
       endif
       if (min(mz(m), 1 - mz(m)) < cutoff) ok = .false.
    enddo
    if (.not. ok) then
       dok = .false.
       nlo2_ncut = nlo2_ncut + 1
    endif
  end subroutine nlo2_real_kin

  ! weights of the real event (wr) and the six counterevents (wd, with the
  ! minus sign of the subtraction) of the point set up by nlo2_real_kin,
  ! at the factorisation scale muf of the line; fo: number densities of
  ! the other line at its Born momentum fraction. Includes the phase-space
  ! weight of gen_four.
  subroutine nlo2_real_me(muf, fo, wr, wd)
    real(dp), intent(in) :: muf, fo(-6:6)
    real(dp), intent(out) :: wr, wd(6)
    real(dp) :: xi4, f4(-6:6), pdf, me, d(6)
    integer :: ig, i
    wr = 0; wd = 0
    ncache = 0
    xi4 = cur_xb(cur_line)*cur_pa(0)/cur_pb(0,cur_line)
    call hoppetEval(xi4, muf, f4)
    f4 = f4/xi4
    do ig = 1, ngrp
       if (grp(ig)%line /= cur_line) cycle
       pdf = 0
       do i = 1, grp(ig)%n
          pdf = pdf + f4(grp(ig)%fa(i))*fo(grp(ig)%fo(i))
       enddo
       if (pdf == 0) cycle
       if (memo_on) then
          call group_me_memo(ig, me, d)
       else
          call group_me(ig, cur_pa, cur_k, me, d, .true.)
       endif
       wr = wr + me*pdf
       wd = wd - d*pdf
    enddo
    wr = wr*cur_w4
    wd = wd*cur_w4
  end subroutine nlo2_real_me

  ! group_me(ig, cur_pa, cur_k, me, d, .true.) through the per-point cache
  ! (scale variations): the real and its dipoles do not depend on the scales
  subroutine group_me_memo(ig, me, d)
    integer, intent(in) :: ig
    real(dp), intent(out) :: me, d(6)
    integer :: key(memo_nkey), slot
    real(dp) :: mom(16), val(7)
    key = 0; key(1) = 5; key(2) = ig; key(3) = cur_line
    mom(1:4) = cur_pa; mom(5:16) = reshape(cur_k, [12])
    if (memo_get(key, mom, 16, val, 7, slot)) then
       me = val(1); d = val(2:7)
    else
       call group_me(ig, cur_pa, cur_k, me, d, .true.)
       call memo_put(slot, key, mom, 16, [me, d], 7)
    endif
  end subroutine group_me_memo

  ! real matrix element of group ig at the line configuration (pa, k) of
  ! the current Born point, and its dipoles summed per mapped
  ! configuration (only those with mok, if usemok)
  subroutine group_me(ig, pa, k, me, d, usemok)
    integer, intent(in) :: ig
    real(dp), intent(in) :: pa(0:3), k(0:3,3)
    real(dp), intent(out) :: me, d(6)
    logical, intent(in) :: usemok
    real(dp) :: p7(0:3,7), p6(0:3,6), b, bmn(0:3,0:3), h, pre, zq, zi, u, v(0:3)
    real(dp) :: pi_(0:3), pj(0:3)
    integer :: bfl(6), line, o, s, id, m, key(10), i
    type(dipole) :: dp_
    type(real_group) :: g
    g = grp(ig)
    line = g%line; o = 3 - line; s = g%s
    if (s == 3 .and. g%fl(4) == g%fl(1)) s = 5
    p7(:,g%leg(1)) = pa
    p7(:,g%leg(2)) = cur_pb(:,o)
    p7(:,g%leg(3)) = cur_pb(:,3)
    p7(:,g%leg(4)) = k(:,1); p7(:,g%leg(5)) = k(:,2); p7(:,g%leg(6)) = k(:,3)
    p7(:,g%leg(7)) = cur_pb(:,3+o)
    d = 0
    call cs_real_me(g%rep, p7, me)
    nlo2_ncount(1) = nlo2_ncount(1) + 1
    do id = 1, ndip(s)
       dp_ = dip(id, s)
       m = map_index(dp_%i, dp_%j, dp_%typ)
       if (usemok .and. .not. mok(m)) cycle
       ! Born flavours and momenta (H+3j order, extra parton on the line)
       bfl(o) = g%fl(2); bfl(3+o) = g%fl(7); bfl(3) = g%fl(3)
       if (dp_%bin == 1) then
          bfl(line) = g%fl(1)
       elseif (dp_%bin == 0) then
          bfl(line) = 0
       else
          bfl(line) = -g%fl(3 - dp_%bin)
       endif
       bfl(3+line) = g%fl(3 + dp_%bout)
       bfl(6) = 0
       if (dp_%bext > 0) bfl(6) = g%fl(3 + dp_%bext)
       p6(:,o) = cur_pb(:,o); p6(:,3+o) = cur_pb(:,3+o); p6(:,3) = cur_pb(:,3)
       p6(:,line) = ma(:,m)
       p6(:,3+line) = pick(dp_%mout)
       p6(:,6) = pick(dp_%mext)
       key = [m, bfl, dp_%mout, dp_%mext, line]
       call born_cached(key, p6, bfl, line, b, bmn)
       ! kernel
       pi_ = k(:,dp_%i); pj = k(:,dp_%j)
       select case (dp_%typ)
       case (1)
          pre = 8*pi/(2*mdot(pi_, pj))
       case (2)
          pre = 8*pi/(2*mdot(pi_, pj)*mvar(m))
       case default
          pre = 8*pi/(2*mdot(pa, pi_)*mvar(m))
       end select
       ! z of role i in the map of the pair (FF: relative to the spectator;
       ! FI: relative to pa); mz(m) is that of the lower role of the pair
       zi = mz(m)
       if (dp_%i > dp_%j) zi = 1 - mz(m)
       select case (dp_%kern)
       case (1)
          zq = zi
          if (dp_%typ == 1) then
             h = h_ff_qg(zq, mvar(m), b)
          else
             h = h_fi_qg(zq, mvar(m), b)
          endif
       case (2)
          if (dp_%typ == 1) then
             h = h_ff_gg(pi_, pj, zi, 1 - zi, mvar(m), b, bmn)
          else
             h = h_fi_gg(pi_, pj, zi, 1 - zi, mvar(m), b, bmn)
          endif
       case (3)
          h = h_qqb(pi_, pj, zi, 1 - zi, b, bmn)
       case (4)
          h = h_if_qg(mvar(m), zi, b)
       case (5)
          h = h_if_gq(mvar(m), b)
       case (6)
          h = h_if_qq(pi_, pj, mvar(m), zi, b, bmn)
       case (7)
          h = h_if_gg(pi_, pj, mvar(m), zi, b, bmn)
       end select
       d(m) = d(m) + symf(s)*dp_%col*pre*h
    enddo
  contains
    function pick(code) result(p)
      integer, intent(in) :: code
      real(dp) :: p(0:3)
      if (code == -1) then
         p = mt(:,m)
      elseif (code == -2) then
         p = mk(:,m)
      else
         p = k(:,code)
      endif
    end function pick
  end subroutine group_me

  subroutine born_cached(key, p6, bfl, line, b, bmn)
    integer, intent(in) :: key(10), bfl(6), line
    real(dp), intent(in) :: p6(0:3,6)
    real(dp), intent(out) :: b, bmn(0:3,0:3)
    integer :: ic
    do ic = 1, ncache
       if (all(ckey(:,ic) == key)) then
          b = cb(ic); bmn = cbmn(:,:,ic)
          return
       endif
    enddo
    call cs_hjjj_born_line(p6, bfl, line, b, bmn)
    nlo2_ncount(2) = nlo2_ncount(2) + 1
    if (ncache < maxcache) then
       ncache = ncache + 1
       ckey(:,ncache) = key; cb(ncache) = b; cbmn(:,:,ncache) = bmn
    endif
  end subroutine born_cached

  ! fingerprint of a real entry: its matrix element and dipoles at two
  ! test configurations
  subroutine fingerprint(j, line, s, leg, fl, pbtest, xbtest, cutoff, fp)
    integer, intent(in) :: j, line, s, leg(7), fl(7)
    real(dp), intent(in) :: pbtest(0:3,5,2), xbtest(2,2), cutoff
    real(dp), intent(out) :: fp(8)
    real(dp) :: r(7), pr(0:3,7), pd(0:3,6,6), me, d(6)
    logical :: dok(6), ok
    integer :: it, save_n
    save_n = ngrp
    ! a temporary group in the last slot
    grp(maxgrp)%line = line; grp(maxgrp)%s = s; grp(maxgrp)%rep = j
    grp(maxgrp)%leg = leg; grp(maxgrp)%fl = fl
    do it = 1, 2
       r = [0.31_dp, 0.62_dp, 0.17_dp, (it - 0.5_dp)/2.0_dp*0.93_dp + 0.02_dp, 0.41_dp, 0.73_dp, 0.29_dp]
       ! fixed test points: matrix elements only, no technical cut
       call nlo2_real_kin(line, pbtest(:,:,it), xbtest(:,it), r, cutoff, pr, pd, dok, ok, nocut=.true.)
       if (.not. ok) stop 'cs_nlo2: test point not generated'
       ncache = 0
       call group_me(maxgrp, cur_pa, cur_k, me, d, .false.)
       fp(4*it-3) = me
       fp(4*it-2) = sum(d(1:3))
       fp(4*it-1) = sum(d(4:6))
       fp(4*it) = d(1) + 2*d(2) + 3*d(3) + 4*d(4) + 5*d(5) + 6*d(6)
    enddo
    grp(maxgrp) = real_group()
    if (ngrp /= save_n) stop 'cs_nlo2: fingerprint'
  end subroutine fingerprint

  !-------------------------------------------------------------------
  ! finite part of the I operator for a Born line with incoming parton
  ! pa and final partons p1, p2 (btype 1: q -> Q(p1) g(p2); btype 2:
  ! g -> Q(p1) Qbar(p2)), with the virtual's normalisation (CDR,
  ! (4 pi)^eps/Gamma(1-eps) (mu^2)^eps removed), at mu^2 = mur2:
  !   -sum_I sum_{J/=I} T_I.T_J/T_I^2 [T_I^2 L^2/2 + gamma_I L + gamma_I
  !                                    + K_I - T_I^2 pi^2/3],
  ! L = ln(mur2/(2 pI.pJ)). The line's contribution is alpha_s/(2 pi) I B.
  real(dp) function nlo2_ifin(btype, pa, p1, p2, mur2) result(res)
    integer, intent(in) :: btype
    real(dp), intent(in) :: pa(0:3), p1(0:3), p2(0:3), mur2
    real(dp) :: p(0:3,3), t2(3), gm(3), kk(3), tt(3,3), l
    integer :: i, j
    p(:,1) = pa; p(:,2) = p1; p(:,3) = p2
    tt = 0
    if (btype == 1) then
       t2 = [CF, CF, CA]; gm = [gam_q, gam_q, gam_g]; kk = [K_q, K_q, K_g]
       tt(1,2) = CA/2 - CF; tt(1,3) = -CA/2; tt(2,3) = -CA/2
    else
       t2 = [CA, CF, CF]; gm = [gam_g, gam_q, gam_q]; kk = [K_g, K_q, K_q]
       tt(2,3) = CA/2 - CF; tt(1,2) = -CA/2; tt(1,3) = -CA/2
    endif
    tt = tt + transpose(tt)
    res = 0
    do i = 1, 3
       do j = 1, 3
          if (i == j) cycle
          l = log(mur2/(2*mdot(p(:,i), p(:,j))))
          res = res - tt(i,j)/t2(i)*(t2(i)*l**2/2 + gm(i)*l + gm(i) + kk(i) - t2(i)*pi**2/3)
       enddo
    enddo
  end function nlo2_ifin

  !-------------------------------------------------------------------
  ! K + P for a three-parton line event with incoming momentum fraction
  ! xi, incoming pa and final partons p1, p2 (the same momenta for both
  ! Born types: q -> Q(p1) g(p2) and g -> Q(p1) Qbar(p2)), at the
  ! factorisation scale muf: effective number densities
  !   fq(a') = sum_a [(K + P)^{a a'} (x) f_a](xi)  for a quark Born with
  !            incoming a' (a' /= 0),
  !   fg     = sum_a [(K + P)^{a g} (x) f_a](xi)   for a gluon Born,
  ! with (g (x) f)(xi) = int_xi^1 dz/z g(z) f(xi/z), the plus distributions
  ! on [0,1] (DISENT's KPFUNS, MSbar, the xmin terms giving the integral
  ! from 0 to xi). The line's contribution is alpha_s/(2 pi) B_a' fq(a')
  ! (or B_g fg) in place of B_a' f_a'(xi).
  subroutine nlo2_kp(xi, pa, p1, p2, muf, fq, fg, fgq)
    real(dp), intent(in) :: xi, pa(0:3), p1(0:3), p2(0:3), muf
    real(dp), intent(out) :: fq(-6:6), fg
    ! fgq (optional): the quark part of fg, sum_q [(K + P)^{q g} (x) f_q](xi)
    real(dp), intent(out), optional :: fgq
    real(dp) :: kqf, kgf, lsc, lsg, s1, s2
    ! colour factors (as DISENT's COLFOR and VIRTHR)
    kqf = (1.5_dp*(CF - CA/2) + 0.5_dp*gam_g)/CF
    kgf = 1.5_dp
    s1 = log(muf**2/(2*mdot(pa, p1)))
    s2 = log(muf**2/(2*mdot(pa, p2)))
    ! log(scale) + PQF = -sum_I T_I.T_a/T_a^2 ln(muf^2/(2 pa.pI))
    lsc = -((CA/2 - CF)/CF*s1 - CA/(2*CF)*s2)
    lsg = 0.5_dp*(s1 + s2)
    call kp_core(xi, muf, kqf, kgf, lsc, lsg, fq, fg, fgq)
  end subroutine nlo2_kp

  ! K + P for a DIS line (stage 3): Born with incoming quark pa and
  ! outgoing quark pb (colour: T_a.T_b = -CF), so kqf = -sum_I
  ! T_I.T_a gamma_I/(T_I^2 T_a^2) = 3/2 and lsc = ln(muf^2/(2 pa.pb)).
  ! fq(a') as in nlo2_kp (the line's contribution is alpha_s/(2 pi) B
  ! fq(a') in place of B f_a'(xi)).
  subroutine nlo2_kp_dis(xi, pa, pb, muf, fq)
    real(dp), intent(in) :: xi, pa(0:3), pb(0:3), muf
    real(dp), intent(out) :: fq(-6:6)
    real(dp) :: lsc, fg
    lsc = log(muf**2/(2*mdot(pa, pb)))
    call kp_core(xi, muf, 1.5_dp, 1.5_dp, lsc, lsc, fq, fg)
  end subroutine nlo2_kp_dis

  ! finite part of the I operator of a DIS line (incoming quark pa,
  ! outgoing quark pb), normalised as nlo2_ifin:
  !   2 [CF L^2/2 + gamma_q L + gamma_q + K_q - CF pi^2/3], L = ln(mur2/(2 pa.pb))
  real(dp) function nlo2_ifin_dis(pa, pb, mur2) result(res)
    real(dp), intent(in) :: pa(0:3), pb(0:3), mur2
    real(dp) :: l
    l = log(mur2/(2*mdot(pa, pb)))
    res = 2*(CF*l**2/2 + gam_q*l + gam_q + K_q - CF*pi**2/3)
  end function nlo2_ifin_dis

  ! the convolutions of nlo2_kp for given colour-dependent factors
  subroutine kp_core(xi, muf, kqf, kgf, lsc, lsg, fq, fg, fgq)
    real(dp), intent(in) :: xi, muf, kqf, kgf, lsc, lsg
    real(dp), intent(out) :: fq(-6:6), fg
    real(dp), intent(out), optional :: fgq
    real(dp) :: lm, dl, zm, z, wz, l, f(-6:6), f1(-6:6), fsum
    real(dp) :: qqp, qqr, gqr, qgr, ggp, ggr, qqd, ggd
    integer :: ip, iq, a
    ! the delta terms (with the integral of the plus distributions from 0 to xi)
    lm = log(1 - xi)
    dl = li2(1 - xi)
    call hoppetEval(xi, muf, f1)
    f1 = f1/xi
    qqd = -CF*(5 - pi**2 + kqf + pi**2/3 - lm**2 - 2*dl + kqf*lm + (2*lm + xi + xi**2/2)*lsc)
    ggd = -CA*(50.0_dp/9 - pi**2 + kgf + pi**2/3 - lm**2 - 2*dl + kgf*lm + 2*lm*lsg) &
         & + TR*nf*16.0_dp/9 - gam_g*lsg
    fq = qqd*f1
    fg = ggd*f1(0)
    if (present(fgq)) fgq = 0
    ! the z integral: [xi, zm] logarithmic, [zm, 1] with 1-z = (1-zm) v^2
    zm = (1 + xi)/2
    do ip = 1, 2
       do iq = 1, nq
          if (ip == 1) then
             z = xi*(zm/xi)**gx(iq)
             wz = gw(iq)*z*log(zm/xi)
          else
             z = 1 - (1 - zm)*gx(iq)**2
             wz = gw(iq)*2*(1 - zm)*gx(iq)
          endif
          call hoppetEval(xi/z, muf, f)
          f = f/xi               ! f_a(xi/z)/z
          l = log((1 - z)/z)
          qqp = CF*2/(1 - z)*(l - kqf/2) - CF*(1 + z**2)/(1 - z)*lsc
          ggp = CA*2/(1 - z)*(l - kgf/2) - CA*2/(1 - z)*lsg
          qqr = CF*(-(1 + z)*l + (1 - z))
          gqr = TR*((z**2 + (1 - z)**2)*(l - lsc) + 2*z*(1 - z))
          qgr = CF*((1 + (1 - z)**2)/z*(l - lsg) + z)
          ggr = CA*((1 - z)/z - 1 + z*(1 - z))*2*(l - lsg)
          fsum = 0
          do a = -nf, nf
             if (a /= 0) fsum = fsum + f(a)
          enddo
          do a = -6, 6
             if (a == 0) cycle
             fq(a) = fq(a) + wz*(qqp*(f(a) - f1(a)) + qqr*f(a) + gqr*f(0))
          enddo
          fg = fg + wz*(ggp*(f(0) - f1(0)) + ggr*f(0) + qgr*fsum)
          if (present(fgq)) fgq = fgq + wz*qgr*fsum
       enddo
    enddo
    fq(0) = 0
  end subroutine kp_core

  ! dilogarithm Li2(x), 0 <= x <= 1
  real(dp) function li2(x)
    real(dp), intent(in) :: x
    real(dp) :: y, s, t
    integer :: k
    if (x > 0.5_dp) then
       ! Li2(x) = pi^2/6 - ln(x) ln(1-x) - Li2(1-x)
       y = 1 - x
       s = 0; t = 1
       do k = 1, 200
          t = t*y
          s = s + t/k**2
          if (t < 1d-18) exit
       enddo
       li2 = pi**2/6 - log(x)*log(max(y, 1d-300)) - s
    else
       s = 0; t = 1
       do k = 1, 200
          t = t*x
          s = s + t/k**2
          if (t < 1d-18) exit
       enddo
       li2 = s
    endif
  end function li2

  ! Gauss-Legendre nodes and weights on [0,1]
  subroutine gauleg(n, x, w)
    integer, intent(in) :: n
    real(dp), intent(out) :: x(n), w(n)
    integer :: i, j, it
    real(dp) :: z, z1, p1, p2, p3, pp
    do i = 1, (n + 1)/2
       z = cos(pi*(i - 0.25_dp)/(n + 0.5_dp))
       do it = 1, 100
          p1 = 1; p2 = 0
          do j = 1, n
             p3 = p2; p2 = p1
             p1 = ((2*j - 1)*z*p2 - (j - 1)*p3)/j
          enddo
          pp = n*(z*p1 - p2)/(z*z - 1)
          z1 = z
          z = z1 - p1/pp
          if (abs(z - z1) < 1d-15) exit
       enddo
       x(i) = (1 - z)/2; x(n+1-i) = (1 + z)/2
       w(i) = 1/((1 - z*z)*pp*pp); w(n+1-i) = w(i)
    enddo
  end subroutine gauleg

end module cs_nlo2
