!----------------------------------------------------------------------
! Diagnostic channel split (notes/2026-10-cs-production, 8 Oct 2026).
!
! Initial state: flags per beam on the parton type of every PDF
! evaluation, f(0) (gluon) times chg(beam), f(+-1..6) (quarks and
! antiquarks) times chq(beam). Applied after every hoppetEval of the
! exclusive part (Born-level, three-parton and four-parton momentum
! fractions, and the x/z of the K + P convolutions), so that a real and
! its dipoles, and the K + P terms, are classified by the parton taken
! from the proton. qq + qg + gq + gg = all, exactly.
!
! Boson: NC (both lines neutral, Z exchange) or CC (W exchange), i.e. the
! class pairs selected by compatible() in cs_exclusive and the real
! groups of cs_nlo2 (NC: the Born-level line keeps its flavour).
!
! Environment:
!   CHAN_BOSON = NC | CC | ALL (default ALL)
!   CHAN_INIT  = ALL (default) or two letters, beam 1 then beam 2, each
!                q (quark or antiquark), g (gluon) or a (any)
!   CHAN_MULTI = 1: nine weights per event in one run, W1 = all,
!                W2-W5 = NC qq, qg, gq, gg, W6-W9 = CC qq, qg, gq, gg
!                (CHAN_BOSON and CHAN_INIT must then be unset or ALL)
!----------------------------------------------------------------------
module cs_chan
  use types, only: dp
  implicit none
  private
  public :: chan_init, chan_set, chan_pdf, chan_w_ok, chan_multi, chan_trivial
  public :: kp_beam
  real(dp), save :: chq(2) = 1, chg(2) = 1
  logical, save :: ch_nc = .true., ch_cc = .true.
  ! the fixed selection from the environment
  real(dp), save :: env_q(2) = 1, env_g(2) = 1
  logical, save :: env_nc = .true., env_cc = .true.
  logical, save :: chan_multi = .false.
  ! the beam of the K + P convolution that follows (set by the caller of
  ! nlo2_kp / nlo2_kp_dis; 0: not set)
  integer, save :: kp_beam = 0

contains

  subroutine chan_init()
    character(len=64) :: s
    integer :: l, b
    call get_environment_variable('CHAN_BOSON', s, l)
    s = upcase(s)
    select case (trim(s))
    case ('', 'ALL')
    case ('NC'); env_cc = .false.
    case ('CC'); env_nc = .false.
    case default; stop 'cs_chan: CHAN_BOSON must be NC, CC or ALL'
    end select
    call get_environment_variable('CHAN_INIT', s, l)
    s = upcase(s)
    if (trim(s) /= '' .and. trim(s) /= 'ALL') then
       if (len_trim(s) /= 2) stop 'cs_chan: CHAN_INIT must be ALL or two of q, g, a'
       do b = 1, 2
          select case (s(b:b))
          case ('Q'); env_g(b) = 0
          case ('G'); env_q(b) = 0
          case ('A')
          case default; stop 'cs_chan: CHAN_INIT must be ALL or two of q, g, a'
          end select
       enddo
    endif
    call get_environment_variable('CHAN_MULTI', s, l)
    chan_multi = trim(s) == '1'
    if (chan_multi .and. .not. (env_nc .and. env_cc .and. all(env_q == 1) .and. all(env_g == 1))) &
         & stop 'cs_chan: CHAN_MULTI with CHAN_BOSON or CHAN_INIT'
    call chan_set(1)
    if (chan_multi) then
       write(6,'(a)') ' cs_chan: CHAN_MULTI, weights W1 all, W2-W5 NC qq qg gq gg, W6-W9 CC qq qg gq gg'
    else
       write(6,'(a,2l2,a,2f4.1,a,2f4.1)') ' cs_chan: NC, CC', ch_nc, ch_cc, '  quark flags (beam 1, 2)', chq, &
            & '  gluon flags', chg
    endif
  end subroutine chan_init

  ! the selection of weight k (CHAN_MULTI), or the environment's (k = 1)
  subroutine chan_set(k)
    integer, intent(in) :: k
    integer :: j
    if (.not. chan_multi .or. k == 1) then
       chq = env_q; chg = env_g; ch_nc = env_nc; ch_cc = env_cc
       return
    endif
    if (k < 2 .or. k > 9) stop 'cs_chan: weight out of range'
    ch_nc = k <= 5; ch_cc = k >= 6
    j = mod(k - 2, 4)          ! 0 qq, 1 qg, 2 gq, 3 gg
    chq = 0; chg = 0
    if (j == 0 .or. j == 1) then
       chq(1) = 1
    else
       chg(1) = 1
    endif
    if (j == 0 .or. j == 2) then
       chq(2) = 1
    else
       chg(2) = 1
    endif
  end subroutine chan_set

  logical function chan_trivial()
    chan_trivial = ch_nc .and. ch_cc .and. all(chq == 1) .and. all(chg == 1)
  end function chan_trivial

  ! number densities f(-6:6) of beam b: apply the parton-type flags
  subroutine chan_pdf(b, f)
    integer, intent(in) :: b
    real(dp), intent(inout) :: f(-6:6)
    if (b /= 1 .and. b /= 2) stop 'cs_chan: chan_pdf without a beam'
    if (chq(b) == 1 .and. chg(b) == 1) return
    f(0) = f(0)*chg(b)
    f(-6:-1) = f(-6:-1)*chq(b)
    f(1:6) = f(1:6)*chq(b)
  end subroutine chan_pdf

  ! a class pair of W charge w (of either line; compatible pairs only)
  logical function chan_w_ok(w)
    integer, intent(in) :: w
    chan_w_ok = (w == 0 .and. ch_nc) .or. (w /= 0 .and. ch_cc)
  end function chan_w_ok

  function upcase(s) result(u)
    character(len=*), intent(in) :: s
    character(len=len(s)) :: u
    integer :: i
    u = s
    do i = 1, len(s)
       if (s(i:i) >= 'a' .and. s(i:i) <= 'z') u(i:i) = achar(iachar(s(i:i)) - 32)
    enddo
  end function upcase

end module cs_chan
