! Independent explicit-Dirac-matrix (Weyl representation) reference for the tree
! amplitudes of  (line a)--g--(line b),  one V on a or b,  V--H--V,  spectator line z.
! Written for the w4q-spinor check; uses nothing from VBFNLO / POWHEG except the
! numerical values of couplings, masses and widths read from common blocks.
module ref
  implicit none
  integer, parameter :: dp = kind(1d0)
  complex(dp), save :: gam(4,4,0:3), gam0(4,4), PL(4,4), PR(4,4)
  real(dp), parameter :: gmet(0:3) = [1d0,-1d0,-1d0,-1d0]
  ! V couplings of the quark types (3 = up-type, 4 = down-type), boson 2 = Z, 3 = W,
  ! chirality (-1 = L, +1 = R): set by the caller from the VBFNLO common block
  real(dp), save :: cpl(3:4,2:3,-1:1)
  logical, save :: ref_inited = .false.
  logical, save :: skip_charge = .false.   ! diagnostic: drop the W+W- charge check
  real(dp), parameter :: colsq = 6d0, colint = -2d0
contains
  subroutine ref_init()
    complex(dp) :: s(2,2,0:3)
    integer :: mu, i, j
    s = 0
    s(1,1,0) = 1; s(2,2,0) = 1
    s(1,2,1) = 1; s(2,1,1) = 1
    s(1,2,2) = (0d0,-1d0); s(2,1,2) = (0d0,1d0)
    s(1,1,3) = 1; s(2,2,3) = -1
    gam = 0
    do mu = 0, 3
       ! gamma^mu = [[0, sigma^mu],[sigmabar^mu, 0]], sigma = (1,sigma_i), sigmabar = (1,-sigma_i)
       gam(1:2,3:4,mu) = s(:,:,mu)
       if (mu == 0) then
          gam(3:4,1:2,mu) = s(:,:,mu)
       else
          gam(3:4,1:2,mu) = -s(:,:,mu)
       end if
    end do
    gam0 = gam(:,:,0)
    PL = 0; PR = 0
    do i = 1, 2
       PL(i,i) = 1
       PR(i+2,i+2) = 1
    end do
    ref_inited = .true.
  end subroutine ref_init

  function slash(p) result(m)
    complex(dp), intent(in) :: p(0:3)
    complex(dp) :: m(4,4)
    integer :: mu
    m = 0
    do mu = 0, 3
       m = m + gmet(mu)*p(mu)*gam(:,:,mu)
    end do
  end function slash

  function minkdot(a,b) result(r)
    real(dp), intent(in) :: a(0:3), b(0:3)
    real(dp) :: r
    r = a(0)*b(0)-a(1)*b(1)-a(2)*b(2)-a(3)*b(3)
  end function minkdot

  ! massless chirality eigenspinor u_c(p), c=-1: left, +1: right; p = physical momentum E>0
  function spinor(p, c) result(u)
    real(dp), intent(in) :: p(0:3)
    integer, intent(in) :: c
    complex(dp) :: u(4)
    real(dp) :: E, ct, ch, sh, phi
    complex(dp) :: chip(2), chim(2)
    E = sqrt(p(1)**2+p(2)**2+p(3)**2)
    ct = p(3)/E
    ch = sqrt(max(0d0,(1d0+ct)/2d0)); sh = sqrt(max(0d0,(1d0-ct)/2d0))
    phi = atan2(p(2),p(1))
    chip = [cmplx(ch,0d0,dp), cmplx(cos(phi),sin(phi),dp)*sh]
    chim = [-cmplx(cos(phi),-sin(phi),dp)*sh, cmplx(ch,0d0,dp)]
    u = 0
    if (c < 0) then
       u(1:2) = sqrt(2d0*E)*chim
    else
       u(3:4) = sqrt(2d0*E)*chip
    end if
  end function spinor

  function bar(u) result(ub)
    complex(dp), intent(in) :: u(4)
    complex(dp) :: ub(4)
    ub = matmul(conjg(u), gam0)
  end function bar

  function cpm(flv, bos) result(m)   ! coupling matrix cL PL + cR PR
    integer, intent(in) :: flv, bos
    complex(dp) :: m(4,4)
    m = cpl(flv,bos,-1)*PL + cpl(flv,bos,1)*PR
  end function cpm

  function curr(ub, u, flv, bos) result(j)   ! J^mu = ubar gamma^mu (cL PL + cR PR) u  (bos=0: pure vector)
    complex(dp), intent(in) :: ub(4), u(4)
    integer, intent(in) :: flv, bos
    complex(dp) :: j(0:3), c(4,4)
    integer :: mu, i
    if (bos == 0) then
       c = 0
       do i = 1, 4
          c(i,i) = 1
       end do
    else
       c = cpm(flv,bos)
    end if
    do mu = 0, 3
       j(mu) = sum(ub*matmul(matmul(gam(:,:,mu),c), u))
    end do
  end function curr

  function chain3(ub, A, B, C, u) result(r)
    complex(dp), intent(in) :: ub(4), A(4,4), B(4,4), C(4,4), u(4)
    complex(dp) :: r
    r = sum(ub*matmul(A, matmul(B, matmul(C, u))))
  end function chain3

  logical function lineok(fs, fe, cc)
    integer, intent(in) :: fs, fe
    logical, intent(in) :: cc
    integer :: part
    lineok = (fs == fe)
    if (cc .and. fs >= 1 .and. fs <= 4) then
       if (mod(fs,2) == 1) then
          part = fs + 1
       else
          part = fs - 1
       end if
       if (fe == part) lineok = .true.
    end if
  end function lineok

  integer function qc(f)   ! charge in units of 1/3, f = 2 up-type, 1 down-type
    integer, intent(in) :: f
    if (f == 2) then
       qc = 2
    else
       qc = -1
    end if
  end function qc

  ! flavour/charge consistency of a pairing: returns okv(1:2) (V on line a / b)
  subroutine pair_flavours(sa, ea, sb, eb, sz, ez, ft, fn, cc, okv)
    integer, intent(in) :: sa, ea, sb, eb, sz, ez, ft(6), fn(6)
    logical, intent(in) :: cc
    logical, intent(out) :: okv(2)
    logical :: chg_a, chg_b, chg_z
    integer :: dq_a, dq_b, dq_z
    okv = .false.
    if (.not. lineok(fn(sa),fn(ea),cc) .or. .not. lineok(fn(sb),fn(eb),cc) .or. .not. lineok(fn(sz),fn(ez),cc)) return
    chg_a = fn(sa) /= fn(ea); chg_b = fn(sb) /= fn(eb); chg_z = fn(sz) /= fn(ez)
    dq_a = qc(ft(ea)) - qc(ft(sa)); dq_b = qc(ft(eb)) - qc(ft(sb)); dq_z = qc(ft(ez)) - qc(ft(sz))
    if (.not. cc) then
       if (.not.(chg_a .or. chg_b .or. chg_z)) okv = .true.
    else
       if (chg_z .and. chg_a .and. .not.chg_b .and. (skip_charge .or. dq_z + dq_a == 0)) okv(1) = .true.
       if (chg_z .and. chg_b .and. .not.chg_a .and. (skip_charge .or. dq_z + dq_b == 0)) okv(2) = .true.
    end if
  end subroutine pair_flavours

  ! Amplitude (without Fermi sign and colour) for the cluster lines a = (sa -> ea), b = (sb -> eb)
  ! connected by a gluon, one V (Z or W, to the spectator line z) on a or b.
  ! kk = arrow momentum of each leg, uu/ubb = (barred) spinors, ft = field flavour type.
  function amp_pair(sa, ea, sb, eb, sz, ez, kk, uu, ubb, ft, fn, cc, valid) result(M)
    integer, intent(in) :: sa, ea, sb, eb, sz, ez, ft(6), fn(6)
    real(dp), intent(in) :: kk(0:3,6)
    complex(dp), intent(in) :: uu(4,6), ubb(4,6)
    logical, intent(in) :: cc
    logical, intent(out) :: valid
    complex(dp) :: M
    integer :: bos, place, fl_a, fl_z, s(2), e(2)
    logical :: okv(2)
    complex(dp) :: Jz(0:3), Jb(0:3), Vm(4,4), Gm(4,4), Sm(4,4)
    real(dp) :: g(0:3), kin(0:3), g2

    bos = 2
    if (cc) bos = 3
    M = 0; valid = .false.
    call pair_flavours(sa, ea, sb, eb, sz, ez, ft, fn, cc, okv)
    if (.not.(okv(1) .or. okv(2))) return
    valid = .true.
    fl_z = 3 + (2 - ft(sz))        ! up-type -> 3, down-type -> 4
    Jz = curr(ubb(:,ez), uu(:,sz), fl_z, bos)
    do place = 1, 2     ! V on line a (1) or b (2)
       if (.not. okv(place)) cycle
       if (place == 1) then
          s = [sa, sb]; e = [ea, eb]
       else
          s = [sb, sa]; e = [eb, ea]
       end if
       ! V-line (s(1)->e(1)); gluon-only line (s(2)->e(2))
       Jb = curr(ubb(:,e(2)), uu(:,s(2)), 3, 0)
       g = kk(:,s(2)) - kk(:,e(2))       ! gluon momentum flowing into the V-line
       g2 = g(0)**2 - g(1)**2 - g(2)**2 - g(3)**2
       fl_a = 3 + (2 - ft(s(1)))
       Vm = matmul(slash(Jz), cpm(fl_a, bos))
       Gm = slash(Jb)
       ! ordering 1: gluon next to the start spinor, then V
       kin = kk(:,s(1)) + g
       Sm = slash(cmplx(kin,0d0,dp))/(kin(0)**2-kin(1)**2-kin(2)**2-kin(3)**2)
       M = M + chain3(ubb(:,e(1)), Vm, Sm, Gm, uu(:,s(1)))/g2
       ! ordering 2: V first, then gluon
       kin = kk(:,e(1)) - g
       Sm = slash(cmplx(kin,0d0,dp))/(kin(0)**2-kin(1)**2-kin(2)**2-kin(3)**2)
       M = M + chain3(ubb(:,e(1)), Gm, Sm, Vm, uu(:,s(1)))/g2
    end do
  end function amp_pair

  ! Entry-level truth.  fl(7), tg(7): flavours and tags of a real-list entry (leg order 1,2 in, 3 H,
  ! 4..7 out), pin(0:3,7) physical momenta.  Returns colour-summed, spin-summed (no coupling/propagator
  ! factors):  Nn = normal pairing squared, Ff = flipped pairing squared, Ii = 2 colint Re(Mn Mf*) (Fermi sign in Mf),
  ! plus the V virtualities pc2 (cluster), ps2 (spectator) and flags.
  subroutine ref_entry(fl, tg, tagfac, pin, cc, Nn, Ff, Ii, pc2, ps2, flipvalid, ident, ok, sigN, sigF, wpair)
    integer, intent(in) :: fl(7), tg(7), tagfac
    real(dp), intent(in) :: pin(0:3,7)
    logical, intent(in) :: cc
    real(dp), intent(out) :: Nn, Ff, Ii, pc2, ps2
    logical, intent(out) :: flipvalid, ident, ok, wpair
    integer, intent(out) :: sigN(3), sigF(3)
    logical :: okn(2), okf(2)
    integer :: dd(6)
    integer, parameter :: idxleg(6) = [1,2,4,5,6,7]
    integer :: i, l, ft(6), fn(6), role(6), tgf(6), lp, c(6), ic, nsp, sz, ez
    integer :: sp, ep, so, eo, np_, nS, nE, flg(6), pflag(6)
    real(dp) :: kk(0:3,6), pph(0:3,6), pc(0:3), ps(0:3), amp2n, amp2f
    complex(dp) :: uu(4,6), ubb(4,6), Mn, Mf
    logical :: isspec(6), v1, v2
    integer :: slist(6), elist(6)

    Nn = 0; Ff = 0; Ii = 0; pc2 = 0; ps2 = 0; flipvalid = .false.; ident = .false.; ok = .false.
    wpair = .false.; sigN = 0; sigF = 0
    do i = 1, 6
       l = idxleg(i)
       flg(i) = fl(l)
       pph(:,i) = pin(:,l)
       ft(i) = 2 - mod(abs(fl(l)),2)
       fn(i) = abs(fl(l))
       if ((l <= 2 .and. fl(l) > 0) .or. (l > 2 .and. fl(l) < 0)) then
          role(i) = +1    ! start of arrow
       else
          role(i) = -1    ! end of arrow
       end if
       if (fl(l) > 0) then
          kk(:,i) = pin(:,l)
       else
          kk(:,i) = -pin(:,l)
       end if
       tgf(i) = mod(tg(l), tagfac)
       pflag(i) = tg(l)/tagfac
    end do
    lp = mod(tg(7), tagfac)
    ! spectator legs: tag_f /= line of the pair
    nsp = 0
    isspec = .false.
    do i = 1, 6
       if (tgf(i) /= lp) then
          isspec(i) = .true.; nsp = nsp + 1
       end if
    end do
    if (nsp /= 2) return
    sz = 0; ez = 0
    do i = 1, 6
       if (isspec(i)) then
          if (role(i) == +1) sz = i
          if (role(i) == -1) ez = i
       end if
    end do
    if (sz == 0 .or. ez == 0) return
    sp = 0; ep = 0; so = 0; eo = 0
    do i = 1, 6
       if (isspec(i)) cycle
       if (pflag(i) >= 1) then
          if (role(i) == +1) sp = i
          if (role(i) == -1) ep = i
       else
          if (role(i) == +1) so = i
          if (role(i) == -1) eo = i
       end if
    end do
    if (min(sp,ep,so,eo) == 0) return
    ok = .true.
    do i = 1, 6
       l = idxleg(i)
       dd(i) = (fl(l)+10)*10
       if (l <= 2) dd(i) = dd(i) + l
    end do
    call pair_flavours(sp,ep,so,eo,sz,ez,ft,fn,cc,okn)
    call pair_flavours(sp,eo,so,ep,sz,ez,ft,fn,cc,okf)
    flipvalid = any(okf)
    call mksig(dd(sp)*1000+dd(ep), dd(so)*1000+dd(eo), dd(sz)*1000+dd(ez), sigN)
    call mksig(dd(sp)*1000+dd(eo), dd(so)*1000+dd(ep), dd(sz)*1000+dd(ez), sigF)
    wpair = (fn(sp) /= fn(eo))      ! the flipped line starting at the pair antiquark carries the W
    ident = (flg(ep) == flg(eo)) .and. ((idxleg(ep) <= 2) .eqv. (idxleg(eo) <= 2))
    ! V virtualities
    pc = 0
    do i = 1, 6
       if (isspec(i)) cycle
       pc = pc - role(i)*kk(:,i)      ! sum_E k - sum_S k
    end do
    pc2 = pc(0)**2-pc(1)**2-pc(2)**2-pc(3)**2
    ps = -role(sz)*kk(:,sz) - role(ez)*kk(:,ez)
    ps2 = ps(0)**2-ps(1)**2-ps(2)**2-ps(3)**2

    call ref_core(sp,ep,so,eo,sz,ez,kk,pph,ft,fn,cc,Nn,Ff,Ii)
  end subroutine ref_entry


  ! colour- and spin-summed |normal|^2, |flipped|^2 and interference for explicit lines
  subroutine ref_core(sp,ep,so,eo,sz,ez,kk,pph,ft,fn,cc,Nn,Ff,Ii)
    integer, intent(in) :: sp,ep,so,eo,sz,ez,ft(6),fn(6)
    real(dp), intent(in) :: kk(0:3,6), pph(0:3,6)
    logical, intent(in) :: cc
    real(dp), intent(out) :: Nn, Ff, Ii
    integer :: ic, i, c(6)
    complex(dp) :: uu(4,6), ubb(4,6), Mn, Mf
    logical :: v1, v2
    Nn = 0; Ff = 0; Ii = 0
    do ic = 0, 63
       do i = 1, 6
          if (btest(ic,i-1)) then
             c(i) = +1
          else
             c(i) = -1
          end if
          uu(:,i) = spinor(pph(:,i), c(i))
          ubb(:,i) = bar(uu(:,i))
       end do
       Mn = amp_pair(sp,ep,so,eo,sz,ez,kk,uu,ubb,ft,fn,cc,v1)
       Mf = -amp_pair(sp,eo,so,ep,sz,ez,kk,uu,ubb,ft,fn,cc,v2)
       Nn = Nn + colsq*abs(Mn)**2
       Ff = Ff + colsq*abs(Mf)**2
       Ii = Ii + 2d0*colint*real(Mn*conjg(Mf))
    end do
  end subroutine ref_core

  ! like ref_core, but the flipped pairing uses its own end-flavour assignment (ftF, fnF); the interference is
  ! formed from the two amplitudes irrespective of whether they describe the same external state
  subroutine ref_core_mixed(sp,ep,so,eo,sz,ez,kk,pph,ft,fn,ftF,fnF,cc,Nn,Ff,Ii)
    integer, intent(in) :: sp,ep,so,eo,sz,ez,ft(6),fn(6),ftF(6),fnF(6)
    real(dp), intent(in) :: kk(0:3,6), pph(0:3,6)
    logical, intent(in) :: cc
    real(dp), intent(out) :: Nn, Ff, Ii
    integer :: ic, i, c(6)
    complex(dp) :: uu(4,6), ubb(4,6), Mn, Mf
    logical :: v1, v2
    Nn = 0; Ff = 0; Ii = 0
    do ic = 0, 63
       do i = 1, 6
          if (btest(ic,i-1)) then
             c(i) = +1
          else
             c(i) = -1
          end if
          uu(:,i) = spinor(pph(:,i), c(i))
          ubb(:,i) = bar(uu(:,i))
       end do
       Mn = amp_pair(sp,ep,so,eo,sz,ez,kk,uu,ubb,ft,fn,cc,v1)
       Mf = -amp_pair(sp,eo,so,ep,sz,ez,kk,uu,ubb,ftF,fnF,cc,v2)
       Nn = Nn + colsq*abs(Mn)**2
       Ff = Ff + colsq*abs(Mf)**2
       Ii = Ii + 2d0*colint*real(Mn*conjg(Mf))
    end do
  end subroutine ref_core_mixed

  ! amplitude for one explicit chirality assignment c(1:6) (external spinor chiralities) of a pairing
  function amp_chir(sa,ea,sb,eb,sz,ez,kk,pph,ft,fn,cc,c) result(M)
    integer, intent(in) :: sa,ea,sb,eb,sz,ez,ft(6),fn(6),c(6)
    real(dp), intent(in) :: kk(0:3,6), pph(0:3,6)
    logical, intent(in) :: cc
    complex(dp) :: M, uu(4,6), ubb(4,6)
    logical :: v
    integer :: i
    do i = 1, 6
       uu(:,i) = spinor(pph(:,i), c(i))
       ubb(:,i) = bar(uu(:,i))
    end do
    M = amp_pair(sa,ea,sb,eb,sz,ez,kk,uu,ubb,ft,fn,cc,v)
  end function amp_chir

  subroutine mksig(c1, c2, cz, sig)
    integer, intent(in) :: c1, c2, cz
    integer, intent(out) :: sig(3)
    sig(1) = min(c1,c2); sig(2) = max(c1,c2); sig(3) = cz
  end subroutine mksig
end module ref
