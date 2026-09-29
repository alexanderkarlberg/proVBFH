!----------------------------------------------------------------------
! Catani-Seymour dipole maps for one VBF quark line (DIS-like: one
! incoming parton pa, final-state partons; the line's momentum transfer
! q = sum(final) - pa is kept fixed by all maps, and so are the Higgs and
! the other line). Massless, d = 4.
!
! Maps from n+1 to n final partons of the line (i emitted, merged with j
! into ij, spectator k or the incoming parton):
!   FF (ij,k):  y = pi.pj/(pi.pj + pi.pk + pj.pk),  z = pi.pk/(pi.pk + pj.pk)
!               ptk = pk/(1-y),  ptij = pi + pj - y/(1-y) pk
!   FI (ij;a):  x = 1 - pi.pj/((pi+pj).pa),  z = pi.pa/((pi+pj).pa)
!               ptij = pi + pj - (1-x) pa,  pta = x pa
!   IF (ai;k):  x = 1 - pi.pk/((pi+pk).pa),  u = pi.pa/((pi+pk).pa)
!               ptk = pk + pi - (1-x) pa,  pta = x pa  (a -> ai emits i)
! and their inverses (the splittings used to generate n+1 from n
! partons). Phase-space factorisation (CS, d = 4), with the line measure
! (dxi/xi) dPhi_n unchanged under pa -> x pa:
!   FF: dPhi_{n+1} = dPhi_n (2 ptij.ptk)/(16 pi^2) (1-y) dy dz dphi/(2 pi)
!   FI: (dxi/xi) dPhi_{n+1} = (dxi~/xi~) dPhi_n (2 ptij.pa)/(16 pi^2) dx dz dphi/(2 pi)
!   IF: (dxi/xi) dPhi_{n+1} = (dxi~/xi~) dPhi_n (2 ptk.pa)/(16 pi^2) dx du dphi/(2 pi)
! (pa the incoming momentum of the n+1 configuration).
!
! gen_four generates the line's four-parton configuration (incoming pa,
! three final partons k) from its Born line, DISENT-style: a
! three-parton configuration with line_radiation (the FI splitting of the
! Born), then one more FF or FI splitting of one of its two final
! partons, the outputs in a random order. four_weight is the resulting
! multichannel weight 1/density, from the 24 ways (splitting pair, FF or
! FI, which three-parton final parton split, order) of reaching the
! point, so that
!   (dxi/xi) dPhi_{H+4} = (dx1/x1) dPhi_{H+2} * w * d^7 r
! as line_radiation's wrad for three partons. All samplings are
! logarithmic down to the cutoff (y, 1-x, min(z,1-z) >= cutoff), and the
! seven random numbers must be uniform (VEGAS dimensions that are not
! adapted).
!
! The kernels: dipole = C * 8 pi alpha_s/(2 pi.pj) [/x for FI] * H, or
! C * 8 pi alpha_s/(2 pa.pi x) * H for IF, with C = -T_k.T_ij (numbers
! for a line: its partons form a colour singlet with at most one gluon
! at Born level) and H the CS splitting function V (d = 4) divided by
! 8 pi alpha_s T_ij^2 and contracted with the Born (B, and B^{mu nu} for
! a gluon emitter, -g_{mu nu} B^{mu nu} = B).
! Momenta (E, px, py, pz), index 0:3.
!----------------------------------------------------------------------
module cs_dipoles
  use cs_kinematics, only: mdot
  implicit none
  private
  public :: map_ff, map_fi, map_if, split_ff, split_fi, split_if, dip_basis
  public :: gen_four, four_weight, contract
  public :: h_ff_qg, h_fi_qg, h_ff_gg, h_fi_gg, h_qqb, h_if_qg, h_if_gq, h_if_qq, h_if_gg

  integer, parameter :: dp = kind(1d0)
  real(dp), parameter :: pi = 3.141592653589793238462643383279502884197_dp
  real(dp), parameter :: CF = 4.0_dp/3.0_dp, CA = 3.0_dp, TR = 0.5_dp

  ! diagnostic output of four_weight (per path)
  logical, public, save :: fw_verbose = .false.
  ! diagnostic: azimuthally averaged spin correlations (not for physics runs)
  logical, public, save :: spin_avg = .false.

contains

  !-------------------------------------------------------------------
  ! forward maps: (pi, pj, pk / pa) -> tilde momenta and variables
  subroutine map_ff(pi_, pj, pk, ptij, ptk, y, z)
    real(dp), intent(in)  :: pi_(0:3), pj(0:3), pk(0:3)
    real(dp), intent(out) :: ptij(0:3), ptk(0:3), y, z
    real(dp) :: dij, dik, djk
    dij = mdot(pi_, pj); dik = mdot(pi_, pk); djk = mdot(pj, pk)
    y = dij/(dij + dik + djk)
    z = dik/(dik + djk)
    ptk = pk/(1 - y)
    ptij = pi_ + pj - y/(1 - y)*pk
  end subroutine map_ff

  subroutine map_fi(pi_, pj, pa, ptij, pta, x, z)
    real(dp), intent(in)  :: pi_(0:3), pj(0:3), pa(0:3)
    real(dp), intent(out) :: ptij(0:3), pta(0:3), x, z
    real(dp) :: dij, dia, dja
    dij = mdot(pi_, pj); dia = mdot(pi_, pa); dja = mdot(pj, pa)
    x = 1 - dij/(dia + dja)
    z = dia/(dia + dja)
    ptij = pi_ + pj - (1 - x)*pa
    pta = x*pa
  end subroutine map_fi

  subroutine map_if(pi_, pk, pa, ptk, pta, x, u)
    real(dp), intent(in)  :: pi_(0:3), pk(0:3), pa(0:3)
    real(dp), intent(out) :: ptk(0:3), pta(0:3), x, u
    real(dp) :: dik, dia, dka
    dik = mdot(pi_, pk); dia = mdot(pi_, pa); dka = mdot(pk, pa)
    x = 1 - dik/(dia + dka)
    u = dia/(dia + dka)
    ptk = pk + pi_ - (1 - x)*pa
    pta = x*pa
  end subroutine map_if

  !-------------------------------------------------------------------
  ! inverse maps (splittings), with the azimuth phi of pi around the
  ! (ptij, ptk) or (ptij, pa) or (ptk, pa) plane, relative to the basis
  ! of dip_basis(v1, v2, e1, e2)
  subroutine split_ff(ptij, ptk, y, z, phi, pi_, pj, pk)
    real(dp), intent(in)  :: ptij(0:3), ptk(0:3), y, z, phi
    real(dp), intent(out) :: pi_(0:3), pj(0:3), pk(0:3)
    real(dp) :: e1(0:3), e2(0:3), s, kt
    ! pi = z ptij + y(1-z) ptk + kt, pj = (1-z) ptij + y z ptk - kt,
    ! pk = (1-y) ptk, kt^2 = -2 ptij.ptk y z (1-z)
    call dip_basis(ptij, ptk, e1, e2)
    s = 2*mdot(ptij, ptk)
    kt = sqrt(max(0.0_dp, s*y*z*(1 - z)))
    pi_ = z*ptij + y*(1 - z)*ptk + kt*(cos(phi)*e1 + sin(phi)*e2)
    pj = (1 - z)*ptij + y*z*ptk - kt*(cos(phi)*e1 + sin(phi)*e2)
    pk = (1 - y)*ptk
  end subroutine split_ff

  subroutine split_fi(ptij, pta, x, z, phi, pi_, pj, pa)
    real(dp), intent(in)  :: ptij(0:3), pta(0:3), x, z, phi
    real(dp), intent(out) :: pi_(0:3), pj(0:3), pa(0:3)
    real(dp) :: e1(0:3), e2(0:3), s, kt
    ! pa = pta/x; from pi + pj = ptij + (1-x) pa, pi.pa/(pi+pj).pa = z,
    ! pi^2 = pj^2 = 0: pi = z ptij + (1-z)(1-x) pa + kt,
    ! kt^2 = -2 ptij.pa (1-x) z (1-z)
    call dip_basis(ptij, pta, e1, e2)
    pa = pta/x
    s = 2*mdot(ptij, pa)
    kt = sqrt(max(0.0_dp, s*(1 - x)*z*(1 - z)))
    pi_ = z*ptij + (1 - z)*(1 - x)*pa + kt*(cos(phi)*e1 + sin(phi)*e2)
    pj = (1 - z)*ptij + z*(1 - x)*pa - kt*(cos(phi)*e1 + sin(phi)*e2)
  end subroutine split_fi

  subroutine split_if(ptk, pta, x, u, phi, pi_, pk, pa)
    real(dp), intent(in)  :: ptk(0:3), pta(0:3), x, u, phi
    real(dp), intent(out) :: pi_(0:3), pk(0:3), pa(0:3)
    real(dp) :: e1(0:3), e2(0:3), s, kt
    ! pa = pta/x; pi + pk = ptk + (1-x) pa, pi.pa/(pi+pk).pa = u:
    ! pi = u ptk + (1-u)(1-x) pa + kt (the same map as FI with z -> u)
    call dip_basis(ptk, pta, e1, e2)
    pa = pta/x
    s = 2*mdot(ptk, pa)
    kt = sqrt(max(0.0_dp, s*(1 - x)*u*(1 - u)))
    pi_ = u*ptk + (1 - u)*(1 - x)*pa + kt*(cos(phi)*e1 + sin(phi)*e2)
    pk = u*(1 - x)*pa + (1 - u)*ptk - kt*(cos(phi)*e1 + sin(phi)*e2)
  end subroutine split_if

  !-------------------------------------------------------------------
  ! unit spacelike e1, e2 orthogonal to the massless v1, v2 (and to each
  ! other), from the lab axes
  subroutine dip_basis(v1, v2, e1, e2)
    real(dp), intent(in)  :: v1(0:3), v2(0:3)
    real(dp), intent(out) :: e1(0:3), e2(0:3)
    real(dp) :: ref(0:3), n2, d
    integer :: i, got
    d = mdot(v1, v2)
    got = 0
    do i = 1, 3
       ref = 0; ref(i) = 1
       if (got == 0) then
          e1 = ref - mdot(ref, v2)/d*v1 - mdot(ref, v1)/d*v2
          n2 = -mdot(e1, e1)
          if (n2 > 1d-6) then
             e1 = e1/sqrt(n2); got = 1
          endif
       else
          e2 = ref - mdot(ref, v2)/d*v1 - mdot(ref, v1)/d*v2
          e2 = e2 + mdot(e2, e1)*e1
          n2 = -mdot(e2, e2)
          if (n2 > 1d-6) then
             e2 = e2/sqrt(n2); return
          endif
       endif
    enddo
    stop 'cs_dipoles: no transverse basis'
  end subroutine dip_basis

  !-------------------------------------------------------------------
  ! four-parton line configuration from the Born line (pB incoming, pOB
  ! outgoing, momentum fraction xB); see the header
  subroutine gen_four(pB, pOB, xB, r, cutoff, pa, k, w, ok)
    use cs_kinematics, only: line_radiation
    real(dp), intent(in)  :: pB(0:3), pOB(0:3), xB, r(7), cutoff
    real(dp), intent(out) :: pa(0:3), k(0:3,3), w
    logical,  intent(out) :: ok
    real(dp) :: pin(0:3), a(0:3), b(0:3), xp, z3, wrad, e(0:3), o(0:3)
    real(dp) :: y, x, z, phi, xi3, q(0:3,3)
    integer :: ipath, iperm, islot, ityp
    integer, parameter :: perm(3,6) = reshape([1,2,3, 1,3,2, 2,1,3, 2,3,1, 3,1,2, 3,2,1], [3,6])
    ok = .false.; pa = 0; k = 0; w = 0
    call line_radiation(pB, pOB, xB, r(1:3), 0, cutoff, pin, a, b, xp, z3, wrad, ok)
    if (.not. ok) return
    ok = .false.
    ipath = min(int(24*r(4)), 23)
    iperm = mod(ipath, 6) + 1
    islot = mod(ipath/6, 2)
    ityp = ipath/12
    if (islot == 0) then
       e = a; o = b
    else
       e = b; o = a
    endif
    z = symlog(r(6), cutoff)
    phi = 2*pi*r(7)
    if (ityp == 0) then
       y = cutoff**(1 - r(5))
       call split_ff(e, o, y, z, phi, q(:,1), q(:,2), q(:,3))
       pa = pin
    else
       xi3 = xB/xp
       if (1 - xi3 <= cutoff) return
       x = 1 - cutoff*exp(log((1 - xi3)/cutoff)*r(5))
       call split_fi(e, pin, x, z, phi, q(:,1), q(:,2), pa)
       q(:,3) = o
    endif
    k(:,perm(1,iperm)) = q(:,1)
    k(:,perm(2,iperm)) = q(:,2)
    k(:,perm(3,iperm)) = q(:,3)
    w = four_weight(pB, pOB, xB, pa, k, cutoff)
    ok = w > 0
  contains
    real(dp) function symlog(rr, c)
      real(dp), intent(in) :: rr, c
      if (rr < 0.5_dp) then
         symlog = c*exp(log(0.5_dp/c)*2*rr)
      else
         symlog = 1 - c*exp(log(0.5_dp/c)*2*(1 - rr))
      endif
    end function symlog
  end subroutine gen_four

  ! multichannel weight 1/density of gen_four at the line configuration
  ! (pa, k(:,1:3)) with Born line (pB, pOB, xB); 0 if no path reaches it
  real(dp) function four_weight(pB, pOB, xB, pa, k, cutoff) result(w)
    real(dp), intent(in) :: pB(0:3), pOB(0:3), xB, pa(0:3), k(0:3,3), cutoff
    real(dp) :: ptij(0:3), ptk(0:3), pta(0:3), y, z, x, xi3, xp3, z3, wr, ws, dens, Q2
    real(dp) :: lc, lh
    integer :: ip, i, j, l, ityp
    integer, parameter :: pairs(3,3) = reshape([1,2,3, 1,3,2, 2,3,1], [3,3])
    Q2 = 2*mdot(pB, pOB)
    lc = log(1/cutoff); lh = log(0.5_dp/cutoff)
    dens = 0
    do ip = 1, 3
       i = pairs(1,ip); j = pairs(2,ip); l = pairs(3,ip)
       do ityp = 0, 1
          if (ityp == 0) then
             call map_ff(k(:,i), k(:,j), k(:,l), ptij, ptk, y, z)
             if (y < cutoff .or. min(z, 1 - z) < cutoff) cycle
             ws = 2*mdot(ptij, ptk)/(16*pi**2)*(1 - y)*(y*lc)*(2*min(z, 1 - z)*lh)
             pta = pa
          else
             call map_fi(k(:,i), k(:,j), pa, ptij, pta, x, z)
             ptk = k(:,l)
             xi3 = xB*pta(0)/pB(0)
             if (1 - x < cutoff .or. min(z, 1 - z) < cutoff .or. 1 - xi3 <= cutoff) cycle
             ws = 2*mdot(ptij, pa)/(16*pi**2)*((1 - x)*log((1 - xi3)/cutoff))*(2*min(z, 1 - z)*lh)
          endif
          ! three-parton configuration (pta; ptij, ptk) from the Born by
          ! line_radiation (logarithmic sampling)
          xp3 = 1 - mdot(ptij, ptk)/mdot(pta, ptij + ptk)
          z3 = mdot(pta, ptij)/mdot(pta, ptij + ptk)
          if (1 - xp3 < cutoff .or. min(z3, 1 - z3) < cutoff .or. 1 - xB <= cutoff) cycle
          wr = Q2/(16*pi**2)/xp3*((1 - xp3)*log((1 - xB)/cutoff))*(2*min(z3, 1 - z3)*lh)
          ! 2 orders of (i,j) x 2 slots of the split parton, each 1/24
          dens = dens + 4/(24*wr*ws)
          if (fw_verbose) write(6,'(a,i2,a,i2,a,es10.3,a,es10.3,a,es10.3,a,2es10.2)') '     path pair', ip, ' typ', ityp, &
               & '  1/density =', 24*wr*ws/4, '  wr =', wr, '  ws =', ws, '  seed 1-xp3, z3 =', 1 - xp3, z3
       enddo
    enddo
    w = 0
    if (dens > 0) w = 1/dens
  end function four_weight

  !-------------------------------------------------------------------
  ! v_mu v_nu B^{mu nu}
  real(dp) function contract(v, bmn)
    real(dp), intent(in) :: v(0:3), bmn(0:3,0:3)
    real(dp) :: vl(0:3)
    integer :: m, n
    vl = [v(0), -v(1), -v(2), -v(3)]
    contract = 0
    do m = 0, 3
       do n = 0, 3
          contract = contract + vl(m)*vl(n)*bmn(m,n)
       enddo
    enddo
  end function contract

  ! kernels H (see the header). zq: momentum fraction of the quark
  real(dp) function h_ff_qg(zq, y, b) result(h)
    real(dp), intent(in) :: zq, y, b
    h = (2/(1 - zq*(1 - y)) - (1 + zq))*b
  end function h_ff_qg

  real(dp) function h_fi_qg(zq, x, b) result(h)
    real(dp), intent(in) :: zq, x, b
    h = (2/(2 - zq - x) - (1 + zq))*b
  end function h_fi_qg

  ! g -> g_i g_j, FF: v = zi pi - zj pj
  real(dp) function h_ff_gg(pi_, pj, zi, zj, y, b, bmn) result(h)
    real(dp), intent(in) :: pi_(0:3), pj(0:3), zi, zj, y, b, bmn(0:3,0:3)
    h = 2*((1/(1 - zi*(1 - y)) + 1/(1 - zj*(1 - y)) - 2)*b &
         & + sc(zi, zj, pi_, pj, b, bmn))
  end function h_ff_gg

  ! contract(zi pi - zj pj, B)/(pi.pj), or its azimuthal average zi zj B
  real(dp) function sc(zi, zj, pi_, pj, b, bmn)
    real(dp), intent(in) :: zi, zj, pi_(0:3), pj(0:3), b, bmn(0:3,0:3)
    if (spin_avg) then
       sc = zi*zj*b
    else
       sc = contract(zi*pi_ - zj*pj, bmn)/mdot(pi_, pj)
    endif
  end function sc

  real(dp) function h_fi_gg(pi_, pj, zi, zj, x, b, bmn) result(h)
    real(dp), intent(in) :: pi_(0:3), pj(0:3), zi, zj, x, b, bmn(0:3,0:3)
    h = 2*((1/(2 - zi - x) + 1/(2 - zj - x) - 2)*b &
         & + sc(zi, zj, pi_, pj, b, bmn))
  end function h_fi_gg

  ! g -> q_i qbar_j (FF and FI)
  real(dp) function h_qqb(pi_, pj, zi, zj, b, bmn) result(h)
    real(dp), intent(in) :: pi_(0:3), pj(0:3), zi, zj, b, bmn(0:3,0:3)
    h = TR/CA*(b - 2*sc(zi, zj, pi_, pj, b, bmn))
  end function h_qqb

  ! IF: incoming quark emits a gluon (the Born has the quark)
  real(dp) function h_if_qg(x, u, b) result(h)
    real(dp), intent(in) :: x, u, b
    h = (2/(1 - x + u) - (1 + x))*b
  end function h_if_qg

  ! IF: incoming gluon emits a quark or antiquark (the Born has the
  ! antiquark or quark)
  real(dp) function h_if_gq(x, b) result(h)
    real(dp), intent(in) :: x, b
    h = TR/CF*(1 - 2*x*(1 - x))*b
  end function h_if_gq

  ! IF: incoming quark emits the same quark (the Born has a gluon);
  ! w = pi/u - pk/(1-u)
  real(dp) function h_if_qq(pi_, pk, x, u, b, bmn) result(h)
    real(dp), intent(in) :: pi_(0:3), pk(0:3), x, u, b, bmn(0:3,0:3)
    h = CF/CA*(x*b + (1 - x)/x*2*u*(1 - u)/mdot(pi_, pk)*contract(pi_/u - pk/(1 - u), bmn))
  end function h_if_qq

  ! IF: incoming gluon emits a gluon
  real(dp) function h_if_gg(pi_, pk, x, u, b, bmn) result(h)
    real(dp), intent(in) :: pi_(0:3), pk(0:3), x, u, b, bmn(0:3,0:3)
    h = 2*((1/(1 - x + u) - 1 + x*(1 - x))*b &
         & + (1 - x)/x*u*(1 - u)/mdot(pi_, pk)*contract(pi_/u - pk/(1 - u), bmn))
  end function h_if_gg

end module cs_dipoles
