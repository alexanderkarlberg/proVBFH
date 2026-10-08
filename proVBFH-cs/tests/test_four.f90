!----------------------------------------------------------------------
! Unit test of cs_dipoles: the maps and splittings, and the multichannel
! four-parton generator gen_four with its weight four_weight.
!   1. split_ff/fi/if invert map_ff/fi/if, and all keep q fixed;
!   2. at every generated point the weight from gen_four is the
!      multichannel weight of the point, and the momenta are massless
!      with q fixed;
!   3. integrals of test functions (with cuts on all invariants, inside
!      the generator's support) agree with an independent flat sampling
!      of the line's phase space,
!        (dxi/xi) dPhi_{H+4} = (dx1/x1) dPhi_{H+2} (Q^2/(2 pi)) (du/u)
!                              dPhi_3(W),
!      pa = pB/u, W^2 = Q^2 (1-u)/u, dPhi_3 flat (RAMBO) in the rest
!      frame of pa + q.
! Run with (four_hard, four_hard2, four_hard2mode) = (0, 0, 0), (0.3, 0, 0),
! (0.3, 0.3, 0), (0.3, 0.5, 1): the hard channels of the first and second
! step, and the second step's hard channel only after the first's.
! Exits with status 1 if a check fails.
!----------------------------------------------------------------------
program test_four
  use cs_kinematics, only: mdot
  use cs_dipoles
  implicit none
  integer, parameter :: dp = kind(1d0)
  real(dp), parameter :: pi = 3.141592653589793238462643383279502884197_dp
  real(dp) :: pB(0:3), pOB(0:3), q(0:3), xB, Q2, E, th, ph, r(7), rr(3)
  real(dp) :: pa(0:3), k(0:3,3), w, cutoff, errmax(4), t(0:3,3), a(0:3), b(0:3), c(0:3)
  real(dp) :: y, z, x, u, ptij(0:3), ptk(0:3), pta(0:3)
  real(dp) :: s1(2), s2(2), f(2), m1(2), m2(2), sig1(2), sig2(2), pull
  logical :: ok, fail
  integer :: ipt, n, i, nev, ifn, ih
  fail = .false.; errmax = 0
  cutoff = 1d-4
  call random_seed()
  do ih = 1, 4
  four_hard = merge(0.0_dp, 0.3_dp, ih == 1)
  four_hard2 = merge(0.3_dp, 0.0_dp, ih == 3)
  if (ih == 4) four_hard2 = 0.5_dp
  four_hard2mode = merge(1, 0, ih == 4)
  write(*,'(a,2f5.1,i2)') ' four_hard, four_hard2, four_hard2mode =', four_hard, four_hard2, four_hard2mode
  do ipt = 1, 6
     call random_number(rr)
     xB = 0.02_dp + 0.5_dp*rr(1)
     E = 1000*xB
     pB = [E, 0.0_dp, 0.0_dp, merge(E, -E, ipt <= 3)]
     call random_number(rr)
     th = pi*rr(1); ph = 2*pi*rr(2)
     pOB = (20 + 300*rr(3))*[1.0_dp, sin(th)*cos(ph), sin(th)*sin(ph), cos(th)]
     Q2 = 2*mdot(pB, pOB)
     q = pOB - pB
     ! 1. maps and splittings
     do i = 1, 100
        call random_number(r)
        call rambo3(pB/(xB + (1 - xB)*r(1)), q, r(2:7), t)
        a = t(:,1); b = t(:,2); c = t(:,3)
        pta = pB/(xB + (1 - xB)*r(1))
        call map_ff(a, b, c, ptij, ptk, y, z)
        call split_ff(ptij, ptk, y, z, azim(a, ptij, ptk), t(:,1), t(:,2), t(:,3))
        errmax(1) = max(errmax(1), maxval(abs(t(:,1) - a)), maxval(abs(t(:,2) - b)), &
             & maxval(abs(t(:,3) - c)))
        call map_fi(a, b, pta, ptij, ptk, x, z)
        errmax(2) = max(errmax(2), maxval(abs(ptij - ptk + c - q - (a + b + c - pta - q))))
        call split_fi(ptij, ptk, x, z, azim(a, ptij, ptk), t(:,1), t(:,2), t(:,3))
        errmax(1) = max(errmax(1), maxval(abs(t(:,1) - a)), maxval(abs(t(:,2) - b)), &
             & maxval(abs(t(:,3) - pta)))
        call map_if(a, b, pta, ptk, ptij, x, u)
        call split_if(ptk, ptij, x, u, azim(a, ptk, ptij), t(:,1), t(:,2), t(:,3))
        errmax(1) = max(errmax(1), maxval(abs(t(:,1) - a)), maxval(abs(t(:,2) - b)), &
             & maxval(abs(t(:,3) - pta)))
     enddo
     errmax(1) = errmax(1)/pB(0)
     ! 2. and 3.
     nev = 2000000
     s1 = 0; s2 = 0; m1 = 0; m2 = 0
     do i = 1, nev
        call random_number(r)
        call gen_four(pB, pOB, xB, r, cutoff, pa, k, w, ok)
        if (ok) then
           errmax(3) = max(errmax(3), abs(w/four_weight(pB, pOB, xB, pa, k, cutoff) - 1))
           errmax(4) = max(errmax(4), maxval(abs(k(:,1) + k(:,2) + k(:,3) - pa - q))/pB(0), &
                & abs(mdot(k(:,1),k(:,1)))/Q2, abs(mdot(k(:,2),k(:,2)))/Q2, &
                & abs(mdot(k(:,3),k(:,3)))/Q2, abs(mdot(pa,pa))/Q2)
           call testf(pa, k, Q2, f)
           s1 = s1 + w*f; m1 = m1 + (w*f)**2
        endif
        ! flat reference
        call random_number(r)
        u = xB**r(1)
        call rambo3(pB/u, q, r(2:7), t)
        call testf(pB/u, t, Q2, f)
        w = Q2/(2*pi)*log(1/xB)*Q2*(1 - u)/u/(256*pi**3)
        s2 = s2 + w*f; m2 = m2 + (w*f)**2
     enddo
     do ifn = 1, 2
        s1(ifn) = s1(ifn)/nev; s2(ifn) = s2(ifn)/nev
        sig1(ifn) = sqrt((m1(ifn)/nev - s1(ifn)**2)/nev)
        sig2(ifn) = sqrt((m2(ifn)/nev - s2(ifn)**2)/nev)
        pull = (s1(ifn) - s2(ifn))/sqrt(sig1(ifn)**2 + sig2(ifn)**2)
        write(*,'(a,i2,a,i2,a,es12.4,a,es9.2,a,es12.4,a,es9.2,a,f6.2)') ' point', ipt, ' f', ifn, ': gen_four ', &
             & s1(ifn), ' +- ', sig1(ifn), '  flat ', s2(ifn), ' +- ', sig2(ifn), '  pull ', pull
        if (abs(pull) > 4) fail = .true.
     enddo
  enddo
  enddo
  write(*,'(a,4es10.2)') ' max errors (maps, q, weight, momenta): ', errmax
  if (maxval(errmax) > 1d-8) fail = .true.
  if (fail) then
     write(*,*) 'test_four: FAILED'
     call exit(1)
  endif
  write(*,*) 'test_four: passed'
contains
  ! test functions with cuts on all invariants (relative to Q^2)
  subroutine testf(pa, k, Q2, f)
    real(dp), intent(in) :: pa(0:3), k(0:3,3), Q2
    real(dp), intent(out) :: f(2)
    real(dp) :: smin, delta = 2d-2
    smin = min(mdot(k(:,1),k(:,2)), mdot(k(:,1),k(:,3)), mdot(k(:,2),k(:,3)), &
         & mdot(pa,k(:,1)), mdot(pa,k(:,2)), mdot(pa,k(:,3)))*2/Q2
    f = 0
    if (smin < delta) return
    f(1) = 1
    f(2) = mdot(pa, k(:,1))/mdot(pa, k(:,1) + k(:,2) + k(:,3)) + 2*mdot(k(:,2),k(:,3))/Q2
  end subroutine testf

  ! azimuth of a around the (v1, v2) plane in the basis of dip_basis
  real(dp) function azim(a, v1, v2)
    real(dp), intent(in) :: a(0:3), v1(0:3), v2(0:3)
    real(dp) :: e1(0:3), e2(0:3)
    call dip_basis(v1, v2, e1, e2)
    azim = atan2(-mdot(a, e2), -mdot(a, e1))
  end function azim

  ! three massless momenta, flat in phase space, in the rest frame of
  ! P = pin + q, boosted to the frame of pin, q (RAMBO)
  subroutine rambo3(pin, q, r, t)
    real(dp), intent(in) :: pin(0:3), q(0:3), r(6)
    real(dp), intent(out) :: t(0:3,3)
    real(dp) :: ptot(0:3), W, qq(0:3,3), qsum(0:3), msum, bvec(3), gam, xx, aa, c, s, f, en, bq
    integer :: i
    ptot = pin + q
    W = sqrt(mdot(ptot, ptot))
    do i = 1, 3
       c = 2*r(2*i-1) - 1; s = sqrt(1 - c*c); f = 2*pi*r(2*i)
       ! energies from -log of a product: use two extra uniforms per particle
       en = -log(max(ranu(), 1d-300)*max(ranu(), 1d-300))
       qq(:,i) = en*[1.0_dp, s*cos(f), s*sin(f), c]
    enddo
    qsum = qq(:,1) + qq(:,2) + qq(:,3)
    msum = sqrt(mdot(qsum, qsum))
    bvec = -qsum(1:3)/msum
    gam = qsum(0)/msum
    aa = 1/(1 + gam)
    xx = W/msum
    do i = 1, 3
       bq = dot_product(bvec, qq(1:3,i))
       t(0,i) = xx*(gam*qq(0,i) + bq)
       t(1:3,i) = xx*(qq(1:3,i) + bvec*qq(0,i) + aa*bq*bvec)
    enddo
    ! boost from the P rest frame to the frame of P
    do i = 1, 3
       t(:,i) = boost(t(:,i), ptot)
    enddo
  end subroutine rambo3

  real(dp) function ranu()
    call random_number(ranu)
  end function ranu

  ! boost p (given in the rest frame of P) to the frame where P has its
  ! given momentum
  function boost(p, ptot) result(pb)
    real(dp), intent(in) :: p(0:3), ptot(0:3)
    real(dp) :: pb(0:3), m, bp
    m = sqrt(mdot(ptot, ptot))
    bp = dot_product(ptot(1:3), p(1:3))
    pb(0) = (ptot(0)*p(0) + bp)/m
    pb(1:3) = p(1:3) + ptot(1:3)*(bp/(m*(ptot(0) + m)) + p(0)/m)
  end function boost
end program test_four
