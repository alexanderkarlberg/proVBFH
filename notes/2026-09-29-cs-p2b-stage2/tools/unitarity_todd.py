#!/usr/bin/env python3
"""Absorptive (T-odd) part of the one-loop VBF H+3j line amplitude from
unitarity, independent of VBFNLO and NNLOJET.

For the radiating line q(pa) + J -> Q(pQ) + g(pg), J the current of the other
line, the only timelike channel is s_Qg; at one loop

    Abs A1 = 1/2 sum_{h', lam', colours} int dPhi_2(p', k')
             A_tree(q J -> Q'(p') g'(k')) A_tree(Q'(p') g'(k') -> Q g),

and the reflection-odd part of the virtual is
    T = 8 pi^2 * sum 2 Re(A0^* i Abs) / sum |A0|^2     [units alpha_s/(2 pi)]
(g_s = 1 in the amplitudes: A1/A0 ~ g_s^2 = 4 pi alpha_s).

The IR (collinear, Q' || Q) divergence of the cut is proportional to A0 with a
real coefficient and drops out of T; it is regulated by |t| > delta s and the
delta dependence is checked. The overall sign of the convention is fixed by
the LL helicities (where VBFNLO and NNLOJET agree); LR decides.

Input: heli_points.dat (harness_virt: p1, p2, k1, k2, k3 per point),
heli_vbfnlo.dat (VBFNLO V/B per placement and helicity).
"""
import numpy as np

# ---------------------------------------------------------------- Dirac algebra
s0 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
Z2 = np.zeros((2, 2), dtype=complex)
sig = [s0, sx, sy, sz]
sigb = [s0, -sx, -sy, -sz]
# Weyl basis: gamma^mu = [[0, sigma^mu], [sigmabar^mu, 0]], gamma5 = diag(-1,-1,1,1)
G = [np.block([[Z2, sig[m]], [sigb[m], Z2]]) for m in range(4)]
METRIC = np.array([1.0, -1.0, -1.0, -1.0])


def mdot(a, b):
    return a[0]*b[0] - a[1]*b[1] - a[2]*b[2] - a[3]*b[3]


def slash(v):
    """v^mu gamma_mu (v contravariant, possibly complex)"""
    return sum(METRIC[m]*v[m]*G[m] for m in range(4))


def spinor(p, h):
    """massless u(p) with helicity h = +1 (R) or -1 (L), positive energy"""
    E = p[0]
    n = np.array(p[1:])/np.linalg.norm(p[1:])
    th = np.arccos(np.clip(n[2], -1, 1))
    ph = np.arctan2(n[1], n[0])
    if h > 0:
        xi = np.array([np.cos(th/2), np.exp(1j*ph)*np.sin(th/2)])
        return np.sqrt(2*E)*np.concatenate([np.zeros(2), xi])
    xi = np.array([-np.exp(-1j*ph)*np.sin(th/2), np.cos(th/2)])
    return np.sqrt(2*E)*np.concatenate([xi, np.zeros(2)])


def ubar(u):
    return u.conj() @ G[0]


def polvecs(k):
    """two circular polarisation vectors (contravariant) for momentum k"""
    n = np.array(k[1:])/np.linalg.norm(k[1:])
    a = np.array([1.0, 0, 0]) if abs(n[0]) < 0.9 else np.array([0, 1.0, 0])
    e1 = a - n*np.dot(a, n); e1 /= np.linalg.norm(e1)
    e2 = np.cross(n, e1)
    return [np.concatenate([[0], (e1 + 1j*lam*e2)/np.sqrt(2)]) for lam in (+1, -1)]


# ---------------------------------------------------------------- colour
lam = np.zeros((8, 3, 3), dtype=complex)
lam[0] = [[0, 1, 0], [1, 0, 0], [0, 0, 0]]
lam[1] = [[0, -1j, 0], [1j, 0, 0], [0, 0, 0]]
lam[2] = [[1, 0, 0], [0, -1, 0], [0, 0, 0]]
lam[3] = [[0, 0, 1], [0, 0, 0], [1, 0, 0]]
lam[4] = [[0, 0, -1j], [0, 0, 0], [1j, 0, 0]]
lam[5] = [[0, 0, 0], [0, 0, 1], [0, 1, 0]]
lam[6] = [[0, 0, 0], [0, 0, -1j], [0, 1j, 0]]
lam[7] = np.diag([1, 1, -2])/np.sqrt(3)
T = lam/2
# f^{abc} = -2 i tr([T^a, T^b] T^c)
F = np.zeros((8, 8, 8))
for a in range(8):
    for b in range(8):
        comm = T[a] @ T[b] - T[b] @ T[a]
        for c in range(8):
            F[a, b, c] = np.real(-2j*np.trace(comm @ T[c]))
FT = np.einsum('abe,eij->abij', F, T)          # f^{a b e} T^e_{ij}
TT = np.einsum('aij,bjk->abik', T, T)          # (T^c T^c')_{ik}
TTr = np.einsum('bij,ajk->abik', T, T)         # (T^c' T^c)_{ik}
FTr = np.transpose(FT, (1, 0, 2, 3))            # f^{c' c e} T^e at [c, c']
TSIGN = -1.0   # M_t = -i f^{c' c e} T^e Nt/t from the Feynman rules; the gauge check tests it


# ---------------------------------------------------------------- amplitudes
def amp_qJ(pa, hq, pout, k, eps_star, J):
    """q(pa) + J -> Q(pout) + g(k, eps*), colour-stripped (colour T^c_{ij}):
    M = -ubar(pout)[eps*/ (pout+k)/ J/ / (pout+k)^2 + J/ (pa-k)/ eps*/ / (pa-k)^2] u(pa)"""
    u = spinor(pa, hq)
    ub = ubar(spinor(pout, hq))
    es, Js = slash(eps_star), slash(J)
    t1 = es @ slash(pout + k) @ Js/mdot(pout + k, pout + k)
    t2 = Js @ slash(pa - k) @ es/mdot(pa - k, pa - k)
    return -(ub @ (t1 + t2) @ u)


def amp_compton(pin, hq, kin, eps_in, pout, kout, eps_out_star):
    """Q(pin) + g'(kin, eps_in, c') -> Q(pout) + g(kout, eps*, c), full colour:
    returns A[c, c', i, i'] (i: outgoing quark colour, i': incoming)"""
    u = spinor(pin, hq)
    ub = ubar(spinor(pout, hq))
    ei, eo = slash(eps_in), slash(eps_out_star)
    s = mdot(pin + kin, pin + kin)
    uu = mdot(pin - kout, pin - kout)
    t = mdot(kin - kout, kin - kout)
    Ns = ub @ eo @ slash(pin + kin) @ ei @ u
    Nu = ub @ ei @ slash(pin - kout) @ eo @ u
    # three-gluon vertex, all momenta incoming: (c', nu, k1 = kin),
    # (c, mu, k2 = -kout), (e, sigma, k3 = kout - kin)
    k1, k2, k3 = kin, -kout, kout - kin
    # V^sigma = g^{nu mu}(k1-k2)^sigma + g^{mu sigma}(k2-k3)^nu + g^{sigma nu}(k3-k1)^mu
    V = (mdot(eps_in, eps_out_star)*(k1 - k2) + mdot(eps_in, k2 - k3)*eps_out_star
         + mdot(eps_out_star, k3 - k1)*eps_in)       # contravariant sigma index
    Nt = ub @ slash(V) @ u
    return -Ns/s*TT - Nu/uu*TTr + TSIGN*1j*Nt/t*FTr


def current(pin, pout, h):
    """J^mu = ubar(pout) gamma^mu u(pin) (contravariant components)"""
    ub, u = ubar(spinor(pout, h)), spinor(pin, h)
    return np.array([ub @ G[m] @ u for m in range(4)])


# ---------------------------------------------------------------- checks
def gauge_checks(pa, pQ, pg, J):
    """tree amplitudes vanish for eps -> momentum"""
    r1 = abs(amp_qJ(pa, -1, pQ, pg, pg.astype(complex), J))/abs(amp_qJ(pa, -1, pQ, pg, polvecs(pg)[0].conj(), J))
    # Compton with the outgoing gluon's eps* -> its momentum
    P = pQ + pg
    # a physical 2 -> 2 configuration: rotate in the rest frame
    pin, kk = cut_momenta(P, pQ, 0.7, 1.3)
    A1 = amp_compton(pin, -1, kk, polvecs(kk)[0], pQ, pg, pg.astype(complex))
    A0 = amp_compton(pin, -1, kk, polvecs(kk)[0], pQ, pg, polvecs(pg)[0].conj())
    r2 = np.max(np.abs(A1))/np.max(np.abs(A0))
    return r1, r2


def compton_check(P, pQ):
    """sum |M|^2/96 for q g -> q g against -(4/9)(s^2+u^2)/(s u) + (s^2+u^2)/t^2"""
    pin, kin = cut_momenta(P, pQ, 0.3, 2.1)
    pg = P - pQ
    tot = 0.0
    for h in (-1, 1):
        for ei in polvecs(kin):
            for eo in polvecs(pg):
                A = amp_compton(pin, h, kin, ei, pQ, pg, eo.conj())
                tot += np.sum(np.abs(A)**2)
    s = mdot(P, P); t = mdot(kin - pg, kin - pg); u = mdot(pin - pg, pin - pg)
    return tot/96, -(4/9)*(s*s + u*u)/(s*u) + (s*s + u*u)/t**2


# ---------------------------------------------------------------- cut kinematics
def boost(p, b):
    """boost p by velocity vector b"""
    b2 = np.dot(b, b)
    if b2 == 0:
        return p.copy()
    g = 1/np.sqrt(1 - b2)
    bp = np.dot(b, p[1:])
    E = g*(p[0] + bp)
    v = p[1:] + ((g - 1)*bp/b2 + g*p[0])*b
    return np.concatenate([[E], v])


def cut_momenta(P, pQ, cth, phi):
    """p', k' with p' + k' = P; direction of p' in the rest frame of P at
    polar angle arccos(cth) from pQ's direction, azimuth phi"""
    b = P[1:]/P[0]
    M = np.sqrt(mdot(P, P))
    q = boost(pQ, -b)
    z = q[1:]/np.linalg.norm(q[1:])
    a = np.array([1.0, 0, 0]) if abs(z[0]) < 0.9 else np.array([0, 1.0, 0])
    x = a - z*np.dot(a, z); x /= np.linalg.norm(x)
    y = np.cross(z, x)
    sth = np.sqrt(max(0.0, 1 - cth**2))
    n = cth*z + sth*(np.cos(phi)*x + np.sin(phi)*y)
    pp = np.concatenate([[M/2], M/2*n])
    kp = np.concatenate([[M/2], -M/2*n])
    return boost(pp, b), boost(kp, b)


# ---------------------------------------------------------------- T-odd part
def todd(pa, pQ, pg, J, hq, delta, ncos=48, nphi=32):
    """T = 8 pi^2 sum 2 Re(A0^* i Abs) / sum |A0|^2 for the radiating line with
    helicity hq and external current J; |t| > delta s"""
    P = pQ + pg
    s = mdot(P, P)
    # tree, summed over the gluon's polarisations and colours
    eg = polvecs(pg)
    A0 = [amp_qJ(pa, hq, pQ, pg, e.conj(), J) for e in eg]     # colour T^c_{ij}
    # colour sum of |A0|^2: tr(T^c T^c) = 4
    B = sum(abs(a)**2 for a in A0)*4.0
    # cut integral: t = -s (1 - cth)/2 >= delta s  ->  1 - cth >= 2 delta.
    # Panels in cth: [-1, 0] with 1 + cth = w^2 (u-channel), [0, 1 - 2 delta]
    # with ln(1 - cth) uniform (t-channel)
    xg, wg = np.polynomial.legendre.leggauss(ncos)
    nodes, weights = [], []
    for x, w in zip(xg, wg):
        v = (x + 1)/2                                   # [0,1]
        nodes.append(v**2 - 1); weights.append(w/2*2*v)  # cth = w^2 - 1, w = v
        lo, hi = np.log(1.0), np.log(2*delta)            # 1 - cth from 1 to 2 delta
        l = lo + (hi - lo)*v
        c = 1 - np.exp(l)
        nodes.append(c); weights.append(w/2*abs(hi - lo)*np.exp(l))
    phis = 2*np.pi*(np.arange(nphi) + 0.5)/nphi
    num = 0.0
    for cth, wc in zip(nodes, weights):
        for ph in phis:
            pp, kp = cut_momenta(P, pQ, cth, ph)
            dphi2 = wc*(2*np.pi/nphi)/(32*np.pi**2)      # dPhi_2 = dOmega/(32 pi^2)
            ekp = polvecs(kp)
            for lamg, e_out in enumerate(eg):
                abs_amp = np.zeros((8, 3, 3), dtype=complex)   # [c, i, j]
                for e_int in ekp:
                    AL = amp_qJ(pa, hq, pp, kp, e_int.conj(), J)          # colour T^{c'}_{i'j}
                    AR = amp_compton(pp, hq, kp, e_int, pQ, pg, e_out.conj())   # [c, c', i, i']
                    abs_amp += 0.5*AL*np.einsum('xyik,ykj->xij', AR, T)
                # interference with A0: sum_c,i,j conj(A0 T^c_ij) abs[c,i,j]
                num += dphi2*2*np.real(np.conj(A0[lamg])*1j*np.einsum('cij,cij->', T.conj(), abs_amp))
    return 8*np.pi**2*num/B


# ---------------------------------------------------------------- two-chart version
# The cut integrand has, besides the t-channel collinear log at p' || pQ, an
# integrable point singularity ~ e^{i phi}/theta where the cut gluon k' becomes
# collinear to the incoming quark pa (AL ~ 1/sqrt(pa.k')). A tensor Gauss rule
# around the pQ axis converges slowly there, so split the sphere with a
# partition of unity: chart 1 (polar axis pQ) carries pi1 = da/(da + dQ), which
# vanishes like theta_a^2 at the ISR point; chart 2 (polar axis at the ISR point,
# Gauss in theta) carries pi2 = dQ/(da + dQ), which vanishes like sin^2 theta_Q
# on the pQ axis (both ends), so needs no cutoff.
def _frame(z):
    a = np.array([1.0, 0, 0]) if abs(z[0]) < 0.9 else np.array([0, 1.0, 0])
    x = a - z*np.dot(a, z); x /= np.linalg.norm(x)
    return x, np.cross(z, x)


def todd2(pa, pQ, pg, J, hq, delta, n1=(48, 32), n2=(48, 32)):
    P = pQ + pg
    b = P[1:]/P[0]
    M = np.sqrt(mdot(P, P))
    zQ = boost(pQ, -b)[1:]; zQ /= np.linalg.norm(zQ)
    za = -boost(pa, -b)[1:]; za /= np.linalg.norm(za)     # p' direction with k' || pa
    eg = polvecs(pg)
    A0 = [amp_qJ(pa, hq, pQ, pg, e.conj(), J) for e in eg]
    B = sum(abs(a)**2 for a in A0)*4.0
    TC = T.conj()

    def integrand(n):
        pp = boost(np.concatenate([[M/2], M/2*n]), b)
        kp = boost(np.concatenate([[M/2], -M/2*n]), b)
        ekp = polvecs(kp)
        ALs = [[amp_qJ(pa, hq, pp, kp, e.conj(), J) for e in ekp]]
        val = 0.0
        for lamg, e_out in enumerate(eg):
            ab = np.zeros((8, 3, 3), dtype=complex)
            for ii, e_int in enumerate(ekp):
                AR = amp_compton(pp, hq, kp, e_int, pQ, pg, e_out.conj())
                ab += 0.5*ALs[0][ii]*np.einsum('xyik,ykj->xij', AR, T)
            val += 2*np.real(np.conj(A0[lamg])*1j*np.einsum('cij,cij->', TC, ab))
        return val

    def pis(n):
        da = 1 - np.dot(n, za)
        dQ = 1 - np.dot(n, zQ)**2
        return da/(da + dQ), dQ/(da + dQ)

    xg, wg = np.polynomial.legendre.leggauss(n1[0])
    tot = 0.0
    # chart 1: axis zQ, panels as in todd
    x1, y1 = _frame(zQ)
    for x, w in zip(xg, wg):
        v = (x + 1)/2
        lo, hi = 0.0, np.log(2*delta)
        l = lo + (hi - lo)*v
        for cth, wc in ((v**2 - 1, w/2*2*v), (1 - np.exp(l), w/2*abs(hi - lo)*np.exp(l))):
            sth = np.sqrt(max(0.0, 1 - cth**2))
            for ph in 2*np.pi*(np.arange(n1[1]) + 0.5)/n1[1]:
                n = cth*zQ + sth*(np.cos(ph)*x1 + np.sin(ph)*y1)
                tot += wc*(2*np.pi/n1[1])*pis(n)[0]*integrand(n)
    # chart 2: axis za, Gauss in theta
    x2, y2 = _frame(za)
    xg2, wg2 = np.polynomial.legendre.leggauss(n2[0])
    for x, w in zip(xg2, wg2):
        th = np.pi*(x + 1)/2
        wt = w*np.pi/2*np.sin(th)
        for ph in 2*np.pi*(np.arange(n2[1]) + 0.5)/n2[1]:
            n = np.cos(th)*za + np.sin(th)*(np.cos(ph)*x2 + np.sin(ph)*y2)
            tot += wt*(2*np.pi/n2[1])*pis(n)[1]*integrand(n)
    return 8*np.pi**2*tot/(32*np.pi**2)/B


if __name__ == '__main__':
    pts = np.loadtxt('heli_points.dat')
    vb = np.loadtxt('heli_vbfnlo.dat')
    for row in pts[:2]:
        ipt = int(row[0])
        p1, p2, k1, k2, k3 = (row[1+4*i:5+4*i] for i in range(5))
        print(f'point {ipt}: s_Qg line 1 = {mdot(k1+k3,k1+k3):.4e}, momentum conservation '
              f'{np.max(np.abs(p1 + p2 - k1 - k2 - k3)):.2e} (the Higgs carries the rest)')
        # gauge checks
        J2 = current(p2, k2, -1)
        r1, r2 = gauge_checks(p1, k1, k3, J2)
        print(f'  gauge checks (should be ~1e-15): q J -> Q g {r1:.1e}, Compton {r2:.1e}')
        print('  q g -> q g |M|^2 check (code, formula):', compton_check(k1 + k3, k1))
        for hname, (h1, h2) in (('LL', (-1, -1)), ('LR', (-1, +1)), ('RL', (+1, -1)), ('RR', (+1, +1))):
            # placement 1: gluon on line 1 (p1 -> k1), current of line 2
            J2 = current(p2, k2, h2)
            J1 = current(p1, k1, h1)
            res = []
            for delta in (1e-5, 1e-7):
                t1 = todd2(p1, k1, k3, J2, h1, delta)
                t2 = todd2(p2, k2, k3, J1, h2, delta)
                res.append((t1, t2))
            # VBFNLO: T = [V(h1,h2) - V(-h1,-h2)]/2 per placement
            ih1, ih2 = (1 if h1 < 0 else 2), (1 if h2 < 0 else 2)
            v = vb[(vb[:, 0] == ipt) & (vb[:, 1] == ih1) & (vb[:, 2] == ih2)][0]
            vf = vb[(vb[:, 0] == ipt) & (vb[:, 1] == 3 - ih1) & (vb[:, 2] == 3 - ih2)][0]
            tv1, tv2 = (v[3] - vf[3])/2, (v[4] - vf[4])/2
            print(f'  {hname}: unitarity T (placement 1, 2) delta 1e-5: {res[0][0]:+.5f} {res[0][1]:+.5f}'
                  f'  delta 1e-7: {res[1][0]:+.5f} {res[1][1]:+.5f}   VBFNLO: {tv1:+.5f} {tv2:+.5f}')
