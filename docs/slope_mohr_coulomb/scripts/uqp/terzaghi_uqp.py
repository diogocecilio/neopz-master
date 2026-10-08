#!/usr/bin/env python3
"""
1D sanity check of the three-field (u-q-p) Biot formulation against the two-field (u-p) one
(the formulation of Material/Plasticity/TPZMatPoroElastoPlasticUP.cpp) and against the Terzaghi series.

Conventions (as in the repository):
  sigma = sigma' - alpha p,  tension positive, p positive in compression, backward Euler,
  mobility k, Darcy  q = -k (grad p - rho_w g),  mass balance  alpha div u_t + (1/M) p_t + div q = 0.

u-p  (code):   [ K    -Qc        ] [U]   [ f_t                         ]
               [ Qc^T  S+dt H    ] [P] = [ Qc^T U_n + S P_n + dt f_g   ]      (p = p_D essential)

u-q-p (new):   [ K     0    -Qc  ] [U]   [ f_t                         ]
               [ 0     A    -B^T ] [Q] = [ g_q - h_pD                   ]      (q.n = q_N essential, p_D natural)
               [ Qc^T  dt B  S   ] [P]   [ Qc^T U_n + S P_n             ]

  K=int E N_u' N_u', Qc=int alpha N_u' N_p, S=int (1/M) N_p N_p, H=int k N_p' N_p', f_g=int k N_p' (rho_w g),
  A=int (1/k) N_q N_q, B_ij = int N_p,i N_q,j', g_q=int rho_w g N_q, h_pD = p_D N_q(L) (n=+1 at x=L).

Column 0<x<1: u(0)=0, impermeable bottom (q(0)=0), top x=1: traction -sigma0, p_D = 0 (drained).
In 1D H(div) = H1, so the Raviart-Thomas-like space of order k is continuous P_k.

Run:  python3 terzaghi_uqp.py
"""
import numpy as np

# ----------------------------------------------------------------------------------------------
# tiny 1D finite element machinery
# ----------------------------------------------------------------------------------------------
XG, WG = np.polynomial.legendre.leggauss(6)
XG = 0.5 * (XG + 1.0)
WG = 0.5 * WG


def lagrange(order, xi):
    """values and derivatives (w.r.t. xi in [0,1]) of the equispaced Lagrange basis of given order"""
    nodes = np.linspace(0.0, 1.0, order + 1)
    n = order + 1
    phi = np.ones(n)
    dphi = np.zeros(n)
    for a in range(n):
        for b in range(n):
            if b != a:
                phi[a] *= (xi - nodes[b]) / (nodes[a] - nodes[b])
        for c in range(n):
            if c == a:
                continue
            term = 1.0 / (nodes[a] - nodes[c])
            for b in range(n):
                if b != a and b != c:
                    term *= (xi - nodes[b]) / (nodes[a] - nodes[b])
            dphi[a] += term
    return phi, dphi


class Space:
    """continuous Lagrange space of order k (kind 'cg') or piecewise constant (kind 'dg0') on n elements"""

    def __init__(self, n, kind, order):
        self.n, self.kind, self.order = n, kind, order
        if kind == 'cg':
            self.ndof = n * order + 1
            self.loc = [np.arange(order + 1) + e * order for e in range(n)]
        else:
            self.ndof = n
            self.order = 0
            self.loc = [np.array([e]) for e in range(n)]

    def eval(self, e, xi, h):
        """phi, dphi/dx at the reference point xi of element e (size h)"""
        if self.kind == 'cg':
            phi, dphi = lagrange(self.order, xi)
            return phi, dphi / h
        return np.array([1.0]), np.array([0.0])

    def coords(self, h):
        if self.kind == 'cg':
            return np.linspace(0.0, 1.0, self.ndof)
        return (np.arange(self.n) + 0.5) * h


def bilinear(sa, da, sb, db, coef, n):
    """M[i,j] = int coef * d^da(N_a,i) d^db(N_b,j) ,  d=0 value, d=1 derivative"""
    h = 1.0 / n
    M = np.zeros((sa.ndof, sb.ndof))
    for e in range(n):
        Me = np.zeros((len(sa.loc[e]), len(sb.loc[e])))
        for xi, w in zip(XG, WG):
            pa, dpa = sa.eval(e, xi, h)
            pb, dpb = sb.eval(e, xi, h)
            fa = pa if da == 0 else dpa
            fb = pb if db == 0 else dpb
            Me += w * h * coef * np.outer(fa, fb)
        M[np.ix_(sa.loc[e], sb.loc[e])] += Me
    return M


def linear(sa, da, coef, n):
    h = 1.0 / n
    v = np.zeros(sa.ndof)
    for e in range(n):
        for xi, w in zip(XG, WG):
            pa, dpa = sa.eval(e, xi, h)
            v[sa.loc[e]] += w * h * coef * (pa if da == 0 else dpa)
    return v


# ----------------------------------------------------------------------------------------------
# analytical solution
# ----------------------------------------------------------------------------------------------
def cv_of(par):
    return par['k'] / (par['invM'] + par['alpha'] ** 2 / par['E'])


def p0_of(par):
    a, iM, E = par['alpha'], par['invM'], par['E']
    return a * par['sig0'] / (a * a + E * iM)


def terzaghi(x, T, par, nterms=40000):
    """p, q, u_top of the Terzaghi column (x from the impermeable bottom, drained top at x=1)"""
    m = np.arange(nterms)
    Mm = (2 * m + 1) * np.pi / 2
    p0 = p0_of(par)
    ex = np.exp(-Mm ** 2 * T)
    x = np.atleast_1d(x)
    z = 1.0 - x
    p = (2 * p0 / Mm * ex)[None, :] * np.sin(np.outer(z, Mm))
    q = par['k'] * (2 * p0 * ex)[None, :] * np.cos(np.outer(z, Mm))  # q = -k dp/dx
    pbar = np.sum(2 * p0 / Mm ** 2 * ex)
    utop = (-par['sig0'] + par['alpha'] * pbar) / par['E']
    return p.sum(1), q.sum(1), utop


# ----------------------------------------------------------------------------------------------
# the two discretisations
# ----------------------------------------------------------------------------------------------
def build(n, par, uord=2, qord=2, pord=1, pkind='cg'):
    su = Space(n, 'cg', uord)
    sp = Space(n, pkind, pord)
    sq = Space(n, 'cg', qord)
    E, al, iM, k = par['E'], par['alpha'], par['invM'], par['k']
    rg = -par['gammaw']  # rho_w g (x points upwards)
    ops = {'su': su, 'sp': sp, 'sq': sq}
    ops['K'] = bilinear(su, 1, su, 1, E, n)
    ops['Qc'] = bilinear(su, 1, sp, 0, al, n)
    ops['S'] = bilinear(sp, 0, sp, 0, iM, n)
    ops['A'] = bilinear(sq, 0, sq, 0, 1.0 / k, n)
    ops['B'] = bilinear(sp, 0, sq, 1, 1.0, n)            # (np x nq):  int N_p div phi_q
    ops['gq'] = linear(sq, 0, rg, n)                     # int rho_w g . phi_q
    if pkind == 'cg':
        ops['H'] = bilinear(sp, 1, sp, 1, k, n)
        ops['fg'] = linear(sp, 1, k * rg, n)             # int k grad N_p . rho_w g
    ft = np.zeros(su.ndof)
    ft[-1] = -par['sig0']                                # traction t = sigma n,  sigma = -sig0,  n=+1
    ops['ft'] = ft
    return ops


def solve_bc(Kf, rhs, cons):
    """solve Kf x = rhs with x[i]=val for i in cons (dict)"""
    N = Kf.shape[0]
    c = np.array(sorted(cons.keys()), dtype=int)
    f = np.setdiff1d(np.arange(N), c)
    x = np.zeros(N)
    x[c] = [cons[i] for i in c]
    x[f] = np.linalg.solve(Kf[np.ix_(f, f)], rhs[f] - Kf[np.ix_(f, c)] @ x[c])
    return x


def system(ops, form, dt, par):
    """returns (K, F, Hist, cons, slices) of one backward-Euler step"""
    nu, nq, npp = ops['su'].ndof, ops['sq'].ndof, ops['sp'].ndof
    K, Qc, S = ops['K'], ops['Qc'], ops['S']
    pD = par['pD']
    if form == 'up':
        N = nu + npp
        M = np.zeros((N, N))
        M[:nu, :nu] = K
        M[:nu, nu:] = -Qc
        M[nu:, :nu] = Qc.T
        M[nu:, nu:] = S + dt * ops['H']
        F = np.zeros(N)
        F[:nu] = ops['ft']
        F[nu:] = dt * ops['fg']
        Hist = np.zeros((N, N))
        Hist[nu:, :nu] = Qc.T
        Hist[nu:, nu:] = S
        cons = {0: 0.0, nu + npp - 1: pD}                # u(0)=0, p(L)=p_D  (essential)
        sl = (slice(0, nu), None, slice(nu, N))
    else:
        N = nu + nq + npp
        M = np.zeros((N, N))
        M[:nu, :nu] = K
        M[:nu, nu + nq:] = -Qc
        M[nu:nu + nq, nu:nu + nq] = ops['A']
        M[nu:nu + nq, nu + nq:] = -ops['B'].T
        M[nu + nq:, :nu] = Qc.T
        M[nu + nq:, nu:nu + nq] = dt * ops['B']
        M[nu + nq:, nu + nq:] = S
        F = np.zeros(N)
        F[:nu] = ops['ft']
        hpD = np.zeros(nq)
        hpD[-1] = pD                                     # int_{Gamma_p} p_D r.n,  n=+1 at x=1
        F[nu:nu + nq] = ops['gq'] - hpD
        Hist = np.zeros((N, N))
        Hist[nu + nq:, :nu] = Qc.T
        Hist[nu + nq:, nu + nq:] = S
        cons = {0: 0.0, nu: 0.0}                         # u(0)=0, q(0)=0 (q.n=0 essential); p_D is NATURAL
        sl = (slice(0, nu), slice(nu, nu + nq), slice(nu + nq, N))
    return M, F, Hist, cons, sl


def march(n, par, form, schedule, report, uord=2, qord=2, pord=1, pkind='cg'):
    """schedule = [(dt, nsteps), ...]; report = list of times at which results are stored"""
    ops = build(n, par, uord, qord, pord, pkind)
    nu, nq, npp = ops['su'].ndof, ops['sq'].ndof, ops['sp'].ndof
    # undrained step (dt = 0): instantaneous load, as the driver does
    M, F, Hist, cons, sl = system(ops, form, 0.0, par)
    N = M.shape[0]
    # history of the undrained step = the state BEFORE the load: U_n = 0, P_n = 0 (S P_n matters only for S>0)
    Xn = np.zeros(N)
    X = solve_bc(M, F + Hist @ Xn, cons)
    t = 0.0
    out = {}
    rep = sorted(report)
    ir = 0
    for dt, ns in schedule:
        M, F, Hist, cons, sl = system(ops, form, dt, par)
        f = np.setdiff1d(np.arange(N), np.array(sorted(cons.keys()), dtype=int))
        Minv = np.linalg.inv(M[np.ix_(f, f)])
        c = np.array(sorted(cons.keys()), dtype=int)
        xc = np.array([cons[i] for i in c])
        for _ in range(ns):
            rhs = F + Hist @ X
            Xnew = np.zeros(N)
            Xnew[c] = xc
            Xnew[f] = Minv @ (rhs[f] - M[np.ix_(f, c)] @ xc)
            Xprev = X
            X = Xnew
            t += dt
            while ir < len(rep) and abs(t - rep[ir]) < 0.5 * dt * 1.0001:
                out[rep[ir]] = (X.copy(), ops, sl, form, Xprev.copy(), dt)
                ir += 1
    return out


def midpoint_q(X, ops, sl, form, par):
    """flux at the element midpoints: -k p' (u-p, post-processed) or q_h (u-q-p)"""
    n = ops['su'].n
    h = 1.0 / n
    res = np.zeros(n)
    sp, sq = ops['sp'], ops['sq']
    for e in range(n):
        if form == 'up':
            ph, dph = sp.eval(e, 0.5, h)
            dp = dph @ X[sl[2]][sp.loc[e]]
            res[e] = -par['k'] * (dp - (-par['gammaw']))      # q = -k (p' - rho_w g)
        else:
            ph, dph = sq.eval(e, 0.5, h)
            res[e] = ph @ X[sl[1]][sq.loc[e]]
    return res


def errors(out, par, T, n):
    cv = cv_of(par)
    t = T / cv
    key = min(out.keys(), key=lambda s: abs(s - t))
    X, ops, sl, form, Xprev, dtl = out[key]
    sp = ops['sp']
    xp = sp.coords(1.0 / n)
    pex, _, utop_ex = terzaghi(xp, T, par)
    ph = X[sl[2]]
    xm = (np.arange(n) + 0.5) / n
    _, qex, _ = terzaghi(xm, T, par)
    qh = midpoint_q(X, ops, sl, form, par)
    ep = np.max(np.abs(ph - pex)) / p0_of(par)
    eu = abs(X[sl[0]][-1] - utop_ex) / abs(utop_ex) if abs(utop_ex) > 1e-12 else abs(X[sl[0]][-1] - utop_ex)
    eq = np.max(np.abs(qh - qex)) / np.max(np.abs(qex))
    # discharge through the drained face versus the rate of volume change (global mass balance, alpha div u_t + (1/M) p_t + div q = 0
    # integrated over the column, q(0)=0):  q(1) = -alpha d/dt u(1) - (1/M) d/dt int p
    du = (X[sl[0]][-1] - Xprev[sl[0]][-1]) / dtl
    pint = lambda Y: np.sum(linear(ops['sp'], 0, 1.0, n) * Y[sl[2]])
    rate = -par['alpha'] * du - par['invM'] * (pint(X) - pint(Xprev)) / dtl
    if form == 'uqp':
        qtop = X[sl[1]][-1]                     # the flux dof at x=1 (normal trace)
    else:
        qtop = midpoint_q(X, ops, sl, form, par)[-1]   # u-p: only -k grad p_h of the last element exists
    eqtop = abs(qtop - rate) / abs(rate)
    return ep, eu, eq, eqtop, ph, pex


def main():
    base = dict(E=1.0, k=1.0, alpha=1.0, invM=0.0, sig0=1.0, gammaw=0.0, pD=0.0)
    gen = dict(E=1.0, k=1.0, alpha=0.9, invM=0.5, sig0=1.0, gammaw=0.0, pD=0.0)
    results = []
    for name, par in (('alpha=1, 1/M=0', base), ('alpha=0.9, 1/M=0.5', gen)):
        cv = cv_of(par)
        print('=' * 100)
        print('Parameter set: %s  (cv = %.5f, p0 = %.5f)' % (name, cv, p0_of(par)))
        for n in (10, 20, 40):
            dt = 2.0e-4 / cv if n == 40 else 4.0e-4 / cv
            Ts = (0.1, 1.0)
            t_end = max(Ts) / cv
            ns = int(round(t_end / dt))
            sched = [(dt, ns)]
            rep = [T / cv for T in Ts]
            outs = {}
            configs = [('u-p  (u P2, p P1)', 'up', dict(uord=2, pord=1)),
                       ('u-q-p (u P2, q P2, p P1)  [div Q = P1 >= P_h]', 'uqp', dict(uord=2, qord=2, pord=1)),
                       ('u-q-p (u P2, q P1, p P1)  [div Q = P0 < P_h]', 'uqp', dict(uord=2, qord=1, pord=1)),
                       ('u-q-p (u P2, q P1, p P0)  [RT0-like]', 'uqp', dict(uord=2, qord=1, pord=0, pkind='dg0')),
                       ]
            print('-' * 100)
            print('n = %d elements, dt = %.3e (T-step %.2e)' % (n, dt, dt * cv))
            print('%-52s %-5s %10s %10s %10s %12s' % ('scheme', 'T', 'max|dp|/p0', '|du_top|rel', 'max|dq|rel', 'massbal.defect'))
            for label, form, kw in configs:
                out = march(n, par, form, sched, rep, **kw)
                for T in Ts:
                    ep, eu, eq, eqt, ph, pex = errors(out, par, T, n)
                    print('%-52s %-5g %10.3e %10.3e %10.3e %12.3e' % (label, T, ep, eu, eq, eqt))
                    results.append((name, n, label, T, ep, eu, eq, eqt))
    return results


def schur_checks(n=12):
    """(a) algebraic elimination of Q from the u-q-p system (Schur complement) reproduces the full solution:
           S + dt H_mix,  H_mix = B A^-1 B^T ,  f_g,mix = B A^-1 (h_pD - g_q)    [max |diff| of the solutions]
       (b) H_up - H_mix (on p(L)=0) is positive semidefinite: the continuous-normal-flux space cannot represent
           the discontinuous -k grad p_h exactly, so H_mix <= H_up (energy of the L2-projection of -k grad p_h)."""
    par = dict(E=1.0, k=0.7, alpha=0.9, invM=0.3, sig0=1.0, gammaw=1.5, pD=0.4)
    dt = 0.01
    ops = build(n, par, 2, 2, 1)
    M, F, Hist, cons, sl = system(ops, 'uqp', dt, par)
    N = M.shape[0]
    Xn = np.random.RandomState(1).rand(N)
    full = solve_bc(M, F + Hist @ Xn, cons)
    nu, nq, npp = ops['su'].ndof, ops['sq'].ndof, ops['sp'].ndof
    fq = np.arange(1, nq)
    B = ops['B'][:, fq]
    A = ops['A'][np.ix_(fq, fq)]
    Hm = B @ np.linalg.solve(A, B.T)
    hpD = np.zeros(nq)
    hpD[-1] = par['pD']
    fgm = B @ np.linalg.solve(A, (hpD - ops['gq'])[fq])
    # reduced (U,P) system: [K -Qc ; Qc^T  S + dt H_mix] [U;P] = [f_t ; Qc^T U_n + S P_n + dt f_g,mix]
    Kr = np.zeros((nu + npp, nu + npp))
    Kr[:nu, :nu] = ops['K']
    Kr[:nu, nu:] = -ops['Qc']
    Kr[nu:, :nu] = ops['Qc'].T
    Kr[nu:, nu:] = ops['S'] + dt * Hm
    Fr = np.concatenate([ops['ft'], dt * fgm + ops['Qc'].T @ Xn[:nu] + ops['S'] @ Xn[nu + nq:]])
    red = solve_bc(Kr, Fr, {0: 0.0})
    d1 = max(np.max(np.abs(red[:nu] - full[:nu])), np.max(np.abs(red[nu:] - full[nu + nq:])))
    # (b)
    Hup = ops['H']
    pf = np.arange(0, npp - 1)
    ev = np.linalg.eigvalsh(Hup[np.ix_(pf, pf)] - Hm[np.ix_(pf, pf)])
    rel = np.linalg.norm(Hup[np.ix_(pf, pf)] - Hm[np.ix_(pf, pf)]) / np.linalg.norm(Hup[np.ix_(pf, pf)])
    return d1, ev.min(), rel


def hydrostatic(n=20):
    """gravity: after a long time both forms give p = gamma_w (1-x), q = 0 (checks the sign of rho_w g)"""
    par = dict(E=1.0, k=1.0, alpha=1.0, invM=0.0, sig0=1.0, gammaw=2.0, pD=0.0)
    sched = [(1e-3, 300), (1e-2, 200), (1.0, 40)]
    tend = 1e-3 * 300 + 1e-2 * 200 + 40.0
    res = {}
    for form in ('up', 'uqp'):
        out = march(n, par, form, sched, [tend])
        X, ops, sl, fm = out[tend][:4]
        xp = ops['sp'].coords(1.0 / n)
        pex = par['gammaw'] * (1.0 - xp)
        ep = np.max(np.abs(X[sl[2]] - pex))
        qm = np.max(np.abs(midpoint_q(X, ops, sl, fm, par)))
        res[form] = (ep, qm)
    return res


def inhomogeneous_pD(n=20):
    """p_D = 0.5 on the top (natural BC in the u-q-p form): the steady state is p = p_D, q = 0 (no gravity)"""
    par = dict(E=1.0, k=1.0, alpha=1.0, invM=0.0, sig0=1.0, gammaw=0.0, pD=0.5)
    sched = [(1e-3, 300), (1e-2, 200), (1.0, 40)]
    tend = 1e-3 * 300 + 1e-2 * 200 + 40.0
    res = {}
    for form in ('up', 'uqp'):
        out = march(n, par, form, sched, [tend])
        X, ops, sl, fm = out[tend][:4]
        ep = np.max(np.abs(X[sl[2]] - 0.5))
        res[form] = ep
    return res


def early_time(n=20):
    """first 10 steps after the load (T = 10 dt), dt well below h^2/(6 cv) (Vermeer-Verruijt limit): oscillations and weak p_D"""
    par = dict(E=1.0, k=1.0, alpha=1.0, invM=0.0, sig0=1.0, gammaw=0.0, pD=0.0)
    print('###  early time (n=%d, h^2/(6 cv) = %.1e): min p, max p, p at the drained node, max|p-p_exact|' % (n, (1.0 / n) ** 2 / 6))
    for dt in (1e-4, 1e-5):
        for label, form, kw in [('u-p  (u2,p1)', 'up', dict()), ('u-q-p (u2,q2,p1)', 'uqp', dict(qord=2)),
                                ('u-q-p (u2,q1,p1)', 'uqp', dict(qord=1)),
                                ('u-q-p (u2,q1,p0)', 'uqp', dict(qord=1, pord=0, pkind='dg0'))]:
            T = 10 * dt
            out = march(n, par, form, [(dt, 10)], [T], **kw)
            X, ops, sl, fm = out[min(out)][:4]
            p = X[sl[2]]
            pex, _, _ = terzaghi(ops['sp'].coords(1.0 / n), T, par)
            print('   dt=%.0e %-18s min=%6.3f max=%6.3f p(top)=%6.3f err=%.2e' % (dt, label, p.min(), p.max(), p[-1], np.max(np.abs(p - pex))))


# ----------------------------------------------------------------------------------------------
# plastic 1D check of the consistent tangent (3-field Newton)
# ----------------------------------------------------------------------------------------------
def plastic_newton(n=10, consistent=True, dt=0.05, nsteps=6, verbose=True):
    """1D elastoplastic column, J2-like yield |s'| <= sy(kappa), sy = sinf-(sinf-sy0) exp(-beta kappa).
    Full 3-field Newton (U,Q,P), memory at the Gauss points updated after convergence."""
    E, k, al, iM = 1.0, 1.0, 1.0, 0.0
    sy0, sinf, beta = 0.4, 0.9, 6.0
    sig0 = 1.0
    su, sq, sp = Space(n, 'cg', 2), Space(n, 'cg', 2), Space(n, 'cg', 1)
    nu, nq, npp = su.ndof, sq.ndof, sp.ndof
    N = nu + nq + npp
    h = 1.0 / n
    Qc = bilinear(su, 1, sp, 0, al, n)
    A = bilinear(sq, 0, sq, 0, 1.0 / k, n)
    B = bilinear(sp, 0, sq, 1, 1.0, n)
    S = bilinear(sp, 0, sp, 0, iM, n)
    ft = np.zeros(nu)
    ft[-1] = -sig0
    ng = len(XG)
    epn = np.zeros((n, ng))
    kapn = np.zeros((n, ng))

    def sy(kap):
        return sinf - (sinf - sy0) * np.exp(-beta * kap)

    def dsy(kap):
        return beta * (sinf - sy0) * np.exp(-beta * kap)

    def stress(eps, ep_n, kap_n):
        st = E * (eps - ep_n)
        f = abs(st) - sy(kap_n)
        if f <= 0.0:
            return st, E, ep_n, kap_n
        s = np.sign(st)
        dg = 0.0
        for _ in range(50):
            g = abs(st) - E * dg - sy(kap_n + dg)
            dgn = dg + g / (E + dsy(kap_n + dg))
            if abs(dgn - dg) < 1e-15:
                dg = dgn
                break
            dg = dgn
        sig = s * (abs(st) - E * dg)
        Dep = E * dsy(kap_n + dg) / (E + dsy(kap_n + dg))
        return sig, Dep, ep_n + s * dg, kap_n + dg

    def assemble(X, Xn, dt, upd=False):
        U, Q, P = X[:nu], X[nu:nu + nq], X[nu + nq:]
        Un = Xn[:nu]
        KT = np.zeros((nu, nu))
        Fint = np.zeros(nu)
        new_ep, new_kap = epn.copy(), kapn.copy()
        for e in range(n):
            for ig, (xi, w) in enumerate(zip(XG, WG)):
                pu, dpu = su.eval(e, xi, h)
                eps = dpu @ U[su.loc[e]]
                sig, Dep, ep1, kap1 = stress(eps, epn[e, ig], kapn[e, ig])
                new_ep[e, ig], new_kap[e, ig] = ep1, kap1
                if not consistent:
                    Dep = E                        # elastic (inconsistent) tangent
                loc = su.loc[e]
                Fint[loc] += w * h * dpu * sig
                KT[np.ix_(loc, loc)] += w * h * Dep * np.outer(dpu, dpu)
        if upd:
            return new_ep, new_kap
        R = np.zeros(N)
        R[:nu] = Fint - Qc @ P - ft
        R[nu:nu + nq] = A @ Q - B.T @ P
        R[nu + nq:] = Qc.T @ (U - Un) + S @ (P - Xn[nu + nq:]) + dt * (B @ Q)
        T = np.zeros((N, N))
        T[:nu, :nu] = KT
        T[:nu, nu + nq:] = -Qc
        T[nu:nu + nq, nu:nu + nq] = A
        T[nu:nu + nq, nu + nq:] = -B.T
        T[nu + nq:, :nu] = Qc.T
        T[nu + nq:, nu:nu + nq] = dt * B
        T[nu + nq:, nu + nq:] = S
        return R, T

    cons = np.array([0, nu])                       # u(0)=0, q(0)=0  (p_D = 0 natural)
    free = np.setdiff1d(np.arange(N), cons)
    X = np.zeros(N)
    X[nu + nq:] = p0_of(dict(E=E, k=k, alpha=al, invM=iM, sig0=sig0))
    # undrained step (dt=0) in the elastic range (sigma' = 0), then steps
    allres = []
    Xn = X.copy()
    for step in range(nsteps):
        X = Xn.copy()
        res = []
        for it in range(30):
            R, T = assemble(X, Xn, dt)
            r = np.linalg.norm(R[free])
            res.append(r)
            if r < 1e-13:
                break
            dX = np.zeros(N)
            dX[free] = np.linalg.solve(T[np.ix_(free, free)], -R[free])
            X = X + dX
        allres.append(res)
        epn, kapn = assemble(X, Xn, dt, upd=True)
        Xn = X.copy()
    if verbose:
        for s, res in enumerate(allres):
            print('  step %d: residual norms:' % (s + 1), ' '.join('%.2e' % r for r in res))
    return allres, epn, kapn


if __name__ == '__main__':
    np.set_printoptions(linewidth=160)
    print('###  (5a) Terzaghi: u-p vs u-q-p  (errors with respect to the analytical series)')
    main()
    d1, evmin, rel = schur_checks()
    print('=' * 100)
    print('###  Schur complement: reduced (U,P) system with S + dt B A^-1 B^T reproduces the full u-q-p step: max|diff| = %.3e' % d1)
    print('     H_up - H_mix (p(L)=0): smallest eigenvalue = %.3e (>= 0 up to round-off), relative norm of the difference = %.3e' % (evmin, rel))
    hs = hydrostatic()
    print('###  hydrostatic steady state (gamma_w = 2): max|p-p_hs| , max|q|  :', {k: ('%.2e' % v[0], '%.2e' % v[1]) for k, v in hs.items()})
    ip = inhomogeneous_pD()
    print('###  steady state with p_D = 0.5 at the top: max|p - 0.5| :', {k: '%.2e' % v for k, v in ip.items()})
    early_time()
    print('=' * 100)
    print('###  1D elastoplastic (nonlinear hardening) column, 3-field Newton')
    print(' consistent tangent  K_T = int B^T D_ep B :')
    plastic_newton(consistent=True)
    print(' elastic (inconsistent) tangent K_T = int B^T E B :')
    plastic_newton(consistent=False)
