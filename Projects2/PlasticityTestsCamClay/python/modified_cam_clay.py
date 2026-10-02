"""
modified_cam_clay.py
====================

Integração implícita (return mapping) do modelo Cam-Clay modificado.

Referência do algoritmo: de Souza Neto, Perić & Owen (2008), *Computational Methods
for Plasticity*, Seção 10.1, eqs. (10.1)-(10.22): sistema reduzido nas incógnitas
{Δγ, α_{n+1}} resolvido por Newton-Raphson com a jacobiana (10.20).

Extensões em relação ao texto do livro (mesmas do pacote ModifiedCamClay.m):
  * endurecimento de estado crítico  p_c(α) = p_c0 exp(v0 α/(λ-κ))  (RS2, eq. 8.8, v constante);
  * elasticidade "linear" (livro) ou "pressure_dependent": p = p_ini exp(-(v0/κ) ε_v^e),
    i.e. K = -v0 p/κ (RS2, eq. 8.7), com G constante ou ν constante (forma secante);
  * tensão inicial σ0;
  * operador tangente consistente 6x6.

Convenções (iguais às do código Mohr-Coulomb / HWTools.m):
  * tração positiva (p < 0 em compressão);
  * Voigt {xx, xy, xz, yy, yz, zz}; deformações com distorção de engenharia (γ = 2ε);
  * α = -ε_v^p (deformação volumétrica plástica, compressão positiva);
  * ε^e medida a partir do estado inicial (σ = σ0 quando ε^e = 0).
"""
import numpy as np

XX, XY, XZ, YY, YZ, ZZ = range(6)
M_VEC = np.array([1., 0., 0., 1., 0., 1.])                 # identidade em Voigt
SHEAR_FAC = np.array([1., .5, .5, 1., .5, 1.])             # engenharia -> componentes tensoriais
I_DEV = np.diag(SHEAR_FAC) - np.outer(M_VEC, M_VEC) / 3.0  # ε (Voigt, eng.) -> desviador (tensor)


def tnorm(t):
    """Norma de um tensor simétrico guardado em Voigt (componentes tensoriais)."""
    return np.sqrt(t[XX]**2 + t[YY]**2 + t[ZZ]**2 + 2.0 * (t[XY]**2 + t[XZ]**2 + t[YZ]**2))


def invariants(sig):
    """{p, q} de uma tensão em Voigt."""
    p = (sig[XX] + sig[YY] + sig[ZZ]) / 3.0
    return p, np.sqrt(1.5) * tnorm(sig - p * M_VEC)


# ----------------------------------------------------------------------------------------------
# Parâmetros
# ----------------------------------------------------------------------------------------------
def mcc_parameters(M=1.2, lam=0.077, kap=0.0066, N=1.788, v0=None, pc0=200., p0=200.,
                   pt=0., beta=1., elasticity='linear', shear='constant_nu', G=20000., nu=0.3,
                   sigma0=None):
    """Monta o dicionário de parâmetros.

    p0, pc0: valores positivos (compressão), como na Tabela 8.1 do RS2.
    v0=None -> v0 = N - λ ln(pc0) + κ ln(pc0/p0)  (linha virgem + linha de descarregamento).
    elasticity: 'linear' (K0 = v0 p0/κ constante) ou 'pressure_dependent' (K = -v0 p/κ).
    shear: 'constant_G' (usa G), 'constant_nu' (G = 3(1-2ν)/(2(1+ν)) K, forma secante s = 2G(p) e^e)
           ou 'hypo_nu' (hipoelástico: ds = 2 G_n de^e com G_n = 3(1-2ν)/(2(1+ν)) K(p_n), como no
           FLAC3D; exige sig_n e eps_n no return mapping).
    """
    if v0 is None:
        v0 = N - lam * np.log(pc0) + kap * np.log(pc0 / p0)
    if sigma0 is None:
        sigma0 = [-p0, 0., 0., -p0, 0., -p0]
    P = dict(M=M, lam=lam, kap=kap, N=N, v0=float(v0), pc0=pc0, p0=p0, pt=pt, beta=beta,
             elasticity=elasticity, shear=shear, G=G, nu=nu, sigma0=np.array(sigma0, float))
    P['pini'] = P['sigma0'][[XX, YY, ZZ]].sum() / 3.0            # p inicial (negativo)
    P['K0'] = -P['v0'] * P['pini'] / kap                          # K no estado inicial
    P['gfac'] = 3.0 * (1.0 - 2.0 * nu) / (2.0 * (1.0 + nu))
    P['G0'] = G if shear == 'constant_G' else P['gfac'] * P['K0']  # G no estado inicial
    P['e0'] = (P['sigma0'] - P['pini'] * M_VEC) / (2.0 * P['G0'])  # desviador inicial / 2G0
    return P


# ----------------------------------------------------------------------------------------------
# Lei elástica, endurecimento, função de escoamento
# ----------------------------------------------------------------------------------------------
def pressure(P, ee):
    """p(ε_v^e), K = dp/dε_v^e e dK/dε_v^e."""
    if P['elasticity'] == 'linear':
        return P['pini'] + P['K0'] * ee, P['K0'], 0.0
    c = P['v0'] / P['kap']
    p = P['pini'] * np.exp(-c * ee)
    K = -c * p
    return p, K, -c * K


def shear_modulus(P, ee):
    """G(ε_v^e) e dG/dε_v^e."""
    if P['shear'] == 'constant_G':
        return P['G'], 0.0
    if P['elasticity'] == 'linear':
        return P['gfac'] * P['K0'], 0.0
    _, K, dK = pressure(P, ee)
    return P['gfac'] * K, P['gfac'] * dK


def hardening(P, al):
    """a(α) (eq. 10.9) e H = da/dα (eq. 10.21) com p_c(α) = p_c0 exp(v0 α/(λ-κ))."""
    c = P['v0'] / (P['lam'] - P['kap'])
    pc = P['pc0'] * np.exp(c * al)
    return (pc + P['pt']) / (1.0 + P['beta']), c * pc / (1.0 + P['beta'])


def pc_of_alpha(P, al):
    return P['pc0'] * np.exp(P['v0'] * al / (P['lam'] - P['kap']))


def yield_function(P, p, q, a):
    """Eqs. (10.1)-(10.2) (p < 0 em compressão)."""
    b = 1.0 if p >= P['pt'] - a else P['beta']
    return (p - P['pt'] + a)**2 / b**2 + (q / P['M'])**2 - a**2


# ----------------------------------------------------------------------------------------------
# Return mapping
# ----------------------------------------------------------------------------------------------
class ReturnMappingError(RuntimeError):
    pass


def return_mapping(P, eps, epsp_n, alpha_n, tol=1e-11, maxit=50, verbose=False, sig_n=None, eps_n=None):
    """Return mapping implícito.

    eps: ε_{n+1} total; epsp_n: ε^p_n; alpha_n: α_n.  (Voigt, engenharia)
    sig_n, eps_n: tensão e deformação total do passo anterior (só para shear = 'hypo_nu').
    Retorna dict com 'stress', 'Dep' (tangente consistente), 'alpha', 'dgamma',
    'elastic_strain', 'plastic_strain', 'plastic', 'iterations', 'b', 'p', 'q', 'pc'.
    """
    M = P['M']; pt = P['pt']
    eps = np.asarray(eps, float); epsp_n = np.asarray(epsp_n, float)
    hypo = P['shear'] == 'hypo_nu'
    # ---- estado tentativa elástico (10.10)
    epse_tr = eps - epsp_n
    xv = epse_tr @ M_VEC                              # ε_v^{e,trial}
    if hypo:
        # G congelado no início do passo; s^trial = s_n + 2 G_n Δe  ->  ê = s_n/(2G_n) + Δe
        if sig_n is None:
            sig_n, eps_n = P['sigma0'], np.zeros(6)
        p_n = (sig_n[XX] + sig_n[YY] + sig_n[ZZ]) / 3.0
        Kn = P['K0'] if P['elasticity'] == 'linear' else -(P['v0'] / P['kap']) * p_n
        Gh = P['gfac'] * Kn
        ehat = (np.asarray(sig_n) - p_n * M_VEC) / (2.0 * Gh) + I_DEV @ (eps - np.asarray(eps_n))
    else:
        ehat = I_DEV @ epse_tr + P['e0']              # ε_d^{e,trial} (+ desviador inicial)
    nrm = tnorm(ehat)
    eq_tr = np.sqrt(2.0 / 3.0) * nrm
    n = ehat / nrm if nrm > 1e-14 else np.zeros(6)
    p_tr, K_tr, _ = pressure(P, xv)
    G_tr, dG_tr = (Gh, 0.0) if hypo else shear_modulus(P, xv)
    q_tr = 3.0 * G_tr * eq_tr
    a_n, _ = hardening(P, alpha_n)
    if yield_function(P, p_tr, q_tr, a_n) <= tol * a_n**2:
        # ---- passo elástico
        sig = 2.0 * G_tr * ehat + p_tr * M_VEC
        D = K_tr * np.outer(M_VEC, M_VEC) + 2.0 * G_tr * I_DEV + 2.0 * dG_tr * np.outer(ehat, M_VEC)
        if hypo:
            return dict(stress=sig, Dep=D, alpha=alpha_n, dgamma=0.0, elastic_strain=eps - epsp_n,
                        plastic_strain=epsp_n.copy(), plastic=False, iterations=0, b=1.0,
                        p=p_tr, q=q_tr, pc=pc_of_alpha(P, alpha_n))
        epse = (ehat - P['e0'] + xv / 3.0 * M_VEC) / SHEAR_FAC
        return dict(stress=sig, Dep=D, alpha=alpha_n, dgamma=0.0, elastic_strain=epse,
                    plastic_strain=eps - epse, plastic=False, iterations=0, b=1.0,
                    p=p_tr, q=q_tr, pc=pc_of_alpha(P, alpha_n))

    # ---- corretor plástico: Newton-Raphson no sistema reduzido (10.17)
    if P['beta'] == 1.0:
        blist = [1.0]
    else:
        blist = [1.0, P['beta']] if p_tr >= pt - a_n else [P['beta'], 1.0]
    done = False
    for b in blist:
        dg, al = 0.0, alpha_n
        conv = False
        for it in range(1, maxit + 1):
            ee = xv + al - alpha_n                      # ε_v^e_{n+1}
            p, K, _ = pressure(P, ee)                   # p(α)  (10.16)
            G, dG = (Gh, 0.0) if hypo else shear_modulus(P, ee)
            f = M**2 / (M**2 + 6.0 * G * dg)
            q = 3.0 * G * f * eq_tr                     # q(Δγ)  (10.15)
            a, H = hardening(P, al)
            pb = p - pt + a                             # p̄  (10.22)
            R1 = pb**2 / b**2 + (q / M)**2 - a**2
            R2 = al - alpha_n + dg * 2.0 * pb / b**2
            if verbose:
                print(f"   iter {it:2d}   |R1|/a^2 = {abs(R1)/a**2:.3e}   |R2| = {abs(R2):.3e}"
                      f"   dgamma = {dg:.6e}   alpha = {al:.6e}")
            if abs(R1) <= tol * a**2 and abs(R2) <= tol:
                conv = True
                break
            J = np.array([[-12.0 * G * f * q**2 / M**4,
                           2.0 * pb / b**2 * (K + H) + 2.0 * q / M**2 * (q * f / G) * dG - 2.0 * a * H],
                          [2.0 * pb / b**2, 1.0 + 2.0 * dg / b**2 * (K + H)]])   # (10.20)
            d = np.linalg.solve(J, -np.array([R1, R2]))
            dg += d[0]; al += d[1]
        if conv and (P['beta'] == 1.0 or (b == 1.0 and pb >= 0.0) or (b != 1.0 and pb < 0.0)):
            done = True
            break
    if not done or dg < -1e-12:
        raise ReturnMappingError(f'return mapping não convergiu (eps = {eps})')

    # ---- atualização (10.14), (10.18)
    sig = 2.0 * G * f * ehat + p * M_VEC
    if hypo:
        depsp = ((1.0 - f) * ehat - (al - alpha_n) / 3.0 * M_VEC) / SHEAR_FAC
        epse = eps - (epsp_n + depsp)
    else:
        epse = (f * ehat - P['e0'] + ee / 3.0 * M_VEC) / SHEAR_FAC
    # ---- tangente consistente
    dqdG = q * f / G
    dqdg = -6.0 * G * f * q / M**2
    J = np.array([[-12.0 * G * f * q**2 / M**4,
                   2.0 * pb / b**2 * (K + H) + 2.0 * q / M**2 * dqdG * dG - 2.0 * a * H],
                  [2.0 * pb / b**2, 1.0 + 2.0 * dg / b**2 * (K + H)]])
    dR1dx = 2.0 * pb / b**2 * K + 2.0 * q / M**2 * dqdG * dG     # ∂R1/∂ε_v^trial
    dR1deq = 2.0 * q / M**2 * 3.0 * G * f                        # ∂R1/∂ε_q^trial
    dR2dx = 2.0 * dg / b**2 * K                                  # ∂R2/∂ε_v^trial
    deq_de = np.sqrt(2.0 / 3.0) * n
    rhs = np.outer([dR1dx, dR2dx], M_VEC) + np.outer([dR1deq, 0.0], deq_de)
    ddg_de, dal_de = -np.linalg.solve(J, rhs)
    dq_de = dqdg * ddg_de + dqdG * dG * (M_VEC + dal_de) + 3.0 * G * f * deq_de
    dp_de = K * (M_VEC + dal_de)
    D = (2.0 * G * f * (I_DEV - np.outer(n, n)) + np.sqrt(2.0 / 3.0) * np.outer(n, dq_de)
         + np.outer(M_VEC, dp_de))
    return dict(stress=sig, Dep=D, alpha=al, dgamma=dg, elastic_strain=epse,
                plastic_strain=eps - epse, plastic=True, iterations=it, b=b,
                p=p, q=q, pc=pc_of_alpha(P, al))


def apply_strain_compute_sigma_dep(epst, epsp, alpha_n, P):
    """Mesma saída de ApplyStrainComputeSigmaDep do código MC: (σ, Dep, ε^e, α, tipo)."""
    r = return_mapping(P, epst, epsp, alpha_n)
    return r['stress'], r['Dep'], r['elastic_strain'], r['alpha'], int(r['plastic'])


def tangent_fd(P, eps, epsp, alpha_n, h=1e-8):
    """Tangente por diferenças finitas centradas (verificação)."""
    D = np.zeros((6, 6))
    for j in range(6):
        e1 = np.array(eps, float); e1[j] += h
        e2 = np.array(eps, float); e2[j] -= h
        D[:, j] = (return_mapping(P, e1, epsp, alpha_n)['stress']
                   - return_mapping(P, e2, epsp, alpha_n)['stress']) / (2 * h)
    return D


# ----------------------------------------------------------------------------------------------
# Ensaio triaxial drenado (ponto material, controle misto)
# ----------------------------------------------------------------------------------------------
def _triaxial_step(P, state, dea, tol, maxit):
    eps, epsp, al, ratio = state
    sc = P['sigma0'][XX]
    e = eps.copy()
    e[ZZ] -= dea
    e[XX] -= ratio * dea; e[YY] -= ratio * dea          # preditor das deformações laterais
    for it in range(maxit + 1):
        try:
            r = return_mapping(P, e, epsp, al)
        except ReturnMappingError:
            return None
        res = np.array([r['stress'][XX] - sc, r['stress'][YY] - sc])
        if np.abs(res).max() <= tol * abs(sc):
            new = (e, r['plastic_strain'], r['alpha'], (e[XX] - eps[XX]) / (-dea))
            return new, r, it
        Dl = r['Dep'][np.ix_([XX, YY], [XX, YY])]
        d = np.linalg.solve(Dl, -res)
        e[XX] += d[0]; e[YY] += d[1]
    return None


def _triaxial_advance(P, state, dea, level, tol, maxit, maxcut, counter):
    s = _triaxial_step(P, state, dea, tol, maxit)
    if s is not None:
        return s
    if level >= maxcut:
        return None
    counter[0] += 1
    s1 = _triaxial_advance(P, state, dea / 2, level + 1, tol, maxit, maxcut, counter)
    if s1 is None:
        return None
    s2 = _triaxial_advance(P, s1[0], dea / 2, level + 1, tol, maxit, maxcut, counter)
    if s2 is None:
        return None
    return s2[0], s2[1], s1[2] + s2[2]


def triaxial_drained(P, ea_max=0.2, nsteps=400, tol=1e-10, maxit=30, maxcut=8, verbose=False):
    """ε_zz prescrita (compressão), σ_xx = σ_yy = σ0_xx, distorções nulas.

    Retorna array (nsteps+1, 7) com colunas (compressão positiva):
    ε_a, p', q, ε_v, ε_q, σ_a, p_c
    """
    K0, G0 = P['K0'], P['G0']
    nu0 = (3 * K0 - 2 * G0) / (2 * (3 * K0 + G0))
    state = (np.zeros(6), np.zeros(6), 0.0, -nu0)
    p, q = invariants(P['sigma0'])
    out = [(0.0, -p, q, 0.0, 0.0, -P['sigma0'][ZZ], P['pc0'])]
    iters, counter = [], [0]
    dea = ea_max / nsteps
    for k in range(nsteps):
        s = _triaxial_advance(P, state, dea, 0, tol, maxit, maxcut, counter)
        if s is None:
            print(f'triaxial_drained: falha no passo {k + 1}')
            break
        state, r, it = s
        iters.append(it)
        eps, sig = state[0], r['stress']
        p, q = invariants(sig)
        out.append((-eps[ZZ], -p, q, -(eps[XX] + eps[YY] + eps[ZZ]),
                    2.0 / 3.0 * abs(eps[ZZ] - eps[XX]), -sig[ZZ], r['pc']))
    if verbose:
        print(f'passos: {len(iters)} | iterações globais por passo: máx = {max(iters)}, '
              f'média = {np.mean(iters):.2f} | subdivisões: {counter[0]}')
    return np.array(out)


# ----------------------------------------------------------------------------------------------
# Solução analítica (forma fechada) do ensaio triaxial drenado convencional
# ----------------------------------------------------------------------------------------------
def triaxial_analytical(P, npts=400, ea_max=np.inf):
    """Trajetória q = 3(p' - p'0) a partir de estado isotrópico, v = v0 constante.

    ε_v^p = (λ-κ)/v0 ln(p_c/p_c0),  p_c = p'(1 + η²/M²),  p' = 3p'0/(3-η)
    ε_q^p = (λ-κ)/v0 [F(η) - F(η_y)]  (integral de dε_q^p/dε_v^p = 2η/(M²-η²))
    F(x) = (1/M) ln|(M+x)/(M-x)| - (2/M) atan(x/M) - ln|M-x|/(3-M) - ln(M+x)/(3+M) + 6 ln(3-x)/(9-M²)
    Parcela elástica com a mesma lei do return mapping. Válida para p_t = 0 e β = 1.
    Retorna array com colunas: ε_a, p', q, ε_v, ε_q, σ_a  (compressão positiva).
    """
    if P['pt'] != 0.0 or P['beta'] != 1.0:
        print('triaxial_analytical: solução fechada válida apenas para pt = 0 e beta = 1.')
    M, lam, kap, v0 = P['M'], P['lam'], P['kap'], P['v0']
    p0, pc0 = -P['pini'], P['pc0']

    def F(x):
        return ((1 / M) * np.log(abs((M + x) / (M - x))) - (2 / M) * np.arctan(x / M)
                - np.log(abs(M - x)) / (3 - M) - np.log(M + x) / (3 + M) + 6 * np.log(3 - x) / (9 - M**2))

    A = 9 + M**2; B = -(18 * p0 + M**2 * pc0); C = 9 * p0**2
    disc = B * B - 4 * A * C
    py = min(r for r in ((-B + np.sqrt(disc)) / (2 * A), (-B - np.sqrt(disc)) / (2 * A)) if r >= p0 - 1e-9)
    qy = 3 * (py - p0); etay = qy / py

    def elastic(pp, q):
        if P['elasticity'] == 'linear':
            ev = (pp - p0) / P['K0']; K = P['K0']
        else:
            ev = (kap / v0) * np.log(pp / p0); K = v0 * pp / kap
        G = P['G'] if P['shear'] == 'constant_G' else P['gfac'] * K
        return ev, q / (3 * G)

    rows = []
    if py - p0 > 1e-9:
        for pp in np.linspace(p0, py, 41):
            q = 3 * (pp - p0)
            rows.append((pp, q) + elastic(pp, q))
    else:
        rows.append((p0, 0.0, 0.0, 0.0))
    g1 = etay + (M - etay) * np.arange(1, npts) / npts
    g2 = M - (M - etay) * np.exp(-12.0 * np.arange(1, npts + 1) / npts)
    etas = np.unique(np.r_[g1, g2])
    if etay > M:
        etas = etas[::-1]
    for eta in etas:
        pp = 3 * p0 / (3 - eta); q = eta * pp; pc = pp * (1 + eta**2 / M**2)
        ev, eq = elastic(pp, q)
        ev += (lam - kap) / v0 * np.log(pc / pc0)
        eq += (lam - kap) / v0 * (F(eta) - F(etay))
        rows.append((pp, q, ev, eq))
    R = np.array(rows)
    out = np.c_[R[:, 3] + R[:, 2] / 3, R[:, 0], R[:, 1], R[:, 2], R[:, 3], R[:, 0] + 2 * R[:, 1] / 3]
    return out[out[:, 0] <= ea_max]


def interp(data, x):
    """Interpolação linear em x (data: array Nx2)."""
    d = np.asarray(data)
    d = d[np.argsort(d[:, 0])]
    return float(np.interp(x, d[:, 0], d[:, 1]))
