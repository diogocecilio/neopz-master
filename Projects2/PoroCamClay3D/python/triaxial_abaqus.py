"""
triaxial_abaqus.py
==================

Abaqus Benchmarks Manual (v2016), 1.15.2 "Consolidation of a triaxial test specimen", reproduzido com o
FE u-p 3D (fe3d_up.py) e hexaedros de 20 nós (u quadrático, p linear nos vértices: Q2-Q1, o mesmo par do
C3D20P/C3D20RP do Abaqus).

Problema (Fig. 1.15.2-1): corpo de prova cilíndrico, altura/diâmetro = 3, metade superior modelada
(simetria no plano médio), H = 60 mm, r = 20 mm. Aqui em 3D: um quarto do cilindro (simetria em x = 0 e
y = 0), malha "O-grid" de Hex20 com os nós de meio de aresta da superfície lateral sobre o arco.
Cam-Clay modificado com elasticidade porosa: ν = 0.3, κ = 0.026, λ = 0.174, M = 1.0, a0 = 58.3 kPa
(p_c0 = 2 a0), e0 = 1.08 (v0 = 2.08); k = 1.728e-4 m/dia, γ_w = 10 kN/m³; tensão inicial efetiva
isotrópica de 100 kPa, pressão confinante P = 100 kPa constante na face lateral; a placa desce até
δ/H = 0.6 em 400 dias (34.56e6 s), com drenagem livre (p = 0) no topo. Placa lisa: só u_z prescrito no
topo; placa rugosa: também u_x = u_y = 0 no topo. Pequenas deformações (sem NLGEOM, como no Abaqus).

Unidades: m, kPa, kN, s.  Resultados no ponto A (r = 5 mm, z = 7.5 mm: centróide do elemento do Abaqus
mostrado na Fig. 1.15.2-1), interpolados dos pontos de Gauss.
"""
import json
import os
import time
import numpy as np
import fe3d_up as fe
import modified_cam_clay as mcc

# ------------------------------------------------------------------------------------------- dados
R, H = 0.020, 0.060                      # raio e meia altura (m)
P_CONF = 100.0                           # pressão confinante = tensão efetiva inicial (kPa)
E0 = 1.08
V0 = 1.0 + E0
KAP, LAM, M_CS, A0, NU = 0.026, 0.174, 1.0, 58.3, 0.3
K_HYD = 1.728e-4 / 86400.0               # condutividade hidráulica (m/s)
GAM_W = 10.0                             # peso específico da água (kN/m³)
MOBILITY = K_HYD / GAM_W                 # k/γ_w (m²/(kPa s))
T_END = 34.56e6                          # 400 dias (s)
DH_END = 0.6                             # δ/H no fim (corpo de prova a 40% da altura)
E_LIN = 15.0e3                           # E da variante com elasticidade linear (kPa): 15 MPa (o texto diz 15 GPa)
POINT_A = np.array([0.005, 0.0, 0.0075])  # r = 5 mm, z = 7.5 mm (no plano y = 0)

# hexaedro de 20 nós com integração reduzida 2x2x2 (análogo do C3D20RP / CAX8RP do Abaqus)
fe.ELEMENTS.setdefault('hex20r', dict(fe.ELEMENTS['hex20'], rule=lambda: fe.gauss_hex(2)))

# marcadores das faces de contorno do quarto de cilindro
X0, LATERAL, Y0, BASE, TOPO = 1, 2, 3, 5, 6


def parameters(elastic='porous'):
    """Cam-Clay modificado do benchmark. 'porous': p = p0 exp(-v0 ε_v^e/κ), G de ν na forma incremental
    (ds = 2G(p_n) de, a elasticidade porosa do Abaqus); 'linear': E = 15 MPa, ν = 0.3."""
    s0 = [-P_CONF, 0.0, 0.0, -P_CONF, 0.0, -P_CONF]
    if elastic == 'porous':
        return mcc.mcc_parameters(M=M_CS, lam=LAM, kap=KAP, v0=V0, pc0=2 * A0, p0=P_CONF,
                                  elasticity='pressure_dependent', shear='hypo_nu', nu=NU, sigma0=s0)
    G = E_LIN / (2 * (1 + NU))
    P = mcc.mcc_parameters(M=M_CS, lam=LAM, kap=KAP, v0=V0, pc0=2 * A0, p0=P_CONF,
                           elasticity='linear', shear='constant_G', G=G, nu=NU, sigma0=s0)
    P['K0'] = E_LIN / (3 * (1 - 2 * NU))
    P['G0'] = G
    return P


# ------------------------------------------------------------------------------------------- malha
def quarter_disk(r, nc, nr, a_ratio=0.5):
    """Malha 'O-grid' de quadriláteros de um quarto de círculo: quadrado central [0,a]² (nc x nc) e dois
    blocos externos (nc ao longo do arco de 45°, nr na direção radial). Retorna (coords2d, quads CCW)."""
    a = a_ratio * r
    index, pts = {}, []

    def node(p):
        key = (round(p[0] / r, 9), round(p[1] / r, 9))
        if key not in index:
            index[key] = len(pts)
            pts.append((float(p[0]), float(p[1])))
        return index[key]

    quads = []
    S = [[node((a * i / nc, a * j / nc)) for j in range(nc + 1)] for i in range(nc + 1)]
    for i in range(nc):
        for j in range(nc):
            quads.append([S[i][j], S[i + 1][j], S[i + 1][j + 1], S[i][j + 1]])
    for block in ('leste', 'norte'):
        G = []
        for j in range(nc + 1):
            if block == 'leste':
                inner, th = np.array([a, a * j / nc]), 0.25 * np.pi * j / nc
            else:
                inner, th = np.array([a * j / nc, a]), 0.5 * np.pi - 0.25 * np.pi * j / nc
            outer = r * np.array([np.cos(th), np.sin(th)])
            G.append([node(inner + (k / nr) * (outer - inner)) for k in range(nr + 1)])
        for j in range(nc):
            for k in range(nr):
                quads.append([G[j][k], G[j][k + 1], G[j + 1][k + 1], G[j + 1][k]])
    C = np.array(pts)
    for q in quads:                                    # orientação anti-horária (normal +z)
        x, y = C[q, 0], C[q, 1]
        if 0.5 * np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y) < 0:
            q[1], q[3] = q[3], q[1]
    return C, quads


def cylinder_mesh(r=R, h=H, nc=2, nr=2, nz=4, etype='hex20', a_ratio=0.5):
    """Quarto de cilindro 0 <= z <= h com hexaedros (ordem dos nós do NDSolve`FEM`: vértices e depois os nós
    de meio de aresta). Os nós de meio de aresta da superfície lateral são levados ao arco r = R.
    Marcadores: 1 x = 0, 2 lateral, 3 y = 0, 5 base (z = 0), 6 topo (z = h)."""
    C2, quads = quarter_disk(r, nc, nr, a_ratio)
    n2 = len(C2)
    coords = [[x, y, h * k / nz] for k in range(nz + 1) for x, y in C2]
    cells = [[k * n2 + q[0], k * n2 + q[1], k * n2 + q[2], k * n2 + q[3],
              (k + 1) * n2 + q[0], (k + 1) * n2 + q[1], (k + 1) * n2 + q[2], (k + 1) * n2 + q[3]]
             for k in range(nz) for q in quads]
    nvert = len(coords)
    edge_node = {}
    tol = 1e-9 * r

    def on_lateral(i):
        return abs(np.hypot(coords[i][0], coords[i][1]) - r) < tol

    def mid(a, b):
        key = (min(a, b), max(a, b))
        if key not in edge_node:
            p = 0.5 * (np.array(coords[a]) + np.array(coords[b]))
            if on_lateral(a) and on_lateral(b):           # aresta na superfície lateral: nó no arco
                p[:2] *= r / np.hypot(p[0], p[1])
            edge_node[key] = len(coords)
            coords.append(list(p))
        return edge_node[key]

    if etype == 'hex20':
        elements = [c + [mid(c[a], c[b]) for a, b in fe.HEX20_EDGES] for c in cells]
    else:
        elements = [list(c) for c in cells]
    coords = np.array(coords)
    lf = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
    count = {}
    for e in elements:
        for f in lf:
            count.setdefault(tuple(sorted(e[i] for i in f)), []).append([e[i] for i in f])
    faces = []
    for lst in count.values():
        if len(lst) != 1:
            continue
        f = lst[0]
        X = coords[f]
        if np.all(abs(X[:, 0]) < tol): m = X0
        elif np.all(abs(X[:, 1]) < tol): m = Y0
        elif np.all(abs(X[:, 2]) < tol): m = BASE
        elif np.all(abs(X[:, 2] - h) < tol): m = TOPO
        else: m = LATERAL
        if etype == 'hex20':
            f = f + [edge_node[(min(f[a], f[b]), max(f[a], f[b]))] for a, b in fe.QUAD8_EDGES]
        faces.append((f, m))
    pnodes = sorted({v for e in elements for v in e[:8]})
    return dict(coords=coords, elements=[np.array(e) for e in elements], faces=faces, etype=etype,
                pnodes=np.array(pnodes), nvert=nvert, L=(r, r, h), n=(nc, nr, nz))


def face_pressure(mesh, select, pressure, outward):
    """Forças nodais consistentes de uma pressão (positiva comprimindo) nas faces com select(f, m):
    F_a = -∫ N_a p n dA, com a normal n orientada por outward(x) (n · outward > 0)."""
    el = fe.ELEMENTS[mesh['etype']]
    C = mesh['coords']
    F = np.zeros(3 * len(C))
    pts, wts = el['frule']()
    for f, m in mesh['faces']:
        if not select(f, m):
            continue
        X = C[f]
        for xi, w in zip(pts, wts):
            N, dN = el['shf'](xi)
            J = dN @ X
            n = np.cross(J[0], J[1])                      # |n| = dA/(dξ dη)
            if n @ outward(N @ X) < 0:
                n = -n
            for a, nd in enumerate(f):
                F[3 * nd:3 * nd + 3] -= pressure * N[a] * n * w
    return F


# ------------------------------------------------------------------------------------------- pós-processamento
def locate(mesh, x, tol=1e-9):
    """Elemento que contém o ponto x e as coordenadas paramétricas (Newton no mapeamento isoparamétrico)."""
    el = fe.ELEMENTS[mesh['etype']]
    C = mesh['coords']
    for k, e in enumerate(mesh['elements']):
        X = C[e]
        if np.any(x < X.min(axis=0) - 1e-12) or np.any(x > X.max(axis=0) + 1e-12):
            continue
        xi = np.zeros(3)
        for _ in range(30):
            N, dN = el['shu'](xi)
            r = N @ X - x
            if np.linalg.norm(r) < 1e-14:
                break
            xi -= np.linalg.solve((dN @ X).T, r)
        if np.all(np.abs(xi) <= 1 + 1e-7):
            return k, xi
    raise ValueError('ponto fora da malha')


def gauss_interpolator(mesh, xi):
    """Pesos de interpolação dos valores nos pontos de Gauss (n x n x n) para o ponto paramétrico xi
    (polinômios de Lagrange nas abscissas de Gauss: exato para campos de grau n-1 por direção)."""
    pts, _ = fe.ELEMENTS[mesh['etype']]['rule']()
    g = np.unique(np.round(pts[:, 0], 14))

    def lag(t):
        return np.array([np.prod([(t - g[m]) / (g[l] - g[m]) for m in range(len(g)) if m != l])
                         for l in range(len(g))])

    Lx, Ly, Lz = lag(xi[0]), lag(xi[1]), lag(xi[2])
    w = []
    for p in pts:
        i, j, k = (int(np.argmin(abs(g - p[c]))) for c in range(3))
        w.append(Lx[i] * Ly[j] * Lz[k])
    return np.array(w)


# ------------------------------------------------------------------------------------------- modelo e passos
def build(platen='rough', elastic='porous', nc=2, nr=2, nz=4, etype='hex20', a_ratio=0.5):
    mesh = cylinder_mesh(R, H, nc, nr, nz, 'hex20' if etype.startswith('hex20') else etype, a_ratio)
    mesh['etype'] = etype
    par = parameters(elastic)
    model = fe.PoroField(mesh, lambda x: par, alpha=1.0, biot_modulus=np.inf, perm=MOBILITY, gam_w=0.0,
                         body=(0.0, 0.0, 0.0), p0_of_x=lambda x: 0.0, init='gp')
    C = mesh['coords']
    tol = 1e-9 * R
    sym_x = np.where(abs(C[:, 0]) < tol)[0]
    sym_y = np.where(abs(C[:, 1]) < tol)[0]
    base = np.where(abs(C[:, 2]) < tol)[0]
    top = np.where(abs(C[:, 2] - H) < tol)[0]
    fixed0 = {}
    for nd in sym_x: fixed0[3 * nd] = 0.0
    for nd in sym_y: fixed0[3 * nd + 1] = 0.0
    for nd in base: fixed0[3 * nd + 2] = 0.0
    if platen == 'rough':
        for nd in top:
            fixed0[3 * nd] = 0.0
            fixed0[3 * nd + 1] = 0.0
    fixed_p = {int(model.pmap[nd]): 0.0 for nd in top if model.pmap[nd] >= 0}
    Fp = face_pressure(mesh, lambda f, m: m == LATERAL, P_CONF, lambda x: np.array([x[0], x[1], 0.0]))
    kA, xiA = locate(mesh, POINT_A)
    wA = gauss_interpolator(mesh, xiA)
    return dict(mesh=mesh, model=model, par=par, fixed0=fixed0, top=top, fixed_p=fixed_p, F=Fp,
                kA=kA, wA=wA, platen=platen, elastic=elastic)


def fixed_at(case, dh):
    fx = dict(case['fixed0'])
    for nd in case['top']:
        fx[3 * nd + 2] = -dh * H
    return fx


def stress_at_A(case):
    st = case['model'].state[case['kA']]
    sig = np.array([s['sig'] for s in st])                 # tensões efetivas nos pontos de Gauss
    return case['wA'] @ sig


def platen_stress(case):
    """tensão axial média sob a placa: -Σ R_z(topo)/(π R²/4) (compressão positiva)."""
    m = case['model']
    Rv = m.Fint - m.Q @ m.P - case['F']
    return -Rv[3 * case['top'] + 2].sum() / (0.25 * np.pi * R * R)


def advance(case, dh0, dh1, rate, level=0, maxlevel=8):
    """Avança a placa de δ/H = dh0 a dh1 (Δt proporcional). Preditor: U_n + rate (dh1 - dh0), com rate o
    último incremento convergido por unidade de δ/H; se o Newton ou o return mapping falharem, divide o
    incremento ao meio. Retorna (iterações, rate, subdivisões)."""
    m = case['model']
    dt = T_END * (dh1 - dh0) / DH_END
    try:
        U0 = m.U.copy()
        it = m.step(case['F'], fixed_at(case, dh1), case['fixed_p'], dt=dt, flow=True,
                    U_guess=m.U + rate * (dh1 - dh0))
        return it, (m.U - U0) / (dh1 - dh0), 0
    except (mcc.ReturnMappingError, RuntimeError, np.linalg.LinAlgError):
        if level >= maxlevel:
            raise
        dm = 0.5 * (dh0 + dh1)
        it1, rate, c1 = advance(case, dh0, dm, rate, level + 1, maxlevel)
        it2, rate, c2 = advance(case, dm, dh1, rate, level + 1, maxlevel)
        return it1 + it2, rate, c1 + c2 + 1


def run(platen='rough', elastic='porous', nc=2, nr=2, nz=4, nsteps=150, etype='hex20', verbose=True):
    case = build(platen, elastic, nc, nr, nz, etype)
    m = case['model']
    C = case['mesh']['coords']
    t0 = time.time()
    # etapa geostática: verifica o equilíbrio do estado inicial com a pressão confinante
    F0, _, _ = m.assemble(m.U)
    m.Fint = F0
    free = np.setdiff1d(np.arange(m.nu), list(case['fixed0']) + [3 * nd + 2 for nd in case['top']])
    res0 = np.abs((F0 - case['F'])[free]).max()
    it = m.step(case['F'], fixed_at(case, 0.0), case['fixed_p'], dt=1.0, flow=False)
    hist = dict(dh=[0.0], p=[], q=[], sa=[], pmax=[0.0], it=[it])
    pA, qA = mcc.invariants(stress_at_A(case))
    hist['p'].append(-pA); hist['q'].append(qA); hist['sa'].append(platen_stress(case))
    # preditor do 1º passo: campo homogêneo elástico (u_z = -z δ/H, expansão lateral ν)
    rate = np.zeros(m.nu)
    rate[2::3] = -C[:, 2]
    rate[0::3] = NU * C[:, 0]; rate[1::3] = NU * C[:, 1]
    cuts = 0
    for k in range(1, nsteps + 1):
        dh0, dh = DH_END * (k - 1) / nsteps, DH_END * k / nsteps
        it, rate, c = advance(case, dh0, dh, rate)
        cuts += c
        pA, qA = mcc.invariants(stress_at_A(case))
        hist['dh'].append(dh); hist['p'].append(-pA); hist['q'].append(qA)
        hist['sa'].append(platen_stress(case)); hist['pmax'].append(np.abs(m.P).max()); hist['it'].append(it)
    hist = {k: np.array(v) for k, v in hist.items()}
    if verbose:
        print(f"{etype} {platen:6s} {elastic:6s} malha nc={nc} nr={nr} nz={nz}: {len(C)} nós, "
              f"{len(case['mesh']['elements'])} elementos; resíduo inicial {res0:.1e} kN; "
              f"Newton {hist['it'][1:].mean():.1f} it/passo (máx {hist['it'].max()}), subdivisões {cuts}; "
              f"|p_poro| máx {hist['pmax'].max():.4f} kPa; q(A) final {hist['q'][-1]:.2f} kPa; "
              f"{time.time() - t0:.0f} s")
    case['hist'] = hist
    case['res0'] = res0
    case['cuts'] = cuts
    return case


# ------------------------------------------------------------------------------------------- referência homogênea
def homogeneous(elastic='porous', nsteps=600, dh_end=DH_END):
    """Placa lisa = estado homogêneo: ε_zz = -δ/H prescrita, σ_xx = σ_yy = -P, distorções nulas
    (Newton nas deformações laterais, com σ_n e ε_n do passo anterior para a lei incremental)."""
    P = parameters(elastic)
    XX, YY, ZZ = mcc.XX, mcc.YY, mcc.ZZ
    eps, epsp, al = np.zeros(6), np.zeros(6), 0.0
    sig = P['sigma0'].copy()
    out = [(0.0, P_CONF, 0.0, 0.0)]
    ratio = NU / (1 - NU)
    for k in range(1, nsteps + 1):
        dea = dh_end / nsteps
        e = eps.copy()
        e[ZZ] -= dea
        e[XX] += ratio * dea; e[YY] += ratio * dea
        for _ in range(50):
            r = mcc.return_mapping(P, e, epsp, al, sig_n=sig, eps_n=eps)
            res = np.array([r['stress'][XX] + P_CONF, r['stress'][YY] + P_CONF])
            if np.abs(res).max() < 1e-10 * P_CONF:
                break
            d = np.linalg.solve(r['Dep'][np.ix_([XX, YY], [XX, YY])], -res)
            e[XX] += d[0]; e[YY] += d[1]
        ratio = (e[XX] - eps[XX]) / dea
        eps, epsp, al, sig = e, r['plastic_strain'], r['alpha'], r['stress']
        p, q = mcc.invariants(sig)
        out.append((dh_end * k / nsteps, -p, q, -(e[XX] + e[YY] + e[ZZ])))
    return np.array(out)


def load_abaqus(path=None):
    path = path or os.path.join(os.path.dirname(os.path.abspath(__file__)), 'abaqus_1_15_2_digitalizado.json')
    return json.load(open(path))


if __name__ == '__main__':
    import sys
    t = time.time()
    hom = homogeneous()
    print('referência homogênea: q em δ/H = 0.03, 0.1, 0.2, 0.3, 0.5, 0.6:',
          np.round(np.interp([0.03, 0.1, 0.2, 0.3, 0.5, 0.6], hom[:, 0], hom[:, 2]), 2), f'({time.time() - t:.1f} s)')
    for platen in (sys.argv[1:] or ['smooth', 'rough']):
        run(platen)
