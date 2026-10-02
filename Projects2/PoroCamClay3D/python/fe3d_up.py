"""
fe3d_up.py
==========

Elementos finitos acoplados u-p (Biot) 3D com Cam-Clay modificado, generalizado para quatro tipos de elemento:

    'hex8'   hexaedro de 8 nós   (u e p trilineares, Q1-Q1: mesma ordem, viola inf-sup no limite não drenado)
    'hex20'  hexaedro de 20 nós  (u serendipity quadrático, p trilinear nos 8 vértices: Q2-Q1, estável)
    'tet4'   tetraedro de 4 nós  (P1-P1, mesma ordem)
    'tet10'  tetraedro de 10 nós (u quadrático, p linear nos 4 vértices: P2-P1, Taylor-Hood, estável)

Ordem dos nós igual à do NDSolve`FEM` do Mathematica (vértices primeiro, depois os nós de meio de aresta),
de modo que este código serve de referência numérica para o pacote Mathematica.

Formulação (formulacao.tex):  K u - G p = f_u ;  Q u' + S p' + H p = f_p  (Q = Gᵀ), Euler implícito
multiplicado por Δt:   [[K_T, -G], [Gᵀ, S + Δt H]] {Δu, Δp} = -{R_u, R_p} (Newton), com
R_u = F_int(σ') - G P - F_ext,   R_p = Gᵀ(U - U_n) + S(P - P_n) + Δt (H P - f_g),
G = α ∫ B_uᵀ m N_p dΩ,  S = ∫ (1/M) N_pᵀ N_p dΩ,  H = ∫ k ∇N_pᵀ ∇N_p dΩ,  f_g = ∫ k ∇N_pᵀ (ρ_w g⃗) dΩ.
Convenções: tração positiva, p > 0 em compressão, σ = σ' - α p m; Voigt do FE {xx, yy, zz, γyz, γxz, γxy}.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import modified_cam_clay as mcc

# ======================================================================================= elementos de referência
HEX8_REF = np.array([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
                     [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]], float)
HEX20_EDGES = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4), (0, 4), (1, 5), (2, 6), (3, 7)]
HEX20_REF = np.vstack([HEX8_REF, [(HEX8_REF[a] + HEX8_REF[b]) / 2 for a, b in HEX20_EDGES]])
TET4_REF = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], float)
TET10_EDGES = [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)]
QUAD4_REF = np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]], float)
QUAD8_EDGES = [(0, 1), (1, 2), (2, 3), (3, 0)]
QUAD8_REF = np.vstack([QUAD4_REF, [(QUAD4_REF[a] + QUAD4_REF[b]) / 2 for a, b in QUAD8_EDGES]])
TRI6_EDGES = [(0, 1), (1, 2), (2, 0)]


def shape_hex8(x):
    s = HEX8_REF
    a, b, c = 1 + s[:, 0] * x[0], 1 + s[:, 1] * x[1], 1 + s[:, 2] * x[2]
    N = a * b * c / 8
    dN = np.array([s[:, 0] * b * c, s[:, 1] * a * c, s[:, 2] * a * b]) / 8
    return N, dN


def shape_hex20(x):
    N = np.zeros(20); dN = np.zeros((3, 20))
    for i, (xi, yi, zi) in enumerate(HEX20_REF):
        a, b, c = 1 + xi * x[0], 1 + yi * x[1], 1 + zi * x[2]
        if i < 8:
            g = xi * x[0] + yi * x[1] + zi * x[2] - 2
            N[i] = a * b * c * g / 8
            dN[0, i] = xi * b * c * (g + a) / 8
            dN[1, i] = yi * a * c * (g + b) / 8
            dN[2, i] = zi * a * b * (g + c) / 8
        elif xi == 0:
            N[i] = (1 - x[0] ** 2) * b * c / 4
            dN[:, i] = [-2 * x[0] * b * c / 4, (1 - x[0] ** 2) * yi * c / 4, (1 - x[0] ** 2) * b * zi / 4]
        elif yi == 0:
            N[i] = a * (1 - x[1] ** 2) * c / 4
            dN[:, i] = [xi * (1 - x[1] ** 2) * c / 4, -2 * x[1] * a * c / 4, a * (1 - x[1] ** 2) * zi / 4]
        else:
            N[i] = a * b * (1 - x[2] ** 2) / 4
            dN[:, i] = [xi * b * (1 - x[2] ** 2) / 4, a * yi * (1 - x[2] ** 2) / 4, -2 * x[2] * a * b / 4]
    return N, dN


def _bary(x):
    L = np.array([1 - x[0] - x[1] - x[2], x[0], x[1], x[2]])
    dL = np.array([[-1, 1, 0, 0], [-1, 0, 1, 0], [-1, 0, 0, 1]], float)
    return L, dL


def shape_tet4(x):
    return _bary(x)


def shape_tet10(x):
    L, dL = _bary(x)
    N = np.zeros(10); dN = np.zeros((3, 10))
    for i in range(4):
        N[i] = L[i] * (2 * L[i] - 1); dN[:, i] = (4 * L[i] - 1) * dL[:, i]
    for k, (i, j) in enumerate(TET10_EDGES):
        N[4 + k] = 4 * L[i] * L[j]; dN[:, 4 + k] = 4 * (L[i] * dL[:, j] + L[j] * dL[:, i])
    return N, dN


def shape_quad4(x):
    s = QUAD4_REF
    a, b = 1 + s[:, 0] * x[0], 1 + s[:, 1] * x[1]
    return a * b / 4, np.array([s[:, 0] * b, s[:, 1] * a]) / 4


def shape_quad8(x):
    N = np.zeros(8); dN = np.zeros((2, 8))
    for i, (xi, yi) in enumerate(QUAD8_REF):
        a, b = 1 + xi * x[0], 1 + yi * x[1]
        if i < 4:
            g = xi * x[0] + yi * x[1] - 1
            N[i] = a * b * g / 4; dN[:, i] = [xi * b * (g + a) / 4, yi * a * (g + b) / 4]
        elif xi == 0:
            N[i] = (1 - x[0] ** 2) * b / 2; dN[:, i] = [-x[0] * b, (1 - x[0] ** 2) * yi / 2]
        else:
            N[i] = a * (1 - x[1] ** 2) / 2; dN[:, i] = [xi * (1 - x[1] ** 2) / 2, -x[1] * a]
    return N, dN


def shape_tri3(x):
    return np.array([1 - x[0] - x[1], x[0], x[1]]), np.array([[-1, 1, 0], [-1, 0, 1]], float)


def shape_tri6(x):
    L = np.array([1 - x[0] - x[1], x[0], x[1]]); dL = np.array([[-1, 1, 0], [-1, 0, 1]], float)
    N = np.zeros(6); dN = np.zeros((2, 6))
    for i in range(3):
        N[i] = L[i] * (2 * L[i] - 1); dN[:, i] = (4 * L[i] - 1) * dL[:, i]
    for k, (i, j) in enumerate(TRI6_EDGES):
        N[3 + k] = 4 * L[i] * L[j]; dN[:, 3 + k] = 4 * (L[i] * dL[:, j] + L[j] * dL[:, i])
    return N, dN


def gauss_line(n):
    if n == 1: return np.array([0.0]), np.array([2.0])
    if n == 2: return np.array([-1, 1]) / np.sqrt(3.0), np.array([1.0, 1.0])
    if n == 3: return np.array([-np.sqrt(0.6), 0.0, np.sqrt(0.6)]), np.array([5, 8, 5]) / 9.0
    raise ValueError(n)


def gauss_hex(n):
    p, w = gauss_line(n)
    pts = [(p[i], p[j], p[k]) for i in range(n) for j in range(n) for k in range(n)]
    wts = [w[i] * w[j] * w[k] for i in range(n) for j in range(n) for k in range(n)]
    return np.array(pts), np.array(wts)


def gauss_quad(n):
    p, w = gauss_line(n)
    return np.array([(p[i], p[j]) for i in range(n) for j in range(n)]), \
        np.array([w[i] * w[j] for i in range(n) for j in range(n)])


def rule_tet4pt():
    a, b = 0.5854101966249685, 0.1381966011250105
    pts = np.array([[b, b, b], [a, b, b], [b, a, b], [b, b, a]])
    return pts, np.full(4, 1.0 / 24.0)


def rule_tri3pt():
    return np.array([[1 / 6, 1 / 6], [2 / 3, 1 / 6], [1 / 6, 2 / 3]]), np.full(3, 1.0 / 6.0)


# tipo: (função de forma de u, nº de nós de u, função de forma de p, nº de nós de p, regra, face, regra da face)
ELEMENTS = {
    'hex8': dict(shu=shape_hex8, nu=8, shp=shape_hex8, np_=8, rule=lambda: gauss_hex(2),
                 shf=shape_quad4, nf=4, frule=lambda: gauss_quad(2)),
    'hex20': dict(shu=shape_hex20, nu=20, shp=shape_hex8, np_=8, rule=lambda: gauss_hex(3),
                  shf=shape_quad8, nf=8, frule=lambda: gauss_quad(3)),
    'tet4': dict(shu=shape_tet4, nu=4, shp=shape_tet4, np_=4, rule=rule_tet4pt,
                 shf=shape_tri3, nf=3, frule=rule_tri3pt),
    'tet10': dict(shu=shape_tet10, nu=10, shp=shape_tet4, np_=4, rule=rule_tet4pt,
                  shf=shape_tri6, nf=6, frule=rule_tri3pt),
}

# Voigt do FE {xx, yy, zz, γyz, γxz, γxy}  <->  Voigt do Cam-Clay {xx, xy, xz, yy, yz, zz}
FE2MCC = [0, 5, 4, 1, 3, 2]
MCC2FE = [0, 3, 5, 4, 2, 1]
M_FE = np.array([1., 1., 1., 0., 0., 0.])


def compute_B(gradphi):
    """ComputeBN (dim = 3) do poro-elastic.nb: linhas {xx, yy, zz, γyz, γxz, γxy}."""
    n = gradphi.shape[1]
    B = np.zeros((6, 3 * n))
    for i in range(n):
        dx, dy, dz = gradphi[:, i]
        B[0, 3 * i] = dx; B[1, 3 * i + 1] = dy; B[2, 3 * i + 2] = dz
        B[3, 3 * i + 1] = dz; B[3, 3 * i + 2] = dy
        B[4, 3 * i] = dz; B[4, 3 * i + 2] = dx
        B[5, 3 * i] = dy; B[5, 3 * i + 1] = dx
    return B


# ======================================================================================= malha estruturada
def box_mesh(L, n, etype='hex8'):
    """Caixa [0,Lx]×[0,Ly]×[0,Lz] com nx×ny×nz células. Marcadores de face: 1 x=0, 2 x=Lx, 3 y=0, 4 y=Ly,
    5 z=0, 6 z=Lz.  Retorna dict com coords, elements (ordem do NDSolve`FEM`), faces [(nós, marcador)],
    pnodes (nós de pressão = vértices), cells (8 vértices de cada célula, para valores de 'zona')."""
    Lx, Ly, Lz = L; nx, ny, nz = n
    xs, ys, zs = (np.linspace(0, Lx, nx + 1), np.linspace(0, Ly, ny + 1), np.linspace(0, Lz, nz + 1))
    idx = lambda i, j, k: k * (ny + 1) * (nx + 1) + j * (nx + 1) + i
    coords = [[xs[i], ys[j], zs[k]] for k in range(nz + 1) for j in range(ny + 1) for i in range(nx + 1)]
    cells = []
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                cells.append([idx(i, j, k), idx(i + 1, j, k), idx(i + 1, j + 1, k), idx(i, j + 1, k),
                              idx(i, j, k + 1), idx(i + 1, j, k + 1), idx(i + 1, j + 1, k + 1), idx(i, j + 1, k + 1)])
    nvert = len(coords)
    edge_node = {}

    def mid(a, b):
        key = (min(a, b), max(a, b))
        if key not in edge_node:
            edge_node[key] = len(coords)
            coords.append(list((np.array(coords[a]) + np.array(coords[b])) / 2))
        return edge_node[key]

    if etype in ('hex8', 'hex20'):
        elements = [list(c) for c in cells]
        if etype == 'hex20':
            elements = [c + [mid(c[a], c[b]) for a, b in HEX20_EDGES] for c in cells]
    else:
        # divisão de Kuhn: 6 tetraedros em torno da diagonal 0-6 de cada célula (malha conforme)
        paths = [(1, 2), (1, 5), (4, 5), (4, 7), (3, 7), (3, 2)]
        C = np.array(coords)
        elements = []
        for c in cells:
            for a, b in paths:
                t = [c[0], c[a], c[b], c[6]]
                if np.linalg.det(np.array([C[t[1]] - C[t[0]], C[t[2]] - C[t[0]], C[t[3]] - C[t[0]]])) < 0:
                    t = [t[0], t[2], t[1], t[3]]
                elements.append(t)
        if etype == 'tet10':
            elements = [t + [mid(t[a], t[b]) for a, b in TET10_EDGES] for t in elements]
    coords = np.array(coords)
    # faces de contorno: faces de vértices que aparecem uma única vez
    if etype in ('hex8', 'hex20'):
        lf = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
    else:
        lf = [(0, 2, 1), (0, 1, 3), (1, 2, 3), (0, 3, 2)]
    count = {}
    for e in elements:
        for f in lf:
            key = tuple(sorted(e[i] for i in f))
            count.setdefault(key, []).append([e[i] for i in f])
    faces = []
    tol = 1e-9 * max(L)
    for key, lst in count.items():
        if len(lst) != 1:
            continue
        f = lst[0]; X = coords[f]
        if np.all(abs(X[:, 0]) < tol): m = 1
        elif np.all(abs(X[:, 0] - Lx) < tol): m = 2
        elif np.all(abs(X[:, 1]) < tol): m = 3
        elif np.all(abs(X[:, 1] - Ly) < tol): m = 4
        elif np.all(abs(X[:, 2]) < tol): m = 5
        else: m = 6
        if etype == 'hex20':
            f = f + [edge_node[(min(f[a], f[b]), max(f[a], f[b]))] for a, b in QUAD8_EDGES]
        elif etype == 'tet10':
            f = f + [edge_node[(min(f[a], f[b]), max(f[a], f[b]))] for a, b in TRI6_EDGES]
        faces.append((f, m))
    pnodes = sorted({v for e in elements for v in e[:ELEMENTS[etype]['np_']]})
    return dict(coords=coords, elements=[np.array(e) for e in elements], faces=faces, etype=etype,
                pnodes=np.array(pnodes), nvert=nvert, cells=[np.array(c) for c in cells], L=L, n=n)


def nodes_on(mesh, marker):
    return sorted({v for f, m in mesh['faces'] if m == marker for v in f})


def face_load(mesh, select, t):
    """forças nodais consistentes de uma tração uniforme t nas faces f com select(f, marcador) verdadeiro."""
    el = ELEMENTS[mesh['etype']]
    F = np.zeros(3 * len(mesh['coords']))
    pts, wts = el['frule']()
    for f, m in mesh['faces']:
        if not select(f, m):
            continue
        X = mesh['coords'][f]
        for xi, w in zip(pts, wts):
            N, dN = el['shf'](xi)
            J = dN @ X
            da = np.linalg.norm(np.cross(J[0], J[1]))
            for a, nd in enumerate(f):
                F[3 * nd:3 * nd + 3] += N[a] * np.asarray(t) * w * da
    return F


def element_geometry(mesh, e):
    """por ponto de Gauss: Nu, Np, B, ∇Np, w·detJ e a posição."""
    el = ELEMENTS[mesh['etype']]
    X = mesh['coords'][e]
    Xp = X[:el['np_']]
    out = []
    for xi, w in zip(*el['rule']()):
        Nu, dNu = el['shu'](xi)
        Np, dNp = el['shp'](xi)
        J = dNu @ X
        detJ = np.linalg.det(J)
        if detJ <= 0:
            raise ValueError('jacobiano não positivo')
        Ji = np.linalg.inv(J)
        Gu = Ji @ dNu
        Gp = Ji @ dNp
        out.append(dict(Nu=Nu, Np=Np, B=compute_B(Gu), Gp=Gp, wdJ=w * detJ, x=Nu @ X))
    return out


# ======================================================================================= materiais
def elastic_parameters(E, nu, sigma0=None):
    """material elástico linear (para verificações, p.ex. Terzaghi), no formato de parâmetros do Cam-Clay."""
    lam, G = E * nu / ((1 + nu) * (1 - 2 * nu)), E / (2 * (1 + nu))
    D = np.zeros((6, 6))                                   # Voigt do Cam-Clay {xx,xy,xz,yy,yz,zz}, γ = 2ε
    for i in (0, 3, 5):
        for j in (0, 3, 5):
            D[i, j] = lam + (2 * G if i == j else 0.0)
    for i in (1, 2, 4):
        D[i, i] = G
    s0 = np.zeros(6) if sigma0 is None else np.asarray(sigma0, float)
    return dict(model='elastic', E=E, nu=nu, D=D, sigma0=s0, K0=lam + 2 * G / 3, G0=G)


def material(P, eps, s_n):
    """tensão efetiva, tangente e variáveis internas no ponto de Gauss."""
    if P.get('model') == 'elastic':
        return dict(stress=P['sigma0'] + P['D'] @ eps, Dep=P['D'], plastic_strain=np.zeros(6), alpha=0.0,
                    plastic=False)
    return mcc.return_mapping(P, eps, s_n['epsp'], s_n['al'], sig_n=s_n['sig'], eps_n=s_n['eps'])


# ======================================================================================= modelo acoplado
class PoroField:
    """FE u-p com Cam-Clay por ponto de Gauss, estado inicial por elemento (centróide), peso próprio,
    poropressão inicial e termo gravitacional no fluxo."""

    def __init__(self, mesh, par_of_x, alpha=1.0, biot_modulus=np.inf, perm=0.0, gam_w=0.0, body=(0, 0, 0),
                 p0_of_x=lambda x: 0.0, init='gp'):
        """init = 'gp': estado inicial (σ'0, v0, K0) avaliado em cada ponto de Gauss (equilíbrio inicial exato
        para qualquer elemento); 'centroid': um estado por elemento, como as zonas do FLAC3D."""
        self.mesh = mesh
        el = ELEMENTS[mesh['etype']]
        self.nn = len(mesh['coords']); self.nu = 3 * self.nn
        self.pmap = -np.ones(self.nn, int); self.pmap[mesh['pnodes']] = np.arange(len(mesh['pnodes']))
        self.npd = len(mesh['pnodes'])
        self.geo = [element_geometry(mesh, e) for e in mesh['elements']]
        if init == 'centroid':
            pc = [par_of_x(mesh['coords'][e[:el['np_']]].mean(axis=0)) for e in mesh['elements']]
            self.par = [[P] * len(gps) for P, gps in zip(pc, self.geo)]
        else:
            self.par = [[par_of_x(g['x']) for g in gps] for gps in self.geo]
        self.state = [[dict(epsp=np.zeros(6), al=0.0, sig=P['sigma0'].copy(), eps=np.zeros(6)) for P in pars]
                      for pars in self.par]
        invM = 0.0 if np.isinf(biot_modulus) else 1.0 / biot_modulus
        rq, cq, vq, rs, cs, vs, vh = [], [], [], [], [], [], []
        self.fg = np.zeros(self.npd); self.Fb = np.zeros(self.nu)
        for e, gps in zip(mesh['elements'], self.geo):
            du = np.array([[3 * a, 3 * a + 1, 3 * a + 2] for a in e]).ravel()
            dp = self.pmap[e[:el['np_']]]
            Qe = sum(alpha * np.outer(g['B'].T @ M_FE, g['Np']) * g['wdJ'] for g in gps)
            Se = sum(invM * np.outer(g['Np'], g['Np']) * g['wdJ'] for g in gps)
            He = sum(perm * g['Gp'].T @ g['Gp'] * g['wdJ'] for g in gps)
            for g in gps:
                self.fg[dp] += perm * (g['Gp'].T @ np.array([0.0, 0.0, -gam_w])) * g['wdJ']
                self.Fb[du] += np.kron(g['Nu'], np.asarray(body, float)) * g['wdJ']
            rq.append(np.repeat(du, len(dp))); cq.append(np.tile(dp, len(du))); vq.append(Qe.ravel())
            rs.append(np.repeat(dp, len(dp))); cs.append(np.tile(dp, len(dp))); vs.append(Se.ravel()); vh.append(He.ravel())
        self.Q = sp.csr_matrix((np.concatenate(vq), (np.concatenate(rq), np.concatenate(cq))), shape=(self.nu, self.npd))
        self.S = sp.csr_matrix((np.concatenate(vs), (np.concatenate(rs), np.concatenate(cs))), shape=(self.npd, self.npd))
        self.H = sp.csr_matrix((np.concatenate(vh), (np.concatenate(rs), np.concatenate(cs))), shape=(self.npd, self.npd))
        self.U = np.zeros(self.nu)
        self.P = np.array([p0_of_x(mesh['coords'][v]) for v in mesh['pnodes']], float)
        self.trial = None

    def assemble(self, U):
        F = np.zeros(self.nu); rows, cols, vals = [], [], []; trial = []
        for e, gps, st, pars in zip(self.mesh['elements'], self.geo, self.state, self.par):
            du = np.array([[3 * a, 3 * a + 1, 3 * a + 2] for a in e]).ravel()
            ue = U[du]; n = len(du); fe = np.zeros(n); ke = np.zeros((n, n)); tr = []
            for g, s_n, P in zip(gps, st, pars):
                eps = (g['B'] @ ue)[FE2MCC]
                r = material(P, eps, s_n)
                r['eps'] = eps
                fe += g['B'].T @ r['stress'][MCC2FE] * g['wdJ']
                ke += g['B'].T @ r['Dep'][np.ix_(MCC2FE, MCC2FE)] @ g['B'] * g['wdJ']
                tr.append(r)
            F[du] += fe
            rows.append(np.repeat(du, n)); cols.append(np.tile(du, n)); vals.append(ke.ravel()); trial.append(tr)
        K = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(self.nu, self.nu))
        return F, K, trial

    def step(self, Fext, fixed_u, fixed_p, dt, flow, tol=1e-8, maxit=30, verbose=False, U_guess=None):
        """Um passo de Euler implícito com Newton. U_guess: chute inicial dos deslocamentos (preditor);
        o estado convergido do passo anterior (self.U, self.P) continua sendo o U_n, P_n do passo."""
        nu, npd = self.nu, self.npd
        U = (self.U if U_guess is None else np.asarray(U_guess, float)).copy(); P = self.P.copy()
        for i, v in fixed_u.items(): U[i] = v
        for i, v in fixed_p.items(): P[i] = v
        fixed = np.array(sorted(list(fixed_u) + [nu + i for i in fixed_p]), int)
        free = np.setdiff1d(np.arange(nu + npd), fixed)
        Hs = self.H if flow else 0.0 * self.H
        fgs = self.fg if flow else 0.0 * self.fg
        ref = max(np.linalg.norm(Fext), 1.0)
        for it in range(1, maxit + 1):
            F, K, trial = self.assemble(U)
            Ru = F - self.Q @ P - Fext
            Rp = self.Q.T @ (U - self.U) + self.S @ (P - self.P) + dt * (Hs @ P - fgs)
            R = np.r_[Ru, Rp]
            nr = np.linalg.norm(R[free]) / ref
            if verbose:
                print(f'      it {it}: |R|/|F| = {nr:.3e}')
            if nr < tol:
                break
            J = sp.bmat([[K, -self.Q], [self.Q.T, self.S + dt * Hs]], format='csc')
            dx = np.zeros(nu + npd)
            dx[free] = spla.spsolve(J[free][:, free], -R[free])
            U += dx[:nu]; P += dx[nu:]
        else:
            raise RuntimeError(f'Newton não convergiu (|R|/|F| = {nr:.3e})')
        self.U, self.P = U, P
        self.state = [[dict(epsp=r['plastic_strain'], al=r['alpha'], sig=r['stress'], eps=r['eps']) for r in tr]
                      for tr in trial]
        self.trial = trial
        self.Fint = F
        return it

    def reactions(self, Fext):
        """R = F_int(σ') - G P - F_ext: zero nos gdl livres; nas restrições, as reações de apoio."""
        F, _, _ = self.assemble(self.U)
        return F - self.Q @ self.P - Fext

    def zone_pp(self, x, z):
        """poropressão de 'zona' (média dos 8 vértices da célula da malha estruturada que contém (x, z))."""
        C = self.mesh['coords']
        for c in self.mesh['cells']:
            X = C[c]
            if X[:, 0].min() <= x <= X[:, 0].max() and X[:, 2].min() <= z <= X[:, 2].max():
                return self.P[self.pmap[c]].mean()
