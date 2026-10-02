"""
poro_camclay_fem.py
===================

Elementos finitos acoplados u-p (Biot) em hexaedros de 8 nós com o modelo Cam-Clay modificado
(modified_cam_clay.py) como lei constitutiva da tensão efetiva.

É a adaptação, para plasticidade, da estrutura do código poroelástico (poro-elastic.nb):
  * mesmas matrizes de acoplamento  Q = α ∫ Bᵀ m Np dΩ,  S = ∫ (1/M_b) Npᵀ Np dΩ,
    H = ∫ (k/μ) ∇Npᵀ ∇Np dΩ  (ContributePorous);
  * mesma matriz B (ComputeBN, ordem {xx, yy, zz, γyz, γxz, γxy}) e interpolação de mesma ordem
    para u e p;
  * mesmo sistema de Euler implícito  [[K, −Q], [Qᵀ, S + Δt H]];
  mas agora não linear: K → tangente consistente, K u → forças internas ∫ Bᵀ σ' dΩ, com
  Newton-Raphson em cada passo e estado (ε^p, α) guardado em cada ponto de Gauss.

Convenções: tração positiva; pressão de poros p > 0 em compressão; σ = σ' − α p m.
"""
import numpy as np
import modified_cam_clay as mcc

# ------------------------------------------------------------------------------------------
# Elemento hexaédrico de 8 nós (ordem de nós do HexahedronElement do Mathematica)
# ------------------------------------------------------------------------------------------
HEX8_NODES = np.array([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
                       [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]], float)
_g = 1.0 / np.sqrt(3.0)
GAUSS_HEX8 = [(np.array([a, b, c]) * _g, 1.0) for c in (-1, 1) for b in (-1, 1) for a in (-1, 1)]
GAUSS_QUAD4 = [(np.array([a, b]) * _g, 1.0) for b in (-1, 1) for a in (-1, 1)]

# FE (ComputeBN): {xx, yy, zz, yz, xz, xy}  <->  MCC/HWTools: {xx, xy, xz, yy, yz, zz}
FE2MCC = [0, 5, 4, 1, 3, 2]      # eps_mcc = eps_fe[FE2MCC]
MCC2FE = [0, 3, 5, 4, 2, 1]      # sig_fe  = sig_mcc[MCC2FE]
M_FE = np.array([1., 1., 1., 0., 0., 0.])


def shape_hex8(xi):
    """N (8,) e dN/dξ (3, 8) no ponto ξ = (ξ, η, ζ)."""
    s = HEX8_NODES
    N = 0.125 * (1 + s[:, 0] * xi[0]) * (1 + s[:, 1] * xi[1]) * (1 + s[:, 2] * xi[2])
    dN = np.array([0.125 * s[:, 0] * (1 + s[:, 1] * xi[1]) * (1 + s[:, 2] * xi[2]),
                   0.125 * s[:, 1] * (1 + s[:, 0] * xi[0]) * (1 + s[:, 2] * xi[2]),
                   0.125 * s[:, 2] * (1 + s[:, 0] * xi[0]) * (1 + s[:, 1] * xi[1])])
    return N, dN


def shape_quad4(xi):
    s = np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]], float)
    N = 0.25 * (1 + s[:, 0] * xi[0]) * (1 + s[:, 1] * xi[1])
    dN = np.array([0.25 * s[:, 0] * (1 + s[:, 1] * xi[1]), 0.25 * s[:, 1] * (1 + s[:, 0] * xi[0])])
    return N, dN


def compute_B(gradphi):
    """ComputeBN (dim = 3) do código original: linhas {xx, yy, zz, γyz, γxz, γxy}."""
    n = gradphi.shape[1]
    B = np.zeros((6, 3 * n))
    for i in range(n):
        dx, dy, dz = gradphi[:, i]
        B[0, 3 * i] = dx
        B[1, 3 * i + 1] = dy
        B[2, 3 * i + 2] = dz
        B[3, 3 * i + 1] = dz; B[3, 3 * i + 2] = dy
        B[4, 3 * i] = dz;     B[4, 3 * i + 2] = dx
        B[5, 3 * i] = dy;     B[5, 3 * i + 1] = dx
    return B


# ------------------------------------------------------------------------------------------
# Malha estruturada de uma caixa (marcadores: 1 x=0, 2 x=Lx, 3 y=0, 4 y=Ly, 5 z=0, 6 z=Lz)
# ------------------------------------------------------------------------------------------
def box_mesh(L=(1., 1., 1.), n=(1, 1, 1)):
    Lx, Ly, Lz = L; nx, ny, nz = n
    xs, ys, zs = np.linspace(0, Lx, nx + 1), np.linspace(0, Ly, ny + 1), np.linspace(0, Lz, nz + 1)
    idx = lambda i, j, k: k * (ny + 1) * (nx + 1) + j * (nx + 1) + i
    coords = np.array([[xs[i], ys[j], zs[k]] for k in range(nz + 1) for j in range(ny + 1) for i in range(nx + 1)])
    els = []
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                els.append([idx(i, j, k), idx(i + 1, j, k), idx(i + 1, j + 1, k), idx(i, j + 1, k),
                            idx(i, j, k + 1), idx(i + 1, j, k + 1), idx(i + 1, j + 1, k + 1), idx(i, j + 1, k + 1)])
    faces = []
    for k in range(nz):
        for j in range(ny):
            faces.append(([idx(0, j, k), idx(0, j, k + 1), idx(0, j + 1, k + 1), idx(0, j + 1, k)], 1))
            faces.append(([idx(nx, j, k), idx(nx, j + 1, k), idx(nx, j + 1, k + 1), idx(nx, j, k + 1)], 2))
    for k in range(nz):
        for i in range(nx):
            faces.append(([idx(i, 0, k), idx(i + 1, 0, k), idx(i + 1, 0, k + 1), idx(i, 0, k + 1)], 3))
            faces.append(([idx(i, ny, k), idx(i, ny, k + 1), idx(i + 1, ny, k + 1), idx(i + 1, ny, k)], 4))
    for j in range(ny):
        for i in range(nx):
            faces.append(([idx(i, j, 0), idx(i, j + 1, 0), idx(i + 1, j + 1, 0), idx(i + 1, j, 0)], 5))
            faces.append(([idx(i, j, nz), idx(i + 1, j, nz), idx(i + 1, j + 1, nz), idx(i, j + 1, nz)], 6))
    return dict(coords=coords, elements=np.array(els), faces=faces)


def nodes_on(mesh, marker):
    return sorted({n for f, m in mesh['faces'] if m == marker for n in f})


def face_load(mesh, marker, t):
    """Forças nodais consistentes de uma tração uniforme t (3,) nas faces com o marcador."""
    F = np.zeros(3 * len(mesh['coords']))
    for f, m in mesh['faces']:
        if m != marker:
            continue
        X = mesh['coords'][f]
        for xi, w in GAUSS_QUAD4:
            N, dN = shape_quad4(xi)
            J = dN @ X                               # (2, 3)
            dA = np.linalg.norm(np.cross(J[0], J[1]))
            for a, node in enumerate(f):
                F[3 * node:3 * node + 3] += N[a] * np.asarray(t) * w * dA
    return F


# ------------------------------------------------------------------------------------------
# Contribuições de elemento
# ------------------------------------------------------------------------------------------
def element_geometry(X):
    """Pré-calcula, por ponto de Gauss: N, B, ∇N e peso·detJ."""
    gps = []
    for xi, w in GAUSS_HEX8:
        N, dN = shape_hex8(xi)
        J = dN @ X                                   # (3, 3)
        detJ = np.linalg.det(J)
        gradphi = np.linalg.solve(J, dN)             # (3, 8)
        gps.append(dict(N=N, B=compute_B(gradphi), G=gradphi, wdJ=w * abs(detJ)))
    return gps


def contribute_porous(gps, alpha, inv_mb, perm_mu):
    """Q, S, H do elemento (mesmas expressões de ContributePorous)."""
    n = len(gps[0]['N'])
    Q = np.zeros((3 * n, n)); S = np.zeros((n, n)); H = np.zeros((n, n))
    for g in gps:
        Q += alpha * np.outer(g['B'].T @ M_FE, g['N']) * g['wdJ']
        S += inv_mb * np.outer(g['N'], g['N']) * g['wdJ']
        H += perm_mu * (g['G'].T @ g['G']) * g['wdJ']
    return Q, S, H


def contribute_plasticity(gps, ue, state_n, P):
    """Forças internas efetivas e tangente consistente do elemento (Cam-Clay modificado)."""
    n = len(gps[0]['N'])
    fe = np.zeros(3 * n); ke = np.zeros((3 * n, 3 * n)); trial = []
    for g, (epsp_n, al_n) in zip(gps, state_n):
        eps_fe = g['B'] @ ue
        r = mcc.return_mapping(P, eps_fe[FE2MCC], epsp_n, al_n)
        sig = r['stress'][MCC2FE]
        D = r['Dep'][np.ix_(MCC2FE, MCC2FE)]
        fe += g['B'].T @ sig * g['wdJ']
        ke += g['B'].T @ D @ g['B'] * g['wdJ']
        trial.append((r['plastic_strain'], r['alpha'], r))
    return fe, ke, trial


# ------------------------------------------------------------------------------------------
# Problema acoplado
# ------------------------------------------------------------------------------------------
class PoroCamClay:
    def __init__(self, mesh, P, alpha=1.0, biot_modulus=np.inf, perm_mu=0.0):
        self.mesh, self.P = mesh, P
        self.alpha, self.perm_mu = alpha, perm_mu
        self.inv_mb = 0.0 if np.isinf(biot_modulus) else 1.0 / biot_modulus
        self.nn = len(mesh['coords']); self.ndof_u = 3 * self.nn
        self.geo = [element_geometry(mesh['coords'][e]) for e in mesh['elements']]
        self.state = [[(np.zeros(6), 0.0) for _ in g] for g in self.geo]
        self.U = np.zeros(self.ndof_u); self.Pp = np.zeros(self.nn)
        # Q, S, H não dependem da plasticidade (pequenas deformações): montar uma vez
        nt = self.ndof_u + self.nn
        self.Qg = np.zeros((self.ndof_u, self.nn)); self.Sg = np.zeros((self.nn, self.nn))
        self.Hg = np.zeros((self.nn, self.nn))
        for e, gps in zip(mesh['elements'], self.geo):
            Q, S, H = contribute_porous(gps, alpha, self.inv_mb, perm_mu)
            du = np.array([[3 * a, 3 * a + 1, 3 * a + 2] for a in e]).ravel()
            self.Qg[np.ix_(du, e)] += Q; self.Sg[np.ix_(e, e)] += S; self.Hg[np.ix_(e, e)] += H

    def assemble(self, U):
        F = np.zeros(self.ndof_u); K = np.zeros((self.ndof_u, self.ndof_u)); trial = []
        for e, gps, st in zip(self.mesh['elements'], self.geo, self.state):
            du = np.array([[3 * a, 3 * a + 1, 3 * a + 2] for a in e]).ravel()
            fe, ke, tr = contribute_plasticity(gps, U[du], st, self.P)
            F[du] += fe; K[np.ix_(du, du)] += ke; trial.append(tr)
        return F, K, trial

    def step(self, Fext, fixed_u, fixed_p, dt=1.0, tol=1e-10, maxit=25, verbose=False):
        """Um passo de Euler implícito. fixed_u: {dof: valor total}, fixed_p: {nó: valor}."""
        nu, nn = self.ndof_u, self.nn
        U = self.U.copy(); Pp = self.Pp.copy()
        for d, v in fixed_u.items(): U[d] = v
        for i, v in fixed_p.items(): Pp[i] = v
        fixed = np.array(sorted(list(fixed_u.keys()) + [nu + i for i in fixed_p.keys()]), int)
        free = np.setdiff1d(np.arange(nu + nn), fixed)
        ref = max(np.linalg.norm(Fext), 1.0)
        for it in range(1, maxit + 1):
            Fint, K, trial = self.assemble(U)
            Ru = Fint - self.Qg @ Pp - Fext
            Rp = self.Qg.T @ (U - self.U) + self.Sg @ (Pp - self.Pp) + dt * self.Hg @ Pp
            R = np.r_[Ru, Rp]
            nr = np.linalg.norm(R[free]) / ref
            if verbose:
                print(f'      it {it}: |R|/|F| = {nr:.3e}')
            if nr < tol:
                break
            J = np.block([[K, -self.Qg], [self.Qg.T, self.Sg + dt * self.Hg]])
            dx = np.zeros(nu + nn)
            dx[free] = np.linalg.solve(J[np.ix_(free, free)], -R[free])
            U += dx[:nu]; Pp += dx[nu:]
        else:
            raise RuntimeError('Newton global não convergiu')
        self.U, self.Pp = U, Pp
        self.state = [[(tr[0], tr[1]) for tr in tre] for tre in trial]
        return it, trial


# ------------------------------------------------------------------------------------------
# Ensaio triaxial (benchmark Itasca): 1 elemento ou malha nx×ny×nz
# ------------------------------------------------------------------------------------------
def triaxial_test(P, drained=True, ea_max=0.5, nsteps=500, n=(1, 1, 1), L=(1., 1., 1.),
                  biot_modulus=None, sigma_cell=None, verbose=False):
    mesh = box_mesh(L, n)
    if sigma_cell is None:
        sigma_cell = P['sigma0'][0]                          # tensão total lateral (tração +)
    if drained:
        model = PoroCamClay(mesh, P, alpha=1.0, biot_modulus=np.inf, perm_mu=0.0)
    else:
        model = PoroCamClay(mesh, P, alpha=1.0, biot_modulus=biot_modulus, perm_mu=0.0)
    Fext = face_load(mesh, 2, [sigma_cell, 0, 0]) + face_load(mesh, 4, [0, sigma_cell, 0])
    fix0 = {}
    for nd in nodes_on(mesh, 1): fix0[3 * nd] = 0.0          # ux = 0 em x = 0 (simetria)
    for nd in nodes_on(mesh, 3): fix0[3 * nd + 1] = 0.0      # uy = 0 em y = 0
    for nd in nodes_on(mesh, 5): fix0[3 * nd + 2] = 0.0      # uz = 0 na base
    top = nodes_on(mesh, 6)
    fixed_p = {i: 0.0 for i in range(model.nn)} if drained else {}
    Lz = L[2]; v0 = P['v0']; nel = len(mesh['elements'])
    hist = []

    def record(ea):
        # médias sobre todos os pontos de Gauss
        sig = np.mean([np.array(r[2]['stress']) for tre in last for r in tre], axis=0) if last else P['sigma0']
        p, q = mcc.invariants(sig)
        eps = np.mean([g['B'] @ model.U[np.array([[3 * a, 3 * a + 1, 3 * a + 2] for a in e]).ravel()]
                       for e, gps in zip(mesh['elements'], model.geo) for g in gps], axis=0)
        ev = eps[0] + eps[1] + eps[2]
        u = float(np.mean(model.Pp))
        pc = np.mean([r[2]['pc'] for tre in last for r in tre]) if last else P['pc0']
        hist.append(dict(ea=ea, p=-p, q=q, v=v0 * (1 + ev), u=u, ptot=-p + u, pc=pc, ev=-ev))

    last = None
    record(0.0)
    iters = []
    for k in range(1, nsteps + 1):
        ea = ea_max * k / nsteps
        fu = dict(fix0)
        for nd in top: fu[3 * nd + 2] = -ea * Lz
        it, last = model.step(Fext, fu, fixed_p, verbose=False)
        iters.append(it)
        record(ea)
    if verbose:
        print(f'   {nel} elemento(s), {nsteps} passos, iterações de Newton por passo: máx {max(iters)}, '
              f'média {np.mean(iters):.2f}')
    return {k: np.array([h[k] for h in hist]) for k in hist[0]}, iters
