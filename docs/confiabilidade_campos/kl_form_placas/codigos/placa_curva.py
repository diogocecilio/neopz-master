"""Placa curva (Ghanem & Spanos, Fig. 5.45) com modulo de elasticidade aleatorio.

Compara tres bases KL para o campo E(x) = 1 + sigma * sum_k sqrt(lam_k) f_k(x) xi_k:
  - Nystrom (quadratura nodal, como no codigo Mathematica do usuario)
  - Galerkin (metodo do livro, eq. 5.61-5.68)
  - KL analitico do quadrado envolvente [0,1.25]^2 restrito a placa
e calcula sigma de u_x(A), Pf por FORM (HLRF) e por Monte Carlo.
"""
import time
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import brentq
from scipy.stats import norm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

t0 = time.time()
LP = 1.25          # lados retos
LC = 1.0           # comprimento de correlacao (x e y)
NT, NS = 24, 4     # elementos ao longo do arco e na espessura


# ---------------------------------------------------------------- malha
def make_mesh(nt, ns):
    """Malha Q4 mapeada: interpola entre o arco interno e o contorno externo em L."""
    xy = np.zeros(((nt + 1) * (ns + 1), 2))
    for i in range(nt + 1):
        t = i / nt
        th = t * np.pi / 2
        pin = np.array([LP - np.cos(th), np.sin(th)])            # arco de raio 1, centro (1.25, 0)
        d = 2 * LP * t                                            # comprimento ao longo do L externo
        pout = np.array([0.0, d]) if d <= LP else np.array([d - LP, LP])
        for j in range(ns + 1):
            s = j / ns
            xy[i * (ns + 1) + j] = (1 - s) * pin + s * pout
    nid = lambda i, j: i * (ns + 1) + j
    els = []
    for i in range(nt):
        for j in range(ns):
            e = [nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1)]
            x = xy[e]
            area2 = np.sum(x[:, 0] * np.roll(x[:, 1], -1) - np.roll(x[:, 0], -1) * x[:, 1])
            els.append(e if area2 > 0 else e[::-1])               # sentido anti-horario
    base = [nid(0, j) for j in range(ns + 1)]                     # t=0: y=0, 0<=x<=0.25
    right = [nid(nt, j) for j in range(ns + 1)]                   # t=1: x=1.25, 1<=y<=1.25
    nodeA = nid(nt, ns)                                           # (1.25, 1.25)
    return xy, np.array(els), base, right, nodeA


# ---------------------------------------------------------------- FEM Q4
gpt = np.array([-1, 1]) / np.sqrt(3)
GAUSS = [(a, b) for a in gpt for b in gpt]


def shape(r, s):
    N = 0.25 * np.array([(1 - r) * (1 - s), (1 + r) * (1 - s), (1 + r) * (1 + s), (1 - r) * (1 + s)])
    dN = 0.25 * np.array([[-(1 - s), (1 - s), (1 + s), -(1 + s)],
                          [-(1 - r), -(1 + r), (1 + r), (1 - r)]])
    return N, dN


def fem_setup(xy, els):
    """Retorna G[n] = rigidez da 'funcao chapeu' do no n (K(E) = sum_n E_n G[n]) e os pontos de Gauss."""
    nn = len(xy)
    nd = 2 * nn
    D1 = np.diag([1.0, 1.0, 0.5])               # E=1, nu=0 (plano de tensao = plano de deformacao)
    G = np.zeros((nn, nd, nd))
    gx, gw, gN = [], [], []
    for e in els:
        x = xy[e]
        dofs = np.ravel([[2 * n, 2 * n + 1] for n in e])
        for (r, s) in GAUSS:
            N, dN = shape(r, s)
            J = dN @ x
            dJ = np.linalg.det(J)
            dNx = np.linalg.solve(J, dN)
            B = np.zeros((3, 8))
            B[0, 0::2] = dNx[0]
            B[1, 1::2] = dNx[1]
            B[2, 0::2] = dNx[1]
            B[2, 1::2] = dNx[0]
            kg = B.T @ D1 @ B * dJ
            for a in range(4):
                G[e[a]][np.ix_(dofs, dofs)] += N[a] * kg
            gx.append(N @ x)
            gw.append(dJ)
            row = np.zeros(nn)
            row[e] = N
            gN.append(row)
    return G, np.array(gx), np.array(gw), np.array(gN)


def kernel(p, q):
    return np.exp(-np.abs(p[:, None, 0] - q[None, :, 0]) / LC - np.abs(p[:, None, 1] - q[None, :, 1]) / LC)


xy, els, base, right, nodeA = make_mesh(NT, NS)
nn = len(xy)
G, gx, gw, gN = fem_setup(xy, els)
fixed = np.ravel([[2 * n, 2 * n + 1] for n in base])
free = np.setdiff1d(np.arange(2 * nn), fixed)
Gf = G[:, free][:, :, free]                       # (nn, nf, nf)
nf = len(free)

# forca total P = 1 em x, distribuida uniformemente na borda direita
F = np.zeros(2 * nn)
yr = xy[right, 1]
order = np.argsort(yr)
rn = np.array(right)[order]
Ltot = yr.max() - yr.min()
for k in range(len(rn) - 1):
    L = xy[rn[k + 1], 1] - xy[rn[k], 1]
    F[2 * rn[k]] += 0.5 * L / Ltot
    F[2 * rn[k + 1]] += 0.5 * L / Ltot
Ff = F[free]
posA = int(np.where(free == 2 * nodeA)[0][0])     # u_x em A (longitudinal)


# padrao de esparsidade: K(E)[I,J] = Gv @ E  (Gv: nnz x nn) -> montagem barata para muitas amostras
I, J = np.nonzero(np.abs(Gf).sum(0) > 0)
Gv = Gf[:, I, J].T.copy()                         # (nnz, nn)


def K_of_E(E):
    K = np.zeros((nf, nf))
    K[I, J] = Gv @ E
    return K


uA_det = np.linalg.solve(K_of_E(np.ones(nn)), Ff)[posA]
area = gw.sum()
print(f"nos={nn} elementos={len(els)} gdl livres={nf}  area={area:.4f} (exata {LP**2 - np.pi / 4:.4f})")
print(f"u_x(A) deterministico = {uA_det:.3f}")

# convergencia de malha do deterministico
for (a, b) in [(48, 8)]:
    xy2, els2, base2, right2, A2 = make_mesh(a, b)
    G2, _, _, _ = fem_setup(xy2, els2)
    fx2 = np.ravel([[2 * n, 2 * n + 1] for n in base2])
    fr2 = np.setdiff1d(np.arange(2 * len(xy2)), fx2)
    F2 = np.zeros(2 * len(xy2))
    r2 = np.array(right2)[np.argsort(xy2[right2, 1])]
    for k in range(len(r2) - 1):
        L = xy2[r2[k + 1], 1] - xy2[r2[k], 1]
        F2[2 * r2[k]] += 0.5 * L / 0.25
        F2[2 * r2[k + 1]] += 0.5 * L / 0.25
    K2 = np.tensordot(np.ones(len(xy2)), G2[:, fr2][:, :, fr2], axes=1)
    u2 = np.linalg.solve(K2, F2[fr2])
    print(f"   malha {a}x{b}: u_x(A) = {u2[np.where(fr2 == 2 * A2)[0][0]]:.3f}")

# ---------------------------------------------------------------- 1. Nystrom
w = np.zeros(nn)
for e in els:
    x = xy[e]
    ar = 0.5 * abs(np.sum(x[:, 0] * np.roll(x[:, 1], -1) - np.roll(x[:, 0], -1) * x[:, 1]))
    w[e] += ar / 4
C = kernel(xy, xy)
sw = np.sqrt(w)
lamN, V = np.linalg.eigh(sw[:, None] * C * sw[None, :])
idx = np.argsort(lamN)[::-1]
lamN, phiN = lamN[idx], (V[:, idx] / sw[:, None]).T              # phiN[k, n]

# ---------------------------------------------------------------- 2. Galerkin
Cq = kernel(gx, gx)
Cg = gN.T @ (gw[:, None] * Cq * gw[None, :]) @ gN
Bg = gN.T @ (gw[:, None] * gN)
lamG, Vg = eigh(Cg, Bg)                                          # d^T B d = 1
idx = np.argsort(lamG)[::-1]
lamG, phiG = lamG[idx], Vg[:, idx].T

# ---------------------------------------------------------------- 3. analitico do quadrado [0,1.25]^2
a = LP / 2
c = 1 / LC
NW = 25
om = []
for i in range(1, NW + 1):
    lo, hi = (i - 1) * np.pi / (2 * a) + 1e-9, i * np.pi / (2 * a) - 1e-9
    f = (lambda x: c - x * np.tan(x * a)) if i % 2 == 1 else (lambda x: x + c * np.tan(x * a))
    om.append(brentq(f, lo, hi))
om = np.array(om)
lam1 = 2 * c / (om ** 2 + c ** 2)


def f1(i, x):                       # i base 0; x centrado em [-a, a]
    o = om[i]
    if i % 2 == 0:
        return np.cos(o * x) / np.sqrt(a + np.sin(2 * o * a) / (2 * o))
    return np.sin(o * x) / np.sqrt(a - np.sin(2 * o * a) / (2 * o))


pairs = sorted(((i, j) for i in range(NW) for j in range(NW)), key=lambda p: -lam1[p[0]] * lam1[p[1]])
pairs = pairs[:nn]
lamQ = np.array([lam1[i] * lam1[j] for i, j in pairs])
phiQ = np.array([f1(i, xy[:, 0] - a) * f1(j, xy[:, 1] - a) for i, j in pairs])

bases = {"Nystrom": (lamN, phiN), "Galerkin": (lamG, phiG), "Analitico quadrado": (lamQ, phiQ)}


def Amat(name, M):
    lam, phi = bases[name]
    return (np.sqrt(np.maximum(lam[:M], 0))[:, None] * phi[:M]).T    # (nn, M)


print("\n--- autovalores (8 primeiros) ---")
for k, (lam, _) in bases.items():
    print(f"{k:20s}", np.round(lam[:8], 4), f" soma={lam.sum():.3f}")
print(f"(area da placa = {area:.3f}; o quadrado tem area {LP**2:.3f})")

# ortogonalidade das autofuncoes NA PLACA (Gram com pesos de Nystrom)
for k, (lam, phi) in bases.items():
    Gm = (phi[:6] * w) @ phi[:6].T
    print(f"{k:20s} max |<f_i,f_j>_placa - delta_ij| (6 modos) = {np.abs(Gm - np.eye(6)).max():.3f}")

# ---------------------------------------------------------------- qualidade da representacao
Ms = [2, 4, 6, 10, 20, 40, nn]
print("\n--- variancia media de E na placa capturada por M modos (exato = 1) ---")
print("M     " + "".join(f"{k:>20s}" for k in bases))
rep = {}
for M in Ms:
    row = []
    for k in bases:
        A = Amat(k, M)
        var = (A ** 2).sum(1)
        CM = A @ A.T
        errC = np.linalg.norm(CM - C) / np.linalg.norm(C)
        row.append((np.sum(w * var) / area, var[nodeA], errC))
    rep[M] = row
    print(f"{M:<6d}" + "".join(f"   {r[0]:.3f} A={r[1]:.2f} eC={r[2]:.3f}" for r in row))


# ---------------------------------------------------------------- MC e FORM
def sample_uA(sig, A, xis, batch=1000):
    out = np.empty(len(xis))
    Emin = np.empty(len(xis))
    for b in range(0, len(xis), batch):
        E = 1 + sig * xis[b:b + batch] @ A.T                          # (B, nn)
        K = np.zeros((len(E), nf, nf))
        K[:, I, J] = E @ Gv.T                                         # (B, nf, nf)
        u = np.linalg.solve(K, np.broadcast_to(Ff, (len(E), nf))[..., None])[..., 0]
        out[b:b + batch] = u[:, posA]
        Emin[b:b + batch] = E.min(1)
    return out, Emin


def u_grad(xi, sig, A):
    E = 1 + sig * A @ xi
    K = K_of_E(E)
    u = np.linalg.solve(K, Ff)
    lamb = np.linalg.solve(K, np.eye(nf)[posA])                    # adjunto (K simetrica)
    dudE = -(lamb[I] * u[J]) @ Gv
    return u[posA], sig * A.T @ dudE


def form_hlrf(delta, sig, A, maxit=30, tol=1e-8):
    """HLRF no espaco dos xi (normais padrao independentes): g = delta - u_A."""
    y = np.zeros(A.shape[1])
    for it in range(maxit):
        uA, du = u_grad(y, sig, A)
        g, dg = delta - uA, -du
        ynew = (dg @ y - g) / (dg @ dg) * dg
        if np.linalg.norm(ynew - y) < tol:
            y = ynew
            break
        y = ynew
    beta = np.linalg.norm(y)
    alpha2 = (dg / np.linalg.norm(dg)) ** 2
    return beta, norm.cdf(-beta), y, alpha2, it + 1


# checagem do gradiente adjunto
A4 = Amat("Nystrom", 4)
x0 = np.array([0.3, -0.5, 0.8, 0.1])
gadj = u_grad(x0, 0.2, A4)[1]
gfd = np.array([(u_grad(x0 + 1e-6 * np.eye(4)[k], 0.2, A4)[0] - u_grad(x0, 0.2, A4)[0]) / 1e-6 for k in range(4)])
print("\ngradiente adjunto:", np.round(gadj, 5), "\ndiferenca finita :", np.round(gfd, 5))

rng = np.random.default_rng(2026)

# --- sigma_u x sigma_E (livro Fig. 5.66, M=4, sigma_max=19.4)
book = np.array([[0.02, 0.04], [0.04, 0.085], [0.06, 0.125], [0.08, 0.17], [0.10, 0.215], [0.12, 0.26],
                 [0.14, 0.31], [0.16, 0.36], [0.18, 0.415], [0.20, 0.475], [0.22, 0.55], [0.24, 0.63],
                 [0.26, 0.72], [0.28, 0.84]])
xiS = rng.standard_normal((5000, nn))
sigs = np.array([0.04, 0.08, 0.12, 0.16, 0.20, 0.24, 0.28])
print("\n--- sigma de u_x(A) / u_det (MC 5000) ---")
print("sigma_E   Nys M=4  Gal M=4  Quad M=4  Nys todos  Quad todos")
stdtab = []
for s in sigs:
    row = [s]
    for k, M in [("Nystrom", 4), ("Galerkin", 4), ("Analitico quadrado", 4), ("Nystrom", nn), ("Analitico quadrado", nn)]:
        u, _ = sample_uA(s, Amat(k, M), xiS[:, :M])
        row.append(u.std() / uA_det)
    stdtab.append(row)
    print("  ".join(f"{v:8.3f}" for v in row))
stdtab = np.array(stdtab)

# --- Pf: FORM x MC
sigE = 0.2
facs = [1.3, 1.5, 1.7]
NMC = 20000
xiP = rng.standard_normal((NMC, nn))
print(f"\n--- Pf = P[u_x(A) > fator*u_det], sigma_E={sigE}, MC N={NMC} ---")
print(f"{'base':20s} {'M':>4s} " + "".join(f"  FORM({f})  MC({f})" for f in facs) + "   beta(1.5) iters")
pf_rows = []
for k in bases:
    for M in [2, 4, 10, 20, nn]:
        A = Amat(k, M)
        u, Emin = sample_uA(sigE, A, xiP[:, :M])
        line, rec = f"{k:20s} {M:4d} ", {"base": k, "M": M}
        for f in facs:
            beta, pf, y, a2, it = form_hlrf(f * uA_det, sigE, A)
            pmc = np.mean((u >= f * uA_det) | (Emin <= 0))
            line += f"  {pf:9.4f} {pmc:7.4f}"
            rec[f] = (pf, pmc)
            if f == 1.5:
                b15, it15 = beta, it
        print(line + f"   {b15:.3f} {it15}")
        pf_rows.append(rec)

beta, pf, y, a2, it = form_hlrf(1.5 * uA_det, sigE, Amat("Nystrom", nn))
print("\nalpha^2 (Nystrom, todos, fator 1.5), 8 maiores modos:", np.round(a2[:8], 3), " soma 8:", round(a2[:8].sum(), 3))
beta, pf, y, a2, it = form_hlrf(1.5 * uA_det, sigE, Amat("Analitico quadrado", nn))
print("alpha^2 (Quadrado, todos, fator 1.5), 8 primeiros modos:", np.round(a2[:8], 3), " soma 8:", round(a2[:8].sum(), 3))

# ---------------------------------------------------------------- figuras
tri = mtri.Triangulation(xy[:, 0], xy[:, 1], np.vstack([els[:, [0, 1, 2]], els[:, [0, 2, 3]]]))
fig, ax = plt.subplots(3, 4, figsize=(15, 11))
for r, k in enumerate(bases):
    lam, phi = bases[k]
    for m in range(4):
        sgn = np.sign(phi[m][np.argmax(np.abs(phi[m]))])
        tc = ax[r, m].tripcolor(tri, sgn * phi[m], shading="gouraud", cmap="RdBu_r")
        ax[r, m].set_aspect("equal")
        ax[r, m].set_title(f"{k}: f{m + 1}, λ={lam[m]:.3f}", fontsize=9)
        plt.colorbar(tc, ax=ax[r, m], shrink=0.7)
plt.tight_layout()
plt.savefig("curva_modos.png", dpi=110)

fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))
for i, k in enumerate(bases):
    var = (Amat(k, 10) ** 2).sum(1)
    tc = ax[i].tripcolor(tri, var, shading="gouraud", vmin=0, vmax=1, cmap="viridis")
    ax[i].set_aspect("equal")
    ax[i].set_title(f"{k}: variância de E com M=10")
    plt.colorbar(tc, ax=ax[i])
plt.tight_layout()
plt.savefig("curva_variancia.png", dpi=110)

fig, ax = plt.subplots(1, 2, figsize=(14, 5))
for k, mk in zip(bases, ["o", "s", "^"]):
    ax[0].plot(Ms, [rep[M][list(bases).index(k)][0] for M in Ms], mk + "-", label=k)
ax[0].set_xscale("log"); ax[0].set_xlabel("M (número de variáveis aleatórias)")
ax[0].set_ylabel("variância média de E capturada na placa"); ax[0].grid(True, which="both"); ax[0].legend()
for k, mk in zip(bases, ["o", "s", "^"]):
    rows = [r for r in pf_rows if r["base"] == k]
    ax[1].plot([r["M"] for r in rows], [r[1.5][0] for r in rows], mk + "-", label=f"FORM {k}")
    ax[1].plot([r["M"] for r in rows], [r[1.5][1] for r in rows], mk + "--", mfc="none", label=f"MC {k}")
ax[1].set_xscale("log"); ax[1].set_xlabel("M"); ax[1].set_ylabel("Pf  (δ = 1.5 u_det)")
ax[1].set_title(f"σ_E = {sigE}"); ax[1].grid(True, which="both"); ax[1].legend(fontsize=8)
plt.tight_layout()
plt.savefig("curva_convergencia.png", dpi=110)

fig, ax = plt.subplots(figsize=(7, 5))
ax.plot(book[:, 0], book[:, 1] * 19.4 / uA_det, "k*", ms=10, label="Livro MCS 500 (Fig. 5.66) × 19.4 / u_det")
for j, lab in enumerate(["Nystrom M=4", "Galerkin M=4", "Quadrado M=4", "Nystrom todos", "Quadrado todos"], start=1):
    ax.plot(stdtab[:, 0], stdtab[:, j], "o-" if j < 4 else "s--", label=lab)
ax.set_xlabel("σ_E"); ax.set_ylabel("σ de u_x(A) / u_det"); ax.grid(True); ax.legend(fontsize=8)
plt.tight_layout()
plt.savefig("curva_sigma.png", dpi=110)
print(f"\ntempo total {time.time() - t0:.1f} s")
