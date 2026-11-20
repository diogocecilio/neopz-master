"""Percolacao permanente saturada (Eqs. 19-21 do artigo) em P1, para
diagnostico: resolve div(K grad u) = 0 na mesma malha da exportacao, com
K = k (e_y e_y + alfa (1 - e_y e_y)) e k nodal (multiplicador).

Contornos (sistema do artigo, y para baixo, nivel final y = hw):
  * crista (y = 0, x <= 0): u = 0;
  * face exposta (0 <= y <= hw): u = -gamma_w y;
  * face submersa, terreno a frente do pe e lateral direita abaixo do
    nivel: u = -gamma_w hw;
  * lateral esquerda: u = 0 (lateral='dirichlet', nivel inicial mantido,
    como no transiente do NeoPZ) ou fluxo nulo (lateral='neumann');
  * base: fluxo nulo.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl

from cinematico import GAMMA_W, CampoMalha


def resolver(xy, tri, talude, hw=None, alfa=1.0, k=None, lateral="dirichlet", tol=1e-7):
    xy = np.asarray(xy, float)
    tri = np.asarray(tri, int)
    n = len(xy)
    hw = talude.H if hw is None else hw
    k = np.ones(n) if k is None else np.asarray(k, float)
    P = xy[tri]
    B = np.stack([P[:, 1] - P[:, 0], P[:, 2] - P[:, 0]], axis=2)    # (nt,2,2) colunas
    det = B[:, 0, 0] * B[:, 1, 1] - B[:, 0, 1] * B[:, 1, 0]
    A = 0.5 * np.abs(det)
    Binv = np.linalg.inv(B)                                          # (nt,2,2)
    Gref = np.array([[-1.0, 1.0, 0.0], [-1.0, 0.0, 1.0]])
    G = np.einsum("kji,jl->kil", Binv, Gref)                         # grad das bases (nt,2,3)
    Kt = np.diag([alfa, 1.0])                                         # (x: horizontal, y: vertical)
    ke = k[tri].mean(axis=1)
    Ke = np.einsum("k,kia,ij,kjb->kab", A * ke, G, Kt, G)
    rows = np.repeat(tri, 3, axis=1).ravel()
    cols = np.tile(tri, (1, 3)).ravel()
    K = sp.csr_matrix((Ke.ravel(), (rows, cols)), shape=(n, n))
    x, y = xy[:, 0], xy[:, 1]
    H, tb = talude.H, np.tan(talude.beta)
    crista = (np.abs(y) < tol) & (x <= tol)
    face = (x > -tol) & (x < talude.xT + tol) & (np.abs(y - x * tb) < tol)
    pe = (np.abs(y - H) < tol) & (x >= talude.xT - tol)
    direita = (np.abs(x - talude.x_dir) < tol) & (y >= H - tol)
    esquerda = np.abs(x - talude.x_esq) < tol
    ud = np.full(n, np.nan)
    ud[crista] = 0.0
    ud[face] = np.where(y[face] <= hw, -GAMMA_W * y[face], -GAMMA_W * hw)
    ud[pe | direita] = -GAMMA_W * hw
    if lateral == "dirichlet":
        ud[esquerda & np.isnan(ud)] = 0.0
    D = ~np.isnan(ud)
    F = ~D
    u = np.zeros(n)
    u[D] = ud[D]
    u[F] = spl.spsolve(K[F][:, F].tocsc(), -K[F][:, D] @ u[D])
    return u


def campo_permanente(campo_base: CampoMalha, talude, **kw) -> CampoMalha:
    """Copia c, phi da exportacao e troca u pela solucao permanente."""
    u = resolver(campo_base.xy, campo_base.tri, talude, **kw)
    return CampoMalha(campo_base.xy, campo_base.tri, campo_base.cn, campo_base.phin, u,
                      campo_base.gamma_ef)
