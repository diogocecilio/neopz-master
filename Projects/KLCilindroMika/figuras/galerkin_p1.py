"""galerkin_p1.py: uma KL de Galerkin pequena em numpy/scipy, com elementos P1 numa malha de triangulos, so para mostrar
que, na KL numerica, os modos (e o campo para os mesmos xi) mudam com a malha.

    C v = lambda B v,  B_ij = int N_i N_j,  C_ij = int int N_i(x) C(x, x') N_j(x') dx dx'

com a regra de Dunavant de 6 pontos (grau 4) em cada triangulo nas duas integrais (o produto pelo nucleo e feito em
blocos de pontos, sem montar a matriz do nucleo inteira). Modos B-ortonormais (int phi^2 = 1)."""
import numpy as np
import scipy.linalg as sla
import scipy.sparse as sp

# Dunavant, grau 4: coordenadas baricentricas e pesos (somam 1)
_A, _WA = 0.445948490915965, 0.223381589678011
_B, _WB = 0.091576213509771, 0.109951743655322
BARY = np.array([[_A, _A, 1 - 2 * _A], [_A, 1 - 2 * _A, _A], [1 - 2 * _A, _A, _A],
                 [_B, _B, 1 - 2 * _B], [_B, 1 - 2 * _B, _B], [1 - 2 * _B, _B, _B]])
PESOS = np.array([_WA] * 3 + [_WB] * 3)


def le_malha(arquivo):
    """nos (N x 2) e triangulos (nt x 3, indices 0..N-1) de uma malha do Gmsh (so os triangulos)"""
    import gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try:
        gmsh.open(arquivo)
        tags, coords, _ = gmsh.model.mesh.getNodes()
        _, nos = gmsh.model.mesh.getElementsByType(2)
    finally:
        gmsh.finalize()
    X = coords.reshape(-1, 3)[:, :2]
    ind = np.full(int(tags.max()) + 1, -1)
    ind[tags.astype(int)] = np.arange(len(tags))
    T = ind[nos.astype(int)].reshape(-1, 3)
    usados = np.unique(T)
    novo = np.full(len(X), -1)
    novo[usados] = np.arange(len(usados))
    return X[usados], novo[T]


def quadratura(X, T):
    """pontos (nq x 2), pesos x area (nq) e a matriz das funcoes de forma nos pontos (nq x N, esparsa)"""
    P = X[T]  # (nt, 3, 2)
    area = 0.5 * np.abs((P[:, 1, 0] - P[:, 0, 0]) * (P[:, 2, 1] - P[:, 0, 1]) -
                        (P[:, 2, 0] - P[:, 0, 0]) * (P[:, 1, 1] - P[:, 0, 1]))
    nt, nr = len(T), len(PESOS)
    xq = np.einsum("rk,tkd->trd", BARY, P).reshape(-1, 2)
    wq = (area[:, None] * PESOS[None, :]).ravel()
    lin = np.repeat(np.arange(nt * nr), 3)
    col = np.repeat(T, nr, axis=0).ravel()
    val = np.tile(BARY, (nt, 1)).ravel()
    Phi = sp.csr_matrix((val, (lin, col)), shape=(nt * nr, len(X)))
    return xq, wq, Phi


def kl_p1(X, T, M, Lx=20., Ly=2., bloco=1500):
    """os M maiores autopares (lambda, V) da KL P1; V (N x M) B-ortonormal"""
    xq, wq, Phi = quadratura(X, T)
    WPhi = sp.diags(wq) @ Phi
    B = (Phi.T @ WPhi).toarray()
    N = len(X)
    C = np.zeros((N, N))
    WPhiT = WPhi.T.tocsr()
    for a in range(0, len(wq), bloco):
        b = min(len(wq), a + bloco)
        K = np.exp(-np.abs(xq[a:b, 0][:, None] - xq[None, :, 0]) / Lx - np.abs(xq[a:b, 1][:, None] - xq[None, :, 1]) / Ly)
        G = (WPhiT @ K.T).T  # (bloco x N) = K @ WPhi
        C += WPhi[a:b].T @ G
    C = 0.5 * (C + C.T)
    lam, V = sla.eigh(C, B, subset_by_index=[N - M, N - 1])
    o = np.argsort(-lam)
    lam, V = lam[o], V[:, o]
    # sinal: o maior coeficiente em modulo positivo (arbitrario, como em qualquer solver)
    for k in range(M):
        i = np.argmax(np.abs(V[:, k]))
        if V[i, k] < 0:
            V[:, k] *= -1
    return lam, V, B


def avalia_p1(X, T, valores_nos, pontos):
    """interpolacao P1 (exata para campos P1) dos valores nodais em pontos quaisquer; valores_nos (N,) ou (N, m)"""
    from matplotlib.tri import Triangulation, TriFinder  # noqa: F401
    tri = Triangulation(X[:, 0], X[:, 1], T)
    finder = tri.get_trifinder()
    t = finder(pontos[:, 0], pontos[:, 1])
    if np.any(t < 0):
        raise ValueError("%d pontos fora da malha" % np.sum(t < 0))
    P = X[T[t]]
    # coordenadas baricentricas
    d = (P[:, 1, 1] - P[:, 2, 1]) * (P[:, 0, 0] - P[:, 2, 0]) + (P[:, 2, 0] - P[:, 1, 0]) * (P[:, 0, 1] - P[:, 2, 1])
    l0 = ((P[:, 1, 1] - P[:, 2, 1]) * (pontos[:, 0] - P[:, 2, 0]) + (P[:, 2, 0] - P[:, 1, 0]) * (pontos[:, 1] - P[:, 2, 1])) / d
    l1 = ((P[:, 2, 1] - P[:, 0, 1]) * (pontos[:, 0] - P[:, 2, 0]) + (P[:, 0, 0] - P[:, 2, 0]) * (pontos[:, 1] - P[:, 2, 1])) / d
    l2 = 1 - l0 - l1
    v = np.asarray(valores_nos)
    nos = T[t]
    if v.ndim == 1:
        return l0 * v[nos[:, 0]] + l1 * v[nos[:, 1]] + l2 * v[nos[:, 2]]
    return l0[:, None] * v[nos[:, 0]] + l1[:, None] * v[nos[:, 1]] + l2[:, None] * v[nos[:, 2]]
