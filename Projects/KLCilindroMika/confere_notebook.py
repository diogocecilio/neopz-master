# confere_notebook.py - conferencia independente (numpy/scipy) da KL de Galerkin na malha do notebook randomfied-v2.nb:
# meio cilindro 8 <= r <= 10, 0 <= theta <= pi, 0 <= z <= 15, nucleo exp(-|x - x'| / 5), hexaedros retos com os nos em
# (r cos theta, r sin theta, z) (mapa trilinear). Espaco q2 (27 funcoes por hexaedro, o do NeoPZ com P=2) ou s20
# (serendipity de 20 nos, o hexaedro de segunda ordem do NDSolve`FEM`). B com Gauss 6 x 6 x 6 (exata nesta geometria),
# C com Gauss n x n x n em cada elemento de cada par. Escreve os 20 maiores autovalores e a soma de todos.
#   python3 confere_notebook.py q2 5            (malha 15x1x3 do notebook; o NeoPZ da o mesmo com GEOMETRIA=reta NG=5)
#   python3 confere_notebook.py s20 5 15x1x3    (malha nt x nr x nz)
import sys, itertools
import numpy as np, scipy.linalg as sla

espaco = sys.argv[1]
ngC = int(sys.argv[2])
nt, nr, nz = (int(v) for v in (sys.argv[3].split('x') if len(sys.argv) > 3 else ('15', '1', '3')))
ri, re, H, l = 8., 10., 15., 5.
th = np.linspace(0., np.pi, nt + 1); rr = np.linspace(ri, re, nr + 1); zz = np.linspace(0., H, nz + 1)

def X(i, j, k):
    return np.array([rr[j] * np.cos(th[i]), rr[j] * np.sin(th[i]), zz[k]])

cantos = [(-1, -1, -1), (1, -1, -1), (1, 1, -1), (-1, 1, -1), (-1, -1, 1), (1, -1, 1), (1, 1, 1), (-1, 1, 1)]
# nos locais: (a, b, c) em {-1, 0, 1}^3; o serendipity nao tem os centros das faces nem o do elemento
locais = [t for t in itertools.product([-1, 0, 1], repeat=3) if espaco == 'q2' or sum(s == 0 for s in t) <= 1]
L = {-1: lambda s: 0.5 * s * (s - 1), 0: lambda s: 1 - s * s, 1: lambda s: 0.5 * s * (s + 1)}

def forma(t, xi):
    a, b, c = t
    if espaco == 'q2':  # Lagrange quadratico tensorial
        return L[a](xi[:, 0]) * L[b](xi[:, 1]) * L[c](xi[:, 2])
    if 0 not in t:  # canto do serendipity
        return 0.125 * (1 + a * xi[:, 0]) * (1 + b * xi[:, 1]) * (1 + c * xi[:, 2]) * (a * xi[:, 0] + b * xi[:, 1] + c * xi[:, 2] - 2)
    f = np.ones(len(xi))  # meio de aresta
    for d, s in enumerate(t):
        f = f * ((1 - xi[:, d] ** 2) if s == 0 else 0.5 * (1 + s * xi[:, d]))
    return f

def regra(n):
    x, w = np.polynomial.legendre.leggauss(n)
    return np.array(list(itertools.product(x, x, x))), np.array([a * b * c for a, b, c in itertools.product(w, w, w)])

def mapa(cel, xi):  # mapa trilinear dos 8 nos e |det J|
    i, j, k = cel
    Xc = np.array([X(i + (a + 1) // 2, j + (b + 1) // 2, k + (c + 1) // 2) for a, b, c in cantos])
    N = np.array([0.125 * (1 + a * xi[:, 0]) * (1 + b * xi[:, 1]) * (1 + c * xi[:, 2]) for a, b, c in cantos]).T
    dN = np.zeros((len(xi), 3, 8))
    for m, (a, b, c) in enumerate(cantos):
        dN[:, 0, m] = 0.125 * a * (1 + b * xi[:, 1]) * (1 + c * xi[:, 2])
        dN[:, 1, m] = 0.125 * b * (1 + a * xi[:, 0]) * (1 + c * xi[:, 2])
        dN[:, 2, m] = 0.125 * c * (1 + a * xi[:, 0]) * (1 + b * xi[:, 1])
    return N @ Xc, np.abs(np.linalg.det(dN @ Xc))

# numeracao global: o no local (a, b, c) da celula (i, j, k) e o ponto (2i+1+a, 2j+1+b, 2k+1+c) da grade de passo 1/2
ids, celulas, glob = {}, [(i, j, k) for k in range(nz) for j in range(nr) for i in range(nt)], []
for (i, j, k) in celulas:
    glob.append([ids.setdefault((2 * i + 1 + t[0], 2 * j + 1 + t[1], 2 * k + 1 + t[2]), len(ids)) for t in locais])
neq = len(ids)

def monta(n):  # pontos, pesos x |det J| e matriz (pontos x equacoes) das funcoes de forma
    P, W = regra(n)
    Phi = np.array([forma(t, P) for t in locais]).T
    Xs, Ws, S = [], [], np.zeros((len(celulas) * len(P), neq))
    for e, cel in enumerate(celulas):
        x, dj = mapa(cel, P)
        Xs.append(x); Ws.append(W * dj)
        S[e * len(P):(e + 1) * len(P), glob[e]] += Phi
    return np.vstack(Xs), np.concatenate(Ws), S

xB, wB, SB = monta(6)
B = SB.T @ (wB[:, None] * SB)
xC, wC, SC = monta(ngC)
K = np.exp(-np.sqrt(((xC[:, None, :] - xC[None, :, :]) ** 2).sum(-1)) / l)
SW = wC[:, None] * SC
lam = sla.eigh(SW.T @ K @ SW, B, eigvals_only=True)[::-1]
print("espaco %s, malha %dx%dx%d, %d equacoes, C com %d^3 pontos por elemento, volume %.10f" % (espaco, nt, nr, nz, neq, ngC, wB.sum()))
print(" ".join("%.7f" % v for v in lam[:20]))
print("soma de todos os autovalores %.7f" % lam.sum())
