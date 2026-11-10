"""Verificacao independente da percolacao permanente saturada do exemplo do Ceron (Eqs. 19-21 de Vargas Ceron et al.
2025) na malha de uma rodada do SlopeRandomField: elementos P1 nos mesmos triangulos (triangulos.csv), carga total
h = y + p/gamma_w com div(k grad h) = 0 e k constante, e o contorno da rodada:
  superficie (superficie=1 em pw.csv): p = gamma_w max(nivel - y, 0)  (p = 0 acima do nivel, SUPERFICIE_P0=1);
  lateral esquerda (x = xl), y <= hl: p = gamma_w (hl - y);
  lateral direita (x = xr), y <= nivel: p = gamma_w (nivel - y);
  base: fluxo nulo;
  com o ultimo argumento 'fluxo0' (LATERAIS_SEM_FLUXO=1), as laterais tambem sem fluxo.
Compara p nos nos com o pw_kPa de pw.csv no instante pedido.
Uso: python3 -I laplace_ceron.py <pasta da rodada> [t_h=0.01] [nivel=5] [hl=10] [xl=0] [xr=25] [fluxo0]"""
import csv
import sys

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

pasta = sys.argv[1]
t_h = float(sys.argv[2]) if len(sys.argv) > 2 else 0.01
nivel = float(sys.argv[3]) if len(sys.argv) > 3 else 5.
hl = float(sys.argv[4]) if len(sys.argv) > 4 else 10.
xl = float(sys.argv[5]) if len(sys.argv) > 5 else 0.
xr = float(sys.argv[6]) if len(sys.argv) > 6 else 25.
lateraisSemFluxo = len(sys.argv) > 7 and sys.argv[7] == 'fluxo0'
gw = 9.81

nos = {}
for r in csv.DictReader(open(pasta + '/pw.csv')):
    if abs(float(r['t_h']) - t_h) > 1e-9:
        continue
    nos[int(r['no'])] = (float(r['x']), float(r['y']), float(r['pw_kPa']), int(r['superficie']))
ids = sorted(nos)
ind = {n: i for i, n in enumerate(ids)}
N = len(ids)
X = np.array([nos[n][0] for n in ids])
Y = np.array([nos[n][1] for n in ids])
Pfe = np.array([nos[n][2] for n in ids])
sup = np.array([nos[n][3] for n in ids], dtype=bool)

tri = np.array([[ind[int(r['n0'])], ind[int(r['n1'])], ind[int(r['n2'])]]
                for r in csv.DictReader(open(pasta + '/triangulos.csv'))])
lin, col, val = [], [], []
for t in tri:
    x, y = X[t], Y[t]
    area2 = (x[1] - x[0]) * (y[2] - y[0]) - (x[2] - x[0]) * (y[1] - y[0])
    b = np.array([y[1] - y[2], y[2] - y[0], y[0] - y[1]])
    c = np.array([x[2] - x[1], x[0] - x[2], x[1] - x[0]])
    K = (np.outer(b, b) + np.outer(c, c)) / (2. * abs(area2))
    for i in range(3):
        for j in range(3):
            lin.append(t[i]); col.append(t[j]); val.append(K[i, j])
A = coo_matrix((val, (lin, col)), shape=(N, N)).tocsr()

tol = 1e-9
fix = {}
for i in range(N):
    if sup[i]:
        fix[i] = Y[i] + max(nivel - Y[i], 0.)        # h = y + p/gw
    if lateraisSemFluxo:
        continue
    if abs(X[i] - xl) < tol and Y[i] <= hl + tol:
        fix[i] = hl
    if abs(X[i] - xr) < tol and Y[i] <= nivel + tol:
        fix[i] = nivel
livres = np.array([i for i in range(N) if i not in fix])
fixos = np.array(sorted(fix))
hfix = np.array([fix[i] for i in fixos])
h = np.zeros(N)
h[fixos] = hfix
Aff = A[livres][:, livres]
rhs = -A[livres][:, fixos] @ hfix
h[livres] = spsolve(Aff.tocsc(), rhs)
P = gw * (h - Y)

d = P - Pfe
print('nos %d (fixos %d, livres %d), triangulos %d' % (N, len(fixos), len(livres), len(tri)))
print('p FE: min %.6g max %.6g; p Laplace: min %.6g max %.6g' % (Pfe.min(), Pfe.max(), P.min(), P.max()))
print('diferenca max |p_Laplace - p_FE| = %.3e kPa (relativa ao max: %.3e); media |dif| = %.3e kPa'
      % (np.abs(d).max(), np.abs(d).max() / np.abs(P).max(), np.abs(d).mean()))
i = int(np.argmax(np.abs(d)))
print('maior diferenca no no %d em (%.3f, %.3f): Laplace %.6g, FE %.6g' % (ids[i], X[i], Y[i], P[i], Pfe[i]))
for (xq, yq) in [(5., 5.), (10., 5.), (12.5, 7.5), (15., 2.5), (20., 2.5), (10., 0.)]:
    k = int(np.argmin((X - xq) ** 2 + (Y - yq) ** 2))
    print('  no mais perto de (%.1f, %.1f): (%.3f, %.3f) p Laplace %.4f p FE %.4f' % (xq, yq, X[k], Y[k], P[k], Pfe[k]))
