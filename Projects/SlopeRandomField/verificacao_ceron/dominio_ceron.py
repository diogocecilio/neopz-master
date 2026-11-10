"""Efeito do dominio hidraulico no Gamma do metodo do Ceron (espiral_ceron.py), no rebaixamento total (nivel no pe):
gera a malha do talude (crista em (10, 10), pe em (15, 5)) com o gmsh, com <atras> m de crista antes da face, <frente> m
de terreno depois do pe e a base <abaixo> m abaixo do pe; resolve a percolacao permanente saturada P1 (a mesma do
material com SATURADO=1 SUPERFICIE_P0=1 KW=1e30, verificada por laplace_ceron.py), com p = 0 na superficie e:
  contorno 'carga': lateral esquerda com a carga inicial (h = 10), direita com a do reservatorio (h = 5), base sem fluxo
                    (o do SlopeRandomField);
  contorno 'fluxo0': as duas laterais e a base sem fluxo (a zona de fluxo do Ceron, com v = 0 alem de R_e);
e calcula o Gamma pela espiral logaritmica com essa pressao.
  contorno 'esq,base,dir' (cada um f, h10 ou h5): por exemplo h10,h10,h5 (u = 0 na esquerda e na base).
Uso: python3 -I dominio_ceron.py <atras> <abaixo> <frente> <carga|fluxo0|esq,base,dir> [hFoco=0.25] [c=10] [phi=30]
     [gamma=20]"""
import importlib.util
import os
import sys

import gmsh
import matplotlib.tri as mtri
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

spec = importlib.util.spec_from_file_location(
    'espiral_ceron', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'espiral_ceron.py'))
esp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(esp)

atras, abaixo, frente = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
contorno = sys.argv[4]
spec_cont = {'carga': 'h10,f,h5', 'fluxo0': 'f,f,f'}.get(contorno, contorno).split(',')
assert len(spec_cont) == 3 and all(p in ('f', 'h10', 'h5') for p in spec_cont), contorno
hFoco = float(sys.argv[5]) if len(sys.argv) > 5 else 0.25
c = float(sys.argv[6]) if len(sys.argv) > 6 else 10.
phi = np.radians(float(sys.argv[7]) if len(sys.argv) > 7 else 30.)
g = float(sys.argv[8]) if len(sys.argv) > 8 else 20.
gw, nivel, hl = 9.81, 5., 10.
xl, xr, yb = 10. - atras, 15. + frente, 5. - abaixo
esp.XMAX = xr

gmsh.initialize()
gmsh.option.setNumber('General.Terminal', 0)
gmsh.model.add('t')
pts = [(xl, yb), (xr, yb), (xr, 5.), (15., 5.), (10., 10.), (xl, 10.)]
tags = [gmsh.model.geo.addPoint(x, y, 0.) for x, y in pts]
lin = [gmsh.model.geo.addLine(tags[i], tags[(i + 1) % 6]) for i in range(6)]
loop = gmsh.model.geo.addCurveLoop(lin)
gmsh.model.geo.addPlaneSurface([loop])
gmsh.model.geo.synchronize()
hMax = max(1.0, min(6.0, 0.12 * max(atras, abaixo, frente)))
f1 = gmsh.model.mesh.field.add('Box')
for k, v in (('VIn', hFoco), ('VOut', hMax), ('XMin', 7.4), ('XMax', 17.5), ('YMin', 3.5), ('YMax', 10.1),
             ('Thickness', 3.)):
    gmsh.model.mesh.field.setNumber(f1, k, v)
f2 = gmsh.model.mesh.field.add('Distance')
gmsh.model.mesh.field.setNumbers(f2, 'CurvesList', [lin[2], lin[3], lin[4]])
gmsh.model.mesh.field.setNumber(f2, 'Sampling', 2000)
f3 = gmsh.model.mesh.field.add('Threshold')
for k, v in (('InField', f2), ('SizeMin', 0.6), ('SizeMax', hMax), ('DistMin', 1.), ('DistMax', 25.)):
    gmsh.model.mesh.field.setNumber(f3, k, v)
f4 = gmsh.model.mesh.field.add('Min')
gmsh.model.mesh.field.setNumbers(f4, 'FieldsList', [f1, f3])
gmsh.model.mesh.field.setAsBackgroundMesh(f4)
for k, v in (('Mesh.MeshSizeExtendFromBoundary', 0), ('Mesh.MeshSizeFromPoints', 0), ('Mesh.MeshSizeFromCurvature', 0),
             ('Mesh.Algorithm', 6)):
    gmsh.option.setNumber(k, v)
gmsh.model.mesh.generate(2)
ntags, coords, _ = gmsh.model.mesh.getNodes()
X = coords[0::3].copy(); Y = coords[1::3].copy()
ind = {int(t): i for i, t in enumerate(ntags)}
et, _, enodes = gmsh.model.mesh.getElements(2)
T = np.array([ind[int(n)] for n in enodes[0]]).reshape(-1, 3)
gmsh.finalize()
N = len(X)

lin_, col_, val_ = [], [], []
for t in T:
    x, y = X[t], Y[t]
    a2 = (x[1] - x[0]) * (y[2] - y[0]) - (x[2] - x[0]) * (y[1] - y[0])
    b = np.array([y[1] - y[2], y[2] - y[0], y[0] - y[1]])
    cc = np.array([x[2] - x[1], x[0] - x[2], x[1] - x[0]])
    K = (np.outer(b, b) + np.outer(cc, cc)) / (2. * abs(a2))
    for i in range(3):
        for j in range(3):
            lin_.append(t[i]); col_.append(t[j]); val_.append(K[i, j])
A = coo_matrix((val_, (lin_, col_)), shape=(N, N)).tocsr()
tol = 1e-7
fix = {}
for i in range(N):
    x, y = X[i], Y[i]
    na_sup = (abs(y - 10.) < tol and x <= 10. + tol) or (abs(y - 5.) < tol and x >= 15. - tol) or \
             (10. - tol <= x <= 15. + tol and abs((y - 10.) + (x - 10.)) < tol)
    if na_sup:
        fix[i] = max(y, nivel)  # h = y + p/gw, p = gw max(nivel - y, 0)
        continue
    # contorno 'esq,base,dir' com cada parte em f (sem fluxo), h10 (carga inicial, u = 0) ou h5 (carga do reservatorio)
    for parte, cond in zip(spec_cont, (abs(x - xl) < tol, abs(y - yb) < tol, abs(x - xr) < tol)):
        if cond and parte != 'f':
            fix[i] = hl if parte == 'h10' else nivel
livres = np.array([i for i in range(N) if i not in fix])
fixos = np.array(sorted(fix))
hfix = np.array([fix[i] for i in fixos])
h = np.zeros(N); h[fixos] = hfix
h[livres] = spsolve(A[livres][:, livres].tocsc(), -A[livres][:, fixos] @ hfix)
P = gw * (h - Y)


class PressaoMalha:
    def __init__(self):
        self.f = mtri.LinearTriInterpolator(mtri.Triangulation(X, Y, T), P)

    def __call__(self, x, y):
        return np.ma.filled(self.f(x, y), np.nan)


(G, v), nval = esp.busca(c, phi, g, PressaoMalha())
Bx, By = esp.ponto_B(min(max(v[2], -1.), 0.99))
k = int(np.argmin((X - 10.) ** 2 + (Y - 5.) ** 2))
print('dominio atras %g, abaixo do pe %g, frente %g, contorno %s: %d triangulos, p(%.1f, %.1f) = %.2f kPa, p min %.2f; '
      'Gamma = %.4f; O = (%.3f, %.3f), B = (%.3f, %.3f)'
      % (atras, abaixo, frente, contorno, len(T), X[k], Y[k], P[k], P.min(), G, v[0], v[1], Bx, By))
