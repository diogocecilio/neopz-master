"""Malha de triângulos lineares (gmsh) de um talude de superfície poligonal, com camadas paralelas à superfície.

Domínio: a superfície (polilinha da esquerda para a direita), as laterais verticais e a base em y = ybase.
Camadas (opcional): faixas de profundidade vertical [0, d1], [d1, d2], ... abaixo de uma superfície de referência
(--ref; padrão: a própria superfície), recortadas pela superfície do terreno (um corte expõe as camadas profundas);
o material id é 1 na primeira camada, 2 na segunda, ... e N + 1 abaixo da última (sem camadas: id 1).
Grupos físicos das linhas (o programa C++ usa o id -tag): 1 base, 2 lateral esquerda, 3 lateral direita,
5 superfície.

Exemplos:
  python malha_poligono.py --sup 0:40,20:40,45.17:10,70:10 --hsup 0.5 --hmax 3 --saida zhang_A.msh
  python malha_poligono.py --sup 0:28,8:28,9.377:32,34:42,69:42 --ref 0:28.21,34:42,69:42 \\
         --camadas 0.5,1,1.5 --hsup 0.25 --dsup 2 --hmax 4 --saida pradhan_corte.msh
"""
import argparse

import gmsh
import numpy as np


def pontos(s):
    return [tuple(map(float, p.split(":"))) for p in s.split(",")]


def cota(poli, x):
    xs, ys = zip(*poli)
    return float(np.interp(x, xs, ys))


def gerar(sup, saida, ybase=0.0, hsup=0.5, hmax=3.0, dsup=2.0, dcres=None, camadas=(), ref=None):
    ref = ref or sup
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("talude")
    occ = gmsh.model.occ
    x0, x1 = sup[0][0], sup[-1][0]
    contorno = list(sup) + [(x1, ybase), (x0, ybase)]
    pts = [occ.addPoint(x, y, 0.0) for x, y in contorno]
    lin = [occ.addLine(pts[i], pts[(i + 1) % len(pts)]) for i in range(len(pts))]
    dom = occ.addPlaneSurface([occ.addCurveLoop(lin)])
    cortes = []
    for d in camadas:
        # offset of the reference surface, extended beyond the lateral boundaries
        xs = [p[0] for p in ref]
        xs = [x0 - 1.0] + [x for x in xs if x0 < x < x1] + [x1 + 1.0]
        q = [occ.addPoint(x, cota(ref, x) - d, 0.0) for x in xs]
        cortes += [occ.addLine(q[i], q[i + 1]) for i in range(len(q) - 1)]
    if cortes:
        occ.fragment([(2, dom)], [(1, c) for c in cortes])
    occ.synchronize()
    # remove the pieces of the offset lines outside the domain
    for dim, tag in gmsh.model.getEntities(1):
        up, _ = gmsh.model.getAdjacencies(1, tag)
        if len(up) == 0:
            gmsh.model.removeEntities([(1, tag)], recursive=True)
    # layers by the depth of the centre of each surface below the reference surface
    grupos2 = {}
    for _, s in gmsh.model.getEntities(2):
        xc, yc, _ = occ.getCenterOfMass(2, s)
        prof = cota(ref, xc) - yc
        lay = 1 + sum(1 for d in camadas if prof > d)
        grupos2.setdefault(lay, []).append(s)
    for lay, ss in grupos2.items():
        gmsh.model.addPhysicalGroup(2, ss, lay)
    # boundary curves (one adjacent surface): base, laterals, surface
    grupos1 = {1: [], 2: [], 3: [], 5: []}
    tol = 1e-6 * (x1 - x0)
    for _, c in gmsh.model.getEntities(1):
        up, _ = gmsh.model.getAdjacencies(1, c)
        if len(up) != 1:
            continue
        a0, a1 = gmsh.model.getParametrizationBounds(1, c)
        xm, ym = gmsh.model.getValue(1, c, [0.5 * (a0[0] + a1[0])])[:2]
        if abs(ym - ybase) < tol:
            grupos1[1].append(c)
        elif abs(xm - x0) < tol:
            grupos1[2].append(c)
        elif abs(xm - x1) < tol:
            grupos1[3].append(c)
        else:
            grupos1[5].append(c)
    for tag, cs in grupos1.items():
        if cs:
            gmsh.model.addPhysicalGroup(1, cs, tag)
    f = gmsh.model.mesh.field
    d = f.add("Distance")
    f.setNumbers(d, "CurvesList", grupos1[5])
    f.setNumber(d, "Sampling", 600)
    t = f.add("Threshold")
    f.setNumber(t, "InField", d)
    f.setNumber(t, "SizeMin", hsup)
    f.setNumber(t, "SizeMax", hmax)
    f.setNumber(t, "DistMin", dsup)
    f.setNumber(t, "DistMax", dsup + (dcres if dcres else 4 * hmax))
    f.setAsBackgroundMesh(t)
    for o, v in (("Mesh.MeshSizeExtendFromBoundary", 0), ("Mesh.MeshSizeFromPoints", 0),
                 ("Mesh.MeshSizeFromCurvature", 0), ("Mesh.Algorithm", 6), ("Mesh.ElementOrder", 1),
                 ("Mesh.MshFileVersion", 2.2), ("Mesh.Binary", 0)):
        gmsh.option.setNumber(o, v)
    gmsh.model.mesh.generate(2)
    gmsh.write(saida)
    gmsh.finalize()
    return retag(saida, ref, camadas)


def retag(saida, ref, camadas):
    """Material of each triangle by the depth of its centroid below the reference surface (the pieces of a layer
    are not convex, so their centres of mass do not identify them)"""
    linhas = open(saida).read().split("\n")
    i = linhas.index("$Nodes")
    n = int(linhas[i + 1])
    xy = {}
    for l in linhas[i + 2:i + 2 + n]:
        v = l.split()
        xy[int(v[0])] = (float(v[1]), float(v[2]))
    j = linhas.index("$Elements")
    m = int(linhas[j + 1])
    usados, ntri = set(), 0
    for k in range(j + 2, j + 2 + m):
        v = linhas[k].split()
        if v[1] != "2":
            continue
        nt = int(v[2])
        nos = [int(t) for t in v[3 + nt:]]
        xc = sum(xy[q][0] for q in nos) / 3
        yc = sum(xy[q][1] for q in nos) / 3
        lay = 1 + sum(1 for d in camadas if cota(ref, xc) - yc > d)
        v[3] = str(lay)
        linhas[k] = " ".join(v)
        usados.add(lay)
        ntri += 1
    open(saida, "w").write("\n".join(linhas))
    return ntri, sorted(usados)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sup", required=True, help="superfície x:y,x:y,... da esquerda para a direita")
    p.add_argument("--ref", help="superfície de referência das camadas (padrão: --sup)")
    p.add_argument("--camadas", default="", help="profundidades verticais das bases das camadas, p. ex. 0.5,1,1.5")
    p.add_argument("--ybase", type=float, default=0.0)
    p.add_argument("--hsup", type=float, default=0.5, help="tamanho dos elementos até dsup de profundidade")
    p.add_argument("--dsup", type=float, default=2.0)
    p.add_argument("--hmax", type=float, default=3.0)
    p.add_argument("--dcres", type=float, help="distância de transição de hsup a hmax (padrão 4 hmax)")
    p.add_argument("--saida", required=True)
    a = p.parse_args()
    cam = [float(v) for v in a.camadas.split(",") if v]
    n, lays = gerar(pontos(a.sup), a.saida, a.ybase, a.hsup, a.hmax, a.dsup, a.dcres, cam,
                    pontos(a.ref) if a.ref else None)
    print(f"{a.saida}: {n} triângulos, materiais {lays}")
