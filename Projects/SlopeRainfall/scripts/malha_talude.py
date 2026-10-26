"""Malha de triângulos lineares (gmsh) do talude homogêneo de Rahardjo et al. (2007), Fig. 2.

Mesma geometria e mesmo refinamento de talude3f/malha_talude.py (o programa Python de referência), mas com
elementos de primeira ordem: o NeoPZ cria os espaços u (P2) e p_w, p_a (P1) sobre os triângulos.

Geometria: base em y = 0, pé em y = 2 Hs, crista em y = 3 Hs, platôs de 3 Hs na crista e no pé.
Grupos físicos das linhas (o programa C++ usa o id -tag):
  1 base   2 lateral esquerda   3 lateral direita   4 crista   5 face   6 terreno do pé
Saída: .msh 2.2 ASCII (triângulos tipo 2, linhas tipo 1), lida por ReadMsh (SlopeRainfall/RainfallModel.h).

Uso: python malha_talude.py Hs alpha [h_sup] [h_max_fator] [saida.msh]
     (h_sup = 0,75 m e h_max = 0,3 Hs por padrão: a malha 'media' de talude3f)
"""
import sys

import gmsh
import numpy as np


def geometria(hs=10.0, alpha_graus=45.0):
    xb = 3 * hs
    xc = xb + hs / np.tan(np.radians(alpha_graus))
    L = xc + 3 * hs
    return dict(a=(0.0, 3 * hs), b=(xb, 3 * hs), c=(xc, 2 * hs), d=(L, 2 * hs),
                f=(L, 0.0), g=(0.0, 0.0), L=L, xb=xb, xc=xc, hs=hs)


def _dist_seg(p, a, b):
    a, b, p = map(np.asarray, (a, b, p))
    t = np.clip(np.dot(p - a, b - a) / np.dot(b - a, b - a), 0, 1)
    return np.linalg.norm(p - (a + t * (b - a)))


def gerar(arq, hs=10.0, alpha=45.0, h_sup=0.75, h_max=None):
    """h_sup até 2 m de profundidade, crescendo até h_max a 2 m + Hs (a frente de umedecimento tem a mesma
    espessura em qualquer talude)"""
    G = geometria(hs, alpha)
    h_max = h_max or 0.3 * hs
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("talude")
    occ = gmsh.model.occ
    pts = [occ.addPoint(*G[k], 0.0) for k in "gfdcba"]
    lin = [occ.addLine(pts[i], pts[(i + 1) % 6]) for i in range(6)]
    occ.addPlaneSurface([occ.addCurveLoop(lin)])
    occ.synchronize()
    lados = {1: (G["g"], G["f"]), 3: (G["f"], G["d"]), 6: (G["d"], G["c"]), 5: (G["c"], G["b"]),
             4: (G["b"], G["a"]), 2: (G["a"], G["g"])}
    grupos = {}
    for _, c in gmsh.model.getEntities(1):
        a0, a1 = gmsh.model.getParametrizationBounds(1, c)
        pm = gmsh.model.getValue(1, c, [0.5 * (a0[0] + a1[0])])[:2]
        tag = [k for k, (a, b) in lados.items() if _dist_seg(pm, a, b) < 1e-6 * hs][0]
        grupos.setdefault(tag, []).append(c)
    for tag, cs in grupos.items():
        gmsh.model.addPhysicalGroup(1, cs, tag)
    gmsh.model.addPhysicalGroup(2, [s for _, s in gmsh.model.getEntities(2)], 1)
    f = gmsh.model.mesh.field
    d = f.add("Distance")
    f.setNumbers(d, "CurvesList", grupos[4] + grupos[5] + grupos[6])
    f.setNumber(d, "Sampling", 400)
    t = f.add("Threshold")
    f.setNumber(t, "InField", d)
    f.setNumber(t, "SizeMin", h_sup)
    f.setNumber(t, "SizeMax", h_max)
    f.setNumber(t, "DistMin", 2.0)
    f.setNumber(t, "DistMax", 2.0 + hs)
    f.setAsBackgroundMesh(t)
    for o, v in (("Mesh.MeshSizeExtendFromBoundary", 0), ("Mesh.MeshSizeFromPoints", 0),
                 ("Mesh.MeshSizeFromCurvature", 0), ("Mesh.Algorithm", 6), ("Mesh.ElementOrder", 1),
                 ("Mesh.MshFileVersion", 2.2), ("Mesh.Binary", 0)):
        gmsh.option.setNumber(o, v)
    gmsh.model.mesh.generate(2)
    gmsh.write(arq)
    gmsh.finalize()
    return arq


if __name__ == "__main__":
    hs = float(sys.argv[1]) if len(sys.argv) > 1 else 10.0
    alpha = float(sys.argv[2]) if len(sys.argv) > 2 else 45.0
    h_sup = float(sys.argv[3]) if len(sys.argv) > 3 else 0.75
    fmax = float(sys.argv[4]) if len(sys.argv) > 4 else 0.3
    arq = sys.argv[5] if len(sys.argv) > 5 else f"talude_H{hs:g}_a{alpha:g}.msh"
    gerar(arq, hs, alpha, h_sup, fmax * hs)
    txt = open(arq).read()
    print(arq)
