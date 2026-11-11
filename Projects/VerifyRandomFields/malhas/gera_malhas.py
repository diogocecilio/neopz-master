#!/usr/bin/env python3
"""Gera as malhas uniformes do talude (talude_uniforme_h025.msh e talude_uniforme_h05.msh) a partir dos .geo desta
pasta com o modulo gmsh do Python (4.15), numa thread so (malha reprodutivel). Uso: python3 gera_malhas.py [arquivo.geo ...]
"""
import os
import sys

import gmsh

AQUI = os.path.dirname(os.path.abspath(__file__))


def gera(geo):
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.option.setNumber("General.NumThreads", 1)
        gmsh.option.setNumber("Mesh.MaxNumThreads2D", 1)
        gmsh.open(geo)
        gmsh.model.mesh.generate(2)
        msh = os.path.splitext(geo)[0] + ".msh"
        gmsh.write(msh)
        _, coords, _ = gmsh.model.mesh.getNodes()
        tipos, elems, _ = gmsh.model.mesh.getElements(2)
        ntri = sum(len(e) for e in elems)
        print(f"{os.path.basename(msh)}: {len(coords) // 3} nos, {ntri} triangulos (gmsh {gmsh.__version__})")
    finally:
        gmsh.finalize()


if __name__ == "__main__":
    arquivos = sys.argv[1:] or [os.path.join(AQUI, f) for f in ("talude_uniforme_h025.geo", "talude_uniforme_h05.geo")]
    for g in arquivos:
        gera(g)
