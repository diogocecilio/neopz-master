#!/usr/bin/env python3
"""refina_banda.py <rodada grossa> <malha grossa .msh ou .geo> <saida sem extensao> [h] [folga] [limiar]

A segunda etapa do refinamento por realizacao (REFINA_BANDA=1 em rodar_mc.sh): uma malha conforme refinada so em volta
da banda de ruptura que a rodada na malha grossa encontrou, gerada pelo gmsh a partir do .geo da malha grossa.

- A banda: os pontos de integracao com sqrt(J2(eps_p)) >= limiar x maximo (padrao 0,1, o BANDA_LIMIAR do programa) no
  instante do fator (o menor FS de fs.csv), lidos de <rodada grossa>/mecanismo.csv e arredondados numa grade de h/1 m.
- O .geo novo e o da malha grossa (mesma geometria, mesmos grupos fisicos, mesmos campos de tamanho) com mais um campo:
  elementos de tamanho h (padrao 0,125 m) ate 'folga' (padrao 0,375 m) dos pontos da banda, crescendo ate o hMax da malha
  grossa a folga + 3 m; o campo de fundo passa a ser o minimo do antigo e deste. Assim o contorno, e com ele o retangulo
  do Karhunen-Loeve, nao muda, e a realizacao tem o mesmo campo aleatorio nas duas malhas.
- O gmsh: o executavel gmsh, se estiver no PATH; senao o modulo python gmsh (python3 -m pip install --user gmsh).

Escreve <saida>.geo e <saida>.msh e imprime o numero de triangulos."""

import csv
import math
import os
import re
import shutil
import subprocess
import sys


def banda(rodada, limiar):
    fs = os.path.join(rodada, "fs.csv")
    t_fs = None
    if os.path.isfile(fs):
        linhas = [r for r in csv.DictReader(open(fs)) if r.get("FS_SRM")]
        if linhas:
            t_fs = float(min(linhas, key=lambda r: float(r["FS_SRM"]))["t_h"])
    pts = []
    for r in csv.DictReader(open(os.path.join(rodada, "mecanismo.csv"))):
        if t_fs is None or abs(float(r["t_h"]) - t_fs) < 1e-9:
            pts.append((float(r["x"]), float(r["y"]), float(r["epsp_raizJ2"])))
    if not pts:
        raise SystemExit("refina_banda: mecanismo.csv sem pontos no instante do fator em " + rodada)
    mx = max(p[2] for p in pts)
    if not mx > 0:
        raise SystemExit("refina_banda: sem deformacao plastica em " + rodada)
    return [(x, y) for x, y, e in pts if e >= limiar * mx]


def geo_refinado(geo_grosso, pontos, h, folga):
    txt = open(geo_grosso).read()
    m = re.search(r"Background\s+Field\s*=\s*(\d+)\s*;", txt)
    if not m:
        raise SystemExit("refina_banda: " + geo_grosso + " nao tem 'Background Field = N;'")
    fundo = int(m.group(1))
    usados = [int(k) for k in re.findall(r"Field\[(\d+)\]", txt)]
    f1 = max(usados + [fundo]) + 1
    pontos_ids = [int(k) for k in re.findall(r"Point\((\d+)\)", txt)]
    p0 = max(pontos_ids + [0]) + 1
    if "hMax" not in txt:
        raise SystemExit("refina_banda: " + geo_grosso + " nao define hMax")
    novo = ["", "// refina_banda.py: elementos de %g m ate %g m da banda de ruptura (%d pontos) da rodada na malha grossa"
            % (h, folga, len(pontos)), "bandaH = %g; bandaFolga = %g;" % (h, folga)]
    novo += ["Point(%d) = {%.4f, %.4f, 0};" % (p0 + i, x, y) for i, (x, y) in enumerate(pontos)]
    novo += ["Field[%d] = Distance; Field[%d].PointsList = {%d:%d};" % (f1, f1, p0, p0 + len(pontos) - 1),
             "Field[%d] = Threshold; Field[%d].InField = %d; Field[%d].SizeMin = bandaH; Field[%d].SizeMax = hMax;"
             % (f1 + 1, f1 + 1, f1, f1 + 1, f1 + 1),
             "Field[%d].DistMin = bandaFolga; Field[%d].DistMax = bandaFolga + 3.0;" % (f1 + 1, f1 + 1),
             "Field[%d] = Min; Field[%d].FieldsList = {%d, %d};" % (f1 + 2, f1 + 2, fundo, f1 + 1),
             "Background Field = %d;" % (f1 + 2)]
    # o fundo antigo sai; os pontos e campos novos entram no lugar dele (antes das opcoes de malha e dos grupos)
    return txt[:m.start()] + "\n".join(novo) + "\n" + txt[m.end():]


def gera_msh(geo, msh):
    if shutil.which("gmsh"):
        r = subprocess.run(["gmsh", "-2", geo, "-o", msh], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                           errors="replace")
        if r.returncode != 0 or not os.path.isfile(msh):
            raise SystemExit("refina_banda: gmsh falhou:\n" + r.stdout[-2000:])
        return
    try:
        import gmsh
    except ImportError:
        raise SystemExit("refina_banda: nem o executavel gmsh nem o modulo python gmsh "
                         "(python3 -m pip install --user gmsh)")
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.open(geo)
        gmsh.model.mesh.generate(2)
        gmsh.write(msh)
    finally:
        gmsh.finalize()


def triangulos(msh):
    t = open(msh).read()
    i = t.index("$Elements")
    L = t[i:t.index("$EndElements")].split("\n")[2:]
    k = n = 0
    while k < len(L):
        p = L[k].split()
        if len(p) == 4:
            num = int(p[3])
            n += num if int(p[2]) == 2 else 0
            k += num + 1
        else:
            k += 1
    return n


def main():
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    rodada, grossa, saida = sys.argv[1], sys.argv[2], sys.argv[3]
    h = float(sys.argv[4]) if len(sys.argv) > 4 else 0.125
    folga = float(sys.argv[5]) if len(sys.argv) > 5 else 0.375
    limiar = float(sys.argv[6]) if len(sys.argv) > 6 else 0.1
    geo_grosso = os.path.splitext(grossa)[0] + ".geo"
    if not os.path.isfile(geo_grosso):
        sys.exit("refina_banda: falta o .geo da malha grossa: " + geo_grosso)
    celula = h
    pts = sorted({(round(x / celula) * celula, round(y / celula) * celula) for x, y in banda(rodada, limiar)})
    with open(saida + ".geo.tmp", "w") as f:
        f.write(geo_refinado(geo_grosso, pts, h, folga))
    os.replace(saida + ".geo.tmp", saida + ".geo")
    gera_msh(saida + ".geo", saida + ".msh.tmp.msh")
    os.replace(saida + ".msh.tmp.msh", saida + ".msh")
    print("%s.msh: %d triangulos, %d pontos de banda (h = %g, folga = %g, limiar = %g)" % (
        saida, triangulos(saida + ".msh"), len(pts), h, folga, limiar))


if __name__ == "__main__":
    main()
