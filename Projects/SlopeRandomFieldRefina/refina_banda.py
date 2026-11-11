#!/usr/bin/env python3
"""refina_banda.py <rodada grossa> <malha grossa .msh ou .geo> <saida sem extensao> [h] [folga] [limiar]
refina_banda.py --malha <saida sem extensao>: so gera <saida>.msh de novo, a partir de <saida>.geo
refina_banda.py --gmsh: so mostra qual gmsh sera usado (codigo 2 se nao ha nenhum)

A segunda etapa de srf_refinado.sh: uma malha conforme refinada so em volta da banda de ruptura que a rodada na malha
grossa encontrou, gerada pelo gmsh a partir do .geo da malha grossa.

- A banda: os pontos de integracao com sqrt(J2(eps_p)) >= limiar x maximo (padrao 0,1, o BANDA_LIMIAR do programa) no
  instante do fator (o menor FS de fs.csv), lidos de <rodada grossa>/mecanismo.csv e arredondados numa grade de h.
- O .geo novo e o da malha grossa (mesma geometria, mesmos grupos fisicos, mesmos campos de tamanho) com mais um campo:
  elementos de tamanho h (padrao 0,125 m) ate 'folga' (padrao 0,375 m) dos pontos da banda, crescendo ate o hMax da malha
  grossa a folga + 3 m; o campo de fundo passa a ser o minimo do antigo e deste. Assim o contorno, e com ele o retangulo
  do Karhunen-Loeve, nao muda, e a realizacao tem o mesmo campo aleatorio nas duas malhas.
- O gmsh: o modulo python gmsh deste python (python3 -m pip install --user gmsh); sem ele, ou se ele nao carrega, o
  executavel gmsh do PATH, se ele roda (gmsh -version). O executavel gmsh que o pip instala junto com o modulo e um
  script com '#!/usr/bin/env python', que nao roda onde nao ha um 'python' no PATH; por isso o modulo vem primeiro.

Escreve <saida>.geo e <saida>.msh (cada um com troca atomica) e imprime uma linha chave=valor: triangulos,
celulas_banda, fator_grossa (o fator no instante da banda), t_h, h, folga e limiar.
Codigos de saida: 0 ok; 1 erro; 2 uso errado ou sem gmsh; 3 a rodada grossa nao tem banda (sem fator positivo no
instante do menor FS, ou sem deformacao plastica)."""

import csv
import os
import re
import shutil
import subprocess
import sys

SEM_BANDA = 3


class SemBanda(Exception):
    pass


def banda(rodada, limiar):
    """(pontos da banda, fator, instante) da rodada grossa"""
    linhas = [r for r in csv.DictReader(open(os.path.join(rodada, "fs.csv"))) if r.get("FS_SRM")]
    if not linhas:
        raise SemBanda("fs.csv sem fator em " + rodada)
    menor = min(linhas, key=lambda r: float(r["FS_SRM"]))
    fator, t_fs = float(menor["FS_SRM"]), float(menor["t_h"])
    if not fator > 0:
        raise SemBanda("fator %g em t = %g h em %s" % (fator, t_fs, rodada))
    pts = []
    for r in csv.DictReader(open(os.path.join(rodada, "mecanismo.csv"))):
        if abs(float(r["t_h"]) - t_fs) < 1e-9:
            pts.append((float(r["x"]), float(r["y"]), float(r["epsp_raizJ2"])))
    mx = max((p[2] for p in pts), default=0.)
    if not mx > 0:
        raise SemBanda("sem deformacao plastica em t = %g h em %s" % (t_fs, rodada))
    return [(x, y) for x, y, e in pts if e >= limiar * mx], fator, t_fs


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


def modulo_gmsh():
    """(o modulo python gmsh, None), ou (None, motivo) sem o modulo ou sem as bibliotecas de que ele depende"""
    try:
        import gmsh
    except ImportError:
        return None, "sem o modulo python gmsh em %s" % sys.executable
    except OSError as e:  # o libgmsh.so do pip depende de bibliotecas do X e do OpenGL (libGLU, libXft, ...)
        return None, "o modulo python gmsh nao carrega: %s" % e
    return gmsh, None


def executavel_gmsh():
    """(caminho, versao) do executavel gmsh do PATH, se 'gmsh -version' roda; ou (None, motivo)"""
    g = shutil.which("gmsh")
    if not g:
        return None, "sem o executavel gmsh no PATH"
    try:
        r = subprocess.run([g, "-version"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                           errors="replace", timeout=120)
    except (OSError, subprocess.SubprocessError) as e:
        return None, "o executavel %s nao roda: %s" % (g, e)
    if r.returncode != 0:  # a ultima linha da saida: o erro (de um script python, a ultima do traceback)
        fim = r.stdout.strip().splitlines()[-1] if r.stdout.strip() else ""
        return None, "o executavel %s nao roda (codigo %d): %s" % (g, r.returncode, fim)
    return g, r.stdout.strip().splitlines()[-1] if r.stdout.strip() else "-"


_GMSH = []


def qual_gmsh():
    """('modulo', modulo) ou ('executavel', caminho, versao): o modulo python gmsh ou, sem ele, o executavel gmsh do
    PATH que roda; sai com codigo 2 se nao ha nenhum"""
    if not _GMSH:
        mod, motivo_mod = modulo_gmsh()
        if mod:
            _GMSH.append(("modulo", mod))
        else:
            exe, info = executavel_gmsh()
            if not exe:
                print("refina_banda: sem gmsh: %s; %s. Instale o modulo (%s -m pip install --user gmsh) ou ponha no "
                      "PATH um executavel gmsh que rode" % (motivo_mod, info, os.path.basename(sys.executable)),
                      file=sys.stderr)
                sys.exit(2)
            _GMSH.append(("executavel", exe, info))
    return _GMSH[0]


def gera_msh(geo, msh):
    g = qual_gmsh()
    if g[0] == "executavel":
        r = subprocess.run([g[1], "-2", geo, "-o", msh], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                           errors="replace")
        if r.returncode != 0 or not os.path.isfile(msh):
            raise SystemExit("refina_banda: gmsh falhou:\n" + r.stdout[-2000:])
        return
    gmsh = g[1]
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.open(geo)
        gmsh.model.mesh.generate(2)
        gmsh.write(msh)
    finally:
        gmsh.finalize()


def msh_atomico(saida):
    """<saida>.msh a partir de <saida>.geo, com troca atomica (uma queda no meio nao deixa .msh incompleto)"""
    gera_msh(saida + ".geo", saida + ".msh.tmp.msh")
    os.replace(saida + ".msh.tmp.msh", saida + ".msh")


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
    if len(sys.argv) == 2 and sys.argv[1] == "--gmsh":
        g = qual_gmsh()
        if g[0] == "modulo":
            print("modulo %s (%s), em %s" % (g[1].__file__, g[1].__version__, sys.executable))
        else:
            print("executavel %s (%s)" % (g[1], g[2]))
        return
    if len(sys.argv) == 3 and sys.argv[1] == "--malha":
        saida = sys.argv[2]
        msh_atomico(saida)
        print("malha=%s.msh triangulos=%d" % (saida, triangulos(saida + ".msh")))
        return
    if len(sys.argv) < 4:
        print(__doc__, file=sys.stderr)
        sys.exit(2)
    rodada, grossa, saida = sys.argv[1], sys.argv[2], sys.argv[3]
    h = float(sys.argv[4]) if len(sys.argv) > 4 else 0.125
    folga = float(sys.argv[5]) if len(sys.argv) > 5 else 0.375
    limiar = float(sys.argv[6]) if len(sys.argv) > 6 else 0.1
    geo_grosso = os.path.splitext(grossa)[0] + ".geo"
    if not os.path.isfile(geo_grosso):
        sys.exit("refina_banda: falta o .geo da malha grossa: " + geo_grosso)
    try:
        pontos, fator, t_fs = banda(rodada, limiar)
    except SemBanda as e:
        print("refina_banda: sem banda de ruptura: %s" % e, file=sys.stderr)
        sys.exit(SEM_BANDA)
    pts = sorted({(round(x / h) * h, round(y / h) * h) for x, y in pontos})
    with open(saida + ".geo.tmp", "w") as f:
        f.write(geo_refinado(geo_grosso, pts, h, folga))
    os.replace(saida + ".geo.tmp", saida + ".geo")
    msh_atomico(saida)
    print("malha=%s.msh triangulos=%d celulas_banda=%d fator_grossa=%.10g t_h=%g h=%g folga=%g limiar=%g" % (
        saida, triangulos(saida + ".msh"), len(pts), fator, t_fs, h, folga, limiar))


if __name__ == "__main__":
    main()
