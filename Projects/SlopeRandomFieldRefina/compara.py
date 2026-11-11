#!/usr/bin/env python3
"""compara.py: tempo e precisao do metodo de duas etapas (srf_refinado.sh) contra malhas fixas, no exemplo saturado do
Ceron (o das fases do fermi: ../SlopeRandomField/README.md, "Validacao no limite saturado"), GIM com arc-length em 0,01 h.

compara.py rodar <pasta> [--casos det,0-7] [--metodos foco125,refinada,foco1875,foco25,foco375] [--paralelo 2]
    Roda cada caso (det: o deterministico, CAMPOS=0; n: a realizacao n) com cada metodo, em <pasta>/<metodo>/<caso>/.
    - Uma fila so, com PARALELO processos ao mesmo tempo (padrao 2), intercalada: para cada caso, todos os metodos, na
      ordem de --metodos. Assim todo metodo roda ao lado do mesmo numero de processos desta fila, e uma variacao da carga
      da maquina ao longo da rodada pesa sobre todos os metodos.
    - Metodos: o nome de uma malha de ../LearningMohrCoulomb2/cluster/malhas (foco375 e ceron_h5_a45_foco375.msh) roda o
      programa nessa malha; "refinada" roda srf_refinado.sh com a malha grossa MALHA_GROSSA (padrao foco375).
    - tempo.txt de cada caso: parede_s (relogio, do inicio ao fim do processo), cpu_s (usuario + sistema do processo e
      dos seus filhos, de wait4), inicio e fim (s desde 1970), paralelo e a carga media da maquina (/proc/loadavg) no
      inicio e no fim. Um caso com resultado.txt e tempo.txt e pulado; um caso sem tempo.txt comeca do zero (a pasta e
      apagada), para o tempo ser o de uma execucao inteira. Ctrl-C para tudo; rodar de novo continua.
compara.py tabela <pasta> [--ref foco125] [--csv <arquivo>]
    Tabelas em Markdown (decimais com virgula): o fator de cada metodo e a diferenca para a referencia, caso a caso; e, por
    metodo, triangulos, equacoes, diferenca media, maxima e media com sinal, tempo de parede e de CPU por caso, a razao
    para o tempo da referencia e o numero medio de processos da fila rodando junto. --csv grava os dados de cada caso.
Ambiente: SRF_BIN (o executavel, padrao ../SlopeRandomField/bin/SlopeRandomField), SOLVER (passado ao programa se
definido), MALHA_GROSSA e os REFINA_* de srf_refinado.sh."""

import argparse
import collections
import csv
import math
import os
import shutil
import subprocess
import sys
import threading
import time

AQUI = os.path.dirname(os.path.abspath(__file__))
MALHAS = os.path.normpath(os.path.join(AQUI, "..", "LearningMohrCoulomb2", "cluster", "malhas"))
SRF_BIN = os.environ.get("SRF_BIN", os.path.normpath(os.path.join(AQUI, "..", "SlopeRandomField", "bin", "SlopeRandomField")))
REFINADO = os.path.join(AQUI, "srf_refinado.sh")
# o exemplo saturado do Ceron, como nas fases do fermi (../SlopeRandomField/README.md)
MODELO = dict(SOLO="f50", I_MM_H="0", VARIANTE="iii", PHI_DEG="30", HL="10", HR="10", XL="0", XR="25", XTOE="15",
              YTOE="5", XCRISTA="10", YCRISTA="10", NIVEL_FINAL="5", T_REBAIX_H="0", FS_H="0.01", DT_H="0.01",
              FATOR="gim_arco", SATURADO="1", SUPERFICIE_P0="1", KW="1e30", MKL_NUM_THREADS="1", OMP_NUM_THREADS="1")


def malha(nome):
    for f in ("ceron_h5_a45_%s.msh" % nome, nome + ".msh", nome):
        if os.path.isfile(os.path.join(MALHAS, f)):
            return os.path.join(MALHAS, f)
    sys.exit("compara: malha %s nao encontrada em %s" % (nome, MALHAS))


def casos(txt):
    r = []
    for p in txt.split(","):
        if p == "det":
            r.append("det")
        elif "-" in p:
            a, b = p.split("-")
            r += [str(i) for i in range(int(a), int(b) + 1)]
        else:
            r.append(str(int(p)))
    return r


def pasta_caso(base, metodo, caso):
    return os.path.join(base, metodo, "det" if caso == "det" else "r%05d" % int(caso))


def chaves(arq):
    """chave=valor de uma linha de texto"""
    try:
        return dict(p.split("=", 1) for p in open(arq).read().split() if "=" in p)
    except OSError:
        return {}


def linha_csv(arq):
    try:
        return next(csv.DictReader(open(arq)))
    except (OSError, StopIteration):
        return {}


def carga():
    try:
        return open("/proc/loadavg").read().split()[0]
    except OSError:
        return "-"


def versao():
    try:
        return open(os.path.join(os.path.dirname(SRF_BIN), "versao.txt")).read().strip()
    except OSError:
        return "-"


def rodar(a):
    base = os.path.abspath(a.pasta)
    metodos = a.metodos.split(",")
    for m in metodos:
        if m != "refinada":
            malha(m)
    grossa = malha(os.environ.get("MALHA_GROSSA", "foco375"))
    if not os.access(SRF_BIN, os.X_OK):
        sys.exit("compara: executavel %s nao encontrado (SRF_BIN)" % SRF_BIN)
    fila = collections.deque()
    for c in casos(a.casos):
        for m in metodos:
            d = pasta_caso(base, m, c)
            if os.path.isfile(os.path.join(d, "resultado.txt")) and os.path.isfile(os.path.join(d, "tempo.txt")):
                continue
            fila.append((m, c, d))
    print("%s: %d casos a rodar, %d ao mesmo tempo, programa %s (%s)" % (base, len(fila), a.paralelo, SRF_BIN, versao()),
          flush=True)
    trava = threading.Lock()
    parar = threading.Event()

    def um(m, c, d):
        if os.path.isdir(d):
            shutil.rmtree(d)  # sem tempo.txt: do zero, para o tempo ser o de uma execucao inteira
        os.makedirs(d)
        env = dict(os.environ)
        env.update(MODELO)
        env.update(CAMPOS="0" if c == "det" else "1", REAL="0" if c == "det" else c, SRF_VERSAO=versao(), SRF_BIN=SRF_BIN)
        if m == "refinada":
            cmd, env["MESH"] = [REFINADO], grossa
        else:
            cmd, env["MESH"] = [SRF_BIN], malha(m)
        c0, t0 = carga(), time.time()
        with open(os.path.join(d, "run.log"), "a") as out, open(os.path.join(d, "run.err"), "a") as err:
            p = subprocess.Popen(cmd, cwd=d, env=env, stdout=out, stderr=err, stdin=subprocess.DEVNULL)
            _, status, uso = os.wait4(p.pid, 0)
        t1 = time.time()
        p.returncode = os.waitstatus_to_exitcode(status)
        ok = os.path.isfile(os.path.join(d, "resultado.txt"))
        if p.returncode < 0 or p.returncode >= 128:
            parar.set()  # um sinal (Ctrl-C): para a fila
            ok = False
        if ok:
            with open(os.path.join(d, "tempo.txt.tmp"), "w") as f:
                f.write("parede_s=%.3f cpu_s=%.3f inicio=%.3f fim=%.3f paralelo=%d carga_inicio=%s carga_fim=%s\n"
                        % (t1 - t0, uso.ru_utime + uso.ru_stime, t0, t1, a.paralelo, c0, carga()))
            os.replace(os.path.join(d, "tempo.txt.tmp"), os.path.join(d, "tempo.txt"))
        r = chaves(os.path.join(d, "resultado.txt"))
        with trava:
            print("%s %s %s: fator %s, %.1f s (codigo %d)" % (time.strftime("%H:%M:%S"), m, c, r.get("FS_min", "-"),
                                                             t1 - t0, p.returncode), flush=True)

    def trabalhador():
        while not parar.is_set():
            with trava:
                if not fila:
                    return
                m, c, d = fila.popleft()
            try:
                um(m, c, d)
            except Exception as e:  # o caso fica sem tempo.txt e e rodado de novo na proxima vez
                with trava:
                    print("compara: erro em %s %s: %s" % (m, c, e), flush=True)

    ts = [threading.Thread(target=trabalhador, daemon=True) for _ in range(a.paralelo)]
    for t in ts:
        t.start()
    try:
        while any(t.is_alive() for t in ts):
            time.sleep(0.5)
    except KeyboardInterrupt:
        parar.set()
        print("compara: interrompido; os casos em andamento ficam sem tempo.txt e recomecam na proxima vez", flush=True)
        for t in ts:
            t.join()
        sys.exit(130)


def media(v):
    v = [x for x in v if x is not None]
    return sum(v) / len(v) if v else None


def fmt(x, d=4, sinal=False):
    if x is None:
        return "-"
    s = ("%+.*f" if sinal else "%.*f") % (d, x)
    return s.replace(".", ",").replace("-", "−")


def banda(pasta, limiar=0.1):
    """os pontos (x, y) da banda de ruptura de uma rodada: sqrt(J2) >= limiar do maximo no instante do menor fator"""
    fs = [r for r in csv.DictReader(open(os.path.join(pasta, "fs.csv"))) if r.get("FS_SRM")]
    t = float(min(fs, key=lambda r: float(r["FS_SRM"]))["t_h"])
    pts = [(float(r["x"]), float(r["y"]), float(r["epsp_raizJ2"]))
           for r in csv.DictReader(open(os.path.join(pasta, "mecanismo.csv"))) if abs(float(r["t_h"]) - t) < 1e-9]
    mx = max(e for _, _, e in pts)
    return [(x, y) for x, y, e in pts if e >= limiar * mx]


def foco(geo):
    """(pontos do foco, folga) do .geo de uma malha focada (os Point(1000...) do campo de distancia), ou None"""
    import re
    txt = open(geo).read()
    f = re.search(r"folga\s*=\s*([0-9.]+)", txt)
    pts = [(float(x), float(y)) for i, x, y in
           re.findall(r"Point\((\d+)\)\s*=\s*\{\s*([-0-9.eE]+),\s*([-0-9.eE]+),", txt) if int(i) >= 1000]
    return (pts, float(f.group(1))) if f and pts else None


def celulas_refino(geo):
    """(celulas da banda da malha grossa, folga) de um malha_refinada.geo de refina_banda.py"""
    import re
    txt = open(geo).read()
    i = txt.index("// refina_banda.py")
    f = re.search(r"bandaFolga\s*=\s*([0-9.eE+-]+)", txt)
    pts = re.findall(r"Point\((\d+)\)\s*=\s*\{\s*([-0-9.eE]+),\s*([-0-9.eE]+),", txt[i:])
    return [(float(x), float(y)) for _, x, y in pts], float(f.group(1))


def tabela(a):
    base = os.path.abspath(a.pasta)
    metodos = sorted(m for m in os.listdir(base) if os.path.isdir(os.path.join(base, m)))
    if a.ref not in metodos:
        sys.exit("compara: sem a referencia %s em %s" % (a.ref, base))
    dados = {}
    for m in metodos:
        for d in sorted(os.listdir(os.path.join(base, m))):
            p = os.path.join(base, m, d)
            r = chaves(os.path.join(p, "resultado.txt"))
            if not r.get("FS_min"):
                continue
            c = "det" if d == "det" else str(int(d[1:]))
            t = chaves(os.path.join(p, "tempo.txt"))
            info = chaves(os.path.join(p, "realizacao.txt"))
            try:
                tri = sum(1 for _ in open(os.path.join(p, "triangulos.csv"))) - 1
            except OSError:
                tri = None
            ref = linha_csv(os.path.join(p, "refino.csv"))
            num = lambda k, dic: float(dic[k]) if dic.get(k) not in (None, "") else None
            fora = dist = None
            if m == "refinada" and os.path.isfile(os.path.join(p, "malha_refinada.geo")):
                try:
                    cel, folga = celulas_refino(os.path.join(p, "malha_refinada.geo"))
                    d2 = [min(math.hypot(x - u, y - v) for u, v in cel) for x, y in banda(p)]
                    fora, dist = sum(x > folga for x in d2) / len(d2), max(d2)
                except (OSError, ValueError):
                    pass
            dados[m, c] = dict(fora_refino=fora, dist_refino=dist, metodo=m, caso=c, fator=float(r["FS_min"]), seg=float(r["seg"]),
                               parede_s=num("parede_s", t), cpu_s=num("cpu_s", t), inicio=num("inicio", t),
                               fim=num("fim", t), neq=int(info["neq"]) if "neq" in info else None, triangulos=tri,
                               fator_grossa=num("fator_grossa", ref), seg_grossa=num("seg_grossa", ref),
                               seg_malha=num("seg_malha", ref), seg_fina=num("seg_fina", ref),
                               celulas_banda=num("celulas_banda", ref), metodo_refino=ref.get("metodo", ""))
    # processos da fila rodando junto com cada caso (a media no tempo, contando ele mesmo)
    todos = [v for v in dados.values() if v["inicio"] is not None]
    for v in todos:
        dur = v["fim"] - v["inicio"]
        sob = sum(max(0., min(v["fim"], u["fim"]) - max(v["inicio"], u["inicio"])) for u in todos)
        v["simultaneos"] = sob / dur if dur > 0 else None
    cs = sorted({c for (_, c) in dados}, key=lambda c: -1 if c == "det" else int(c))
    cs = [c for c in cs if (a.ref, c) in dados]
    # as malhas fixas da mais grossa a mais fina, o metodo de duas etapas e a referencia
    fixas = sorted((m for m in metodos if m not in (a.ref, "refinada")),
                   key=lambda m: media([dados[m, c]["neq"] for c in cs if (m, c) in dados]) or 0)
    ordem = fixas + (["refinada"] if "refinada" in metodos else []) + [a.ref]
    for v in dados.values():
        r = dados.get((a.ref, v["caso"]))
        v["dif"] = (v["fator"] - r["fator"]) / r["fator"] if r else None
    if a.csv:
        cols = ["metodo", "caso", "fator", "dif", "seg", "parede_s", "cpu_s", "simultaneos", "neq", "triangulos",
                "fator_grossa", "seg_grossa", "seg_malha", "seg_fina", "celulas_banda", "metodo_refino", "fora_refino",
                "dist_refino"]
        with open(a.csv, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(cols)
            for m in ordem:
                for c in cs:
                    if (m, c) in dados:
                        w.writerow(["" if dados[m, c].get(k) is None else dados[m, c].get(k) for k in cols])
    print("Fator (Γ) de cada caso e a diferença para a referência (%s):\n" % a.ref)
    print("| caso | " + " | ".join(ordem) + " |")
    print("|---|" + "---|" * len(ordem))
    for c in cs:
        cel = []
        for m in ordem:
            v = dados.get((m, c))
            if v is None:
                cel.append("-")
            elif m == a.ref:
                cel.append(fmt(v["fator"]))
            else:
                cel.append("%s (%s%%)" % (fmt(v["fator"]), fmt(100 * v["dif"], 2, True)))
        print("| %s | %s |" % (c, " | ".join(cel)))
    print("\nResumo por método (%d casos: %s):\n" % (len(cs), ", ".join(cs)))
    print("| método | triângulos | equações | dif. média | dif. máxima | dif. média com sinal | Γ < 1 | "
          "tempo de parede (s) | tempo de CPU (s) | parede / %s | processos juntos |" % a.ref)
    print("|---|---|---|---|---|---|---|---|---|---|---|")
    tref = media([dados[a.ref, c]["parede_s"] for c in cs])
    for m in ordem:
        vs = [dados[m, c] for c in cs if (m, c) in dados]
        if len(vs) < len(cs):
            continue
        difs = [v["dif"] for v in vs]
        tp, tc = media([v["parede_s"] for v in vs]), media([v["cpu_s"] for v in vs])
        print("| %s | %s | %s | %s%% | %s%% | %s%% | %d | %s | %s | %s | %s |" % (
            m, fmt(media([v["triangulos"] for v in vs]), 0), fmt(media([v["neq"] for v in vs]), 0),
            fmt(100 * media([abs(x) for x in difs]), 2), fmt(100 * max(abs(x) for x in difs), 2),
            fmt(100 * media(difs), 2, True), sum(v["fator"] < 1 for v in vs), fmt(tp, 1), fmt(tc, 1),
            fmt(tp / tref if tp and tref else None, 2), fmt(media([v.get("simultaneos") for v in vs]), 2)))
    # a banda da referencia e a regiao focada das malhas fixas (hFoco ate 'folga' dos pontos do foco)
    geo = os.path.splitext(malha(a.ref))[0] + ".geo"
    fc = foco(geo) if os.path.isfile(geo) else None
    if fc:
        dist = []
        for c in cs:
            try:
                b = banda(pasta_caso(base, a.ref, c))
            except (OSError, ValueError):
                continue
            dist.append(max(min(math.hypot(x - u, y - v) for u, v in fc[0]) for x, y in b))
        if dist:
            print("\nBanda de ruptura da referência (√J2 ≥ 0,1 do máximo) e a região focada de %s: em %d de %d casos, "
                  "algum ponto da banda a mais de %s m dos pontos do foco; distância máxima %s m." % (
                      os.path.basename(geo), sum(d > fc[1] for d in dist), len(dist), fmt(fc[1], 2),
                      fmt(max(dist), 2)))
    vs = [dados["refinada", c] for c in cs if ("refinada", c) in dados]
    if vs:
        print("\nAs etapas do método de duas etapas (médias por caso, tempos do relógio de cada etapa):\n")
        print("| fator na malha grossa: dif. média | células de banda | triângulos refinados | equações refinadas | "
              "malha grossa (s) | gmsh (s) | malha refinada (s) | sem banda | banda final fora da região refinada |")
        print("|---|---|---|---|---|---|---|---|---|")
        dg = [(v["fator_grossa"] - dados[a.ref, v["caso"]]["fator"]) / dados[a.ref, v["caso"]]["fator"] for v in vs
              if v["fator_grossa"] is not None]
        fr = [v for v in vs if v["fora_refino"] is not None]
        print("| %s%% | %s | %s | %s | %s | %s | %s | %d | %s |" % (
            fmt(100 * media([abs(x) for x in dg]), 2), fmt(media([v["celulas_banda"] for v in vs]), 0),
            fmt(media([v["triangulos"] for v in vs]), 0), fmt(media([v["neq"] for v in vs]), 0),
            fmt(media([v["seg_grossa"] for v in vs]), 1), fmt(media([v["seg_malha"] for v in vs]), 2),
            fmt(media([v["seg_fina"] for v in vs]), 1), sum(v["metodo_refino"] == "grossa" for v in vs),
            "%d de %d casos (no pior, %s%% dos pontos; distância máxima %s m)" % (
                sum(v["fora_refino"] > 0 for v in fr), len(fr), fmt(100 * max(v["fora_refino"] for v in fr), 1),
                fmt(max(v["dist_refino"] for v in fr), 2)) if fr else "-"))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("rodar")
    r.add_argument("pasta")
    r.add_argument("--casos", default="det,0-7")
    r.add_argument("--metodos", default="foco125,refinada,foco1875,foco25,foco375")
    r.add_argument("--paralelo", type=int, default=2)
    t = sub.add_parser("tabela")
    t.add_argument("pasta")
    t.add_argument("--ref", default="foco125")
    t.add_argument("--csv")
    a = ap.parse_args()
    if a.cmd == "rodar":
        rodar(a)
    else:
        tabela(a)


if __name__ == "__main__":
    main()
