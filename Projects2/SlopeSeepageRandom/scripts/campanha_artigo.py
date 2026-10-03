#!/usr/bin/env python3
"""Campanha de Monte Carlo de todos os casos do artigo (Vargas Ceron et al., IJNAMG 2025), em blocos retomáveis.

Cada caso é dividido em blocos de amostras (um CSV por bloco: <dir>/<caso>/b<k>.csv); a fila intercala os casos
(rodada r = bloco r de todos os casos), de modo que todos avançam juntos e os resultados parciais são sempre
utilizáveis. Rodar de novo a mesma fila retoma: blocos completos são pulados e o executável pula as amostras já
gravadas num bloco interrompido. Não mude --bloco depois de começar (os blocos são numerados por ele).

  campanha_artigo.py jobs   <executável> <dir> [--alvo artigo|N|covX] [--min 1000] [--bloco 50] [--casos a,b]
                            [--adapt 2] > jobs.txt
  fila.sh jobs.txt P                                   # P processos (um por núcleo físico)
  PROCESSOS=P campanha_artigo.py status <dir>          # amostras por caso, Pf, ritmo e previsão de tempo
  campanha_artigo.py junta  <dir>                      # <dir>/<caso>.csv, .mec, .modo e .param (todos os blocos)
  campanha_artigo.py pacote <dir> [arquivo.tar.gz]     # só os resultados (sem caches): para enviar/arquivar

--alvo artigo: o número de amostras S de cada caso no artigo (critério CoV(Pf) < 5 %, 10 000 a 100 000);
--alvo N: N por caso (limitado ao S do artigo);
--alvo covX (p.ex. cov5): o protocolo do artigo com o Pf deste código: N tal que CoV(Pf) = sqrt((1-Pf)/(N Pf))
  < X %, com Pf estimado das amostras já calculadas (casos com menos de --min amostras vão até --min; casos com
  Pf = 0 usam o S do artigo). Gere a fila de novo de tempos em tempos: o alvo é reavaliado com as amostras novas.
"""
import csv
import glob
import math
import os
import sys
import tarfile

# γw = 9.81 kN/m³: os valores do artigo (Fig. 8 em h_w = 0, Fig. 13) correspondem a 9.81, não a 10
BASE = "caso=percolacao modelo=mc h=1 adapt=2 hkl=1 covc=0.3 covphi=0.1 covk=0.6 Lx=20 Ly=2 gw=9.81"

# nome: (argumentos, S do artigo, referência)
CASOS = {
    "ref": (BASE, 10000, "Tab. 2-6, ref."),
    # Tabela 3: coeficientes de variação
    "covk0": (BASE + " covk=0", 10000, "Tab. 3"),
    "covk75": (BASE + " covk=0.75", 10000, "Tab. 3"),
    "covk90": (BASE + " covk=0.9", 10000, "Tab. 3"),
    "covk100": (BASE + " covk=1.0", 10000, "Tab. 3"),
    "covc10": (BASE + " covc=0.1", 90000, "Tab. 3"),
    "covc50": (BASE + " covc=0.5", 10000, "Tab. 3"),
    "covc70": (BASE + " covc=0.7", 10000, "Tab. 3"),
    "covphi5": (BASE + " covphi=0.05", 10000, "Tab. 3"),
    "covphi15": (BASE + " covphi=0.15", 10000, "Tab. 3"),
    "covphi20": (BASE + " covphi=0.2", 10000, "Tab. 3"),
    # Tabela 4: (Lx, Ly) = s (20 m, 2 m)
    "s1.5": (BASE + " Lx=30 Ly=3", 10000, "Tab. 4"),
    "s2": (BASE + " Lx=40 Ly=4", 10000, "Tab. 4"),
    "s5": (BASE + " Lx=100 Ly=10", 10000, "Tab. 4"),
    "s10": (BASE + " Lx=200 Ly=20", 10000, "Tab. 4"),
    "s20": (BASE + " Lx=400 Ly=40", 10000, "Tab. 4"),
    "s400": (BASE + " Lx=8000 Ly=800", 10000, "Tab. 4"),
    # Tabela 5: anisotropia α = kh/kv
    "alfa2": (BASE + " alpha=2", 20000, "Tab. 5"),
    "alfa3": (BASE + " alpha=3", 30000, "Tab. 5"),
    "alfa4": (BASE + " alpha=4", 40000, "Tab. 5"),
    "alfa5": (BASE + " alpha=5", 70000, "Tab. 5"),
    # Tabela 6: rebaixamento h_w/H (H = 5 m)
    "hw0.5": (BASE + " hw=2.5", 30000, "Tab. 6"),
    "hw0.6": (BASE + " hw=3", 10000, "Tab. 6"),
    "hw0.7": (BASE + " hw=3.5", 10000, "Tab. 6"),
    "hw0.8": (BASE + " hw=4", 10000, "Tab. 6"),
    "hw0.9": (BASE + " hw=4.5", 10000, "Tab. 6"),
    # seção 5.3: Cho (2010), sem percolação
    "cho_coesivo": ("caso=cho_coesivo modelo=mc h=1 adapt=2 hkl=1 covc=0.3 covphi=0 covk=0 Lx=20 Ly=2", 100000,
                    "seç. 5.3.1"),
    "cho_cphi": ("caso=cho_cphi modelo=mc h=2 adapt=2 hkl=1 covc=0.3 covphi=0.2 covk=0 Lx=20 Ly=2", 50000,
                 "seç. 5.3.2"),
}


def opt(argv, key, default):
    if key in argv:
        i = argv.index(key)
        v = argv[i + 1]
        del argv[i:i + 2]
        return v
    return default


def linhas(path):
    try:
        with open(path) as f:
            return max(sum(1 for _ in f) - 1, 0)
    except OSError:
        return 0


def pf_atual(d):
    rows, _ = ler(d) if os.path.isdir(d) else ({}, {})
    n = len(rows)
    return n, (sum(float(r["fator"]) < 1. for r in rows.values()) / n if n else 0.)


def jobs(argv):
    alvo = opt(argv, "--alvo", "artigo")
    bloco = int(opt(argv, "--bloco", "50"))
    nmin = int(opt(argv, "--min", "1000"))
    adapt = opt(argv, "--adapt", "")  # níveis de refinamento da malha (padrão 2; mude só num diretório novo)
    casos = opt(argv, "--casos", "")
    exe, out = os.path.abspath(argv[0]), os.path.abspath(argv[1])
    nomes = [c for c in CASOS if not casos or c in casos.split(",")]
    args = {c: CASOS[c][0] if not adapt else CASOS[c][0].replace("adapt=2", f"adapt={adapt}") for c in nomes}
    alvos = {}
    for c in nomes:
        S = CASOS[c][1]
        if alvo == "artigo":
            alvos[c] = S
        elif alvo.startswith("cov"):
            x = float(alvo[3:]) / 100.
            n, pf = pf_atual(os.path.join(out, c))
            if n < nmin:
                alvos[c] = nmin
            elif pf <= 0.:
                alvos[c] = S
            else:
                alvos[c] = max(nmin, int(math.ceil((1. - pf) / (pf * x * x))))
            alvos[c] = -(-alvos[c] // bloco) * bloco
            print(f"# {c}: {n} amostras, Pf = {pf * 100:.2f} % -> alvo {alvos[c]} (artigo: {S})", file=sys.stderr)
        else:
            alvos[c] = min(int(alvo), S)
    for c in nomes:
        os.makedirs(os.path.join(out, c), exist_ok=True)
    # preparo: malha adaptada e autopares da KL (caches no diretório do caso) antes dos blocos
    for c in nomes:
        d = os.path.join(out, c)
        if not glob.glob(os.path.join(d, "kl_*_v3.bin")):
            print(f"cd {d} && {exe} mc {args[c]} n=0 saida=preparo.csv > preparo.log 2>&1")
    rodadas = max((alvos[c] + bloco - 1) // bloco for c in nomes)
    for r in range(rodadas):
        for c in nomes:
            ini = r * bloco
            if ini >= alvos[c]:
                continue
            n = min(bloco, alvos[c] - ini)
            d = os.path.join(out, c)
            if linhas(os.path.join(d, f"b{r}.csv")) >= n:
                continue
            print(f"cd {d} && {exe} mc {args[c]} inicio={ini} n={n} saida=b{r}.csv >> b{r}.log 2>&1")


def ler(d):
    rows, mec = {}, {}
    for f in glob.glob(os.path.join(d, "b*.csv")):
        with open(f) as fh:
            for r in csv.DictReader(fh):
                try:
                    rows[int(r["amostra"])] = r
                except (ValueError, KeyError):
                    pass
        try:
            with open(f + ".mec") as fh:
                for r in csv.DictReader(fh):
                    mec[int(r["amostra"])] = r
        except OSError:
            pass
    return rows, mec


def status(argv):
    out = argv[0]
    tot, ttot, falta_s = 0, 0.0, 0.0
    print(f"{'caso':12s} {'N':>7s} {'S artigo':>8s} {'μ':>7s} {'σ':>7s} {'Pf %':>7s} {'CoV(Pf)%':>8s} {'s/amostra':>9s}")
    for c, (_, S, _) in CASOS.items():
        d = os.path.join(out, c)
        if not os.path.isdir(d):
            continue
        rows, _ = ler(d)
        n = len(rows)
        if n == 0:
            print(f"{c:12s} {0:7d} {S:8d}")
            continue
        g = [float(r["fator"]) for r in rows.values()]
        t = [float(r["tempo_s"]) for r in rows.values()]
        mu = sum(g) / n
        sd = (sum((x - mu) ** 2 for x in g) / max(n - 1, 1)) ** 0.5
        pf = sum(x < 1 for x in g) / n
        cv = ((1 - pf) / (n * pf)) ** 0.5 * 100 if pf > 0 else float("nan")
        tm = sum(t) / n
        tot += n
        ttot += sum(t)
        falta_s += max(S - n, 0) * tm
        print(f"{c:12s} {n:7d} {S:8d} {mu:7.3f} {sd:7.3f} {pf * 100:7.2f} {cv:8.1f} {tm:9.2f}")
    print(f"total {tot} amostras, {ttot / 3600:.1f} h de CPU; até o S do artigo faltam {falta_s / 3600:.0f} h de CPU")
    # previsão: tempo médio por amostra de cada caso (casos sem amostras: média geral), P processos
    P = int(os.environ.get("PROCESSOS", "4"))
    tg = ttot / tot if tot else 5.
    tc = {}
    for c in CASOS:
        rows, _ = ler(os.path.join(out, c)) if os.path.isdir(os.path.join(out, c)) else ({}, {})
        tc[c] = (len(rows), sum(float(r["tempo_s"]) for r in rows.values()) / len(rows) if rows else tg)
    print(f"previsão com {P} processos (tempo por amostra medido com {P} processos simultâneos):")
    for nome, alvo in (("N = 500", 500), ("N = 1000", 1000), ("N = 2000", 2000), ("S do artigo", None)):
        cpu = sum(max((alvo if alvo else CASOS[c][1]) if alvo is None or alvo <= CASOS[c][1] else CASOS[c][1], 0)
                  * t for c, (n, t) in tc.items())
        falta = sum(max(min(alvo or CASOS[c][1], CASOS[c][1]) - n, 0) * t for c, (n, t) in tc.items())
        print(f"  {nome:12s}: {cpu / 3600:7.1f} h de CPU no total, faltam {falta / 3600:7.1f} h de CPU = "
              f"{falta / 3600 / P:6.1f} h ({falta / 86400 / P:5.2f} dias) com {P} processos")


def junta(argv):
    out = argv[0]
    for c in CASOS:
        d = os.path.join(out, c)
        if not os.path.isdir(d):
            continue
        rows, mec = ler(d)
        if not rows:
            continue
        cols = list(next(iter(rows.values())).keys())
        with open(os.path.join(out, c + ".csv"), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=cols)
            w.writeheader()
            for k in sorted(rows):
                w.writerow(rows[k])
        if mec:
            mcols = list(next(iter(mec.values())).keys())
            with open(os.path.join(out, c + ".mec"), "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=mcols)
                w.writeheader()
                for k in sorted(mec):
                    w.writerow(mec[k])
        modo = {}
        for f in glob.glob(os.path.join(d, "b*.csv.modo")):
            for ln in open(f):
                p = ln.strip().split(",")
                if len(p) == 4 and p[0].isdigit():
                    modo[int(p[0])] = ln.strip()
        if modo:
            with open(os.path.join(out, c + ".modo"), "w") as f:
                f.write("amostra,modo,y_min_banda,dist_pe\n" + "\n".join(modo[k] for k in sorted(modo)) + "\n")
        par = sorted(glob.glob(os.path.join(d, "b*.csv.param")))
        if par:
            with open(par[0]) as fi, open(os.path.join(out, c + ".param"), "w") as fo:
                fo.write(fi.read())
        print(f"{c}: {len(rows)} amostras, {len(mec)} com médias no mecanismo, {len(modo)} com modo de ruptura")


def pacote(argv):
    out = os.path.abspath(argv[0])
    arq = argv[1] if len(argv) > 1 else os.path.join(out, "resultados_campanha.tar.gz")
    n = 0
    with tarfile.open(arq, "w:gz") as tar:
        for c in CASOS:
            d = os.path.join(out, c)
            if not os.path.isdir(d):
                continue
            for f in sorted(glob.glob(os.path.join(d, "b*.csv")) + glob.glob(os.path.join(d, "b*.csv.*")) +
                            glob.glob(os.path.join(d, "preparo.log"))):
                if f.endswith((".lock", ".tmp")):
                    continue
                tar.add(f, arcname=os.path.join("campanha", c, os.path.basename(f)))
                n += 1
    print(f"{arq}: {n} arquivos ({os.path.getsize(arq) / 1e6:.1f} MB)")


if __name__ == "__main__":
    cmds = {"jobs": jobs, "status": status, "junta": junta, "pacote": pacote}
    if len(sys.argv) < 3 or sys.argv[1] not in cmds:
        print(__doc__)
        sys.exit(1)
    cmds[sys.argv[1]](sys.argv[2:])
