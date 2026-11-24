#!/usr/bin/env python3
"""lote_ref.py: o Monte Carlo de referencia em Python (coluna_ref.py mc) de varios casos em paralelo, em blocos de
realizacoes intercalados entre os casos (todos andam juntos; um bloco pronto nao e refeito: da para parar e retomar).

    python3 -I lote_ref.py <casos.txt> <pasta> [processos=8] [bloco=100] [M=100] [opcoes do coluna_ref.py ...]

casos.txt: 'id Pf18 Pf24 N opcoes...' por linha (# comenta; Pf18 e Pf24 sao os valores do artigo, '-' se nao ha), o
formato do casos_cho2014.txt da spec (D.1). Cada caso precisa de lv=; as opcoes do caso vem depois das gerais.
A base de Nystrom de cada lv e <pasta>/bases/base_lv<lv>_M<M>.txt (Projects/Cho2010/nystrom_base.py coluna=H lv=lv
M=M: n = 8M nos, dy = 0.0025 com H = 2 e M = 100), feita aqui se falta. As realizacoes do caso vao em
<pasta>/<id>/<id>_<a>_<b>.csv (as colunas do SlopeRainColumn); o andamento, em <pasta>/status.txt (reescrito a cada
bloco) e <pasta>/lote.log. Com a mesma semente (padrao 2025), os casos ficam pareados (os mesmos xi)."""
import multiprocessing as mp
import os
import subprocess
import sys
import time

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)

import coluna_ref as R  # noqa: E402

NYSTROM = os.path.join(AQUI, "..", "..", "Cho2010", "nystrom_base.py")


def le_casos(arq):
    casos = []
    for linha in open(arq):
        t = linha.split()
        if not t or t[0].startswith("#"):
            continue
        if len(t) < 4:
            sys.exit(f"{arq}: linha curta: {linha.strip()}")
        op = dict(x.split("=", 1) for x in t[4:])
        if "lv" not in op:
            sys.exit(f"{arq}: o caso {t[0]} precisa de lv=")
        casos.append({"id": t[0], "Pf18": t[1], "Pf24": t[2], "N": int(t[3]), "opcoes": t[4:], "lv": op["lv"],
                      "H": op.get("H", "2")})
    return casos


def base_de(pasta, lv, M, H):
    arq = os.path.join(pasta, "bases", f"base_lv{lv}_M{M}.txt")
    if not os.path.exists(arq):
        os.makedirs(os.path.dirname(arq), exist_ok=True)
        with open(arq[:-4] + ".log", "w") as log:
            subprocess.run([sys.executable, "-I", NYSTROM, arq + ".tmp", f"coluna={H}", f"lv={lv}", f"M={M}",
                            "mercer=2000"], stdout=log, stderr=subprocess.STDOUT, check=True)
        os.replace(arq + ".tmp", arq)
    return arq


_CACHE = {}


def tarefa(t):
    """um bloco de realizacoes de um caso (no processo filho)"""
    caso, a, n, args = t
    if caso not in _CACHE:
        o, _ = R.le_opcoes(args, "mc")
        _CACHE.clear()  # um caso por vez na memoria de cada processo
        _CACHE[caso] = (o, R.prepara(o), R.Base(o["nystrom"], o["M"]))
    o, slz, base = _CACHE[caso]
    t0 = time.time()
    try:
        r = R.mc_bloco(o, a, n, so_lei_zs=slz, base=base, quieto=True)
    except Exception as e:  # o bloco fica sem arquivo (refeito na proxima rodada); o lote segue
        return caso, a, n, time.time() - t0, 0., {"erro": f"{type(e).__name__}: {e}"}
    return caso, a, n, time.time() - t0, r["s"], r["falhas"]


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    casos = le_casos(sys.argv[1])
    pasta = os.path.abspath(sys.argv[2])
    g = {"processos": "8", "bloco": "100", "M": "100"}
    gerais = []
    for x in sys.argv[3:]:
        k, v = x.split("=", 1)
        if k in ("processos", "bloco"):
            g[k] = v
        else:
            if k == "M":
                g["M"] = v
            gerais.append(x)
    P, bloco, M = int(g["processos"]), int(g["bloco"]), int(g["M"])
    os.makedirs(pasta, exist_ok=True)
    tarefas, total, feitas0 = [], {}, {}
    nb = max((c["N"] + bloco - 1) // bloco for c in casos)
    for c in casos:
        base = base_de(pasta, c["lv"], M, c["H"])
        os.makedirs(os.path.join(pasta, c["id"]), exist_ok=True)
        c["args"] = gerais + c["opcoes"] + [f"nystrom={base}", f"saida={os.path.join(pasta, c['id'], c['id'])}"]
        R.le_opcoes(c["args"], "mc")  # confere as opcoes antes de lancar
        total[c["id"]], feitas0[c["id"]] = c["N"], 0
    for b in range(nb):  # intercalados: o bloco b de todos os casos, depois o b + 1
        for c in casos:
            a = b * bloco
            if a >= c["N"]:
                continue
            n = min(bloco, c["N"] - a)
            if os.path.exists(os.path.join(pasta, c["id"], f"{c['id']}_{a}_{a + n - 1}.csv")):
                feitas0[c["id"]] += n
                continue
            tarefas.append((c["id"], a, n, c["args"]))
    feitas = dict(feitas0)
    t0, s_cpu, n_novas = time.time(), 0., 0
    log = open(os.path.join(pasta, "lote.log"), "a")
    log.write(f"{time.strftime('%F %T')} inicio: {len(casos)} casos, {len(tarefas)} blocos a fazer ({sum(t[2] for t in tarefas)}"
              f" realizacoes), {P} processos, bloco {bloco}; {sum(feitas0.values())} realizacoes ja feitas\n")
    log.flush()

    def status(fim=False):
        el = time.time() - t0
        faltam = sum(total.values()) - sum(feitas.values())
        taxa = el / n_novas if n_novas else float("nan")
        with open(os.path.join(pasta, "status.txt.tmp"), "w") as s:
            s.write(f"{'FIM' if fim else 'rodando'} {time.strftime('%F %T')} (inicio {time.strftime('%F %T', time.localtime(t0))},"
                    f" pid {os.getpid()}, {P} processos)\n")
            s.write(f"realizacoes: {sum(feitas.values())} de {sum(total.values())}; {el / 60:.1f} min; "
                    f"{s_cpu / max(n_novas, 1):.2f} s por realizacao por processo")
            if not fim and n_novas:
                s.write(f"; fim previsto {time.strftime('%F %T', time.localtime(time.time() + faltam * taxa))}")
            s.write("\n")
            for c in casos:
                s.write(f"  {c['id']:22s} {feitas[c['id']]:6d} / {total[c['id']]}\n")
        os.replace(os.path.join(pasta, "status.txt.tmp"), os.path.join(pasta, "status.txt"))

    status()
    with mp.get_context("fork").Pool(P) as pool:
        for caso, a, n, s, sc, falhas in pool.imap_unordered(tarefa, tarefas):
            if "erro" in falhas:
                log.write(f"{time.strftime('%F %T')} {caso} {a}..{a + n - 1}: ERRO {falhas['erro']}\n")
                log.flush()
                continue
            feitas[caso] += n
            n_novas += n
            s_cpu += sc
            log.write(f"{time.strftime('%F %T')} {caso} {a}..{a + n - 1}: {s:.0f} s; FS < 1: " +
                      " ".join(f"{t:g}h:{k}" for t, k in sorted(falhas.items())) + "\n")
            log.flush()
            status()
    status(fim=True)
    log.write(f"{time.strftime('%F %T')} fim\n")
    log.close()


if __name__ == "__main__":
    main()
