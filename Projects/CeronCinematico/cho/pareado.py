#!/usr/bin/env python3
"""pareado.py: analise limite (LA) x elementos finitos (MEF, Cho2010) nas MESMAS realizacoes dos campos.

    python3 -I pareado.py gerar pasta=<p> r=-0.5 [semente=2010] [primeira=0] [n=20] [campo=analitico] [centro=<ustar>]
                                [base=<base_nys2000.txt>]
        grava <p>/xi/real_<n>.xi (os 2M = 4000 xi de cada realizacao, os M de c' e depois os M de phi', ja com o centro
        da ISMC se houver) e <p>/parametros.json.
    python3 -I pareado.py la pasta=<p> [procs=4] [cache=<psi.npy>] [h=<o do gerar>]
        Gamma e FS da analise limite de todas as realizacoes -> <p>/la.csv (o mesmo formato do mc_la.py), com a busca
        forte em todas (robusto=99; a comparacao valor a valor precisa de Gamma bom tambem longe de 1).
    python3 -I pareado.py mef pasta=<p> exe=<Cho2010> [procs=2] [base=<base no fermi>]
        o Cho2010 em cada xi (n=1 xi=<arq> M=2000 rcphi=r fator=tri2 nref=4 srm_f0=0.8 metodo=srm [nystrom=base]: o mesmo
        fator do MC direto e das ISMC de referencia e o comando do form_kl.comando) -> <p>/mef.csv (real, FS_MEF,
        situacao, t_s). Retoma (pula as ja "ok" no mef.csv; as com erro sao tentadas de novo). Roda com nice,
        PZ_THREADS=0, OMP_NUM_THREADS=1. Sai com codigo 3 se alguma realizacao ficar sem FS (erro do Cho2010: veja
        <p>/mef/real_<n>.log), para o job nao terminar "bem" sem o MEF.
    python3 -I pareado.py juntar pasta=<p>                         (la.csv x mef.csv da pasta)
    python3 -I pareado.py juntar la=<la.csv> mef=<pasta com mc_*.csv do MEF> [saida=<pares.csv>]
        (o MC direto ou uma ISMC do MEF ja rodados: mesma semente e mesmo centro => mesmas realizacoes)
        junta por (real) e imprime: N, correlacoes de Gamma_LA e FS_LA com FS_MEF, FS_LA - FS_MEF (media, dp, min, max),
        P_f de cada um no conjunto pareado e as falhas discordantes; grava os pares em <saida> (padrao <p>/pares.csv)."""
from __future__ import annotations

import concurrent.futures as cf
import csv
import glob
import json
import math
import os
import subprocess
import sys
import time

import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
import campos_cho as cc  # noqa: E402
import la  # noqa: E402
import mc_la  # noqa: E402

FATOR = dict(M=2000, fator="tri2", nref=4, srm_f0=0.8, metodo="srm")   # o das referencias (mc_parametros.txt)


def gerar(a):
    p = a["pasta"]
    r, sem, prim, n = float(a["r"]), int(a.get("semente", 2010)), int(a.get("primeira", 0)), int(a.get("n", 20))
    campo, centro = a.get("campo", "analitico"), a.get("centro", "")
    h = float(a.get("h", la.H_PADRAO))
    os.makedirs(os.path.join(p, "xi"), exist_ok=True)
    mu = cc.le_desloc(centro) if centro else None
    for i in range(prim, prim + n, 100):
        reais = list(range(i, min(i + 100, prim + n)))
        X = cc.xi_cpp(sem, reais, FATOR["M"])
        for real, xi in zip(reais, X):
            if mu is not None:
                xi = xi + mu
            with open(os.path.join(p, "xi", f"real_{real}.xi"), "w") as f:
                f.write("\n".join(repr(float(t)) for t in xi.ravel()) + "\n")
    par = dict(r=r, semente=sem, primeira=prim, n=n, campo=campo, centro=os.path.abspath(centro) if centro else "",
               base=os.path.abspath(a["base"]) if a.get("base") else "", h=h, **FATOR)
    json.dump(par, open(os.path.join(p, "parametros.json"), "w"), indent=1)
    print(f"gerar: {n} arquivos xi em {p}/xi (r = {r}, semente {sem}, reais {prim}..{prim + n - 1}, campo {campo}"
          f"{', centro ' + centro if centro else ''})")


def roda_la(a):
    p = a["pasta"]
    par = json.load(open(os.path.join(p, "parametros.json")))
    argv = [f"saida={os.path.join(p, 'la.csv')}", f"r={par['r']}", f"campo={par['campo']}", f"semente={par['semente']}",
            f"primeira={par['primeira']}", f"n={par['n']}", "cov=0", f"procs={a.get('procs', 4)}", "bloco=1",
            f"fs_ate={par['primeira'] + par['n']}", f"robusto={a.get('robusto', 99)}",
            f"h={a.get('h', par.get('h', la.H_PADRAO))}"]
    if par["centro"]:
        argv.append(f"centro={par['centro']}")
    if par["campo"] == "nystrom":
        argv.append(f"base={a.get('base', par['base'])}")
        if a.get("cache"):
            argv.append(f"cache={a['cache']}")
    mc_la.main(argv)


def comando(exe, par, arq_xi, rotulo, base=""):
    cmd = [exe, "n=1", "primeira=0", f"xi={arq_xi}", f"M={par['M']}", f"rcphi={par['r']}", f"fator={par['fator']}",
           f"nref={par['nref']}", f"srm_f0={par['srm_f0']}", f"metodo={par['metodo']}", f"saida={rotulo}"]
    if par["campo"] == "nystrom":
        cmd.append(f"nystrom={base or par['base']}")
    return cmd


def mef(a):
    p, exe = os.path.abspath(a["pasta"]), os.path.abspath(a["exe"])
    par = json.load(open(os.path.join(p, "parametros.json")))
    procs = int(a.get("procs", 2))
    trab = os.path.join(p, "mef")
    os.makedirs(trab, exist_ok=True)
    arq = os.path.join(p, "mef.csv")
    feitas = set()
    if os.path.exists(arq):
        mc_la.corta_parcial(arq)
        feitas = {int(l["real"]) for l in csv.DictReader(open(arq)) if l.get("situacao") == "ok"}
    reais = [i for i in range(par["primeira"], par["primeira"] + par["n"]) if i not in feitas]
    novo = not os.path.exists(arq) or os.path.getsize(arq) == 0
    f = open(arq, "a", newline="")
    if novo:
        f.write("real,FS_MEF,situacao,t_s\n")
        f.flush()
    env = dict(os.environ, PZ_THREADS="0", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
    print(f"mef: {len(reais)} realizacoes (ja feitas {len(feitas)}), {procs} processos, {exe}", flush=True)

    def uma(real):
        rot = f"real_{real}"
        t0 = time.time()
        cmd = ["nice", "-n", "10"] + comando(exe, par, os.path.join(p, "xi", rot + ".xi"), rot, a.get("base", ""))
        with open(os.path.join(trab, rot + ".log"), "w") as lg:
            subprocess.run(cmd, stdout=lg, stderr=subprocess.STDOUT, env=env, cwd=trab, check=False)
        try:
            r = next(csv.DictReader(open(os.path.join(trab, rot + "_0_0.csv"))))
            fs, sit = float(r["FS_min"]), r["situacao"]
        except (OSError, StopIteration, KeyError, ValueError):
            fs, sit = float("nan"), "erro"
        for suf in ("_mapa_0_0.csv", "_parametros.txt"):
            try:
                os.remove(os.path.join(trab, rot + suf))
            except FileNotFoundError:
                pass
        return real, fs, sit, time.time() - t0

    erros = []
    with cf.ThreadPoolExecutor(procs) as ex:
        for k, (real, fs, sit, t) in enumerate(ex.map(uma, reais)):
            f.write(f"{real},{fs!r},{sit},{t:.1f}\n")
            f.flush()
            if not math.isfinite(fs):
                erros.append(real)
            if k < 5 or not math.isfinite(fs) or (k + 1) % 24 == 0 or k + 1 == len(reais):
                print(f"  real {real}: FS_MEF = {fs:.6f} ({sit}), {t:.0f} s  [{k + 1}/{len(reais)}]", flush=True)
    f.close()
    if erros:
        lg = os.path.join(trab, f"real_{erros[0]}.log")
        cauda = open(lg).read()[-1500:] if os.path.exists(lg) else "(sem log)"
        print(f"ERRO: {len(erros)} de {len(reais)} realizacoes sem FS do MEF (reais {erros[:10]}...). Fim do log de "
              f"{lg}:\n{cauda}", file=sys.stderr, flush=True)
        sys.exit(3)


def le_mef(fonte):
    """real -> FS_MEF: de um mef.csv ou de uma pasta de rodada do Cho2010 (mc_*.csv)"""
    out = {}
    if os.path.isdir(fonte):
        for arq in glob.glob(os.path.join(fonte, "mc_[0-9]*_[0-9]*.csv")):
            for l in csv.DictReader(open(arq)):
                if int(l["real"]) >= 0 and l["situacao"] == "ok":
                    out[int(l["real"])] = (float(l["FS_min"]), float(l["log_peso"]))
    else:
        for l in csv.DictReader(open(fonte)):
            if l["situacao"] == "ok":
                out[int(l["real"])] = (float(l["FS_MEF"]), None)
    return out


def juntar(a):
    p = a.get("pasta", "")
    arq_la = a.get("la", os.path.join(p, "la.csv"))
    fonte = a.get("mef", os.path.join(p, "mef.csv"))
    saida = a.get("saida", os.path.join(p, "pares.csv") if p else os.path.splitext(arq_la)[0] + "_pares.csv")
    L = {int(l["real"]): l for l in mc_la.ler([arq_la])}
    M = le_mef(fonte)
    reais = sorted(set(L) & set(M))
    if not reais:
        sys.exit(f"nenhuma realizacao em comum entre {arq_la} ({len(L)}) e {fonte} ({len(M)})")
    G = np.array([float(L[i]["Gamma"]) for i in reais])
    F = np.array([M[i][0] for i in reais])
    fin = np.isfinite(G)                       # Gamma = inf: nenhum mecanismo admissivel (fora das correlacoes)
    fsl = np.array([float(L[i]["FS"]) if L[i].get("FS") else np.nan for i in reais])
    lw_la = np.array([float(L[i]["log_peso"] or 0) for i in reais])
    lw_m = np.array([M[i][1] if M[i][1] is not None else np.nan for i in reais])
    w = np.exp(lw_la)
    with open(saida, "w", newline="") as f:
        f.write("real,Gamma_LA,FS_LA,FS_MEF,modo,log_peso\n")
        for i, g, fl, fm, lw in zip(reais, G, fsl, F, lw_la):
            f.write(f"{i},{g:.6f},{'' if np.isnan(fl) else f'{fl:.6f}'},{fm:.6f},{L[i]['modo']},{lw!r}\n")
    N = len(reais)
    print(f"pareado: {arq_la} x {fonte}: N = {N} realizacoes em comum")
    if np.any(np.isfinite(lw_m)):
        d = np.nanmax(np.abs(lw_m - lw_la))
        print(f"  log_peso LA x MEF: dif max {d:.1e} (as mesmas realizacoes: {'OK' if d < 1e-5 else 'NAO!'})")
    Gf, Ff = G[fin], F[fin]
    print(f"  corr(Gamma_LA, FS_MEF) = {np.corrcoef(Gf, Ff)[0, 1]:.4f}; corr(ln Gamma_LA, ln FS_MEF) = "
          f"{np.corrcoef(np.log(Gf), np.log(Ff))[0, 1]:.4f}; Spearman = {_spearman(G, F):.4f}"
          + (f" ({int((~fin).sum())} com Gamma = inf fora)" if (~fin).any() else ""))
    q = Gf / Ff
    print(f"  Gamma_LA / FS_MEF: media {q.mean():.4f}, mediana {np.median(q):.4f}; perto da falha (FS_MEF < 1,1, "
          f"{int(np.sum(Ff < 1.1))}): media {q[Ff < 1.1].mean() if np.any(Ff < 1.1) else float('nan'):.4f}")
    ok = np.isfinite(fsl)
    if ok.sum() > 2:
        d = fsl[ok] - F[ok]
        print(f"  FS_LA x FS_MEF ({ok.sum()} com FS): corr {np.corrcoef(fsl[ok], F[ok])[0, 1]:.4f}; FS_LA/FS_MEF media "
              f"{np.mean(fsl[ok] / F[ok]):.4f}; FS_LA - FS_MEF: media "
              f"{d.mean():+.4f}, dp {d.std(ddof=1):.4f}, min {d.min():+.4f}, max {d.max():+.4f}; FS_LA/FS_MEF - 1: media "
              f"{100 * np.mean(fsl[ok] / F[ok] - 1):+.2f} %; FS_LA < FS_MEF em {int(np.sum(d < 0))}")
    fl, fm = G < 1, F < 1
    pl, pm = np.mean(fl * w), np.mean(fm * w)
    print(f"  falhas: LA {int(fl.sum())}, MEF {int(fm.sum())}, as duas {int((fl & fm).sum())}, so LA {int((fl & ~fm).sum())},"
          f" so MEF {int((~fl & fm).sum())}")
    print(f"  P_f no conjunto pareado{' (com o peso da ISMC)' if np.any(lw_la != 0) else ''}: LA {100 * pl:.4f} %, MEF "
          f"{100 * pm:.4f} %, razao LA/MEF {pl / pm if pm > 0 else float('nan'):.3f}")
    if (fl | fm).any():
        dif = fm.astype(float) - fl.astype(float)
        sd = np.std(dif * w, ddof=1) / math.sqrt(N)
        print(f"  P_f(MEF) - P_f(LA) = {100 * (pm - pl):.4f} % +- {100 * 1.96 * sd:.4f} % (IC 95 % pareado)")
    print(f"  pares em {saida}")


def _spearman(a, b):
    ra, rb = np.argsort(np.argsort(a)), np.argsort(np.argsort(b))
    return float(np.corrcoef(ra, rb)[0, 1])


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    f = {"gerar": gerar, "la": roda_la, "mef": mef, "juntar": juntar}.get(cmd)
    if f is None:
        sys.exit(__doc__)
    f(la.le_kv(sys.argv[2:]))
