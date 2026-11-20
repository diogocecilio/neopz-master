#!/usr/bin/env python3
"""ls_la.py: line sampling com a analise limite, na direcao de um FORM (do MEF ou do form_la.py).

    python3 -I ls_la.py direcao=<form.json> pasta=<p> r=-0.5 [n=100] [semente=7] [K=2000] [campo=nystrom]
                        [base=<base_nys2000.txt>] [cache=<psi.npy>] [procs=4] [tol=1e-3] [cov=0] [nmin=30]

  * a mesma construcao do Projects/Cho2010/line_sampling.py: e = u*/|u*| (os M do FORM nos K termos de cada campo,
    completados com zeros), z da linha i = o i-esimo vetor de 2K normais de random.Random(semente).gauss (os MESMOS z
    do line sampling do MEF com a mesma semente e o mesmo K), zp = z - (z.e) e, c_i: Gamma(zp + c e) = 1, p_i = Phi(-c_i),
    P_f = media de p_i. Com a direcao e a semente do MEF, a linha i daqui e a linha i do MEF: c_i(LA) x c_i(MEF).
  * raiz de g(c) = ln Gamma: secante protegida (Illinois) a partir de c = beta da direcao e do passo de Newton com a
    inclinacao medida no ponto de projeto (como no MEF); para quando o passo < tol. Sem raiz em [-8, 12]: p = 0 ou 1.
  * parada: n linhas ou, com cov > 0, o CoV(P_f) <= cov depois de nmin linhas.
  * linha perdida: 25 analises sem cercar a raiz (g nao muda de sinal). Vai para <p>/perdidas.csv, fica fora da P_f
    (como no line_sampling.py do MEF) e o ls.json conta as perdidas; nao e refeita na retomada.
  * saidas: <p>/linhas.csv (linha, c, analises, p), <p>/perdidas.csv (linha, analises, ultimo_c), <p>/analises.csv
    (linha, c, Gamma), <p>/inclinacao.json, <p>/ls.json (o resultado). Retoma pelas linhas ja no linhas.csv ou no
    perdidas.csv (uma ultima linha cortada e apagada antes)."""
from __future__ import annotations

import csv
import json
import math
import multiprocessing as mp
import os
import random
import sys
import time

import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
import campos_cho as cc  # noqa: E402
import la  # noqa: E402
import mc_la  # noqa: E402

FIM = (-8.0, 12.0)
_G = {}


def cauda(c):
    return 0.5 * math.erfc(c / math.sqrt(2.0))


def _gamma(x):
    av = _G["av"]
    c, phi = av.cc.campos(np.asarray(x).reshape(2, -1))
    return min(av.minimizar(c, phi).Gamma, 1e6)      # inf = nenhum mecanismo admissivel (muito seguro)


def _linha(args):
    """c_i da linha i: g(c) = ln Gamma(zp + c e) = 0 (g > 0 seguro, g <= 0 falha, a falha do lado de c grande)"""
    i, zp = args
    e, beta, s, tol = _G["e"], _G["beta"], _G["s"], _G["tol"]
    reg = []

    def g(c):
        G = _gamma(zp + c * e)
        reg.append((c, G))
        return math.log(G)

    s = s if s < 0 else -0.3
    pts = [(beta, g(beta))]
    cn = beta - pts[0][1] / s
    # 1) cercar a raiz: passos de secante (o 1o com a inclinacao s do ponto de projeto), no maximo 3 de cada vez
    while True:
        segs = [q for q in pts if q[1] > 0]
        fals = [q for q in pts if q[1] <= 0]
        if segs and fals:
            break
        if len(reg) >= 25:
            return i, float("nan"), len(reg), float("nan"), reg     # sem cercar a raiz: linha perdida
        cn = min(max(cn, FIM[0]), FIM[1])
        gn = g(cn)
        pts.append((cn, gn))
        if gn > 0 and cn >= FIM[1]:
            return i, FIM[1], len(reg), 0.0, reg
        if gn <= 0 and cn <= FIM[0]:
            return i, FIM[0], len(reg), 1.0, reg
        (c1, g1), (c2, g2) = pts[-2], pts[-1]
        sl = (g2 - g1) / (c2 - c1) if abs(c2 - c1) > 1e-12 else s
        if not sl < 0:
            sl = s
        cn = c2 + max(-3.0, min(3.0, -g2 / sl))
    # 2) Illinois entre o maior c seguro e o menor c na falha
    ac, ag = max(segs)
    bc, bg = min(fals)
    lado = 0
    c = ac - ag * (bc - ac) / (bg - ag)
    while bc - ac > tol and len(reg) < 40:
        c = ac - ag * (bc - ac) / (bg - ag)
        gc = g(c)
        if abs(gc) < tol * abs(s):
            break
        if gc > 0:
            ac, ag = c, gc
            if lado == 1:
                bg *= 0.5
            lado = 1
        else:
            bc, bg = c, gc
            if lado == -1:
                ag *= 0.5
            lado = -1
    return i, c, len(reg), cauda(c), reg


def main(argv):
    a = la.le_kv(argv, n="100", semente="7", K="2000", campo="nystrom", procs="4", tol="1e-3", cov="0", nmin="30",
                 h=str(la.H_PADRAO))
    if "direcao" not in a or "pasta" not in a or "r" not in a:
        sys.exit(__doc__)
    p, K, nmax, procs = a["pasta"], int(a["K"]), int(a["n"]), int(a["procs"])
    alvo, nmin = float(a["cov"]), int(a["nmin"])
    os.makedirs(p, exist_ok=True)
    form = json.load(open(a["direcao"]))
    u = np.array(form["u"], float)
    M = len(u) // 2
    e = np.zeros(2 * K)
    e[:M], e[K:K + M] = u[:M], u[M:]
    beta = float(np.linalg.norm(e))
    e /= beta
    t0 = time.perf_counter()
    campos = cc.CamposCho(float(a["r"]), a["campo"], float(a["h"]), M=K, base=a.get("base"), cache=a.get("cache"))
    _G.update(av=la.Avaliador(campos), e=e, beta=beta, tol=float(a["tol"]))
    arq_s = os.path.join(p, "inclinacao.json")
    if os.path.exists(arq_s):
        s = json.load(open(arq_s))["s"]
        n_s = 0
    else:
        gp, gm = math.log(_gamma((beta + 0.1) * e)), math.log(_gamma((beta - 0.1) * e))
        s = (gp - gm) / 0.2
        n_s = 2
        json.dump({"s": s, "beta": beta, "Gamma_p": math.exp(gp), "Gamma_m": math.exp(gm)}, open(arq_s, "w"))
    _G["s"] = s
    print(f"ls_la: direcao {a['direcao']} (M = {M}), beta = {beta:.4f}, P_f FORM = {100 * cauda(beta):.4f} %, "
          f"d ln Gamma/dc = {s:.4f}; campo {a['campo']} K = {K}, r = {a['r']}, semente {a['semente']}", flush=True)
    arq_l, arq_p, arq_a = (os.path.join(p, x) for x in ("linhas.csv", "perdidas.csv", "analises.csv"))
    for arq in (arq_l, arq_p, arq_a):
        mc_la.corta_parcial(arq)
    feitas, perdidas = {}, {}
    if os.path.exists(arq_l):
        for l in csv.DictReader(open(arq_l)):
            try:
                float(l["c"]), int(l["analises"])
                feitas[int(l["linha"])] = l
            except (TypeError, ValueError, KeyError):
                pass
    if os.path.exists(arq_p):
        for l in csv.DictReader(open(arq_p)):
            try:
                perdidas[int(l["linha"])] = l
            except (TypeError, ValueError, KeyError):
                pass
    vazio = lambda arq: not os.path.exists(arq) or os.path.getsize(arq) == 0  # noqa: E731
    novo_l, novo_p = vazio(arq_l), vazio(arq_p)
    fl = open(arq_l, "a", newline="")
    if novo_l:
        fl.write("linha,c,analises,p\n")
    fp = open(arq_p, "a", newline="")
    if novo_p:
        fp.write("linha,analises,ultimo_c\n")
    fa = open(arq_a, "a")
    rng = random.Random(int(a["semente"]))
    dim = 2 * K

    def zs():
        for i in range(nmax):
            z = np.array([rng.gauss(0.0, 1.0) for _ in range(dim)])
            if i in feitas or i in perdidas:
                continue
            yield i, z - (z @ e) * e

    res = {i: (float(l["c"]), int(l["analises"])) for i, l in feitas.items()}

    def estat():
        ps = np.array([cauda(c) for c, _ in res.values()])
        N = len(ps)
        m = ps.mean() if N else float("nan")
        sd = ps.std(ddof=1) if N > 1 else float("nan")
        return N, m, sd / (math.sqrt(N) * m) if N > 1 and m > 0 else float("inf"), sd

    motivo = ""
    with mp.get_context("fork").Pool(procs) as pool:
        it = zs()
        voo = []
        acabou = False
        while True:
            while not motivo and not acabou and len(voo) < 2 * procs:
                try:
                    voo.append(pool.apply_async(_linha, (next(it),)))
                except StopIteration:
                    acabou = True
            if not voo:
                break
            pr = [v for v in voo if v.ready()]
            if not pr:
                time.sleep(0.05)
                continue
            for v in pr:
                voo.remove(v)
                i, c, na, pi, reg = v.get()
                if not math.isfinite(c):
                    perdidas[i] = {"analises": na}
                    fp.write(f"{i},{na},{reg[-1][0]!r}\n")
                    fp.flush()
                    print(f"  linha {i} perdida: {na} analises sem cercar a raiz", flush=True)
                else:
                    res[i] = (c, na)
                    fl.write(f"{i},{c!r},{na},{pi!r}\n")
                for cc_, G in reg:
                    fa.write(f"{i},{cc_!r},{G!r}\n")
            fl.flush()
            fa.flush()
            N, m, cov, _ = estat()
            if N % 10 == 0:
                print(f"  {N} linhas: P_f = {100 * m:.4f} %, CoV = {100 * cov:.2f} %", flush=True)
            if alvo > 0 and N >= nmin and cov <= alvo and not motivo:
                motivo = f"CoV {100 * cov:.2f} % <= {100 * alvo:g} %"
    N, m, cov, sd = estat()
    na = [x[1] for x in res.values()]
    out = {"direcao": os.path.abspath(a["direcao"]), "beta": beta, "s": s, "P_f": m, "CoV": cov,
           "CoV_linha": sd / m if m > 0 else float("nan"), "linhas": N, "analises_linha": float(np.mean(na)),
           "analises": int(sum(na)) + n_s + sum(int(x["analises"]) for x in perdidas.values()), "analises_form": form.get("avaliacoes"),
           "analises_cov5": (form.get("avaliacoes") or 0) + 2 + float(np.mean(na)) * N * (cov / 0.05) ** 2,
           "perdidas": len(perdidas), "linhas_perdidas": sorted(perdidas),
           "fracao_perdida": len(perdidas) / max(1, N + len(perdidas)), "motivo": motivo or "n linhas", "r": float(a["r"]), "K": K, "semente": int(a["semente"]),
           "campo": a["campo"], "tempo_s": time.perf_counter() - t0}
    json.dump(out, open(os.path.join(p, "ls.json"), "w"), indent=1)
    print(f"ls_la: {N} linhas, P_f = {100 * m:.4f} %, CoV = {100 * cov:.2f} %, CoV por linha {100 * out['CoV_linha']:.0f} %,"
          f" {out['analises_linha']:.2f} analises/linha, {len(perdidas)} perdidas, {out['tempo_s']:.0f} s; {p}/ls.json",
          flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])
