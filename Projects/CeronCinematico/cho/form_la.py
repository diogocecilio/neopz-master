#!/usr/bin/env python3
"""form_la.py: FORM (HL-RF) com a analise limite: g(u) = ln Gamma_LA(u), gradiente por diferencas finitas.

    python3 -I form_la.py pasta=<p> r=-0.5 [M=20] [campo=nystrom] [base=<base_nys2000.txt>] [cache=<psi.npy>]
                          [procs=4] [dif=central] [passo=0.2] [iteracoes=15] [tol_g=1e-3] [tol_beta=1e-2]

  * u = os M primeiros xi de c' e os M primeiros de phi' (os outros termos ficam em 0), no campo dado (padrao Nystrom de
    2000 modos, o mesmo das ISMC e do line sampling do MEF; as FORM do MEF M20/M100 tambem sao de Nystrom).
  * gradiente: diferencas centrais (2 x 2M minimizacoes; dif=frente: 2M + 1) com passo 'passo' em cada u_k, em paralelo.
    A busca do mecanismo tem ruido ~1e-5 em Gamma, desprezivel com passo 0,2.
  * HL-RF: u <- ((grad g . u - g)/|grad g|^2) grad g; para quando |g| < tol_g e |delta beta| < tol_beta.
  * saida: <p>/form.json (o formato do form_kl.py: parametros, beta, P_f, u, alfa, convergiu, iteracoes, avaliacoes,
    historico) e <p>/ustar.txt (os 2M numeros de u*, o centro da ISMC: mc_la.py centro=<p>/ustar.txt)."""
from __future__ import annotations

import json
import math
import multiprocessing as mp
import os
import sys
import time

import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
import campos_cho as cc  # noqa: E402
import la  # noqa: E402

_G = {}


def _gamma(u):
    av, M = _G["av"], _G["M"]
    xi = np.zeros((2, av.cc.M))
    xi[0, :M], xi[1, :M] = u[:M], u[M:]
    c, phi = av.cc.campos(xi)
    return av.minimizar(c, phi).Gamma


def cauda(b):
    return 0.5 * math.erfc(b / math.sqrt(2.0))


def main(argv):
    a = la.le_kv(argv, M="20", campo="nystrom", procs="4", dif="central", passo="0.2", iteracoes="15",
                 tol_g="1e-3", tol_beta="1e-2", h=str(la.H_PADRAO))
    if "pasta" not in a or "r" not in a:
        sys.exit(__doc__)
    p, M, procs = a["pasta"], int(a["M"]), int(a["procs"])
    passo, central = float(a["passo"]), a["dif"] == "central"
    os.makedirs(p, exist_ok=True)
    t0 = time.perf_counter()
    campos = cc.CamposCho(float(a["r"]), a["campo"], float(a["h"]), base=a.get("base"), cache=a.get("cache"))
    _G.update(av=la.Avaliador(campos), M=M)
    n = 2 * M
    u = np.zeros(n)
    hist, nav = [], 0
    beta_ant, conv = None, False
    print(f"form_la: r = {a['r']}, campo {a['campo']}, M = {M} ({n} variaveis), diferencas {a['dif']} passo {passo}, "
          f"{procs} processos", flush=True)
    with mp.get_context("fork").Pool(procs) as pool:
        for it in range(int(a["iteracoes"])):
            E = np.eye(n) * passo
            pts = [u] + [u + e for e in E] + ([u - e for e in E] if central else [])
            G = np.array(pool.map(_gamma, pts, chunksize=1))
            nav += len(pts)
            g = math.log(G[0])
            gp = np.log(G[1:n + 1])
            grad = (gp - np.log(G[n + 1:])) / (2 * passo) if central else (gp - g) / passo
            beta = float(np.linalg.norm(u))
            hist.append(dict(it=it, beta=beta, Gamma=float(G[0]), g=g, norma_grad=float(np.linalg.norm(grad))))
            print(f"  it {it}: beta = {beta:.4f}, Gamma = {G[0]:.5f}, |grad ln Gamma| = {np.linalg.norm(grad):.4f}, "
                  f"{nav} avaliacoes, {time.perf_counter() - t0:.0f} s", flush=True)
            if it > 0 and abs(g) < float(a["tol_g"]) and beta_ant is not None and \
                    abs(beta - beta_ant) < float(a["tol_beta"]):
                conv = True
                break
            unovo = ((grad @ u - g) / (grad @ grad)) * grad
            beta_ant = beta
            u = unovo
    beta = float(np.linalg.norm(u))
    alfa = (u / beta).tolist() if beta > 0 else u.tolist()
    al = np.array(alfa)
    res = {"parametros": {"M": M, "rcphi": float(a["r"]), "campo": a["campo"], "prog": "la", "dif": a["dif"],
                          "passo": passo, "nystrom": a.get("base", ""), "h": float(a["h"]), "metodo": "la"},
           "beta": beta, "P_f": cauda(beta), "convergiu": conv, "iteracoes": len(hist), "u": u.tolist(),
           "alfa": alfa, "alfa2_c": float(np.sum(al[:M] ** 2)), "alfa2_phi": float(np.sum(al[M:] ** 2)),
           "avaliacoes": nav, "historico": hist, "tempo_s": time.perf_counter() - t0}
    json.dump(res, open(os.path.join(p, "form.json"), "w"), indent=1)
    np.savetxt(os.path.join(p, "ustar.txt"), u, fmt="%.17g")
    print(f"form_la: beta = {beta:.4f}, P_f FORM = {100 * cauda(beta):.4f} %, convergiu {conv}, {nav} avaliacoes, "
          f"alfa^2: c {res['alfa2_c']:.3f}, phi {res['alfa2_phi']:.3f}; {p}/form.json", flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])
