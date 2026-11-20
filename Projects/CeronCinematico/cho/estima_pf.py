#!/usr/bin/env python3
"""estima_pf.py: previsao da P_f da analise limite (LA) e do custo do MC direto LA ate CoV 5 %, pelo piloto pareado.

    python3 -I estima_pf.py la=<piloto la.csv com FS> mef=<pasta do MC direto do MEF> [t=<s por realizacao>]

O piloto (mc_la.py com fs_ate=N, semente 2010) da FS_LA e FS_MEF nas mesmas realizacoes. Ajusta
ln FS_LA = a + b ln FS_MEF + e, e ~ N(0, s^2), e aplica a TODAS as realizacoes do MC direto do MEF:
P_f(LA) ~ media de Phi(-(a + b ln FS_MEF)/s). Como Gamma < 1 <=> FS_LA < 1, e a P_f da LA. N para CoV 5 % = (1 - p)/(p 0,05^2).
E uma extrapolacao (o ajuste vem de FS ~ 1,0-1,5, a cauda pesa): serve so para o tamanho e o custo do MC LA."""
from __future__ import annotations

import math
import os
import sys

import numpy as np
from scipy.stats import norm

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
import la  # noqa: E402
import mc_la  # noqa: E402
import pareado  # noqa: E402


def main(argv):
    a = la.le_kv(argv)
    L = {int(l["real"]): l for l in mc_la.ler([a["la"]]) if l.get("FS")}
    M = pareado.le_mef(a["mef"])
    comum = sorted(set(L) & set(M))
    x = np.log([M[i][0] for i in comum])
    y = np.log([float(L[i]["FS"]) for i in comum])
    A = np.c_[np.ones_like(x), x]
    (c0, c1), *_ = np.linalg.lstsq(A, y, rcond=None)
    s = float(np.std(y - A @ [c0, c1], ddof=2))
    todos = np.log(np.array([v[0] for v in M.values()]))
    p_la = float(np.mean(norm.cdf(-(c0 + c1 * todos) / s)))
    p_mef = float(np.mean(todos < 0))
    n_la = (1 - p_la) / (p_la * 0.0025)
    print(f"piloto pareado: {len(comum)} realizacoes; ln FS_LA = {c0:+.4f} + {c1:.4f} ln FS_MEF, residuo dp {s:.4f}; "
          f"FS_LA - FS_MEF media {np.mean(np.exp(y) - np.exp(x)):+.4f}")
    print(f"MC direto do MEF ({a['mef']}): N = {len(todos)}, P_f(MEF) = {100 * p_mef:.3f} %")
    print(f"P_f(LA) prevista = {100 * p_la:.3f} % (razao LA/MEF {p_la / p_mef:.2f}); N para CoV 5 % = {n_la:.0f}")
    if a.get("t"):
        t = float(a["t"])
        print(f"custo com {t:.2f} s por realizacao: {n_la * t / 3600:.1f} h de CPU; PC com 4 processos "
              f"{n_la * t / 4 / 3600:.1f} h")


if __name__ == "__main__":
    main(sys.argv[1:])
