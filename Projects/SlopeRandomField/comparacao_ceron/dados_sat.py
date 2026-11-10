#!/usr/bin/env python3
"""dados_sat.py: le os resultados das tres fases do exemplo de Vargas Ceron et al. (2025, IJNAMG 49:2459, Secao 6) no
limite saturado (SATURADO=1 SUPERFICIE_P0=1 KW=1e30, README "As tres fases no fermi") e calcula, por caso, as grandezas
das Tabelas 3 a 6 e das Figs. 25 a 33 do artigo, do mesmo jeito que ele:

- S (realizacoes usadas), n_f = #(Gamma < 1), P_f = n_f / S e CoV(P_f) = sqrt((1 - P_f) / (S P_f)) (Eqs. 67 a 69);
- o ajuste log-normal (maxima verossimilhanca): media, mediana, desvio e CoV da distribuicao ajustada, que sao os mu,
  sigma e CoV das tabelas do Ceron ("based on the fit of the PDF", Secao 6); as estatisticas amostrais ao lado;
- o teste qui-quadrado dos ajustes normal e log-normal (k = ceil(2 S^0,4) classes equiprovaveis, gl = k - 3);
- a convergencia de P_f com o numero de simulacoes, na ordem das realizacoes (Fig. 25b);
- a correlacao de Spearman parcial de 2a ordem entre Gamma e (c', phi', k) de cada realizacao (Eqs. 70 e 71, Fig. 26),
  com duas definicoes das variaveis: as medias na banda de ruptura e na massa que desliza, e as medias do campo;
- o Gamma deterministico (CAMPOS=0) de cada caso que o tem e a banda de ruptura dele (Fig. 24).

O Gamma de cada realizacao e o FS_min de mc_resumo.csv (com FATOR=gim_arco, o fator do GIM em t = 0,01 h). Uma
realizacao entra se tem situacao ok e Gamma > 0 finito; as outras sao contadas e listadas.

Estrutura lida (pacote do README): <raiz>/<prefixo>_fase<N>/<caso>/mc_resumo.csv e parametros.txt, e
<raiz>/<prefixo>_fase<N>/<caso>_det/r*/ (resumo.csv, fs.csv, mecanismo.csv). Os dados so sao lidos.
Uso: python3 -I dados_sat.py <raiz> <prefixo, por exemplo ceron_sat> <saida.json>"""

import csv
import glob
import json
import math
import os
import sys

import numpy as np
from scipy import stats

# (familia, valor, fase, sufixo da pasta): o valor da referencia de cada familia e VALOR_REF
CASOS = [
    ("ref", None, 1, ""),
    ("covk", 0.0, 2, "_covk0"), ("covk", 0.75, 2, "_covk0.75"), ("covk", 0.9, 2, "_covk0.9"), ("covk", 1.0, 2, "_covk1"),
    ("covc", 0.1, 2, "_covc0.1"), ("covc", 0.5, 2, "_covc0.5"), ("covc", 0.7, 2, "_covc0.7"),
    ("covphi", 0.05, 2, "_covphi0.05"), ("covphi", 0.15, 2, "_covphi0.15"), ("covphi", 0.2, 2, "_covphi0.2"),
    ("s", 1.5, 2, "_lx30_ly3"), ("s", 2.0, 2, "_lx40_ly4"), ("s", 5.0, 2, "_lx100_ly10"), ("s", 10.0, 2, "_lx200_ly20"),
    ("s", 20.0, 2, "_lx400_ly40"), ("s", 400.0, 2, "_lx8000_ly800"),
    ("alfa", 2.0, 2, "_alfa2"), ("alfa", 3.0, 2, "_alfa3"), ("alfa", 4.0, 2, "_alfa4"), ("alfa", 5.0, 2, "_alfa5"),
    ("hw", 0.5, 3, "_nivelfinal7.5"), ("hw", 0.6, 3, "_nivelfinal7"), ("hw", 0.7, 3, "_nivelfinal6.5"),
    ("hw", 0.8, 3, "_nivelfinal6"), ("hw", 0.9, 3, "_nivelfinal5.5"),
]
VALOR_REF = {"covk": 0.6, "covc": 0.3, "covphi": 0.1, "s": 1.0, "alfa": 1.0, "hw": 1.0}
# familias cujo deterministico muda com o valor (rodar_mc.sh: ALFA e NIVEL_FINAL entram no nome do deterministico)
DET_PROPRIO = ("alfa", "hw")
BASE = "f50_I0_iii"
# as variaveis da correlacao de Spearman: (nome, coluna da banda/massa, coluna da media do campo)
VARS_SPEARMAN = (("c", "banda_c", "c_med"), ("phi", "banda_phi_deg", "phi_med"), ("k", "massa_k_mult", "k_mult_med"))


def _f(v):
    try:
        x = float(v)
    except (TypeError, ValueError):
        return math.nan
    return x if math.isfinite(x) else math.nan


def le_csv(arq):
    with open(arq, newline="") as f:
        return list(csv.DictReader(f))


def le_parametros(arq):
    p = {}
    if os.path.isfile(arq):
        for linha in open(arq):
            linha = linha.strip()
            if "=" in linha:
                k, v = linha.split("=", 1)
                p[k] = v
    return p


def postos(x):
    return stats.rankdata(x)


def parcial(C, i, j):
    """correlacao parcial de i e j controlando as outras variaveis de C; com 4 variaveis, a de 2a ordem (Eq. 71)"""
    P = np.linalg.inv(C)
    return float(-P[i, j] / math.sqrt(P[i, i] * P[j, j]))


def ic95(r, n, k):
    if n - 3 - k <= 0 or not math.isfinite(r) or abs(r) >= 1:
        return (math.nan, math.nan)
    z, s = math.atanh(r), 1. / math.sqrt(n - 3 - k)
    return (math.tanh(z - 1.96 * s), math.tanh(z + 1.96 * s))


def qui2(g, dist, k):
    """p-valor do qui-quadrado com k classes equiprovaveis da distribuicao ajustada (gl = k - 3)"""
    bordas = dist.ppf(np.linspace(0., 1., k + 1))
    obs = np.histogram(g, bins=np.concatenate([[-np.inf], bordas[1:-1], [np.inf]]))[0]
    esp = len(g) / k
    x2 = float(((obs - esp) ** 2 / esp).sum())
    gl = k - 3
    return x2, gl, float(stats.chi2.sf(x2, gl))


def estatisticas(g):
    g = np.asarray(g, float)
    S = len(g)
    nf = int((g < 1.).sum())
    pf = nf / S
    covpf = math.sqrt((1. - pf) / (S * pf)) if pf > 0 else math.inf
    ln = np.log(g)
    ml, sl = float(ln.mean()), float(ln.std(ddof=0))  # maxima verossimilhanca da log-normal
    cov = math.sqrt(math.exp(sl ** 2) - 1.)
    media = math.exp(ml + sl ** 2 / 2.)
    k = int(math.ceil(2. * S ** 0.4))
    x2n = qui2(g, stats.norm(g.mean(), g.std(ddof=0)), k)
    x2l = qui2(g, stats.lognorm(s=sl, scale=math.exp(ml)), k)
    return {
        "S": S, "n_f": nf, "P_f": pf, "CoV_Pf": covpf,
        "mu_ln": ml, "sigma_ln": sl, "media": media, "mediana": math.exp(ml), "desvio": media * cov, "CoV": cov,
        "media_amostral": float(g.mean()), "desvio_amostral": float(g.std(ddof=1)),
        "CoV_amostral": float(g.std(ddof=1) / g.mean()), "minimo": float(g.min()), "maximo": float(g.max()),
        "qui2_classes": k, "qui2_normal": {"x2": x2n[0], "gl": x2n[1], "p": x2n[2]},
        "qui2_lognormal": {"x2": x2l[0], "gl": x2l[1], "p": x2l[2]},
    }


def spearman(linhas, g, qual):
    """Spearman parcial de 2a ordem entre Gamma e cada variavel, controlando as outras duas"""
    col = 1 if qual == "banda" else 2
    X = [g] + [np.array([_f(r[v[col]]) for r in linhas]) for v in VARS_SPEARMAN]
    D = np.column_stack(X)
    ok = np.all(np.isfinite(D), axis=1)
    D = D[ok]
    n = len(D)
    if n < 10 or np.any(np.std(D, axis=0) == 0):
        return None  # poucas realizacoes ou uma variavel constante
    C = np.corrcoef(np.column_stack([postos(D[:, j]) for j in range(D.shape[1])]), rowvar=False)
    res = {"n": n}
    for j, v in enumerate(VARS_SPEARMAN, start=1):
        r = parcial(C, 0, j)
        res[v[0]] = {"rho_parcial": r, "ic95": ic95(r, n, 2), "rho_simples": float(C[0, j])}
    return res


def _media_finita(v):
    v = [x for x in v if math.isfinite(x)]
    return float(np.mean(v)) if v else math.nan


def deterministico(pasta):
    arqs = sorted(glob.glob(os.path.join(pasta, "r*", "resumo.csv")))
    if not arqs:
        return None
    d = os.path.dirname(arqs[0])
    r = le_csv(arqs[0])[0]
    det = {"pasta": os.path.relpath(d, os.path.dirname(os.path.dirname(pasta))), "Gamma": _f(r.get("FS_min")),
           "situacao": r.get("situacao"), "versao": r.get("versao")}
    for c in ("banda_area_m2", "banda_xc", "banda_yc", "banda_ymin"):
        det[c] = _f(r.get(c))
    mec = os.path.join(d, "mecanismo.csv")
    if os.path.isfile(mec):
        pts = [(_f(m["x"]), _f(m["y"]), _f(m["epsp_raizJ2"])) for m in le_csv(mec)]
        if pts:
            emax = max(p[2] for p in pts)
            det["banda_pontos"] = [(x, y) for x, y, e in pts if e >= 0.1 * emax]  # BANDA_LIMIAR padrao
    return det


def calcula(raiz, prefixo):
    casos, avisos, versoes = [], [], {}
    for fam, valor, fase, suf in CASOS:
        pasta = os.path.join(raiz, f"{prefixo}_fase{fase}", BASE + suf)
        arq = os.path.join(pasta, "mc_resumo.csv")
        cid = fam if fam == "ref" else f"{fam}={valor:g}"
        if not os.path.isfile(arq):
            avisos.append(f"{cid}: sem {os.path.relpath(arq, raiz)}")
            continue
        linhas = le_csv(arq)
        linhas.sort(key=lambda r: int(_f(r.get("real")) if math.isfinite(_f(r.get("real"))) else 1e9))
        bons, fora = [], []
        for r in linhas:
            g = _f(r.get("FS_min"))
            if r.get("situacao") == "ok" and math.isfinite(g) and g > 0:
                bons.append(r)
            else:
                fora.append({"real": r.get("real"), "situacao": r.get("situacao"), "Gamma": r.get("FS_min")})
        for r in linhas:
            versoes[r.get("versao", "?")] = versoes.get(r.get("versao", "?"), 0) + 1
        g = np.array([_f(r["FS_min"]) for r in bons])
        caso = {"id": cid, "familia": fam, "valor": VALOR_REF.get(fam) if fam == "ref" else valor,
                "pasta": os.path.relpath(pasta, raiz), "parametros": le_parametros(os.path.join(pasta, "parametros.txt")),
                "realizacoes": len(linhas), "excluidas": fora,
                "KL_variancia_media": _media_finita([_f(r.get("KL_variancia")) for r in bons])}
        # os parametros dos campos, conferidos em cada linha
        for col in ("cov_c", "cov_phi", "cov_k", "Lx", "Ly", "alfa"):
            vals = sorted({r.get(col) for r in linhas})
            caso[col] = vals[0] if len(vals) == 1 else vals
        if len(g):
            caso["estat"] = estatisticas(g)
            caso["gamma"] = g.tolist()
            caso["conv_pf"] = (np.cumsum(g < 1.) / np.arange(1, len(g) + 1)).tolist()
            if fam == "ref":
                caso["spearman_banda"] = spearman(bons, g, "banda")
                caso["spearman_campo"] = spearman(bons, g, "campo")
        det_pasta = pasta + "_det" if (fam == "ref" or fam in DET_PROPRIO) else os.path.join(
            raiz, f"{prefixo}_fase1", BASE + "_det")
        caso["det"] = deterministico(det_pasta)
        if caso["det"] is None:
            avisos.append(f"{cid}: sem deterministico em {os.path.relpath(det_pasta, raiz)}")
        casos.append(caso)
    return {"raiz": os.path.abspath(raiz), "prefixo": prefixo, "casos": casos, "versoes": versoes, "avisos": avisos}


def _limpa(x):
    if isinstance(x, float):
        return None if not math.isfinite(x) else x
    if isinstance(x, dict):
        return {k: _limpa(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_limpa(v) for v in x]
    return x


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    est = calcula(sys.argv[1], sys.argv[2])
    with open(sys.argv[3], "w") as f:
        json.dump(_limpa(est), f)
    for c in est["casos"]:
        e = c.get("estat")
        if e:
            print("%-14s S=%5d P_f=%6.2f%% CoV(P_f)=%5.1f%% mu=%.3f sigma=%.3f CoV=%.1f%% det=%s" % (
                c["id"], e["S"], 100 * e["P_f"], 100 * e["CoV_Pf"], e["media"], e["desvio"], 100 * e["CoV"],
                "%.4f" % c["det"]["Gamma"] if c.get("det") else "-"))
    for a in est["avisos"]:
        print("aviso:", a)
    print("versoes:", est["versoes"])
