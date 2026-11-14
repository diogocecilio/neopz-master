#!/usr/bin/env python3
"""dados_chuva.py: le os resultados do Monte Carlo com chuva (SlopeRainMonteCarlo, montecarlo.cpp: talude da Fig. 2 de
Rahardjo et al. 2007, H = 10 m a 45 graus, sob chuva, com campos log-normais de c', phi' e k_s como em Vargas Ceron et
al. 2025) e calcula, por caso, as grandezas das Secoes 6.1 a 6.5 do Ceron do mesmo jeito que o dados_sat.py do
rebaixamento (SlopeRandomField/comparacao_ceron), no instante das familias (t = 24 h, o fim da chuva, se nada for dito):

- S (realizacoes usadas), n_f = #(FS < 1), P_f = n_f / S e CoV(P_f) = sqrt((1 - P_f) / (S P_f)) (Eqs. 67 a 69 do
  Ceron), o intervalo de 95% de Wilson de P_f e o P_f do ajuste log-normal;
- o ajuste log-normal (maxima verossimilhanca): media, mediana, desvio e CoV, que sao os mu, sigma e CoV das tabelas
  do Ceron; as estatisticas amostrais (e os quantis de 5% e 95%) ao lado;
- o qui-quadrado dos ajustes normal e log-normal (k = ceil(2 S^0,4) classes equiprovaveis, gl = k - 3);
- a convergencia de P_f com o numero de simulacoes, na ordem das realizacoes;
- o Spearman parcial de 2a ordem entre FS e (c', phi', k) (Eqs. 70 e 71) com o IC de 95% (Fisher), com duas definicoes
  das variaveis de cada realizacao: (banda_c, banda_phi_deg, massa_k_mult) e as medias do campo (c_med, phi_med,
  k_mult_med);
- o FS deterministico (a linha real = -1, rodada com det=1) e a banda dele (mc_real-1/mecanismo.csv, se gravada); nas
  familias em que o deterministico nao muda com o valor (covk, covc, covphi, s), o da referencia quando o caso nao tem
  o seu.

E os da chuva:
- as mesmas estatisticas em cada instante de cada caso (na referencia, a curva da duracao: o FS em t e o do fim de uma
  chuva de t horas com a mesma intensidade);
- os pares FS(0) e FS(t) de cada realizacao (casos com os dois instantes);
- o mapa da frequencia da banda do mecanismo (mc_mapa_*.csv, contagens somadas nas partes) em cada instante;
- a categoria do nivel d'agua de cada realizacao (Tabela 3 do Rahardjo: I sobe ate a crista, II ate o meio da face,
  III ate o pe, IV sem subida ou subida pequena), pelas cotas nivel_x30, nivel_x35 e nivel_x40 contra as cotas da
  superficie (30, 25 e 20 m), e a subida de cada vertical desde o nivel inicial (a linha deterministica em t = 0).

Estrutura lida: <raiz>/<caso>/, um caso por pasta, com o nome do caso: ref, covk=<v>, covc=<v>, covphi=<v>, s=<v>,
alfa=<v> ou I=<v> (I em mm/h). Em cada pasta: as partes <prefixo>_<a>_<b>.csv (juntadas; uma linha repetida de
(real, t_h) entra uma vez so), <prefixo>_parametros.txt, <prefixo>_mapa_<a>_<b>.csv e <prefixo>_real<n>/. O valor da
referencia em cada familia vem do <prefixo>_parametros.txt da pasta ref. Uma realizacao entra num instante se tem
situacao = ok e FS_min > 0 finito; as outras sao contadas e listadas. Os dados so sao lidos.

Uso: python3 -I dados_chuva.py <raiz> <saida.json> [t=24] [prefixo=mc]"""

import csv
import json
import math
import os
import re
import sys

import numpy as np
from scipy import stats

FAMILIAS = ("covk", "covc", "covphi", "s", "alfa", "I")
# nome da familia na pasta (minusculas) -> familia
NOME_PASTA = {"covk": "covk", "covc": "covc", "covphi": "covphi", "s": "s", "alfa": "alfa", "i": "I"}
# valores de referencia, se a pasta ref nao tiver o mc_parametros.txt
VALOR_REF_PADRAO = {"covk": 0.6, "covc": 0.3, "covphi": 0.1, "s": 1.0, "alfa": 1.0, "I": 33.0}
# a chave do parametros.txt de cada familia (s = LX / 20)
CHAVE_PARAMETRO = {"covk": "COV_K", "covc": "COV_C", "covphi": "COV_PHI", "s": "LX", "alfa": "ALFA", "I": "I"}
# a coluna do csv de cada familia, para conferir o nome da pasta
COLUNA_FAMILIA = {"covk": "cov_k", "covc": "cov_c", "covphi": "cov_phi", "s": "Lx", "alfa": "alfa", "I": "I"}
# familias cujo deterministico muda com o valor (as outras usam o da referencia se o caso nao tem o seu)
DET_PROPRIO = ("alfa", "I")
# as variaveis da correlacao de Spearman: (nome, coluna da banda/massa, coluna da media do campo)
VARS_SPEARMAN = (("c", "banda_c", "c_med"), ("phi", "banda_phi_deg", "phi_med"), ("k", "massa_k_mult", "k_mult_med"))
N_MIN_SPEARMAN = 10
# nivel d'agua: as verticais do montecarlo.cpp (crista, meio da face e pe) e a cota da superficie em cada uma
VERTICAIS = (("x30", "nivel_x30", 30.0), ("x35", "nivel_x35", 25.0), ("x40", "nivel_x40", 20.0))
TOL_NIVEL = 0.15  # NivelAgua sobe de 0,1 em 0,1 m: o nivel na superficie aparece ate 0,1 m abaixo dela
CATEGORIAS = ("I", "II", "III", "IV")
# Rahardjo et al. (2007), solo f50,-5 a 45 graus: F_min(I) da Fig. 7b (Eq. 6, I em mm/h) e F_ini da Tabela 7
RAHARDJO = {"F_ini": 2.118, "sigmoide": {"a": 0.82, "b": 1.83, "c": 10.7, "d": -13.13}}
BANDA_LIMIAR = 0.1  # banda: sqrt(J2(eps_p)) >= 0,1 do maior, como no montecarlo.cpp
TALUDE = [[0, 0], [70, 0], [70, 20], [40, 20], [30, 30], [0, 30]]  # Fig. 2 do Rahardjo (montecarlo.cpp, kTalude)


def rahardjo_fmin(I):
    """F_min(I) da sigmoide da Fig. 7b do Rahardjo (f50,-5, 45 graus): a + b / (1 + exp(-(I - c) / d))"""
    s = RAHARDJO["sigmoide"]
    return s["a"] + s["b"] / (1. + math.exp(-(I - s["c"]) / s["d"]))


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
        with open(arq) as f:
            for linha in f:
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


def wilson(nf, S, z=1.96):
    """intervalo de 95% de Wilson de uma proporcao (nao some quando n_f = 0)"""
    if S <= 0:
        return (math.nan, math.nan)
    p = nf / S
    den = 1. + z * z / S
    c = (p + z * z / (2. * S)) / den
    h = z * math.sqrt(p * (1. - p) / S + z * z / (4. * S * S)) / den
    return (max(0., c - h), min(1., c + h))


def qui2(g, dist, k):
    """p-valor do qui-quadrado com k classes equiprovaveis da distribuicao ajustada (gl = k - 3)"""
    bordas = dist.ppf(np.linspace(0., 1., k + 1))
    obs = np.histogram(g, bins=np.concatenate([[-np.inf], bordas[1:-1], [np.inf]]))[0]
    esp = len(g) / k
    x2 = float(((obs - esp) ** 2 / esp).sum())
    gl = k - 3
    return {"x2": x2, "gl": gl, "p": float(stats.chi2.sf(x2, gl)) if gl > 0 else math.nan}


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
    pf_ln = 0.5 * math.erfc(ml / (sl * math.sqrt(2.))) if sl > 0 else float(ml < 0.)
    e = {
        "S": S, "n_f": nf, "P_f": pf, "CoV_Pf": covpf, "IC95_Pf": wilson(nf, S), "P_f_lognormal": pf_ln,
        "mu_ln": ml, "sigma_ln": sl, "media": media, "mediana": math.exp(ml), "desvio": media * cov, "CoV": cov,
        "media_amostral": float(g.mean()), "desvio_amostral": float(g.std(ddof=1)) if S > 1 else math.nan,
        "CoV_amostral": float(g.std(ddof=1) / g.mean()) if S > 1 else math.nan, "minimo": float(g.min()),
        "maximo": float(g.max()), "q05": float(np.quantile(g, 0.05)), "q95": float(np.quantile(g, 0.95)),
        "qui2_classes": None, "qui2_normal": None, "qui2_lognormal": None,
    }
    k = int(math.ceil(2. * S ** 0.4))
    if S >= 5 and sl > 0 and g.std() > 0 and k > 3:
        e["qui2_classes"] = k
        e["qui2_normal"] = qui2(g, stats.norm(g.mean(), g.std(ddof=0)), k)
        e["qui2_lognormal"] = qui2(g, stats.lognorm(s=sl, scale=math.exp(ml)), k)
    return e


def spearman(linhas, g, qual):
    """Spearman parcial de 2a ordem entre FS e cada variavel, controlando as outras duas"""
    col = 1 if qual == "banda" else 2
    X = [np.asarray(g, float)] + [np.array([_f(r[v[col]]) for r in linhas]) for v in VARS_SPEARMAN]
    D = np.column_stack(X)
    ok = np.all(np.isfinite(D), axis=1)
    D = D[ok]
    n = len(D)
    if n < N_MIN_SPEARMAN or np.any(np.std(D, axis=0) == 0):
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


def _num(s):
    """'covc=0.5' -> 0.5 (None se nao for numero)"""
    try:
        return float(s)
    except ValueError:
        return None


def caso_da_pasta(nome):
    """(familia, valor) do nome da pasta, ou None"""
    if nome == "ref":
        return ("ref", None)
    m = re.fullmatch(r"([A-Za-z]+)=([-+0-9.eE]+)", nome)
    if not m or m.group(1).lower() not in NOME_PASTA or _num(m.group(2)) is None:
        return None
    return (NOME_PASTA[m.group(1).lower()], _num(m.group(2)))


def mesmo_t(a, b):
    return abs(a - b) < 1e-6


def ok(r):
    return r["situacao"] == "ok" and math.isfinite(r["fs"]) and r["fs"] > 0


# ---- leitura de uma pasta
def le_partes(pasta, prefixo, avisos, cid):
    """junta as partes <prefixo>_<a>_<b>.csv; devolve (linhas, nomes das partes)"""
    re_parte = re.compile(r"^%s_(-?\d+)_(-?\d+)\.csv$" % re.escape(prefixo))
    partes = sorted((f for f in os.listdir(pasta) if re_parte.match(f)),
                    key=lambda f: tuple(int(x) for x in re_parte.match(f).groups()))
    linhas, vistas, repetidas = [], {}, 0
    for p in partes:
        for r in le_csv(os.path.join(pasta, p)):
            real, t = _f(r.get("real")), _f(r.get("t_h"))
            if not (math.isfinite(real) and math.isfinite(t)):
                continue
            r["real"], r["t"], r["fs"], r["parte"] = int(real), t, _f(r.get("FS_min")), p
            chave = (r["real"], round(t, 6))
            if chave in vistas:
                repetidas += r["real"] >= 0  # o deterministico em cada parte (det=1 em todas) e normal
                if not (vistas[chave]["fs"] == r["fs"] or (math.isnan(vistas[chave]["fs"]) and math.isnan(r["fs"]))):
                    avisos.append("%s: real %d, t = %g h repetida com FS diferente (%s: %g; %s: %g); fica a primeira" % (
                        cid, r["real"], t, vistas[chave]["parte"], vistas[chave]["fs"], p, r["fs"]))
                continue
            vistas[chave] = r
            linhas.append(r)
    if repetidas:
        avisos.append("%s: %d linhas de realizacoes repetidas entre as partes (partes sobrepostas?), contadas uma vez"
                      % (cid, repetidas))
    return linhas, partes


def le_mapa(pasta, prefixo):
    """soma as contagens das partes do mapa: {t: {(x, y): [n, n_falha]}} e n_real somado (das partes com linhas)"""
    re_mapa = re.compile(r"^%s_mapa_(-?\d+)_(-?\d+)\.csv$" % re.escape(prefixo))
    soma, n_real, arqs = {}, 0, []
    for f in sorted(os.listdir(pasta)):
        if not re_mapa.match(f):
            continue
        arqs.append(f)
        linhas = le_csv(os.path.join(pasta, f))
        if linhas:
            n_real += int(_f(linhas[0]["n_real"]))
        for r in linhas:
            t = round(_f(r["t_h"]), 6)
            c = soma.setdefault(t, {}).setdefault((_f(r["x"]), _f(r["y"])), [0, 0])
            c[0] += int(_f(r["n"]))
            c[1] += int(_f(r["n_falha"]))
    return soma, n_real, arqs


def le_realizacao_figura(pasta):
    """mc_real<n>/: a banda de cada instante (mecanismo.csv) e o p_w nos vertices (vertices, triangles, snaps_pw)"""
    out = {}
    mec = os.path.join(pasta, "mecanismo.csv")
    if os.path.isfile(mec):
        por_t = {}
        for m in le_csv(mec):
            por_t.setdefault(round(_f(m["t_h"]), 6), []).append((_f(m["x"]), _f(m["y"]), _f(m["epsp_raizJ2"])))
        bandas = []
        for t in sorted(por_t):
            pts = por_t[t]
            emax = max(p[2] for p in pts)
            bandas.append({"t": t, "pontos": [[round(x, 3), round(y, 3)] for x, y, e in pts if e >= BANDA_LIMIAR * emax]})
        out["bandas"] = bandas
    arqs = [os.path.join(pasta, n) for n in ("vertices.csv", "triangles.csv", "snaps_pw.csv")]
    if all(os.path.isfile(a) for a in arqs):
        v = le_csv(arqs[0])
        tri = [[int(r["n0"]), int(r["n1"]), int(r["n2"])] for r in le_csv(arqs[1])]
        snaps = []
        with open(arqs[2]) as f:
            for linha in f:
                vals = [_f(x) for x in linha.strip().split(",") if x != ""]
                if len(vals) == len(v) + 1:
                    snaps.append({"t": vals[0], "pw": [round(x, 3) for x in vals[1:]]})
        if snaps and tri:
            out["malha_pw"] = {"x": [_f(r["x"]) for r in v], "y": [_f(r["y"]) for r in v], "tri": tri,
                               "instantes": snaps}
    fs = os.path.join(pasta, "fs.csv")
    if os.path.isfile(fs):
        out["fs"] = [[_f(r["t_h"]), _f(r["FS"])] for r in le_csv(fs)]
    return out or None


# ---- nivel d'agua (Tabela 3 do Rahardjo)
def niveis(r):
    return {k: _f(r.get(col)) for k, col, _ in VERTICAIS}


def categoria(nv):
    """I: o nivel chega a superficie na crista (x = 30 m); II: no meio da face (x = 35 m); III: no pe (x = 40 m);
    IV: nao chega a superficie em nenhuma das tres (sem subida ou subida pequena)"""
    for (k, _, cota), cat in zip(VERTICAIS, CATEGORIAS):
        if math.isfinite(nv[k]) and nv[k] >= cota - TOL_NIVEL:
            return cat
    if all(math.isfinite(nv[k]) for k, _, _ in VERTICAIS):
        return "IV"
    return None


def resumo_niveis(bons, iniciais, det):
    """categorias, subidas e P_f por categoria das realizacoes bons num instante"""
    cont = {c: 0 for c in CATEGORIAS}
    nf = {c: 0 for c in CATEGORIAS}
    soma_fs = {c: 0. for c in CATEGORIAS}
    por_real, subidas = [], {k: [] for k, _, _ in VERTICAIS}
    for r in bons:
        nv = niveis(r)
        cat = categoria(nv)
        por_real.append([r["real"], cat, nv["x30"], nv["x35"], nv["x40"], r["fs"]])
        if cat:
            cont[cat] += 1
            nf[cat] += int(r["fs"] < 1.)
            soma_fs[cat] += r["fs"]
        if iniciais:
            for k, _, _ in VERTICAIS:
                if math.isfinite(nv[k]) and math.isfinite(iniciais[k]):
                    subidas[k].append(nv[k] - iniciais[k])
    S = sum(cont.values())
    out = {"contagem": cont, "fracao": {c: (cont[c] / S if S else math.nan) for c in CATEGORIAS},
           "n_f": nf, "fs_medio": {c: (soma_fs[c] / cont[c] if cont[c] else math.nan) for c in CATEGORIAS},
           "P_f_condicional": {c: (nf[c] / cont[c] if cont[c] else math.nan) for c in CATEGORIAS},
           "nivel_medio": {k: _media_finita([p[i] for p in por_real]) for i, (k, _, _) in enumerate(VERTICAIS, start=2)},
           "subida_media": {k: _media_finita(v) for k, v in subidas.items()},
           "por_real": por_real}
    if det:
        out["det_niveis"] = niveis(det)
        out["det_categoria"] = categoria(out["det_niveis"])
    return out


# ---- um caso
def resumo_instante(linhas, t, iniciais, det):
    """as estatisticas de um instante (realizacoes aleatorias com situacao ok)"""
    bons = sorted((r for r in linhas if mesmo_t(r["t"], t) and r["real"] >= 0 and ok(r)), key=lambda r: r["real"])
    g = np.array([r["fs"] for r in bons])
    out = {"t": t, "S": len(bons), "estat": estatisticas(g) if len(g) else None,
           "det_FS": det["fs"] if det else None, "reais": [r["real"] for r in bons], "fs": g.tolist(),
           "spearman_banda": spearman(bons, g, "banda") if len(g) else None,
           "spearman_campo": spearman(bons, g, "campo") if len(g) else None,
           "nivel": resumo_niveis(bons, iniciais, det)}
    for col in ("banda_prof", "banda_area_m2", "banda_xc", "banda_yc", "banda_pw", "nos_face", "eq_fator", "s_fator"):
        out[col + "_medio"] = _media_finita([_f(r.get(col)) for r in bons])
    return out, bons


def checa_parametros(caso, linhas, ref_param, avisos):
    """os parametros das linhas (unicos?) e o valor do nome da pasta contra a coluna da familia"""
    for col in ("solo", "I", "chuva", "hw", "metodo", "cov_c", "cov_phi", "cov_k", "Lx", "Ly", "alfa"):
        vals = sorted({r.get(col) for r in linhas})
        caso[col] = vals[0] if len(vals) == 1 else vals
        if len(vals) > 1:
            avisos.append("%s: a coluna %s tem valores diferentes entre as linhas: %s" % (caso["id"], col, vals))
    fam, v = caso["familia"], caso["valor"]
    if fam in COLUNA_FAMILIA and not isinstance(caso[COLUNA_FAMILIA[fam]], list):
        lido = _f(caso[COLUNA_FAMILIA[fam]])
        if fam == "s":
            lido /= 20.
        if math.isfinite(lido) and abs(lido - v) > 1e-6 * max(1., abs(v)):
            avisos.append("%s: o nome da pasta diz %s = %g, mas o csv tem %g" % (caso["id"], fam, v, lido))
    # os outros parametros do caso devem ser os da referencia
    if ref_param:
        chaves = {"covk": "COV_K", "covc": "COV_C", "covphi": "COV_PHI", "s": "LX", "alfa": "ALFA", "I": "I"}
        p = caso["parametros"]
        dif = []
        for f2, k in chaves.items():
            if f2 == fam or k not in p or k not in ref_param:
                continue
            if abs(_f(p[k]) - _f(ref_param[k])) > 1e-9 * max(1., abs(_f(ref_param[k]))):
                dif.append("%s = %s (ref. %s)" % (k, p[k], ref_param[k]))
        for k in ("SOLO", "HW", "CHUVA", "METODO", "MALHA_TRANSIENTE", "MALHA_FATOR", "NREF"):
            if k in p and k in ref_param and p[k] != ref_param[k]:
                dif.append("%s = %s (ref. %s)" % (k, p[k], ref_param[k]))
        if dif:
            avisos.append("%s: parametros diferentes da referencia alem da familia: %s" % (caso["id"], "; ".join(dif)))


def calcula_caso(raiz, nome, fam, valor, prefixo, t_fam, ref_info, avisos):
    pasta = os.path.join(raiz, nome)
    cid = nome
    linhas, partes = le_partes(pasta, prefixo, avisos, cid)
    if not linhas:
        avisos.append("%s: sem linhas nas partes %s_<a>_<b>.csv" % (cid, prefixo))
        return None
    param = le_parametros(os.path.join(pasta, prefixo + "_parametros.txt"))
    caso = {"id": cid, "familia": fam, "valor": valor, "pasta": os.path.relpath(pasta, raiz), "parametros": param,
            "partes": partes, "linhas": len(linhas)}
    checa_parametros(caso, linhas, ref_info.get("parametros") if fam != "ref" else None, avisos)
    instantes = sorted({round(r["t"], 6) for r in linhas})
    caso["instantes_h"] = instantes
    # instantes adaptativos em todas as realizacoes (nas das figuras todos sao adaptativos)
    caso["adaptativos_h"] = [t for t in instantes if all(r.get("adaptativo") == "1" for r in linhas
                                                         if r["real"] >= 0 and mesmo_t(r["t"], t) and ok(r))]
    reais = sorted({r["real"] for r in linhas if r["real"] >= 0})
    caso["realizacoes"] = len(reais)
    caso["excluidas"] = [{"real": r["real"], "t_h": r["t"], "situacao": r["situacao"], "FS": r["fs"]}
                         for r in linhas if r["real"] >= 0 and not ok(r)]
    caso["KL_variancia_media"] = _media_finita([_f(r.get("KL_variancia")) for r in linhas])
    # custo medio por realizacao: transiente (acumulado ate o ultimo instante) e fator (soma dos instantes)
    sT, sF = {}, {}
    for r in linhas:
        if r["real"] >= 0:
            sT[r["real"]] = max(sT.get(r["real"], 0.), _f(r.get("s_transiente")) if math.isfinite(
                _f(r.get("s_transiente"))) else 0.)
            sF[r["real"]] = sF.get(r["real"], 0.) + (_f(r.get("s_fator")) if math.isfinite(_f(r.get("s_fator"))) else 0.)
    caso["custo_medio_s"] = {"transiente": _media_finita(list(sT.values())), "fator": _media_finita(list(sF.values()))}

    # o instante das familias
    t = next((x for x in instantes if mesmo_t(x, t_fam)), None)
    if t is None:
        t = instantes[-1]
        avisos.append("%s: sem o instante %g h; as estatisticas do caso sao as de %g h" % (cid, t_fam, t))
    caso["t"] = t

    # deterministico: o do caso (real = -1) ou, nas familias em que ele nao muda com o valor, o da referencia
    det_linhas = {round(r["t"], 6): r for r in linhas if r["real"] == -1 and ok(r)}
    fig_det = le_realizacao_figura(os.path.join(pasta, prefixo + "_real-1"))
    det_origem = "proprio" if det_linhas else None
    caso["_det_linhas"] = det_linhas
    if not det_linhas and fam not in DET_PROPRIO and ref_info.get("det_linhas"):
        # o da referencia (as mesmas medias); a banda e o p_w ficam so na referencia
        det_linhas, fig_det, det_origem = ref_info["det_linhas"], None, "ref"
    if not det_linhas:
        avisos.append("%s: sem deterministico (rodar com det=1)" % cid)

    # niveis iniciais: a linha deterministica em t = 0; senao uma realizacao em t = 0; senao os da referencia
    iniciais, origem = None, None
    if 0.0 in det_linhas:
        iniciais, origem = niveis(det_linhas[0.0]), "deterministico em t = 0"
    else:
        r0 = next((r for r in linhas if mesmo_t(r["t"], 0.) and ok(r)), None)
        if r0:
            iniciais, origem = niveis(r0), "realizacao %d em t = 0" % r0["real"]
        elif ref_info.get("iniciais") and ref_info.get("hw") == param.get("HW", ref_info.get("hw")):
            iniciais, origem = ref_info["iniciais"], "referencia (mesmo nivel inicial)"
    caso["niveis_iniciais"] = iniciais
    caso["niveis_iniciais_origem"] = origem

    caso["instantes"] = []
    bons_t = None
    for ti in instantes:
        d = det_linhas.get(ti)
        res, bons = resumo_instante(linhas, ti, iniciais, d)
        caso["instantes"].append(res)
        if mesmo_t(ti, t):
            bons_t, res_t = bons, res
    caso["estat"] = res_t["estat"]
    caso["fs"] = res_t["fs"]
    caso["reais"] = res_t["reais"]
    g = np.asarray(res_t["fs"])
    caso["conv_pf"] = (np.cumsum(g < 1.) / np.arange(1, len(g) + 1)).tolist() if len(g) else []
    caso["spearman_banda"] = res_t["spearman_banda"]
    caso["spearman_campo"] = res_t["spearman_campo"]
    caso["nivel"] = res_t["nivel"]
    caso["S_excluidas_t"] = sum(1 for e in caso["excluidas"] if mesmo_t(e["t_h"], t))

    # o deterministico no instante das familias e nos outros
    if det_linhas:
        dt = det_linhas.get(t)
        det = {"origem": det_origem, "FS": dt["fs"] if dt else None,
               "FS_t": [[ti, det_linhas[ti]["fs"]] for ti in sorted(det_linhas)]}
        if dt:
            for c in ("banda_area_m2", "banda_prof", "banda_xmin", "banda_xmax", "banda_xc", "banda_yc", "banda_ymin",
                      "banda_pw", "massa_area_m2", "massa_xc", "massa_yc", "nos_face", "c_med", "phi_med", "k_mult_med"):
                det[c] = _f(dt.get(c))
            det["niveis"] = niveis(dt)
            det["categoria"] = categoria(det["niveis"])
        det["niveis_t"] = [[ti, niveis(det_linhas[ti]), categoria(niveis(det_linhas[ti]))] for ti in sorted(det_linhas)]
        if fig_det:
            det.update(fig_det)
        caso["det"] = det
    else:
        caso["det"] = None

    # pares FS(0) e FS(t) de cada realizacao
    if 0.0 in instantes and not mesmo_t(t, 0.):
        fs0 = {r["real"]: r["fs"] for r in linhas if mesmo_t(r["t"], 0.) and r["real"] >= 0 and ok(r)}
        cats = {p[0]: p[1] for p in res_t["nivel"]["por_real"]}
        pares = [(r["real"], fs0[r["real"]], r["fs"], cats.get(r["real"])) for r in bons_t if r["real"] in fs0]
        if pares:
            a0, a1 = np.array([p[1] for p in pares]), np.array([p[2] for p in pares])
            caso["pares"] = {"t0": 0.0, "t1": t, "reais": [p[0] for p in pares], "fs0": a0.tolist(), "fs1": a1.tolist(),
                             "categoria_t1": [p[3] for p in pares],
                             "pearson": float(np.corrcoef(a0, a1)[0, 1]) if len(pares) > 2 and a0.std() > 0 and
                             a1.std() > 0 else math.nan,
                             "spearman": float(stats.spearmanr(a0, a1)[0]) if len(pares) > 2 and a0.std() > 0 and
                             a1.std() > 0 else math.nan,
                             "razao_media": float(np.mean(a1 / a0)), "razao_desvio": float(np.std(a1 / a0, ddof=1))
                             if len(pares) > 1 else math.nan,
                             "det": [det_linhas[0.0]["fs"], det_linhas[t]["fs"]] if (0.0 in det_linhas and
                                                                                   t in det_linhas) else None}

    # mapa da banda: contagens somadas nas partes, por instante; S e n_f (os divisores das fracoes) so das realizacoes
    # na faixa <a>..<b> de algum <prefixo>_mapa_<a>_<b>.csv (o montecarlo.cpp grava o mapa no fim: um processo
    # interrompido deixa linhas sem mapa)
    soma, n_real, arqs = le_mapa(pasta, prefixo)
    re_mapa = re.compile(r"^%s_mapa_(-?\d+)_(-?\d+)\.csv$" % re.escape(prefixo))
    faixas = [tuple(int(x) for x in re_mapa.match(f).groups()) for f in arqs]
    reais_mapa = {x for x in reais if any(a <= x <= b for a, b in faixas)}
    sem_mapa = [x for x in reais if x not in reais_mapa]
    if arqs and sem_mapa:
        avisos.append("%s: %d realizacoes sem arquivo do mapa (%d a %d), fora das fracoes do mapa" % (
            cid, len(sem_mapa), sem_mapa[0], sem_mapa[-1]))
    mapa = {"arquivos": arqs, "n_real_arquivos": n_real, "instantes": []}
    for ti in instantes:
        bons_i = [r for r in linhas if mesmo_t(r["t"], ti) and r["real"] in reais_mapa and ok(r)]
        cel = soma.get(round(ti, 6), {})
        mapa["instantes"].append({"t": ti, "S": len(bons_i), "n_f": sum(1 for r in bons_i if r["fs"] < 1.),
                                  "celulas": [[x, y, n, nf] for (x, y), (n, nf) in sorted(cel.items())]})
    if arqs and n_real and n_real != len(reais):
        avisos.append("%s: o mapa soma %d realizacoes e as partes tem %d" % (cid, n_real, len(reais)))
    caso["mapa"] = mapa
    return caso


def calcula(raiz, t_fam=24.0, prefixo="mc"):
    avisos = []
    pastas = sorted(d for d in os.listdir(raiz) if os.path.isdir(os.path.join(raiz, d)))
    achados = []
    for d in pastas:
        fv = caso_da_pasta(d)
        if fv is None:
            avisos.append("%s: pasta ignorada (o nome nao e ref nem <familia>=<valor>)" % d)
            continue
        achados.append((d,) + fv)
    # a referencia primeiro: o valor de referencia de cada familia, o deterministico e os niveis iniciais dela
    achados.sort(key=lambda a: (a[1] != "ref", FAMILIAS.index(a[1]) if a[1] in FAMILIAS else -1, a[2] or 0.))
    ref_info, casos = {}, []
    valor_ref = dict(VALOR_REF_PADRAO)
    for nome, fam, valor in achados:
        caso = calcula_caso(raiz, nome, fam, valor, prefixo, t_fam, ref_info, avisos)
        if caso is None:
            continue
        if fam == "ref":
            p = caso["parametros"]
            for f2, k in CHAVE_PARAMETRO.items():
                if k in p and math.isfinite(_f(p[k])):
                    valor_ref[f2] = _f(p[k]) / (20. if f2 == "s" else 1.)
            ref_info = {"parametros": p, "hw": p.get("HW"), "iniciais": caso.get("niveis_iniciais"),
                        "det_linhas": caso["_det_linhas"]}
        caso.pop("_det_linhas", None)
        casos.append(caso)
    if not any(c["familia"] == "ref" for c in casos):
        avisos.append("sem a pasta ref: as familias ficam sem o ponto de referencia")
    # um caso com o valor de referencia repete a referencia
    for c in casos:
        if c["familia"] in FAMILIAS and abs(c["valor"] - valor_ref[c["familia"]]) < 1e-9:
            avisos.append("%s: tem o valor de referencia da familia; a figura usa este caso no lugar de ref" % c["id"])
    I_ref = valor_ref["I"]
    return {"raiz": os.path.abspath(raiz), "prefixo": prefixo, "t_familias": t_fam, "valor_ref": valor_ref,
            "rahardjo": dict(RAHARDJO, F_min_ref=rahardjo_fmin(I_ref)), "talude": TALUDE,
            "verticais": [[k, col, cota] for k, col, cota in VERTICAIS], "tol_nivel": TOL_NIVEL,
            "banda_limiar": BANDA_LIMIAR, "casos": casos, "avisos": avisos}


def _limpa(x):
    if isinstance(x, float):
        return None if not math.isfinite(x) else x
    if isinstance(x, (np.floating,)):
        return _limpa(float(x))
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, dict):
        return {k: _limpa(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_limpa(v) for v in x]
    return x


def _pct(v):
    return "   -  " if v is None or not math.isfinite(v) else "%6.2f" % (100. * v)


def imprime(est):
    print("instante das familias: t = %g h; valores de referencia: %s" % (
        est["t_familias"], ", ".join("%s = %g" % kv for kv in est["valor_ref"].items())))
    print("%-12s %5s %5s %8s %8s %8s %7s %7s %6s  %-4s %8s" % ("caso", "t(h)", "S", "P_f(%)", "P_fLN(%)", "mu", "sigma",
                                                             "CoV(%)", "FSdet", "cat", "I/II/III/IV"))
    for c in est["casos"]:
        e = c.get("estat")
        det = c.get("det") or {}
        nv = c.get("nivel") or {}
        cont = nv.get("contagem", {})
        if e:
            print("%-12s %5g %5d %8s %8s %8.4f %7.4f %7.2f %6s  %-4s %8s" % (
                c["id"], c["t"], e["S"], _pct(e["P_f"]), _pct(e["P_f_lognormal"]), e["media"], e["desvio"],
                100 * e["CoV"], "%.4f" % det["FS"] if det.get("FS") else "  -  ", det.get("categoria") or "-",
                "/".join(str(cont.get(k, 0)) for k in CATEGORIAS)))
        else:
            print("%-12s sem realizacoes validas no instante %g h" % (c["id"], c["t"]))
    ref = next((c for c in est["casos"] if c["familia"] == "ref"), None)
    if ref and len(ref["instantes"]) > 1:
        print("\nreferencia por instante (duracao da chuva):")
        print("  t(h)     S   FSdet      mu   sigma   P_f(%)  P_fLN(%)  I/II/III/IV")
        for i in ref["instantes"]:
            e = i["estat"]
            if not e:
                continue
            print("  %4g %5d  %6s  %6.4f  %6.4f  %7s  %8s  %s" % (
                i["t"], e["S"], "%.4f" % i["det_FS"] if i["det_FS"] else "  -  ", e["media"], e["desvio"],
                _pct(e["P_f"]), _pct(e["P_f_lognormal"]),
                "/".join(str(i["nivel"]["contagem"][k]) for k in CATEGORIAS)))
        if ref.get("pares"):
            p = ref["pares"]
            print("  FS(0) x FS(%g): %d pares, Pearson %.3f, Spearman %.3f, FS(%g)/FS(0) medio %.3f" % (
                p["t1"], len(p["fs0"]), p["pearson"], p["spearman"], p["t1"], p["razao_media"]))
        sp = ref.get("spearman_banda")
        if sp:
            print("  Spearman parcial (banda e massa) em %g h: c %.3f, phi %.3f, k %.3f" % (
                ref["t"], sp["c"]["rho_parcial"], sp["phi"]["rho_parcial"], sp["k"]["rho_parcial"]))
    print("Rahardjo (f50, 45 graus): F_ini = %.3f; F_min(I = %g mm/h) = %.3f" % (
        RAHARDJO["F_ini"], est["valor_ref"]["I"], est["rahardjo"]["F_min_ref"]))
    for a in est["avisos"]:
        print("aviso:", a)


def main(argv):
    if len(argv) < 3:
        sys.exit(__doc__)
    raiz, saida = argv[1], argv[2]
    op = {"t": "24", "prefixo": "mc"}
    for a in argv[3:]:
        k, _, v = a.partition("=")
        if k not in op or not v:
            sys.exit("opcao desconhecida: %s\n\n%s" % (a, __doc__))
        op[k] = v
    est = calcula(raiz, float(op["t"]), op["prefixo"])
    d = os.path.dirname(os.path.abspath(saida))
    os.makedirs(d, exist_ok=True)
    with open(saida, "w") as f:
        json.dump(_limpa(est), f)
    imprime(est)
    print("gravado:", saida)


if __name__ == "__main__":
    main(sys.argv)
