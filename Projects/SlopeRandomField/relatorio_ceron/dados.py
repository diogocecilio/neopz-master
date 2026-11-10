#!/usr/bin/env python3
"""dados.py: le as campanhas de Monte Carlo do SlopeRandomField e calcula, por caso, as grandezas da secao 6 de Vargas
Ceron et al. (2025, IJNAMG 49:2459, doi 10.1002/nag.3993):

- S (N realizacoes usadas), n_f = #(FS < 1), P_f = n_f/S (eqs. 67 e 69) e CoV(P_f) = sqrt((1 - P_f)/(S P_f)) (eq. 68);
- o ajuste log-normal de 2 parametros por maxima verossimilhanca (media, mediana, desvio e CoV ajustados, como os mu,
  sigma e CoV das Tabelas 3 a 6 do Ceron) e as estatisticas amostrais;
- o teste qui-quadrado dos ajustes normal e log-normal (k = ceil(2 S^0,4) classes equiprovaveis, gl = k - 3);
- a correlacao de postos de Spearman parcial de 2a ordem entre o FS e (c', phi', k_v) equivalentes (eqs. 70-71, Fig. 26);
- o FS deterministico (rodada com as medias, CAMPOS=0) de cada caso que o tem.

calcula(raiz_dados) -> dict serializavel em JSON (esquema no FILE CONTRACT do relatorio); grava(est, saida) grava
<saida>/estatisticas.json. Os dados sao so lidos.

Uma realizacao conta como concluida se e so se tem resultado.txt nao vazio. Cada realizacao e conferida: o FS de cada
instante em fs.csv e o de resumo.csv, FS_min e t_min_h de resumo.csv e resultado.txt, a banda e a massa do instante do
FS. Pastas fundidas num caso (f100 + f100_fase1, f100 + f100_fase3) precisam de parametros.txt identico, numeros de
realizacao disjuntos e os mesmos parametros em resumo.csv; senao, erro."""

import json
import math
import os
import sys

_AQUI = os.path.dirname(os.path.abspath(__file__))
_PAI = os.path.dirname(_AQUI)  # Projects/SlopeRandomField, onde esta sensibilidade_mc.py
for _p in (_PAI, _AQUI):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import numpy as np  # noqa: E402
from scipy import stats  # noqa: E402

from sensibilidade_mc import ic95, parcial, spearman  # noqa: E402

# ---- constantes do programa (main.cpp, valores padrao; parametros.txt tem "-" nessas chaves em todas as pastas usadas)
K_S = {"f10": 1.0e-4, "f50": 1.0e-5, "f100": 1.0e-6}  # p.ks, m/s
GAMMA_SOLO = 20.0  # p.gamma, kN/m3
GAMMA_W = 9.81  # kN/m3: k (m/s) / gamma_w = k na unidade do Ceron (m2 kPa^-1 s^-1 = m4 kN^-1 s^-1)
POROSIDADE, E_MPA, NU = 0.45, 20.0, 0.3  # p.n, p.E = 2e4 kPa, p.nu
C_MEDIO, PHI_MEDIO = 10.0, 26.0  # p.c (kPa), p.phi (graus)
TCHUVA_H = 24.0
L_REF = (20.0, 2.0)  # (Lx, Ly) de referencia: (Lx, Ly) = s (20 m, 2 m), eq. 72 do Ceron

# ---- os casos de cada solo
SUFIXOS = [  # (familia, sufixo da pasta, valor da variavel da familia)
    ("covk", "covk0", 0.0), ("covk", "covk0.75", 0.75), ("covk", "covk0.9", 0.9), ("covk", "covk1", 1.0),
    ("covc", "covc0.1", 0.1), ("covc", "covc0.5", 0.5), ("covc", "covc0.7", 0.7),
    ("covphi", "covphi0.05", 0.05), ("covphi", "covphi0.15", 0.15), ("covphi", "covphi0.2", 0.2),
    ("s", "lx40_ly4", 2.0), ("s", "lx100_ly10", 5.0), ("s", "lx8000_ly800", 400.0),
    ("alfa", "alfa2", 2.0), ("alfa", "alfa5", 5.0), ("alfa", "alfa10", 10.0),
]
VALOR_REF = {"covk": 0.6, "covc": 0.3, "covphi": 0.1, "s": 1.0, "alfa": 1.0}
# a variavel de cada familia, como campo do CASO
CAMPO_FAMILIA = {"covk": "cov_k", "covc": "cov_c", "covphi": "cov_phi", "s": "s", "alfa": "alfa", "chuva": "I_mm_h"}
FAMILIAS = ("ref", "covk", "covc", "covphi", "s", "alfa", "chuva")

# combinacoes que existem nos dados e nao entram no relatorio (sem analogo no Ceron)
MOTIVO_FTRANS1 = ("F_TRANS = 1 (resistência real também no transiente): sem análogo no Ceron")
MOTIVO_KRAMIN = ("KRA_MIN = 0,01 (piso da permeabilidade relativa ao ar 100 vezes maior que o padrão 0,0001): "
                 "sem análogo no Ceron")
EXCLUIDAS = {
    ("monte_carlo", "f50_I36_i"): MOTIVO_FTRANS1,
    ("monte_carlo", "f50_I36_i_det"): MOTIVO_FTRANS1,
    ("monte_carlo", "f50_I36_i_ftrans0.5_kramin0.01"): MOTIVO_KRAMIN,
    ("monte_carlo", "f50_I36_i_ftrans0.5_kramin0.01_det"): MOTIVO_KRAMIN,
    ("f100", "f100_I36_i_kramin0.01"): MOTIVO_KRAMIN,
    ("f100", "f100_I36_i_kramin0.01_det"): MOTIVO_KRAMIN,
    ("f100_fase3", "f100_I36_i_kramin0.01"): MOTIVO_KRAMIN,
}

# campos de resumo.csv que identificam o modelo e devem ser constantes num caso
CAMPOS_IDENT = ("solo", "I_mm_h", "variante", "f_trans", "kra_min", "ar_laterais", "campos", "semente", "KL_M",
                "KL_variancia", "cov_c", "cov_phi", "cov_k", "Lx", "Ly", "alfa", "area_solo_m2")
# o que vai para "amostras" (do instante do FS de cada realizacao)
CAMPOS_AMOSTRA = ("banda_ymin", "banda_xc", "banda_yc", "banda_c", "banda_phi_deg", "massa_k_mult")
# variaveis equivalentes da correlacao de Spearman (Hu et al.; Ceron, secao 6.1)
VARS_SPEARMAN = (("c", "banda_c"), ("phi", "banda_phi_deg"), ("k", "massa_k_mult"))


def _id_caso(solo, fam, v):
    return f"{solo}_{fam}{v:g}"


def _num_pt(v, casas=None):
    """numero com virgula decimal (texto simples)"""
    s = f"{v:g}" if casas is None else f"{v:.{casas}f}"
    return s.replace(".", ",")


def _lista_pt(itens):
    itens = list(itens)
    return itens[0] if len(itens) == 1 else ", ".join(itens[:-1]) + " e " + itens[-1]


def _rotulo(fam, v):
    if fam == "ref":
        return "Referência"
    if fam == "covk":
        return f"CoV(k_v) = {_num_pt(100 * v)}%"
    if fam == "covc":
        return f"CoV(c') = {_num_pt(100 * v)}%"
    if fam == "covphi":
        return f"CoV(φ') = {_num_pt(100 * v)}%"
    if fam == "s":
        return f"s = {_num_pt(v)}"
    if fam == "alfa":
        return f"α = {_num_pt(v)}"
    return f"I = {_num_pt(v)} mm/h"


def especificacao():
    """{solo: {id: {"pastas": [(campanha, combinacao)], "familia", "valor", "fs_def", "det": (campanha, comb) | None}}}"""
    esp = {"f50": {}, "f100": {}}
    # f50: referencia e casos de sensibilidade da campanha sensibilidade (FS so em 24 h); chuva da campanha monte_carlo
    esp["f50"]["f50_ref"] = {"pastas": [("sensibilidade", "f50_I36_i")], "familia": "ref", "valor": None,
                             "fs_def": "FS_24h", "det": ("monte_carlo", "f50_I36_i_ftrans0.5_det")}
    for fam, suf, v in SUFIXOS:
        esp["f50"][_id_caso("f50", fam, v)] = {"pastas": [("sensibilidade", f"f50_I36_i_{suf}")], "familia": fam,
                                               "valor": v, "fs_def": "FS_24h", "det": None}
    for txt in ("9", "36", "80"):
        esp["f50"][f"f50_I{txt}"] = {"pastas": [("monte_carlo", f"f50_I{txt}_i_ftrans0.5")], "familia": "chuva",
                                     "valor": float(txt), "fs_def": "FS_min",
                                     "det": ("monte_carlo", f"f50_I{txt}_i_ftrans0.5_det")}
    # f100: referencia = f100/f100_I36_i (0 a 299) + f100_fase1/f100_I36_i (300 a 999); FS minimo de 12 a 48 h
    pastas_ref = [("f100", "f100_I36_i"), ("f100_fase1", "f100_I36_i")]
    esp["f100"]["f100_ref"] = {"pastas": pastas_ref, "familia": "ref", "valor": None, "fs_def": "FS_min",
                               "det": ("f100", "f100_I36_i_det")}
    for fam, suf, v in SUFIXOS:
        esp["f100"][_id_caso("f100", fam, v)] = {
            "pastas": [("f100", f"f100_I36_i_{suf}")], "familia": fam, "valor": v, "fs_def": "FS_min",
            "det": ("f100", f"f100_I36_i_{suf}_det") if fam == "alfa" else None}
    for txt in ("3.6", "9", "36", "80"):
        if txt == "36":
            pastas, det = list(pastas_ref), ("f100", "f100_I36_i_det")
        else:
            pastas, det = [("f100", f"f100_I{txt}_i"), ("f100_fase3", f"f100_I{txt}_i")], ("f100", f"f100_I{txt}_i_det")
        esp["f100"][f"f100_I{txt}"] = {"pastas": pastas, "familia": "chuva", "valor": float(txt), "fs_def": "FS_min",
                                       "det": det}
    return esp


# ---------------------------------------------------------------------------------------------------------- leitura
def _f(s):
    """texto -> float; vazio -> None"""
    s = s.strip()
    return float(s) if s != "" else None


def _le_csv(arq):
    with open(arq, encoding="utf-8") as f:
        cab = f.readline().rstrip("\n").split(",")
        linhas = [l.rstrip("\n").split(",") for l in f if l.strip()]
    return cab, linhas


def _le_resumo(arq):
    cab, linhas = _le_csv(arq)
    if len(linhas) != 1:
        raise RuntimeError(f"{arq}: {len(linhas)} linhas de dados (esperada 1)")
    if len(linhas[0]) != len(cab):
        raise RuntimeError(f"{arq}: {len(linhas[0])} valores para {len(cab)} colunas")
    return dict(zip(cab, linhas[0]))


def _le_chave_valor(arq):
    """'a=1 b=2 ...' (resultado.txt, realizacao.txt): tokens sem '=' sao ignorados"""
    d = {}
    with open(arq, encoding="utf-8") as f:
        for tok in f.read().split():
            if "=" in tok:
                k, v = tok.split("=", 1)
                d[k] = v
    return d


def _le_parametros(arq):
    with open(arq, encoding="utf-8") as f:
        texto = f.read()
    d = {}
    for l in texto.splitlines():
        if "=" in l:
            k, v = l.split("=", 1)
            d[k.strip()] = v.strip()
    return texto, d


def _le_campanha(arq):
    d = {}
    if not os.path.exists(arq):
        return d
    with open(arq, encoding="utf-8") as f:
        for l in f:
            if "=" in l:
                k, v = l.strip().split("=", 1)
                d[k] = v.strip().strip('"')
    return d


def _instantes(param):
    """FS_H de parametros.txt; '-' = padrao do programa (0, 12, 24, 36 e 48 h)"""
    v = param.get("FS_H", "-")
    if v in ("", "-"):
        v = "0,12,24,36,48"
    return [float(x) for x in v.split(",") if x.strip()]


def _nome_fs(t):
    return f"FS_{t:g}h"


class Pasta:
    """uma combinacao de uma campanha: as realizacoes concluidas (resultado.txt nao vazio), as com falhou.txt e as sem
    nenhum dos dois"""

    def __init__(self, raiz, campanha, comb):
        self.campanha, self.comb = campanha, comb
        self.nome = f"{campanha}/{comb}"
        self.dir = os.path.join(raiz, campanha, comb)
        arq = os.path.join(self.dir, "parametros.txt")
        self.param_texto, self.param = _le_parametros(arq) if os.path.exists(arq) else ("", {})
        self.instantes = _instantes(self.param)
        self.dirs = {}
        for e in os.scandir(self.dir):
            if e.is_dir() and e.name.startswith("r") and e.name[1:].isdigit():
                self.dirs[int(e.name[1:])] = e.path
        self.concl, self.falhou, self.sem = [], [], []
        for n in sorted(self.dirs):
            d = self.dirs[n]
            res = os.path.join(d, "resultado.txt")
            if os.path.exists(res) and os.path.getsize(res) > 0:
                self.concl.append(n)
            elif os.path.exists(os.path.join(d, "falhou.txt")):
                self.falhou.append(n)
            else:
                self.sem.append(n)
        self._reais = None

    def parametros_nao_padrao(self):
        return {k: v for k, v in self.param.items() if v not in ("-", "") and k != "MESH"}

    def resumos(self):
        """{n: (resultado, resumo)} das realizacoes concluidas (leitura leve, para a proveniencia)"""
        out = {}
        for n in self.concl:
            d = self.dirs[n]
            res = _le_chave_valor(os.path.join(d, "resultado.txt"))
            arq = os.path.join(d, "resumo.csv")
            out[n] = (res, _le_resumo(arq) if os.path.exists(arq) else {})
        return out

    def realizacoes(self):
        """{n: dict} das realizacoes concluidas, conferidas (leitura completa: resultado.txt, resumo.csv e fs.csv)"""
        if self._reais is None:
            self._reais = {n: _le_realizacao(self.dirs[n], n, self.instantes) for n in self.concl}
            for r in self._reais.values():
                r["pasta"] = self.nome
        return self._reais


def _iguais(a, b, tol=1e-9):
    if a is None or b is None:
        return a is None and b is None
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def _linhas_fs(arq, instantes):
    """as linhas integras de fs.csv ({t: linha}) e o numero de linhas corrompidas. Integra: so caracteres ASCII
    imprimiveis, o numero de campos do cabecalho, t_h um dos instantes de FS_H e FS_SRM numerico."""
    with open(arq, "rb") as f:
        texto = f.read().decode("latin-1")
    linhas = texto.split("\n")
    cab = linhas[0].rstrip("\r").split(",")
    boas, ruins = {}, 0
    for l in linhas[1:]:
        l = l.rstrip("\r")
        if not l.strip():
            continue
        campos = l.split(",")
        ok = all(32 <= ord(ch) < 127 for ch in l) and len(campos) == len(cab)
        if ok:
            try:
                t, _ = float(campos[0]), float(campos[1])
            except ValueError:
                ok = False
        if not ok or t not in instantes:
            ruins += 1
            continue
        if t in boas:
            raise RuntimeError(f"{arq}: duas linhas com t_h = {t:g}")
        boas[t] = dict(zip(cab, campos))
    return boas, ruins


def _le_realizacao(d, n, instantes):
    """resultado.txt, resumo.csv e fs.csv de uma realizacao concluida, com as conferencias entre eles. Um instante sem
    linha integra em fs.csv (arquivo corrompido) usa o FS_<t>h de resumo.csv, e a realizacao fica marcada."""
    res = _le_chave_valor(os.path.join(d, "resultado.txt"))
    resumo = _le_resumo(os.path.join(d, "resumo.csv"))
    linhas, ruins = _linhas_fs(os.path.join(d, "fs.csv"), instantes)
    erro = lambda msg: RuntimeError(f"{d}: {msg}")  # noqa: E731
    if int(res.get("real", -1)) != n or int(resumo.get("real", -1)) != n:
        raise erro(f"numero da realizacao: pasta {n}, resultado {res.get('real')}, resumo {resumo.get('real')}")
    fs_t, fs_txt, de_resumo = {}, {}, []
    for t in sorted(instantes):
        col = _nome_fs(t)
        v_res = _f(resumo[col]) if col in resumo else None
        if t in linhas:
            fs_txt[t] = linhas[t]["FS_SRM"]
            fs_t[t] = _f(fs_txt[t])
            if col in resumo and not _iguais(v_res, fs_t[t]):
                raise erro(f"{col}: resumo.csv {resumo[col]} e fs.csv {fs_txt[t]}")
        elif v_res is not None:
            fs_t[t], fs_txt[t] = v_res, resumo[col]
            de_resumo.append(t)
        else:
            raise erro(f"sem FS em t = {t:g} h em fs.csv e em resumo.csv")
    # o minimo e o primeiro instante em que ocorre (como main.cpp: comparacao estrita, em ordem de tempo)
    t_min = None
    for t in sorted(instantes):
        if t_min is None or fs_t[t] < fs_t[t_min]:
            t_min = t
    fs_min = fs_t[t_min]
    for nome, v in (("resumo.csv FS_min", _f(resumo.get("FS_min", ""))),
                    ("resultado.txt FS_min", _f(res.get("FS_min", "")))):
        if not _iguais(v, fs_min):
            raise erro(f"{nome} = {v}, minimo dos instantes = {fs_min}")
    for nome, v in (("resumo.csv t_min_h", _f(resumo.get("t_min_h", ""))),
                    ("resultado.txt t_min_h", _f(res.get("t_min_h", "")))):
        if not _iguais(v, t_min):
            raise erro(f"{nome} = {v}, instante do minimo = {t_min}")
    # banda e massa de resumo.csv sao as do instante do minimo
    if t_min in linhas:
        for c in CAMPOS_AMOSTRA:
            if c in resumo and c in linhas[t_min] and not _iguais(_f(resumo[c]), _f(linhas[t_min][c])):
                raise erro(f"{c}: resumo.csv {resumo[c]} e fs.csv (t = {t_min:g} h) {linhas[t_min][c]}")
    situacao = res.get("situacao", "")
    if resumo.get("situacao", situacao) != situacao:
        raise erro(f"situacao: resultado.txt {situacao}, resumo.csv {resumo.get('situacao')}")
    return {"n": n, "dir": d, "res": res, "resumo": resumo, "linhas": linhas, "fs_t": fs_t, "fs_txt": fs_txt,
            "fs_min": fs_min, "t_min": t_min, "situacao": situacao, "fs_csv_linhas_corrompidas": ruins,
            "fs_de_resumo": de_resumo}


def _fs_do_caso(r, fs_def, instantes):
    """(FS, instante, {campo da amostra: valor} do instante, texto do FS) de uma realizacao, pela definicao do caso"""
    if fs_def == "FS_24h":
        t = 24.0
        if t not in r["fs_t"]:
            raise RuntimeError(f"{r['dir']}: sem FS em t = 24 h")
    else:
        t = None
        for tt in sorted(instantes):  # o primeiro minimo entre os instantes do caso
            if t is None or r["fs_t"][tt] < r["fs_t"][t]:
                t = tt
    if t in r["linhas"]:
        vals = {c: _f(r["linhas"][t].get(c, "")) for c in CAMPOS_AMOSTRA}
    elif t == r["t_min"]:  # linha corrompida em fs.csv: banda e massa do instante do minimo, de resumo.csv
        vals = {c: _f(r["resumo"].get(c, "")) for c in CAMPOS_AMOSTRA}
    else:
        raise RuntimeError(f"{r['dir']}: sem a banda e a massa de t = {t:g} h")
    return r["fs_t"][t], t, vals, r["fs_txt"][t]


# ------------------------------------------------------------------------------------------------------ estatistica
def _qui2(x, cdf, k):
    """qui-quadrado com k classes equiprovaveis sob a distribuicao ajustada: cada valor vai para a classe floor(k F(x))"""
    u = cdf(x)
    idx = np.minimum(np.floor(u * k).astype(int), k - 1)
    obs = np.bincount(idx, minlength=k)
    esp = len(x) / k
    chi2 = float(np.sum((obs - esp) ** 2) / esp)
    gl = k - 3
    return {"chi2": chi2, "gl": int(gl), "p": float(stats.chi2.sf(chi2, gl)), "observados": [int(o) for o in obs],
            "esperado": float(esp)}


H_SRM = 0.00625  # resolucao da bisseccao do SRM (FactorSRMRetry): quase todo FS e multiplo dela


def _qui2_grade(x, cdf, ppf, k, h=H_SRM):
    """variante de _qui2 para FS em grade: as bordas das classes sao os pontos medios da grade (multiplos impares de h/2)
    mais proximos dos quantis j/k da distribuicao ajustada (bordas repetidas se fundem), e a frequencia esperada de
    cada classe e n vezes a probabilidade ajustada dela; um valor sobre uma borda vai para a classe de cima"""
    bordas = sorted({round((math.floor(ppf(j / k) / h) + 0.5) * h, 10) for j in range(1, k)})
    prob = np.diff([0.0] + [float(cdf(b)) for b in bordas] + [1.0])
    idx = np.searchsorted(np.array(bordas), x, side="right")
    obs = np.bincount(idx, minlength=len(prob))
    esp = len(x) * prob
    chi2 = float(np.sum((obs - esp) ** 2 / esp))
    gl = len(prob) - 3
    return {"chi2": chi2, "gl": int(gl), "p": float(stats.chi2.sf(chi2, gl)), "classes": int(len(prob)),
            "bordas": [float(b) for b in bordas], "observados": [int(o) for o in obs],
            "esperados": [float(e) for e in esp]}


def estatisticas(fs):
    """as grandezas do Ceron para uma amostra de FS (FS = -1: sem equilibrio nem com c'/0,6; conta em P_f e fica fora
    dos ajustes e das estatisticas amostrais)"""
    fs = np.asarray(fs, dtype=float)
    N = int(len(fs))
    n_f = int(np.sum(fs < 1.0))
    Pf = n_f / N
    pos = fs[fs > 0]
    m = int(len(pos))
    media = float(np.mean(pos))
    dp = float(np.std(pos, ddof=1))
    amostra = {"n": m, "media": media, "dp": dp, "cov": dp / media, "mediana": float(np.median(pos)),
               "min": float(np.min(pos)), "max": float(np.max(pos)),
               "assimetria": float(stats.skew(pos, bias=False)),
               "assimetria_g1": float(stats.skew(pos, bias=True))}
    ln = np.log(pos)
    mu_ln, s_ln = float(np.mean(ln)), float(np.std(ln, ddof=0))
    m_ln = math.exp(mu_ln + 0.5 * s_ln ** 2)
    cov_ln = math.sqrt(math.expm1(s_ln ** 2))
    lognormal = {"mu_ln": mu_ln, "sigma_ln": s_ln, "media": m_ln, "mediana": math.exp(mu_ln), "dp": m_ln * cov_ln,
                 "cov": cov_ln}
    k = int(math.ceil(2.0 * m ** 0.4))
    mu_n, s_n = float(np.mean(pos)), float(np.std(pos, ddof=0))
    qn = _qui2(pos, lambda x: stats.norm.cdf(x, mu_n, s_n), k)
    qn.update({"mu": mu_n, "sigma": s_n})
    ql = _qui2(pos, lambda x: stats.norm.cdf(np.log(x), mu_ln, s_ln), k)
    menor = "lognormal" if ql["chi2"] < qn["chi2"] else ("normal" if qn["chi2"] < ql["chi2"] else "iguais")
    gn = _qui2_grade(pos, lambda x: stats.norm.cdf(x, mu_n, s_n), lambda u: stats.norm.ppf(u, mu_n, s_n), k)
    gln = _qui2_grade(pos, lambda x: stats.norm.cdf(np.log(x), mu_ln, s_ln),
                      lambda u: math.exp(stats.norm.ppf(u, mu_ln, s_ln)), k)
    menor_g = "lognormal" if gln["chi2"] < gn["chi2"] else ("normal" if gn["chi2"] < gln["chi2"] else "iguais")
    r = pos / H_SRM
    return {"N": N, "n_fs_negativo": int(np.sum(fs < 0)), "n_fs_igual_1": int(np.sum(fs == 1.0)), "n_f": n_f,
            "Pf": Pf, "cov_Pf": math.sqrt((1.0 - Pf) / (N * Pf)) if Pf > 0 else None,
            "n_fs_fora_grade": int(np.sum(np.abs(r - np.round(r)) > 1e-6)),
            "amostra": amostra, "lognormal": lognormal,
            "qui2": {"k": k, "n": m, "normal": qn, "lognormal": ql, "menor": menor},
            "qui2_grade": {"k": k, "n": m, "h": H_SRM, "normal": gn, "lognormal": gln, "menor": menor_g}}


def correlacoes(fs, colunas):
    """Spearman entre o FS e (c, phi, k) equivalentes; so as realizacoes com FS > 0 e as tres variaveis definidas"""
    X = np.column_stack([np.asarray(fs, float)] + [np.array([np.nan if v is None else v for v in colunas[nome]],
                                                            dtype=float) for nome, _ in VARS_SPEARMAN])
    ok = np.all(np.isfinite(X), axis=1) & (X[:, 0] > 0)
    X = X[ok]
    n = int(len(X))
    if n < 5:
        return None
    # uma variavel constante (k_mult = 1 em todo ponto com CoV(k_v) = 0) nao tem correlacao: fica de fora, e as
    # parciais das outras controlam so as variaveis nao constantes (ordem menor que 2)
    var = [j for j in range(1, X.shape[1]) if np.ptp(X[:, j]) > 0]
    constantes = [VARS_SPEARMAN[j - 1][0] for j in range(1, X.shape[1]) if j not in var]
    C = spearman(X[:, [0] + var])
    out = {"n": n, "segunda_ordem": {}, "ic95": {}, "simples": {}, "ic95_simples": {}, "ordem": {},
           "constantes": constantes, "colunas": {nome: col for nome, col in VARS_SPEARMAN}}
    for j, (nome, _) in enumerate(VARS_SPEARMAN):
        if j + 1 not in var:
            for chave in ("segunda_ordem", "ic95", "simples", "ic95_simples", "ordem"):
                out[chave][nome] = None
            continue
        i = var.index(j + 1) + 1
        q = len(var) - 1  # controles da parcial
        rp = float(parcial(C, 0, i))
        r1 = float(C[0, i])
        out["segunda_ordem"][nome] = rp
        out["ic95"][nome] = [float(v) for v in ic95(rp, n, q)]
        out["ordem"][nome] = q
        out["simples"][nome] = r1
        out["ic95_simples"][nome] = [float(v) for v in ic95(r1, n, 0)]
    return out


# ----------------------------------------------------------------------------------------------------------- casos
def _constantes(reais, nome_caso):
    """os campos de identificacao de resumo.csv, que devem ser os mesmos em todas as realizacoes do caso"""
    vals = {}
    for r in reais:
        for c in CAMPOS_IDENT:
            if c in r["resumo"]:
                vals.setdefault(c, set()).add(r["resumo"][c])
    out = {}
    for c, s in vals.items():
        if len(s) > 1:
            raise RuntimeError(f"caso {nome_caso}: {c} varia entre as realizacoes: {sorted(s)}")
        out[c] = next(iter(s))
    return out


def _confere_caso(cid, e, ident):
    """os parametros lidos sao os do mapa de casos"""
    fam, v = e["familia"], e["valor"]
    esperado = {"cov_c": 0.3, "cov_phi": 0.1, "cov_k": 0.6, "Lx": L_REF[0], "Ly": L_REF[1], "alfa": 1.0,
                "I_mm_h": 36.0}
    if fam in ("covk", "covc", "covphi", "alfa"):
        esperado[CAMPO_FAMILIA[fam]] = v
    elif fam == "s":
        esperado["Lx"], esperado["Ly"] = v * L_REF[0], v * L_REF[1]
    elif fam == "chuva":
        esperado["I_mm_h"] = v
    for c, x in esperado.items():
        if not _iguais(float(ident[c]), x):
            raise RuntimeError(f"caso {cid}: {c} = {ident[c]} nos dados, {x} no mapa de casos")
    if ident.get("campos") != "1":
        raise RuntimeError(f"caso {cid}: campos = {ident.get('campos')} (esperado 1)")
    if ident.get("solo") != cid.split("_")[0]:
        raise RuntimeError(f"caso {cid}: solo {ident.get('solo')}")


def _deterministico(pasta, fs_def, instantes, ident_caso):
    if pasta.concl != [0]:
        raise RuntimeError(f"{pasta.nome}: realizacoes concluidas {pasta.concl} (esperada so a 0)")
    r = pasta.realizacoes()[0]
    if r["resumo"].get("campos") != "0":
        raise RuntimeError(f"{pasta.nome}: campos = {r['resumo'].get('campos')} (esperado 0)")
    if r["situacao"] != "ok":
        raise RuntimeError(f"{pasta.nome}: situacao {r['situacao']}")
    for c in ("solo", "I_mm_h", "variante", "f_trans", "alfa"):
        if c in r["resumo"] and c in ident_caso and r["resumo"][c] != ident_caso[c]:
            raise RuntimeError(f"{pasta.nome}: {c} = {r['resumo'][c]}, no caso {ident_caso[c]}")
    for t in instantes:
        if t not in r["fs_t"]:
            raise RuntimeError(f"{pasta.nome}: sem FS em t = {t:g} h")
    fs, t, _, txt = _fs_do_caso(r, fs_def, instantes)
    return {"pasta": pasta.nome, "FS": fs, "FS_txt": txt, "t_h": t, "versao": r["resumo"].get("versao", ""),
            "serie": [[t, r["fs_t"][t], r["fs_txt"][t]] for t in sorted(r["fs_t"])],
            "instantes_h": sorted(r["fs_t"])}


def _texto_fs_def(fs_def, instantes):
    if fs_def == "FS_24h":
        return "FS em t = 24 h"
    return "FS mínimo nos instantes " + _lista_pt([_num_pt(t) for t in instantes]) + " h"


def _caso(cid, e, pastas, solo):
    ps = [pastas[p] for p in e["pastas"]]
    # fusao: parametros.txt identico e realizacoes disjuntas
    fusao = None
    if len(ps) > 1:
        base = ps[0]
        for p in ps[1:]:
            if p.param != base.param or p.param_texto.strip() != base.param_texto.strip():
                raise RuntimeError(f"caso {cid}: parametros.txt de {p.nome} diferente do de {base.nome}; nao junto")
        vistos = {}
        for p in ps:
            for n in p.dirs:
                if n in vistos:
                    raise RuntimeError(f"caso {cid}: realizacao {n} em {vistos[n]} e em {p.nome}; nao junto")
                vistos[n] = p.nome
        fusao = {"caso": cid, "pastas": [p.nome for p in ps], "parametros_identicos": True, "reais_disjuntos": True,
                 "N_por_pasta": [len(p.concl) for p in ps],
                 "faixas": [[min(p.dirs), max(p.dirs)] for p in ps]}
    instantes = ps[0].instantes
    reais = []
    for p in ps:
        reais.extend(p.realizacoes().values())
    reais.sort(key=lambda r: r["n"])
    ident = _constantes(reais, cid)
    _confere_caso(cid, e, ident)
    nao_ok = [r["n"] for r in reais if r["situacao"] != "ok"]
    usadas = [r for r in reais if r["situacao"] == "ok"]
    fs, ts, amos = [], [], {c: [] for c in CAMPOS_AMOSTRA}
    versoes = {}
    for r in usadas:
        v, t, vals, _ = _fs_do_caso(r, e["fs_def"], instantes)
        fs.append(v)
        ts.append(t)
        for c in CAMPOS_AMOSTRA:
            amos[c].append(vals[c])
        ver = r["resumo"].get("versao", "")
        versoes[ver] = versoes.get(ver, 0) + 1
    est = estatisticas(fs)
    tcont = {}
    for t in ts:
        tcont[t] = tcont.get(t, 0) + 1
    fam, v = e["familia"], e["valor"]
    s = float(ident["Lx"]) / L_REF[0]
    caso = {
        "id": cid, "solo": solo, "pastas": [p.nome for p in ps], "familia": fam, "valor": v,
        "valores": ({f: VALOR_REF[f] for f in ("covk", "covc", "covphi", "s", "alfa")} if fam == "ref"
                    else {fam: v}),
        "rotulo": _rotulo(fam, v),
        "cov_c": float(ident["cov_c"]), "cov_phi": float(ident["cov_phi"]), "cov_k": float(ident["cov_k"]),
        "Lx": float(ident["Lx"]), "Ly": float(ident["Ly"]), "s": s, "alfa": float(ident["alfa"]),
        "I_mm_h": float(ident["I_mm_h"]), "variante": ident.get("variante"), "f_trans": float(ident["f_trans"]),
        "kra_min": float(ident["kra_min"]) if "kra_min" in ident else None,
        "KL_M": int(ident["KL_M"]), "KL_variancia": float(ident["KL_variancia"]),
        "semente": int(ident["semente"]),
        "fs_def": e["fs_def"], "fs_def_texto": _texto_fs_def(e["fs_def"], instantes), "instantes_h": instantes,
        "N": est["N"], "n_falhou": sum(len(p.falhou) for p in ps),
        "falhou": [{"pasta": p.nome, "reais": list(p.falhou)} for p in ps if p.falhou],
        "n_sem_resultado": sum(len(p.sem) for p in ps), "n_situacao_nao_ok": len(nao_ok),
        "reais_situacao_nao_ok": nao_ok,
        "fs_csv_corrompido": [{"pasta": r["pasta"], "real": r["n"], "linhas_corrompidas": r["fs_csv_linhas_corrompidas"],
                               "instantes_de_resumo": r["fs_de_resumo"]}
                              for r in reais if r["fs_csv_linhas_corrompidas"] or r["fs_de_resumo"]],
        "n_fs_negativo": est["n_fs_negativo"], "n_fs_igual_1": est["n_fs_igual_1"], "n_f": est["n_f"],
        "Pf": est["Pf"], "cov_Pf": est["cov_Pf"],
        "n_fs_fora_grade": est["n_fs_fora_grade"],
        "amostra": est["amostra"], "lognormal": est["lognormal"], "qui2": est["qui2"], "qui2_grade": est["qui2_grade"],
        "det": None, "versoes": versoes,
        "t_min_h_contagem": [[t, tcont[t]] for t in sorted(tcont)],
        "faixa_reais": [usadas[0]["n"], usadas[-1]["n"]], "fusao": fusao,
        "spearman": correlacoes(fs, {nome: amos[col] for nome, col in VARS_SPEARMAN}),
        "amostras": {"real": [r["n"] for r in usadas], "FS": fs, "t_min_h": ts, **amos},
    }
    if e["det"] is not None:
        caso["det"] = _deterministico(pastas[e["det"]], e["fs_def"], instantes, ident)
    return caso


# ------------------------------------------------------------------------------------------------------ proveniencia
def _proveniencia_pasta(p, uso):
    res = p.resumos()
    versoes, situacoes, f_trans, kra = {}, {}, set(), set()
    sem_fs = []
    for n, (rt, rs) in res.items():
        ver = rs.get("versao", "")
        versoes[ver] = versoes.get(ver, 0) + 1
        sit = rt.get("situacao", "")
        situacoes[sit] = situacoes.get(sit, 0) + 1
        if "f_trans" in rs:
            f_trans.add(rs["f_trans"])
        if "kra_min" in rs:
            kra.add(rs["kra_min"])
        if rt.get("FS_min", "") == "":
            sem_fs.append(n)
    d = {"pasta": p.nome, "campanha": p.campanha, "combinacao": p.comb,
         "parametros": p.parametros_nao_padrao(), "instantes_h": p.instantes,
         "n_pastas_r": len(p.dirs), "faixa_reais": [min(p.dirs), max(p.dirs)] if p.dirs else None,
         "concluidas": len(p.concl), "falhou": list(p.falhou), "sem_resultado": list(p.sem),
         "concluidas_sem_FS": sem_fs, "versoes": versoes, "situacao": situacoes,
         "f_trans": sorted(f_trans), "kra_min": sorted(kra), "deterministica": p.comb.endswith("_det"),
         "uso": sorted(uso.get((p.campanha, p.comb), [])), "excluida": EXCLUIDAS.get((p.campanha, p.comb))}
    if d["excluida"] is None and not d["uso"]:
        d["excluida"] = "não prevista no mapa de casos do relatório"
        print(f"aviso: {p.nome} nao esta no mapa de casos nem na lista de excluidas", file=sys.stderr)
    return d


def _compara_f50(pastas):
    """sensibilidade/f50_I36_i (referencia f50) contra monte_carlo/f50_I36_i_ftrans0.5 nas realizacoes em comum"""
    a, b = pastas.get(("sensibilidade", "f50_I36_i")), pastas.get(("monte_carlo", "f50_I36_i_ftrans0.5"))
    if a is None or b is None:
        return None
    ra, rb = a.realizacoes(), b.realizacoes()
    comuns = sorted(set(ra) & set(rb))
    if not comuns:
        return None
    iguais = {"FS_24h": 0, "c_med": 0, "phi_med": 0, "k_mult_med": 0, "banda_c": 0, "massa_k_mult": 0}
    dif_max = 0.0
    for n in comuns:
        x, y = ra[n], rb[n]
        if _iguais(x["fs_t"][24.0], y["fs_t"][24.0], 0.0):
            iguais["FS_24h"] += 1
        dif_max = max(dif_max, abs(x["fs_t"][24.0] - y["fs_t"][24.0]))
        for c in ("c_med", "phi_med", "k_mult_med", "banda_c", "massa_k_mult"):
            if x["resumo"].get(c) == y["resumo"].get(c):
                iguais[c] += 1
    va = {r["resumo"].get("versao", "") for r in ra.values()}
    vb = {r["resumo"].get("versao", "") for r in rb.values()}
    return {"pastas": [a.nome, b.nome], "versoes": [sorted(va), sorted(vb)], "reais_comuns": len(comuns),
            "faixa_comum": [comuns[0], comuns[-1]], "iguais": iguais, "dif_max_FS_24h": dif_max,
            "descricao": ("realizações em comum das duas pastas: quantas têm o mesmo valor (texto idêntico em "
                          "resumo.csv; FS em t = 24 h pelo valor) de cada grandeza")}


def _limpa(x):
    """numpy -> python, NaN/inf -> None, chaves -> texto, tuplas -> listas"""
    if isinstance(x, dict):
        return {str(k): _limpa(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_limpa(v) for v in x]
    if isinstance(x, np.ndarray):
        return [_limpa(v) for v in x.tolist()]
    if isinstance(x, (np.bool_, bool)):
        return bool(x)
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (float, np.floating)):
        x = float(x)
        return x if math.isfinite(x) else None
    return x


METODOS = {
    "S": "S = número de realizações concluídas (resultado.txt não vazio) com situacao = ok",
    "Pf": "P_f = n_f/S, n_f = número de realizações com FS < 1 (Ceron, eqs. 67 e 69); FS = −1 (sem equilíbrio nem "
          "com c'/0,6 e tan φ'/0,6) conta como FS < 1",
    "cov_Pf": "CoV(P_f) = √((1 − P_f)/(S P_f)) (Ceron, eq. 68)",
    "lognormal": "log-normal de 2 parâmetros por máxima verossimilhança sobre as realizações com FS > 0: "
                 "μ_ln = média de ln FS, σ_ln = desvio-padrão populacional de ln FS; média = exp(μ_ln + σ_ln²/2), "
                 "mediana = exp(μ_ln), CoV = √(exp(σ_ln²) − 1), σ = média × CoV",
    "amostra": "realizações com FS > 0: média, desvio-padrão amostral (divisor n − 1), CoV = desvio/média, mediana, "
               "mínimo, máximo e coeficiente de assimetria amostral ajustado de Fisher-Pearson G1 "
               "(assimetria_g1: o não ajustado)",
    "qui2": "qui-quadrado de aderência das distribuições normal (μ e σ de máxima verossimilhança) e log-normal "
            "ajustadas: k = ⌈2 n^0,4⌉ classes de igual probabilidade sob cada distribuição ajustada (cada FS vai "
            "para a classe ⌊k F(FS)⌋), n = realizações com FS > 0, frequência esperada n/k, gl = k − 3",
    "qui2_grade": "variante do qui-quadrado para FS em grade: as bordas das classes são os pontos médios da grade do "
                  "FS (múltiplos ímpares de 0,003125) mais próximos dos quantis j/k da distribuição ajustada; a "
                  "frequência esperada de cada classe é n vezes a probabilidade ajustada dela; gl = classes − 3",
    "spearman": "correlação de postos de Spearman (postos médios nos empates) entre o FS e as variáveis "
                "equivalentes de cada realização, no instante do FS: c = c' médio na banda de ruptura (banda_c), "
                "φ = φ' médio na banda (banda_phi_deg), k = k/k_s médio na massa que desliza (massa_k_mult); "
                "parcial de 2ª ordem ρ_XFS,YZ = −P_ij/√(P_ii P_jj), P = inversa da matriz 4×4 de correlação de "
                "postos (equivalente às eqs. 70-71 do Ceron); intervalo de 95% pela transformação de Fisher com "
                "erro-padrão 1/√(n − 3 − q), q = 2 controles (q = 0 na correlação simples); só realizações com FS > 0",
    "banda": "banda de ruptura: pontos de integração com √J2 da deformação plástica ≥ 0,1 × máximo, no último estado "
             "convergido do SRM do instante do FS",
    "massa": "massa que desliza: pontos com |Δu| ≥ 0,5 × máximo, Δu o incremento de deslocamento entre os dois "
             "últimos estados convergidos do SRM do instante do FS",
    "FS": "fator de segurança por redução de resistência (SRM, bissecção com resolução 0,00625: os FS são múltiplos "
          "de 0,00625); FS = −1: sem equilíbrio nem com c'/0,6",
    "FS_24h": "FS do SRM em t = 24 h (fim da chuva)",
    "FS_min": "o menor FS do SRM entre os instantes calculados da realização (o primeiro, em caso de empate)",
    "convergencia": "P_f(n) = fração de FS < 1 entre as n primeiras realizações, em ordem crescente do número da "
                    "realização",
}


def calcula(raiz_dados):
    raiz = os.path.abspath(raiz_dados)
    esp = especificacao()
    # todas as pastas de todas as campanhas
    campanhas, pastas = {}, {}
    for camp in sorted(os.listdir(raiz)):
        cd = os.path.join(raiz, camp)
        if not os.path.isdir(cd):
            continue
        campanhas[camp] = _le_campanha(os.path.join(cd, "campanha.txt"))
        for comb in sorted(os.listdir(cd)):
            if os.path.isdir(os.path.join(cd, comb)):
                pastas[(camp, comb)] = Pasta(raiz, camp, comb)
    uso = {}
    for solo, cs in esp.items():
        for cid, e in cs.items():
            for p in e["pastas"] + ([e["det"]] if e["det"] else []):
                if p not in pastas:
                    raise RuntimeError(f"caso {cid}: pasta {p[0]}/{p[1]} nao existe em {raiz}")
                uso.setdefault(p, set()).add(cid)
    solos, fusoes = {}, []
    for solo, cs in esp.items():
        casos = {cid: _caso(cid, e, pastas, solo) for cid, e in cs.items()}
        for c in casos.values():
            if c["fusao"] and not any(f["pastas"] == c["fusao"]["pastas"] for f in fusoes):
                fusoes.append(c["fusao"])
        fams = {f: [] for f in FAMILIAS}
        ref = f"{solo}_ref"
        fams["ref"] = [ref]
        for f in ("covk", "covc", "covphi", "s", "alfa"):
            ids = [cid for cid, c in casos.items() if c["familia"] == f] + [ref]
            fams[f] = sorted(ids, key=lambda i: casos[i][CAMPO_FAMILIA[f]])
        fams["chuva"] = sorted([cid for cid, c in casos.items() if c["familia"] == "chuva"],
                               key=lambda i: casos[i]["I_mm_h"])
        ref_chuva = f"{solo}_I36"
        solos[solo] = {"k_s_m_s": K_S[solo], "k_s_ceron": K_S[solo] / GAMMA_W, "fs_def_ref": casos[ref]["fs_def"],
                       "casos": casos, "familias": fams, "ref": ref, "ref_chuva": ref_chuva,
                       "campo_familia": dict(CAMPO_FAMILIA),
                       "mesmos_dados_ref_chuva": casos[ref]["pastas"] == casos[ref_chuva]["pastas"]}
    prov = {
        "raiz_dados": os.path.basename(raiz),
        "campanhas": campanhas,
        "pastas": [_proveniencia_pasta(p, uso) for p in pastas.values()],
        "fusoes": fusoes,
        "excluidas": [{"pasta": f"{c}/{p}", "motivo": m} for (c, p), m in EXCLUIDAS.items() if (c, p) in pastas],
        "comparacao_f50_ref": _compara_f50(pastas),
        "metodos": METODOS,
        "constantes": {"k_s_m_s": K_S, "gamma_w_kN_m3": GAMMA_W, "gamma_kN_m3": GAMMA_SOLO, "c_medio_kPa": C_MEDIO,
                       "phi_medio_graus": PHI_MEDIO, "porosidade": POROSIDADE, "E_MPa": E_MPA, "nu": NU,
                       "chuva_duracao_h": TCHUVA_H, "L_ref_m": list(L_REF), "fonte": "main.cpp (valores padrão)"},
    }
    return _limpa({"solos": solos, "proveniencia": prov})


def grava(est, saida):
    os.makedirs(saida, exist_ok=True)
    arq = os.path.join(saida, "estatisticas.json")
    tmp = arq + ".tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(est, f, ensure_ascii=False, separators=(",", ":"), allow_nan=False)
    os.replace(tmp, arq)
    return arq


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("uso: python3 -I dados.py <raiz_dados> <saida>")
        sys.exit(2)
    print(grava(calcula(sys.argv[1]), sys.argv[2]))
