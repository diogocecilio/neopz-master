#!/usr/bin/env python3
"""coluna_ref.py: a referencia 1D em Python da coluna de chuva de Cho (2014) para o SlopeRainColumn (coluna.cpp): os
volumes finitos do col1d.py (copia literal do Apendice do cho2014.md, como o mc2.py; aqui so importados), forma mista
(ou forma=h), Newton, troca fluxo <-> face no topo como o RainTransient, e o FS de Lu-Godt de Cho (2014, Eqs. 5, 7 e 9)
nos pontos do perfil.

  python3 -I coluna_ref.py det [opcoes]                            o deterministico (k_s medio, ou kmult=)
  python3 -I coluna_ref.py mc nystrom=<base.txt> [opcoes]          Monte Carlo com a base de Nystrom 1D
  python3 -I coluna_ref.py xi nystrom=<base.txt> <xi...> [opcoes]  as realizacoes dos arquivos de xi do driver C++
  python3 -I coluna_ref.py fs <perfil.csv> [opcoes]                o FS de um perfil de p_w (V1.3)
  python3 -I coluna_ref.py tabela [eps=1e-3] [fecho=reta] [opcoes]  a tabela de k_r do vg_tab (T1)

Opcoes (chave=valor; os padroes sao o Cho 2014 e, onde ha, as opcoes e os nomes do driver C++, spec B.5):
  corrida   n=10 primeira=0 semente=2025 saida=coluna tempos=6,12,18,24 perfis=0 grava_xi=0 kmult=1 localiza=1
            gerador=pz|numpy desloc_xi=<arquivo>
  malha     H=2 N=200 (celulas de H/N) nos=<arquivo com as profundidades z, m> (os pontos do FS: os nos do NeoPZ)
  solo      vg=0.358,3.58e-4,0.149,1.289 (theta_s, theta_r, alfa em 1/kPa, n; S_wr = theta_r/theta_s) ks=3e-6
            Kw=2.2e6 (S_s = theta_s gw/Kw, constante, como o col1d) gw=9.81
            k_r: vg_tab=1e-3,reta (o padrao: o MualemTabelado.h do NeoPZ, fecho reta|c1|ippisch)
                 | vg_tab_psi=psi_c[,fecho] (eps = 1 - S_e(psi_c); 0.01,ippisch e o Ippisch pela tabela, so no k_r)
                 | vg_tab=0 (VG-Mualem pura: o Newton para em ~4,6 h, V2.0)
                 | ippisch=psi_e (o Solo do col1d: retencao S_e/S_c e k_r normalizado de Ippisch et al. 2006)
                 | col1d_tab=eps (o SoloTab do col1d: 200 pontos em 1 - S_e)
  chuva     I=10.8 (mm/h) chuva=24 (h; depois, fluxo 0) succao0=20 (kPa) nivel0=<z_w, m: p_w = gw (z_w - y)> succao_max=1e9
  numerica  dt=0.1 (h) dt_min=1s (h ou Ns; o passo e dividido ao meio ate dt_min) newton=1e-6,40 (|dp| em kPa, iteracoes)
            forma=mista|h (h: armazenamento C(h) (h - h_n), a forma em h da Eq. 1 de Cho, nao conservativa; Newton sem
            o termo C'(h)(h - h_n))
  FS        fs=talude|nenhum beta=40 c=5 phi=32 gamma_d=16 fs_peso=integrado|constante gamma=20 fs_zmin=0
            fs_em=nos|centros (nos: a superficie, os centros e as faces das celulas e a base, z = j dz/2, com a carga
            interpolada e extrapolada pelos contornos; centros: o fs_perfil do col1d, so nos centros)
  campo     covk=0.5 M=100 nystrom=<base.txt> (lv sai do cabecalho da base; lv= so confere)

det: imprime FS_min, z_crit, a frente, p_w na base e no topo e a menor carga de pressao em cada instante, os passos, as
iteracoes de Newton e o tempo; grava <saida>_0_-1.csv (real = -1) e <saida>_perfil_-1.csv. mc e xi: grava
<saida>_<a>_<b>.csv (a = primeira, b = a + n - 1) com as colunas do SlopeRainColumn (spec B.10; o que nao existe aqui
vai com o valor 'vf', 0 ou -1), <saida>_parametros.txt e, se pedidos, <saida>_perfil_<real>.csv (as perfis primeiras) e
<saida>_xi_<real>.txt (grava_xi=1: os M xi, um por linha, o formato do driver). O FS e calculado depois de cada passo:
FS_passagem e o menor ate o instante e t_passagem a primeira travessia de FS = 1 (localiza=1: bissecao no tempo, 7
vezes, ate dt/128, e interpolacao linear no ultimo intervalo; o estado volta ao fim do passo, na grade de dt). Se o
transiente para, os instantes seguintes ficam com situacao=transiente e t_parou.

O campo: ln k_s = mu + sigma Z com sigma = sqrt(ln(1 + covk^2)) e mu = -sigma^2/2 (k_mult de media 1), Z(y) =
sum_k sqrt(lambda_k) xi_k f_k(y) com f_k interpolado por Nystrom nos centros das celulas (o RandomField::At). gerador=pz
(o padrao): xi da realizacao n = klsep::Normal(semente * 0x9E3779B97F4A7C15 + n).Vetor(M) (KLSeparavel.h), os mesmos xi
do RandomField::Sample(semente, n) com um campo so (o k da coluna, spec B.8): o MC daqui fica pareado com o do NeoPZ
com a mesma semente. gerador=numpy: numpy.random.default_rng([semente, n]).

fs <perfil.csv>: le as colunas z (ou y, com H=) e pw (kPa), e t_h e FS se houver; calcula o FS de cada instante nos
pontos dados (o mesmo FS do det) e, com a coluna FS, a maior diferenca. tabela: psi_c, k_c, os pontos, o erro maximo de
k_r em [psi_c, 1e3 kPa], k_r em 0,01/0,1/1/10 kPa contra Mualem e o maior |dk_r/dpsi| (spec A.1); grava=<arq> grava a
tabela (S_w, ln k_r)."""
import csv
import math
import os
import re
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np  # noqa: E402

import col1d  # noqa: E402
from col1d import Coluna, Solo, SoloTab, fs_perfil  # noqa: E402
from scipy.linalg import solve_banded  # noqa: E402

COLUNAS = ("real,t_h,FS_min,situacao,z_crit,base,FS_passagem,t_passagem,frente_m,pw_base,pw_topo,nos_face,passos,newton,"
           "s_transiente,t_parou,k_mult_med,k_mult_min,k_mult_max,log_peso,ny,dt,elementos,massa,lei_k,cov_k,lv,M,"
           "KL_variancia")

PADRAO = {"n": "10", "primeira": "0", "semente": "2025", "saida": "coluna", "tempos": "6,12,18,24", "perfis": "0",
          "grava_xi": "0", "kmult": "1", "localiza": "1", "gerador": "pz", "desloc_xi": "",
          "H": "2", "N": "200", "nos": "",
          "vg": "0.358,3.58e-4,0.149,1.289", "ks": "3e-6", "Kw": "2.2e6", "gw": "9.81",
          "vg_tab": "1e-3,reta", "vg_tab_psi": "", "ippisch": "", "col1d_tab": "",
          "I": "10.8", "chuva": "24", "succao0": "20", "nivel0": "", "succao_max": "1e9",
          "dt": "0.1", "dt_min": "1s", "newton": "1e-6,40", "forma": "mista",
          "fs": "talude", "beta": "40", "c": "5", "phi": "32", "gamma_d": "16", "fs_peso": "integrado", "gamma": "20",
          "fs_zmin": "0", "fs_em": "nos",
          "covk": "0.5", "M": "100", "nystrom": "", "lv": "",
          "eps": "", "fecho": "", "grava": ""}
G_SEMENTE = 0x9E3779B97F4A7C15
# o que para o transiente (o col1d levanta RuntimeError quando o passo cai abaixo de dt_min; o solve_banded, LinAlgError
# com a matriz singular)
PARADAS = (RuntimeError, np.linalg.LinAlgError, FloatingPointError)
MASK = (1 << 64) - 1


def erro(msg):
    sys.exit(f"coluna_ref.py: {msg}")


# ---- o gerador do RandomField (KLSeparavel.h, klsep::Normal): xoshiro256** semeado por splitmix64, Box-Muller
class NormalPZ:
    def __init__(s, semente):
        x = semente & MASK
        s.s = []
        for _ in range(4):
            x = (x + 0x9E3779B97F4A7C15) & MASK
            z = x
            z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & MASK
            z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & MASK
            s.s.append(z ^ (z >> 31))
        s.tem, s.guarda = False, 0.

    def proximo(s):
        s0, s1, s2, s3 = s.s
        r1 = (s1 * 5) & MASK
        r = ((((r1 << 7) | (r1 >> 57)) & MASK) * 9) & MASK
        t = (s1 << 17) & MASK
        s2 ^= s0
        s3 ^= s1
        s1 ^= s2
        s0 ^= s3
        s2 ^= t
        s3 = ((s3 << 45) | (s3 >> 19)) & MASK
        s.s = [s0, s1, s2, s3]
        return r

    def uniforme(s):  # em (0, 1)
        return ((s.proximo() >> 11) + 0.5) * (1.0 / 9007199254740992.0)

    def __call__(s):
        if s.tem:
            s.tem = False
            return s.guarda
        u1, u2 = s.uniforme(), s.uniforme()
        r, t = math.sqrt(-2. * math.log(u1)), 2. * math.pi * u2
        s.guarda, s.tem = r * math.sin(t), True
        return r * math.cos(t)

    def vetor(s, n):
        return np.array([s() for _ in range(n)])


def xi_realizacao(o, n, M):
    """os M xi da realizacao n (RandomField::Sample com um campo, ou o numpy)"""
    if o["gerador"] == "pz":
        return NormalPZ((int(o["semente"]) * G_SEMENTE + int(n)) & MASK).vetor(M)
    return np.random.default_rng([int(o["semente"]), int(n)]).standard_normal(M)


# ---- a tabela de k_r de Mualem do vg_tab (spec A.3, o MualemTabelado.h): nos em x = alfa psi, so depende de n, S_wr e eps
def ln_kr_mualem(x, n):
    """ln k_r de Mualem em x = alfa psi (y = x^n): -(m/2) log1p(y) + 2 log(-expm1(-m log1p(1/y)))"""
    m = 1. - 1. / n
    y = np.asarray(x, dtype=float) ** n
    return -(m / 2.) * np.log1p(y) + 2. * np.log(-np.expm1(-m * np.log1p(1. / y)))


def se_vg(x, n):
    return (1. + np.asarray(x, dtype=float) ** n) ** (-(1. - 1. / n))


def eps_de_psi(a, n, psic):
    """1 - S_e(psi_c), a em 1/kPa"""
    m = 1. - 1. / n
    return -math.expm1(-m * math.log1p((a * psic) ** n))


def tabela_vg(n, swr, eps, fecho="reta", N=400):
    """(S_w, ln k_r, info) da tabela: N + 1 valores de x log-espacados de 1e6 a x_c = ((1 - eps)^(-1/m) - 1)^(1/n), so os
    S_w estritamente crescentes, mais o fecho ate S_w = 1: reta (1, 0); c1 (ln k = L_c + g_c t + A t^2 em 50 pontos, com
    g_c = d ln k_r/dS_w analitico em S_c e A = (-L_c - g_c e')/e'^2, e' = 1 - S_c); ippisch (L_c subtraido de tudo e
    (1, 0))"""
    if not 1e-5 <= eps <= 0.1:
        erro(f"vg_tab: eps = {eps:g} fora de [1e-5, 0.1]")
    if fecho not in ("reta", "c1", "ippisch"):
        erro(f"vg_tab: fecho {fecho}: reta, c1 ou ippisch")
    m = 1. - 1. / n
    xc = math.expm1(-math.log1p(-eps) / m) ** (1. / n)
    x = np.exp(np.linspace(math.log(1e6), math.log(xc), N + 1))
    x[0], x[-1] = 1e6, xc
    y = x ** n
    Sw = swr + (1. - swr) * (1. + y) ** (-m)
    lnk = -(m / 2.) * np.log1p(y) + 2. * np.log(-np.expm1(-m * np.log1p(1. / y)))
    S, L = [Sw[0]], [lnk[0]]
    for s, v in zip(Sw[1:], lnk[1:]):
        if s > S[-1]:
            S.append(s)
            L.append(v)
    S, L = np.array(S), np.array(L)
    Sc, Lc = S[-1], L[-1]
    info = {"xc": xc, "Sc": Sc, "kc": math.exp(Lc), "tabela": len(S)}
    if fecho == "reta":
        S, L = np.append(S, 1.), np.append(L, 0.)
    elif fecho == "ippisch":
        S, L = np.append(S, 1.), np.append(L - Lc, 0.)
    else:
        yc = xc ** n
        Sec = (1. + yc) ** (-m)
        q = 1. / (1. + yc)  # S_e^(1/m)
        f = -math.expm1(-m * math.log1p(1. / yc))  # 1 - (1 - q)^m
        df = math.exp(-(m - 1.) * math.log1p(1. / yc)) * q / Sec  # (1 - q)^(m-1) q / S_e
        gc = (0.5 / Sec + 2. * df / f) / (1. - swr)
        ep = 1. - Sc
        A = (-Lc - gc * ep) / ep ** 2
        t = ep * np.arange(1, 51) / 50.
        Se_, Le_ = Sc + t, Lc + gc * t + A * t * t
        Se_[-1], Le_[-1] = 1., 0.
        S, L = np.append(S, Se_), np.append(L, Le_)
        info["gc"] = gc
    return S, L, info


def kr_tabela(S, St, Lt, kmin=1e-12):
    """k_r e dk_r/dS_w como o ETabulated do TPZUnsaturatedHydraulics.h (ln k linear entre os pontos; constante fora,
    derivada 0; piso fKrwMin = 1e-12)"""
    S = np.asarray(S, dtype=float)
    dentro = (S > St[0]) & (S < St[-1])
    i = np.clip(np.searchsorted(St, S, side="right") - 1, 0, len(St) - 2)
    g = (Lt[i + 1] - Lt[i]) / (St[i + 1] - St[i])
    k = np.where(dentro, np.exp(Lt[i] + g * (S - St[i])), np.exp(np.where(S <= St[0], Lt[0], Lt[-1])))
    dk = np.where(dentro, k * g, 0.)
    baixo = ~(k > kmin)
    return np.where(baixo, kmin, k), np.where(baixo, 0., dk)


class SoloTabVG(Solo):
    """retencao VG pura e k_r da tabela do vg_tab (o material do NeoPZ com o MualemTabelado.h)"""

    def __init__(s, eps=1e-3, fecho="reta", **kw):
        super().__init__(**kw)
        s.swr = s.tr / s.ts
        s.eps, s.fecho = eps, fecho
        s.St, s.Lt, s.info = tabela_vg(s.n, s.swr, eps, fecho)

    def Sw(s, h):
        return s.swr + (1. - s.swr) * s.Se(h)

    def kr(s, h):
        return kr_tabela(s.Sw(h), s.St, s.Lt)[0]

    def dkr(s, h):  # dk_r/dh = (dk_r/dS_w) dS_w/dh, dS_w/dh = (1 - S_wr) dS_e/dh = C/theta_s
        return kr_tabela(s.Sw(h), s.St, s.Lt)[1] * s.C(h) / s.ts


# ---- a coluna: o Coluna do col1d com a contagem dos passos aceitos e a forma em h
class ColunaRef(Coluna):
    forma = "mista"

    def __init__(s, *a, **kw):
        super().__init__(*a, **kw)
        s.stats["passos"] = 0

    def step(s, dt, metodo="newton", tol=1e-6, maxit=40):
        ok, its = super().step(dt, metodo, tol, maxit)
        if ok:
            s.stats["passos"] += 1
        return ok, its

    def solve(s, hn, dt, face, metodo, tol, maxit):
        if s.forma == "mista":
            return super().solve(hn, dt, face, metodo, tol, maxit)
        # forma em h: o mesmo solve do col1d com o armazenamento C(h) (h - h_n) no lugar de theta(h) - theta(h_n)
        so, dz, N = s.so, s.dz, s.N
        h = hn.copy()
        for it in range(1, maxit + 1):
            K = s.K(h)
            Kf = 0.5 * (K[:-1] + K[1:])
            g = (h[1:] - h[:-1]) / dz + 1.0
            R = (so.C(h) * (h - hn) + so.Ss * (h - hn)) * dz / dt
            R[:-1] += -Kf * g
            R[1:] -= -Kf * g
            gtop = -h[-1] / (dz / 2) + 1.0
            R[-1] -= K[-1] * gtop if face else s.r
            main = so.C(h) * dz / dt + so.Ss * dz / dt
            up = np.zeros(N)
            lo = np.zeros(N)
            a = Kf / dz
            main[:-1] += a
            up[1:] -= a
            main[1:] += a
            lo[:-1] -= a
            if face:
                main[-1] += K[-1] / (dz / 2)
            dK = so.ks * s.km * so.dkr(h)
            main[:-1] += -0.5 * dK[:-1] * g
            up[1:] += -0.5 * dK[1:] * g
            main[1:] += 0.5 * dK[1:] * g
            lo[:-1] += 0.5 * dK[:-1] * g
            if face:
                main[-1] -= dK[-1] * gtop
            dh = solve_banded((1, 1), np.vstack([up, main, lo]), -R)
            if not np.all(np.isfinite(dh)):
                return None, it
            h = h + dh
            if np.abs(dh).max() * col1d.GW < tol:
                return h, it
        return None, maxit

    def salva(s):
        return (s.h.copy(), s.face, s.t, s.r, dict(s.stats))

    def carrega(s, e):
        s.h, s.face, s.t, s.r, s.stats = e[0].copy(), e[1], e[2], e[3], dict(e[4])


# ---- opcoes
def le_opcoes(args, sub):
    o, pos = dict(PADRAO), []
    dadas = set()
    for a in args:
        if "=" not in a:
            pos.append(a)
            continue
        k, v = a.split("=", 1)
        if k not in PADRAO:
            erro(f"opcao desconhecida: {a}")
        o[k] = v
        dadas.add(k)
    if sub == "tabela" and ("eps" in dadas or "fecho" in dadas):
        if "vg_tab" in dadas or "vg_tab_psi" in dadas:
            erro("tabela: eps=/fecho= ou vg_tab=/vg_tab_psi=, nao os dois")
        o["vg_tab"] = (o["eps"] or "1e-3") + "," + (o["fecho"] or "reta")
        dadas.add("vg_tab")
    elif o["eps"] or o["fecho"]:
        erro("eps= e fecho= so no tabela (aqui: vg_tab=eps,fecho)")
    leis = [k for k in ("vg_tab_psi", "ippisch", "col1d_tab") if k in dadas] + (["vg_tab"] if "vg_tab" in dadas else [])
    if len(leis) > 1:
        erro(f"lei de k_r: uma so entre vg_tab, vg_tab_psi, ippisch e col1d_tab (dadas: {', '.join(leis)})")
    if leis and leis[0] != "vg_tab":
        o["vg_tab"] = "0"
    f = {}
    for k in ("n", "primeira", "semente", "perfis", "grava_xi", "localiza", "N", "M"):
        f[k] = int(o[k])
    for k in ("kmult", "H", "ks", "Kw", "gw", "I", "chuva", "succao0", "succao_max", "dt", "beta", "c", "phi", "gamma_d",
              "gamma", "fs_zmin", "covk"):
        f[k] = float(o[k])
    f["dt_min"] = float(o["dt_min"][:-1]) / 3600. if o["dt_min"].endswith("s") else float(o["dt_min"])
    f["tempos"] = sorted({float(t) for t in o["tempos"].split(",") if t})
    vg = [float(v) for v in o["vg"].split(",")]
    if len(vg) != 4:
        erro("vg=theta_s,theta_r,alfa,n")
    f["vg"] = vg
    nw = [float(v) for v in o["newton"].split(",")]
    f["tol"], f["maxit"] = nw[0], int(nw[1]) if len(nw) > 1 else 40
    for k in ("saida", "gerador", "desloc_xi", "nos", "vg_tab", "vg_tab_psi", "ippisch", "col1d_tab", "nivel0", "forma",
              "fs", "fs_peso", "fs_em", "nystrom", "lv", "grava"):
        f[k] = o[k]
    if f["gerador"] not in ("pz", "numpy"):
        erro("gerador=pz|numpy")
    if f["forma"] not in ("mista", "h"):
        erro("forma=mista|h")
    if f["fs"] not in ("talude", "nenhum") or f["fs_peso"] not in ("integrado", "constante") or \
            f["fs_em"] not in ("nos", "centros"):
        erro("fs=talude|nenhum, fs_peso=integrado|constante, fs_em=nos|centros")
    if f["fs_em"] == "centros" and (f["fs_peso"] != "integrado" or f["fs_zmin"] > 0. or f["nos"]):
        erro("fs_em=centros (o fs_perfil do col1d): so com fs_peso=integrado, fs_zmin=0 e sem nos=")
    if f["dt"] <= 0. or not f["tempos"] or f["N"] < 2:
        erro("dt > 0, tempos= e N >= 2")
    f["dadas"] = dadas
    return f, pos


def faz_solo(o):
    """o solo da lei de k_r escolhida e o rotulo lei_k (sem virgulas: vai no CSV)"""
    ts, tr, a, n = o["vg"]
    kw = dict(ts=ts, tr=tr, a_kpa=a, n=n, ks=o["ks"], Ss=ts * o["gw"] / o["Kw"])
    if o["ippisch"]:
        pe = float(o["ippisch"])
        return Solo(psi_e_kpa=pe, **kw), f"ippisch_{pe:g}"
    if o["col1d_tab"]:
        e = float(o["col1d_tab"])
        return SoloTab(eps=e, **kw), f"col1dtab_{e:g}"
    if o["vg_tab_psi"]:
        v = o["vg_tab_psi"].split(",")
        psic, fecho = float(v[0]), (v[1] if len(v) > 1 else "reta")
        e = eps_de_psi(a, n, psic)
        return SoloTabVG(eps=e, fecho=fecho, **kw), f"tabpsi_{psic:g}_{fecho}"
    v = o["vg_tab"].split(",")
    e = float(v[0])
    if e == 0.:
        return Solo(psi_e_kpa=0., **kw), "vg"
    fecho = v[1] if len(v) > 1 else "reta"
    return SoloTabVG(eps=e, fecho=fecho, **kw), f"tab_{e:g}_{fecho}"


# ---- a base de Nystrom 1D (nystrom_base.py coluna=H)
class Base:
    def __init__(s, arq, M):
        with open(arq) as a:
            t = a.readline().split()
            n, mArq = int(t[0]), int(t[1])
            cab = dict(x.split("=", 1) for x in t[2:] if "=" in x)
            nos = np.array([[float(v) for v in a.readline().split()] for _ in range(n)])
            lam = np.array([float(v) for v in a.readline().split()])
            f = np.loadtxt(a, max_rows=n, ndmin=2)
        if cab.get("nucleo", "exp") not in ("exp", "gauss") or "lh" not in cab or "lv" not in cab:
            erro(f"{arq}: a base precisa do cabecalho nucleo=, lh=, lv= (nystrom_base.py coluna=H)")
        s.lh, s.lv, s.nucleo = float(cab["lh"]), float(cab["lv"]), cab.get("nucleo", "exp")
        if not s.lh > 1e20:
            erro(f"{arq}: lh = {s.lh:g}: nao e uma base 1D (nystrom_base.py coluna=H grava lh=1e+30)")
        if np.abs(nos[:, 0]).max() > 0. or len(lam) != mArq or f.shape != (n, mArq):
            erro(f"{arq}: base 1D ilegivel (x != 0 ou tamanhos errados)")
        s.M = min(M, mArq)
        s.y, s.w, s.lam, s.f = nos[:, 1], nos[:, 2], lam[:s.M], f[:, :s.M]
        s.variancia = s.lam.sum() / s.w.sum()
        s.arq = arq

    def C(s, d):
        return np.exp(-np.abs(d) / s.lv) if s.nucleo == "exp" else np.exp(-(d / s.lv) ** 2)

    def modos(s, y):
        """sqrt(lambda_k) f_k(y) pela formula de Nystrom (o RandomField::ModeValues): (len(y) x M)"""
        c = s.C(np.asarray(y)[:, None] - s.y[None, :]) * s.w[None, :]
        return (c @ s.f) / np.sqrt(s.lam)[None, :]


# ---- o FS de Lu-Godt nos pontos do perfil (spec B.9)
def perfil(col, zs=None, t0=False):
    """os pontos do perfil de cima para baixo: z (m) e a carga h (m). A superficie: 0 na face, a extrapolacao do fluxo
    (h_N + dz/2 (r/K_N - 1), a do criterio de troca do col1d) fora dela; a base: h_1 + dz/2 (fluxo nulo); entre os
    centros, linear; no estado inicial (t0), a formula do estado inicial em cada ponto. Sem zs: z = j dz/2, j = 0..2N"""
    N, dz, H = col.N, col.dz, col.H
    z = np.arange(2 * N + 1) * (0.5 * dz) if zs is None else zs
    if t0:  # o estado inicial nos proprios pontos (a sua formula: succao uniforme ou nivel)
        return z, col.h0f(H - z)
    hs = 0. if col.face else col.h[-1] + 0.5 * dz * (col.r / col.K(col.h)[-1] - 1.)
    zk = np.concatenate([[0.], H - col.y[::-1], [H]])
    hk = np.concatenate([[hs], col.h[::-1], [col.h[0] + 0.5 * dz]])
    return z, np.interp(z, zk, hk)


def fs_pontos(z, h, so, o):
    """FS = [c' + (W cos^2 b - sigma_s) tan phi'] / (W sin b cos b) em z > max(0, fs_zmin); W = trapezios de
    gamma_d + theta gw desde z = 0 (ou gamma z com fs_peso=constante); sigma_s = S_e p_w se p_w < 0, senao p_w"""
    gw = col1d.GW
    z, h = np.asarray(z, dtype=float), np.asarray(h, dtype=float)
    if o["fs_peso"] == "integrado":
        gt = o["gamma_d"] + so.theta(h) * gw
        W = np.concatenate([[0.], np.cumsum(0.5 * (gt[1:] + gt[:-1]) * np.diff(z))])
        if z[0] > 0.:
            erro("fs: o perfil comeca em z = 0 (o W integra desde a superficie)")
    else:
        W = o["gamma"] * z
    pw = gw * h
    sig = np.where(pw < 0., so.Se(h) * pw, pw)
    b, tp = math.radians(o["beta"]), math.tan(math.radians(o["phi"]))
    ok = z > max(0., o["fs_zmin"]) + 1e-12
    with np.errstate(divide="ignore", invalid="ignore"):
        FS = np.where(ok, (o["c"] + (W * math.cos(b) ** 2 - sig) * tp) / (W * math.sin(b) * math.cos(b)), np.inf)
    return FS, W


def avalia(col, so, o, zs, pw0=None, t0=False):
    """FS_min, z_crit, base, frente, p_w na base e no topo, a menor carga e o perfil (z, h, FS, W)"""
    if o["fs_em"] == "centros":
        z, FS, h = fs_perfil(col, c=o["c"], phi=o["phi"], beta=o["beta"], gd=o["gamma_d"])
        W = np.full(len(z), np.nan)
        zb, hb = perfil(col, None, t0)
        base_lim = col.H - col.dz
    else:
        z, h = perfil(col, zs, t0)
        FS, W = fs_pontos(z, h, so, o) if o["fs"] == "talude" else (np.full(len(z), np.inf), np.zeros(len(z)))
        zb, hb = z, h
        base_lim = col.H - 1e-9
    k = int(np.argmin(FS))
    pw = col1d.GW * hb
    d = {"FS": float(FS[k]) if o["fs"] == "talude" else float("nan"), "z_crit": float(z[k]),
         "base": int(z[k] >= base_lim), "pw_base": float(pw[-1]), "pw_topo": float(pw[0]), "h_min": float(hb.min()),
         "face": int(col.face)}
    if pw0 is not None:
        j = np.flatnonzero(pw - pw0 < 1.)
        d["frente"] = float(zb[j[0]]) if len(j) else float(zb[-1])
    d["perfil"] = (z, h, FS, W)
    return d


# ---- uma realizacao
def grade(o):
    """os fins de passo (h): k dt, os tempos e o fim da chuva"""
    T = max(o["tempos"])
    ts = {round(k * o["dt"], 9) for k in range(1, int(math.ceil(T / o["dt"] - 1e-9)) + 1)}
    ts |= {round(t, 9) for t in o["tempos"] if t > 0.}
    if 0. < o["chuva"] < T:
        ts.add(round(o["chuva"], 9))
    return sorted(t for t in ts if t <= T + 1e-9)


def nova_coluna(o, so, kmult):
    gw = col1d.GW
    col = ColunaRef(so, H=o["H"], N=o["N"], kmult=kmult, r=o["I"] / 1000. / 3600., h0=-o["succao0"] / gw)
    col.forma = o["forma"]
    if o["nivel0"]:
        zw, hmin = float(o["nivel0"]), -o["succao_max"] / gw
        col.h0f = lambda y: np.maximum(zw - np.asarray(y, dtype=float), hmin)
    else:
        h0 = -o["succao0"] / gw
        col.h0f = lambda y: np.full(np.shape(y), h0)
    col.h = col.h0f(col.y)
    return col


def roda(o, so, kmult, zs=None, perfis=None):
    """o transiente de uma realizacao: o FS depois de cada passo; devolve as linhas dos tempos (dicts) e o resumo"""
    t0 = time.perf_counter()
    col = nova_coluna(o, so, kmult)
    q = col.r
    _, h0 = perfil(col, zs, True) if o["fs_em"] == "nos" else perfil(col, None, True)
    pw0 = col1d.GW * h0
    kw = dict(metodo="newton", dtmin=o["dt_min"] * 3600., tol=o["tol"], maxit=o["maxit"])
    linhas, tpass, fsmin_hist, tparou = [], -1., math.inf, -1.
    a = avalia(col, so, o, zs, pw0, t0=True)
    fs_ant = a["FS"]
    if o["fs"] == "talude" and fs_ant < 1.:
        tpass = 0.
    fsmin_hist = min(fsmin_hist, fs_ant) if o["fs"] == "talude" else math.nan

    def linha(t, a, situacao="ok"):
        d = {"t_h": t, "situacao": situacao, "FS_passagem": fsmin_hist, "t_passagem": tpass if 0. <= tpass <= t + 1e-9
             else -1., "passos": col.stats["passos"], "newton": col.stats["it"], "s": time.perf_counter() - t0,
             "t_parou": tparou}
        d.update({k: v for k, v in a.items() if k != "perfil"})
        if perfis is not None and "perfil" in a:
            perfis.append((t, a["perfil"]))
        return d

    if 0. in o["tempos"]:
        linhas.append(linha(0., a))
    tcur = 0.
    for tn in grade(o):
        col.r = q if tcur < o["chuva"] - 1e-9 else 0.
        e0 = col.salva()
        try:
            col.advance_to(tn * 3600., dt=(tn - tcur) * 3600., **kw)
        except PARADAS:
            col.carrega(e0)
            tparou = tcur
            break
        a = avalia(col, so, o, zs, pw0)
        if o["fs"] == "talude":
            if tpass < 0. and a["FS"] < 1. <= fs_ant:
                tL, fL, tR, fR = tcur, fs_ant, tn, a["FS"]
                if o["localiza"]:
                    fim = col.salva()
                    eL = e0
                    for _ in range(7):
                        tm = 0.5 * (tL + tR)
                        col.carrega(eL)
                        try:
                            col.advance_to(tm * 3600., dt=(tm - tL) * 3600., **kw)
                        except PARADAS:
                            break
                        fm = avalia(col, so, o, zs, pw0)["FS"]
                        if fm < 1.:
                            tR, fR = tm, fm
                        else:
                            tL, fL, eL = tm, fm, col.salva()
                    col.carrega(fim)
                tpass = tL + (fL - 1.) / (fL - fR) * (tR - tL) if fL > fR else tR
            fsmin_hist = min(fsmin_hist, a["FS"])
            fs_ant = a["FS"]
        tcur = tn
        if any(abs(tn - t) < 1e-9 for t in o["tempos"]):
            linhas.append(linha(tn, a))
    if tparou >= 0.:
        for t in o["tempos"]:
            if t > tparou + 1e-9:
                d = linha(t, {"FS": math.nan, "z_crit": math.nan, "base": 0, "pw_base": math.nan, "pw_topo": math.nan,
                              "h_min": math.nan, "face": 0, "frente": math.nan}, "transiente")
                linhas.append(d)
    return linhas, {"s": time.perf_counter() - t0, "passos": col.stats["passos"], "newton": col.stats["it"],
                    "subdivisoes": col.stats["fail"], "t_parou": tparou, "t_passagem": tpass}


def fmt(v):
    if isinstance(v, str):
        return v
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return "nan"
    return f"{v:.10g}"


def grava_csv(arq, real, linhas, extra):
    """as linhas de uma realizacao nas colunas do SlopeRainColumn"""
    out = []
    for d in linhas:
        v = {"real": real, "t_h": d["t_h"], "FS_min": d["FS"], "situacao": d["situacao"], "z_crit": d["z_crit"],
             "base": d["base"], "FS_passagem": d["FS_passagem"], "t_passagem": d["t_passagem"],
             "frente_m": d.get("frente", math.nan), "pw_base": d["pw_base"], "pw_topo": d["pw_topo"],
             "nos_face": d["face"], "passos": d["passos"], "newton": d["newton"], "s_transiente": d["s"],
             "t_parou": d["t_parou"]}
        v.update(extra)
        out.append(",".join(fmt(v[c]) for c in COLUNAS.split(",")))
    arq.write("\n".join(out) + "\n")


def grava_perfil(nome, perfis):
    with open(nome, "w") as a:
        a.write("t_h,z,pw,FS,W\n")
        for t, (z, h, FS, W) in perfis:
            for i in range(len(z)):
                a.write(f"{t:g},{z[i]:.15g},{col1d.GW * h[i]:.15g},{FS[i]:.15g},{W[i]:.15g}\n")


def le_nos(arq):
    """as profundidades z (m) da primeira coluna de um arquivo (linhas que nao sao numeros sao puladas)"""
    z = []
    for linha in open(arq):
        t = linha.replace(",", " ").split()
        if t:
            try:
                z.append(float(t[0]))
            except ValueError:
                pass
    z = np.unique(np.array(z))
    if not len(z) or z[0] < 0.:
        erro(f"nos={arq}: profundidades z >= 0 (m)")
    return z if z[0] == 0. else np.concatenate([[0.], z])


def prepara(o):
    """o solo, a gravidade e os nos do perfil"""
    col1d.GW = o["gw"]
    so, lei = faz_solo(o)
    zs = le_nos(o["nos"]) if o["nos"] else None
    if zs is not None and zs[-1] > o["H"] + 1e-9:
        erro(f"nos=: profundidade {zs[-1]} > H = {o['H']}")
    return so, lei, zs


def grava_parametros(o, lei, base, so):
    with open(o["saida"] + "_parametros.txt", "w") as a:
        a.write(f"# coluna_ref.py {' '.join(sys.argv[1:])}\n")
        for k in sorted(o):
            if k != "dadas":
                a.write(f"{k}={o[k]}\n")
        a.write(f"lei_k={lei}\n")
        if isinstance(so, SoloTabVG):
            a.write(f"tabela: psi_c={so.info['xc'] / so.a * col1d.GW:.6g} kPa k_c={so.info['kc']:.6g} "
                    f"pontos={len(so.St)}\n")
        if base is not None:
            a.write(f"base={base.arq} nucleo={base.nucleo} lv={base.lv:g} M={base.M} KL_variancia={base.variancia:.6f}\n")


def det(o):
    so, lei, zs = prepara(o)
    perf = []
    linhas, r = roda(o, so, np.full(o["N"], o["kmult"]), zs, perf)
    extra = {"k_mult_med": o["kmult"], "k_mult_min": o["kmult"], "k_mult_max": o["kmult"], "log_peso": 0.,
             "ny": o["N"], "dt": o["dt"], "elementos": "vf", "massa": o["forma"], "lei_k": lei, "cov_k": 0., "lv": 0.,
             "M": 0, "KL_variancia": 0.}
    with open(o["saida"] + "_0_-1.csv", "w") as a:
        a.write(COLUNAS + "\n")
        grava_csv(a, -1, linhas, extra)
    grava_perfil(o["saida"] + "_perfil_-1.csv", perf)
    grava_parametros(o, lei, None, so)
    print(f"coluna de {o['H']:g} m, {o['N']} celulas (dz = {o['H'] / o['N']:.4g} m), dt = {o['dt']:g} h, k_r {lei}, "
          f"forma {o['forma']}, FS nos {'centros' if o['fs_em'] == 'centros' else 'pontos ' + ('de ' + o['nos'] if zs is not None else 'z = j dz/2')}"
          f", k_s = {o['kmult'] * o['ks']:.4g} m/s, I = {o['I']:g} mm/h")
    print("     t (h)   FS_min  z_crit (m) base  frente (m)  pw_base  pw_topo (kPa)  h_min (m)  face  passos  newton")
    for d in linhas:
        print(f"   {d['t_h']:7.3f}  {d['FS']:7.4f}  {d['z_crit']:9.4f}  {d['base']:3d}  {d.get('frente', math.nan):9.4f}"
              f"  {d['pw_base']:8.3f} {d['pw_topo']:8.3f}       {d['h_min']:8.4f}  {d['face']:3d}  {d['passos']:6d}"
              f"  {d['newton']:6d}  {d['situacao']}")
    print(f"  t_passagem = {r['t_passagem']:.4f} h; {r['s']:.2f} s, {r['passos']} passos, {r['newton']} iteracoes de "
          f"Newton, {r['subdivisoes']} subdivisoes" + (f"; o transiente parou em {r['t_parou']:.4f} h" if r["t_parou"] >= 0
                                                       else ""))
    print(f"  arquivos: {o['saida']}_0_-1.csv {o['saida']}_perfil_-1.csv {o['saida']}_parametros.txt")
    return 2 if r["t_parou"] >= 0. else 0


class Campo:
    """o campo de k_mult nos centros das celulas de uma base 1D"""

    def __init__(s, o, base, y):
        s.base, s.B = base, base.modos(y)
        s.sl = math.sqrt(math.log(1. + o["covk"] ** 2))
        s.mu = -0.5 * s.sl ** 2

    def kmult(s, xi):
        return np.exp(s.mu + s.sl * (s.B @ xi))


def mc_bloco(o, a, n, xis=None, so_lei_zs=None, base=None, quieto=False):
    """as realizacoes a..a+n-1 (ou as dos xi dados: lista de (real, xi)) em <saida>_<a>_<b>.csv (gravado por inteiro
    no fim, por um .tmp); devolve o resumo"""
    so, lei, zs = so_lei_zs if so_lei_zs else prepara(o)
    if base is None:
        if not o["nystrom"]:
            erro("mc/xi: nystrom=<base.txt> (nystrom_base.py coluna=H)")
        base = Base(o["nystrom"], o["M"])
    if o["lv"] and abs(float(o["lv"]) - base.lv) > 1e-12 * base.lv:
        erro(f"lv={o['lv']} mas a base {base.arq} tem lv = {base.lv:g}")
    y = (np.arange(o["N"]) + 0.5) * (o["H"] / o["N"])
    campo = Campo(o, base, y)
    desl = None
    if o["desloc_xi"]:
        desl = np.zeros(base.M)
        v = np.loadtxt(o["desloc_xi"], ndmin=1)
        if len(v) > base.M:
            erro(f"desloc_xi: {len(v)} numeros > M = {base.M}")
        desl[:len(v)] = v
    reais = [(r, None) for r in range(a, a + n)] if xis is None else xis
    nome = f"{o['saida']}_{a}_{a + n - 1}.csv"
    tmp = nome + ".tmp"
    if a == 0 or xis is not None:
        grava_parametros(o, lei, base, so)
    soma_s, falhas = 0., {}
    with open(tmp, "w") as arq:
        arq.write(COLUNAS + "\n")
        for r, xi in reais:
            if xi is None:
                xi = xi_realizacao(o, r, base.M)
                if desl is not None:
                    xi = xi + desl
            logp = float(-xi @ desl + 0.5 * desl @ desl) if desl is not None else 0.
            if o["grava_xi"]:
                np.savetxt(f"{o['saida']}_xi_{r}.txt", xi, fmt="%.17g")
            km = campo.kmult(xi)
            perf = [] if r < o["perfis"] else None
            linhas, res = roda(o, so, km, zs, perf)
            if perf is not None:
                grava_perfil(f"{o['saida']}_perfil_{r}.csv", perf)
            extra = {"k_mult_med": float(km.mean()), "k_mult_min": float(km.min()), "k_mult_max": float(km.max()),
                     "log_peso": logp, "ny": o["N"], "dt": o["dt"], "elementos": "vf", "massa": o["forma"],
                     "lei_k": lei, "cov_k": o["covk"], "lv": base.lv, "M": base.M, "KL_variancia": base.variancia}
            grava_csv(arq, r, linhas, extra)
            soma_s += res["s"]
            for d in linhas:
                if d["situacao"] == "ok" and d["FS"] < 1.:
                    falhas[d["t_h"]] = falhas.get(d["t_h"], 0) + 1
            if not quieto:
                print(f"  realizacao {r}: {res['s']:.2f} s, {res['passos']} passos, {res['newton']} iteracoes de Newton"
                      + (f", t_passagem {res['t_passagem']:.3f} h" if res["t_passagem"] >= 0. else "")
                      + (f", o transiente parou em {res['t_parou']:.3f} h" if res["t_parou"] >= 0. else ""), flush=True)
    os.replace(tmp, nome)
    return {"arquivo": nome, "n": len(reais), "s": soma_s, "falhas": falhas}


def le_xi(arqs, M):
    """os xi dos arquivos do driver (um numero por linha, os M de k_s); a realizacao vem do nome (_xi_<n>) ou da ordem"""
    out = []
    for i, f in enumerate(arqs):
        v = np.loadtxt(f, ndmin=1)
        if len(v) > M:
            erro(f"xi: {f} tem {len(v)} numeros > M = {M}")
        xi = np.zeros(M)
        xi[:len(v)] = v
        m = re.search(r"_xi_(-?\d+)\.txt$", f)
        out.append((int(m.group(1)) if m else i, xi))
    return out


def fs_csv(o, arq):
    """o FS de um perfil gravado (V1.3)"""
    col1d.GW = o["gw"]
    so, _ = faz_solo(o)
    with open(arq) as a:
        rd = csv.DictReader(a)
        cols = rd.fieldnames
        dados = list(rd)
    cz = next((c for c in ("z", "z_m") if c in cols), None)
    cy = next((c for c in ("y", "y_m") if c in cols), None)
    cp = next((c for c in ("pw", "p_w", "pw_kPa") if c in cols), None)
    if cp is None or (cz is None and cy is None):
        erro(f"fs: {arq} precisa das colunas z (ou y) e pw (kPa); tem {cols}")
    grupos = {}
    for d in dados:
        grupos.setdefault(float(d["t_h"]) if "t_h" in cols else 0., []).append(d)
    print(f"{arq}: {len(dados)} linhas, {len(grupos)} instante(s)")
    pior = 0.
    for t, g in sorted(grupos.items()):
        z = np.array([float(d[cz]) for d in g]) if cz else o["H"] - np.array([float(d[cy]) for d in g])
        pw = np.array([float(d[cp]) for d in g])
        ordem = np.argsort(z, kind="stable")
        z, pw = z[ordem], pw[ordem]
        FS, _ = fs_pontos(z, pw / col1d.GW, so, o)
        k = int(np.argmin(FS))
        txt = f"  t = {t:g} h: FS_min = {FS[k]:.10g} em z = {z[k]:.6g} m (base {int(z[k] >= o['H'] - 1e-9)})"
        if "FS" in cols:
            ref = np.array([float(d["FS"]) for d in g])[ordem]
            fin = np.isfinite(ref) & np.isfinite(FS)
            dif = float(np.abs(FS[fin] - ref[fin]).max()) if fin.any() else 0.
            pior = max(pior, dif)
            txt += f"; max |FS - FS do arquivo| = {dif:.2e}"
        print(txt)
    if pior:
        print(f"  maior diferenca: {pior:.2e}")


def tabela(o):
    ts, tr, a, n = o["vg"]
    swr = tr / ts
    if o["vg_tab_psi"]:
        v = o["vg_tab_psi"].split(",")
        eps, fecho = eps_de_psi(a, n, float(v[0])), (v[1] if len(v) > 1 else "reta")
    else:
        v = o["vg_tab"].split(",")
        eps, fecho = float(v[0]), (v[1] if len(v) > 1 else "reta")
    S, L, info = tabela_vg(n, swr, eps, fecho)
    m = 1. - 1. / n
    psic = info["xc"] / a
    # erro de k_r em [psi_c, 1e3 kPa] (2000 pontos log)
    psis = np.exp(np.linspace(math.log(psic), math.log(1e3), 2000))
    Sw = swr + (1. - swr) * se_vg(a * psis, n)
    kt = kr_tabela(Sw, S, L)[0]
    km = np.exp(ln_kr_mualem(a * psis, n))
    if fecho == "ippisch":
        km = km / info["kc"]
    erro_max = float(np.abs(kt / km - 1.).max())
    print(f"tabela vg_tab: n = {n:g}, S_wr = {swr:.6g}, eps = {eps:.6g}, fecho {fecho}: psi_c = {psic:.4g} kPa, "
          f"k_c = {info['kc']:.4f}, {len(S)} pontos ({info['tabela']} da Mualem), erro max de k_r em [psi_c, 1e3 kPa] "
          f"(2000 pontos) = {100 * erro_max:.4f}%")
    ps = np.array([0.01, 0.1, 1., 10.])
    kp = kr_tabela(swr + (1. - swr) * se_vg(a * ps, n), S, L)[0]
    kmp = np.exp(ln_kr_mualem(a * ps, n))
    print("  k_r(0.01/0.1/1/10 kPa) = " + "/".join(f"{x:.4f}" for x in kp) + ", Mualem " +
          "/".join(f"{x:.4f}" for x in kmp))
    print("  k_r(0.01/0.1/1/10 kPa) com 15 digitos: " + " ".join(f"{x:.15g}" for x in kp))
    # |dk_r/dpsi| = k_r g |dS_w/dpsi| (1/kPa), numa grade fina de 1e-7 a 1e3 kPa, e o salto em psi_c
    pg = np.exp(np.linspace(math.log(1e-7), math.log(1e3), 400001))

    def dk(p):
        y = (a * p) ** n
        dS = (1. - swr) * m * n * a ** n * p ** (n - 1.) * (1. + y) ** (-m - 1.)
        return kr_tabela(swr + (1. - swr) * (1. + y) ** (-m), S, L)[1] * dS
    d = dk(pg)
    i = int(np.argmax(d))
    print(f"  max |dk_r/dpsi| = {d[i]:.3g} /kPa em psi = {pg[i]:.3g} kPa; em psi_c: {float(dk(psic * (1 + 1e-9))):.3g} "
          f"(tabela) -> {float(dk(psic * (1 - 1e-9))):.3g} /kPa (fecho)" + (f"; g_c = {info['gc']:.6g}" if "gc" in info
                                                                          else ""))
    if o["grava"]:
        with open(o["grava"], "w") as f:
            f.write("i,S_w,lnk_r,k_r\n")
            for j in range(len(S)):
                f.write(f"{j},{S[j]:.17g},{L[j]:.17g},{math.exp(L[j]):.17g}\n")
        print(f"  gravada: {o['grava']}")


def main():
    if len(sys.argv) < 2 or sys.argv[1] not in ("det", "mc", "xi", "fs", "tabela"):
        sys.exit(__doc__)
    sub = sys.argv[1]
    o, pos = le_opcoes(sys.argv[2:], sub)
    if sub == "det":
        if pos:
            erro(f"det: argumentos sem '=': {pos}")
        sys.exit(det(o))
    if sub == "tabela":
        tabela(o)
        return
    if sub == "fs":
        if len(pos) != 1:
            erro("fs <perfil.csv> [opcoes]")
        fs_csv(o, pos[0])
        return
    if sub == "mc":
        if pos:
            erro(f"mc: argumentos sem '=': {pos}")
        r = mc_bloco(o, o["primeira"], o["n"])
    else:
        if not pos:
            erro("xi nystrom=<base.txt> <arquivos de xi>")
        base = Base(o["nystrom"], o["M"]) if o["nystrom"] else erro("xi: nystrom=<base.txt>")
        xis = le_xi(pos, base.M)
        r = mc_bloco(o, xis[0][0], len(xis), xis=xis, base=base)
    print(f"  {r['n']} realizacoes, {r['s'] / max(r['n'], 1):.2f} s por realizacao; FS < 1 por instante: " +
          ", ".join(f"{t:g} h {k}" for t, k in sorted(r["falhas"].items())))
    print(f"  arquivos: {r['arquivo']} {o['saida']}_parametros.txt")


if __name__ == "__main__":
    main()
