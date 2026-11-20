#!/usr/bin/env python3
"""la.py: a analise limite cinematica do Ceron (2025) no talude seco do Cho (2010, Exemplo 2), realizacao a realizacao.

  * talude: H = 10 m, beta = 45 graus, crista de 10 m e terreno de 10 m a frente do pe, base 5 m abaixo do pe (o dominio
    [0, 30] x [0, 15] do Cho2010). Coordenadas do Ceron: origem na aresta da crista, x para fora, y para baixo:
    kin.Talude(10, 45 graus, x_esq = -10, x_dir = 20, y_base = 15). Solo seco: u = 0, gamma = 20 kN/m3.
  * campos: campos_cho.CamposCho (c' em kPa e phi' em rad na grade de passo h, os mesmos do MEF dados os xi).
  * Gamma: rapido.minimizar_rapido (espiral log por trechos de H/100, phi_i = limite superior de max phi ao longo do
    proprio arco, dissipacao 'exata', mecanismos I (face/pe) e II (abaixo do pe)), com partida quente nos mecanismos
    criticos do deterministico (cinematico.minimizar, o metodo exato; guardado em partidas.json).
  * FS (Eqs. 57-58 do artigo, generalizadas ao campo): o F tal que Gamma(c/F, atan(tan phi/F)) = 1, com c e phi
    reduzidos no por no. Secante em ln Gamma x ln F (Gamma ~ F^(-1/0,32) perto do deterministico), cada passo uma
    minimizacao com partida quente nos mecanismos do passo anterior. Gamma < 1 <=> FS < 1, entao a P_f e a mesma; o FS
    so serve para comparar valores com o FS do SRM do MEF.
  * modo de ruptura (como na Fig. 20 do Ceron): 'pe' (classe I com B no pe), 'face' (classe I acima do pe), 'base'
    (classe II, B a frente do pe)."""
from __future__ import annotations

import json
import math
import os
import sys
import time

import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
import campos_cho as cc  # noqa: E402
import cinematico as kin  # noqa: E402
import rapido as rp  # noqa: E402

TALUDE = kin.Talude(10.0, math.radians(45.0), -10.0, 20.0, 15.0)
GAMMA = 20.0
EXPOENTE_FS = 0.32          # ln FS ~ 0,32 ln Gamma (deterministico: ln 1,203 / ln 1,777)
# A busca rapida padrao do rapido.py (grade 5 x 5 x 4, 1 Nelder-Mead por classe) perde as vezes um mecanismo pela face
# quando o melhor e pelo pe (3 de 30 realizacoes do Cho com +2,5 a +6,6 % em Gamma, todas com Gamma > 1,8). A busca
# forte (grade 7 x 7 x 5 e ate 3 Nelder-Mead por classe, em bacias diferentes) acerta todas, com 2,7 vezes o custo.
# Politica: a rapida em todas e a forte de novo quando Gamma < G_ROBUSTO (a regiao que decide a falha) e em toda
# minimizacao do FS (Gamma perto de 1).
BUSCA_RAPIDA = {}
BUSCA_FORTE = dict(grade=(7, 7, 5), n_nm=3, folga_2a=0.2)
G_ROBUSTO = 1.5
# Passo da grade dos campos. O phi_i (max phi ao longo do arco) le o campo P1 da grade, e o campo de 2000 modos com o
# nucleo exponencial e aspero na escala de 0,1 m: com h = 0,1 o max e subestimado e Gamma (FS) saem ~2,5 % (~1 %)
# abaixo do convergido. Estudo de h (2 realizacoes, Gamma 0,95 e 1,22): h 0,1 -> 0,05 -> 0,025 -> 0,0125 muda o FS
# +0,5 %, +0,2 %, ver README; o custo por realizacao quase nao muda (os arcos tem amostras fixas). Padrao: 0,025.
H_PADRAO = 0.025


def grade(h, c, phi):
    return rp.CampoGrade(-cc.X_CRISTA, h, c, phi, np.zeros_like(c), GAMMA, TALUDE)


class _SemMecanismo:
    """nenhum mecanismo admissivel achado (campo muito resistente: com phi ate ~70 graus as espirais da grade da busca
    nao chegam a superficie, "A fora da crista" ou "nao entra no solo"; 1 em ~150 realizacoes do r07 com h = 0,025, em
    que o h = 0,1 dava Gamma ~ 60). Gamma = inf: nao e falha."""
    classe = "-"
    sB = xA = xB = float("nan")


def nenhum_mecanismo():
    return kin.Resultado(math.inf, _SemMecanismo(), {}, 0, 0.0, 0, 0)


def modo(m):
    if m.classe not in ("I", "II"):
        return "nenhum"
    if m.classe == "II":
        return "base"
    return "pe" if abs(m.sB) < 1e-9 else "face"


def partidas_det(h=0.1, arquivo=None):
    """os mecanismos criticos (por classe) do deterministico (c = 10, phi = 30 uniformes) pelo metodo exato. O campo
    uniforme nao depende da grade: o arquivo e sempre o da grade h = 0,1 (partidas_h0.1.json), qualquer que seja o h
    das realizacoes (so serve de partida para a busca)."""
    h = 0.1
    arquivo = arquivo or os.path.join(AQUI, f"partidas_h{h:g}.json")
    if os.path.exists(arquivo):
        d = json.load(open(arquivo))
        return {k: tuple(v) for k, v in d["partidas"].items()}, d["Gamma"]
    nx, ny = int(round(cc.DX / h)) + 1, int(round(cc.DY / h)) + 1
    cg = grade(h, np.full((nx, ny), cc.MEDIA_C), np.full((nx, ny), math.radians(cc.MEDIA_PHI)))
    r = kin.minimizar(cg, TALUDE)
    part = {cl: (m.sB, m.t2, m.rB) for cl, m in r.por_classe.items()}
    tmp = arquivo + f".{os.uname().nodename}.{os.getpid()}"
    json.dump({"partidas": part, "Gamma": r.Gamma, "por_classe": {cl: m.Gamma for cl, m in r.por_classe.items()}},
              open(tmp, "w"), indent=1)
    os.replace(tmp, arquivo)
    return part, r.Gamma


class Avaliador:
    """Gamma (e FS) da analise limite dados os xi (2, M)."""

    def __init__(self, campos: cc.CamposCho, partidas=None, robusto=G_ROBUSTO):
        self.cc = campos
        self.h = campos.h
        self.part = partidas if partidas is not None else partidas_det(self.h)[0]
        self.robusto = robusto

    def minimizar(self, c, phi, partidas=None, forte=None):
        """forte=None: busca rapida (BUSCA_RAPIDA) e, se Gamma < robusto, de novo com BUSCA_FORTE partindo tambem dos
        mecanismos achados, e fica o menor; forte=False: so a rapida; forte=True: rapida e forte sempre"""
        cg = grade(self.h, c, phi)
        part = partidas or self.part

        def busca(p, opc):
            try:
                return rp.minimizar_rapido(cg, TALUDE, partidas=p, **opc)
            except RuntimeError:            # "nenhum mecanismo valido"
                return None

        r = busca(part, BUSCA_RAPIDA)
        if r is None:                       # a forte tem mais candidatos; sem nenhum: Gamma = inf (nao falha)
            r = busca({**self.part, **part}, BUSCA_FORTE)
            return r if r is not None else nenhum_mecanismo()
        if forte is False or (forte is None and r.Gamma >= self.robusto):
            return r
        part = {**part, **{cl: (m.sB, m.t2, m.rB) for cl, m in r.por_classe.items()}}
        r2 = busca(part, BUSCA_FORTE)
        if r2 is None:
            return r
        a, b = (r, r2) if r.Gamma <= r2.Gamma else (r2, r)
        a.n_aval += b.n_aval
        a.tempo_s += b.tempo_s
        return a

    def fs(self, c, phi, res, tol=2e-4, max_it=14):
        """F com Gamma(c/F, atan(tan phi/F)) = 1; res = o resultado com F = 1 (os seus mecanismos sao a partida).
        1) secante em ln Gamma x ln F com a busca rapida partindo dos mecanismos do passo anterior (o mecanismo critico
           muda pouco com F); 2) no F final, a busca forte confere: se achar ln Gamma < -5e-4 (outro mecanismo), a
           secante continua com a busca forte. Devolve (F, minimizacoes, convergiu); convergiu = False se max_it
           acabou antes do passo < tol (o F e o ultimo iterado)."""
        tphi = np.tan(phi)
        hist = [(0.0, math.log(res.Gamma))] if math.isfinite(res.Gamma) else []   # (ln F, ln Gamma)
        est = {"part": {cl: (m.sB, m.t2, m.rB) for cl, m in res.por_classe.items()}, "n": 0, "conv": True}

        def lgam(lf, forte):
            F = math.exp(lf)
            r = self.minimizar(c / F, np.arctan(tphi / F), est["part"], forte=forte)
            if r.por_classe:
                est["part"] = {cl: (m.sB, m.t2, m.rB) for cl, m in r.por_classe.items()}
            est["n"] += 1
            return math.log(r.Gamma)        # inf se nenhum mecanismo (F ainda pequeno)

        def secante(lf, forte):
            while est["n"] < max_it:
                lg = lgam(lf, forte)
                if not math.isfinite(lg):   # nenhum mecanismo: F maior
                    lf += 0.3
                    continue
                hist.append((lf, lg))
                if abs(lg) < 1e-7:
                    return lf
                if len(hist) < 2:           # so um ponto: a inclinacao tipica d ln Gamma / d ln F = -1/0,32
                    hist.insert(0, (lf - 0.1, lg + 0.1 / EXPOENTE_FS))
                (a0, g0), (a1, g1) = hist[-2], hist[-1]
                novo = a1 - g1 * (a1 - a0) / (g1 - g0) if abs(g1 - g0) > 1e-12 else a1
                novo = min(max(novo, a1 - 0.5), a1 + 0.5) if np.isfinite(novo) else a1
                if abs(novo - a1) < tol:
                    return novo
                lf = novo
            est["conv"] = False
            return lf

        lf = secante(EXPOENTE_FS * math.log(res.Gamma) if hist else math.log(3.0), False)
        lg = lgam(lf, True)
        if lg < -5e-4 and math.isfinite(lg):
            hist.append((lf, lg))
            est["conv"] = True
            lf = secante(lf + lg * EXPOENTE_FS, True)   # ln Gamma ~ ln Gamma0 - ln F / 0,32
        return math.exp(lf), est["n"], est["conv"]

    def avalia(self, xi, fs=False, fs_gamma=0.0):
        """Gamma, modo e tempos; FS se fs=True ou se Gamma < fs_gamma"""
        t0 = time.perf_counter()
        c, phi = self.cc.campos(xi)
        t1 = time.perf_counter()
        r = self.minimizar(c, phi)
        m = r.mecanismo
        out = dict(Gamma=r.Gamma, modo=modo(m), classe=m.classe, sB=m.sB, xA=m.xA, xB=m.xB, n_aval=r.n_aval,
                   descartes=r.n_normalidade + r.n_phi_nao_conv, t_campos=t1 - t0, t_gamma=r.tempo_s)
        if fs or r.Gamma < fs_gamma or not math.isfinite(r.Gamma):   # sem mecanismo: o FS decide (F > 1: seguro)
            t2 = time.perf_counter()
            F, nf, conv = self.fs(c, phi, r)
            out.update(FS=F, n_fs=nf, fs_conv=conv, t_fs=time.perf_counter() - t2)
        out["t"] = time.perf_counter() - t0
        return out


def le_kv(argv, **pad):
    a = dict(pad)
    for x in argv:
        if "=" not in x:
            raise SystemExit(f"argumento sem '=': {x}")
        k, v = x.split("=", 1)
        a[k] = v
    return a
