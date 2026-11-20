#!/usr/bin/env python3
"""mc.py: Monte Carlo do Ceron (2025, Secao 6, configuracao de referencia) com o avaliador rapido (rapido.py).

Cada realizacao i: xi ~ N(0, I) com o gerador default_rng(semente + i) (campos c, phi, k nessa ordem, M = 2000
termos cada), c e phi nos nos da grade dos campos, k por celula da malha do Darcy, Darcy permanente com o contorno
do Ceron, u interpolado na grade e Gamma pela busca rapida (partida quente no mecanismo critico deterministico).

    python3 -I mc.py n=200 primeira=0 semente=2025 saida=resultados/mc_piloto.csv procs=4

  * grava uma linha por realizacao no CSV assim que ela termina (retomada: realizacoes ja gravadas sao puladas);
  * no fim (e com resumo=1 so le o CSV) imprime P_f = P(Gamma < 1), IC de 95% (Wilson), CoV(P_f) = sqrt((1 - P_f)/
    (S P_f)) (Eq. 68), media e desvio de Gamma.
Outras opcoes: h=0.05 (passo da grade dos campos), hD=0.25 (malha do Darcy na janela), resumo=1."""
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import csv  # noqa: E402
import json  # noqa: E402
import math  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402

import numpy as np  # noqa: E402

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
import campos as cp  # noqa: E402
import cinematico as kin  # noqa: E402
import rapido as rp  # noqa: E402

COLUNAS = ["real", "Gamma", "classe", "Gamma_I", "Gamma_II", "xA", "xB", "yB", "y_max", "sB", "t2", "rB",
           "avaliacoes", "tempo", "t_campos", "t_darcy", "t_busca", "n_normalidade", "n_phi_nao_conv", "fallbacks",
           "c_curva", "phi_curva", "c_bloco", "phi_bloco", "k_bloco", "max_varr", "n_max_regra",
           "semente", "h", "hD", "M"]
# parametros que definem a amostra: gravados em cada linha e conferidos na retomada
CHAVES = ("semente", "h", "hD", "M")

_P = None
_PART = None


def _args(argv):
    a = dict(n="200", primeira="0", semente="2025", saida=os.path.join(AQUI, "resultados", "mc_piloto.csv"),
             procs="4", h="0.05", hD="0.25", resumo="0")
    for s in argv:
        k, _, v = s.partition("=")
        if k not in a:
            raise SystemExit(f"opcao desconhecida: {k}")
        a[k] = v
    return a


def _arquivo_partida(h, hD):
    return os.path.join(AQUI, "dados", f"partida_det_h{h}_hD{hD}.json")


def _assinatura(P):
    """O que define o caso deterministico (o mecanismo de partida so vale para ele)."""
    t = P.dom
    return dict(h=P.h, hD=P.hD_janela, M=P.kl.M, hw=P.hw, gamma_ef=P.gamma_ef, H=P.H, beta=P.beta,
                stats={k: list(v) for k, v in P.stats.items()}, dominio=[t.x_esq, t.x_dir, t.y_base],
                janela=[P.tal.x_esq, P.tal.x_dir, P.tal.y_base])


def _confere_partida(P, arq):
    """Partida guardada: confere a assinatura; se nao bater (ou faltar), recalcula e regrava."""
    if os.path.exists(arq):
        with open(arq) as f:
            d = json.load(f)
        if d.get("assinatura") == json.loads(json.dumps(_assinatura(P))):
            return
        print(f"AVISO: {arq} e de outro caso (assinatura diferente): recalculando", flush=True)
        os.remove(arq)
    rp.partidas_deterministicas(P, arq)
    with open(arq) as f:
        d = json.load(f)
    d["assinatura"] = _assinatura(P)
    with open(arq, "w") as f:
        json.dump(d, f, indent=1)


def _inicia(h, hD, arq_partida):
    global _P, _PART
    try:
        os.nice(10)
    except OSError:
        pass
    _P = cp.Problema(h=h, hD_janela=hD)
    _P.hD_janela = hD
    _PART, _ = rp.partidas_deterministicas(_P, arq_partida)


def dentro_do_bloco(m, talude, x, y):
    """Mascara dos pontos (x, y) dentro do bloco do mecanismo m (pela quadratura polar: intervalos de rho)."""
    rho = np.hypot(x - m.Cx, y + m.Cy)
    th = np.arctan2(y + m.Cy, m.Cx - x)
    ray, ra, rb, _, _ = kin._cruzamentos(th, m, talude)
    dentro = np.zeros(len(x), bool)
    ok = (rho[ray] >= ra) & (rho[ray] <= rb)
    dentro[ray[ok]] = True
    return dentro


def uma_realizacao(i, semente):
    P = _P
    t0 = time.perf_counter()
    d = P.realizacao(i, semente)
    cg = rp.CampoGrade(P.xg[0], P.h, d["c"], d["phi"], d["u"], P.gamma_ef, P.tal)
    rp.EST.zera()
    t1 = time.perf_counter()
    lin = dict(real=i, t_campos=d["t_campos"], t_darcy=d["t_darcy"], semente=semente, h=P.h, hD=P.hD_janela,
               M=P.kl.M)
    try:
        r = rp.minimizar_rapido(cg, P.tal, partidas=_PART)
    except RuntimeError:
        lin.update(Gamma=float("nan"), classe="", tempo=time.perf_counter() - t0)
        print(f"AVISO: realizacao {i}: nenhum mecanismo valido (Gamma = NaN)", flush=True)
        return lin
    t2 = time.perf_counter()
    m = r.mecanismo
    X, Y = m.pontos(2)
    lin.update(Gamma=r.Gamma, classe=m.classe, Gamma_I=r.por_classe.get("I", kin.Mecanismo(0, 0, 0)).Gamma,
               Gamma_II=r.por_classe.get("II", kin.Mecanismo(0, 0, 0)).Gamma, xA=m.xA, xB=m.xB, yB=m.yB,
               y_max=float(Y.max()), sB=m.sB, t2=m.t2, rB=m.rB, avaliacoes=r.n_aval, t_busca=t2 - t1,
               n_normalidade=r.n_normalidade, n_phi_nao_conv=r.n_phi_nao_conv, fallbacks=rp.EST.fallback,
               max_varr=rp.EST.max_varr, n_max_regra=rp.EST.n_max_regra)
    # medias de c e phi na linha de ruptura (pontos de Gauss, arcos de mesmo comprimento) e no bloco
    xi, _ = kin._GL[2]
    a, b = m.th[:-1], m.th[1:]
    tg = 0.5 * (a + b)[:, None] + 0.5 * (b - a)[:, None] * xi[None, :]
    rg = m.r[:-1, None] * np.exp((tg - a[:, None]) * m.tp[:, None])
    xg, yg = m.Cx - rg * np.cos(tg), -m.Cy + rg * np.sin(tg)
    lin["c_curva"] = float(cg.c(xg, yg).mean())
    lin["phi_curva"] = float(np.degrees(cg.phi(xg, yg).mean()))
    Xg, Yg = np.meshgrid(P.xg, P.yg, indexing="ij")
    caixa = (Xg >= X.min()) & (Xg <= m.xB) & (Yg <= Y.max())
    xs, ys = Xg[caixa], Yg[caixa]
    den = dentro_do_bloco(m, P.tal, xs, ys)
    if np.any(den):
        lin["c_bloco"] = float(d["c"][caixa][den].mean())
        lin["phi_bloco"] = float(np.degrees(d["phi"][caixa][den].mean()))
    cen = P.xyD[P.triD].mean(axis=1)
    caixa = (cen[:, 0] >= X.min()) & (cen[:, 0] <= m.xB) & (cen[:, 1] <= Y.max())
    den = dentro_do_bloco(m, P.tal, cen[caixa, 0], cen[caixa, 1])
    if np.any(den):
        lin["k_bloco"] = float(d["k"][caixa][den].mean())
    lin["tempo"] = time.perf_counter() - t0
    return lin


def _tarefa(args):
    return uma_realizacao(*args)


def _linha_ok(r):
    if None in r or any(r.get(c) is None for c in COLUNAS[:2]):
        return False
    try:
        int(r["real"])
        float(r["Gamma"])
    except (TypeError, ValueError):
        return False
    return True


def ler(saida):
    """Linhas completas do CSV (uma linha cortada por um processo morto e descartada)."""
    if not os.path.exists(saida):
        return []
    with open(saida, newline="") as f:
        linhas = list(csv.DictReader(f))
    boas = [r for r in linhas if _linha_ok(r)]
    if len(boas) < len(linhas):
        print(f"AVISO: {len(linhas) - len(boas)} linha(s) incompleta(s) descartada(s) em {saida}", flush=True)
    return boas


def _corta_ate_ultima_linha(saida):
    """Antes de reabrir em modo 'a': tira um pedaco de linha sem '\\n' no fim (processo morto no meio da escrita)."""
    if not os.path.exists(saida) or os.path.getsize(saida) == 0:
        return
    with open(saida, "rb+") as f:
        dados = f.read()
        if dados.endswith(b"\n"):
            return
        k = dados.rfind(b"\n")
        f.truncate(k + 1 if k >= 0 else 0)
        print(f"AVISO: linha cortada no fim de {saida} removida", flush=True)


def resumo(linhas, texto=True):
    G = np.array([float(r["Gamma"]) for r in linhas])
    nan = [int(r["real"]) for r, g in zip(linhas, G) if not np.isfinite(g)]
    if nan:
        print(f"AVISO: {len(nan)} realizacao(oes) sem mecanismo valido FORA de S: {nan[:20]}", flush=True)
    G = G[np.isfinite(G)]
    S = len(G)
    nf = int(np.sum(G < 1.0))
    pf = nf / S if S else float("nan")
    z = 1.959963984540054
    den = 1 + z * z / S
    cen = (pf + z * z / (2 * S)) / den
    meia = z * math.sqrt(pf * (1 - pf) / S + z * z / (4 * S * S)) / den
    cov = math.sqrt((1 - pf) / (S * pf)) if nf > 0 else float("inf")
    tempos = np.array([float(r["tempo"]) for r in linhas if r.get("tempo") not in (None, "")])
    res = dict(S=S, n_falhas=nf, Pf=pf, IC95=(cen - meia, cen + meia), CoV_Pf=cov, media_Gamma=float(G.mean()),
               desvio_Gamma=float(G.std(ddof=1)) if S > 1 else float("nan"), tempo_medio=float(tempos.mean()),
               n_nan=len(linhas) - S)
    if texto:
        print(f"S = {S} realizacoes ({res['n_nan']} sem mecanismo valido), falhas (Gamma < 1): {nf}")
        print(f"P_f = {100 * pf:.2f}%   IC 95% (Wilson) = [{100 * res['IC95'][0]:.2f}%, {100 * res['IC95'][1]:.2f}%]"
              f"   CoV(P_f) = {100 * cov:.2f}%")
        print(f"Gamma: media {res['media_Gamma']:.4f}, desvio {res['desvio_Gamma']:.4f}   (artigo: P_f = 11,50%,"
              f" CoV 2,77% com S = 10 000; media 1,353, desvio 0,318)")
        print(f"tempo medio por realizacao (um processo): {res['tempo_medio']:.2f} s")
    return res


def main(argv):
    a = _args(argv)
    saida = os.path.abspath(a["saida"])
    if a["resumo"] == "1":
        resumo(ler(saida))
        return
    n, primeira, semente, procs = int(a["n"]), int(a["primeira"]), int(a["semente"]), int(a["procs"])
    h, hD = float(a["h"]), float(a["hD"])
    os.makedirs(os.path.dirname(saida), exist_ok=True)
    arq_partida = _arquivo_partida(h, hD)
    P0 = cp.Problema(h=h, hD_janela=hD)
    P0.hD_janela = hD
    _confere_partida(P0, arq_partida)
    del P0
    with open(arq_partida) as f:
        print(f"deterministico: Gamma = {json.load(f)['Gamma']:.4f}", flush=True)
    _corta_ate_ultima_linha(saida)
    ant = ler(saida)
    esperado = dict(semente=semente, h=h, hD=hD, M=cp.M_PADRAO)
    for r in ant:
        for c in CHAVES:
            v = r.get(c)
            if v in (None, "") or float(v) != float(esperado[c]):
                raise SystemExit(f"{saida} tem a realizacao {r['real']} com {c} = {v!r}, diferente de "
                                 f"{esperado[c]}: use outro arquivo de saida")
    feitas = {int(r["real"]) for r in ant}
    todo = [(i, semente) for i in range(primeira, primeira + n) if i not in feitas]
    print(f"{len(todo)} realizacoes a fazer ({len(feitas)} ja no arquivo), {procs} processos, semente {semente}",
          flush=True)
    novo = not os.path.exists(saida) or os.path.getsize(saida) == 0
    t0 = time.perf_counter()
    with open(saida, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLUNAS, extrasaction="ignore")
        if novo:
            w.writeheader()
            f.flush()
        if procs <= 1:
            _inicia(h, hD, arq_partida)
            it = map(_tarefa, todo)
            pool = None
        else:
            import multiprocessing as mp
            pool = mp.get_context("fork").Pool(procs, initializer=_inicia, initargs=(h, hD, arq_partida))
            it = pool.imap_unordered(_tarefa, todo)
        for k, lin in enumerate(it, 1):
            w.writerow(lin)
            f.flush()
            if k % 10 == 0 or k == len(todo):
                print(f"  {k}/{len(todo)}  ({time.perf_counter() - t0:.0f} s)", flush=True)
        if pool is not None:
            pool.close()
            pool.join()
    print(f"tempo de parede: {time.perf_counter() - t0:.0f} s", flush=True)
    resumo(ler(saida))


if __name__ == "__main__":
    main(sys.argv[1:])
