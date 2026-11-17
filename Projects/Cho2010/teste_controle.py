#!/usr/bin/env python3
"""teste_controle.py: a variavel de controle funciona com o que se tem na pratica? Um teste com as realizacoes ja feitas,
sem nenhuma analise nova de elementos finitos (r = -0,5, Cho2010).

O diagnostico_ai.py ajusta o substituto (ln FS = a + b . xi, nos 200 primeiros xi de cada campo) e o probit no direto
inteiro (25 000 realizacoes), que na pratica nao existe: se existisse, nao precisariamos da AI. Aqui:

1. O direto em lotes de 1000 realizacoes (o tamanho de uma AI). Em cada lote, tres estimativas:
   - P: a media de I, o Monte Carlo direto do lote;
   - P_cv realista: o substituto ajustado so no proprio lote, com validacao cruzada em 5 dobras (o de cada dobra sai das
     outras 4, com 800 realizacoes), nos K primeiros xi de cada campo (2K + 1 coeficientes);
   - P_cv ideal: o substituto ajustado fora do lote (as outras 24 000 realizacoes), com K = 200.
   O desvio padrao entre os lotes e a dispersao verdadeira de cada estimador; a razao das variancias e o ganho
   verdadeiro. A media dos lotes contra a P_f do direto inteiro mostra o vies.
2. As AI de r = -0,5 (as mesmas do diagnostico_ai.py), com o substituto ajustado so na propria amostra da AI (5 dobras,
   K termos por campo) e, para comparar, no direto inteiro (o ideal). Na AI, o controle de cada realizacao e
   H_j(xi) w - E[H_j], com H_j ajustado sem a dobra j: tem media zero sob h, entao a estimativa continua sem vies.
Uso: python3 -I teste_controle.py  (grava <runs>/cho2010/diagnostico/teste_controle.json e a tabela na tela)"""
import glob
import json
import math
import os
import runpy

import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
D = runpy.run_path(os.path.join(AQUI, "diagnostico_ai.py"), run_name="diagnostico")
ajusta, h_de, media_exata = D["ajusta"], D["h_de"], D["media_exata"]
CACHE, E1 = D["CACHE"], D["E1"]
KS = (10, 20, 50)  # termos por campo no substituto realista
LOTE = 1000


def colunas(K, Ktot=200):
    return list(range(K)) + list(range(Ktot, Ktot + K))


def cv_dobras(X, fs, w, dobras=5, semente=12345):
    """P_cv com o substituto ajustado fora de cada dobra; w = 1 no direto. Devolve (P_cv, desvio, rho)"""
    n = len(fs)
    y = (fs < 1) * w
    dob = np.random.default_rng(semente).permutation(n) % dobras
    c = np.empty(n)  # o controle: H_j(xi) w - E[H_j], media zero
    for j in range(dobras):
        sub = ajusta(X[dob != j], fs[dob != j])
        c[dob == j] = h_de(sub, X[dob == j]) * w[dob == j] - media_exata(sub)
    g = np.cov(y, c)[0, 1] / c.var(ddof=1)
    p = y.mean() - g * c.mean()
    return float(p), float((y - g * c).std(ddof=1) / math.sqrt(n)), float(np.corrcoef(y, c)[0, 1])


def com_sub(sub, X, fs, w):
    """P_cv com um substituto ja ajustado (fora da amostra)"""
    y = (fs < 1) * w
    c = h_de(sub, X) * w - media_exata(sub)
    g = np.cov(y, c)[0, 1] / c.var(ddof=1)
    return float(y.mean() - g * c.mean()), float((y - g * c).std(ddof=1) / math.sqrt(len(y))), \
        float(np.corrcoef(y, c)[0, 1])


def direto_r05():
    R = D["realizacoes"](f"{E1}/r05_direto")
    reais = np.array(sorted(k for k, l in R.items() if l["situacao"] == "ok"))
    fs = np.array([float(R[k]["FS_min"]) for k in reais])
    z = np.load(f"{CACHE}/xi_direto_2010_25000_24999.npz")  # os xi das realizacoes 0 a 24 999 (semente 2010)
    assert reais.max() < 25000
    return reais, fs, z["XA"][reais].reshape(len(reais), -1), z["XN"][reais].reshape(len(reais), -1)


def main():
    reais, fs, XA, XN = direto_r05()
    P = float((fs < 1).mean())
    nl = len(fs) // LOTE
    print(f"direto r = -0,5: {len(fs)} realizacoes, P_f = {100 * P:.3f}%; {nl} lotes de {LOTE}")
    est = {"P": []}
    est.update({f"realista K={K}": [] for K in KS})
    est["ideal K=200"] = []
    um = np.ones(LOTE)
    for b in range(nl):
        i = np.arange(b * LOTE, (b + 1) * LOTE)
        fora = np.setdiff1d(np.arange(len(fs)), i)
        p = float((fs[i] < 1).mean())
        est["P"].append((p, math.sqrt(max(p * (1 - p), 1e-12) / LOTE), 0.))
        for K in KS:
            est[f"realista K={K}"].append(cv_dobras(XA[i][:, colunas(K)], fs[i], um))
        est["ideal K=200"].append(com_sub(ajusta(XA[fora], fs[fora]), XA[i], fs[i], um))
        print(f"  lote {b:2d}: P {100 * p:.2f}%  " + "  ".join(f"{k} {100 * v[-1][0]:.2f}%" for k, v in est.items()
                                                                  if k != "P"), flush=True)
    var0 = np.var([v[0] for v in est["P"]], ddof=1)
    lotes = {}
    print(f"\n  {'estimador':18s} {'media':>7s} {'desvio verdadeiro':>18s} {'desvio estimado':>16s} {'ganho':>6s} "
          f"{'rho':>5s}")
    for k, v in est.items():
        a = np.array(v)
        lotes[k] = dict(media=float(a[:, 0].mean()), desvio=float(a[:, 0].std(ddof=1)),
                        desvio_estimado=float(np.sqrt((a[:, 1] ** 2).mean())), ganho=float(var0 / a[:, 0].var(ddof=1)),
                        rho=float(a[:, 2].mean()), por_lote=a[:, 0].tolist())
        L = lotes[k]
        print(f"  {k:18s} {100 * L['media']:6.3f}% {100 * L['desvio']:17.3f}% {100 * L['desvio_estimado']:15.3f}% "
              f"{L['ganho']:6.2f} {L['rho']:5.2f}")
    print(f"  (P_f do direto inteiro: {100 * P:.3f}%; desvio da media de {nl} lotes de P: "
          f"{100 * math.sqrt(P * (1 - P) / LOTE / nl):.3f}%)")

    # 2. as AI de r = -0,5, com o substituto ajustado na propria amostra
    sub_kl, sub_nys = ajusta(XA, fs), ajusta(XN, fs)
    ais = []
    print(f"\n  AI de r = -0,5 (P_f do direto {100 * P:.3f}%): P, e P_cv com o substituto da propria amostra (K termos) "
          f"e do direto inteiro")
    for r, base, rot, pasta in D["casos"]():
        if r != "r05" or len(glob.glob(f"{pasta}/mc_[0-9]*_[0-9]*.csv")) == 0:
            continue
        if len(D["realizacoes"](pasta)) < 1000:  # so as AI completas (a fase 4 parcial tem ~65 por semente)
            continue
        par = D["parametros"](pasta)
        sem = int(par["SEMENTE"])
        R = D["realizacoes"](pasta)
        rr = np.array(sorted(k for k, l in R.items() if l["situacao"] == "ok"))
        mu = D["desloc"](pasta)
        k = mu.shape[1]
        z = D["xi"](sem, rr, max(200, k))
        z[:, :, :k] += mu
        w = np.exp(np.array([float(R[j]["log_peso"]) for j in rr]))
        fsa = np.array([float(R[j]["FS_min"]) for j in rr])
        Z = z[:, :, :200].reshape(len(rr), -1)
        if base == "nys":  # o substituto ideal da base de Nystrom e o do direto projetado (XN)
            ideal = com_sub(sub_nys, Z, fsa, w)
        else:
            ideal = com_sub(sub_kl, Z, fsa, w)
        y = (fsa < 1) * w
        linha = dict(rotulo=rot, base=base, semente=sem, pasta=pasta, n=len(rr), p=float(y.mean()),
                     sd=float(y.std(ddof=1) / math.sqrt(len(y))), ideal=ideal,
                     realista={K: cv_dobras(Z[:, colunas(K)], fsa, w) for K in KS})
        ais.append(linha)
        print(f"  {base:3s} s{sem:<2d} {rot[:34]:34s} P {100 * linha['p']:.3f} +- {100 * linha['sd']:.3f}   " +
              "  ".join(f"K={K}: {100 * v[0]:.3f} +- {100 * v[1]:.3f}" for K, v in linha["realista"].items()) +
              f"   ideal {100 * ideal[0]:.3f} +- {100 * ideal[1]:.3f}", flush=True)
    with open(f"{CACHE}/teste_controle.json", "w") as f:
        json.dump(dict(P_direto=P, n_direto=len(fs), lote=LOTE, lotes=lotes, ai=ais), f, indent=1)


if __name__ == "__main__":
    main()
