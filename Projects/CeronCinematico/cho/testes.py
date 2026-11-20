#!/usr/bin/env python3
"""testes.py: conferencias da analise limite no Cho.

    python3 -I testes.py campos [mef=testes_mef] [base=<base_nys2000.txt>]
        (1) c e phi de campos_cho contra o campos.csv que o Cho2010 grava (figuras=<n>): 2 realizacoes analiticas r = -0,5
            (semente 2010, reais 0 e 1), 1 analitica r = -0,7 (real 5) e 1 de Nystrom com o centro da ISMC (semente 21,
            real 0, form_r05_M20/ustar.txt); a grade da analise limite (forma tensorial / Psi de Nystrom) contra a formula
            ponto a ponto; media, CoV e correlacao amostrais em 400 realizacoes contra as teoricas.
    python3 -I testes.py det
        (2) deterministico: Gamma e FS (rapido na grade e cinematico.fs_classico) contra 1,777 e 1,203.
    python3 -I testes.py busca [n=30] [procs=4] [r=-0.5] [semente=2010]
        (3) rapido (minimizar_rapido) x exato (cinematico.minimizar na mesma triangulacao) nas realizacoes 0..n-1.
Os comandos do Cho2010 que geraram testes_mef/ estao no README.md."""
from __future__ import annotations

import csv
import math
import os
import sys
import time

import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
import campos_cho as cc  # noqa: E402
import cinematico as kin  # noqa: E402
import la  # noqa: E402
import rapido as rp  # noqa: E402

RUNS = os.path.expanduser("~/projects/runs")
BASE = f"{RUNS}/cho2010/nystrom2000/base_nys2000.txt"
USTAR = f"{RUNS}/fermi_etapa3b/etapa3b/form_r05_M20/ustar.txt"


def le_campos_csv(arq):
    rs = list(csv.DictReader(open(arq)))
    return tuple(np.array([float(r[k]) for r in rs]) for k in ("x", "y", "c", "phi"))


def t_campos(a):
    mef = a.get("mef", os.path.join(AQUI, "testes_mef"))
    base = a.get("base", BASE)
    casos = [("ana05_real0", "analitico", -0.5, 2010, 0, None), ("ana05_real1", "analitico", -0.5, 2010, 1, None),
             ("ana07_real5", "analitico", -0.7, 2010, 5, None), ("nys05_real0", "nystrom", -0.5, 21, 0, USTAR)]
    B = cc.le_nystrom(base)
    ok = True
    print("(1a) campo do Python x campos.csv do Cho2010 (pontos de 0,5 m no talude)")
    for pasta, campo, r, sem, real, centro in casos:
        arq = os.path.join(mef, pasta, "campos.csv")
        if not os.path.exists(arq):
            print(f"  {pasta}: falta {arq}")
            ok = False
            continue
        X, Y, c, phi = le_campos_csv(arq)
        xi = cc.xi_cpp(sem, [real])[0]
        if centro:
            xi = xi + cc.le_desloc(centro)
        cp, pp = cc.campos_em_pontos(xi, X, Y, r, campo, B)
        ec, ep = np.max(np.abs(cp / c - 1)), np.max(np.abs(pp / phi - 1))
        corr = np.corrcoef(np.log(c), np.log(phi))[0, 1]
        bom = ec < 1e-5 and ep < 1e-5
        ok &= bom
        print(f"  {pasta:12s} {campo:9s} r = {r:+.1f}: {len(X)} pontos, erro rel max c {ec:.1e}, phi {ep:.1e} "
              f"(precisao do CSV ~5e-6); corr(ln c, ln phi) espacial {corr:+.2f} -> {'OK' if bom else 'FALHOU'}")
    print("(1b) grade da analise limite (h = 0,1) x formula ponto a ponto, 500 nos de solo ao acaso")
    rng = np.random.default_rng(1)
    for campo in ("analitico", "nystrom"):
        C = cc.CamposCho(-0.7, campo, 0.1, base=B)
        xi = cc.xi_cpp(2010, [3])[0]
        cg, pg = C.campos(xi)
        xs, ys = np.meshgrid(C.x, C.y, indexing="ij")
        solo = (ys >= np.where(xs <= 0, 0, np.where(xs < 10, xs, 10))).ravel()
        idx = rng.choice(np.flatnonzero(solo), 500, replace=False)
        Xc, Yc = (xs + 10).ravel()[idx], (15 - ys).ravel()[idx]
        cp, pp = cc.campos_em_pontos(xi, Xc, Yc, -0.7, campo, B)
        ec = np.max(np.abs(cg.ravel()[idx] / cp - 1))
        ep = np.max(np.abs(np.degrees(pg.ravel()[idx]) / pp - 1))
        bom = ec < 1e-5 and ep < 1e-5
        ok &= bom
        print(f"  {campo:9s}: erro rel max c {ec:.1e}, phi {ep:.1e} -> {'OK' if bom else 'FALHOU'}")
    print("(1c) estatisticas amostrais (400 realizacoes, semente 2010, campo analitico) em 4 pontos")
    for r in (-0.5, -0.7):
        C = cc.CamposCho(r, "analitico", 0.1)
        pts = [(5.0, 10.0), (10.0, 5.0), (15.0, 12.0), (-5.0, 2.0)]   # (x, y) do Ceron
        ii = [(int(round((x + 10) / 0.1)), int(round(y / 0.1))) for x, y in pts]
        X = cc.xi_cpp(2010, range(400))
        vals = np.array([[(c[i, j], np.degrees(p[i, j])) for i, j in ii] for c, p in (C.campos(x) for x in X)])
        P = C.kl.psi(np.array([x + 10 for x, _ in pts]), np.array([15 - y for _, y in pts]))
        var = (P ** 2).sum(1)                            # variancia de Z retida em cada ponto
        for k, (x, y) in enumerate(pts):
            c, p = vals[:, k, 0], vals[:, k, 1]
            cov_t = math.sqrt(math.exp(math.log(1.09) * var[k]) - 1)
            covf_t = math.sqrt(math.exp(math.log(1.04) * var[k]) - 1)
            print(f"  r = {r:+.1f} ({x:+5.1f}, {y:4.1f}): media c {c.mean():6.2f} (teor. "
                  f"{10 / math.sqrt(1.09) * math.exp(0.5 * math.log(1.09) * var[k]):5.2f}), CoV c {c.std() / c.mean():.3f}"
                  f" ({cov_t:.3f}); media phi {p.mean():5.2f}, CoV phi {p.std() / p.mean():.3f} ({covf_t:.3f}); "
                  f"corr(ln c, ln phi) {np.corrcoef(np.log(c), np.log(p))[0, 1]:+.3f} ({r:+.1f})")
    print("TESTE DOS CAMPOS:", "PASSOU" if ok else "FALHOU")
    return ok


def t_det(a):
    part, G = la.partidas_det(0.1)
    C = cc.CamposCho(-0.5, "analitico", la.H_PADRAO)
    av = la.Avaliador(C, part)
    c = np.full(C.shape, 10.0)
    phi = np.full(C.shape, math.radians(30.0))
    t = time.perf_counter()
    r = av.minimizar(c, phi)
    tg = time.perf_counter() - t
    t = time.perf_counter()
    F, nf, _ = av.fs(c, phi, r)
    tf = time.perf_counter() - t
    fsc = kin.fs_classico(la.TALUDE, 10.0, math.radians(30.0), 20.0, classes=("I",), Fmin=1.0, Fmax=1.6)
    ok = abs(r.Gamma - 1.777) < 2e-3 and abs(F - 1.203) < 2e-3
    print(f"(2) deterministico: Gamma exato (cinematico.minimizar) {G:.4f}; rapido {r.Gamma:.4f} ({la.modo(r.mecanismo)},"
          f" {tg:.2f} s); FS rapido {F:.4f} ({nf} minimizacoes, {tf:.2f} s); FS classico (forma fechada) {fsc:.4f}."
          f" Artigo: 1,777 e 1,203 (Cho LEM 1,204; SRM do MEF 1,2044) -> {'PASSOU' if ok else 'FALHOU'}")
    return ok


_B = {}


def _uma(n):
    av = _B["av"]
    xi = cc.xi_cpp(_B["sem"], [n])[0]
    c, phi = av.cc.campos(xi)
    cg = la.grade(av.h, c, phi)
    rp.EST.zera()
    r0 = rp.minimizar_rapido(cg, la.TALUDE, partidas=av.part)
    rr = av.minimizar(c, phi)
    re = kin.minimizar(cg.para_malha(), la.TALUDE)
    return dict(real=n, G0=r0.Gamma, Gr=rr.Gamma, Ge=re.Gamma, rel=rr.Gamma / re.Gamma - 1, mr=la.modo(rr.mecanismo),
                me=la.modo(re.mecanismo), tr=rr.tempo_s, te=re.tempo_s, nr=rr.n_aval, ne=re.n_aval,
                fb=rp.EST.fallback, viol=rr.n_normalidade + rr.n_phi_nao_conv + re.n_normalidade + re.n_phi_nao_conv)


def t_busca(a):
    import multiprocessing as mp
    n, procs = int(a.get("n", 30)), int(a.get("procs", 4))
    r, sem = float(a.get("r", -0.5)), int(a.get("semente", 2010))
    _B.update(av=la.Avaliador(cc.CamposCho(r, "analitico", 0.1)), sem=sem)
    with mp.get_context("fork").Pool(procs) as pool:
        L = sorted(pool.map(_uma, range(n)), key=lambda d: d["real"])
    print(f"(3) rapido x exato, r = {r}, semente {sem}, realizacoes 0..{n - 1}")
    print("G rapido: la.Avaliador.minimizar (a busca rapida e, se Gamma < 1,5, a forte); G r0: so a busca rapida padrao")
    print(f"{'real':>4s} {'G r0':>9s} {'G rapido':>9s} {'G exato':>9s} {'dif rel':>9s} {'modo r/e':>10s} {'t r':>5s} "
          f"{'t e':>5s}")
    for d in L:
        print(f"{d['real']:4d} {d['G0']:9.5f} {d['Gr']:9.5f} {d['Ge']:9.5f} {d['rel']:+9.1e} {d['mr']:>4s}/{d['me']:<5s} "
              f"{d['tr']:5.2f} {d['te']:5.1f}")
    rel0 = np.array([d["G0"] / d["Ge"] - 1 for d in L])
    print(f"so a busca rapida padrao: dif rel max {rel0.max():+.1e}, acima de 0,1 %: {int(np.sum(np.abs(rel0) > 1e-3))}")
    rel = np.array([d["rel"] for d in L])
    print(f"dif rel rapido/exato - 1: max {rel.max():+.1e}, min {rel.min():+.1e}, media {rel.mean():+.1e}; acima de 0,1 %:"
          f" {int(np.sum(np.abs(rel) > 1e-3))} de {n}; falha classificada diferente: "
          f"{sum((d['Gr'] < 1) != (d['Ge'] < 1) for d in L)}; modo diferente: {sum(d['mr'] != d['me'] for d in L)}")
    print(f"tempo medio: rapido {np.mean([d['tr'] for d in L]):.2f} s ({np.mean([d['nr'] for d in L]):.0f} avaliacoes),"
          f" exato {np.mean([d['te'] for d in L]):.1f} s ({np.mean([d['ne'] for d in L]):.0f}); marchas na sequencial "
          f"{sum(d['fb'] for d in L)}; descartes (normalidade/phi_i) {sum(d['viol'] for d in L)}")
    perto = np.array([d["Ge"] < 1.5 for d in L])
    ok = np.all(np.abs(rel[perto]) < 2e-3) and np.all(rel > -2e-3)
    print(f"criterio: |dif| < 0,2 % onde Gamma < 1,5 ({int(perto.sum())} realizacoes) e o rapido nunca 0,2 % abaixo")
    print("TESTE RAPIDO x EXATO:", "PASSOU" if ok else "FALHOU")
    with open(os.path.join(AQUI, "resultados", f"busca_r{abs(r):g}.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(L[0]))
        w.writeheader()
        w.writerows(L)
    return ok


if __name__ == "__main__":
    os.makedirs(os.path.join(AQUI, "resultados"), exist_ok=True)
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    a = la.le_kv(sys.argv[2:])
    f = {"campos": t_campos, "det": t_det, "busca": t_busca}.get(cmd)
    if f is None:
        sys.exit(__doc__)
    sys.exit(0 if f(a) else 1)
