#!/usr/bin/env python3
"""conferencia.py: conferencias do Monte Carlo rapido (campos.py, rapido.py, mc.py).

    python3 -I conferencia.py campos [n=2000]       KL: media, CoV, correlacao amostrais x teoria; variancia retida
    python3 -I conferencia.py det                   deterministico: busca rapida x cinematico.minimizar (exato)
    python3 -I conferencia.py busca [n=30] [procs=4] [semente=2025] [primeira=0] [saida=...csv]
        em n realizacoes: busca rapida x cinematico.minimizar exato com os MESMOS campos (CampoMalha da mesma
        triangulacao), tabela em resultados/conferencia_busca.csv
"""
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import csv  # noqa: E402
import math  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402

import numpy as np  # noqa: E402

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
import campos as cp  # noqa: E402
import cinematico as kin  # noqa: E402
import rapido as rp  # noqa: E402
import mc  # noqa: E402


def _kw(argv, **pad):
    for s in argv:
        k, _, v = s.partition("=")
        pad[k] = v
    return pad


def conf_campos(argv):
    a = _kw(argv, n="2000", semente="2025")
    n, semente = int(a["n"]), int(a["semente"])
    P = cp.Problema()
    kl = P.kl
    print(f"KL: retangulo {kl.Dx:g} x {kl.Dy:g} m, Lx = {kl.Lx:g}, Ly = {kl.Ly:g}, M = {kl.M} termos, "
          f"{len(kl.ux)} modos em x e {len(kl.uy)} em y; fracao de variancia retida (media no retangulo) = "
          f"{100 * kl.fracao:.2f}%")
    pts = np.array([(-1.0, 1.0), (2.5, 3.0), (0.0, 5.0), (5.0, 6.0), (-10.0, 3.0), (-30.0, 20.0),
                    (-1.0 + 2.0, 1.0), (-1.0, 1.0 + 0.5)])
    base = cp.PontosKL(kl, pts[:, 0], pts[:, 1])
    var_trunc = kl.variancia(base.Fx, base.Gy)
    Z = np.empty((3, n, len(pts)))
    for i in range(n):
        xi = cp.sortear_xi(semente, i, kl.M)
        for c in range(3):
            Z[c, i] = base.Z(xi[c])
    nomes = ("c", "phi", "k")
    print(f"\n{n} realizacoes (semente {semente}); pontos (x, y) com y para baixo a partir da crista")
    print(f"{'ponto':>14s} {'var KL':>7s} | " + " | ".join(f"{nm + ' media':>9s} {'(alvo)':>6s} {'CoV':>6s} {'alvo':>5s}"
                                                         f" {'trunc':>5s}" for nm in nomes))
    for p in range(len(pts)):
        cols = []
        for c, nm in enumerate(nomes):
            v = P.lognormal[nm](Z[c, :, p])
            media, cov = P.stats[nm]
            s = P.lognormal[nm].sigma_ln
            # CoV esperado com a variancia truncada: sqrt(exp(s^2 var) - 1); a media cai para exp(mu + s^2 var/2)
            cov_t = math.sqrt(math.exp(s * s * var_trunc[p]) - 1)
            cols.append(f"{v.mean():9.3f} {media:6.2f} {v.std(ddof=1) / v.mean():6.3f} {cov:5.2f} {cov_t:5.3f}")
        print(f"({pts[p, 0]:5.1f},{pts[p, 1]:5.1f}) {var_trunc[p]:7.4f} | " + " | ".join(cols))
    print("\ncorrelacao do campo gaussiano Z (o nucleo e o de Z; Eq. 64: exp(-|dx|/Lx - |dy|/Ly)):")
    for (p, q, rot) in ((0, 6, "dx = 2 m"), (0, 7, "dy = 0,5 m")):
        dx, dy = abs(pts[p, 0] - pts[q, 0]), abs(pts[p, 1] - pts[q, 1])
        teo = math.exp(-dx / kl.Lx - dy / kl.Ly)
        # o que o KL truncado da: sum lambda f(p) f(q) / sqrt(var_p var_q)
        fp = base.Fx[p, kl.cx] * base.Gy[p, kl.cy]
        fq = base.Fx[q, kl.cx] * base.Gy[q, kl.cy]
        trunc = float((kl.lam * fp * fq).sum() / math.sqrt(var_trunc[p] * var_trunc[q]))
        amostra = [np.corrcoef(Z[c, :, p], Z[c, :, q])[0, 1] for c in range(3)]
        lnam = [np.corrcoef(P.lognormal[nm](Z[c, :, p]), P.lognormal[nm](Z[c, :, q]))[0, 1]
                for c, nm in enumerate(nomes)]
        print(f"  {rot:>10s}: teoria {teo:.4f}, KL truncado {trunc:.4f}, amostral (c, phi, k) "
              + ", ".join(f"{v:.4f}" for v in amostra) + "   [lognormais: " + ", ".join(f"{v:.4f}" for v in lnam)
              + "]")
    print("  (erro padrao amostral de uma correlacao ~ (1 - rho^2)/sqrt(n) =",
          f"{(1 - 0.905 ** 2) / math.sqrt(n):.4f} e {(1 - 0.779 ** 2) / math.sqrt(n):.4f})")
    # correlacao cruzada entre campos (deve ser ~0)
    cc = [np.corrcoef(Z[0, :, 0], Z[1, :, 0])[0, 1], np.corrcoef(Z[0, :, 0], Z[2, :, 0])[0, 1],
          np.corrcoef(Z[1, :, 0], Z[2, :, 0])[0, 1]]
    print("  correlacao cruzada no ponto 1 (c-phi, c-k, phi-k):", ", ".join(f"{v:+.4f}" for v in cc))
    # variancia truncada na janela dos mecanismos
    gk = P.gradeKL_G
    vj = kl.variancia(np.repeat(gk.Fx, len(P.yg), axis=0), np.tile(gk.Gy, (len(P.xg), 1)))
    print(f"\nvariancia truncada na grade dos campos (janela x em [{P.xg[0]:g}, {P.xg[-1]:g}], y em [0, {P.yg[-1]:g}]):"
          f" media {vj.mean():.4f}, min {vj.min():.4f}, max {vj.max():.4f}")
    t = []
    for i in range(20):
        t0 = time.perf_counter()
        xi = cp.sortear_xi(semente, i, kl.M)
        P.gradeKL_G.Z(xi[0]); P.gradeKL_G.Z(xi[1]); P.gradeKL_D.Z(xi[2])
        t.append(time.perf_counter() - t0)
    print(f"tempo dos 3 campos (grade {P.forma_grade[0]} x {P.forma_grade[1]} nos + {len(P.gradeKL_D.Fx)} x "
          f"{len(P.gradeKL_D.Gy)} celulas do Darcy): {1e3 * np.median(t):.1f} ms")


def conf_det(argv):
    P = cp.Problema()
    d = P.campos_medios()
    cg = rp.CampoGrade(P.xg[0], P.h, d["c"], d["phi"], d["u"], P.gamma_ef, P.tal)
    t0 = time.perf_counter()
    cm = cg.para_malha()
    re = kin.minimizar(cm, P.tal)
    te = time.perf_counter() - t0
    rp.EST.zera()
    rr = rp.minimizar_rapido(cg, P.tal)               # sem partida quente (ela e o proprio deterministico)
    # Darcy direto na sua malha (sem passar pela grade)
    cD = kin.CampoMalha(P.xyD, P.triD, np.full(len(P.xyD), 10.0), np.full(len(P.xyD), math.radians(30.0)),
                        d["uD"], P.gamma_ef)
    rD = kin.minimizar(cD, P.dom)
    print(f"Darcy: {len(P.xyD)} nos, {len(P.triD)} triangulos; grade dos campos {P.forma_grade} (h = {P.h})")
    print(f"exato  (cinematico.minimizar, CampoMalha da grade): Gamma = {re.Gamma:.5f} ({re.mecanismo.classe}, "
          f"xA = {re.mecanismo.xA:.3f}), {re.n_aval} avaliacoes, {te:.1f} s")
    print(f"rapido (minimizar_rapido, sem partida quente):     Gamma = {rr.Gamma:.5f} ({rr.mecanismo.classe}, "
          f"xA = {rr.mecanismo.xA:.3f}), {rr.n_aval} avaliacoes, {rr.tempo_s:.2f} s")
    print(f"exato direto na malha do Darcy (sem a grade):      Gamma = {rD.Gamma:.5f}")


def _uma(args):
    i, semente = args
    P = mc._P
    d = P.realizacao(i, semente)
    cg = rp.CampoGrade(P.xg[0], P.h, d["c"], d["phi"], d["u"], P.gamma_ef, P.tal)
    rp.EST.zera()
    rr = rp.minimizar_rapido(cg, P.tal, partidas=mc._PART)
    fb, nmax = rp.EST.fallback, rp.EST.n_max_regra
    t0 = time.perf_counter()
    cm = cg.para_malha()
    tm = time.perf_counter() - t0
    re = kin.minimizar(cm, P.tal)
    return dict(real=i, G_rapido=rr.Gamma, G_exato=re.Gamma, rel=rr.Gamma / re.Gamma - 1,
                cl_rapido=rr.mecanismo.classe, cl_exato=re.mecanismo.classe,
                xA_rapido=rr.mecanismo.xA, xA_exato=re.mecanismo.xA, xB_rapido=rr.mecanismo.xB,
                xB_exato=re.mecanismo.xB, nav_rapido=rr.n_aval, nav_exato=re.n_aval,
                t_campos=d["t_campos"], t_darcy=d["t_darcy"], t_rapido=rr.tempo_s, t_exato=re.tempo_s,
                t_malha=tm, fallbacks=fb, n_max_regra=nmax, viol_rapido=rr.n_normalidade + rr.n_phi_nao_conv,
                viol_exato=re.n_normalidade + re.n_phi_nao_conv)


def conf_busca(argv):
    a = _kw(argv, n="30", procs="4", semente="2025", primeira="0",
            saida=os.path.join(AQUI, "resultados", "conferencia_busca.csv"))
    n, procs, semente, primeira = int(a["n"]), int(a["procs"]), int(a["semente"]), int(a["primeira"])
    h, hD = 0.05, 0.25
    arq = mc._arquivo_partida(h, hD)
    P0 = cp.Problema(h=h, hD_janela=hD)
    P0.hD_janela = hD
    mc._confere_partida(P0, arq)
    del P0
    import multiprocessing as mp
    with mp.get_context("fork").Pool(procs, initializer=mc._inicia, initargs=(h, hD, arq)) as pool:
        linhas = sorted(pool.imap_unordered(_uma, [(i, semente) for i in range(primeira, primeira + n)]),
                        key=lambda r: r["real"])
    saida = os.path.abspath(a["saida"])
    os.makedirs(os.path.dirname(saida), exist_ok=True)
    with open(saida, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(linhas[0]))
        w.writeheader()
        w.writerows(linhas)
    print(f"{'real':>4s} {'G rapido':>9s} {'G exato':>9s} {'dif rel':>9s} {'cl':>5s} {'xA r':>6s} {'xA e':>6s}"
          f" {'nav r':>5s} {'nav e':>5s} {'t r (s)':>7s} {'t e (s)':>7s}")
    for r in linhas:
        print(f"{r['real']:4d} {r['G_rapido']:9.5f} {r['G_exato']:9.5f} {r['rel']:+9.2e} {r['cl_rapido']:>2s}/"
              f"{r['cl_exato']:<2s} {r['xA_rapido']:6.2f} {r['xA_exato']:6.2f} {r['nav_rapido']:5d}"
              f" {r['nav_exato']:5d} {r['t_rapido']:7.2f} {r['t_exato']:7.2f}")
    rel = np.array([r["rel"] for r in linhas])
    print(f"\ndiferenca relativa rapido/exato - 1: max {rel.max():+.2e}, min {rel.min():+.2e}, media {rel.mean():+.2e},"
          f" media |.| {np.abs(rel).mean():.2e}")
    print(f"rapido acima do exato por mais de 0,5% (minimo perdido): {int(np.sum(rel > 0.005))} de {len(rel)};"
          f" mais de 0,1%: {int(np.sum(rel > 0.001))}")
    print(f"falha (Gamma < 1) classificada diferente: {int(sum((r['G_rapido'] < 1) != (r['G_exato'] < 1) for r in linhas))}")
    med = lambda k: np.mean([r[k] for r in linhas])  # noqa: E731
    print(f"tempos medios: campos {1e3 * med('t_campos'):.0f} ms, Darcy {1e3 * med('t_darcy'):.0f} ms, busca rapida"
          f" {med('t_rapido'):.2f} s ({med('nav_rapido'):.0f} avaliacoes, {1e3 * med('t_rapido') / med('nav_rapido'):.2f}"
          f" ms cada), exato {med('t_exato'):.1f} s ({med('nav_exato'):.0f} avaliacoes,"
          f" {1e3 * med('t_exato') / med('nav_exato'):.1f} ms cada)")
    print(f"marchas que cairam na sequencial: {sum(r['fallbacks'] for r in linhas)}; marchas com a regra 'so sobe'"
          f" (> 15 varreduras): {sum(r['n_max_regra'] for r in linhas)}; descartes por normalidade/phi_i:"
          f" rapido {sum(r['viol_rapido'] for r in linhas)}, exato {sum(r['viol_exato'] for r in linhas)}")
    print("salvo em", saida)


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    {"campos": conf_campos, "det": conf_det, "busca": conf_busca}.get(
        cmd, lambda a: print(__doc__))(sys.argv[2:])
