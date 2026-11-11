#!/usr/bin/env python3
"""malhas_campos.py: (A8) o mesmo campo (os mesmos xi) em duas malhas diferentes.

KL analitica (a do SlopeRandomField): o programa, com SO_CAMPOS=1 CAMPOS=1 REAL=0 e o modelo de compara.py, grava em
campos.csv c', phi' e k nos pontos de integracao de cada malha (ceron_h5_a45_foco375, 720 triangulos, e a malha
refinada da realizacao 0 do metodo de duas etapas, gerada pelo gmsh de malha_refinada.geo, 2544 triangulos). Aqui o
mesmo campo e refeito em Python (klfechada.py: forma fechada, xi do xoshiro256** do programa e transformacao lognormal)
nos pontos de cada malha, e os pontos de integracao que as duas malhas tem em comum sao comparados.
KL numerica (contraste): KL de Galerkin P1 em numpy (galerkin_p1.py) em cada uma das duas malhas, M = 200 termos, os
mesmos xi de c'; os campos sao comparados nos mesmos pontos (os pontos de integracao da malha refinada), antes e depois
de alinhar o sinal de cada modo.

Uso: malhas_campos.py [dados] [figura]   (padrao: os dois; "dados" precisa das rodadas em $KL_SCRATCH/kl_rel_tal)
Dados (../../SlopeRandomFieldRefina/relatorio/dados): kl_malha_<malha>_nos.csv, kl_malha_<malha>_tri.csv,
kl_malha_<malha>_campos.csv (pontos de integracao: x, y, area, c' do programa, c' da forma fechada em Python, c' da
KL P1 da propria malha e, na refinada, c' da KL P1 da malha foco375 nos mesmos pontos), kl_malha_produto_interno.csv
e kl_malha_resumo.csv."""
import os
import sys
import time

import numpy as np

import estilo as E
from estilo import plt

MALHAS = {"foco375": os.path.join(E.MALHAS, "ceron_h5_a45_foco375.msh"),
          "refinada": os.path.join(E.SCRATCH, "kl_rel_tal", "malhas", "refinada_r00000.msh")}
CAMPOS = {m: os.path.join(E.SCRATCH, "kl_rel_tal", "campos", m, "campos.csv") for m in MALHAS}
# as mesmas realizacoes rodadas pelo compara.py (a simulacao inteira, com o campo usado nela)
CAMPOS_MC = {"foco375": os.path.join(E.SCRATCH, "refina", "comparacao", "foco375", "r00000", "campos.csv"),
             "refinada": os.path.join(E.SCRATCH, "refina", "comparacao", "refinada", "r00000", "campos.csv")}
M_NUM = 200
MEDIA_C, COV_C = 10., 0.3  # c' do modelo (kPa)


def le_campos(arq):
    a = np.loadtxt(arq, delimiter=",", skiprows=1)
    return dict(ip=a[:, 0].astype(int), x=a[:, 1], y=a[:, 2], area=a[:, 3], c=a[:, 4], phi=a[:, 5], k=a[:, 6])


def dados():
    from galerkin_p1 import avalia_p1, kl_p1, le_malha, quadratura
    from klfechada import KL2D, lognormal, xi_slope
    t0 = time.time()
    kl = KL2D(2000)
    xi = [xi_slope(c, real=0) for c in range(3)]
    mu_c, s_c = lognormal(MEDIA_C, COV_C)
    mu_p, s_p = lognormal(30., 0.1)
    mu_k, s_k = lognormal(1., 0.6)
    res = {}
    linhas = []
    for m in MALHAS:
        d = le_campos(CAMPOS[m])
        # forma fechada em Python nos mesmos pontos
        F, G = kl.fx(d["x"]), kl.gy(d["y"])
        z = [np.zeros(len(d["x"])) for _ in range(3)]
        for k in range(kl.M):
            fg = F[kl.ix[k]] * G[kl.iy[k]] * np.sqrt(kl.lam[k])
            for c in range(3):
                z[c] += xi[c][k] * fg
        d["c_py"] = np.exp(mu_c + s_c * z[0])
        d["phi_py"] = np.exp(mu_p + s_p * z[1])
        d["k_py"] = np.exp(mu_k + s_k * z[2])
        d["z_c"] = z[0]
        er = {v: np.max(np.abs(d[v] - d[v + "_py"]) / np.abs(d[v + "_py"])) for v in ("c", "phi", "k")}
        # o campo da simulacao inteira (compara.py) na mesma malha
        mc = le_campos(CAMPOS_MC[m]) if os.path.isfile(CAMPOS_MC[m]) else None
        dif_mc = (np.max(np.abs(mc["c"] - d["c"])) if mc is not None and len(mc["c"]) == len(d["c"]) else np.nan)
        X, T = le_malha(MALHAS[m])
        res[m] = (d, X, T)
        linhas.append(("pontos_%s" % m, len(d["x"]), "triangulos", len(T), "nos", len(X)))
        linhas.append(("analitica_erro_rel_max_programa_python_%s" % m, er["c"], er["phi"], er["k"], "dif_c_simulacao_compara",
                       dif_mc))
        print("%s: %d pontos; programa contra a forma fechada em Python: erro relativo maximo c' %.1e, phi' %.1e, k %.1e;"
              " contra o campo da simulacao do compara.py: %s" % (m, len(d["x"]), er["c"], er["phi"], er["k"], dif_mc))
    # pontos de integracao em comum
    from scipy.spatial import cKDTree
    a, b = res["foco375"][0], res["refinada"][0]
    dist, idx = cKDTree(np.c_[a["x"], a["y"]]).query(np.c_[b["x"], b["y"]])
    comum = dist < 1e-6
    dc = np.abs(b["c"][comum] - a["c"][idx[comum]])
    linhas.append(("analitica_pontos_comuns", int(comum.sum()), "max_dif_c_kPa", dc.max() if comum.any() else np.nan,
                   "max_dif_phi", np.max(np.abs(b["phi"][comum] - a["phi"][idx[comum]])) if comum.any() else np.nan))
    print("pontos de integracao comuns as duas malhas: %d; maior diferenca de c' %.3g kPa" % (comum.sum(), dc.max()))
    # KL numerica P1 em cada malha, mesmos xi (os de c'), M_NUM termos
    num = {}
    for m in MALHAS:
        d, X, T = res[m]
        t = time.time()
        lam, V, B = kl_p1(X, T, M_NUM)
        num[m] = (lam, V)
        print("KL P1 %s: %d nos, %.1f s; lambda_1..3 %s; 1 - sum_{k<=%d} lambda_k / area = %.4f" %
              (m, len(X), time.time() - t, np.round(lam[:3], 4), M_NUM, 1 - lam.sum() / 187.5))
        linhas.append(("p1_%s" % m, len(X), "erro_mq_M%d" % M_NUM, 1 - lam.sum() / 187.5, "lambda1", lam[0]))
    # produto interno dos modos das duas malhas, nos pontos de Dunavant da malha refinada
    Xb, Tb = res["refinada"][1], res["refinada"][2]
    xq, wq, Phi = quadratura(Xb, Tb)
    Xa, Ta = res["foco375"][1], res["foco375"][2]
    phiA = avalia_p1(Xa, Ta, num["foco375"][1], xq)  # (nq, M)
    phiB = Phi @ num["refinada"][1]
    pi = np.einsum("q,qk,qk->k", wq, phiA, phiB)
    sinal = np.where(pi < 0, -1., 1.)
    with open(os.path.join(E.DADOS, "kl_malha_produto_interno.csv"), "w") as f:
        f.write("k,lambda_foco375,lambda_refinada,produto_interno\n")
        for k in range(M_NUM):
            f.write("%d,%.6f,%.6f,%.6f\n" % (k + 1, num["foco375"][0][k], num["refinada"][0][k], pi[k]))
    # campos numericos nos pontos de integracao de cada malha e, na refinada, o da foco375 nos mesmos pontos
    xi_c = xi[0][:M_NUM]
    for m in MALHAS:
        d, X, T = res[m]
        lam, V = num[m]
        P = np.c_[d["x"], d["y"]]
        d["z_num"] = avalia_p1(X, T, V @ (np.sqrt(lam) * xi_c), P)
    d = res["refinada"][0]
    lamA, VA = num["foco375"]
    P = np.c_[d["x"], d["y"]]
    d["z_numA"] = avalia_p1(Xa, Ta, VA @ (np.sqrt(lamA) * xi_c), P)
    d["z_numA_al"] = avalia_p1(Xa, Ta, VA @ (np.sqrt(lamA) * xi_c * sinal), P)
    w = d["area"]

    def rms(v):
        return np.sqrt(np.sum(w * v * v) / np.sum(w))
    r_sem = rms(d["z_num"] - d["z_numA"]) / rms(d["z_num"])
    r_com = rms(d["z_num"] - d["z_numA_al"]) / rms(d["z_num"])
    corr = np.corrcoef(d["z_num"], d["z_numA_al"])[0, 1]
    corr0 = np.corrcoef(d["z_num"], d["z_numA"])[0, 1]
    linhas.append(("p1_rms_dif_relativa", r_sem, "alinhado", r_com, "correlacao_alinhado", corr))
    linhas.append(("p1_correlacao_sem_alinhar", corr0, "modos_com_sinal_trocado", int((pi < 0).sum()),
                   "produto_interno_min_abs_k1a20", np.abs(pi[:20]).min()))
    print("KL P1, mesmos xi: RMS(Z_ref - Z_foco375) / RMS(Z_ref) = %.3f sem alinhar, %.3f com o sinal alinhado "
          "(correlacao %.3f; %d modos com sinal trocado)" % (r_sem, r_com, corr, (pi < 0).sum()))
    for m in MALHAS:
        d, X, T = res[m]
        np.savetxt(os.path.join(E.DADOS, "kl_malha_%s_nos.csv" % m), X, delimiter=",", fmt="%.6f", header="x,y",
                   comments="")
        np.savetxt(os.path.join(E.DADOS, "kl_malha_%s_tri.csv" % m), T, delimiter=",", fmt="%d", header="n0,n1,n2",
                   comments="")
        col = [d["x"], d["y"], d["area"], d["c"], d["c_py"], np.exp(mu_c + s_c * d["z_num"])]
        cab = "x,y,area,c_programa,c_python,c_p1_propria"
        fmt = ["%.6f", "%.6f", "%.6e", "%.7g", "%.10g", "%.6g"]
        if m == "refinada":
            col += [np.exp(mu_c + s_c * d["z_numA"]), np.exp(mu_c + s_c * d["z_numA_al"])]
            cab += ",c_p1_foco375,c_p1_foco375_alinhado"
            fmt += ["%.6g", "%.6g"]
        np.savetxt(os.path.join(E.DADOS, "kl_malha_%s_campos.csv" % m), np.c_[tuple(col)], delimiter=",", fmt=fmt,
                   header=cab, comments="")
    with open(os.path.join(E.DADOS, "kl_malha_resumo.csv"), "w") as f:
        for l_ in linhas:
            f.write(",".join(str(v) for v in l_) + "\n")
    print("dados: %.1f s" % (time.time() - t0))


def media_por_triangulo(X, T, x, y, v):
    """media de v nos pontos (x, y) de cada triangulo"""
    from matplotlib.tri import Triangulation
    t = Triangulation(X[:, 0], X[:, 1], T).get_trifinder()(x, y)
    s = np.bincount(t, weights=v, minlength=len(T))
    n = np.bincount(t, minlength=len(T))
    return s / np.maximum(n, 1)


def figura():
    from matplotlib.tri import Triangulation
    dd = {}
    for m in MALHAS:
        X = np.loadtxt(os.path.join(E.DADOS, "kl_malha_%s_nos.csv" % m), delimiter=",", skiprows=1)
        T = np.loadtxt(os.path.join(E.DADOS, "kl_malha_%s_tri.csv" % m), delimiter=",", skiprows=1).astype(int)
        c = np.genfromtxt(os.path.join(E.DADOS, "kl_malha_%s_campos.csv" % m), delimiter=",", names=True)
        dd[m] = (X, T, c)
    pi = np.genfromtxt(os.path.join(E.DADOS, "kl_malha_produto_interno.csv"), delimiter=",", names=True)
    todos = np.concatenate([dd[m][2]["c_programa"] for m in MALHAS] + [dd[m][2]["c_p1_propria"] for m in MALHAS])
    vmin, vmax = np.percentile(todos, [0.5, 99.5])
    fig = plt.figure(figsize=(E.LARGURA, 6.3))
    gs = fig.add_gridspec(3, 6, height_ratios=[1, 1, 1.45], hspace=0.55, wspace=1.1)
    nomes = {"foco375": "foco375 (720 triângulos)", "refinada": "refinada, realização 0 (2544 triângulos)"}
    for lin, (col, titulo) in enumerate((("c_programa", "KL analítica (SlopeRandomField, M = 2000)"),
                                          ("c_p1_propria", "KL numérica P1 na própria malha (M = %d)" % M_NUM))):
        for j, m in enumerate(MALHAS):
            X, T, c = dd[m]
            ax = fig.add_subplot(gs[lin, 3 * j:3 * j + 3])
            v = media_por_triangulo(X, T, c["x"], c["y"], c[col])
            tri = Triangulation(X[:, 0], X[:, 1], T)
            im = ax.tripcolor(tri, facecolors=v, cmap="viridis", vmin=vmin, vmax=vmax, edgecolors="white",
                              linewidth=0.08, rasterized=True)
            E.contorno_talude(ax, lw=0.6)
            ax.set_aspect("equal")
            ax.grid(False)
            ax.set_xlim(-0.2, 25.2)
            ax.set_ylim(-0.2, 10.2)
            ax.set_xticks([0, 5, 10, 15, 20, 25])
            ax.set_yticks([0, 5, 10])
            ax.tick_params(labelsize=6.5, length=2)
            ax.set_title("(%s) %s: %s" % ("abcd"[2 * lin + j], titulo.split(" (")[0], nomes[m]), fontsize=7, loc="left",
                         pad=2)
        if lin == 0:
            cax = fig.add_axes([0.92, 0.555, 0.012, 0.36])
            cb = fig.colorbar(im, cax=cax)
            cb.set_label("c' (kPa)", fontsize=7)
            cb.ax.tick_params(labelsize=6.5)
            cb.ax.yaxis.set_major_formatter(E.fmt_virgula())
    # (e) programa contra a forma fechada em Python, nos pontos de cada malha
    ax = fig.add_subplot(gs[2, 0:2])
    for j, m in enumerate(MALHAS):
        c = dd[m][2]
        ax.plot(c["c_python"], c["c_programa"], "o", ms=1.6 if m == "refinada" else 2.4, color=E.PALETA[j], alpha=0.8,
                mew=0, label=m, rasterized=True)
    lim = [min(vmin, 4), max(vmax, 18)]
    ax.plot(lim, lim, color=E.TINTA2, lw=0.6)
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_aspect("equal")
    ax.set_xlabel("c' da forma fechada no ponto (kPa)", fontsize=7)
    ax.set_ylabel("c' do programa (kPa)", fontsize=7)
    er = max(np.max(np.abs(dd[m][2]["c_programa"] - dd[m][2]["c_python"])) for m in MALHAS)
    ax.set_title("(e) analítica: programa × função", fontsize=7, loc="left", pad=2)
    ax.text(0.04, 0.96, "máx. |dif.| = %s kPa" % E.cientifico(er, 1), transform=ax.transAxes, fontsize=6.5, va="top")
    ax.legend(loc="lower right", fontsize=6.3, markerscale=2.5, handletextpad=0.2, borderaxespad=0.2)
    E.virgula(ax)
    ax.tick_params(labelsize=6.5)
    # (f) KL P1: os campos das duas malhas nos mesmos pontos (os da malha refinada)
    ax = fig.add_subplot(gs[2, 2:4])
    c = dd["refinada"][2]
    ax.plot(c["c_p1_foco375"], c["c_p1_propria"], "o", ms=1.4, color=E.NEUTRO, alpha=0.5, mew=0, rasterized=True,
            label="sinal do solver")
    ax.plot(c["c_p1_foco375_alinhado"], c["c_p1_propria"], "o", ms=1.4, color=E.PALETA[1], alpha=0.6, mew=0,
            rasterized=True, label="sinal alinhado")
    ax.plot(lim, lim, color=E.TINTA2, lw=0.6)
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_aspect("equal")
    ax.set_xlabel("c' da KL P1 da foco375 (kPa)", fontsize=7)
    ax.set_ylabel("c' da KL P1 da refinada (kPa)", fontsize=7)
    ax.set_title("(f) numérica: mesmos pontos e mesmos ξ", fontsize=7, loc="left", pad=2)
    ax.legend(loc="lower right", fontsize=6.3, markerscale=2.5, handletextpad=0.2, borderaxespad=0.2)
    E.virgula(ax)
    ax.tick_params(labelsize=6.5)
    # (g) produto interno dos modos de mesmo indice das duas malhas
    ax = fig.add_subplot(gs[2, 4:6])
    ax.plot(pi["k"], np.abs(pi["produto_interno"]), "o", ms=2.2, color=E.PALETA[0], mew=0)
    ax.set_xscale("log")
    E.eixo_log(ax, "x")
    ax.set_ylim(0, 1.05)
    ax.set_xlim(0.8, M_NUM * 1.2)
    ax.yaxis.set_major_formatter(E.fmt_virgula())
    ax.set_xlabel("índice do modo k", fontsize=7)
    ax.set_ylabel(r"$|\int f_k^{foco375}\, f_k^{ref}\, d\Omega|$", fontsize=7)
    ax.set_title("(g) numérica: modos das duas malhas", fontsize=7, loc="left", pad=2)
    ax.tick_params(labelsize=6.5)
    E.salva(fig, "kl_malha_independencia", png=True)


if __name__ == "__main__":
    etapas = sys.argv[1:] or ["dados", "figura"]
    if "dados" in etapas:
        dados()
    if "figura" in etapas:
        figura()
