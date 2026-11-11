#!/usr/bin/env python3
"""cilindro.py: figuras, tabelas e dados da secao do meio cilindro de Mika et al. (2021), Secao 6.2, Exemplo 1-1.

Le os resultados do KLCilindroMika (CSV de uma linha por malha), a referencia de Nystrom da verificacao independente
(verif_kl/cil_referencia.npy e reta15_referencia.npy: linha 0 os 20 autovalores, linha 1 a dispersao absoluta entre
as quatro extrapolacoes), os autovalores guardados no notebook randomfied-v2.nb e as contas em numpy de
../confere_notebook.py (q2 com 3 e 5 pontos, s20 com 5). Grava em ../../SlopeRandomFieldRefina/relatorio:
  dados/kl_cil_neopz.csv          uma linha por malha e regra (malha exata) e as duas contas da malha do notebook
  dados/kl_cil_autovalores.csv    por modo: Mika Galerkin 1 e 2, Nystrom, NeoPZ 24x3x12 (5^3 pontos)
  dados/kl_cil_notebook.csv       por modo, malha 15x1x3: Nystrom no semi-anel poligonal, NeoPZ, numpy, notebook
  figuras/kl_cil_convergencia.pdf, kl_cil_espectro.pdf, kl_cil_geometria.pdf
  tabelas/kl_cil_autovalores.tex, kl_cil_convergencia.tex, kl_cil_notebook.tex
Uso: python3 cilindro.py   (KL_SCRATCH: diretorio com klmika/res, verif_kl, nb e kl_cil)
"""
import csv
import os

import numpy as np

import estilo as E
from estilo import plt

RES = os.path.join(E.SCRATCH, "klmika", "res")
VERIF = os.path.join(E.SCRATCH, "verif_kl")
VOLUME = 270 * np.pi

# Mika et al. (2021), Tabela 4: Galerkin isogeometrico com o nucleo exato (NURBS p = 2), 1050 e 6642 graus de liberdade
G1 = np.array([162.7991539, 91.43063070, 57.56702901, 51.09017918, 38.79740931, 27.90386143, 25.05611145, 19.36893419,
               16.15700369, 15.79530798, 15.14622914, 11.21328690, 10.17896939, 9.693564320, 8.054783080, 7.578085990,
               6.722378350, 6.445576690, 6.177345410, 5.763786920])
G2 = np.array([162.7965791, 91.42804062, 57.56447741, 51.08762278, 38.79483423, 27.90128438, 25.05356161, 19.36659412,
               16.15443088, 15.79273209, 15.14385016, 11.21090778, 10.17639037, 9.690982310, 8.052352020, 7.576621410,
               6.719925970, 6.443915210, 6.174771170, 5.761319370])
GL_G1, GL_G2 = 1050, 6642


def le_csv(nome):
    linhas = []
    for r in csv.DictReader(open(os.path.join(RES, nome))):
        linhas.append(dict(geometria=r["geometria"], malha="%sx%sx%s" % (r["nt"], r["nr"], r["nz"]),
                           elementos=int(r["elementos"]), equacoes=int(r["equacoes"]), h=float(r["h"]),
                           ng=int(r["ng"]), lam=np.array([float(r["lambda%d" % k]) for k in range(1, 21)]),
                           soma=float(r["soma"]), volume=float(r["volume_dominio"]), t_C=float(r["t_C"]),
                           t_aut=float(r["t_autovalores"]), t_total=float(r["t_total"]),
                           memoria=float(r["memoria_MB"])))
    return linhas


def numeros(caminho, linha=None):
    """os 20 numeros de uma linha de texto (a linha que tem 20 numeros, se linha for None)"""
    for s in open(caminho):
        p = s.split()
        if len(p) == 20:
            try:
                return np.array([float(v) for v in p])
            except ValueError:
                pass
    raise ValueError(caminho)


def gera_dados():
    """junta as fontes e grava os CSV pequenos em dados/"""
    cil = []
    for nome in ("exata_p2.csv", "exata_p2_24.csv", "exata_p2_ng5.csv", "exata_p2_24_ng5.csv"):
        cil += le_csv(nome)
    nb = le_csv("notebook_ng3.csv") + le_csv("notebook_ng5.csv")
    ref, disp = np.load(os.path.join(VERIF, "cil_referencia.npy"))
    ref_reta, disp_reta = np.load(os.path.join(VERIF, "reta15_referencia.npy"))
    guardado = np.array([float(s.split()[1]) for s in open(os.path.join(E.SCRATCH, "nb",
                                                                        "autovalores_notebook_cilindro.txt"))
                         if s.strip() and not s.startswith("#")])
    py_q2_3 = numeros(os.path.join(E.SCRATCH, "kl_cil", "q2_3.txt"))
    py_q2_5 = numeros(os.path.join(E.SCRATCH, "kl_cil", "q2_5.txt"))
    py_s20_5 = numeros(os.path.join(E.SCRATCH, "klmika", "py", "s20_5.txt"))

    cab = ["geometria", "malha", "elementos", "equacoes", "h", "ng"] + ["lambda%d" % k for k in range(1, 21)] + \
          ["soma", "volume_dominio", "t_C", "t_autovalores", "t_total", "memoria_MB"]
    with open(os.path.join(E.DADOS, "kl_cil_neopz.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(cab)
        for r in cil + nb:
            w.writerow([r["geometria"], r["malha"], r["elementos"], r["equacoes"], r["h"], r["ng"]] +
                       ["%.9f" % v for v in r["lam"]] +
                       ["%.9f" % r["soma"], "%.9f" % r["volume"], r["t_C"], r["t_aut"], r["t_total"], r["memoria"]])
    fino = [r for r in cil if r["malha"] == "24x3x12" and r["ng"] == 5][0]
    with open(os.path.join(E.DADOS, "kl_cil_autovalores.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["modo", "mika_galerkin1", "mika_galerkin2", "nystrom", "nystrom_dispersao", "neopz_24x3x12_ng5",
                    "notebook_guardado"])
        for k in range(20):
            w.writerow([k + 1, "%.10g" % G1[k], "%.10g" % G2[k], "%.8f" % ref[k], "%.2e" % disp[k],
                        "%.9f" % fino["lam"][k], "%.12g" % guardado[k]])
    nb3 = [r for r in nb if r["ng"] == 3][0]
    nb5 = [r for r in nb if r["ng"] == 5][0]
    with open(os.path.join(E.DADOS, "kl_cil_notebook.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["modo", "nystrom_semianel_poligonal", "nystrom_dispersao", "neopz_q2_ng3", "neopz_q2_ng5",
                    "numpy_q2_ng3", "numpy_q2_ng5", "numpy_s20_ng5", "notebook_guardado"])
        for k in range(20):
            w.writerow([k + 1, "%.8f" % ref_reta[k], "%.2e" % disp_reta[k], "%.9f" % nb3["lam"][k],
                        "%.9f" % nb5["lam"][k], "%.7f" % py_q2_3[k], "%.7f" % py_q2_5[k], "%.7f" % py_s20_5[k],
                        "%.12g" % guardado[k]])
    return dict(cil=cil, nb3=nb3, nb5=nb5, ref=ref, disp=disp, ref_reta=ref_reta, disp_reta=disp_reta,
                guardado=guardado, py_q2_3=py_q2_3, py_q2_5=py_q2_5, py_s20_5=py_s20_5, fino=fino)


def erros(lam, ref):
    return (lam - ref) / ref


def fig_convergencia(d):
    ref = d["ref"]
    fig, ax = plt.subplots(figsize=(E.LARGURA, 3.0))
    for ng, cor in ((3, E.PALETA[0]), (5, E.PALETA[1])):
        rs = sorted([r for r in d["cil"] if r["ng"] == ng], key=lambda r: r["equacoes"])
        n = np.array([r["equacoes"] for r in rs])
        emax = np.array([np.max(np.abs(erros(r["lam"], ref))) for r in rs])
        e1 = np.array([abs(erros(r["lam"], ref)[0]) for r in rs])
        ax.plot(n, emax, "-o", color=cor, ms=3.5, label=r"NeoPZ, $%d^3$ pontos: máx. em $k \leq 20$" % ng)
        ax.plot(n, e1, "--s", color=cor, ms=3, mfc="white", label=r"NeoPZ, $%d^3$ pontos: $\lambda_1$" % ng)
    for nome, g, gl, cor in (("Galerkin 1", G1, GL_G1, E.PALETA[2]), ("Galerkin 2", G2, GL_G2, E.PALETA[6])):
        e = np.max(np.abs(erros(g, ref)))
        ax.axhline(e, color=cor, lw=1.0, ls=":")
        ax.plot([gl], [e], marker="*", ms=8, color=cor, ls=":",
                label="Mika, %s (%d g.l.): máx. em $k \\leq 20$" % (nome, gl))
    # inclinacao h^4 (N ~ h^-3 => N^-4/3)
    n0 = np.array([1500., 8000.])
    ax.plot(n0, 3e-2 * (n0 / n0[0]) ** (-4. / 3.), color=E.TINTA2, lw=1.0, label=r"inclinação $h^4 \propto N^{-4/3}$")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(400, 10000)
    ax.set_ylim(5e-7, 6e-2)
    E.eixo_log(ax)
    ax.set_xlabel(r"número de equações $N$")
    ax.set_ylabel(r"erro relativo $(\lambda_k^h - \lambda_k)/\lambda_k$")
    ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), fontsize=7)
    E.salva(fig, "kl_cil_convergencia")


def fig_espectro(d):
    ref, k = d["ref"], np.arange(1, 21)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(E.LARGURA, 2.8))
    series = [("Nyström (referência)", ref, E.TINTA, "o", 6.5, "none"),
              ("Mika, Galerkin 1", G1, E.PALETA[2], "^", 4.5, None),
              ("Mika, Galerkin 2", G2, E.PALETA[6], "v", 4.5, None),
              (r"NeoPZ 24x3x12, $5^3$ pontos", d["fino"]["lam"], E.PALETA[1], "s", 3.2, None),
              ("notebook (guardado)", d["guardado"], E.PALETA[7], "D", 3.2, None)]
    for nome, v, cor, mk, ms, mfc in series:
        a1.plot(k, v, marker=mk, ms=ms, ls="none", color=cor, mfc=mfc if mfc else cor, mew=0.9, label=nome)
    a1.set_yscale("log")
    from matplotlib.ticker import FixedLocator, NullLocator
    a1.yaxis.set_major_locator(FixedLocator([5, 10, 20, 50, 100, 200]))
    a1.yaxis.set_major_formatter(E.fmt_virgula())
    a1.yaxis.set_minor_locator(NullLocator())
    a1.set_ylim(4.5, 220)
    a1.set_xticks([1, 5, 10, 15, 20])
    a1.set_xlabel(r"modo $k$")
    a1.set_ylabel(r"$\lambda_k$")
    a1.legend(loc="upper right", fontsize=6.5)
    for nome, v, cor, mk, ms, mfc in series[1:]:
        a2.plot(k, np.abs(erros(v, ref)), marker=mk, ms=ms, color=cor, lw=0.9, label=nome)
    a2.set_yscale("log")
    E.eixo_log(a2, "y")
    a2.set_xticks([1, 5, 10, 15, 20])
    a2.set_xlabel(r"modo $k$")
    a2.set_ylabel(r"$|\lambda_k^h - \lambda_k|/\lambda_k$ (contra Nyström)")
    a2.set_ylim(3e-7, 2)
    fig.tight_layout(w_pad=1.5)
    E.salva(fig, "kl_cil_espectro")


def fig_geometria():
    """meio cilindro 8 <= rho <= 10, 0 <= theta <= pi, 0 <= z <= 15 e as arestas da malha 12x2x6 na superficie"""
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
    nt, nr, nz = 12, 2, 6
    th, rr, zz = np.linspace(0, np.pi, nt + 1), np.linspace(8, 10, nr + 1), np.linspace(0, 15, nz + 1)
    tf = np.linspace(0, np.pi, 121)
    fig = plt.figure(figsize=(E.LARGURA * 0.62, 3.1))
    ax = fig.add_subplot(111, projection="3d")
    ax.grid(False)
    # superficies (externa, interna, topo) com transparencia
    T, Z = np.meshgrid(tf, [0, 15])
    for r, cor, alfa in ((10, E.PALETA[0], 0.16), (8, E.PALETA[0], 0.08)):
        ax.plot_surface(r * np.cos(T), r * np.sin(T), Z, color=cor, alpha=alfa, linewidth=0, shade=False)
    T2, R2 = np.meshgrid(tf, [8, 10])
    ax.plot_surface(R2 * np.cos(T2), R2 * np.sin(T2), 15 + 0 * T2, color=E.PALETA[0], alpha=0.12, linewidth=0,
                    shade=False)
    opt = dict(color=E.TINTA2, lw=0.45)
    for z in zz:  # arcos nas superficies rho = 8 e rho = 10
        for r in (8, 10):
            ax.plot(r * np.cos(tf), r * np.sin(tf), z + 0 * tf, **opt)
    for z in (0, 15):  # arcos e raios no fundo e no topo
        for r in rr:
            ax.plot(r * np.cos(tf), r * np.sin(tf), z + 0 * tf, **opt)
        for t in th:
            ax.plot([8 * np.cos(t), 10 * np.cos(t)], [8 * np.sin(t), 10 * np.sin(t)], [z, z], **opt)
    for t in th:  # verticais nas superficies rho = 8 e 10
        for r in (8, 10):
            ax.plot([r * np.cos(t)] * 2, [r * np.sin(t)] * 2, [0, 15], **opt)
    for t in (0, np.pi):  # faces planas theta = 0 e pi
        for r in rr:
            ax.plot([r * np.cos(t)] * 2, [0, 0], [0, 15], **opt)
        for z in zz:
            ax.plot([8 * np.cos(t), 10 * np.cos(t)], [0, 0], [z, z], **opt)
    ax.set_box_aspect((20, 10, 15))
    ax.view_init(elev=22, azim=-62)
    ax.set_xlabel(r"$x$ (m)", labelpad=-2)
    ax.set_ylabel(r"$y$ (m)", labelpad=-2)
    ax.set_zlabel(r"$z$ (m)", labelpad=-4)
    ax.set_xticks([-10, -5, 0, 5, 10])
    ax.set_yticks([0, 5, 10])
    ax.set_zticks([0, 5, 10, 15])
    ax.tick_params(pad=-2, labelsize=6.5)
    for a in (ax.xaxis, ax.yaxis, ax.zaxis):
        a.set_major_formatter(E.fmt_virgula())
        a.pane.fill = False
        a.pane.set_edgecolor(E.GRADE)
    E.salva(fig, "kl_cil_geometria")


def tex_num(v, casas):
    return ("%.*f" % (casas, v)).replace(".", ",")


def tab_autovalores(d):
    ref, fino = d["ref"], d["fino"]["lam"]
    L = [r"\begin{tabular}{r rrrr rrr}", r"\toprule",
         r" & \multicolumn{4}{c}{$\lambda_k$} & \multicolumn{3}{c}{$(\lambda_k^h - \lambda_k)/\lambda_k$} \\",
         r"\cmidrule(lr){2-5}\cmidrule(lr){6-8}",
         r"$k$ & Galerkin 1 & Galerkin 2 & NeoPZ & Nyström & Galerkin 1 & Galerkin 2 & NeoPZ \\", r"\midrule"]
    for k in range(20):
        L.append(r"%d & %s & %s & %s & %s & %s & %s & %s \\" % (
            k + 1, tex_num(G1[k], 6), tex_num(G2[k], 6), tex_num(fino[k], 6), tex_num(ref[k], 6),
            E.cientifico_tex(erros(G1, ref)[k]), E.cientifico_tex(erros(G2, ref)[k]),
            E.cientifico_tex(erros(fino, ref)[k])))
    L += [r"\midrule", r"g.l. & %d & %d & %d & -- & & & \\" % (GL_G1, GL_G2, d["fino"]["equacoes"]),
          r"\bottomrule", r"\end{tabular}"]
    open(os.path.join(E.TABELAS, "kl_cil_autovalores.tex"), "w").write("\n".join(L) + "\n")
    print("  tabelas/kl_cil_autovalores.tex")


def tab_convergencia(d):
    ref = d["ref"]
    L = [r"\begin{tabular}{c r r c c c c r r r}", r"\toprule",
         r"malha & elem. & $N$ & pontos & máx. erro & erro & $\sum\lambda/(270\pi)$ & $C$ (s) & autoval. (s) & "
         r"memória (MB) \\",
         r"$n_\theta\times n_\rho\times n_z$ & & & & $k \le 20$ & $\lambda_1$ & & & & \\", r"\midrule"]
    rs = sorted(d["cil"], key=lambda r: (r["ng"], r["equacoes"]))
    for i, r in enumerate(rs):
        if i == 4:
            L.append(r"\midrule")
        e = erros(r["lam"], ref)
        L.append(r"%s & %d & %d & $%d^3$ & %s & %s & %s & %s & %s & %s \\" % (
            r["malha"].replace("x", r"$\times$"), r["elementos"], r["equacoes"], r["ng"],
            E.cientifico_tex(np.max(np.abs(e))), E.cientifico_tex(e[0]), tex_num(r["soma"] / VOLUME, 4),
            tex_num(r["t_C"], 2), tex_num(r["t_aut"], 2), tex_num(r["memoria"], 0)))
    L += [r"\bottomrule", r"\end{tabular}"]
    open(os.path.join(E.TABELAS, "kl_cil_convergencia.tex"), "w").write("\n".join(L) + "\n")
    print("  tabelas/kl_cil_convergencia.tex")


def tab_notebook(d):
    ex = d["ref_reta"]
    nb3, nb5 = d["nb3"]["lam"], d["nb5"]["lam"]
    L = [r"\begin{tabular}{r r rrr r r rr}", r"\toprule",
         r" & Nyström & \multicolumn{3}{c}{Q2 (27 nós)} & serendipity & notebook & "
         r"\multicolumn{2}{c}{$(\lambda_k^h - \lambda_k)/\lambda_k$} \\",
         r"\cmidrule(lr){3-5}\cmidrule(lr){8-9}",
         r"$k$ & (semi-anel) & NeoPZ $3^3$ & NeoPZ $5^3$ & numpy $5^3$ & 20 nós, $5^3$ & (guardado) & NeoPZ $5^3$ & "
         r"notebook \\", r"\midrule"]
    for k in range(20):
        L.append(r"%d & %s & %s & %s & %s & %s & %s & %s & %s \\" % (
            k + 1, tex_num(ex[k], 4), tex_num(nb3[k], 4), tex_num(nb5[k], 4), tex_num(d["py_q2_5"][k], 4),
            tex_num(d["py_s20_5"][k], 4), tex_num(d["guardado"][k], 4), E.cientifico_tex(erros(nb5, ex)[k]),
            E.cientifico_tex(erros(d["guardado"], ex)[k])))
    L += [r"\midrule",
          r"$N$ & -- & %d & %d & %d & 408 & -- & & \\" % (d["nb3"]["equacoes"], d["nb5"]["equacoes"],
                                                         d["nb5"]["equacoes"]),
          r"\bottomrule", r"\end{tabular}"]
    open(os.path.join(E.TABELAS, "kl_cil_notebook.tex"), "w").write("\n".join(L) + "\n")
    print("  tabelas/kl_cil_notebook.tex")


def resumo(d):
    ref = d["ref"]
    print("erros contra Nystrom (lambda1, max 1-20) e taxas em N:")
    for ng in (3, 5):
        rs = sorted([r for r in d["cil"] if r["ng"] == ng], key=lambda r: r["equacoes"])
        ant = None
        for r in rs:
            e = erros(r["lam"], ref)
            s = "  %-8s %d^3 N=%5d e1=%.2e emax=%.2e min=%.2e soma/270pi=%.4f" % (
                r["malha"], ng, r["equacoes"], e[0], np.max(np.abs(e)), np.min(e), r["soma"] / VOLUME)
            if ant is not None:
                s += "  ordem em N (lambda1) %.2f, em h: queda %.2f" % (
                    np.log(ant[1] / e[0]) / np.log(r["equacoes"] / ant[0]), ant[1] / e[0])
            print(s)
            ant = (r["equacoes"], e[0])
    for nome, g in (("G1", G1), ("G2", G2)):
        e = erros(g, ref)
        print("  %s: e1=%.2e max=%.2e min=%.2e" % (nome, e[0], np.max(np.abs(e)), np.min(e)))
    print("  numpy q2 5^3 - NeoPZ 5^3 (notebook): max |dif| = %.2e ; 3^3: %.2e" % (
        np.max(np.abs(d["py_q2_5"] - d["nb5"]["lam"])), np.max(np.abs(d["py_q2_3"] - d["nb3"]["lam"]))))
    print("  s20 vs q2 (5^3): max rel %.2e" % np.max(np.abs(d["py_s20_5"] - d["nb5"]["lam"]) / d["nb5"]["lam"]))
    eg = erros(d["guardado"], d["ref_reta"])
    print("  notebook vs semi-anel: %.2e .. %.2e ; vs cilindro: %.2e .. %.2e" % (
        eg.min(), eg.max(), erros(d["guardado"], ref).min(), erros(d["guardado"], ref).max()))
    e5 = erros(d["nb5"]["lam"], d["ref_reta"])
    print("  NeoPZ 5^3 vs semi-anel: %.2e .. %.2e" % (e5.min(), e5.max()))


if __name__ == "__main__":
    d = gera_dados()
    resumo(d)
    tab_autovalores(d)
    tab_convergencia(d)
    tab_notebook(d)
    fig_convergencia(d)
    fig_espectro(d)
    fig_geometria()
