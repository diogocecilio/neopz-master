#!/usr/bin/env python3
"""Figuras do VerifyRandomFields (rotulos em portugues, virgula decimal) a partir das saidas do programa.

Uso: python3 plota.py [SAIDA] [SAIDA_h05]
  SAIDA      diretorio com pontos.csv, modos.csv e varredura.csv da malha padrao (h = 0,25 m);
             padrao: saida_talude_uniforme_h025 no diretorio atual
  SAIDA_h05  (opcional) o mesmo para a malha h = 0,5 m: entra nas figuras da varredura e dos modos
As figuras (PNG) vao para figuras/, ao lado deste arquivo. Se existir, o espectro de referencia do talude (Nystrom
extrapolado, ../SlopeRandomFieldRefina/relatorio/dados/kl_espectros.csv) entra na figura dos autovalores.

Mapas: os valores nos pontos de integracao da malha de Galerkin (triangulacao de Delaunay dos pontos, sem os triangulos fora
do talude), com o contorno do talude. Cores: a paleta categorica fixa #2a78d6, #eb6834, #1baf7a, #eda100, #e87ba4, #008300,
#6250d6, #e34948 (nessa ordem); escala sequencial azul para a coesao; divergente azul - cinza - vermelho para os erros.
"""
import csv
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker as mticker  # noqa: E402
import matplotlib.tri as mtri  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm  # noqa: E402

AQUI = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(AQUI, "figuras")
REF_ESPECTRO = os.path.join(AQUI, "..", "SlopeRandomFieldRefina", "relatorio", "dados", "kl_espectros.csv")

# paleta categorica (ordem fixa) e tinta
AZUL, LARANJA, AQUA, AMARELO, MAGENTA, VERDE, VIOLETA, VERMELHO = (
    "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948")
TINTA, TINTA2, MUDO, GRADE, EIXO, FUNDO = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb"
# sequencial: rampa azul 100 -> 700; divergente: azul <- cinza -> vermelho (bracos com o mesmo numero de passos)
SEQ = LinearSegmentedColormap.from_list("azul", ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#2a78d6", "#1c5cab",
                                                 "#104281", "#0d366b"])
DIV = LinearSegmentedColormap.from_list("azul_vermelho", ["#0d366b", "#2a78d6", "#9ec5f4", "#f0efec", "#f3b3b2",
                                                          "#e34948", "#8c1c1c"])
ROTULO = {"padrao": "malha padrão", "segunda": "segunda malha"}  # trocados pelo h das malhas (resumo.txt)
PARAM = {"media": 10.0, "cov": 0.3, "Lx": 20.0, "Ly": 2.0, "real": 0}  # lidos de resumo.txt (os padroes do programa)
TALUDE = np.array([[0, 0], [25, 0], [25, 5], [15, 5], [10, 10], [0, 10], [0, 0]], float)

plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
    "text.color": TINTA, "axes.labelcolor": TINTA2, "axes.edgecolor": EIXO, "axes.linewidth": 0.8,
    "xtick.color": TINTA2, "ytick.color": TINTA2, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "grid.color": GRADE, "grid.linewidth": 0.6, "grid.linestyle": "-",
    "figure.facecolor": FUNDO, "axes.facecolor": FUNDO, "savefig.facecolor": FUNDO,
    "legend.frameon": False, "legend.fontsize": 8, "axes.spines.top": False, "axes.spines.right": False,
    "axes.unicode_minus": True,
})


def virgula(s):
    return s.replace(".", ",")


def num(v, nd=2):
    return virgula(f"{v:.{nd}f}").replace("-", "\u2212")


def sig(v, n=3):
    """v com n algarismos significativos e virgula decimal (0,00549; 0,591; 13,0)"""
    s = f"{v:#.{n}g}".rstrip(".")
    if "e" in s:
        s = f"{v:.{n - 1}e}"
    return virgula(s).replace("-", "\u2212")


class FormatoVirgula(mticker.ScalarFormatter):
    """ScalarFormatter com virgula decimal (rotulos e deslocamento)."""

    def __call__(self, x, pos=None):
        return virgula(super().__call__(x, pos))

    def get_offset(self):
        return virgula(super().get_offset())


def formata_eixos(ax, x=True, y=True):
    if x and ax.get_xscale() == "linear":
        ax.xaxis.set_major_formatter(FormatoVirgula())
    if y and ax.get_yscale() == "linear":
        ax.yaxis.set_major_formatter(FormatoVirgula())


def le_csv(caminho):
    with open(caminho) as f:
        r = csv.reader(f)
        cab = next(r)
        dados = np.array([[float(v) for v in linha] for linha in r])
    return {c: dados[:, i] for i, c in enumerate(cab)}


def triangulacao(x, y):
    t = mtri.Triangulation(x, y)
    xc, yc = x[t.triangles].mean(axis=1), y[t.triangles].mean(axis=1)
    fora = (yc > 5.0) & (xc + yc > 20.0)  # acima do terreno do pe e a direita da face
    # triangulos degenerados ao longo do contorno (pontos de integracao nao chegam ao contorno)
    lado = np.max(np.hypot(np.diff(x[t.triangles][:, [0, 1, 2, 0]], axis=1), np.diff(y[t.triangles][:, [0, 1, 2, 0]], axis=1)),
                  axis=1)
    t.set_mask(fora | (lado > 1.5))
    return t


def mapa(ax, tri, v, cmap, norm, titulo, ref=None):
    im = ax.tripcolor(tri, v, cmap=cmap, norm=norm, shading="gouraud", rasterized=True)
    ax.plot(TALUDE[:, 0], TALUDE[:, 1], color=TINTA2, lw=0.8)
    if ref is not None:
        ax.plot(ref[0], ref[1], marker="o", ms=6, mfc="white", mec=TINTA, mew=1.2, ls="none")
    ax.set_aspect("equal")
    ax.set_xlim(-0.3, 25.3)
    ax.set_ylim(-0.3, 10.3)
    ax.set_title(titulo, loc="left")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    formata_eixos(ax)
    return im


def barra(fig, im, ax, rotulo, aspecto=18):
    cb = fig.colorbar(im, ax=ax, shrink=0.92, pad=0.02, aspect=aspecto)
    cb.set_label(rotulo)
    cb.ax.yaxis.set_major_formatter(FormatoVirgula())
    cb.outline.set_edgecolor(EIXO)
    cb.outline.set_linewidth(0.6)
    return cb


def norma_simetrica(*vs, q=99.9):
    lim = max(np.percentile(np.abs(v), q) for v in vs)
    return TwoSlopeNorm(vcenter=0.0, vmin=-lim, vmax=lim)


def salva(fig, nome):
    os.makedirs(FIG, exist_ok=True)
    caminho = os.path.join(FIG, nome)
    fig.savefig(caminho, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print("figura:", caminho)


def rms(w, v):
    return float(np.sqrt(np.sum(w * v * v) / np.sum(w)))


def figuras_mapas(P, M):
    tri = triangulacao(P["x"], P["y"])
    w = P["peso"]
    # 1. coesao nas tres formas, escala comum
    campos = [("cG", "Galerkin (G), M = %d" % M), ("cN", "Nyström (N), M = %d" % M),
              ("cA", "forma fechada do retângulo (A), M = %d" % M), ("cA2000", "forma fechada, 2000 termos (produção)")]
    vmin = min(np.percentile(P[c], 0.1) for c, _ in campos)
    vmax = max(np.percentile(P[c], 99.9) for c, _ in campos)
    fig, axs = plt.subplots(2, 2, figsize=(11.5, 5.6), constrained_layout=True)
    norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    for ax, (c, t) in zip(axs.flat, campos):
        im = mapa(ax, tri, P[c], SEQ, norm, t)
    barra(fig, im, axs[:, 1], "coesão c (kPa)", aspecto=30)
    fig.suptitle(r"Realização REAL = %d com os mesmos $\xi_k$ nos três métodos (média %s kPa, CoV %s; "
                 r"$L_x$ = %s m, $L_y$ = %s m; %s)" % (PARAM["real"], virgula("%g" % PARAM["media"]),
                                                     virgula("%g" % PARAM["cov"]), virgula("%g" % PARAM["Lx"]),
                                                     virgula("%g" % PARAM["Ly"]), ROTULO["padrao"]),
                 x=0.0, ha="left", fontsize=10)
    salva(fig, "campos_coesao.png")

    # 2. diferencas na realizacao (escalas proprias: G - N e muito menor que G - A). G - N com M_VERIF termos (colunas
    # dGN_M<m>, modos que as duas discretizacoes resolvem: a verificacao ponto a ponto) e com M termos (dominada pelos modos
    # altos, que nenhuma das duas resolve); G - A e N - A: bases diferentes
    mv = sorted((int(c[5:]), c) for c in P if c.startswith("dGN_M"))
    difs = [(c, r"$c_G - c_N$, M = %d: modos resolvidos por G e N (verificação)" % m) for m, c in mv]
    difs += [("dGN", r"$c_G - c_N$, M = %d: diferença entre as discretizações G e N" % M),
             ("dGA", r"$c_G - c_A$, M = %d (bases diferentes: outra realização, não é erro)" % M),
             ("dNA", r"$c_N - c_A$, M = %d (bases diferentes: outra realização, não é erro)" % M)]
    ncol = 2 if len(difs) % 2 == 0 else 1
    nlin = len(difs) // ncol
    fig, axs = plt.subplots(nlin, ncol, figsize=(7.2 * ncol, 3.0 * nlin), constrained_layout=True, squeeze=False)
    for ax, (c, t) in zip(axs.flat, difs):
        im = mapa(ax, tri, P[c], DIV, norma_simetrica(P[c]), t + "\nRMS " + sig(rms(w, P[c])) + " kPa, máx. |.| "
                  + sig(np.abs(P[c]).max()) + " kPa")
        barra(fig, im, ax, "kPa")
    salva(fig, "erros_realizacao.png")

    # 3. deficit da variancia truncada, 1 - sigma_M^2 (alvo 0), escala comum: sequencial se nao ha sigma_M^2 > 1
    fig, axs = plt.subplots(3, 1, figsize=(7.2, 9.0), constrained_layout=True)
    ev = {m: 1.0 - P["s2" + m] for m in "GNA"}
    nomes = {"G": "Galerkin", "N": "Nyström", "A": "forma fechada do retângulo"}
    if min(v.min() for v in ev.values()) >= 0.0:  # escala comum nos percentis 0,5 e 99,5 dos tres campos
        cmap, norm = SEQ, matplotlib.colors.Normalize(vmin=min(np.percentile(v, 0.5) for v in ev.values()),
                                                      vmax=max(np.percentile(v, 99.5) for v in ev.values()))
    else:
        cmap, norm = DIV, norma_simetrica(*ev.values(), q=100)
    for ax, m in zip(axs, "GNA"):
        im = mapa(ax, tri, ev[m], cmap, norm, r"$1 - \sigma_M^2$, %s (M = %d): média %s, mín. %s, máx. %s" % (
            nomes[m], M, num(np.sum(w * ev[m]) / np.sum(w), 4), num(ev[m].min(), 3), num(ev[m].max(), 3)))
        barra(fig, im, ax, r"$1 - \sigma_M^2$")
    salva(fig, "erro_variancia.png")

    # 4. Monte Carlo: media - 10 e desvio - 3, escala comum por linha
    if "mcmG" in P and np.any(P["mcmG"] != 0.0):  # R = 0: sem Monte Carlo (colunas zeradas)
        fig, axs = plt.subplots(3, 2, figsize=(13.0, 8.4), constrained_layout=True)
        media, desvio = PARAM["media"], PARAM["media"] * PARAM["cov"]
        em = {m: P["mcm" + m] - media for m in "GNA"}
        ed = {m: P["mcd" + m] - desvio for m in "GNA"}
        nm, nd = norma_simetrica(*em.values(), q=100), norma_simetrica(*ed.values(), q=100)
        for i, m in enumerate("GNA"):
            im1 = mapa(axs[i, 0], tri, em[m], DIV, nm, "%s: média de c − %s kPa (RMS %s kPa)" % (
                nomes[m], virgula("%g" % media), num(rms(w, em[m]), 3)))
            im2 = mapa(axs[i, 1], tri, ed[m], DIV, nd, "%s: desvio de c − %s kPa (RMS %s kPa)" % (
                nomes[m], virgula("%g" % desvio), num(rms(w, ed[m]), 3)))
        barra(fig, im1, axs[:, 0], "kPa", aspecto=40)
        barra(fig, im2, axs[:, 1], "kPa", aspecto=40)
        salva(fig, "monte_carlo.png")

    # 5. covariancia com os pontos de referencia: erro contra o nucleo
    refs = sorted({int(c[3:c.index("_")]) for c in P if c.startswith("cov") and c.endswith("_alvo")})
    if refs:
        fig, axs = plt.subplots(len(refs), 3, figsize=(15.0, 2.45 * len(refs) + 0.3), constrained_layout=True)
        axs = np.atleast_2d(axs)
        ec = {(r, m): P["cov%d_%s" % (r, m)] - P["cov%d_alvo" % r] for r in refs for m in "GNA"}
        norm = norma_simetrica(*ec.values(), q=99.9)
        for i, r in enumerate(refs):
            alvo = P["cov%d_alvo" % r]
            k = int(np.argmax(alvo))
            for j, m in enumerate("GNA"):
                im = mapa(axs[i, j], tri, ec[(r, m)], DIV, norm, "%s, ponto %d: RMS %s" % (
                    nomes[m], r, virgula("%.1e" % rms(w, ec[(r, m)]))), ref=(P["x"][k], P["y"][k]))
        barra(fig, im, axs[:, 2], r"$\sum_{k<M} \lambda_k f_k(x) f_k(x_r) - C(x, x_r)$", aspecto=40)
        salva(fig, "covariancia.png")


def figura_varredura(V, V05):
    fig, axs = plt.subplots(1, 3, figsize=(13.5, 3.9), constrained_layout=True)
    M = V["M"]
    series = [("dGN", r"$c_G - c_N$", AZUL, "o"), ("dGA", r"$c_G - c_A$", LARANJA, "s"), ("dNA", r"$c_N - c_A$", AQUA, "^")]
    ax = axs[0]
    for c, rot, cor, mk in series:
        ax.plot(M, V["rms_" + c], color=cor, marker=mk, ms=5, lw=1.6, label=rot + " (%s)" % ROTULO["padrao"])
    if V05 is not None:
        ax.plot(V05["M"], V05["rms_dGN"], color=VIOLETA, marker="D", ms=5, lw=1.6,
                label=r"$c_G - c_N$ (%s)" % ROTULO["segunda"])
    ax.plot(M, V["desvio_cG"], color=MUDO, lw=1.0, label=r"desvio de $c_G$ na área")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("M (termos)")
    ax.set_ylabel("RMS na área (kPa)")
    ax.set_title("(a) diferenças na realização REAL = 0", loc="left")
    ax.grid(True, which="major")
    ax.legend(loc="lower right")
    ax = axs[1]
    for c, rot, cor, mk in [("GN", r"$Z_G - Z_N$", AZUL, "o"), ("GA", r"$Z_G - Z_A$", LARANJA, "s"),
                            ("NA", r"$Z_N - Z_A$", AQUA, "^")]:
        ax.plot(M, V["ez_" + c], color=cor, marker=mk, ms=5, lw=1.6, label=rot + " (%s)" % ROTULO["padrao"])
    if V05 is not None:
        ax.plot(V05["M"], V05["ez_GN"], color=VIOLETA, marker="D", ms=5, lw=1.6,
                label=r"$Z_G - Z_N$ (%s)" % ROTULO["segunda"])
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("M (termos)")
    ax.set_ylabel(r"RMS na área de $\sqrt{E[(Z_1 - Z_2)^2]}$")
    ax.set_title("(b) diferença esperada (todas as realizações)", loc="left")
    ax.grid(True, which="major")
    ax.legend(loc="lower right")
    ax = axs[2]
    for c, rot, cor, mk in [("G", "Galerkin", AZUL, "o"), ("N", "Nyström", LARANJA, "s"), ("A", "forma fechada", AQUA, "^")]:
        ax.plot(M, V["erro_var_" + c], color=cor, marker=mk, ms=5, lw=1.6, label=rot)
    if "erro_var_lambda_N" in V:  # o mesmo pelos autovalores (em G, igual ao dos pontos: modos de norma 1)
        ax.plot(M, V["erro_var_lambda_N"], color=LARANJA, marker="s", mfc=FUNDO, ms=5, lw=1.2, ls="--",
                label=r"Nyström, 1 $-$ $\Sigma\,\lambda_k$ / 187,5")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("M (termos)")
    ax.set_ylabel(r"1 $-$ média de $\sigma_M^2$ no talude")
    ax.set_title("(c) erro médio de truncamento", loc="left")
    ax.grid(True, which="major")
    ax.legend(loc="lower left")
    salva(fig, "varredura.png")


def figura_modos(Mo, Mo05, ref):
    fig, axs = plt.subplots(1, 3, figsize=(13.5, 3.9), constrained_layout=True)
    k = Mo["k"]
    ax = axs[0]
    ax.plot(k, Mo["lambda_G"], color=AZUL, lw=1.6, label="Galerkin (talude)")
    ax.plot(k, Mo["lambda_N"], color=LARANJA, lw=1.6, label="Nyström (talude)")
    ax.plot(k, Mo["lambda_A"], color=AQUA, lw=1.6, label="forma fechada (retângulo)")
    if ref is not None:
        ax.plot(np.arange(1, len(ref) + 1), ref, ls="none", marker="o", ms=4, mfc="none", mec=TINTA2, mew=0.8,
                label="referência do talude (Nyström extrapolado)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("k")
    ax.set_ylabel(r"$\lambda_k$ (m²)")
    ax.set_title("(a) autovalores", loc="left")
    ax.grid(True, which="major")
    ax.legend(loc="lower left")
    ax = axs[1]
    if ref is not None:
        n = min(len(ref), 60, len(k))
        ax.plot(k[:n], np.abs(Mo["lambda_G"][:n] / ref[:n] - 1), color=AZUL, marker="o", ms=4, lw=1.4,
                label=r"Galerkin, %s" % ROTULO["padrao"])
        ax.plot(k[:n], np.abs(Mo["lambda_N"][:n] / ref[:n] - 1), color=LARANJA, marker="s", ms=4, lw=1.4,
                label="Nyström")
        if Mo05 is not None:
            ax.plot(k[:n], np.abs(Mo05["lambda_G"][:n] / ref[:n] - 1), color=VIOLETA, marker="D", ms=4, lw=1.4,
                    label=r"Galerkin, %s" % ROTULO["segunda"])
        ax.set_title(r"(b) $|\lambda_k / \lambda_k^{ref} - 1|$ (k ≤ %d)" % n, loc="left")
    else:
        n = min(60, len(k))
        ax.plot(k[:n], np.abs(Mo["dif_rel_NG"][:n]), color=AZUL, marker="o", ms=4, lw=1.4,
                label=r"$|\lambda_N / \lambda_G - 1|$")
        ax.set_title("(b) autovalores: Nyström contra Galerkin", loc="left")
    ax.set_yscale("log")
    ax.set_xlabel("k")
    ax.set_ylabel("diferença relativa")
    ax.grid(True, which="major")
    ax.legend(loc="lower right")
    formata_eixos(ax, y=False)
    ax = axs[2]
    ax.plot(k, np.maximum(1.0 - Mo["acordo"], 1e-12), color=AZUL, lw=1.4, label=ROTULO["padrao"])
    if Mo05 is not None:
        ax.plot(Mo05["k"], np.maximum(1.0 - Mo05["acordo"], 1e-12), color=VIOLETA, lw=1.4, label=ROTULO["segunda"])
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("k")
    ax.set_ylabel(r"$1 - |\cos(f_k^G, f_k^N)|$")
    ax.set_title("(c) concordância dos modos de G e N", loc="left")
    ax.grid(True, which="major")
    ax.legend(loc="lower right")
    salva(fig, "modos.png")


def rotulo_malha(saida, padrao):
    """'h = 0,25 m' a partir do nome da malha em resumo.txt (talude_uniforme_h025.msh), ou o padrao."""
    try:
        for linha in open(os.path.join(saida, "resumo.txt")):
            if linha.strip().startswith("malha "):
                nome = os.path.basename(linha.split()[-1])
                if "_h0" in nome:
                    d = nome.split("_h0")[1].split(".")[0]
                    return "h = 0,%s m" % d
                return nome
    except OSError:
        pass
    return padrao


def main():
    saida = sys.argv[1] if len(sys.argv) > 1 else "saida_talude_uniforme_h025"
    saida05 = sys.argv[2] if len(sys.argv) > 2 else None
    ROTULO["padrao"] = rotulo_malha(saida, ROTULO["padrao"])
    if saida05:
        ROTULO["segunda"] = rotulo_malha(saida05, ROTULO["segunda"])
    P = le_csv(os.path.join(saida, "pontos.csv"))
    V = le_csv(os.path.join(saida, "varredura.csv"))
    Mo = le_csv(os.path.join(saida, "modos.csv"))
    V05 = le_csv(os.path.join(saida05, "varredura.csv")) if saida05 else None
    Mo05 = le_csv(os.path.join(saida05, "modos.csv")) if saida05 else None
    ref = None
    if os.path.exists(REF_ESPECTRO):
        with open(REF_ESPECTRO) as f:
            ref = np.array([float(r["nystrom_talude"]) for r in csv.DictReader(f) if r["nystrom_talude"] not in ("", "nan")])
    M = 200
    resumo = os.path.join(saida, "resumo.txt")
    if os.path.exists(resumo):
        for linha in open(resumo):
            if linha.startswith("Realizacao REAL") and "M = " in linha:
                M = int(linha.split("M = ")[1].split()[0])
                PARAM["real"] = int(linha.split("REAL = ")[1].split()[0])
            if "nucleo C = " in linha and "media " in linha:  # Lx = 20 m, Ly = 2 m; coesao lognormal, media 10 kPa, CoV 0.3
                PARAM["Lx"] = float(linha.split("Lx = ")[1].split()[0])
                PARAM["Ly"] = float(linha.split("Ly = ")[1].split()[0].rstrip(";"))
                PARAM["media"] = float(linha.split("media ")[1].split()[0])
                PARAM["cov"] = float(linha.split("CoV ")[1].split()[0])
    figuras_mapas(P, M)
    figura_varredura(V, V05)
    figura_modos(Mo, Mo05, ref)


if __name__ == "__main__":
    main()
