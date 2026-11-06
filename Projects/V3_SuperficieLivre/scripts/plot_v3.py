#!/usr/bin/env python3
# plot_v3.py: compara as linhas freaticas do V3 com as digitalizadas de Sharma et al. (2021) e imprime a tabela
# "artigo x nosso". Uso: python3 plot_v3.py <pasta com uma subpasta por execucao> [--saida figura.png]
# Cada subpasta tem os CSVs do executavel (resumo.csv, linha_freatica.csv, face.csv). As execucoes de um mesmo caso
# (retangular, retangular_a2, ...) vao para o mesmo painel.
import argparse, csv, os, sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

AQUI = os.path.dirname(os.path.abspath(__file__))
DADOS = os.path.join(AQUI, "..", "dados")
GEOM = {  # contorno da barragem: A, B, D, E
    "retangular": [(0, 0), (16, 0), (16, 24), (0, 24), (0, 0)],
    "inclinada": [(0, 0), (7, 0), (2, 5), (0, 5), (0, 0)],
    "dois_blocos": [(0, 0), (5, 0), (5, 10), (0, 10), (0, 0)],
}
REF = {"retangular": "sharma2021_fig7_retangular.csv", "inclinada": "sharma2021_fig14_inclinada.csv",
       "dois_blocos": "sharma2021_fig17_dois_blocos.csv"}
SAIDA_ARTIGO = {"retangular": 12.9, "inclinada": 3.4, "dois_blocos": 2.46}
NIVEIS = {"retangular": (24, 4), "inclinada": (5, 1), "dois_blocos": (10, 2)}


def le_csv(arq):
    with open(arq) as f:
        linhas = [l for l in f if l.strip() and not l.startswith("#")]
    return list(csv.DictReader(linhas))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pasta")
    ap.add_argument("--saida", default="v3_comparacao.png")
    a = ap.parse_args()
    execs = sorted(d for d in os.listdir(a.pasta) if os.path.isfile(os.path.join(a.pasta, d, "resumo.csv")))
    if not execs:
        print("nenhuma execucao em", a.pasta)
        return 1
    casos = [c for c in GEOM if any(le_csv(os.path.join(a.pasta, e, "resumo.csv"))[0]["caso"] == c for e in execs)]
    fig, axs = plt.subplots(1, len(casos), figsize=(5.5 * len(casos), 6), squeeze=False)
    print("execucao | caso | malha | a (1/kPa) | saida nosso (m) | saida artigo (m) | q / q_DC - 1 | "
          "q / q_DC_corr - 1 | desvio medio / max da linha (m) | passos | s")
    for e in execs:
        r = le_csv(os.path.join(a.pasta, e, "resumo.csv"))[0]
        c = r["caso"]
        ax = axs[0][casos.index(c)]
        lf = le_csv(os.path.join(a.pasta, e, "linha_freatica.csv"))
        x = np.array([float(l["x_m"]) for l in lf])
        y = np.array([float(l["y_nosso_m"]) for l in lf])
        m = y >= 0
        ax.plot(x[m], y[m], lw=1.4, label=f"{e} (a = {float(r['aG']):g}/kPa, {r['nx']}x{r['ny']})")
        ys, ya = float(r["y_saida_m"]), float(r["y_acima_m"])
        qdc, qdcc, q = float(r["q_DC_m3s"]), float(r["q_DC_corr_m3s"]), float(r["q_esq_m3s"])
        rel = f"{q / qdc - 1:+.4f}" if qdc > 0 else "-"
        relc = f"{q / qdcc - 1:+.5f}" if qdcc > 0 else "-"
        print(f"{e} | {c} | {r['nx']}x{r['ny']} | {float(r['aG']):g} | {ys:.3f} a {ya:.3f} | "
              f"{float(r['y_saida_artigo_m']):.2f} | {rel} | {relc} | {float(r['desvio_linha_med_m']):.3f} / "
              f"{float(r['desvio_linha_max_m']):.3f} | {r['passos']} | {float(r['seg']):.0f}")
    for c in casos:
        ax = axs[0][casos.index(c)]
        g = np.array(GEOM[c])
        ax.plot(g[:, 0], g[:, 1], "k-", lw=0.8)
        h1, h2 = NIVEIS[c]
        if c in REF:
            pts = le_csv(os.path.join(DADOS, REF[c]))
            ax.plot([float(p["x_m"]) for p in pts], [float(p["y_m"]) for p in pts], "ko", ms=3.5, mfc="none",
                    label=f"Sharma et al. (2021), digitalizado")
        ax.axhline(h2, color="tab:blue", lw=0.6, ls=":")
        ax.plot([], [], " ", label=f"saida no artigo: y = {SAIDA_ARTIGO[c]} m")
        ax.set_aspect("equal")
        ax.set_xlabel("x (m)")
        ax.set_ylabel("y (m)")
        ax.set_title(c)
        ax.legend(fontsize=7, loc="lower left")
        ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(a.saida, dpi=130)
    print("figura:", a.saida)
    return 0


if __name__ == "__main__":
    sys.exit(main())
