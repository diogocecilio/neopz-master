#!/usr/bin/env python3
"""sensibilidade_mc.py <pasta> [<pasta> ...] [--t <horas>] [--saida <prefixo>]: analise de sensibilidade do Monte Carlo
como na secao 6 de Vargas Ceron et al. (2025). Requer numpy e matplotlib.

Para cada pasta de combinacao (le mc_resumo.csv ou, com --t, mc_fs.csv; junta_mc.sh e chamado antes):
- correlacao de postos de Spearman entre o FS e as propriedades equivalentes de cada realizacao, como em Hu et al. e no
  Ceron: c' e phi' medios ao longo da superficie de ruptura (a banda: banda_c, banda_phi_deg) e k medio da massa acima
  dela (a massa que desliza: massa_k_mult; sem a coluna, de uma versao anterior, o k medio da banda). Tambem com as medias
  no solo inteiro (c_med, phi_med, k_mult_med), para comparar;
- de 1a ordem (rho_XY) e parcial de 2a ordem, controlando as outras duas propriedades (rho_cF,phik, rho_phiF,ck,
  rho_kF,cphi; eqs. 70-71 do Ceron), com o intervalo de 95% (Fisher). Vai para <pasta>/sensibilidade.csv e .png.
O FS e o FS minimo de cada realizacao (FS_min, com a banda e a massa do instante do minimo) ou, com --t 24, o FS em
t = 24 h (com a banda e a massa desse instante). Realizacoes sem mecanismo (FS = -1, sem equilibrio nem com c'/0.6) ficam
fora das correlacoes, mas contam como ruptura em P_f.

Com varias pastas (por exemplo os casos de CoV, de distancias de autocorrelacao ou de anisotropia, como nas tabelas 3 a 5
do Ceron), tambem: a tabela por caso (N, media, desvio e CoV do FS, P_f = P(FS < 1) e CoV(P_f) = sqrt((1 - P_f)/(N P_f)))
em <prefixo>_casos.csv e a densidade do FS de cada caso em <prefixo>_casos.png (prefixo padrao: sensibilidade, na pasta
corrente)."""
import os
import subprocess
import sys

import numpy as np

AZUL, LARANJA, AQUA, AMARELO, MAGENTA, VERDE, VIOLETA, VERMELHO = (
    "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948")
CINZA = "#8a8984"


def le_csv(arq):
    """colunas -> arrays (numeros; vazio = nan), lidas sem numpy.genfromtxt (cabecalhos com '.')"""
    with open(arq) as f:
        cab = f.readline().strip().split(",")
        linhas = [l.rstrip("\n").split(",") for l in f if l.strip()]
    col = {}
    for j, nome in enumerate(cab):
        v = [l[j] if j < len(l) else "" for l in linhas]
        try:
            col[nome] = np.array([float(x) if x != "" else np.nan for x in v])
        except ValueError:
            col[nome] = np.array(v, dtype=object)
    return col


def postos(x):
    """postos de 1 a n, com a media nos empates"""
    o = np.argsort(x, kind="mergesort")
    r = np.empty(len(x))
    xs = x[o]
    i = 0
    while i < len(x):
        j = i
        while j + 1 < len(x) and xs[j + 1] == xs[i]:
            j += 1
        r[o[i:j + 1]] = 0.5 * (i + j) + 1.
        i = j + 1
    return r


def spearman(dados):
    """matriz de correlacao de postos (Spearman) das colunas de dados (n x m)"""
    R = np.column_stack([postos(dados[:, k]) for k in range(dados.shape[1])])
    return np.corrcoef(R, rowvar=False)


def parcial(C, i, j):
    """correlacao parcial de i e j controlando todas as outras variaveis de C: -P_ij / sqrt(P_ii P_jj), P = inv(C).
    Com 4 variaveis (FS, c, phi, k), e a de 2a ordem das eqs. 70-71 do Ceron"""
    P = np.linalg.inv(C)
    return -P[i, j] / np.sqrt(P[i, i] * P[j, j])


def ic95(r, n, k):
    """intervalo de 95% de uma correlacao (parcial com k controles) pela transformacao de Fisher"""
    if n - 3 - k <= 0 or not np.isfinite(r) or abs(r) >= 1:
        return (np.nan, np.nan)
    z, s = np.arctanh(r), 1. / np.sqrt(n - 3 - k)
    return (np.tanh(z - 1.96 * s), np.tanh(z + 1.96 * s))


def amostra(D, t):
    """(FS, colunas das propriedades, numero de realizacoes) de uma pasta"""
    junta = os.path.join(os.path.dirname(os.path.abspath(__file__)), "junta_mc.sh")
    if os.path.exists(junta):
        subprocess.run([junta, D], check=False, stdout=subprocess.DEVNULL)
    if t is None:
        c = le_csv(os.path.join(D, "mc_resumo.csv"))
        fs = c["FS_min"]
    else:
        c = le_csv(os.path.join(D, "mc_fs.csv"))
        sel = np.isclose(c["t_h"], t)
        c = {k: v[sel] for k, v in c.items()}
        fs = c["FS_SRM"]
    return fs, c


GRUPOS = [  # (rotulo, colunas c, phi, k): a definicao do Ceron e as medias no solo inteiro
    ("superficie de ruptura e massa acima dela", ("banda_c", "banda_phi_deg", "massa_k_mult")),
    ("medias no solo inteiro", ("c_med", "phi_med", "k_mult_med")),
]
NOMES = ("c'", "phi'", "k_v")


def analisa(D, t, saida_png):
    fs, c = amostra(D, t)
    linhas = []
    resultados = []
    for rotulo, cols in GRUPOS:
        cols = list(cols)
        if cols[2] not in c or np.all(np.isnan(np.asarray(c.get(cols[2], [np.nan]), dtype=float))):
            if cols[2] == "massa_k_mult" and "banda_k_mult" in c:
                cols[2] = "banda_k_mult"  # versao anterior, sem a massa: o k da banda
            else:
                continue
        if not all(k in c for k in cols):
            continue
        X = np.column_stack([fs] + [np.asarray(c[k], dtype=float) for k in cols])
        ok = np.all(np.isfinite(X), axis=1) & (fs > 0)
        X = X[ok]
        n = len(X)
        if n < 5:
            print(f"{os.path.basename(os.path.normpath(D))}: {rotulo}: so {n} realizacoes com mecanismo; nada calculado")
            continue
        C = spearman(X)
        for k in range(3):
            r1 = C[0, k + 1]
            r2 = parcial(C, 0, k + 1)
            a1, b1 = ic95(r1, n, 0)
            a2, b2 = ic95(r2, n, 2)
            linhas.append((rotulo, NOMES[k], cols[k], n, r1, a1, b1, r2, a2, b2))
        resultados.append((rotulo, [parcial(C, 0, k + 1) for k in range(3)], [ic95(parcial(C, 0, k + 1), n, 2) for k in range(3)], n))
    nome = os.path.basename(os.path.normpath(D))
    print(f"== {nome}: FS = {'FS_min' if t is None else f'FS em t = {t:g} h'}; {len(fs)} realizacoes, "
          f"{int(np.sum(fs <= 0))} sem mecanismo (FS = -1)")
    print(f"{'propriedades':42s} {'var':5s} {'N':>5s} {'rho_XY':>8s} {'IC 95%':>16s} {'rho 2a ord.':>11s} {'IC 95%':>16s}")
    for (rot, nm, col, n, r1, a1, b1, r2, a2, b2) in linhas:
        print(f"{rot:42s} {nm:5s} {n:5d} {r1:8.3f} [{a1:6.3f}, {b1:6.3f}] {r2:11.3f} [{a2:6.3f}, {b2:6.3f}]")
    with open(os.path.join(D, "sensibilidade.csv"), "w") as f:
        f.write("propriedades,variavel,coluna,N,rho_spearman,ic95_inf,ic95_sup,rho_parcial_2a,ic95_parcial_inf,ic95_parcial_sup\n")
        for l in linhas:
            f.write(",".join(str(x) for x in l) + "\n")
    if resultados:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(7, 3.6), constrained_layout=True)
        y = np.arange(3)
        h = 0.8 / len(resultados)
        for g, ((rot, rr, ics, n), cor) in enumerate(zip(resultados, (AZUL, LARANJA))):
            pos = y + (g - (len(resultados) - 1) / 2) * h
            err = np.array([[r - a, b - r] for r, (a, b) in zip(rr, ics)]).T
            ax.barh(pos, rr, height=h * 0.9, color=cor, label=f"{rot} (N = {n})", xerr=err, ecolor=CINZA, capsize=3)
        ax.axvline(0, color="#0b0b0b", lw=0.8)
        ax.set_yticks(y, NOMES)
        ax.set_xlim(-1, 1)
        ax.set_xlabel("correlacao parcial de postos (Spearman, 2a ordem) com o FS")
        ax.grid(axis="x", color="0.9", lw=0.6)
        ax.set_axisbelow(True)
        ax.legend(fontsize=8, loc="lower left")
        ax.set_title(f"{nome}: {'FS minimo' if t is None else f'FS em t = {t:g} h'} (barras: IC 95%)", fontsize=10)
        fig.savefig(saida_png, dpi=110)
        print(saida_png)


def casos(pastas, t, prefixo):
    tab = []
    for D in pastas:
        fs, _ = amostra(D, t)
        fs = fs[np.isfinite(fs)]
        n = len(fs)
        if n == 0:
            continue
        pos = fs[fs > 0]
        pf = float(np.mean(fs < 1))
        tab.append((os.path.basename(os.path.normpath(D)), n, np.mean(pos), np.std(pos, ddof=1) if len(pos) > 1 else np.nan,
                    pf, np.sqrt((1 - pf) / (n * pf)) if pf > 0 else np.nan, pos))
    print(f"{'caso':32s} {'N':>5s} {'media':>8s} {'desvio':>8s} {'CoV':>7s} {'P_f':>8s} {'CoV(P_f)':>9s}")
    with open(prefixo + "_casos.csv", "w") as f:
        f.write("caso,N,media_FS,desvio_FS,CoV_FS,P_f,CoV_P_f\n")
        for (nm, n, mu, dp, pf, cv, _) in tab:
            print(f"{nm:32s} {n:5d} {mu:8.4f} {dp:8.4f} {dp / mu:7.3f} {pf:8.4f} {cv:9.3f}")
            f.write(f"{nm},{n},{mu},{dp},{dp / mu},{pf},{cv}\n")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    cores = (AZUL, LARANJA, AQUA, AMARELO, MAGENTA, VERDE, VIOLETA, VERMELHO)
    fig, ax = plt.subplots(1, 2, figsize=(12, 4), constrained_layout=True)
    xs = np.linspace(min(np.min(p) for *_, p in tab) * 0.9, max(np.max(p) for *_, p in tab) * 1.1, 300)
    for k, (nm, n, mu, dp, pf, cv, p) in enumerate(tab[:8]):
        bw = 1.06 * np.std(p) * len(p) ** -0.2 if len(p) > 1 else 0.05  # Silverman
        dens = np.exp(-0.5 * ((xs[:, None] - p[None, :]) / bw) ** 2).sum(axis=1) / (len(p) * bw * np.sqrt(2 * np.pi))
        ax[0].plot(xs, dens, color=cores[k], lw=2, label=nm)
    if len(tab) > 8:
        print(f"(a densidade mostra so os 8 primeiros casos de {len(tab)})")
    ax[0].axvline(1, color=CINZA, ls="--", lw=1)
    ax[0].set_xlabel("FS minimo" if t is None else f"FS em t = {t:g} h")
    ax[0].set_ylabel("densidade de probabilidade")
    ax[0].legend(fontsize=8)
    i = np.arange(len(tab))
    pf = np.array([x[4] for x in tab])
    cv = np.array([x[5] for x in tab])
    ax[1].errorbar(i, pf, yerr=1.96 * cv * pf, fmt="o", color="#0b0b0b", ecolor=CINZA, capsize=3)
    ax[1].set_xticks(i, [x[0] for x in tab], rotation=30, ha="right", fontsize=8)
    ax[1].set_ylabel("P_f = P(FS < 1) (barras: 95%)")
    ax[1].set_ylim(bottom=0)
    for a in ax:
        a.grid(color="0.9", lw=0.6)
        a.set_axisbelow(True)
    fig.savefig(prefixo + "_casos.png", dpi=110)
    print(prefixo + "_casos.csv", prefixo + "_casos.png")


def main():
    args = sys.argv[1:]
    t, prefixo, pastas = None, "sensibilidade", []
    while args:
        a = args.pop(0)
        if a == "--t":
            t = float(args.pop(0))
        elif a == "--saida":
            prefixo = args.pop(0)
        else:
            pastas.append(a)
    if not pastas:
        print(__doc__)
        sys.exit(2)
    for D in pastas:
        analisa(D, t, os.path.join(D, "sensibilidade.png"))
    if len(pastas) > 1:
        casos(pastas, t, prefixo)


if __name__ == "__main__":
    main()
