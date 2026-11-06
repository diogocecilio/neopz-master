#!/usr/bin/env python3
"""sensibilidade_mc.py <pasta> [<pasta> ...] [--t <horas>] [--saida <prefixo>]: analise de sensibilidade do Monte Carlo
como na secao 6 de Vargas Ceron et al. (2025), para as pastas da fila do cluster ($SAIDA/<projeto>/<caso>[_trocas]).
Requer numpy e matplotlib.

Para cada pasta (le mc_resumo.csv e mc_fs.csv; junta_mc.sh e chamado antes):
- correlacao de postos de Spearman entre o FS e as propriedades equivalentes de cada realizacao que o programa grava,
  em grupos: c' e phi' medios na superficie de ruptura (a banda: banda_c, banda_phi_deg) com o k da massa acima dela
  (massa_k_mult, ou o da banda), como em Hu et al. e no Ceron; c', phi' e k medios da massa que desliza (massa_c,
  massa_phi_deg, massa_k_mult); e as medias no solo inteiro (c_med, phi_med, k_mult_med). Num grupo entram as
  propriedades que existem e variam (com c_u so aleatorio, so c');
- de 1a ordem (rho_XY) e parcial de ordem maior, controlando as outras propriedades do grupo (com tres, a de 2a ordem
  das eqs. 70-71 do Ceron), com o intervalo de 95% (Fisher). Vai para <pasta>/sensibilidade.csv e .png.
O FS e o FS minimo de cada realizacao (FS_min) ou, com --t 24 (so com varios instantes do FS por realizacao), o FS em
t = 24 h. Realizacoes sem equilibrio (FS <= 0, ou uma situacao em SITUACOES_RUPTURA, padrao sem_equilibrio_geostatico e
fs_abaixo_de_F0) ficam fora das correlacoes, mas contam como ruptura em P_f; as de outra situacao != ok (o calculo nao
terminou) ficam fora de tudo.

Com varias pastas (por exemplo os casos de CoV, de distancias de autocorrelacao ou de anisotropia, como nas tabelas 3 a 5
do Ceron), tambem: a tabela por caso (N, media, desvio e CoV do FS, P_f = P(FS < 1) e CoV(P_f) = sqrt((1 - P_f)/(N P_f)))
em <prefixo>_casos.csv e a densidade do FS de cada caso em <prefixo>_casos.png (prefixo padrao: sensibilidade, na pasta
corrente). LIMITE_FS (padrao 1) troca o limiar de P_f."""
import os
import subprocess
import sys

import numpy as np

AZUL, LARANJA, AQUA, AMARELO, MAGENTA, VERDE, VIOLETA, VERMELHO = (
    "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948")
CINZA = "#8a8984"
# situacoes de uma realizacao concluida que sao rupturas (sem equilibrio nem com a resistencia minima da busca do FS)
SITUACOES_RUPTURA = {s for s in os.environ.get("SITUACOES_RUPTURA", "sem_equilibrio_geostatico,fs_abaixo_de_F0").split(",")
                     if s}


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


def junta(D):
    """refaz mc_resumo.csv e mc_fs.csv da pasta (junta_mc.sh), se o script estiver ao lado"""
    j = os.path.join(os.path.dirname(os.path.abspath(__file__)), "junta_mc.sh")
    if os.path.exists(j) and any(n.startswith("r") and n[1:].isdigit() for n in os.listdir(D)):
        subprocess.run([j, D], check=False, stdout=subprocess.DEVNULL)


def num(c, k):
    """coluna numerica (None se faltar ou nao for numerica)"""
    if k not in c:
        return None
    try:
        return np.asarray(c[k], dtype=float)
    except (TypeError, ValueError):
        return None


def validas(res):
    """(entra, ruptura): as realizacoes que entram nas estatisticas (situacao ok, ou uma das SITUACOES_RUPTURA) e as que
    sao rupturas pela situacao"""
    n = len(res["real"])
    s = [str(x) for x in res.get("situacao", np.array(["ok"] * n, dtype=object))]
    rup = np.array([x in SITUACOES_RUPTURA for x in s], dtype=bool)
    ok = np.array([x in ("ok", "", "nan") for x in s], dtype=bool) | rup
    return ok, rup


def fs_minimo(res):
    """o FS_min de cada realizacao, -1 nas que romperam pela situacao (sem equilibrio)"""
    F = num(res, "FS_min")
    if F is None:
        F = num(res, "FS_SRM")
    F = np.array(F, dtype=float)
    _, rup = validas(res)
    F[rup & ~(F <= 0)] = -1.
    return F


def instante_unico(fs):
    """True se cada realizacao tem um instante so do FS em mc_fs.csv (um talude seco, um regime permanente): entao o
    FS_min e o FS, e o t_h (que pode variar de uma realizacao para outra) nao define instantes comuns"""
    r = num(fs, "real")
    if r is None or len(r) == 0:
        return True
    _, cnt = np.unique(r[np.isfinite(r)], return_counts=True)
    return len(cnt) == 0 or cnt.max() <= 1


def tabela(D):
    """(res, fs): mc_resumo.csv e mc_fs.csv da pasta; com um instante so por realizacao, as colunas de fs.csv que nao
    estao no resumo (a massa que desliza, por exemplo) entram no res, pela realizacao"""
    res = le_csv(os.path.join(D, "mc_resumo.csv"))
    arq = os.path.join(D, "mc_fs.csv")
    fs = le_csv(arq) if os.path.exists(arq) and os.path.getsize(arq) > 0 else {"real": np.array([])}
    if instante_unico(fs) and len(fs.get("real", [])):
        pos = {int(r): i for i, r in enumerate(fs["real"]) if np.isfinite(r)}
        idx = [pos.get(int(r), -1) if np.isfinite(r) else -1 for r in res["real"]]
        for k, v in fs.items():
            if k in res or k in ("real", "t_h"):
                continue
            if v.dtype == object:
                res[k] = np.array([v[i] if i >= 0 else "" for i in idx], dtype=object)
            else:
                res[k] = np.array([v[i] if i >= 0 else np.nan for i in idx])
    return res, fs


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
    return np.atleast_2d(np.corrcoef(R, rowvar=False))


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


def correlacoes(F, cols, nomes, minimo=5):
    """[(nome, n, rho1, ic1, rho_parcial, ic_parcial)] do FS com cada coluna (None ou sem variacao: fora), a parcial
    controlando as outras colunas do grupo; so realizacoes com FS > 0 e todas as colunas finitas"""
    usa = [(c, n) for c, n in zip(cols, nomes) if c is not None and np.any(np.isfinite(c))]
    if not usa:
        return []
    X = np.column_stack([F] + [c for c, _ in usa])
    ok = np.all(np.isfinite(X), axis=1) & (F > 0)
    X = X[ok]
    if len(X) < minimo or np.ptp(X[:, 0]) == 0:
        return []
    varia = [k for k in range(len(usa)) if np.ptp(X[:, k + 1]) > 0]
    if not varia:
        return []
    X = X[:, [0] + [k + 1 for k in varia]]
    nm = [usa[k][1] for k in varia]
    n = len(X)
    C = spearman(X)
    m = len(nm)
    out = []
    for k in range(m):
        r1 = C[0, k + 1]
        r2 = parcial(C, 0, k + 1) if m > 1 else r1
        out.append((nm[k], n, r1, ic95(r1, n, 0), r2, ic95(r2, n, m - 1)))
    return out


GRUPOS = [  # (rotulo, colunas c, phi, k): a definicao do Ceron, a massa que desliza e as medias no solo inteiro
    ("superficie de ruptura e massa acima dela", ("banda_c", "banda_phi_deg", "massa_k_mult")),
    ("massa que desliza", ("massa_c", "massa_phi_deg", "massa_k_mult")),
    ("medias no solo inteiro", ("c_med", "phi_med", "k_mult_med")),
]
NOMES = ("c'", "phi'", "k_v")


def colunas_grupo(c, cols):
    """(valores, nomes) das colunas de um grupo (None nas que faltam; (None, None) sem nenhuma); no grupo da banda, sem
    massa_k_mult (versao anterior), o k da banda; o grupo da banda exige a banda"""
    cols = list(cols)
    if cols[0] == "banda_c":
        if "banda_c" not in c:
            return None, None
        mk = num(c, "massa_k_mult")
        if mk is None or not np.any(np.isfinite(mk)):
            cols[2] = "banda_k_mult"
    v = [num(c, k) for k in cols]
    return (v, cols) if any(x is not None for x in v) else (None, None)


def amostra(D, t):
    """(FS, colunas, entra): o FS de cada realizacao (FS_min, ou o FS no instante t) e as colunas das propriedades"""
    junta(D)
    res, fs = tabela(D)
    if t is None:
        ok, _ = validas(res)
        return fs_minimo(res)[ok], {k: v[ok] for k, v in res.items()}
    if instante_unico(fs):
        sys.exit(f"{D}: um instante so do FS por realizacao: use sem --t (o FS_min e o FS)")
    sel = np.isclose(num(fs, "t_h"), t)
    c = {k: v[sel] for k, v in fs.items()}
    return num(c, "FS_SRM"), c


def analisa(D, t, saida_png):
    fs, c = amostra(D, t)
    linhas, resultados = [], []
    for rotulo, cols in GRUPOS:
        v, cn = colunas_grupo(c, cols)
        if v is None:
            continue
        p = correlacoes(fs, v, NOMES)
        if not p:
            continue
        nomes_col = dict(zip(NOMES, cn))
        for (nm, n, r1, i1, r2, i2) in p:
            linhas.append((rotulo, nm, nomes_col[nm], n, r1, i1[0], i1[1], r2, i2[0], i2[1]))
        resultados.append((rotulo, p))
    nome = os.path.basename(os.path.normpath(D))
    print(f"== {nome}: FS = {'FS_min' if t is None else f'FS em t = {t:g} h'}; {len(fs)} realizacoes, "
          f"{int(np.sum(fs <= 0))} sem equilibrio (FS <= 0)")
    if not linhas:
        print("  nenhuma propriedade que varie com pelo menos 5 realizacoes com FS > 0: nada calculado")
    else:
        print(f"{'propriedades':42s} {'var':5s} {'N':>5s} {'rho_XY':>8s} {'IC 95%':>16s} {'rho parcial':>11s} {'IC 95%':>16s}")
    for (rot, nm, col, n, r1, a1, b1, r2, a2, b2) in linhas:
        print(f"{rot:42s} {nm:5s} {n:5d} {r1:8.3f} [{a1:6.3f}, {b1:6.3f}] {r2:11.3f} [{a2:6.3f}, {b2:6.3f}]")
    with open(os.path.join(D, "sensibilidade.csv"), "w") as f:
        f.write("propriedades,variavel,coluna,N,rho_spearman,ic95_inf,ic95_sup,rho_parcial,ic95_parcial_inf,ic95_parcial_sup\n")
        for l in linhas:
            f.write(",".join(str(x) for x in l) + "\n")
    if resultados:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(7, 3.6), constrained_layout=True)
        h = 0.8 / len(resultados)
        for g, ((rot, p), cor) in enumerate(zip(resultados, (AZUL, LARANJA, AQUA))):
            yy = np.array([NOMES.index(q[0]) for q in p]) + (g - (len(resultados) - 1) / 2) * h
            rr = np.array([q[4] for q in p])
            err = np.array([[q[4] - q[5][0], q[5][1] - q[4]] for q in p]).T
            ax.barh(yy, rr, height=h * 0.9, color=cor, label=f"{rot} (N = {p[0][1]})", xerr=err, ecolor=CINZA, capsize=3)
        ax.axvline(0, color="#0b0b0b", lw=0.8)
        ax.set_yticks(np.arange(3), NOMES)
        ax.set_xlim(-1, 1)
        ax.set_xlabel("correlacao parcial de postos (Spearman) com o FS")
        ax.grid(axis="x", color="0.9", lw=0.6)
        ax.set_axisbelow(True)
        ax.legend(fontsize=8, loc="lower left")
        ax.set_title(f"{nome}: {'FS minimo' if t is None else f'FS em t = {t:g} h'} (barras: IC 95%)", fontsize=10)
        fig.savefig(saida_png, dpi=110)
        plt.close(fig)
        print(saida_png)


def casos(pastas, t, prefixo):
    lim = float(os.environ.get("LIMITE_FS", "1"))
    tab = []
    for D in pastas:
        fs, _ = amostra(D, t)
        fs = fs[np.isfinite(fs)]
        n = len(fs)
        if n == 0:
            continue
        pos = fs[fs > 0]
        pf = float(np.mean(fs < lim))
        tab.append((os.path.basename(os.path.normpath(D)), n, np.mean(pos) if len(pos) else np.nan,
                    np.std(pos, ddof=1) if len(pos) > 1 else np.nan, pf, np.sqrt((1 - pf) / (n * pf)) if pf > 0 else np.nan,
                    pos))
    if not tab:
        print("nenhuma realizacao concluida nas pastas")
        return
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
    com = [x for x in tab if len(x[6])]
    if com:
        xs = np.linspace(min(np.min(p) for *_, p in com) * 0.9, max(np.max(p) for *_, p in com) * 1.1, 300)
        for k, (nm, n, mu, dp, pf, cv, p) in enumerate(com[:8]):
            bw = 1.06 * np.std(p) * len(p) ** -0.2 if len(p) > 1 and np.std(p) > 0 else 0.05  # Silverman
            dens = np.exp(-0.5 * ((xs[:, None] - p[None, :]) / bw) ** 2).sum(axis=1) / (len(p) * bw * np.sqrt(2 * np.pi))
            ax[0].plot(xs, dens, color=cores[k], lw=2, label=nm)
    if len(com) > 8:
        print(f"(a densidade mostra so os 8 primeiros casos de {len(com)})")
    ax[0].axvline(lim, color=CINZA, ls="--", lw=1)
    ax[0].set_xlabel("FS minimo" if t is None else f"FS em t = {t:g} h")
    ax[0].set_ylabel("densidade de probabilidade")
    ax[0].legend(fontsize=8)
    i = np.arange(len(tab))
    pf = np.array([x[4] for x in tab])
    cv = np.array([x[5] for x in tab])
    ax[1].errorbar(i, pf, yerr=np.nan_to_num(1.96 * cv * pf), fmt="o", color="#0b0b0b", ecolor=CINZA, capsize=3)
    ax[1].set_xticks(i, [x[0] for x in tab], rotation=30, ha="right", fontsize=8)
    ax[1].set_ylabel(f"P_f = P(FS < {lim:g}) (barras: 95%)")
    ax[1].set_ylim(bottom=0)
    for a in ax:
        a.grid(color="0.9", lw=0.6)
        a.set_axisbelow(True)
    fig.savefig(prefixo + "_casos.png", dpi=110)
    plt.close(fig)
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
