#!/usr/bin/env python3
"""compara_mc.py <pasta de referencia> <pasta> [<pasta> ...] [--rotulos A,B,...] [--saida <pasta>]: compara Monte Carlos
realizacao a realizacao. A mesma realizacao (o mesmo numero, com a mesma SEMENTE) tem os mesmos numeros aleatorios em
todas as pastas de um projeto, entao a diferenca entre elas e so a da troca (o nucleo da covariancia, o angulo do
acamamento, a correlacao cruzada, ...): uma comparacao pareada, com muito menos ruido que a de dois P_f. Junta cada pasta
(junta_mc.sh), le mc_resumo.csv e mc_fs.csv e usa so as realizacoes concluidas em todas (situacao ok, ou uma ruptura de
SITUACOES_RUPTURA, que entra com FS = -1). Exemplos, no PC, depois do rsync:
  python3 compara_mc.py ~/projects/runs/monte_carlo/V2_TaludesCho/cho1 ~/projects/runs/monte_carlo/V2_TaludesCho/cho1_nucleoexp
  python3 compara_mc.py .../V4_TaludeLiu/s3 .../V4_TaludeLiu/s3_angulo90 .../V4_TaludeLiu/s3_rho-0.5 \\
      --rotulos "alfa = 0,alfa = 90,rho = -0,5"
Gera em --saida (padrao: comparacao/ ao lado da primeira pasta):
  pareado.csv      uma linha por realizacao comum: FS_min e t_min (num transiente) de cada pasta (colunas com o rotulo),
                   as medias dos campos no solo
  comparacao.png   FS_min de cada pasta contra o da referencia (quadrantes de ruptura), a diferenca, FS(t) e P_f(t)
                   medios nas realizacoes comuns (com um instante so do FS, a distribuicao acumulada do FS_min) e a
                   diferenca contra as medias dos campos no solo
  relatorio.md     por pasta: N, media, desvio e P(FS_min < limite) (1 e, se for outro, LIMITE_FS); contra a referencia:
                   a diferenca media do FS_min e o seu desvio, a correlacao de postos e a tabela de rupturas pareada (as
                   duas, so uma, nenhuma) com o teste exato de McNemar para a diferenca de P_f
Requer numpy e matplotlib."""
import math
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sensibilidade_mc as sm  # noqa: E402  (le_csv, spearman, cores)

CORES = [sm.AZUL, sm.LARANJA, sm.AQUA, sm.MAGENTA, sm.AMARELO, sm.VIOLETA, sm.VERDE, sm.VERMELHO]
TINTA, CINZA = "#0b0b0b", "#8a8984"
LIMITES = tuple(sorted({1.0, float(os.environ.get("LIMITE_FS", "1"))}))


def estilo(a, x=None, y=None, titulo=None):
    if x:
        a.set_xlabel(x)
    if y:
        a.set_ylabel(y)
    if titulo:
        a.set_title(titulo, fontsize=10)
    a.grid(color="0.9", lw=0.6)
    a.set_axisbelow(True)


def rho(x, y):
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 3 or np.ptp(x[ok]) == 0 or np.ptp(y[ok]) == 0:
        return np.nan
    return sm.spearman(np.column_stack([x[ok], y[ok]]))[0, 1]


def mcnemar(b, c):
    """p-valor exato (bilateral) de McNemar: b e c sao os pares discordantes"""
    n = b + c
    if n == 0:
        return 1.0
    k = min(b, c)
    p = sum(math.comb(n, i) for i in range(k + 1)) / 2 ** n
    return min(1.0, 2 * p)


def carrega(D):
    """{real: (FS_min, t_min)}, {real: {t: FS}}, colunas de mc_resumo, um instante so? (das que entram: ok ou ruptura).
    Com um instante so do FS por realizacao (talude seco, regime permanente: o t_h pode variar), a serie e {0: FS_min}"""
    res, fs = sm.tabela(D)
    real = np.asarray(res["real"], dtype=float)
    ok, _ = sm.validas(res)
    fsmin = sm.fs_minimo(res)
    tmin = np.asarray(res["t_min_h"], dtype=float) if sm.num(res, "t_min_h") is not None else np.full(len(real), np.nan)
    resumo = {int(r): (f, t) for r, f, t, k in zip(real, fsmin, tmin, ok) if k and np.isfinite(f) and np.isfinite(r)}
    unico = sm.instante_unico(fs)
    serie = {}
    if unico:
        serie = {r: {0.0: v[0]} for r, v in resumo.items()}
    else:
        for r, t, F in zip(fs["real"], fs["t_h"], fs["FS_SRM"]):
            if np.isfinite(r) and int(r) in resumo:
                serie.setdefault(int(r), {})[round(float(t), 6)] = float(F)
    return resumo, serie, res, unico


def main():
    args = sys.argv[1:]
    pastas, rotulos, saida = [], None, None
    while args:
        a = args.pop(0)
        if a == "--rotulos":
            rotulos = [x.strip() for x in args.pop(0).split(",")]
        elif a == "--saida":
            saida = args.pop(0)
        else:
            pastas.append(os.path.normpath(a))
    if len(pastas) < 2:
        print(__doc__)
        sys.exit(2)
    rotulos = rotulos if rotulos and len(rotulos) == len(pastas) else [os.path.basename(p) for p in pastas]
    saida = saida or os.path.join(os.path.dirname(pastas[0]), "comparacao")
    os.makedirs(saida, exist_ok=True)
    for D in pastas:
        sm.junta(D)
        if not os.path.exists(os.path.join(D, "mc_resumo.csv")):
            sys.exit(f"{D}: sem mc_resumo.csv (nenhuma realizacao concluida?)")
    dados = [carrega(D) for D in pastas]
    comuns = sorted(set.intersection(*[set(d[0]) for d in dados]))
    n = len(comuns)
    if n < 3:
        sys.exit(f"so {n} realizacoes concluidas em comum entre as pastas")
    F = np.array([[d[0][r][0] for r in comuns] for d in dados])  # pasta x realizacao
    T = np.array([[d[0][r][1] for r in comuns] for d in dados])
    # medias dos campos no solo (iguais em todas as pastas para a mesma realizacao): da primeira que as tiver
    campos = {}
    for k in ("c_med", "phi_med", "k_mult_med"):
        for _, _, res, _ in dados:
            if k in res:
                m = {int(r): v for r, v in zip(np.asarray(res["real"], dtype=float), np.asarray(res[k], dtype=float))
                     if np.isfinite(r)}
                campos[k] = np.array([m.get(r, np.nan) for r in comuns])
                break

    # ---------------- relatorio
    L = []
    L.append(f"# Comparacao pareada: {', '.join(rotulos)}\n")
    L.append(f"- Realizacoes concluidas em comum: {n} (os mesmos campos em todas as pastas).")
    L.append(f"- Referencia: {rotulos[0]} ({pastas[0]}).\n")
    L.append("| pasta | N | FS_min medio | desvio | " + " | ".join(f"P(FS_min < {lim:g})" for lim in LIMITES)
             + " | t_min mais comum |")
    L.append("|---|---|---|---|" + "---|" * len(LIMITES) + "---|")
    for j, r in enumerate(rotulos):
        f = F[j]
        p = f[f > 0]
        tm = T[j][np.isfinite(T[j])]
        moda = f"{np.bincount(np.round(tm).astype(int)).argmax()} h" if len(tm) and np.all(tm >= 0) else "-"
        L.append(f"| {r} | {n} | {np.mean(p) if len(p) else np.nan:.4f} | "
                 f"{np.std(p, ddof=1) if len(p) > 1 else np.nan:.4f} | "
                 + " | ".join(f"{100 * np.mean(f < lim):.1f}%" for lim in LIMITES) + f" | {moda} |")
    L.append("\nContra a referencia (d = FS_min da pasta - FS_min da referencia, realizacao a realizacao):\n")
    L.append("| pasta | d medio | desvio de d | d min | d max | rho (postos) | rompem nas duas | so na ref. "
             "| so na pasta | em nenhuma | McNemar p | limiar |")
    L.append("|---|---|---|---|---|---|---|---|---|---|---|---|")
    linhas_csv = []
    for j in range(1, len(pastas)):
        d = F[j] - F[0]
        for lim in LIMITES:
            a, b = F[0] < lim, F[j] < lim
            both, so_ref, so_p, nenhuma = int(np.sum(a & b)), int(np.sum(a & ~b)), int(np.sum(~a & b)), int(np.sum(~a & ~b))
            L.append(f"| {rotulos[j]} | {np.mean(d):+.4f} | {np.std(d, ddof=1):.4f} | {np.min(d):+.4f} | {np.max(d):+.4f} "
                     f"| {rho(F[0], F[j]):.3f} | {both} | {so_ref} | {so_p} | {nenhuma} | {mcnemar(so_ref, so_p):.2g} "
                     f"| FS < {lim:g} |")
    # diferenca contra as medias dos campos no solo
    if campos:
        L.append("\nCorrelacao de postos de d com as medias dos campos no solo (o efeito do modelo depende do campo?):\n")
        L.append("| pasta | " + " | ".join(campos) + " |")
        L.append("|---|" + "---|" * len(campos))
        for j in range(1, len(pastas)):
            d = F[j] - F[0]
            L.append(f"| {rotulos[j]} | " + " | ".join(f"{rho(v, d):+.3f}" for v in campos.values()) + " |")

    # FS(t) medio nas realizacoes comuns, nos instantes de cada pasta
    series = []
    for (_, serie, _, _), r in zip(dados, rotulos):
        ts = sorted(set.intersection(*[set(serie.get(k, {})) for k in comuns]))
        M = np.array([[serie[k][t] for t in ts] for k in comuns]) if ts else np.zeros((n, 0))
        series.append((np.array(ts), M))
    ts_comuns = sorted(set.intersection(*[set(s[0].tolist()) for s in series]))
    unicos = all(d[3] for d in dados)
    if ts_comuns and not unicos:
        L.append("\nFS(t) medio nas realizacoes comuns e P(FS(t) < 1):\n")
        L.append("| t (h) | " + " | ".join(rotulos) + " |")
        L.append("|---|" + "---|" * len(rotulos))
        for t in ts_comuns:
            cel = []
            for ts, M in series:
                c = M[:, list(ts).index(t)]
                cel.append(f"{np.mean(c):.4f} ({100 * np.mean(c < 1):.1f}%)")
            L.append(f"| {t:g} | " + " | ".join(cel) + " |")
    with open(os.path.join(saida, "relatorio.md"), "w") as f:
        f.write("\n".join(L) + "\n\n![comparacao.png](comparacao.png)\n")

    # ---------------- pareado.csv
    with open(os.path.join(saida, "pareado.csv"), "w") as f:
        nomes = [re.sub(r"[^0-9A-Za-z]+", "_", r).strip("_") or str(j) for j, r in enumerate(rotulos)]
        cab = ["real"] + [f"FS_min_{x}" for x in nomes] + [f"t_min_h_{x}" for x in nomes] + list(campos)
        f.write(",".join(cab) + "\n")
        for i, r in enumerate(comuns):
            v = [str(r)] + [f"{F[j, i]:.6g}" for j in range(len(pastas))] + [f"{T[j, i]:.6g}" for j in range(len(pastas))]
            v += [f"{c[i]:.6g}" for c in campos.values()]
            f.write(",".join(v) + "\n")

    # ---------------- figuras
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 2, figsize=(13, 9.5))
    a = ax[0, 0]
    # as rupturas sem equilibrio (FS = -1) ficam na borda de baixo, sem esticar a escala
    Fv = np.where(F > 0, F, np.nan)
    lo, hi = np.nanmin(Fv) - 0.02, np.nanmax(Fv) + 0.02
    a.plot([lo, hi], [lo, hi], color=CINZA, lw=1, label="1:1")
    a.axvline(1, color=TINTA, ls="--", lw=0.8)
    a.axhline(1, color=TINTA, ls="--", lw=0.8)
    for j in range(1, len(pastas)):
        a.scatter(np.clip(F[0], lo, None), np.clip(F[j], lo, None), s=12, color=CORES[(j - 1) % len(CORES)], alpha=0.7,
                  label=rotulos[j])
    a.set_xlim(lo, hi)
    a.set_ylim(lo, hi)
    estilo(a, f"FS_min: {rotulos[0]}", "FS_min", "cada realizacao (os mesmos campos)")
    a.legend(fontsize=8)
    a = ax[0, 1]
    for j in range(1, len(pastas)):
        d = F[j] - F[0]
        a.hist(d, bins=max(8, int(np.sqrt(n))), color=CORES[(j - 1) % len(CORES)], alpha=0.55,
               label=f"{rotulos[j]}: media {np.mean(d):+.3f}")
    a.axvline(0, color=TINTA, lw=0.8)
    estilo(a, f"d = FS_min - FS_min ({rotulos[0]})", "realizacoes", "diferenca do FS minimo")
    a.legend(fontsize=8)
    a = ax[1, 0]
    if unicos:
        # um instante so do FS: a distribuicao acumulada do FS_min de cada pasta nas realizacoes comuns
        for j, r in enumerate(rotulos):
            cor = TINTA if j == 0 else CORES[(j - 1) % len(CORES)]
            x = np.sort(F[j])
            a.step(x, np.arange(1, n + 1) / n, where="post", color=cor, lw=2,
                   label=f"{r}: P(FS_min < 1) = {100 * np.mean(F[j] < 1):.1f}%")
        a.axvline(1, color=TINTA, ls=":", lw=0.8)
        estilo(a, "FS_min", "P(FS_min <= x)", f"distribuicao acumulada nas {n} realizacoes comuns")
        a.legend(fontsize=8, loc="lower right")
    a2 = a.twinx() if not unicos else None
    for j, ((ts, M), r) in enumerate(zip(series, rotulos) if not unicos else []):
        if len(ts) == 0:
            continue
        cor = TINTA if j == 0 else CORES[(j - 1) % len(CORES)]
        a.plot(ts, M.mean(axis=0), "-o", ms=3, color=cor, label=f"{r}: media")
        a.fill_between(ts, np.percentile(M, 25, axis=0), np.percentile(M, 75, axis=0), color=cor, alpha=0.12)
        a2.plot(ts, 100 * np.mean(M < 1, axis=0), "--", color=cor, lw=1)
    if not unicos:
        a.axhline(1, color=TINTA, ls=":", lw=0.8)
        estilo(a, "t (h)", "FS (media e faixa 25-75%)", f"FS(t) nas {n} realizacoes comuns; tracejado: P(FS < 1)")
        a2.set_ylabel("P(FS(t) < 1) (%)")
        a2.set_ylim(0, 100)
        a.legend(fontsize=8, loc="upper right")
    a = ax[1, 1]
    ck = campos.get("k_mult_med")
    cc = campos.get("c_med")
    ck = ck if ck is not None and np.nanstd(ck) > 0 else None
    cc = cc if cc is not None and np.nanstd(cc) > 0 else None
    if ck is not None or cc is not None:
        x, nome = (ck, "k medio no solo (multiplicador de k_s)") if ck is not None else (cc, "c' medio no solo")
        for j in range(1, len(pastas)):
            d = F[j] - F[0]
            a.scatter(x, d, s=12, color=CORES[(j - 1) % len(CORES)], alpha=0.7,
                      label=f"{rotulos[j]} (rho = {rho(x, d):+.2f})")
        a.axhline(0, color=TINTA, lw=0.8)
        estilo(a, nome, f"d = FS_min - FS_min ({rotulos[0]})", "a diferenca depende do campo?")
        a.legend(fontsize=8)
    else:
        a.axis("off")
    fig.suptitle(f"Comparacao pareada ({n} realizacoes): " + " x ".join(rotulos), fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(saida, "comparacao.png"), dpi=130)
    print("\n".join(L))
    print(f"\nem {saida}: relatorio.md, pareado.csv, comparacao.png")


if __name__ == "__main__":
    main()
