#!/usr/bin/env python3
"""posprocessa_mc.py <pasta> [--det <pasta>] [--mapas] [--saida <pasta>] [--limite <FS>]: o pos-processamento de uma
pasta do Monte Carlo da fila do cluster ($SAIDA/<projeto>/<caso>[_trocas], com as realizacoes em r<nnnnn>/). Le
mc_resumo.csv e mc_fs.csv (refeitos antes por junta_mc.sh, se as realizacoes estiverem la) e, com --mapas, os CSV de
cada realizacao. Serve para qualquer projeto que siga o contrato (V2_TaludesCho, V4_TaludeLiu, ...): o que depende de
colunas que o programa nao grava fica de fora. Requer numpy e matplotlib. Gera em <pasta>/posproc/ (ou --saida):
  estatisticas.csv    para o FS minimo de cada realizacao e, num transiente (varios instantes do FS por realizacao), por
                      instante: N, media, desvio, CoV, assimetria, curtose, minimo, percentis 5, 25, 50, 75 e 95, maximo,
                      P_f = P(FS < 1) com o IC de 95% (Wilson) e CoV(P_f)
  convergencia.png    P_f, CoV(P_f), media e desvio do FS minimo contra o numero de realizacoes
  pdf.png             histograma e densidade (kernel) do FS minimo, a distribuicao acumulada e, num transiente, a
                      densidade em cada instante
  evolucao.png        (transiente) FS(t): mediana, faixas de 5-95% e 25-75%, media e o deterministico; P_f(t)
  spearman.csv/.png   correlacoes de postos de Spearman de 1a ordem e parciais (Vargas Ceron et al. 2025, eqs. 70-71,
                      Fig. 26) do FS minimo com as propriedades de cada realizacao que existem e variam: c', phi' e k na
                      superficie de ruptura, na massa que desliza e as medias no solo inteiro
  spearman_tempo.png  (transiente) as parciais em cada instante do FS
  dispersao.png       FS minimo contra essas propriedades
  mecanismo.png       onde rompe: a banda de ruptura (centroide, profundidade e area) ou a massa que desliza (centroide,
                      area e modo de ruptura)
  hidraulica.png      (transiente) infiltracao, face de percolacao e p_w maximo ao longo do tempo; FS contra a infiltracao
  diagnostico.png     tempo de calculo, iteracoes, cortes de passo e retomadas
  mapas.png           (--mapas) frequencia de cada ponto na banda de ruptura (ou na massa que desliza, de mecanismo.csv)
                      no instante do FS minimo; p_w medio e desvio (de pw.csv)
  relatorio.md        os numeros principais
O caso deterministico entra sozinho: a pasta de deterministico.txt (gravado por rodar_mc.sh) ou <pasta>_det. FS <= 0
(sem equilibrio nem com a resistencia minima da busca) e uma situacao em SITUACOES_RUPTURA (padrao
sem_equilibrio_geostatico e fs_abaixo_de_F0) contam como ruptura e ficam fora dos momentos; outra situacao != ok (o
calculo nao terminou) fica fora de tudo. Com --mapas, sem as colunas da massa (versao anterior do SlopeRandomField), o k
medio da massa e estimado pela geometria: os pontos acima da linha media da banda."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sensibilidade_mc as sm  # noqa: E402  (le_csv, tabela, validas, fs_minimo, correlacoes, cores)

AZUL, LARANJA, AQUA, AMARELO, MAGENTA, VERDE, VIOLETA, VERMELHO = sm.AZUL, sm.LARANJA, sm.AQUA, sm.AMARELO, sm.MAGENTA, \
    sm.VERDE, sm.VIOLETA, sm.VERMELHO
TINTA, CINZA = "#0b0b0b", "#8a8984"
num = sm.num
NOMES = ("c'", "phi'", "k_v")


def wilson(k, n, z=1.96):
    if n == 0:
        return (np.nan, np.nan)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0., c - h), min(1., c + h))


def estatisticas(v, limite):
    """v: FS (finitos); FS <= 0 (sem equilibrio) conta como ruptura e fica fora dos momentos"""
    v = v[np.isfinite(v)]
    n = len(v)
    pos = v[v > 0]
    k = int(np.sum(v < limite))
    pf = k / n if n else np.nan
    a, b = wilson(k, n)
    d = {"N": n, "rupturas": k, "sem_equilibrio": int(np.sum(v <= 0)), "P_f": pf, "P_f_ic95_inf": a, "P_f_ic95_sup": b,
         "CoV_P_f": np.sqrt((1 - pf) / (n * pf)) if n and pf > 0 else np.nan}
    if len(pos) > 1:
        m, s = np.mean(pos), np.std(pos, ddof=1)
        z = (pos - m) / s if s > 0 else pos * 0
        d.update({"media": m, "desvio": s, "CoV": s / m, "assimetria": np.mean(z ** 3), "curtose_excesso": np.mean(z ** 4) - 3,
                  "minimo": np.min(pos), "p5": np.percentile(pos, 5), "p25": np.percentile(pos, 25),
                  "mediana": np.percentile(pos, 50), "p75": np.percentile(pos, 75), "p95": np.percentile(pos, 95),
                  "maximo": np.max(pos)})
    return d


def kde(p, xs):
    p = p[np.isfinite(p)]
    if len(p) < 2:
        return None
    bw = 1.06 * np.std(p) * len(p) ** -0.2
    if bw <= 0:
        return None
    return np.exp(-0.5 * ((xs[:, None] - p[None, :]) / bw) ** 2).sum(axis=1) / (len(p) * bw * np.sqrt(2 * np.pi))


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


def varia(x):
    return x is not None and np.sum(np.isfinite(x)) >= 2 and np.nanmax(x) > np.nanmin(x)


def barras_spearman(ax, tab, rots):
    """barras das parciais (com o IC de 95%) por variavel, um grupo de cada cor"""
    h = 0.8 / max(1, len(rots))
    for g, (rot, cor) in enumerate(zip(rots, (AZUL, LARANJA, AQUA))):
        lin = [t for t in tab if t[0] == rot]
        yy = np.array([NOMES.index(t[1]) for t in lin]) + (g - (len(rots) - 1) / 2) * h
        rr = [t[5] for t in lin]
        err = np.array([[t[5] - t[6][0], t[6][1] - t[5]] for t in lin]).T
        ax.barh(yy, rr, height=h * 0.9, color=cor, xerr=np.nan_to_num(err), ecolor=CINZA, capsize=3,
                label=f"{rot} (N = {lin[0][2]})")
    ax.axvline(0, color=TINTA, lw=0.8)
    ax.set_yticks(np.arange(3), NOMES)
    ax.set_xlim(-1, 1)
    ax.legend(fontsize=8, loc="lower left")


def pasta_det(D):
    """a pasta do deterministico: a de deterministico.txt (rodar_mc.sh) ou <pasta>_det"""
    arq = os.path.join(D, "deterministico.txt")
    if os.path.exists(arq):
        with open(arq) as f:
            n = f.read().strip()
        c = os.path.join(os.path.dirname(D), n)
        if n and os.path.isdir(c):
            return c
    return D + "_det" if os.path.isdir(D + "_det") else None


def main():
    args = sys.argv[1:]
    # o limiar de P_f: --limite, ou LIMITE_FS (como em resumo_mc.sh e sensibilidade_mc.py), ou 1
    D, det, saida, mapas, limite = None, None, None, False, float(os.environ.get("LIMITE_FS", "1"))
    while args:
        a = args.pop(0)
        if a == "--det":
            det = args.pop(0)
        elif a == "--saida":
            saida = args.pop(0)
        elif a == "--mapas":
            mapas = True
        elif a == "--limite":
            limite = float(args.pop(0))
        else:
            D = a
    if not D:
        print(__doc__)
        sys.exit(2)
    D = os.path.normpath(D)
    det = det or pasta_det(D)
    saida = saida or os.path.join(D, "posproc")
    for d in [D] + ([det] if det else []):
        sm.junta(d)
    if not os.path.exists(os.path.join(D, "mc_resumo.csv")):
        sys.exit(f"{D}: sem mc_resumo.csv (nenhuma realizacao concluida?)")
    os.makedirs(saida, exist_ok=True)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.tri as mtri

    nome = os.path.basename(D)
    res, fsr = sm.tabela(D)
    real = num(res, "real")
    ok, rup = sm.validas(res)
    fsmin = sm.fs_minimo(res)
    unico = sm.instante_unico(fsr)
    figs, linhas_rel = [], []
    rel = linhas_rel.append

    # deterministico
    det_fsmin, det_t, det_F = None, None, None
    if det and os.path.exists(os.path.join(det, "mc_fs.csv")):
        dr = sm.le_csv(os.path.join(det, "mc_fs.csv"))
        F_det = num(dr, "FS_SRM")
        if F_det is not None and len(F_det):
            t_det = num(dr, "t_h")
            o = np.argsort(t_det) if t_det is not None else np.arange(len(F_det))
            det_t, det_F = (t_det[o] if t_det is not None else None), F_det[o]
            det_fsmin = float(np.nanmin(det_F))

    # ---------------- estatisticas
    o = np.argsort(real)
    F = fsmin[o][ok[o]]
    rl = real[o][ok[o]]
    est = {"FS_min": estatisticas(F, limite)}
    ts = np.array([]) if unico or "t_h" not in fsr else np.unique(fsr["t_h"][np.isfinite(fsr["t_h"])])
    porT = {}
    for t in ts:
        sel = np.isclose(fsr["t_h"], t)
        porT[t] = sel
        est[f"FS({t:g} h)"] = estatisticas(fsr["FS_SRM"][sel], limite)
    chaves = ["N", "rupturas", "sem_equilibrio", "media", "desvio", "CoV", "assimetria", "curtose_excesso", "minimo", "p5",
              "p25", "mediana", "p75", "p95", "maximo", "P_f", "P_f_ic95_inf", "P_f_ic95_sup", "CoV_P_f"]
    with open(os.path.join(saida, "estatisticas.csv"), "w") as f:
        f.write("grandeza," + ",".join(chaves) + "\n")
        for g, d in est.items():
            f.write(g + "," + ",".join(str(d.get(k, "")) for k in chaves) + "\n")
    e = est["FS_min"]
    rel(f"# {nome}: pos-processamento do Monte Carlo\n")
    nsit = int(np.sum(rup))
    rel(f"- Realizacoes: {len(real)} ({int(np.sum(~ok))} com o calculo incompleto, fora das estatisticas; "
        f"{e['sem_equilibrio']} sem equilibrio (FS <= 0{f', {nsit} pela situacao' if nsit else ''}), contadas como "
        f"ruptura).")
    if "media" in e:
        rel(f"- FS minimo: media {e['media']:.4f}, desvio {e['desvio']:.4f}, CoV {100 * e['CoV']:.1f}%, mediana "
            f"{e['mediana']:.4f}, percentis 5-95% [{e['p5']:.4f}, {e['p95']:.4f}], assimetria {e['assimetria']:.2f}.")
    cvtxt = "-" if not np.isfinite(e["CoV_P_f"]) else f"{100 * e['CoV_P_f']:.1f}%"
    rel(f"- P_f = P(FS_min < {limite:g}) = {100 * e['P_f']:.2f}% ({e['rupturas']} de {e['N']}), IC 95% "
        f"[{100 * e['P_f_ic95_inf']:.2f}%, {100 * e['P_f_ic95_sup']:.2f}%], CoV(P_f) = {cvtxt}.")
    if det_fsmin is not None:
        rel(f"- Deterministico (medias, {os.path.basename(det)}): FS minimo {det_fsmin:.4f}.")
    tm = num(res, "t_min_h")
    if tm is not None and len(ts):
        vals, cnt = np.unique(tm[ok][np.isfinite(tm[ok])], return_counts=True)
        rel("- Instante do FS minimo: " + ", ".join(f"{v:g} h em {c}" for v, c in zip(vals, cnt)) + ".")
    if len(ts):
        rel("\n| t (h) | N | media | desvio | CoV | P_f | IC 95% | CoV(P_f) |\n|---|---|---|---|---|---|---|---|")
        for t in ts:
            d = est[f"FS({t:g} h)"]
            cvtxt = "-" if not np.isfinite(d["CoV_P_f"]) else f"{100 * d['CoV_P_f']:.1f}%"
            rel(f"| {t:g} | {d['N']} | {d.get('media', np.nan):.4f} | {d.get('desvio', np.nan):.4f} | "
                f"{d.get('CoV', np.nan):.3f} | {100 * d['P_f']:.2f}% | [{100 * d['P_f_ic95_inf']:.2f}, "
                f"{100 * d['P_f_ic95_sup']:.2f}]% | {cvtxt} |")

    # ---------------- convergencia do Monte Carlo (na ordem das realizacoes)
    if len(F) >= 2:
        n = np.arange(1, len(F) + 1)
        falha = (F < limite).astype(float)
        pf = np.cumsum(falha) / n
        lo, hi = np.array([wilson(k, m) for k, m in zip(np.cumsum(falha), n)]).T
        cv = np.where(pf > 0, np.sqrt((1 - pf) / (n * np.where(pf > 0, pf, 1))), np.nan)
        pos = np.where(F > 0, F, np.nan)
        cnt = np.cumsum(np.isfinite(pos))
        s1 = np.nancumsum(pos)
        s2 = np.nancumsum(pos ** 2)
        mu = s1 / np.maximum(cnt, 1)
        sd = np.sqrt(np.maximum(s2 / np.maximum(cnt, 1) - mu ** 2, 0) * cnt / np.maximum(cnt - 1, 1))
        fig, ax = plt.subplots(2, 2, figsize=(12, 7.5), constrained_layout=True)
        a = ax[0, 0]
        a.fill_between(n, lo, hi, color=AZUL, alpha=0.18, lw=0, label="IC 95% (Wilson)")
        a.plot(n, pf, color=AZUL, lw=2, label="P_f")
        a.legend(fontsize=8)
        estilo(a, "numero de realizacoes", f"P_f = P(FS_min < {limite:g})", "convergencia de P_f")
        a = ax[0, 1]
        if np.any(np.isfinite(cv)):
            a.plot(n, 100 * cv, color=AZUL, lw=2, label="CoV(P_f)")
            a.set_ylim(0, min(100, np.nanmax(100 * cv) * 1.1))
        else:
            a.text(0.5, 0.5, "nenhuma ruptura ainda: CoV(P_f) indefinido", transform=a.transAxes, ha="center", fontsize=9)
            a.set_ylim(0, 100)
        a.set_xlim(1, max(2, n[-1]))
        for lim in (5, 10):
            a.axhline(lim, color=CINZA, ls="--", lw=1)
            a.text(0.99, lim, f"{lim}%", transform=a.get_yaxis_transform(), va="bottom", ha="right", fontsize=8, color=CINZA)
        estilo(a, "numero de realizacoes", "CoV(P_f) (%)", "CoV(P_f) = sqrt((1 - P_f) / (N P_f)) (Ceron, eq. 68)")
        a = ax[1, 0]
        a.plot(n, mu, color=AZUL, lw=2, label="media")
        a.legend(fontsize=8)
        estilo(a, "numero de realizacoes", "media do FS minimo", "convergencia da media")
        a = ax[1, 1]
        a.plot(n, sd, color=LARANJA, lw=2, label="desvio")
        a.legend(fontsize=8)
        estilo(a, "numero de realizacoes", "desvio do FS minimo", "convergencia do desvio")
        fig.suptitle(f"{nome}: convergencia do Monte Carlo ({len(F)} realizacoes)")
        fig.savefig(os.path.join(saida, "convergencia.png"), dpi=110)
        plt.close(fig)
        figs.append("convergencia.png")
        # quantas realizacoes para CoV(P_f) = 10% e 5%, com o P_f atual
        p = e["P_f"]
        if p > 0:
            rel(f"- Para CoV(P_f) = 10%: N = {int(np.ceil((1 - p) / (p * 0.01)))}; para 5%: N = "
                f"{int(np.ceil((1 - p) / (p * 0.0025)))} (com o P_f atual).")

    # ---------------- densidades e distribuicao acumulada
    pos = F[F > 0]
    if len(pos) >= 2:
        npan = 3 if len(ts) >= 2 else 2
        fig, ax = plt.subplots(1, npan, figsize=(5.4 * npan, 4.4), constrained_layout=True)
        xs = np.linspace(max(0, pos.min() * 0.9), pos.max() * 1.1, 400)
        a = ax[0]
        a.hist(pos, bins=max(8, int(np.sqrt(len(pos)))), density=True, color="#c7ddf5", edgecolor="white", label="histograma")
        dk = kde(pos, xs)
        if dk is not None:
            a.plot(xs, dk, color=AZUL, lw=2, label="densidade (kernel)")
        a.axvline(limite, color=TINTA, ls="--", lw=1, label=f"FS = {limite:g}")
        if det_fsmin is not None:
            a.axvline(det_fsmin, color=LARANJA, lw=2, label=f"deterministico ({det_fsmin:.3f})")
        a.legend(fontsize=8)
        estilo(a, "FS minimo", "densidade de probabilidade", f"FS minimo: P_f = {100 * e['P_f']:.2f}%")
        a = ax[1]
        xsrt = np.sort(F)
        a.step(xsrt, np.arange(1, len(xsrt) + 1) / len(xsrt), where="post", color=AZUL, lw=2, label="empirica")
        a.axvline(limite, color=TINTA, ls="--", lw=1)
        a.axhline(e["P_f"], color=CINZA, lw=1)
        a.text(a.get_xlim()[0], e["P_f"], f" P_f = {100 * e['P_f']:.2f}%", va="bottom", fontsize=8)
        estilo(a, "FS minimo", "P(FS_min <= x)", "distribuicao acumulada")
        if npan == 3:
            a = ax[2]
            cmap = plt.get_cmap("Blues")
            todos = fsr["FS_SRM"][np.isfinite(fsr["FS_SRM"]) & (fsr["FS_SRM"] > 0)]
            xs2 = np.linspace(max(0, np.min(todos) * 0.9), np.max(todos) * 1.05, 400)
            for k, t in enumerate(ts):
                v = fsr["FS_SRM"][porT[t]]
                dk = kde(v[np.isfinite(v) & (v > 0)], xs2)
                if dk is not None:
                    a.plot(xs2, dk, color=cmap(0.35 + 0.65 * k / max(1, len(ts) - 1)), lw=2, label=f"t = {t:g} h")
            a.axvline(limite, color=TINTA, ls="--", lw=1)
            if a.get_legend_handles_labels()[0]:
                a.legend(fontsize=7, ncol=2)
            estilo(a, "FS", "densidade de probabilidade", "densidade do FS em cada instante")
        fig.suptitle(nome)
        fig.savefig(os.path.join(saida, "pdf.png"), dpi=110)
        plt.close(fig)
        figs.append("pdf.png")

    # ---------------- evolucao no tempo (transiente)
    if len(ts) >= 2:
        q = {k: [] for k in ("p5", "p25", "mediana", "p75", "p95", "media", "P_f", "P_f_ic95_inf", "P_f_ic95_sup")}
        for t in ts:
            d = est[f"FS({t:g} h)"]
            for k in q:
                q[k].append(d.get(k, np.nan))
        q = {k: np.array(v) for k, v in q.items()}
        fig, ax = plt.subplots(1, 2, figsize=(13, 4.4), constrained_layout=True)
        a = ax[0]
        a.fill_between(ts, q["p5"], q["p95"], color=AZUL, alpha=0.15, lw=0, label="5-95%")
        a.fill_between(ts, q["p25"], q["p75"], color=AZUL, alpha=0.3, lw=0, label="25-75%")
        a.plot(ts, q["mediana"], color=AZUL, lw=2, label="mediana")
        a.plot(ts, q["media"], color=AZUL, lw=1.5, ls="--", label="media")
        if det_t is not None and len(det_t) >= 2:
            a.plot(det_t, det_F, "-o", color=LARANJA, lw=2, ms=4, label="deterministico")
        a.axhline(limite, color=TINTA, ls="--", lw=1)
        a.legend(fontsize=8)
        estilo(a, "t (h)", "FS", "FS ao longo do tempo")
        a = ax[1]
        a.fill_between(ts, 100 * q["P_f_ic95_inf"], 100 * q["P_f_ic95_sup"], color=AZUL, alpha=0.18, lw=0, label="IC 95%")
        a.plot(ts, 100 * q["P_f"], "-o", color=AZUL, lw=2, ms=4, label="P_f")
        a.legend(fontsize=8)
        a.set_ylim(bottom=0)
        estilo(a, "t (h)", "P_f (%)", f"P_f(t) = P(FS(t) < {limite:g})")
        fig.suptitle(nome)
        fig.savefig(os.path.join(saida, "evolucao.png"), dpi=110)
        plt.close(fig)
        figs.append("evolucao.png")

    # ---------------- massa aproximada (versao anterior) e mapas, com --mapas
    massa_aprox = None
    if mapas:
        massa_aprox, mapa_fig = mapas_e_massa(D, res, ok, saida, plt, mtri, nome)
        if mapa_fig:
            figs.append(mapa_fig)

    # ---------------- Spearman (FS minimo; propriedades no instante do minimo)
    def col(k):
        v = num(res, k)
        return v[o][ok[o]] if v is not None else None

    kcol, krot, knome = col("massa_k_mult"), "k da massa (Delta u do SRM)", "massa_k_mult"
    if kcol is None or not np.any(np.isfinite(kcol)):
        if massa_aprox is not None:
            kcol, krot, knome = np.array([massa_aprox.get(int(r), np.nan) for r in rl]), "k da massa (estimada pela banda)", \
                "massa_k_aprox"
        else:
            kcol, krot, knome = col("banda_k_mult"), "k da banda (versao anterior, sem a massa)", "banda_k_mult"
    grupos = []
    if col("banda_c") is not None:
        grupos.append((f"superficie de ruptura e {krot}", [col("banda_c"), col("banda_phi_deg"), kcol]))
    grupos.append(("massa que desliza", [col("massa_c"), col("massa_phi_deg"), col("massa_k_mult")]))
    grupos.append(("medias no solo inteiro", [col("c_med"), col("phi_med"), col("k_mult_med")]))
    tab = []
    for rot, cols in grupos:
        for (nm, nn, r1, i1, r2, i2) in sm.correlacoes(F, cols, NOMES):
            tab.append((rot, nm, nn, r1, i1, r2, i2))
    if tab:
        with open(os.path.join(saida, "spearman.csv"), "w") as f:
            f.write("propriedades,variavel,N,rho_spearman,ic95_inf,ic95_sup,rho_parcial,ic95_parcial_inf,ic95_parcial_sup\n")
            for (rot, nm, nn, r1, i1, r2, i2) in tab:
                f.write(f"{rot},{nm},{nn},{r1},{i1[0]},{i1[1]},{r2},{i2[0]},{i2[1]}\n")
        rel("\n## Correlacao de postos de Spearman com o FS minimo\n")
        rel("Parcial: controlando as outras propriedades do mesmo grupo (com tres, a de 2a ordem do Ceron).\n")
        rel("| propriedades | var | N | rho (1a ordem) | rho parcial | IC 95% |\n|---|---|---|---|---|---|")
        for (rot, nm, nn, r1, i1, r2, i2) in tab:
            rel(f"| {rot} | {nm} | {nn} | {r1:.3f} | {r2:.3f} | [{i2[0]:.3f}, {i2[1]:.3f}] |")
        fig, ax = plt.subplots(figsize=(8, 3.8), constrained_layout=True)
        barras_spearman(ax, tab, list(dict.fromkeys(t[0] for t in tab)))
        estilo(ax, "correlacao parcial de postos (Spearman) com o FS minimo", None,
               f"{nome}: sensibilidade (Ceron, Fig. 26; barras: IC 95%)")
        fig.savefig(os.path.join(saida, "spearman.png"), dpi=110)
        plt.close(fig)
        figs.append("spearman.png")

    # ---------------- Spearman em cada instante (transiente)
    if len(ts) >= 2 and "banda_c" in fsr:
        kk = "massa_k_mult" if "massa_k_mult" in fsr and np.any(np.isfinite(num(fsr, "massa_k_mult"))) else "banda_k_mult"
        serie = {nm: ([], [], [], []) for nm in NOMES}
        for t in ts:
            sel = porT[t]
            cols = [num(fsr, "banda_c")[sel], num(fsr, "banda_phi_deg")[sel],
                    num(fsr, kk)[sel] if num(fsr, kk) is not None else None]
            for (nm, nn, r1, i1, r2, i2) in sm.correlacoes(fsr["FS_SRM"][sel], cols, NOMES):
                for lista, v in zip(serie[nm], (t, r2, i2[0], i2[1])):
                    lista.append(v)
        if any(len(v[0]) >= 2 for v in serie.values()):
            fig, ax = plt.subplots(figsize=(8, 4), constrained_layout=True)
            for nm, cor in zip(NOMES, (AZUL, LARANJA, AQUA)):
                tt, r2, a1, b1 = (np.array(v) for v in serie[nm])
                if len(tt) < 2:
                    continue
                ax.fill_between(tt, a1, b1, color=cor, alpha=0.15, lw=0)
                ax.plot(tt, r2, "-o", color=cor, lw=2, ms=4,
                        label=nm + (" (massa)" if nm == "k_v" and kk == "massa_k_mult" else " (banda)"))
            ax.axhline(0, color=TINTA, lw=0.8)
            ax.set_ylim(-1, 1)
            ax.legend(fontsize=8)
            estilo(ax, "t (h)", "rho parcial com o FS(t)", f"{nome}: sensibilidade ao longo do tempo (IC 95%)")
            fig.savefig(os.path.join(saida, "spearman_tempo.png"), dpi=110)
            plt.close(fig)
            figs.append("spearman_tempo.png")

    # ---------------- dispersao: FS minimo contra cada propriedade que varia
    cand = [("banda_c", "c' na banda"), ("banda_phi_deg", "phi' na banda (graus)"), (knome, krot),
            ("massa_c", "c' na massa que desliza"), ("massa_phi_deg", "phi' na massa (graus)"),
            ("massa_k_mult", "k/k_s na massa"), ("c_med", "c' medio no solo"), ("phi_med", "phi' medio no solo (graus)"),
            ("k_mult_med", "k/k_s medio no solo")]
    pares, vistos = [], set()
    for k, rot in cand:
        if k in vistos:
            continue
        x = kcol if k == knome else col(k)
        if varia(x):
            pares.append((x, rot))
            vistos.add(k)
    if pares:
        nc = min(3, len(pares))
        nl = int(np.ceil(len(pares) / nc))
        fig, ax = plt.subplots(nl, nc, figsize=(4.6 * nc, 3.9 * nl), constrained_layout=True, squeeze=False)
        rupt = F < limite
        for a, (x, rot) in zip(ax.flat, pares):
            a.scatter(x[~rupt], F[~rupt], s=10, color=AZUL, label=f"FS >= {limite:g}")
            if np.any(rupt):
                a.scatter(x[rupt], np.where(F[rupt] > 0, F[rupt], np.nan), s=14, color=LARANJA, label=f"FS < {limite:g}")
            a.axhline(limite, color=TINTA, ls="--", lw=1)
            a.legend(fontsize=7)
            estilo(a, rot, "FS minimo", f"rho de Spearman = {rho(x, F):.2f}")
        for a in list(ax.flat)[len(pares):]:
            a.axis("off")
        fig.suptitle(nome)
        fig.savefig(os.path.join(saida, "dispersao.png"), dpi=110)
        plt.close(fig)
        figs.append("dispersao.png")

    # ---------------- mecanismo
    if col("banda_yc") is not None:
        # a banda de ruptura (SlopeRandomField: o talude de Rahardjo, 45 graus, pe em y = 20 m)
        fig, ax = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
        a = ax[0, 0]
        sc = a.scatter(col("banda_xc"), col("banda_yc"), c=np.where(F > 0, F, np.nan), s=12, cmap="Blues_r")
        fig.colorbar(sc, ax=a, label="FS minimo")
        if col("massa_xc") is not None and np.any(np.isfinite(col("massa_xc"))):
            a.scatter(col("massa_xc"), col("massa_yc"), s=8, marker="x", color=LARANJA, label="centroide da massa")
            a.legend(fontsize=8)
        a.plot([0, 30, 40, 70], [30, 30, 20, 20], color=TINTA, lw=1)
        a.set_xlim(15, 55)
        a.set_ylim(12, 32)
        a.set_aspect("equal")
        estilo(a, "x (m)", "y (m)", "centroide da banda de ruptura (cor: FS minimo)")
        a = ax[0, 1]
        ym = col("banda_ymin")
        a.hist(ym[np.isfinite(ym)], bins=20, color="#c7ddf5", edgecolor="white")
        a.axvline(20, color=CINZA, ls="--", lw=1)
        a.text(20, 0.97, "cota do pe ", transform=a.get_xaxis_transform(), fontsize=8, color=CINZA, va="top", ha="right")
        estilo(a, "ponto mais baixo da banda (m)", "realizacoes", "profundidade da superficie de ruptura")
        a = ax[1, 0]
        x = col("banda_area_m2")
        a.scatter(x, F, s=10, color=AZUL)
        a.axhline(limite, color=TINTA, ls="--", lw=1)
        estilo(a, "area da banda (m2)", "FS minimo", f"rho = {rho(x, F):.2f}")
        a = ax[1, 1]
        x = col("massa_area_m2")
        if x is not None and np.any(np.isfinite(x)):
            a.scatter(x, F, s=10, color=AZUL)
            estilo(a, "area da massa que desliza (m2)", "FS minimo", f"rho = {rho(x, F):.2f}")
        else:
            a.scatter(ym, F, s=10, color=AZUL)
            estilo(a, "ponto mais baixo da banda (m)", "FS minimo", f"rho = {rho(ym, F):.2f}")
        a.axhline(limite, color=TINTA, ls="--", lw=1)
        fig.suptitle(f"{nome}: mecanismo de ruptura no instante do FS minimo")
        fig.savefig(os.path.join(saida, "mecanismo.png"), dpi=110)
        plt.close(fig)
        figs.append("mecanismo.png")
        rel(f"\n- Banda de ruptura: ponto mais baixo em y = {np.nanmedian(ym):.2f} m (mediana; o pe esta em 20 m), "
            f"{100 * np.nanmean(ym < 19.5):.0f}% das realizacoes com a superficie abaixo do pe.")
    elif col("massa_area_m2") is not None and np.sum(np.isfinite(col("massa_area_m2"))) >= 2:
        # a massa que desliza (os pontos com o maior incremento de deslocamento no fim da busca do SRM)
        area = col("massa_area_m2")
        modo = res.get("modo")
        modo = np.asarray(modo, dtype=object)[o][ok[o]] if modo is not None else None
        fig, ax = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
        a = ax[0, 0]
        if varia(col("massa_xc")) and col("massa_yc") is not None:
            sc = a.scatter(col("massa_xc"), col("massa_yc"), c=np.where(F > 0, F, np.nan), s=12, cmap="Blues_r")
            fig.colorbar(sc, ax=a, label="FS minimo")
            a.set_aspect("equal", adjustable="datalim")
            estilo(a, "x (m)", "y (m)", "centroide da massa que desliza (cor: FS minimo)")
        else:
            a.axis("off")
        a = ax[0, 1]
        a.hist(area[np.isfinite(area)], bins=20, color="#c7ddf5", edgecolor="white")
        estilo(a, "area da massa que desliza (m2)", "realizacoes",
               f"area: media {np.nanmean(area):.1f} m2, desvio {np.nanstd(area, ddof=1):.1f} m2")
        a = ax[1, 0]
        a.scatter(area, F, s=10, color=AZUL)
        a.axhline(limite, color=TINTA, ls="--", lw=1)
        estilo(a, "area da massa que desliza (m2)", "FS minimo", f"rho = {rho(area, F):.2f}")
        a = ax[1, 1]
        if modo is not None:
            vals, cnt = np.unique([str(m) for m in modo if str(m) not in ("", "nan")], return_counts=True)
            if len(vals):
                a.bar(range(len(vals)), cnt, color=AZUL)
                a.set_xticks(range(len(vals)), vals)
                estilo(a, "modo de ruptura", "realizacoes", "modo de ruptura")
                rel("- Modo de ruptura: " + ", ".join(f"{v} em {c}" for v, c in zip(vals, cnt)) + ".")
            else:
                a.axis("off")
        else:
            a.axis("off")
        fig.suptitle(f"{nome}: a massa que desliza")
        fig.savefig(os.path.join(saida, "mecanismo.png"), dpi=110)
        plt.close(fig)
        figs.append("mecanismo.png")
        rel(f"- Massa que desliza: area media {np.nanmean(area):.1f} m2, desvio {np.nanstd(area, ddof=1):.1f} m2.")

    # ---------------- hidraulica (transiente)
    hcols = [("inf_acum_m3", "infiltracao acumulada (m3/m)"), ("face_nos", "nos na face de percolacao"),
             ("pw_max_kPa", "p_w maximo (kPa)")]
    hcols = [(k, r) for (k, r) in hcols if k in fsr and num(fsr, k) is not None and np.any(np.isfinite(num(fsr, k)))]
    if hcols and len(ts) >= 2:
        fig, ax = plt.subplots(1, len(hcols) + 1, figsize=(5 * (len(hcols) + 1), 4), constrained_layout=True)
        for a, (k, r) in zip(ax, hcols):
            v = num(fsr, k)
            p5 = [np.nanpercentile(v[porT[t]], 5) for t in ts]
            p50 = [np.nanpercentile(v[porT[t]], 50) for t in ts]
            p95 = [np.nanpercentile(v[porT[t]], 95) for t in ts]
            a.fill_between(ts, p5, p95, color=AZUL, alpha=0.18, lw=0, label="5-95%")
            a.plot(ts, p50, "-o", color=AZUL, lw=2, ms=4, label="mediana")
            a.legend(fontsize=8)
            estilo(a, "t (h)", r, r)
        a = ax[-1]
        x = col("inf_total_m3")
        if x is not None and np.any(np.isfinite(x)):
            a.scatter(x, F, s=10, color=AZUL)
            a.axhline(limite, color=TINTA, ls="--", lw=1)
            estilo(a, "infiltracao total (m3/m)", "FS minimo", f"rho = {rho(x, F):.2f}")
        else:
            a.axis("off")
        fig.suptitle(f"{nome}: hidraulica")
        fig.savefig(os.path.join(saida, "hidraulica.png"), dpi=110)
        plt.close(fig)
        figs.append("hidraulica.png")

    # ---------------- diagnostico numerico
    seg, passos, newton, cortes = col("seg"), col("passos"), col("newton"), col("cortes")
    if seg is not None and np.any(np.isfinite(seg)):
        fig, ax = plt.subplots(1, 3, figsize=(15, 3.8), constrained_layout=True)
        a = ax[0]
        a.hist(seg[np.isfinite(seg)] / 60, bins=20, color="#c7ddf5", edgecolor="white")
        estilo(a, "minutos por realizacao", "realizacoes", f"tempo: media {np.nanmean(seg) / 60:.1f} min")
        a = ax[1]
        if passos is not None and newton is not None and np.any(passos > 0):
            v = newton / np.where(passos > 0, passos, np.nan)
            a.hist(v[np.isfinite(v)], bins=20, color="#c7ddf5", edgecolor="white")
            estilo(a, "iteracoes de Newton por passo", "realizacoes", "convergencia do transiente")
        elif newton is not None and np.any(np.isfinite(newton)):
            a.hist(newton[np.isfinite(newton)], bins=20, color="#c7ddf5", edgecolor="white")
            estilo(a, "iteracoes de Newton por realizacao", "realizacoes", "iteracoes")
        else:
            a.axis("off")
        a = ax[2]
        if cortes is not None and np.any(np.isfinite(cortes)):
            vals, cnt = np.unique(cortes[np.isfinite(cortes)].astype(int), return_counts=True)
            a.bar(vals, cnt, color=AZUL)
            estilo(a, "cortes de passo por realizacao", "realizacoes", "passos divididos")
        else:
            a.axis("off")
        fig.suptitle(f"{nome}: diagnostico numerico")
        fig.savefig(os.path.join(saida, "diagnostico.png"), dpi=110)
        plt.close(fig)
        figs.append("diagnostico.png")
        rel(f"\n- Tempo: {np.nanmean(seg) / 60:.1f} min por realizacao em media"
            + (f"; cortes de passo em {int(np.nansum(cortes > 0))} realizacoes." if cortes is not None else "."))

    rel("\n## Figuras\n")
    for f in figs:
        rel(f"![{f}]({f})")
    with open(os.path.join(saida, "relatorio.md"), "w") as f:
        f.write("\n".join(linhas_rel) + "\n")
    print("\n".join(l for l in linhas_rel if not l.startswith("![")))
    print(f"\nem {saida}: estatisticas.csv, {'spearman.csv, ' if tab else ''}relatorio.md, " + ", ".join(figs))


def le_tab(arq):
    """um CSV com cabecalho (numpy estruturado, ao menos 1 linha) ou None"""
    if not os.path.exists(arq) or os.path.getsize(arq) == 0:
        return None
    try:
        m = np.atleast_1d(np.genfromtxt(arq, delimiter=",", names=True))
    except (ValueError, IndexError):
        return None
    return m if m.size and m.dtype.names else None


def triangulacao(mtri, X, Y, tri=None, idx=None):
    """a triangulacao dos nos: a da malha (triangulos.csv), ou a de Delaunay sem os triangulos de lado longo (o fecho
    convexo cobre o vazio acima da face de um talude)"""
    if tri is not None and idx is not None:
        tris = np.column_stack([[idx.get(int(v), -1) for v in tri[k]] for k in ("n0", "n1", "n2")])
        tris = tris[np.all(tris >= 0, axis=1)]
        if len(tris):
            return mtri.Triangulation(X, Y, tris)
    T = mtri.Triangulation(X, Y)
    xt, yt = X[T.triangles], Y[T.triangles]
    lados = np.max(np.hypot(xt - np.roll(xt, 1, axis=1), yt - np.roll(yt, 1, axis=1)), axis=1)
    T.set_mask(lados > 3 * np.median(lados))
    return T


def mapas_e_massa(D, res, ok, saida, plt, mtri, nome):
    """frequencia na banda de ruptura (sqrt(J2) da deformacao plastica >= 0,1 do maximo) ou na massa que desliza
    (|du| >= MASSA_LIMIAR do maximo, padrao 0,5) no instante do FS minimo (mecanismo.csv), p_w medio e desvio nesse
    instante (pw.csv) e, sem as colunas da massa e com a banda, o k medio da massa estimada pela banda (os pontos acima
    da linha media da banda)"""
    real = num(res, "real")
    tmin = num(res, "t_min_h")
    if tmin is None:
        tmin = np.full(len(real), np.nan)
    precisa_massa = num(res, "massa_k_mult") is None or not np.any(np.isfinite(num(res, "massa_k_mult")))
    limiar_du = float(os.environ.get("MASSA_LIMIAR", "0.5"))
    cont, xy = {}, {}
    pwsum, pw2, npw, nos_xy, tri = {}, {}, 0, {}, None
    massa = {}
    n, tipo = 0, "banda de ruptura"
    for r, t, k in zip(real, tmin, ok):
        if not k:
            continue
        pasta = os.path.join(D, f"r{int(r):05d}")
        m = le_tab(os.path.join(pasta, "mecanismo.csv"))
        if m is not None:
            if "t_h" in m.dtype.names and np.isfinite(t):
                m = m[np.isclose(m["t_h"], t)]
            b = None
            if m.size and "epsp_raizJ2" in m.dtype.names:
                b = m[m["epsp_raizJ2"] >= 0.1 * m["epsp_raizJ2"].max()]
            elif m.size and "du_m" in m.dtype.names and m["du_m"].max() > 0:
                b, tipo = m[m["du_m"] >= limiar_du * m["du_m"].max()], "massa que desliza"
            if b is not None and b.size:
                n += 1
                for ip, x, y in zip(b["ip"], b["x"], b["y"]):
                    cont[int(ip)] = cont.get(int(ip), 0) + 1
                    xy[int(ip)] = (x, y)
                c = le_tab(os.path.join(pasta, "campos.csv")) if precisa_massa and tipo == "banda de ruptura" else None
                if c is not None and "k_mult" in c.dtype.names:
                    # a superficie de ruptura: a linha media (mediana de y) da banda em faixas de 1 m de x; a massa sao os
                    # pontos acima dela, entre as pontas da banda
                    x0, x1 = b["x"].min(), b["x"].max()
                    xm, ym = [], []
                    for a0 in np.arange(np.floor(x0), np.ceil(x1) + 1):
                        s = (b["x"] >= a0) & (b["x"] < a0 + 1)
                        if s.any():
                            xm.append(a0 + 0.5)
                            ym.append(np.median(b["y"][s]))
                    if len(xm) >= 2:
                        ys = np.interp(c["x"], xm, ym)
                        sel = (c["x"] >= x0) & (c["x"] <= x1) & (c["y"] > ys)
                        if sel.any():
                            w = c["area"][sel] if "area" in c.dtype.names else np.ones(sel.sum())
                            massa[int(r)] = float(np.sum(w * c["k_mult"][sel]) / np.sum(w))
        p = le_tab(os.path.join(pasta, "pw.csv"))
        if p is not None and all(k2 in p.dtype.names for k2 in ("no", "x", "y", "pw_kPa")):
            if "t_h" in p.dtype.names and np.isfinite(t):
                p = p[np.isclose(p["t_h"], t)]
            for no, x, y, v in zip(p["no"], p["x"], p["y"], p["pw_kPa"]):
                no = int(no)
                pwsum[no] = pwsum.get(no, 0.) + v
                pw2[no] = pw2.get(no, 0.) + v * v
                nos_xy[no] = (x, y)
            npw += 1
            if tri is None:
                tri = le_tab(os.path.join(pasta, "triangulos.csv"))
    if n == 0 and npw == 0:
        print("--mapas: nenhuma realizacao com mecanismo.csv (MEC_CSV=1) ou pw.csv: sem mapas")
        return (massa or None), None
    paineis = (1 if n else 0) + (2 if npw else 0)
    fig, ax = plt.subplots(paineis, 1, figsize=(10, 4.2 * paineis), constrained_layout=True)
    ax = list(np.atleast_1d(ax))
    if n:
        a = ax.pop(0)
        ips = sorted(cont)
        X = np.array([xy[i][0] for i in ips])
        Y = np.array([xy[i][1] for i in ips])
        fr = np.array([cont[i] / n for i in ips])
        oo = np.argsort(fr)
        sc = a.scatter(X[oo], Y[oo], c=fr[oo], s=6, cmap="Purples", vmin=0, vmax=max(fr.max(), 1e-9))
        fig.colorbar(sc, ax=a, label="fracao das realizacoes")
        a.set_aspect("equal", adjustable="datalim")
        estilo(a, "x (m)", "y (m)", f"frequencia na {tipo} ({n} realizacoes)")
    if npw:
        nos = sorted(nos_xy)
        idx = {no: i for i, no in enumerate(nos)}
        XN = np.array([nos_xy[no][0] for no in nos])
        YN = np.array([nos_xy[no][1] for no in nos])
        mu = np.array([pwsum[no] / npw for no in nos])
        sd = np.sqrt(np.maximum(np.array([pw2[no] / npw for no in nos]) - mu ** 2, 0))
        T = triangulacao(mtri, XN, YN, tri, idx)
        import matplotlib.colors as mcolors
        lim = max(abs(mu.min()), abs(mu.max()), 1e-9)
        a = ax.pop(0)
        cf = a.tricontourf(T, mu, levels=30, cmap="RdBu", norm=mcolors.TwoSlopeNorm(0., -lim, lim))
        fig.colorbar(cf, ax=a, label="p_w medio (kPa)")
        if mu.min() < 0 < mu.max():
            a.tricontour(T, mu, levels=[0.], colors=TINTA, linewidths=1.5)
        a.set_aspect("equal")
        estilo(a, "x (m)", "y (m)", f"p_w medio no instante do FS minimo ({npw} realizacoes; linha: p_w = 0)")
        a = ax.pop(0)
        if np.ptp(sd) > 0:
            cf = a.tricontourf(T, sd, levels=30, cmap="Greens")
            fig.colorbar(cf, ax=a, label="desvio de p_w (kPa)")
        else:
            a.text(0.5, 0.5, f"p_w igual em todas as realizacoes (desvio {sd.max():.2g} kPa)", transform=a.transAxes,
                   ha="center")
        a.set_aspect("equal")
        estilo(a, "x (m)", "y (m)", "desvio de p_w entre as realizacoes")
    fig.suptitle(nome)
    fig.savefig(os.path.join(saida, "mapas.png"), dpi=110)
    plt.close(fig)
    return (massa or None), "mapas.png"


if __name__ == "__main__":
    main()
