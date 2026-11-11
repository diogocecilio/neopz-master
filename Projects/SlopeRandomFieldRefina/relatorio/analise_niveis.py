#!/usr/bin/env python3
"""analise_niveis.py: a analise do estudo dos niveis de refinamento (niveis.py) para o relatorio.

Uso: python3 analise_niveis.py [<pasta niveis> [<pasta comparacao>]]
  (padrao: as pastas da campanha no scratchpad; NIVEIS e COMPARACAO no ambiente tambem servem)

Le as rodadas de <niveis>/<inicial>/<config>/<caso> (inicial foco375 ou g15; config inicial, um_passo_h0.25,
um_passo_h0.1875, um_passo_h0.125) e as da comparacao (<comparacao>/{foco375,foco25,foco1875,foco125}/<caso>, as malhas
fixas, e <comparacao>/refinada/<caso>, a cadeia foco375 -> h = 0,125 m). A referencia de cada caso e a foco125.
Escreve, relativo a esta pasta:
  dados/niveis.csv            uma linha por configuracao e caso (fator, e, tempos da cadeia, contencao da banda)
  dados/niveis_resumo.csv     o resumo por configuracao
  dados/niveis_extrapolacao.csv  Gamma_0 e p por caso e cadeia
  tabelas/niv_resumo.tex, niv_casos.tex, niv_pf.tex, niv_extrap.tex, niv_etapa1.tex
  secoes/refino_niveis_valores.tex  os numeros citados no texto (macros \\niv...)
  figuras/niv_erro_tempo.pdf, niv_convergencia.pdf, niv_etapas.pdf, niv_contencao.pdf
O custo de CPU e comparado dentro de cada campanha: R = CPU medio da cadeia / CPU medio da rodada foco375 (etapa 1) da
mesma campanha (a nova, PARALELO=3, ou a comparacao, PARALELO=2); o 'CPU equivalente' e R vezes o CPU medio da
foco375 na campanha nova."""

import csv
import math
import os
import sys

import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(AQUI))
import compara  # noqa: E402  banda() e celulas_refino()

SCR = "/tmp/claude-0/-home-user-neopz-master/674eb771-9e25-5098-bcbd-eda87ce4d4e9/scratchpad"
NIV = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else os.environ.get("NIVEIS", os.path.join(SCR, "niveis")))
COMP = os.path.abspath(sys.argv[2] if len(sys.argv) > 2 else
                       os.environ.get("COMPARACAO", os.path.join(SCR, "refina", "comparacao")))
CASOS = ["det"] + ["r%05d" % i for i in range(8)]
ROT_CASO = {c: ("det" if c == "det" else "r%d" % int(c[1:])) for c in CASOS}
LIMIAR = 0.1
FOLGA = 0.375
CORES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
CINZA = "#7a7a7a"
# fase 1 no fermi (foco375, 1000 realizacoes): ajuste lognormal e P_f por contagem
FERMI_MEDIA, FERMI_DP, FERMI_PF_CONT = 1.457, 0.355, 0.0760
FERMI_SEG_REAL, FERMI_SIMULT = 17.4, 96


def chaves(arq):
    try:
        txt = open(arq).read().split()
    except OSError:
        return {}
    return dict(t.split("=", 1) for t in txt if "=" in t)


def num(d, k):
    try:
        return float(d[k])
    except (KeyError, TypeError, ValueError):
        return None


def linha_csv(arq):
    try:
        with open(arq) as f:
            return next(csv.DictReader(f))
    except (OSError, StopIteration):
        return {}


def contencao(pasta, geo):
    cel, folga = compara.celulas_refino(os.path.join(pasta, geo))
    b = np.array(compara.banda(pasta, LIMIAR))
    c = np.array(cel)
    d = np.sqrt(((b[:, None, :] - c[None, :, :]) ** 2).sum(-1)).min(1)
    return float(d.max()), float((d > folga + 1e-9).mean()), len(b), len(c)


# ---------------------------------------------------------------------------------------------------------------------
# configuracoes: (chave, rotulo LaTeX, rotulo curto do grafico, h, campanha, familia)
CONFIGS = [
    ("foco375", "foco375 (fixa)", "foco375", 0.375, "comp", "fixa"),
    ("foco25", "foco25 (fixa)", "foco25", 0.25, "comp", "fixa"),
    ("foco1875", "foco1875 (fixa)", "foco1875", 0.1875, "comp", "fixa"),
    ("foco125", "foco125 (fixa, ref.)", "foco125", 0.125, "comp", "fixa"),
    ("foco375_h0.25", r"foco375 $\to$ 0,25", "foco375→0,25", 0.25, "nova", "foco375"),
    ("foco375_h0.1875", r"foco375 $\to$ 0,1875", "foco375→0,1875", 0.1875, "nova", "foco375"),
    ("foco375_h0.125", r"foco375 $\to$ 0,125$^{*}$", "foco375→0,125", 0.125, "comp", "foco375"),
    ("g15", "g15 (só etapa 1)", "g15", None, "nova", "g15i"),
    ("g15_h0.25", r"g15 $\to$ 0,25", "g15→0,25", 0.25, "nova", "g15"),
    ("g15_h0.1875", r"g15 $\to$ 0,1875", "g15→0,1875", 0.1875, "nova", "g15"),
    ("g15_h0.125", r"g15 $\to$ 0,125", "g15→0,125", 0.125, "nova", "g15"),
]
ROT = {c[0]: c[1] for c in CONFIGS}
CURTO = {c[0]: c[2] for c in CONFIGS}
HCFG = {c[0]: c[3] for c in CONFIGS}
CAMP = {c[0]: c[4] for c in CONFIGS}
FAM = {c[0]: c[5] for c in CONFIGS}


def le():
    """dados[(config, caso)] = dict(fator, cpu1, cpu_malha, cpu_run, cpu, parede, triangulos, dist, fora, ...)"""
    dados = {}
    # as malhas fixas da comparacao: triangulos da tabela.csv dela
    tri_fixa = {}
    try:
        with open(os.path.join(COMP, "tabela.csv")) as f:
            for r in csv.DictReader(f):
                c = "det" if r["caso"] == "det" else "r%05d" % int(r["caso"])
                if r.get("triangulos"):
                    tri_fixa[r["metodo"], c] = int(float(r["triangulos"]))
    except OSError:
        pass
    for m in ("foco375", "foco25", "foco1875", "foco125"):
        for c in CASOS:
            p = os.path.join(COMP, m, c)
            r, t = chaves(os.path.join(p, "resultado.txt")), chaves(os.path.join(p, "tempo.txt"))
            if num(r, "FS_min") is None or num(t, "cpu_s") is None:
                continue
            dados[m, c] = dict(fator=num(r, "FS_min"), cpu1=num(t, "cpu_s"), cpu_malha=0., cpu_run=0.,
                               cpu=num(t, "cpu_s"), parede=num(t, "parede_s"), triangulos=tri_fixa.get((m, c)),
                               paralelo=num(t, "paralelo"))
    # a refinada da comparacao: foco375 -> 0,125 (o tempo e o da cadeia inteira; as etapas pelos tempos do programa)
    for c in CASOS:
        p = os.path.join(COMP, "refinada", c)
        r, t, rf = chaves(os.path.join(p, "resultado.txt")), chaves(os.path.join(p, "tempo.txt")), \
            linha_csv(os.path.join(p, "refino.csv"))
        if num(r, "FS_min") is None or num(t, "cpu_s") is None or num(rf, "seg_fina") is None:
            continue
        sg, sm, sf = num(rf, "seg_grossa"), num(rf, "seg_malha"), num(rf, "seg_fina")
        cpu = num(t, "cpu_s")
        tot = sg + sm + sf
        v = dict(fator=r and num(r, "FS_min"), cpu=cpu, parede=num(t, "parede_s"), cpu1=cpu * sg / tot,
                 cpu_malha=cpu * sm / tot, cpu_run=cpu * sf / tot, triangulos=int(num(rf, "triangulos")),
                 celulas_banda=int(num(rf, "celulas_banda")), paralelo=num(t, "paralelo"), fator_pai=num(rf, "fator_grossa"))
        try:
            v["dist"], v["fora"], v["nbanda"], v["ncel"] = contencao(p, "malha_refinada.geo")
        except (OSError, ValueError, IndexError) as e:
            print("analise_niveis: contencao refinada %s: %s" % (c, e), file=sys.stderr)
        dados["foco375_h0.125", c] = v
    # a campanha nova
    for ini in ("foco375", "g15"):
        for c in CASOS:
            p1 = os.path.join(NIV, ini, "inicial", c)
            r1, t1 = chaves(os.path.join(p1, "resultado.txt")), chaves(os.path.join(p1, "tempo.txt"))
            if num(r1, "FS_min") is None or num(t1, "cpu_s") is None:
                continue
            m1 = chaves(os.path.join(p1, "malha.txt"))
            base = dict(fator=num(r1, "FS_min"), cpu1=num(t1, "cpu_s"), cpu_malha=0., cpu_run=0., cpu=num(t1, "cpu_s"),
                        parede=num(t1, "parede_s"), parede1=num(t1, "parede_s"), triangulos=int(num(m1, "triangulos")),
                        paralelo=num(t1, "paralelo"))
            dados[(ini + "_ini") if ini == "foco375" else ini, c] = base
            for h in ("0.25", "0.1875", "0.125"):
                p = os.path.join(NIV, ini, "um_passo_h" + h, c)
                r, t = chaves(os.path.join(p, "resultado.txt")), chaves(os.path.join(p, "tempo.txt"))
                if num(r, "FS_min") is None or num(t, "cpu_s") is None or not os.path.isfile(os.path.join(p, "tempo.txt")):
                    continue
                if ini == "foco375" and h == "0.125":
                    continue  # a de h = 0,125 da foco375 e a refinada da comparacao
                mt = chaves(os.path.join(p, "malha.txt"))
                v = dict(fator=num(r, "FS_min"), cpu1=base["cpu1"], cpu_malha=num(t, "cpu_malha_s") or 0.,
                         cpu_run=num(t, "cpu_s"), parede=base["parede1"] + (num(t, "seg_malha") or 0.) + num(t, "parede_s"),
                         triangulos=int(num(mt, "triangulos")), celulas_banda=int(num(mt, "celulas_banda")),
                         paralelo=num(t, "paralelo"), fator_pai=num(mt, "fator_grossa"))
                v["cpu"] = v["cpu1"] + v["cpu_malha"] + v["cpu_run"]
                try:
                    v["dist"], v["fora"], v["nbanda"], v["ncel"] = contencao(p, "malha.geo")
                except (OSError, ValueError, IndexError) as e:
                    print("analise_niveis: contencao %s %s %s: %s" % (ini, h, c, e), file=sys.stderr)
                dados["%s_h%s" % (ini, h), c] = v
    for (k, c), v in dados.items():
        ref = dados.get(("foco125", c))
        v["e"] = v["fator"] / ref["fator"] - 1. if ref else None
    return dados


# ---------------------------------------------------------------------------------------------------------------------
# formatos

def f(x, d=2, sinal=False):
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return "--"
    s = ("%+.*f" if sinal else "%.*f") % (d, x)
    s = s.replace(".", "{,}")
    return s.replace("-", "$-$").replace("+", "$+$") if sinal else s.replace("-", "$-$")


def fp(x, d=2, sinal=False):  # em porcentagem
    return f(None if x is None else 100. * x, d, sinal)


def ft(x, d=2):  # texto corrido (sem chaves do siunitx)
    return f(x, d).replace("{,}", ",")


def virg(x, d):
    return ("%.*f" % (d, x)).replace(".", ",").replace("-", "−")


def media(v):
    v = [x for x in v if x is not None]
    return sum(v) / len(v) if v else None


# ---------------------------------------------------------------------------------------------------------------------
# a extrapolacao Gamma(h) = Gamma_0 + C h^p com tres niveis

def ajuste_p(h, g):
    """(p, Gamma_0, situacao): p de (g1-g2)/(g2-g3) = (h1^p-h2^p)/(h2^p-h3^p), por bissecao em (0,05; 12)"""
    (h1, h2, h3), (g1, g2, g3) = h, g
    d12, d23 = g1 - g2, g2 - g3
    if d23 == 0 or d12 == 0 or (d12 > 0) != (d23 > 0):
        return None, None, "nao monotona"
    r = d12 / d23
    q = lambda p: (h1 ** p - h2 ** p) / (h2 ** p - h3 ** p)
    a, b = 0.05, 12.
    if not (q(a) < r < q(b)):
        return None, None, ("p < 0,05" if r <= q(a) else "p > 12")
    for _ in range(200):
        mid = 0.5 * (a + b)
        if q(mid) < r:
            a = mid
        else:
            b = mid
    p = 0.5 * (a + b)
    C = d23 / (h2 ** p - h3 ** p)
    return p, g3 - C * h3 ** p, "ok"


def extrap_p(h, g, p):
    """Gamma_0 com p dado, pelos dois niveis mais finos"""
    (h2, h3), (g2, g3) = h[1:], g[1:]
    C = (g2 - g3) / (h2 ** p - h3 ** p)
    return g3 - C * h3 ** p


# ---------------------------------------------------------------------------------------------------------------------

def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter
    from scipy.stats import norm

    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42,
                         "axes.grid": True, "grid.color": "#e4e4e4", "grid.linewidth": 0.6, "legend.frameon": False})
    vir = lambda d: FuncFormatter(lambda x, _pos: virg(x, d))

    dados = le()
    for d in ("dados", "tabelas", "figuras", "secoes"):
        os.makedirs(os.path.join(AQUI, d), exist_ok=True)
    feitos = {k: [c for c in CASOS if (k, c) in dados] for k in [x[0] for x in CONFIGS] + ["foco375_ini"]}
    # os casos com todas as configuracoes
    completos = [c for c in CASOS if all((k, c) in dados for k in [x[0] for x in CONFIGS] + ["foco375_ini"])]
    print("analise_niveis: casos completos: %s" % ", ".join(ROT_CASO[c] for c in completos))
    for k, v in feitos.items():
        if len(v) < len(CASOS):
            print("  %s: %d casos (%s)" % (k, len(v), ", ".join(ROT_CASO[c] for c in v)))
    casos = completos

    # ----- custo: base foco375 de cada campanha (nos casos usados)
    cpu375 = {"nova": media([dados["foco375_ini", c]["cpu"] for c in casos]),
              "comp": media([dados["foco375", c]["cpu"] for c in casos])}
    par375 = {"nova": media([dados["foco375_ini", c]["parede"] for c in casos]),
              "comp": media([dados["foco375", c]["parede"] for c in casos])}
    razao_etapa1 = [dados["foco375_ini", c]["cpu"] / dados["foco375", c]["cpu"] for c in casos]

    resumo = {}
    for k, *_ in CONFIGS:
        vs = [dados[k, c] for c in casos]
        es = [v["e"] for v in vs]
        R = media([v["cpu"] for v in vs]) / cpu375[CAMP[k]]
        resumo[k] = dict(
            n=len(vs), triangulos=media([v["triangulos"] for v in vs]), e_abs=media([abs(e) for e in es]),
            e_max=max(abs(e) for e in es), e_med=media(es), cpu=media([v["cpu"] for v in vs]),
            parede=media([v["parede"] for v in vs]), cpu1=media([v["cpu1"] for v in vs]),
            cpu_malha=media([v["cpu_malha"] for v in vs]), cpu_run=media([v["cpu_run"] for v in vs]), R=R,
            dist=max([v["dist"] for v in vs if v.get("dist") is not None], default=None),
            dist_med=media([v.get("dist") for v in vs]),
            fora=media([v.get("fora") for v in vs]), fora_max=max([v["fora"] for v in vs if v.get("fora") is not None],
                                                                   default=None),
            casos_fora=sum(1 for v in vs if v.get("fora")), paralelo=media([v.get("paralelo") for v in vs]))
    for k in resumo:
        resumo[k]["R125"] = resumo[k]["R"] / resumo["foco125"]["R"]
        resumo[k]["cpu_eq"] = resumo[k]["R"] * cpu375["nova"]
        resumo[k]["cpu1_eq"] = resumo[k]["cpu1"] / cpu375[CAMP[k]] * cpu375["nova"]
        resumo[k]["cpu_malha_eq"] = resumo[k]["cpu_malha"] / cpu375[CAMP[k]] * cpu375["nova"]
        resumo[k]["cpu_run_eq"] = resumo[k]["cpu_run"] / cpu375[CAMP[k]] * cpu375["nova"]

    # ----- csv por caso
    with open(os.path.join(AQUI, "dados", "niveis.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["config", "campanha", "h", "caso", "fator", "e", "triangulos", "celulas_banda", "fator_pai",
                    "cpu_etapa1_s", "cpu_malha_s", "cpu_refinada_s", "cpu_cadeia_s", "parede_cadeia_s", "paralelo",
                    "dist_max_m", "fora", "pontos_banda", "pontos_malha"])
        for k, *_ in CONFIGS + [("foco375_ini",)]:
            for c in CASOS:
                v = dados.get((k, c))
                if v is None:
                    continue
                w.writerow([k, CAMP.get(k, "nova"), "" if HCFG.get(k) is None else HCFG[k], ROT_CASO[c],
                            "%.9f" % v["fator"], "" if v["e"] is None else "%.6f" % v["e"], v.get("triangulos", ""),
                            v.get("celulas_banda", ""), v.get("fator_pai", "") or "", "%.3f" % v["cpu1"],
                            "%.3f" % v["cpu_malha"], "%.3f" % v["cpu_run"], "%.3f" % v["cpu"], "%.3f" % v["parede"],
                            "" if v.get("paralelo") is None else "%g" % v["paralelo"],
                            "" if v.get("dist") is None else "%.4f" % v["dist"],
                            "" if v.get("fora") is None else "%.4f" % v["fora"], v.get("nbanda", ""), v.get("ncel", "")])
    with open(os.path.join(AQUI, "dados", "niveis_resumo.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        cols = ["n", "triangulos", "e_abs", "e_max", "e_med", "cpu", "parede", "cpu1", "cpu_malha", "cpu_run", "R",
                "R125", "cpu_eq", "dist", "fora", "fora_max", "casos_fora"]
        w.writerow(["config", "campanha"] + cols)
        for k, *_ in CONFIGS:
            w.writerow([k, CAMP[k]] + ["" if resumo[k][x] is None else "%.6g" % resumo[k][x] for x in cols])

    # ----- extrapolacao
    H3 = (0.25, 0.1875, 0.125)
    cadeias = {"foco375": ["foco375_h0.25", "foco375_h0.1875", "foco375_h0.125"],
               "g15": ["g15_h0.25", "g15_h0.1875", "g15_h0.125"],
               "fixas": ["foco25", "foco1875", "foco125"]}
    ext = {}
    for nome, ks in cadeias.items():
        for c in casos:
            g = [dados[k, c]["fator"] for k in ks]
            p, g0, sit = ajuste_p(H3, g)
            g0_2 = extrap_p(H3, g, 2.)
            ext[nome, c] = dict(p=p, g0=g0, sit=sit, g0_2=g0_2, g=g)
    with open(os.path.join(AQUI, "dados", "niveis_extrapolacao.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["cadeia", "caso", "fator_h0.25", "fator_h0.1875", "fator_h0.125", "p", "gamma0_p", "situacao",
                    "gamma0_p2", "gamma0_p2_sobre_foco125_menos1"])
        for (nome, c), x in ext.items():
            w.writerow([nome, ROT_CASO[c]] + ["%.9f" % y for y in x["g"]] +
                       ["" if x["p"] is None else "%.4f" % x["p"], "" if x["g0"] is None else "%.9f" % x["g0"], x["sit"],
                        "%.9f" % x["g0_2"], "%.6f" % (x["g0_2"] / dados["foco125", c]["fator"] - 1.)])
    ps = {nome: [ext[nome, c]["p"] for c in casos if ext[nome, c]["p"] is not None] for nome in cadeias}
    # Gamma_0 (p = 2) em relacao a foco125: o vies estimado da propria referencia
    b0 = {nome: [ext[nome, c]["g0_2"] / dados["foco125", c]["fator"] - 1. for c in casos] for nome in cadeias}

    # ----- P_f
    cov = FERMI_DP / FERMI_MEDIA
    s_ln = math.sqrt(math.log(1. + cov ** 2))
    m_ln = math.log(FERMI_MEDIA) - 0.5 * s_ln ** 2
    b375 = resumo["foco375"]["e_med"]

    def pf_config(bx):  # P_f com a configuracao de vies bx (relativo a foco125), a partir do ajuste da foco375
        return norm.cdf((math.log(1. + b375) - math.log(1. + bx) - m_ln) / s_ln)

    pf_ref = pf_config(0.)
    pf_ln375 = norm.cdf(-m_ln / s_ln)
    cv = lambda pf, n: math.sqrt((1. - pf) / (n * pf))
    pf_rows = []
    for k, *_ in CONFIGS:
        bx = resumo[k]["e_med"]
        pf = pf_config(bx)
        pf_rows.append((k, bx, pf, pf - pf_ref, pf / pf_ref - 1.))
    # Gamma_0 estimado (p = 2, cadeia das fixas) como referencia alternativa: vies da foco125 em relacao a ele
    b125_0 = media(b0["fixas"])

    # ----- fermi
    horas = lambda R, n: n * FERMI_SEG_REAL * R / FERMI_SIMULT / 3600.

    # ----- tabelas
    # niv_resumo
    L = [r"\begin{tabular}{@{}l r r r r r r r r r@{}}", r"\toprule",
         r" & & \multicolumn{3}{c}{$e$ (\%)} & \multicolumn{3}{c}{custo da cadeia por caso} & "
         r"\multicolumn{2}{c}{contenção da banda} \\",
         r"\cmidrule(lr){3-5}\cmidrule(lr){6-8}\cmidrule(l){9-10}",
         r"Configuração & triâng. & $\overline{|e|}$ & $\max|e|$ & $\bar e$ & CPU (s) & parede (s) & "
         r"$R/R_{125}$ & $d_{\max}$ (m) & fora (\%) \\", r"\midrule"]
    for k, *_ in CONFIGS:
        x = resumo[k]
        if k in ("foco375_h0.25", "g15"):
            L.append(r"\midrule")
        ref = k == "foco125"
        L.append(" & ".join([ROT[k], "%d" % round(x["triangulos"]),
                             "--" if ref else fp(x["e_abs"], 2), "--" if ref else fp(x["e_max"], 2),
                             "--" if ref else fp(x["e_med"], 2, True), f(x["cpu"], 1), f(x["parede"], 1),
                             f(x["R125"], 3), f(x["dist"], 3) if x["dist"] is not None else "--",
                             fp(x["fora"], 2) if x["fora"] is not None else "--"]) + r" \\")
    L += [r"\bottomrule", r"\end{tabular}"]
    open(os.path.join(AQUI, "tabelas", "niv_resumo.tex"), "w").write("\n".join(L) + "\n")

    # niv_casos: linhas = configuracoes, colunas = casos
    L = [r"\begin{tabular}{@{}l " + "r" * len(casos) + r"@{}}", r"\toprule",
         "Configuração & " + " & ".join(r"\multicolumn{1}{c}{%s}" % ROT_CASO[c] for c in casos) + r" \\", r"\midrule"]
    for k, *_ in CONFIGS:
        if k in ("foco375_h0.25", "g15"):
            L.append(r"\midrule")
        L.append(ROT[k] + " & " + " & ".join(f(dados[k, c]["fator"], 4) for c in casos) + r" \\")
    L.append(r"\midrule")
    L.append(r"\multicolumn{%d}{@{}l}{\emph{Diferença para a foco125, $e$ (\%%)}} \\" % (len(casos) + 1))
    for k, *_ in CONFIGS:
        if k == "foco125":
            continue
        L.append(ROT[k] + " & " + " & ".join(fp(dados[k, c]["e"], 2, True) for c in casos) + r" \\")
    L += [r"\bottomrule", r"\end{tabular}"]
    open(os.path.join(AQUI, "tabelas", "niv_casos.tex"), "w").write("\n".join(L) + "\n")

    # niv_pf
    L = [r"\begin{tabular}{@{}l r r r r@{}}", r"\toprule",
         r"Configuração & $\bar e$ (\%) & $P_f$ (\%) & $\Delta P_f$ (p.p.) & $\Delta P_f/P_f$ (\%) \\", r"\midrule"]
    for k, bx, pf, dpf, rel in pf_rows:
        if k in ("foco375_h0.25", "g15"):
            L.append(r"\midrule")
        L.append(" & ".join([ROT[k], fp(bx, 2, True), fp(pf, 2), fp(dpf, 2, True), fp(rel, 1, True)]) + r" \\")
    pf0 = pf_config(b125_0)
    L += [r"\midrule", " & ".join([r"$\Gamma_0$ estimado ($p = 2$, fixas)", fp(b125_0, 2, True), fp(pf0, 2),
                                    fp(pf0 - pf_ref, 2, True), fp(pf0 / pf_ref - 1., 1, True)]) + r" \\"]
    L += [r"\midrule",
          r"\multicolumn{5}{@{}l}{\emph{Coeficiente de variação de Monte Carlo de $P_f$ (com $P_f$ = %s\,\%%)}} \\"
          % fp(pf_ref, 2)]
    for n in (1000, 5000):
        L.append(r"$N = %s$ & \multicolumn{4}{r@{}}{%s\,\%%} \\" % ("\\num{%d}" % n, fp(cv(pf_ref, n), 1)))
    L += [r"\bottomrule", r"\end{tabular}"]
    open(os.path.join(AQUI, "tabelas", "niv_pf.tex"), "w").write("\n".join(L) + "\n")

    # niv_extrap: por caso, p e Gamma_0/foco125 - 1 de cada cadeia
    L = [r"\begin{tabular}{@{}l r r r r r r@{}}", r"\toprule",
         r" & \multicolumn{2}{c}{cadeia foco375} & \multicolumn{2}{c}{cadeia g15} & "
         r"\multicolumn{2}{c}{malhas fixas} \\", r"\cmidrule(lr){2-3}\cmidrule(lr){4-5}\cmidrule(l){6-7}",
         r"Caso & $p$ & $\Gamma_0^{(2)}$ (\%) & $p$ & $\Gamma_0^{(2)}$ (\%) & $p$ & $\Gamma_0^{(2)}$ (\%) \\",
         r"\midrule"]
    def pcel(x):
        if x["p"] is not None:
            return f(x["p"], 2)
        return {"nao monotona": r"n.\,m.", "p < 0,05": r"$<0{,}05$", "p > 12": r"$>12$"}[x["sit"]]
    for c in casos:
        cel = [ROT_CASO[c]]
        for nome in ("foco375", "g15", "fixas"):
            x = ext[nome, c]
            cel += [pcel(x), fp(x["g0_2"] / dados["foco125", c]["fator"] - 1., 2, True)]
        L.append(" & ".join(cel) + r" \\")
    L.append(r"\midrule")
    cel = ["média"]
    for nome in ("foco375", "g15", "fixas"):
        cel += [(f(media(ps[nome]), 2) + (r"$^{\,(%d)}$" % len(ps[nome]) if len(ps[nome]) < len(casos) else ""))
                if ps[nome] else "--", fp(media(b0[nome]), 2, True)]
    L.append(" & ".join(cel) + r" \\")
    L += [r"\bottomrule", r"\end{tabular}"]
    open(os.path.join(AQUI, "tabelas", "niv_extrap.tex"), "w").write("\n".join(L) + "\n")

    # niv_etapa1: a conferencia do CPU das duas rodadas foco375 (etapa 1)
    L = [r"\begin{tabular}{@{}l " + "r" * len(casos) + r" r@{}}", r"\toprule",
         "CPU (s) & " + " & ".join(r"\multicolumn{1}{c}{%s}" % ROT_CASO[c] for c in casos) + r" & média \\",
         r"\midrule",
         r"nova (\texttt{PARALELO=3}) & " + " & ".join(f(dados["foco375_ini", c]["cpu"], 2) for c in casos) +
         " & " + f(cpu375["nova"], 2) + r" \\",
         r"comparação (\texttt{PARALELO=2}) & " + " & ".join(f(dados["foco375", c]["cpu"], 2) for c in casos) +
         " & " + f(cpu375["comp"], 2) + r" \\",
         r"razão nova/comparação & " + " & ".join(f(x, 3) for x in razao_etapa1) + " & " +
         f(cpu375["nova"] / cpu375["comp"], 3) + r" \\",
         r"\bottomrule", r"\end{tabular}"]
    open(os.path.join(AQUI, "tabelas", "niv_etapa1.tex"), "w").write("\n".join(L) + "\n")

    # niv_fermi: horas no fermi
    L = [r"\begin{tabular}{@{}l r r r@{}}", r"\toprule",
         r"Configuração & $R$ & \num{1000} realizações (h) & \num{5000} realizações (h) \\", r"\midrule"]
    for k in ("foco375", "foco375_h0.1875", "foco375_h0.125", "g15_h0.125", "foco125"):
        L.append(" & ".join([ROT[k], f(resumo[k]["R"], 2), f(horas(resumo[k]["R"], 1000), 2),
                             f(horas(resumo[k]["R"], 5000), 2)]) + r" \\")
    L += [r"\bottomrule", r"\end{tabular}"]
    open(os.path.join(AQUI, "tabelas", "niv_fermi.tex"), "w").write("\n".join(L) + "\n")

    # ----- macros para o texto
    M = {}
    dif1 = [abs(dados["foco375_ini", c]["fator"] / dados["foco375", c]["fator"] - 1.) for c in CASOS
            if ("foco375_ini", c) in dados and ("foco375", c) in dados]
    M["nivDifEtapaUm"] = ("%.1e" % max(dif1)).replace(".", ",") if max(dif1) > 0 else "0"
    def carga(pasta):
        t = chaves(os.path.join(pasta, "tempo.txt"))
        v = [x for x in (num(t, "carga_inicio"), num(t, "carga_fim")) if x is not None]
        return media(v)
    cg_n = [carga(os.path.join(NIV, ini, cfg, c)) for ini in ("foco375", "g15")
            for cfg in ("inicial", "um_passo_h0.25", "um_passo_h0.1875", "um_passo_h0.125") for c in CASOS
            if os.path.isfile(os.path.join(NIV, ini, cfg, c, "tempo.txt"))]
    cg_c = [carga(os.path.join(COMP, m, c)) for m in ("foco375", "foco25", "foco1875", "foco125", "refinada")
            for c in CASOS if os.path.isfile(os.path.join(COMP, m, c, "tempo.txt"))]
    M["nivCargaNova"] = ft(media([x for x in cg_n if x is not None]), 1)
    M["nivCargaComp"] = ft(media([x for x in cg_c if x is not None]), 1)
    M["nivNcasosEtapaUm"] = "%d" % len(dif1)
    pf0 = pf_config(b125_0)
    M["nivPfZero"] = ft(100 * pf0, 2)
    M["nivDpfZero"] = ft(100 * (pf0 / pf_ref - 1.), 1)
    M["nivPfCont"] = ft(100 * FERMI_PF_CONT, 2)
    M["nivCvMilCont"] = ft(100 * cv(FERMI_PF_CONT, 1000), 1)
    M["nivCvCincoCont"] = ft(100 * cv(FERMI_PF_CONT, 5000), 1)
    def area_banda(pasta):
        try:
            with open(os.path.join(pasta, "fs.csv")) as fh:
                fs = [r for r in csv.DictReader(fh) if r.get("FS_SRM")]
            return float(min(fs, key=lambda r: float(r["FS_SRM"]))["banda_area_m2"])
        except (OSError, ValueError, KeyError):
            return None
    M["nivAreaF"] = ft(media([area_banda(os.path.join(NIV, "foco375", "inicial", c)) for c in casos]), 1)
    M["nivAreaG"] = ft(media([area_banda(os.path.join(NIV, "g15", "inicial", c)) for c in casos]), 1)
    M["nivCpuGI"] = ft(resumo["g15"]["cpu"], 1)
    M["nivRazEtapaGF"] = ft(cpu375["nova"] / resumo["g15"]["cpu"], 0)
    for k, nm in (("foco375_h0.125", "FC"), ("g15_h0.125", "GC"), ("foco375_h0.1875", "FB"), ("g15_h0.25", "GA"),
                  ("foco375_h0.25", "FA"), ("g15_h0.1875", "GB")):
        M["niv%sCel" % nm] = "%d" % round(media([dados[k, c].get("celulas_banda") for c in casos]))
    M["nivNcasos"] = "%d" % len(casos)
    M["nivCasos"] = ", ".join(ROT_CASO[c] for c in casos)
    M["nivCpuNova"] = ft(cpu375["nova"], 1)
    M["nivCpuComp"] = ft(cpu375["comp"], 1)
    M["nivParNova"] = ft(par375["nova"], 1)
    M["nivParComp"] = ft(par375["comp"], 1)
    M["nivRazEtapa"] = ft(cpu375["nova"] / cpu375["comp"], 2)
    M["nivRazEtapaMin"] = ft(min(razao_etapa1), 2)
    M["nivRazEtapaMax"] = ft(max(razao_etapa1), 2)
    nomes = {"foco375": "FfixA", "foco25": "FfixB", "foco1875": "FfixC", "foco125": "FfixD",
             "foco375_h0.25": "FA", "foco375_h0.1875": "FB", "foco375_h0.125": "FC", "g15": "GI",
             "g15_h0.25": "GA", "g15_h0.1875": "GB", "g15_h0.125": "GC"}
    for k, nm in nomes.items():
        x = resumo[k]
        M["niv%sEabs" % nm] = ft(100 * x["e_abs"], 2)
        M["niv%sEmax" % nm] = ft(100 * x["e_max"], 2)
        M["niv%sEmed" % nm] = ft(100 * x["e_med"], 2).replace("−", "$-$") if x["e_med"] < 0 else \
            "$+$" + ft(100 * x["e_med"], 2)
        M["niv%sCpu" % nm] = ft(x["cpu"], 1)
        M["niv%sCpuEq" % nm] = ft(x["cpu_eq"], 1)
        M["niv%sPar" % nm] = ft(x["parede"], 1)
        M["niv%sR" % nm] = ft(x["R"], 2)
        M["niv%sRcent" % nm] = ft(x["R125"], 3)
        M["niv%sRpct" % nm] = ft(100 * x["R125"], 1)
        M["niv%sTri" % nm] = "%d" % round(x["triangulos"])
        M["niv%sDist" % nm] = ft(x["dist"], 2) if x["dist"] is not None else "--"
        M["niv%sFora" % nm] = ft(100 * x["fora"], 2) if x["fora"] is not None else "--"
        M["niv%sForaMax" % nm] = ft(100 * x["fora_max"], 1) if x["fora_max"] is not None else "--"
        M["niv%sCasosFora" % nm] = "%d" % x["casos_fora"]
        M["niv%sHmil" % nm] = ft(horas(x["R"], 1000), 2)
        M["niv%sHcinco" % nm] = ft(horas(x["R"], 5000), 2)
        M["niv%sMinmil" % nm] = ft(60 * horas(x["R"], 1000), 0)
        M["niv%sMincinco" % nm] = ft(60 * horas(x["R"], 5000), 0)
        bx = x["e_med"]
        pf = pf_config(bx)
        M["niv%sPf" % nm] = ft(100 * pf, 2)
        M["niv%sDpf" % nm] = ft(100 * (pf / pf_ref - 1.), 1)
        M["niv%sDpfAbs" % nm] = ft(100 * abs(pf / pf_ref - 1.), 1)
    M["nivSln"] = ft(s_ln, 4)
    M["nivMln"] = ft(m_ln, 4)
    M["nivPfLn"] = ft(100 * pf_ln375, 2)
    M["nivPfRef"] = ft(100 * pf_ref, 2)
    M["nivCvMil"] = ft(100 * cv(pf_ref, 1000), 1)
    M["nivCvCinco"] = ft(100 * cv(pf_ref, 5000), 1)
    M["nivBzero"] = ft(100 * abs(b125_0), 2)
    M["nivBzeroSinal"] = ("$-$" if b125_0 < 0 else "$+$") + ft(100 * abs(b125_0), 2)
    for nome, nm in (("foco375", "F"), ("g15", "G"), ("fixas", "X")):
        M["nivP%sMed" % nm] = ft(media(ps[nome]), 2) if ps[nome] else "--"
        M["nivP%sN" % nm] = "%d" % len(ps[nome])
        M["nivP%sMin" % nm] = ft(min(ps[nome]), 2) if ps[nome] else "--"
        M["nivP%sMax" % nm] = ft(max(ps[nome]), 2) if ps[nome] else "--"
        bb = b0[nome]
        M["nivB%sMed" % nm] = ("$-$" if media(bb) < 0 else "$+$") + ft(100 * abs(media(bb)), 2)
        M["nivB%sMin" % nm] = ft(100 * min(bb), 2).replace("−", "$-$")
        M["nivB%sMax" % nm] = ft(100 * max(bb), 2).replace("−", "$-$")
    # diferenca entre as duas cadeias em h = 0,125 (mesmo h, malhas iniciais diferentes)
    dgc = [dados["g15_h0.125", c]["fator"] / dados["foco375_h0.125", c]["fator"] - 1. for c in casos]
    M["nivDgcMed"] = ft(100 * media(dgc), 2).replace("−", "$-$")
    M["nivDgcMax"] = ft(100 * max(abs(x) for x in dgc), 2)
    dgb = [dados["g15_h0.1875", c]["fator"] / dados["foco375_h0.1875", c]["fator"] - 1. for c in casos]
    M["nivDgbMax"] = ft(100 * max(abs(x) for x in dgb), 2)
    dga = [dados["g15_h0.25", c]["fator"] / dados["foco375_h0.25", c]["fator"] - 1. for c in casos]
    M["nivDgaMax"] = ft(100 * max(abs(x) for x in dga), 2)
    # derivada de P_f no vies: dP_f/P_f por 1 % de vies
    eps = 1e-4
    dpf = (pf_config(eps) - pf_config(-eps)) / (2 * eps) / pf_ref
    M["nivSensPf"] = ft(abs(dpf), 1)
    M["nivFermiFoco"] = ft(FERMI_SEG_REAL * 1000 / FERMI_SIMULT / 60., 1)
    with open(os.path.join(AQUI, "secoes", "refino_niveis_valores.tex"), "w") as fh:
        fh.write("%% gerado por analise_niveis.py: os numeros citados em refino_niveis.tex\n")
        for k in sorted(M):
            fh.write("\\def\\%s{%s}\n" % (k, M[k]))

    # ----- figuras
    ordem_cor = {"foco375_h0.25": CORES[0], "foco375_h0.1875": CORES[0], "foco375_h0.125": CORES[0],
                 "g15_h0.25": CORES[1], "g15_h0.1875": CORES[1], "g15_h0.125": CORES[1], "g15": CORES[1]}
    marc = {0.25: "o", 0.1875: "s", 0.125: "D", None: "^", 0.375: "v"}

    # niv_erro_tempo: |e| medio e maximo x CPU equivalente da cadeia
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.2), sharex=True)
    for ax, chave, tit in ((axs[0], "e_abs", r"$\overline{|e|}$ (%)"), (axs[1], "e_max", r"$\max|e|$ (%)")):
        for fam, cor, rot in (("fixa", CINZA, "malhas fixas"), ("foco375", CORES[0], "foco375 → h"),
                              ("g15", CORES[1], "g15 → h")):
            ks = [k for k, *_ in CONFIGS if FAM[k] == fam and k != "foco125"]
            xs = [resumo[k]["cpu_eq"] for k in ks]
            ys = [100 * resumo[k][chave] for k in ks]
            ax.plot(xs, ys, "-", color=cor, lw=1.2, alpha=0.8, zorder=2)
            for k, x, y in zip(ks, xs, ys):
                ax.plot([x], [y], marker=marc[HCFG[k]], color=cor, ms=5.5, ls="none", zorder=3,
                        mfc=cor if fam != "fixa" else "white")
            ax.plot([], [], "-", color=cor, label=rot)
        k = "g15"
        ax.plot([resumo[k]["cpu_eq"]], [100 * resumo[k][chave]], marker="^", color=CORES[1], mfc="white", ms=5.5)
        ax.axvline(resumo["foco125"]["cpu_eq"], color=CINZA, ls=":", lw=1)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_ylabel(tit)
        ax.set_xlabel("CPU da cadeia por caso (s, escala da campanha nova)")
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _p: virg(x, 0) if x >= 1 else virg(x, 1)))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _p: virg(y, 2 if y < 0.1 else (1 if y < 1 else 0))))
        ax.yaxis.set_minor_formatter(FuncFormatter(lambda y, _p: ""))
        for k, *_ in CONFIGS:
            if k == "foco125":
                continue
            x, y = resumo[k]["cpu_eq"], 100 * resumo[k][chave]
            txt = CURTO[k].replace("foco375→", "").replace("g15→", "") if FAM[k] in ("foco375", "g15") else CURTO[k]
            dy = {"foco375": -9, "g15": 5, "fixa": 5, "g15i": 5}[FAM[k]]
            ax.annotate(txt, (x, y), textcoords="offset points", xytext=(4, dy), fontsize=6.5,
                        color={"fixa": "#555555", "foco375": CORES[0], "g15": CORES[1], "g15i": CORES[1]}[FAM[k]])
        ax.text(resumo["foco125"]["cpu_eq"], ax.get_ylim()[1], " foco125", fontsize=6.5, color="#555555", va="top")
    axs[0].legend(loc="lower left", fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(AQUI, "figuras", "niv_erro_tempo.pdf"))
    plt.close(fig)

    # niv_convergencia: e (%) x h por caso, as duas cadeias e as fixas
    n = len(casos)
    nc = 3
    nl = int(math.ceil(n / nc))
    fig, axs = plt.subplots(nl, nc, figsize=(7.2, 2.0 * nl + 0.4), sharex=True)
    axs = np.atleast_1d(axs).ravel()
    for i, c in enumerate(casos):
        ax = axs[i]
        hs = [0.375, 0.25, 0.1875, 0.125]
        ax.plot(hs, [100 * dados[k, c]["e"] for k in ("foco375", "foco25", "foco1875", "foco125")], "-v",
                color=CINZA, mfc="white", ms=4, lw=1, label="malhas fixas")
        for nome, cor, rot in (("foco375", CORES[0], "foco375 → h"), ("g15", CORES[1], "g15 → h")):
            ks = cadeias[nome]
            ax.plot(H3, [100 * dados[k, c]["e"] for k in ks], "-o", color=cor, ms=3.5, lw=1.2, label=rot)
            g0 = ext[nome, c]["g0_2"]
            ax.plot([0.], [100 * (g0 / dados["foco125", c]["fator"] - 1.)], marker="o", mfc="white", color=cor,
                    ms=4, ls="none")
            ax.plot([0., 0.125], [100 * (g0 / dados["foco125", c]["fator"] - 1.), 100 * dados[ks[2], c]["e"]], ":",
                    color=cor, lw=0.9)
        g0 = ext["fixas", c]["g0_2"]
        ax.plot([0.], [100 * (g0 / dados["foco125", c]["fator"] - 1.)], marker="v", mfc="white", color=CINZA, ms=4)
        ax.plot([0., 0.125], [100 * (g0 / dados["foco125", c]["fator"] - 1.), 0.], ":", color=CINZA, lw=0.9)
        ax.axhline(0, color="#999999", lw=0.6)
        ax.set_title(ROT_CASO[c] + "   (Γ$_{125}$ = %s)" % virg(dados["foco125", c]["fator"], 3), fontsize=8)
        todos = [100 * dados[k, c]["e"] for k in cadeias["foco375"] + cadeias["g15"] + ["foco25", "foco1875"]] + \
            [100 * (ext[nm, c]["g0_2"] / dados["foco125", c]["fator"] - 1.) for nm in cadeias]
        lo, hi = min(todos + [-0.3]), max(todos + [0.3])
        mg = 0.15 * (hi - lo)
        ax.set_ylim(lo - mg, hi + mg)
        if dados["foco375", c]["e"] * 100 > hi + mg:
            ax.annotate("foco375: %s%%" % virg(100 * dados["foco375", c]["e"], 1), (0.375, hi + mg),
                        xytext=(0.30, hi + mg - 0.1 * (hi - lo)), fontsize=6, color="#555555", ha="center",
                        arrowprops=dict(arrowstyle="->", color="#999999", lw=0.6))
        ax.set_xlim(-0.015, 0.39)
        ax.set_xticks([0, 0.125, 0.1875, 0.25, 0.375])
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _p: virg(x, 4).rstrip("0").rstrip(",") if x else "0"))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _p: virg(y, 1)))
        ax.tick_params(labelsize=7)
        if i % nc == 0:
            ax.set_ylabel("Γ/Γ$_{125}$ − 1 (%)", fontsize=8)
        if i >= n - nc:
            ax.set_xlabel("h (m)", fontsize=8)
            ax.tick_params(axis="x", labelrotation=45)
    for j in range(n, len(axs)):
        axs[j].axis("off")
    h_, l_ = axs[0].get_legend_handles_labels()
    fig.legend(h_, l_, loc="lower center", ncol=3, fontsize=7.5, bbox_to_anchor=(0.5, 0.0))
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(os.path.join(AQUI, "figuras", "niv_convergencia.pdf"))
    plt.close(fig)

    # niv_etapas: barras horizontais empilhadas do CPU medio por etapa (escala da campanha nova)
    ks = [k for k, *_ in CONFIGS][::-1]
    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    ys = np.arange(len(ks))
    a = np.array([resumo[k]["cpu1_eq"] for k in ks])
    b = np.array([resumo[k]["cpu_malha_eq"] for k in ks])
    cc = np.array([resumo[k]["cpu_run_eq"] for k in ks])
    ax.barh(ys, a, color=[CINZA if FAM[k] == "fixa" else CORES[2] for k in ks], height=0.68)
    ax.barh(ys, b, left=a, color=CORES[3], height=0.68)
    ax.barh(ys, cc, left=a + b, color=[ordem_cor.get(k, CINZA) for k in ks], height=0.68)
    for y, k in zip(ys, ks):
        tot = resumo[k]["cpu_eq"]
        ax.text(tot + 3, y, virg(tot, 1), ha="left", va="center", fontsize=7)
    ax.set_yticks(ys)
    ax.set_yticklabels([CURTO[k] for k in ks], fontsize=7.5)
    ax.set_xlabel("CPU médio por caso (s, escala da campanha nova)")
    ax.set_xlim(0, max(resumo[k]["cpu_eq"] for k in ks) * 1.12)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _p: virg(x, 0)))
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=CINZA, label="malha fixa"), Patch(color=CORES[2], label="etapa 1 (malha inicial)"),
                       Patch(color=CORES[3], label="refina_banda.py + gmsh"),
                       Patch(color=CORES[0], label="rodada refinada (foco375 → h)"),
                       Patch(color=CORES[1], label="rodada refinada (g15 → h)")], fontsize=7, loc="lower right")
    ax.grid(axis="y", visible=False)
    fig.tight_layout()
    fig.savefig(os.path.join(AQUI, "figuras", "niv_etapas.pdf"))
    plt.close(fig)

    # niv_contencao: d_max por caso e configuracao refinada
    kr = [k for k, *_ in CONFIGS if FAM[k] in ("foco375", "g15")]
    fig, ax = plt.subplots(figsize=(7.2, 2.6))
    larg = 0.8 / len(kr)
    for j, k in enumerate(kr):
        ys = [dados[k, c].get("dist") for c in casos]
        cor = CORES[0] if FAM[k] == "foco375" else CORES[1]
        alfa = {0.25: 0.45, 0.1875: 0.7, 0.125: 1.0}[HCFG[k]]
        ax.bar(np.arange(len(casos)) - 0.4 + (j + 0.5) * larg, [y if y is not None else 0 for y in ys], width=larg,
               color=cor, alpha=alfa, label=CURTO[k])
    ax.axhline(FOLGA, color=CORES[7], ls="--", lw=1)
    ax.text(-0.45, FOLGA + 0.005, "folga 0,375 m", color=CORES[7], fontsize=7, ha="left", va="bottom")
    ax.set_xticks(np.arange(len(casos)))
    ax.set_xticklabels([ROT_CASO[c] for c in casos])
    ax.set_ylabel("$d_{\\max}$ (m)")
    ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _p: virg(y, 2)))
    ax.legend(fontsize=6.5, ncol=6, loc="upper center", bbox_to_anchor=(0.5, 1.13), handlelength=1.2,
              columnspacing=0.8)
    ax.set_ylim(0, max(FOLGA * 1.25, max(dados[k, c].get("dist") or 0 for k in kr for c in casos) * 1.35))
    ax.grid(axis="x", visible=False)
    fig.tight_layout()
    fig.savefig(os.path.join(AQUI, "figuras", "niv_contencao.pdf"))
    plt.close(fig)

    # ----- impressao
    print("\nconfig                 tri    |e|med  |e|max  e_med   CPU     CPUeq   parede  R      R/R125  dist   fora")
    for k, *_ in CONFIGS:
        x = resumo[k]
        print("%-20s %6.0f  %6.3f  %6.3f  %+6.3f  %7.1f %7.1f %7.1f  %6.2f %7.4f  %5s  %5s" % (
            k, x["triangulos"], 100 * x["e_abs"], 100 * x["e_max"], 100 * x["e_med"], x["cpu"], x["cpu_eq"],
            x["parede"], x["R"], x["R125"], "--" if x["dist"] is None else "%.3f" % x["dist"],
            "--" if x["fora"] is None else "%.3f" % (100 * x["fora"])))
    print("\nCPU etapa 1 foco375: nova %.2f s, comparacao %.2f s, razao %.3f (%.3f a %.3f)" % (
        cpu375["nova"], cpu375["comp"], cpu375["nova"] / cpu375["comp"], min(razao_etapa1), max(razao_etapa1)))
    print("parede etapa 1 foco375: nova %.2f s, comparacao %.2f s" % (par375["nova"], par375["comp"]))
    for nome in cadeias:
        print("p %s: %s ; Gamma0(p=2)/foco125-1 %%: %s" % (
            nome, ", ".join("%s:%s" % (ROT_CASO[c], "%.2f" % ext[nome, c]["p"] if ext[nome, c]["p"] else ext[nome, c]["sit"])
                            for c in casos), ", ".join("%.3f" % (100 * x) for x in b0[nome])))
    print("lognormal: mu_ln %.4f sigma_ln %.4f P_f(ln, foco375) %.4f P_f(ref) %.4f ; CV N=1000 %.4f N=5000 %.4f" % (
        m_ln, s_ln, pf_ln375, pf_ref, cv(pf_ref, 1000), cv(pf_ref, 5000)))
    for k, bx, pf, dpf_, rel in pf_rows:
        print("  %-18s b %+.4f%%  P_f %.4f%%  dP_f %+.4f pp  rel %+.2f%%  fermi 1000: %.2f h 5000: %.2f h" % (
            k, 100 * bx, 100 * pf, 100 * dpf_, 100 * rel, horas(resumo[k]["R"], 1000), horas(resumo[k]["R"], 5000)))
    print("g15 vs foco375 em h=0,125: media %.3f%% max %.3f%%; 0,1875 max %.3f%%; 0,25 max %.3f%%" % (
        100 * media(dgc), 100 * max(abs(x) for x in dgc), 100 * max(abs(x) for x in dgb), 100 * max(abs(x) for x in dga)))


if __name__ == "__main__":
    main()
