#!/usr/bin/env python3
"""gera_dados.py: as figuras e tabelas do relatorio do Cho2020 (relatorio_cho2020.tex), com os dados de
~/projects/runs/cho2020/passo1 (eps2_p<p>_iii_fs.csv e _perfis.csv; o que faltar fica de fora):
  figs/fig3.pdf      (a) a frente de molhamento e (b) o FS contra o tempo, p = 0,2 a 1,2 k_s, com o artigo;
  figs/fig4.pdf      p/k_s contra o tempo ate a ruptura (FS < 1), NeoPZ e artigo (Fig. 4);
  figs/perfis.pdf    os perfis de p_w e S_w em p = 0,2 k_s;
  tab/tempos.tex     o tempo em que a frente chega a rocha e o FS antes e depois, NeoPZ e artigo."""
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

AQUI = os.path.dirname(os.path.abspath(__file__))
R = os.path.expanduser("~/projects/runs/cho2020/passo1")
TINTA, TINTA2, GRADE = "#0b0b0b", "#5c5c58", "#e4e3df"
PS = [0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0, 1.2]
COR = dict(zip(PS, plt.cm.viridis([i / (len(PS) - 1) for i in range(len(PS))])))
os.makedirs(os.path.join(AQUI, "figs"), exist_ok=True)
os.makedirs(os.path.join(AQUI, "tab"), exist_ok=True)
v = lambda x, d=2: f"{x:.{d}f}".replace(".", "{,}")

lido = [r for r in csv.reader(open(os.path.join(AQUI, "..", "data", "fig3_fig4_lido.csv"))) if r and not r[0].startswith("#")]
lido = {float(r[0]): float(r[1]) for r in lido[1:]}


def le(p):
    arq = os.path.join(R, f"eps2_p{p:.2f}_iii_fs.csv")
    return list(csv.DictReader(open(arq))) if os.path.exists(arq) else []


def eixo(ax):
    ax.grid(True, color=GRADE, lw=0.7)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(labelsize=8, colors=TINTA2)


def t_ruptura(Rr):
    for r in Rr:
        if float(r["FS"]) < 1.:
            return float(r["t_h"])
    return None


# ---- Fig. 3
fig, (a, b) = plt.subplots(1, 2, figsize=(7.6, 3.6))
for ax in (a, b):
    eixo(ax)
for p in PS:
    Rr = le(p)
    if not Rr or p in (0.1, 0.3):
        continue
    t = [float(r["t_h"]) for r in Rr]
    a.plot(t, [float(r["frente_m"]) for r in Rr], color=COR[p], lw=1.6, label=f"p = {p:g} $k_s$")
    b.plot(t, [float(r["FS"]) for r in Rr], color=COR[p], lw=1.6)
    if p in lido:
        a.plot([lido[p]], [2.], "o", ms=6, mfc="none", mec=COR[p], mew=1.3)
        b.plot([lido[p], lido[p]], [1.50, 0.92], ls="none", marker="o", ms=5, mfc="none", mec=COR[p], mew=1.2)
a.plot([], [], "o", ms=6, mfc="none", mec=TINTA2, label="Cho (2020), Fig. 3")
a.set_xlim(0, 16), a.set_ylim(0, 2.5), b.set_xlim(0, 16), b.set_ylim(0.8, 1.6)
a.set_xlabel("t (h)", fontsize=9), b.set_xlabel("t (h)", fontsize=9)
a.set_ylabel("frente de molhamento, profundidade vertical (m)", fontsize=8)
b.set_ylabel("fator de segurança (Eq. 21)", fontsize=8)
a.set_title("(a)", fontsize=9, loc="left"), b.set_title("(b)", fontsize=9, loc="left")
a.legend(fontsize=7, frameon=False, loc="lower right")
fig.tight_layout()
fig.savefig(os.path.join(AQUI, "figs", "fig3.pdf"))
plt.close(fig)

# ---- Fig. 4
fig, ax = plt.subplots(figsize=(5.2, 3.8))
eixo(ax)
ax.set_xscale("log"), ax.set_yscale("log")
nos = [(t_ruptura(le(p)), p) for p in PS if le(p)]
nos = [(t, p) for t, p in nos if t is not None]
ax.plot([t for t, _ in nos], [p for _, p in nos], "-", color="#2a78d6", lw=1.6, marker="o", ms=4,
        label="NeoPZ, coluna normal")
ax.plot(list(lido.values()), list(lido.keys()), "s", ms=5, mfc="none", mec="#e34948", mew=1.2,
        label="Cho (2020), Figs. 3 e 4 (lido)")
ax.set_xlim(1, 100), ax.set_ylim(0.08, 1.5)
ax.set_xlabel("tempo até a ruptura (h)", fontsize=9)
ax.set_ylabel("$p/k_s$", fontsize=9)
ax.legend(fontsize=7, frameon=False)
fig.tight_layout()
fig.savefig(os.path.join(AQUI, "figs", "fig4.pdf"))
plt.close(fig)

# ---- perfis em p = 0,2 k_s
arq = os.path.join(R, "eps2_p0.20_iii_perfis.csv")
if os.path.exists(arq):
    P = list(csv.DictReader(open(arq)))
    horas = [0, 3, 6, 9, 12, 13, 14]
    cor = plt.cm.viridis([i / (len(horas) - 1) for i in range(len(horas))])
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.0, 3.8), sharey=True)
    for ax in (a, b):
        eixo(ax)
    for h, c in zip(horas, cor):
        Q = [r for r in P if abs(float(r["t_h"]) - h) < 1e-6]
        if not Q:
            continue
        z = [float(r["z_m"]) for r in Q]
        a.plot([float(r["pw_kPa"]) for r in Q], z, color=c, lw=1.4, label=f"{h} h")
        b.plot([float(r["Sw"]) for r in Q], z, color=c, lw=1.4)
    a.set_ylim(2, 0)
    a.set_xlabel("$p_w$ (kPa)", fontsize=9), b.set_xlabel("$S_w$", fontsize=9)
    a.set_ylabel("profundidade vertical (m)", fontsize=9)
    a.legend(fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(AQUI, "figs", "perfis.pdf"))
    plt.close(fig)

# ---- tabela
with open(os.path.join(AQUI, "tab", "tempos.tex"), "w") as f:
    f.write("\\begin{tabular}{c|cc|cc|c}\n\\toprule\n")
    f.write("$p/k_s$ & \\multicolumn{2}{c|}{frente na rocha (h)} & FS inicial & FS final & passos\\\\\n")
    f.write(" & NeoPZ & Cho & \\multicolumn{2}{c|}{NeoPZ (Cho: 1,50 e 0,92)} & \\\\\n\\midrule\n")
    for p in PS:
        Rr = le(p)
        if not Rr:
            f.write(f"{v(p, 1)} & (rodando) & {v(lido[p]) if p in lido else '--'} & & & \\\\\n")
            continue
        tc = next((float(r["t_h"]) for r in Rr if float(r["frente_m"]) >= 2. - 1e-9), None)
        f.write(f"{v(p, 1)} & {v(tc) if tc is not None else '--'} & {v(lido[p]) if p in lido else '--'} & "
                f"{v(float(Rr[0]['FS']), 3)} & {v(float(Rr[-1]['FS']), 3)} & {Rr[-1]['passos']}\\\\\n")
    f.write("\\bottomrule\n\\end{tabular}\n")
print("pronto")
