#!/usr/bin/env python3
"""figuras_ls.py: as figuras do line sampling da formulacao (figs_ls/ls_esquema.pdf e figs_ls/ls_exemplos.pdf), em
Python. Os exemplos sao os mesmos do exemplo_ls.wl (Mathematica), que da os numeros do texto; aqui os sorteios sao outros,
entao as curvas mudam um pouco, mas nao o que elas mostram.
Uso: python3 -I figuras_ls.py   (~1 min). Tambem copia as duas para Projects/PlanoDeTrabalho/comparacao_metodos/figs."""
import os
import shutil

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Circle, FancyArrowPatch, Polygon  # noqa: E402
from scipy.special import ndtr  # noqa: E402
from scipy.stats import norm  # noqa: E402

AQUI = os.path.dirname(os.path.abspath(__file__))
FIGS = os.path.join(AQUI, "figs_ls")
os.makedirs(FIGS, exist_ok=True)
plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10.5, "legend.fontsize": 8.5,
                     "xtick.labelsize": 9, "ytick.labelsize": 9, "mathtext.fontset": "cm"})
LARG = 7.2
AZUL, VERM, VERDE, LAR = "#2166ac", "#c0392b", "#1a7f37", "#e08214"


def cauda(c):
    return ndtr(-np.asarray(c, dtype=float))


# ------------------------------------------------------------------ 1. o esquema
th = np.radians(55.)
beta, kap = 2.0, 0.13
e = np.array([np.cos(th), np.sin(th)])            # a direcao (para a falha)
t = np.array([-np.sin(th), np.cos(th)])           # o plano ortogonal


def ponto(s, r):
    return np.outer(np.atleast_1d(s), t).squeeze() + np.outer(np.atleast_1d(r), e).squeeze()


nu = lambda s: beta + kap * s ** 2                # a superficie de falha: nu = beta + kap s^2
ss = np.array([-2.3, -1.25, -0.35, 0.55, 1.45, 2.3])
rr = np.array([0.7, -0.9, 1.1, -0.35, 0.45, -1.3])
zs = np.array([ponto(s, r) for s, r in zip(ss, rr)])
i1 = 4                                            # a linha com os rotulos

fig, (a, b) = plt.subplots(1, 2, figsize=(LARG, 3.75), gridspec_kw={"width_ratios": [1.05, 1]})
sgrid = np.linspace(-6, 6, 400)
sup = np.array([ponto(s, nu(s)) for s in sgrid])
a.add_patch(Polygon(np.vstack([sup, ponto(6, 12), ponto(-6, 12)]), closed=True, color="#fbe3e0", lw=0))
a.plot(sup[:, 0], sup[:, 1], color=VERM, lw=2)
for r_ in (1, 2, 3):
    a.add_patch(Circle((0, 0), r_, fill=False, color="0.8", lw=0.8))
pl = np.array([ponto(-5, 0), ponto(5, 0)])
a.plot(pl[:, 0], pl[:, 1], "k--", lw=0.9)
for i, (s, z) in enumerate(zip(ss, zs)):
    zp = s * t
    xi = ponto(s, nu(s))
    ini = zp - 1.2 * e
    a.plot([ini[0], xi[0]], [ini[1], xi[1]], color="0.6", lw=0.7)
    fim = zp + 6 * e
    a.plot([xi[0], fim[0]], [xi[1], fim[1]], color=VERM, lw=2.2)
    a.plot([z[0], zp[0]], [z[1], zp[1]], ":", color=AZUL, lw=1)
    a.plot(*z, "o", color=AZUL, ms=5)
    a.plot(*zp, "o", mfc="white", mec=AZUL, ms=5.5)
    a.plot(*xi, "o", color=VERM, ms=5.5)
a.add_patch(FancyArrowPatch((0, 0), tuple(beta * e), arrowstyle="-|>", mutation_scale=14, color=VERDE, lw=2.2))
a.plot(0, 0, "ko", ms=5)
a.plot(*(beta * e), "ko", ms=5)
a.annotate(r"$\xi^* = \beta\,\mathbf{e}$", beta * e, xytext=(10, -14), textcoords="offset points", color=VERDE, fontsize=10.5,
           bbox=dict(fc="white", ec="none", pad=0.5))
a.text(*(0.95 * e + 0.28 * t), r"$\mathbf{e}$", color=VERDE, fontsize=12)
a.text(*(ponto(-3.0, -0.32)), "plano ortogonal a $\\mathbf{e}$", fontsize=8.5, style="italic", rotation=np.degrees(th) - 90,
       ha="center", va="center", bbox=dict(fc="white", ec="none", pad=0.3))
a.text(2.2, 3.7, "falha\n$g(\\xi) \\leq 0$", color=VERM, fontsize=10, ha="center", bbox=dict(fc="white", ec="none", pad=1))
z1, zp1, x1 = zs[i1], ss[i1] * t, ponto(ss[i1], nu(ss[i1]))
a.annotate(r"$\xi_\ell \sim N(\mathbf{0}, I)$", z1, xytext=(-62, 14), textcoords="offset points", color=AZUL, fontsize=9.5,
           bbox=dict(fc="white", ec="none", pad=0.3), arrowprops=dict(arrowstyle="-", color=AZUL, lw=0.6))
a.annotate(r"$\xi^\perp_\ell$", zp1, xytext=(-20, -12), textcoords="offset points", color=AZUL, fontsize=10.5)
a.annotate(r"$\nu_\ell$", x1, xytext=(-18, 3), textcoords="offset points", color=VERM, fontsize=11)
a.set_xlim(-3.6, 3.6)
a.set_ylim(-3.2, 4.3)
a.set_aspect("equal")
a.set_xlabel(r"$\xi_1$")
a.set_ylabel(r"$\xi_2$")
a.set_title("(a) cada amostra vira uma linha paralela a $\\mathbf{e}$", fontsize=10)

c = np.linspace(-4, 5.6, 500)
cl = nu(ss[i1])
b.plot(c, norm.pdf(c), "k", lw=1.8)
cc = c[c >= cl]
b.fill_between(cc, norm.pdf(cc), color=VERM, alpha=0.45, lw=0)
b.axvline(beta, color=VERDE, ls="--", lw=1)
b.plot([cl, cl], [0, 0.33], color=VERM, lw=2)
b.text(beta, 0.432, r"$\beta$ (FORM)", color=VERDE, ha="center", fontsize=10)
b.text(cl + 0.1, 0.34, r"$\nu_\ell$", color=VERM, fontsize=11)
b.annotate(r"$p_\ell = \Phi(-\nu_\ell)$", (cl + 0.35, 0.012), xytext=(cl + 0.4, 0.11), color=VERM, fontsize=10.5,
           arrowprops=dict(arrowstyle="-", color=VERM, lw=0.8))
b.text(-9.6, 0.41, "ao longo da linha,\n$\\nu \\sim N(0,1)$", fontsize=9.5, style="italic", va="top")
b.text(-9.6, 0.27, r"$\hat P_f = \frac{1}{N}\sum_{\ell=1}^{N} \Phi(-\nu_\ell)$", fontsize=11.5, va="center")
b.text(-9.6, 0.13, "só a variação dos $\\nu_\\ell$\nentre as linhas entra\nna variância", fontsize=9, va="center")
b.set_ylim(0, 0.46)
b.set_xlim(-9.8, 5.6)
b.set_xticks([-2, 0, 2, 4])
b.set_xlabel(r"$\nu$ (coordenada ao longo de $\mathbf{e}$)")
b.set_ylabel("densidade")
b.set_title("(b) cada linha dá a sua probabilidade exata", fontsize=10)
fig.tight_layout()
fig.savefig(os.path.join(FIGS, "ls_esquema.pdf"))
plt.close(fig)

# ------------------------------------------------------------------ 2. os exemplos
rng = np.random.default_rng(7)
# (a)(b) superficie curva: falha xi2 > 2,5 + 0,15 xi1^2; LS com e = (0, 1) (2 analises por linha) x AI em N(xi*, I)
b2, k2 = 2.5, 0.15
y = np.linspace(-12, 12, 20001)
pfx = np.trapezoid(norm.pdf(y) * cauda(b2 + k2 * y ** 2), y)
pls = cauda(b2 + k2 * rng.standard_normal(500) ** 2)
xa = rng.standard_normal((1000, 2)) + np.array([0., b2])
iai = (xa[:, 1] > b2 + k2 * xa[:, 0] ** 2) * np.exp(-b2 * xa[:, 1] + b2 ** 2 / 2)


def acum(v, custo):
    k = np.arange(1, len(v) + 1)
    m = np.cumsum(v) / k
    var = (np.cumsum(v ** 2) - k * m ** 2) / np.maximum(k - 1, 1)
    cv = np.where((k > 4) & (m > 0), np.sqrt(np.maximum(var, 0) / k) / np.where(m > 0, m, 1), np.nan)
    return custo * k, m, cv


fig, ax = plt.subplots(2, 2, figsize=(LARG, 6.9))
xl, ml, cl_ = acum(pls, 2)
xi_, mi, ci = acum(iai, 1)
A = ax[0, 0]
A.semilogx(xl, 1e3 * ml, color=AZUL, lw=1.6, label="line sampling (2 análises por linha)")
A.semilogx(xi_, 1e3 * mi, color=LAR, lw=1.4, label="AI no ponto de projeto")
A.axhline(1e3 * pfx, color="k", ls="--", lw=1, label="exata")
A.set_ylim(0, 2.5e3 * pfx)
A.set_xlabel("análises")
A.set_ylabel(r"$\hat P_f \times 10^3$")
A.set_title("(a) superfície curva: a estimativa", fontsize=10)
A.legend(loc="upper right")
A = ax[0, 1]
A.loglog(xl, cl_, color=AZUL, lw=1.6, label="line sampling")
A.loglog(xi_, ci, color=LAR, lw=1.4, label="AI no ponto de projeto")
A.axhline(0.05, color="k", ls=":", lw=1)
A.set_ylim(0.01, 1)
A.set_xlabel("análises")
A.set_ylabel("CoV estimado")
A.set_title("(b) o CoV (pontilhado: 5%)", fontsize=10)
A.legend(loc="upper right")

# (c) n = 2000, g = 2,5 - a.xi, direcao truncada em m termos: o CoV por linha exato
n5 = 2000
j = np.arange(1, n5 + 1)
a5 = np.sqrt(np.exp(-j / 60.) + 0.02 * np.exp(-j / 600.))
a5 /= np.linalg.norm(a5)
v = np.linspace(-9, 9, 6001)
ms_, cvs, vis = [], [], []
for m in (5, 10, 20, 50, 100, 200, 500, 1000):
    sa = np.linalg.norm(a5[:m])                    # inclinacao ao longo de e_m (|a| = 1)
    sr = np.sqrt(max(1 - sa ** 2, 0.))
    p = cauda((b2 - sr * v) / sa)
    m1 = np.trapezoid(norm.pdf(v) * p, v)
    m2 = np.trapezoid(norm.pdf(v) * p ** 2, v)
    ms_.append(m)
    cvs.append(np.sqrt(m2 - m1 ** 2) / m1)
    vis.append(sa ** 2)
A = ax[1, 0]
A.loglog(ms_, cvs, "o-", color=AZUL, lw=1.6, ms=5)
for m, cv, f in zip(ms_, cvs, vis):
    A.annotate(f"{100 * f:.0f}%", (m, cv), xytext=(4, 5), textcoords="offset points", fontsize=8, color="0.35")
A.set_xlim(3, 2500)
A.set_xlabel("m (termos da direção do FORM)")
A.set_ylabel("CoV por linha")
A.set_title("(c) n = 2000, direção truncada em m termos\n(rótulos: fração de $|\\mathbf{a}|^2$ vista; m = 2000: CoV = 0)",
            fontsize=9.5)

# (d) com o tempo: dois episodios (12 e 36 h) em direcoes diferentes; LS com a direcao de 36 h
tt = np.arange(0, 49)
bT = 3.0 - 0.55 * np.exp(-((tt - 12) / 5.) ** 2) - 0.7 * np.exp(-((tt - 36) / 5.) ** 2)
thT = np.radians(15 + 70 / (1 + np.exp(-(tt - 24) / 2.)))
D = np.stack([np.cos(thT), np.sin(thT)], axis=1)  # a normal de cada hora
e6 = D[36]
z = rng.standard_normal((1000, 2))
zp = z - np.outer(z @ e6, e6)
nut = (bT[None, :] - zp @ D.T) / (D @ e6)[None, :]  # a raiz de cada hora em cada linha
pfT = cauda(np.minimum.accumulate(nut, axis=1)).mean(axis=0)
pinst = cauda(nut).mean(axis=0)
xm = rng.standard_normal((200000, 2))
gm = bT[None, :] - xm @ D.T
pfMC = (np.minimum.accumulate(gm, axis=1) < 0).mean(axis=0)
A = ax[1, 1]
A.semilogy(tt, pfMC, "k", lw=1.8, label="acumulada até T, MC de 200 000")
A.semilogy(tt[::2], pfT[::2], "s", color=AZUL, ms=4.5, label="acumulada, line sampling (1000 linhas)")
A.semilogy(tt, cauda(bT), "--", color="0.5", lw=1.3, label="instantânea, exata")
A.semilogy(tt[::2], pinst[::2], "o", mfc="none", mec=VERM, ms=6, label="instantânea, line sampling")
A.set_ylim(1e-4, 0.1)
A.set_xlabel("T (h)")
A.set_ylabel("probabilidade de falha")
A.set_title("(d) com o tempo: $P_f(T)$ das mesmas linhas", fontsize=10)
A.legend(loc="lower right", fontsize=7.8)
fig.tight_layout()
fig.savefig(os.path.join(FIGS, "ls_exemplos.pdf"))
plt.close(fig)

dest = os.path.join(AQUI, "..", "..", "..", "Projects", "PlanoDeTrabalho", "comparacao_metodos", "figs")
if os.path.isdir(dest):
    for f in ("ls_esquema.pdf", "ls_exemplos.pdf"):
        shutil.copy(os.path.join(FIGS, f), dest)
print(f"pronto: figs_ls/ls_esquema.pdf, figs_ls/ls_exemplos.pdf (P_f exata da parte 2 = {pfx:.4e}; "
      f"P_f(48 h): LS {pfT[-1]:.4f}, MC {pfMC[-1]:.4f})")
