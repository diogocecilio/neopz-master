"""figuras_gstools.py: figuras da comparacao entre a KL em forma fechada do retangulo (a do SlopeRandomField) e o
metodo de randomizacao (Hesse et al. 2014, o do GSTools) no talude de Vargas Ceron et al. (2025):
  kl_gs_campos.png      uma realizacao de cada gerador no talude, e a covariancia de cada um a partir de um ponto;
  kl_gs_covariancia.pdf a covariancia contra a distancia em x, em y e na diagonal (alvo, KL, randomizacao separavel,
                        GSTools Exponential), e o erro da covariancia contra o numero de realizacoes;
  kl_gs_variancia.png   a variancia pontual (exata para a KL, empirica para a randomizacao).
Geradores (ver compara_gstools.py): KL (M = 2000), RM-sep (randomizacao com o nucleo do Ceron, espectro produto de
duas Cauchy), GS-exp (gs.Exponential, len_scale = [20, 2]: distancia euclidiana escalada, outro nucleo).
Uso: python3 figuras_gstools.py [R=2000]
"""
import sys

import numpy as np

import klfechada as kf
from estilo import DIVERGENTE, FIGURAS, LARGURA, PALETA, plt, num
from matplotlib.ticker import FuncFormatter

R = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
N = 1000
rng = np.random.default_rng(7)
virg = FuncFormatter(lambda v, _: num(round(v, 6)))

h = 0.25
xs, ys = np.arange(h / 2, 25, h), np.arange(h / 2, 10, h)
xg, yg = np.meshgrid(xs, ys)
dentro = kf.no_talude(xg, yg)
X, Y = xg[dentro], yg[dentro]
n = len(X)

kl = kf.KL2D(2000)
Phi = (kl.fx(X)[kl.ix] * kl.gy(Y)[kl.iy]) * np.sqrt(kl.lam)[:, None]   # (M, n)

import gstools as gs  # noqa: E402
srf = gs.SRF(gs.Exponential(dim=2, var=1., len_scale=[20., 2.]), mode_no=N)


def rm_sep():
    kx, ky = rng.standard_cauchy(N) / 20., rng.standard_cauchy(N) / 2.
    ph = np.outer(kx, X) + np.outer(ky, Y)
    return np.sqrt(1. / N) * (rng.standard_normal(N) @ np.cos(ph) + rng.standard_normal(N) @ np.sin(ph))


def gs_exp(s):
    return srf((X, Y), seed=s)


def mapa(v):
    m = np.full(xg.shape, np.nan)
    m[dentro] = v
    return m


def contorno(ax):
    t = np.vstack([kf.TALUDE, kf.TALUDE[:1]])
    ax.plot(t[:, 0], t[:, 1], color="#0b0b0b", lw=0.8)
    ax.set_aspect("equal")
    ax.set_xlim(-0.3, 25.3)
    ax.set_ylim(-0.3, 10.3)
    ax.xaxis.set_major_formatter(virg)
    ax.yaxis.set_major_formatter(virg)
    ax.grid(False)


# ---- ponto de referencia e linhas de pares ----
p0 = np.argmin((X - 5.125) ** 2 + (Y - 2.625) ** 2)
lx = np.where(np.abs(Y - Y[p0]) < 1e-9)[0]; lx = lx[X[lx] >= X[p0]]
ly = np.where(np.abs(X - X[p0]) < 1e-9)[0]; ly = ly[Y[ly] >= Y[p0]]
ld = np.where(np.abs((X - X[p0]) / 20 - (Y - Y[p0]) / 2 * 0.1 * 1) < 1e-9)[0]  # placeholder, refeito abaixo
# diagonal em distancia escalada: dx/20 = dy/2  ->  dy = dx/10; na grade, passos de 2,5 m em x e 0,25 m em y
ld = np.array([np.argmin((X - (X[p0] + 2.5 * i)) ** 2 + (Y - (Y[p0] + 0.25 * i)) ** 2) for i in range(0, 8)])
linhas = {"x": lx, "y": ly, "d": ld}
pares_a = np.concatenate([np.full(len(v), p0) for v in linhas.values()])
pares_b = np.concatenate(list(linhas.values()))
# pares aleatorios para o erro contra R
pa, pb = rng.integers(0, n, 600), rng.integers(0, n, 600)
alvo_ab = np.exp(-np.abs(X[pa] - X[pb]) / 20 - np.abs(Y[pa] - Y[pb]) / 2)

# ---- ensembles ----
marcos = np.unique(np.round(np.logspace(1, np.log10(R), 12)).astype(int))
res = {}
for nome, ger in (("RM-sep", lambda s: rm_sep()), ("GS-exp", gs_exp)):
    Sl = np.zeros(len(pares_b)); S2 = np.zeros(n); Sab = np.zeros(len(pa)); Sp0 = np.zeros(n); err = []
    exemplo = None
    for s in range(R):
        u = ger(s)
        if exemplo is None:
            exemplo = u
        Sl += u[pares_a] * u[pares_b]; S2 += u ** 2; Sab += u[pa] * u[pb]; Sp0 += u[p0] * u
        if s + 1 in marcos:
            err.append(np.abs(Sab / (s + 1) - alvo_ab).mean())
    res[nome] = dict(lin=Sl / R, var=S2 / R, cov0=Sp0 / R, err=np.array(err), ex=exemplo)
    print(nome, "pronto")
xi = rng.standard_normal(2000)
z_kl = xi @ Phi
var_kl = (Phi ** 2).sum(0)
cov0_kl = Phi[:, p0] @ Phi
lin_kl = np.einsum("kp,kp->p", Phi[:, pares_a], Phi[:, pares_b])
err_kl = np.abs(np.einsum("kp,kp->p", Phi[:, pa], Phi[:, pb]) - alvo_ab).mean()
alvo0 = np.exp(-np.abs(X - X[p0]) / 20 - np.abs(Y - Y[p0]) / 2)

# ---- figura 1: campos e covariancia a partir de p0 ----
fig, axs = plt.subplots(3, 2, figsize=(LARGURA, 5.4), sharex=True, sharey=True, layout="constrained")
nomes = [("KL forma fechada, M = 2000", z_kl, cov0_kl), ("Randomização, núcleo do Ceron", res["RM-sep"]["ex"],
                                                            res["RM-sep"]["cov0"]),
         ("GSTools Exponential [20, 2]", res["GS-exp"]["ex"], res["GS-exp"]["cov0"])]
for i, (t, z, c0) in enumerate(nomes):
    a = axs[i, 0]
    im = a.pcolormesh(xs, ys, mapa(z), cmap="viridis", vmin=-2.5, vmax=2.5, shading="auto", rasterized=True)
    contorno(a); a.set_title("%s\numa realização" % t, loc="left", fontsize=7.5)
    b = axs[i, 1]
    im2 = b.pcolormesh(xs, ys, mapa(c0 - alvo0), cmap=DIVERGENTE, vmin=-0.15, vmax=0.15, shading="auto",
                       rasterized=True)
    b.plot(X[p0], Y[p0], "k+", ms=6)
    contorno(b); b.set_title("covariância com o ponto +, menos o alvo\n" + ("exata" if i == 0 else "%d realizações" % R),
                             loc="left", fontsize=7.5)
for a in axs[-1]:
    a.set_xlabel("x (m)")
for a in axs[:, 0]:
    a.set_ylabel("y (m)")
fig.colorbar(im, ax=axs[:, 0], shrink=0.5, location="bottom", label="Z (normal padrão)", format=virg)
fig.colorbar(im2, ax=axs[:, 1], shrink=0.5, location="bottom", label="C − C_alvo", format=virg)
fig.savefig(FIGURAS + "/kl_gs_campos.png", bbox_inches="tight")
plt.close(fig)

# ---- figura 2: covariancia contra a distancia e erro contra R ----
fig, axs = plt.subplots(1, 4, figsize=(LARGURA, 2.4))
pos = 0
dd = {"x": np.abs(X[lx] - X[p0]), "y": np.abs(Y[ly] - Y[p0]), "d": np.abs(X[ld] - X[p0])}
tit = {"x": "ao longo de x", "y": "ao longo de y", "d": "diagonal (dy = dx/10)"}
for j, k in enumerate(("x", "y", "d")):
    a = axs[j]; m = len(linhas[k]); sl = slice(pos, pos + m); pos += m
    d = dd[k]; o = np.argsort(d)
    hx = np.abs(X[linhas[k]] - X[p0]); hy = np.abs(Y[linhas[k]] - Y[p0])
    a.plot(d[o], np.exp(-hx / 20 - hy / 2)[o], color="#0b0b0b", lw=2.2, label="alvo (Ceron)")
    a.plot(d[o], np.exp(-np.sqrt((hx / 20) ** 2 + (hy / 2) ** 2))[o], color=PALETA[1], lw=1, ls="--",
           label="núcleo do GSTools")
    a.plot(d[o], lin_kl[sl][o], color=PALETA[0], lw=1.2, label="KL (exata)")
    a.plot(d[o], res["RM-sep"]["lin"][sl][o], "o", color=PALETA[2], ms=2.5, label="randomização Ceron")
    a.plot(d[o], res["GS-exp"]["lin"][sl][o], "s", color=PALETA[1], ms=2.5, label="GSTools")
    a.set_title(tit[k], loc="left"); a.set_xlabel("dx (m)" if k != "y" else "dy (m)")
    a.xaxis.set_major_formatter(virg); a.yaxis.set_major_formatter(virg)
axs[0].set_ylabel("covariância")
axs[0].legend(loc="lower left", fontsize=6)
a = axs[3]
a.loglog(marcos, res["RM-sep"]["err"], "o-", color=PALETA[2], ms=3, label="randomização Ceron")
a.loglog(marcos, res["GS-exp"]["err"], "s-", color=PALETA[1], ms=3, label="GSTools")
a.axhline(err_kl, color=PALETA[0], lw=1.2, label="KL (exata)")
a.loglog(marcos, 0.8 / np.sqrt(marcos), ":", color="#8a8984", lw=1, label="∝ 1/√R")
a.set_xlabel("realizações R"); a.set_ylabel("erro médio da covariância"); a.set_title("600 pares", loc="left")
a.legend(fontsize=6)
fig.tight_layout()
fig.savefig(FIGURAS + "/kl_gs_covariancia.pdf", bbox_inches="tight")
plt.close(fig)

# ---- figura 3: variancia pontual ----
fig, axs = plt.subplots(1, 3, figsize=(LARGURA, 2.2), sharey=True, layout="constrained")
for a, (t, v) in zip(axs, (("KL (exata)", var_kl), ("randomização Ceron, %d real." % R, res["RM-sep"]["var"]),
                            ("GSTools, %d real." % R, res["GS-exp"]["var"]))):
    im = a.pcolormesh(xs, ys, mapa(v), cmap=DIVERGENTE, vmin=0.85, vmax=1.15, shading="auto", rasterized=True)
    contorno(a); a.set_title("%s\nde %s a %s" % (t, num(v.min(), 3), num(v.max(), 3)), loc="left", fontsize=7)
    a.set_xlabel("x (m)")
axs[0].set_ylabel("y (m)")
fig.colorbar(im, ax=axs, shrink=0.9, label="variância", format=virg)
fig.savefig(FIGURAS + "/kl_gs_variancia.png", bbox_inches="tight")
plt.close(fig)
print("erro medio KL %.5f; RM-sep(R) %s; GS-exp(R) %s" % (err_kl, np.round(res["RM-sep"]["err"], 4),
                                                         np.round(res["GS-exp"]["err"], 4)))
