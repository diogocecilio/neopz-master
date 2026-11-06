#!/usr/bin/env python3
"""plot_realizacao.py <pasta da realizacao> [t_h] [saida.png]: uma realizacao em 6 paineis, a partir dos CSV dela:
FS(t) (fs.csv); chuva, infiltracao e fluxos laterais por hora (passos.csv); infiltracao e armazenamento acumulados
(passos.csv); p_w no instante t_h (padrao: o do FS minimo), com o lencol (p_w = 0), a face de percolacao e a banda de
ruptura (pw.csv, triangulos.csv, mecanismo.csv); c' com a banda (campos.csv); e o mecanismo, sqrt(J2) da deformacao
plastica no ultimo estado convergido do SRM (mecanismo.csv). Requer numpy e matplotlib."""
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np

R = sys.argv[1]
args = sys.argv[2:]
saida = next((a for a in args if a.endswith(".png")), os.path.join(R, "realizacao.png"))
th = next((float(a) for a in args if not a.endswith(".png")), None)
# cores categoricas em ordem fixa (paleta validada do guia de visualizacao) e tinta neutra para texto e referencias
AZUL, LARANJA, AQUA, AMARELO, VERDE, VIOLETA = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#008300", "#6250d6"
TINTA, CINZA = "#0b0b0b", "#8a8984"


def le(nome):
    arq = os.path.join(R, nome)
    if not os.path.exists(arq) or os.path.getsize(arq) == 0:
        return None
    d = np.genfromtxt(arq, delimiter=",", names=True, dtype=None, encoding=None)
    return np.atleast_1d(d) if d.size else None


fs, pas, pw, tri, mec, cam = (le(n) for n in ("fs.csv", "passos.csv", "pw.csv", "triangulos.csv", "mecanismo.csv",
                                              "campos.csv"))
if th is None and fs is not None:
    th = float(fs["t_h"][np.argmin(fs["FS_SRM"])])
banda = None
if mec is not None and th is not None:
    m = mec[np.isclose(mec["t_h"], th)]
    if m.size:
        banda = m[m["epsp_raizJ2"] >= 0.1 * m["epsp_raizJ2"].max()] # o mesmo criterio de BANDA_LIMIAR = 0.1
fig, ax = plt.subplots(3, 2, figsize=(14, 12), constrained_layout=True)


def estilo(a, x, y):
    a.set_xlabel(x)
    a.set_ylabel(y)
    a.grid(color="0.9", lw=0.6)
    a.set_axisbelow(True)


a = ax[0, 0]
if fs is not None:
    o = np.argsort(fs["t_h"])
    a.plot(fs["t_h"][o], fs["FS_SRM"][o], "-o", color=TINTA, lw=2, ms=6)
    a.axhline(1, color=CINZA, ls="--", lw=1)
    if th is not None:
        a.axvline(th, color=CINZA, lw=1)
a.set_title("FS (SRM)")
estilo(a, "t (h)", "FS")

if pas is not None:
    t, dt = pas["t_h"], pas["dt_h"]
    a = ax[0, 1]
    # a chuva por cima e tracejada: sem escoamento superficial, coincide com a infiltracao
    for col, rot, cor, ls, z in (("chuva_m3", "chuva", AZUL, "--", 3), (None, "infiltracao (fluxo + face)", LARANJA, "-", 2),
                                 ("esq_m3", "lateral esquerda", AQUA, "-", 2), ("dir_m3", "lateral direita", AMARELO, "-", 2)):
        v = pas["inf_fluxo_m3"] + pas["inf_face_m3"] if col is None else pas[col]
        a.plot(t, v / dt, color=cor, lw=2, ls=ls, zorder=z, label=rot)
    a.set_title("entrada de agua (positiva) por hora")
    a.legend(fontsize=8)
    estilo(a, "t (h)", "m3 por metro de talude, por hora")
    a = ax[1, 0]
    a.plot(t, pas["inf_acum_m3"], color=LARANJA, lw=2, label="infiltracao acumulada")
    a.plot(t, np.cumsum(pas["inf_fluxo_m3"] + pas["inf_face_m3"] + pas["esq_m3"] + pas["dir_m3"]), color=VERDE, lw=2,
           label="entradas acumuladas (superficie + laterais)")
    a.plot(t, np.cumsum(pas["dV_m3"]), color=VIOLETA, lw=2, ls="--", label="variacao do armazenamento acumulada")
    a.set_title("acumulados (entradas = armazenamento: balanco fechado)")
    a.legend(fontsize=8)
    estilo(a, "t (h)", "m3 por metro de talude")

a = ax[1, 1]
if pw is not None and tri is not None and th is not None:
    s = pw[np.isclose(pw["t_h"], th)]
    n = int(max(s["no"].max(), tri["n0"].max(), tri["n1"].max(), tri["n2"].max())) + 1
    X, Y, P = (np.full(n, np.nan) for _ in range(3))
    X[s["no"]], Y[s["no"]], P[s["no"]] = s["x"], s["y"], s["pw_kPa"]
    T = mtri.Triangulation(X, Y, np.column_stack([tri["n0"], tri["n1"], tri["n2"]]))
    lim = max(abs(np.nanmin(P)), abs(np.nanmax(P)))
    cf = a.tricontourf(T, P, levels=30, cmap="RdBu", norm=mcolors.TwoSlopeNorm(0., -lim, lim))
    fig.colorbar(cf, ax=a, label="p_w (kPa)")
    a.tricontour(T, P, levels=[0.], colors=TINTA, linewidths=1.5)
    f = s[s["face"] == 1]
    if f.size:
        a.plot(f["x"], f["y"], "s", color=TINTA, ms=4, label="face de percolacao")
    if banda is not None:
        a.plot(banda["x"], banda["y"], ".", color=TINTA, ms=2, alpha=0.5, label="banda de ruptura")
    if f.size or banda is not None:
        a.legend(fontsize=8, loc="lower left")
    a.set_title(f"p_w em t = {th:g} h (linha: p_w = 0)")
    a.set_aspect("equal")
    estilo(a, "x (m)", "y (m)")

a = ax[2, 0]
if cam is not None:
    T = mtri.Triangulation(cam["x"], cam["y"])
    xt, yt = cam["x"][T.triangles], cam["y"][T.triangles]
    T.set_mask(np.max(np.hypot(xt - np.roll(xt, 1, axis=1), yt - np.roll(yt, 1, axis=1)), axis=1) > 3.0)
    cf = a.tricontourf(T, cam["c"], levels=30, cmap="Greens")
    fig.colorbar(cf, ax=a, label="c' (kPa)")
    if banda is not None:
        a.plot(banda["x"], banda["y"], ".", color=TINTA, ms=2, alpha=0.5, label="banda de ruptura")
        a.legend(fontsize=8, loc="lower left")
    a.set_title("coesao c' (campo da realizacao)")
    a.set_aspect("equal")
    estilo(a, "x (m)", "y (m)")
else:
    a.axis("off")

a = ax[2, 1]
if mec is not None and th is not None:
    m = mec[np.isclose(mec["t_h"], th)]
    o = np.argsort(m["epsp_raizJ2"])
    sc = a.scatter(m["x"][o], m["y"][o], c=m["epsp_raizJ2"][o], s=6, cmap="Purples")
    fig.colorbar(sc, ax=a, label="sqrt(J2) da deformacao plastica")
    a.set_title(f"mecanismo em t = {th:g} h (FS = {m['FS'][0]:g})" if m.size else "mecanismo")
    a.set_aspect("equal")
    estilo(a, "x (m)", "y (m)")
else:
    a.axis("off")
fig.suptitle(os.path.basename(os.path.normpath(R)))
fig.savefig(saida, dpi=110)
print(saida)
