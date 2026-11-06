#!/usr/bin/env python3
# digitaliza_fig15.py <yuan2021.pdf> <pasta de trabalho>: digitaliza a Fig. 15 de Yuan et al. (2021) ("Least Squares
# Circle", pagina 8): calibra os eixos pelas marcas, acha o centro (ponto vermelho) e os nos pretos da banda de
# deformacao plastica e ajusta circulos a eles. Precisa de pdftoppm (poppler), numpy e Pillow. Gerou
# dados/yuan2021_fig15_circulo.csv (digitalizado, aproximado). As constantes de recorte e de cor valem para a
# renderizacao a 400 dpi do PDF do artigo.
import os
import subprocess
import sys

import numpy as np
from PIL import Image


def main():
    pdf, pasta = sys.argv[1], sys.argv[2]
    os.makedirs(pasta, exist_ok=True)
    base = os.path.join(pasta, "p8")
    subprocess.run(["pdftoppm", "-f", "8", "-l", "8", "-r", "400", "-png", pdf, base], check=True)
    im = Image.open(base + "-08.png")
    w, h = im.size
    im = np.asarray(im.crop((int(0.06 * w), int(0.60 * h), int(0.50 * w), int(0.86 * h))).convert("RGB")).astype(int)
    R, G, B = im[..., 0], im[..., 1], im[..., 2]
    escuro = (R < 90) & (G < 90) & (B < 90)
    # marcas dos eixos (abaixo da moldura de baixo e a esquerda da moldura da esquerda); valores das marcas da figura
    linhas = [i for i in range(im.shape[0]) if escuro[i].sum() > 0.5 * im.shape[1]]
    colunas = [j for j in range(im.shape[1]) if escuro[:, j].sum() > 0.4 * im.shape[0]]
    yb, xe = max(linhas), min(colunas)
    sub = escuro[yb + 2:yb + 10, :].sum(axis=0)
    mx = [j for j in range(xe + 5, im.shape[1] - 5) if sub[j] >= 7]
    gx = np.split(np.array(mx), np.where(np.diff(mx) > 3)[0] + 1)
    xt = np.array([g.mean() for g in gx])
    sub = escuro[:, xe - 12:xe - 1].sum(axis=1)
    my = [i for i in range(min(linhas) + 2, yb - 2) if sub[i] >= 9]
    gy = np.split(np.array(my), np.where(np.diff(my) > 3)[0] + 1)
    yt = np.array([g.mean() for g in gy])
    xv = np.array([-75, -50, -25, 0, 25, 50, 75, 100])[:len(xt)]
    yv = np.array([120, 100, 80, 60, 40, 20, 0])[:len(yt)]
    ax, bx = np.polyfit(xv, xt, 1)
    ay, by = np.polyfit(yv, yt, 1)
    tox = lambda px: (px - bx) / ax  # noqa: E731
    toy = lambda py: (py - by) / ay  # noqa: E731
    print(f"escala: {ax:.4f} px/m em x, {-ay:.4f} px/m em y")
    verm = (R > 200) & (G < 80) & (B < 80)
    ys, xs = np.nonzero(verm)
    a, b = tox(xs.mean()), toy(ys.mean())
    H = 50.
    r = np.hypot(a, b)
    c = np.sqrt(r * r - (H - b) ** 2) + a
    print(f"centro ({a:.2f}, {b:.2f}) m, r = {r:.2f} m, saida na crista em x = {c:.2f} m (Eq. 6)")
    # nos pretos: componentes conexas pequenas de pixels muito escuros na regiao da banda
    preto = (R < 60) & (G < 60) & (B < 60)
    x0, x1 = int(bx + 0.3 * ax), int(bx + 61 * ax)
    y0, y1 = int(by + 52 * ay), int(by - 0.5 * ay)
    S = preto[y0:y1, x0:x1]
    lab = np.zeros(S.shape, int)
    pts = []
    for i in range(S.shape[0]):
        for j in range(S.shape[1]):
            if S[i, j] and not lab[i, j]:
                pilha, comp = [(i, j)], []
                lab[i, j] = 1
                while pilha:
                    p, q = pilha.pop()
                    comp.append((p, q))
                    for dp, dq in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                        u, v = p + dp, q + dq
                        if 0 <= u < S.shape[0] and 0 <= v < S.shape[1] and S[u, v] and not lab[u, v]:
                            lab[u, v] = 1
                            pilha.append((u, v))
                comp = np.array(comp)
                if 4 <= len(comp) <= 40 and np.ptp(comp[:, 0]) < 8 and np.ptp(comp[:, 1]) < 8:
                    pts.append((tox(comp[:, 1].mean() + x0), toy(comp[:, 0].mean() + y0)))
    pts = np.array(sorted(pts, key=lambda t: t[1]))
    X, Y = pts[:, 0], pts[:, 1]
    s = np.linalg.lstsq(np.c_[2 * X, 2 * Y], X * X + Y * Y, rcond=None)[0]
    print(f"{len(pts)} nos; circulo pelos nos forcado a passar pelo pe: centro ({s[0]:.2f}, {s[1]:.2f}), r = {np.hypot(*s):.2f}")
    saida = os.path.join(pasta, "yuan2021_fig15_circulo.csv")
    with open(saida, "w") as f:
        f.write("tipo,x_m,y_m\n")
        f.write(f"centro,{a:.2f},{b:.2f}\npe,0.00,0.00\nsaida,{c:.2f},{H:.2f}\n")
        for x, y in pts:
            f.write(f"no,{x:.3f},{y:.3f}\n")
    print("gravado", saida)


if __name__ == "__main__":
    main()
