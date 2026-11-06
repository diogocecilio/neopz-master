#!/usr/bin/env python3
# digitaliza_fig3.py: os marcadores da Fig. 3 de Liu et al. (2024) (quadrados pretos: Jiang et al. 2014; estrelas
# vermelhas: "This study"), a cada 0,5 m de cota. Uso, com os recortes da pagina 5 renderizada a 300 dpi
# (pdftoppm -r 300 -f 5 -l 5; convert -crop 1150x900+150+180 para (a) e 1150x900+1250+180 para (b)):
#   python3 -I digitaliza_fig3.py fig3a.png 221.5 1042.5 64 662.5 10     (pixels de x = -20 e 50 kPa, z = 10 e 0 m)
#   python3 -I digitaliza_fig3.py fig3b.png 104.5 936.5 62.5 669 7.5
# Imprime z, p_w (quadrado) e p_w (estrela); as leituras espurias (legenda, eixos) foram descartadas a mao ao montar
# dados/fig3_secao_*.csv, que usa a estrela ("This study") e, onde ela fica escondida, o quadrado (as duas coincidem).
import sys
import numpy as np
import matplotlib.image as mpimg
arq, x0, x1, y10, y0, zmax = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6])
img = mpimg.imread(arq)[:, :, :3]
r, g, b = img[:, :, 0], img[:, :, 1], img[:, :, 2]
preto = (r < 0.4) & (g < 0.4) & (b < 0.4)
verm = (r > 0.6) & (g < 0.45) & (b < 0.45)
kx = 70. / (x1 - x0)
ky = (y0 - y10) / 10.
zero = x0 + 20. / kx
def px2p(px): return (px - x0) * kx - 20.
for z in np.arange(0., zmax + 1e-9, 0.5):
    py = y0 - z * ky
    lin = slice(int(round(py - 7)), int(round(py + 8)))
    out = []
    for nome, m in (("quad", preto), ("estr", verm)):
        sub = m[lin, :].copy()
        for c in (x0, x1, zero):
            sub[:, int(c) - 3:int(c) + 4] = False
        cols = np.where(sub.sum(axis=0) > 0)[0]
        cols = cols[(cols > x0 + 3) & (cols < x1 - 3)]
        if len(cols) == 0:
            out.append(float('nan')); continue
        # agrupa colunas contiguas e escolhe o grupo com mais pixels
        grupos, atual = [], [cols[0]]
        for c in cols[1:]:
            if c - atual[-1] <= 2: atual.append(c)
            else: grupos.append(atual); atual = [c]
        grupos.append(atual)
        best = max(grupos, key=lambda gg: sub[:, gg].sum())
        w = sub[:, best].sum(axis=0)
        out.append(px2p(np.sum(np.array(best) * w) / w.sum()))
    print(f"{z:4.1f} {out[0]:7.2f} {out[1]:7.2f}")
