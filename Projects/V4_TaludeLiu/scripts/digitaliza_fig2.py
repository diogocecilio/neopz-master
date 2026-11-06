#!/usr/bin/env python3
# digitaliza_fig2.py: a superficie critica (Bishop) da Fig. 2 de Liu et al. (2024) e a area da massa que desliza
# entre ela e o terreno. Uso: python3 -I digitaliza_fig2.py <fig2.png> (recorte da pagina 4 renderizada a 300 dpi com
# pdftoppm -r 300; os pixels dos eixos abaixo sao os desse recorte). Digitalizado, aproximado.
import sys
import numpy as np
import matplotlib.image as mpimg
img = mpimg.imread(sys.argv[1])[:, :, :3]
preto = (img[:, :, 0] < 0.3) & (img[:, :, 1] < 0.3) & (img[:, :, 2] < 0.3)
# eixos (pixels do recorte): x = 0 e 30 m, z = 0 e 10 m
cols = preto.sum(axis=0); rows = preto.sum(axis=1)
X0, X30 = 136., 1035.
Z0, Z10 = 417., 118.
def m2px(x, z): return X0 + x * (X30 - X0) / 30., Z0 - z * (Z0 - Z10) / 10.
def px2m(px, py): return (px - X0) * 30. / (X30 - X0), (Z0 - py) * 10. / (Z0 - Z10)
def terreno(x): return 10. if x <= 10 else (5. if x >= 20 else 10 - 0.5 * (x - 10))
pts = []
zant = 10.
for px in range(int(m2px(6.0, 0)[0]), int(m2px(21.6, 0)[0])):
    x = px2m(px, 0)[0]
    pyt = m2px(x, terreno(x))[1]
    ys = np.where(preto[:, px])[0]
    ys = ys[(ys > pyt + 7) & (ys < Z0 - 3)]  # abaixo da linha do terreno
    if len(ys) == 0: continue
    # grupos de pixels contiguos; a superficie critica e o grupo mais proximo do ponto anterior (a seta do rotulo
    # "Critical slip surface" cruza a regiao e e descartada por continuidade)
    grupos, g = [], [ys[0]]
    for y in ys[1:]:
        if y - g[-1] <= 2: g.append(y)
        else: grupos.append(g); g = [y]
    grupos.append(g)
    zs = [px2m(px, np.mean(gg))[1] for gg in grupos]
    z = min(zs, key=lambda v: abs(v - zant))
    if abs(z - zant) > 0.6 and pts: continue
    zant = z
    pts.append((x, z))
pts = np.array(pts)
# corta onde a linha encontra o terreno (entrada e saida)
dz = np.array([terreno(x) - z for x, z in pts])
ok = dz > 0.05
pts = pts[ok]
xe, xs = pts[0, 0], pts[-1, 0]
i = np.argmin(pts[:, 1])
area = np.trapz([terreno(x) - z for x, z in pts], pts[:, 0]) if hasattr(np, 'trapz') else np.trapezoid([terreno(x) - z for x, z in pts], pts[:, 0])
print(f"entrada x = {xe:.2f} m, saida x = {xs:.2f} m, ponto mais baixo ({pts[i,0]:.2f}, {pts[i,1]:.2f}) m, D_f = {10-pts[i,1]:.2f} m")
print(f"area entre o terreno e a superficie critica: {area:.1f} m2")
for x in np.arange(7, 21.1, 1.):
    j = np.argmin(abs(pts[:, 0] - x)); print(f"  x = {pts[j,0]:5.2f}  z = {pts[j,1]:5.2f}")
if len(sys.argv) > 2:  # grava a superficie a cada 0,5 m
    with open(sys.argv[2], 'w') as f:
        f.write("# Liu et al. (2024), Fig. 2: superficie critica (Bishop, SLOPE/W) do caso deterministico. Digitalizado, "
                "aproximado (scripts/digitaliza_fig2.py, pagina 4 a 300 dpi).\n")
        f.write(f"# entrada x = {xe:.2f} m, saida x = {xs:.2f} m (o fim da linha encontra o terreno em x ~ 21 m), ponto mais "
                f"baixo ({pts[i,0]:.2f}, {pts[i,1]:.2f}) m, D_f = {10-pts[i,1]:.2f} m, area da massa {area:.1f} m2\n")
        f.write("x_m,z_m\n")
        for x in np.arange(7.0, 20.51, 0.5):
            j = np.argmin(abs(pts[:, 0] - x)); f.write(f"{pts[j,0]:.2f},{pts[j,1]:.2f}\n")
