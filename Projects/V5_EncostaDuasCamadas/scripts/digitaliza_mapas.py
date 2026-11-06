#!/usr/bin/env python3
# V5_EncostaDuasCamadas: digitaliza os mapas de cores de Abbasov et al. (2024): Fig. 4 (theta), Fig. 5 (tensao de
# succao) e Fig. 6 (LFS), nos instantes 0, 7,5 e 15 h, na zona de investigacao (x_zona = x - 15 m, z_zona = z - 17 m).
# Cada pixel e convertido no valor da cor mais proxima da barra de cores da figura (distancia RGB; pixels brancos,
# pretos ou cinza, como as linhas da interface e do contorno LFS = 1, sao descartados). DIGITALIZADO, APROXIMADO:
# resolucao de ~0,02 m por pixel a 250 dpi; a barra de cores satura nos extremos (theta < 0,22 ou > 0,40; LFS < 0,5 ou
# > 3; tensao de succao > 8 kPa ficam no extremo).
#
# Uso: python3 digitaliza_mapas.py <abaso2024.pdf> [--saida ../dados]
# Gera fig4_theta.csv, fig5_sigs.csv, fig6_lfs.csv e fig9_lfs.csv (Fig. 9, painel COMSOL, 15 h) (grade de 0,25 m:
# t_h, x_zona, z_zona, valor) e os perfis verticais em x_zona = 2, 10, 17 e 20,5 m (passo de 0,05 m) em
# perfis_fig4_5_6.csv (a variavel da Fig. 9 e lfs_fig9).
import argparse
import os
import subprocess
import tempfile
import numpy as np
from PIL import Image

# geometria dos paineis a 250 dpi (pixels da pagina, medida nos centros dos rotulos dos eixos): x0 e o pixel de
# x_zona = 0 e escx a escala (px/m); para cada instante (de cima para baixo), o pixel de z_zona = 0 e a escala em z;
# barra de cores: coluna e faixa de linhas, e dois rotulos (linha, valor)
FIGS = {
    'fig4': dict(pagina=5, x0=359.5, escx=57.69, z0=[(704.5, 57.94), (1281.5, 58.19), (1857.0, 58.19)],
                 barra=(1743, 825, 1340), rot=((860., 0.40), (1333., 0.22)), nome='theta'),
    'fig5': dict(pagina=6, x0=600.0, escx=32.78, z0=[(548.5, 32.90), (939.0, 33.05), (1328.0, 33.0)],
                 barra=(1473, 605, 1055), rot=((622., 8.0), (1042., 0.0)), nome='sigs_kPa'),
    'fig6': dict(pagina=7, x0=421.5, escx=52.18, z0=[(657.0, 52.56), (1180.0, 52.81), (1704.5, 52.63)],
                 barra=(1692, 785, 1222), rot=((793., 3.0), (1215., 0.5)), nome='lfs'),
    # Fig. 9, painel de cima (COMSOL), 15 h, K_s(l) = 0,06 cm/h, nu = 0,35, phi' = 35 graus, c' = 5 kPa (caso fig9)
    'fig9': dict(pagina=11, x0=318.0, escx=63.05, z0=[(742.0, 63.06)], tempos=[15.], barra=(1819, 560, 985),
                 rot=((563., 3.0), (976., 0.5)), nome='lfs_fig9'),
}
TEMPOS = [0., 7.5, 15.]


def superficie(x):  # z_zona da superficie
    return 8. if x <= 4. else (2. if x >= 19. else 8. - 6. * (x - 4.) / 15.)


def barra_de_cores(im, conf):
    xb, ya, yb = conf['barra']
    (y1, v1), (y2, v2) = conf['rot']
    cores, vals = [], []
    for y in range(ya, yb):
        c = im[y, xb - 3:xb + 4].mean(axis=0)
        if c.max() - c.min() < 40:
            continue  # moldura ou fundo: so os pixels coloridos da barra
        cores.append(c)
        vals.append(v1 + (y - y1) * (v2 - v1) / (y2 - y1))
    vals = np.clip(np.array(vals), min(v1, v2), max(v1, v2))
    return np.array(cores), vals


def valor(im, cores, vals, px, py, r=2):
    bloco = im[int(py) - r:int(py) + r + 1, int(px) - r:int(px) + r + 1].reshape(-1, 3)
    out = []
    for c in bloco:
        if c.max() - c.min() < 40 and (c.mean() > 200 or c.mean() < 120):
            continue  # branco, preto ou cinza: fundo, linhas, texto
        d = ((cores - c) ** 2).sum(axis=1)
        k = int(np.argmin(d))
        if d[k] < 60. ** 2:
            out.append(vals[k])
    return float(np.median(out)) if out else np.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('pdf')
    ap.add_argument('--saida', default=os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'dados'))
    args = ap.parse_args()
    tmp = tempfile.mkdtemp()
    perfis = []
    for fig, conf in FIGS.items():
        pref = os.path.join(tmp, fig)
        subprocess.run(['pdftoppm', '-r', '250', '-f', str(conf['pagina']), '-l', str(conf['pagina']), '-png',
                        args.pdf, pref], check=True)
        arq = [f for f in os.listdir(tmp) if f.startswith(fig) and f.endswith('.png')][0]
        im = np.asarray(Image.open(os.path.join(tmp, arq)).convert('RGB')).astype(float)
        cores, vals = barra_de_cores(im, conf)
        linhas = []
        for t, (y0, escz) in zip(conf.get('tempos', TEMPOS), conf['z0']):
            for x in np.arange(0.125, 22., 0.25):
                for z in np.arange(0.125, 8., 0.25):
                    if z > superficie(x) - 0.08:
                        continue
                    v = valor(im, cores, vals, conf['x0'] + x * conf['escx'], y0 - z * escz)
                    if np.isfinite(v):
                        linhas.append((t, x, z, v))
            for x in (2., 10., 17., 20.5):
                zs = superficie(x)
                for d in np.arange(0.03, zs, 0.05):
                    v = valor(im, cores, vals, conf['x0'] + x * conf['escx'], y0 - (zs - d) * escz, r=1)
                    perfis.append((fig, conf['nome'], t, x, zs - d, d, v))
        with open(os.path.join(args.saida, f"{fig}_{conf['nome'].split('_')[0]}.csv"), 'w') as f:
            f.write(f"# Abbasov et al. (2024), {fig.replace('fig', 'Fig. ')}: {conf['nome']} na zona de investigacao "
                    "(x_zona = x - 15 m, z_zona = z - 17 m); digitalizado, aproximado (scripts/digitaliza_mapas.py)\n")
            f.write(f"t_h,x_zona,z_zona,{conf['nome']}\n")
            for l in linhas:
                f.write(f"{l[0]},{l[1]:.3f},{l[2]:.3f},{l[3]:.4f}\n")
    with open(os.path.join(args.saida, 'perfis_fig4_5_6.csv'), 'w') as f:
        f.write("# Abbasov et al. (2024), Figs. 4 a 6: perfis verticais nos mapas de cores; digitalizado, aproximado "
                "(scripts/digitaliza_mapas.py); nan = linha preta ou fora da barra de cores\n")
        f.write("figura,variavel,t_h,x_zona,z_zona,prof,valor\n")
        for p in perfis:
            f.write(f"{p[0]},{p[1]},{p[2]},{p[3]},{p[4]:.3f},{p[5]:.3f},{p[6]:.4f}\n")


if __name__ == '__main__':
    main()
