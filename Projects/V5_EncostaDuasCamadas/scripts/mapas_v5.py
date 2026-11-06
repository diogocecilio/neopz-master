#!/usr/bin/env python3
# V5_EncostaDuasCamadas: mapas de theta, tensao de succao S_e |p_w| e LFS na zona de investigacao (x_zona = x - 15 m,
# z_zona = z - 17 m) de uma execucao (nos.csv, elementos.csv, campos_t<h>.csv), lado a lado com os mapas digitalizados
# de Abbasov et al. (2024), Figs. 4 a 6 (dados/fig4_theta.csv, fig5_sigs.csv, fig6_lfs.csv), nas mesmas escalas de cor.
# As figuras do artigo sao do caso tabela1 (Tabela 1); a comparacao de outro caso e so qualitativa.
#
# Uso: python3 mapas_v5.py <pasta da execucao> [--saida mapas.png] [--tempos 0,7.5,15] [--sem-artigo]
import argparse
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

DADOS = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'dados')
X0, Z0 = 15., 17.
VARS = [('theta', 'theta [-]', 0.22, 0.40, 'jet_r', 'fig4_theta.csv'),
        ('sigs', 'tensao de succao S_e|p_w| [kPa]', 0., 8., 'jet', 'fig5_sigs.csv'),
        ('lfs', 'LFS [-]', 0.5, 3., 'jet_r', 'fig6_lfs.csv')]


def le(arq):
    # os CSV de dados/ tem uma linha de comentario antes do cabecalho
    with open(arq) as f:
        pula = 1 if f.readline().startswith('#') else 0
    return np.genfromtxt(arq, delimiter=',', names=True, skip_header=pula)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('pasta')
    ap.add_argument('--saida', default='')
    ap.add_argument('--tempos', default='0,7.5,15')
    ap.add_argument('--sem-artigo', action='store_true')
    args = ap.parse_args()
    tempos = [float(t) for t in args.tempos.split(',')]
    nos = le(os.path.join(args.pasta, 'nos.csv'))
    els = le(os.path.join(args.pasta, 'elementos.csv'))
    zona = {int(e['el']): e for e in els if e['zona'] == 1}
    ncol = 1 if args.sem_artigo else 2
    fig, ax = plt.subplots(len(tempos), 3 * ncol, figsize=(7.2 * ncol * 3 / 2, 2.3 * len(tempos)), squeeze=False)
    for i, t in enumerate(tempos):
        arq = os.path.join(args.pasta, f'campos_t{t:g}.csv')
        if not os.path.exists(arq):
            continue
        c = le(arq)
        sel = [r for r in c if int(r['el']) in zona]
        # triangulacao com os vertices repetidos por elemento (theta e succao com a camada de cada elemento)
        xs, zs, tri = [], [], []
        val = {'theta': [], 'sigs': []}
        lfs = []
        for r in sel:
            e = zona[int(r['el'])]
            k = len(xs)
            for j in range(3):
                n = int(e[f'n{j}'])
                xs.append(nos['x'][n] - X0)
                zs.append(nos['z'][n] - Z0)
                val['theta'].append(r[f'theta{j}'])
                val['sigs'].append(r[f'sigs{j}'])
            tri.append([k, k + 1, k + 2])
            lfs.append(r['lfs'])
        T = mtri.Triangulation(np.array(xs), np.array(zs), np.array(tri))
        for j, (v, rot, vmin, vmax, cm, fdig) in enumerate(VARS):
            a = ax[i, ncol * j]
            if v == 'lfs':
                h = a.tripcolor(T, facecolors=np.array(lfs), cmap=cm, vmin=vmin, vmax=vmax)
                a.tricontour(T, np.repeat(np.array(lfs), 3), levels=[1.], colors='k', linewidths=0.6)
            else:
                h = a.tripcolor(T, np.array(val[v]), shading='gouraud', cmap=cm, vmin=vmin, vmax=vmax)
            a.set_title(f'nosso: {v}, t = {t:g} h', fontsize=8)
            plt.colorbar(h, ax=a, fraction=0.025)
            if not args.sem_artigo:
                b = ax[i, ncol * j + 1]
                d = le(os.path.join(DADOS, fdig))
                nome = d.dtype.names[3]
                s = d[d['t_h'] == t]
                h2 = b.scatter(s['x_zona'], s['z_zona'], c=s[nome], s=5, marker='s', cmap=cm, vmin=vmin, vmax=vmax)
                b.set_title(f'artigo ({fdig.split("_")[0]}, digitalizado): t = {t:g} h', fontsize=8)
                plt.colorbar(h2, ax=b, fraction=0.025)
    for a in ax.flat:
        a.set_aspect('equal')
        a.set_xlim(0, 22)
        a.set_ylim(0, 8.2)
        a.plot([0, 4, 19, 22], [6, 6, 0, 0], 'k-', lw=0.5)
        a.tick_params(labelsize=6)
    fig.tight_layout()
    fig.savefig(args.saida or os.path.join(args.pasta, 'mapas.png'), dpi=120)


if __name__ == '__main__':
    main()
