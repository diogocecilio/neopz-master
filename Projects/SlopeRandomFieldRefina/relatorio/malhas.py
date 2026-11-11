#!/usr/bin/env python3
"""malhas.py: malhas gmsh v4.1 (triangulos lineares) do talude de Vargas Ceron et al. (2025), Secao 6.

Funcoes reutilizaveis:
  ler_msh(caminho)                -> (xy[N,2], tri[M,3]) com indices a partir de 0
  estatisticas(xy, tri)           -> dict (nos, triangulos, area, arestas, h equivalente)
  pontos_banda_geo(caminho_geo)   -> pontos Point(...) que o refina_banda.py acrescentou ao .geo
  banda(mecanismo_csv, fs_csv)    -> (x, y) dos pontos com sqrt(J2) >= limiar * maximo no instante do fator
  desenhar(ax, xy, tri, ...)      -> malha (arestas finas) + contorno + banda + zoom
  gerar_msh(geo, msh)             -> malha 2D com o modulo python do gmsh (1 thread)
  fmt(v, nd)                      -> numero com virgula decimal

CLI:
  python3 malhas.py estat malha.msh [outra.msh ...]
  python3 malhas.py gerar malha.geo saida.msh
  python3 malhas.py figura malha.msh -o saida.png [--mecanismo mecanismo.csv --fs fs.csv]
                    [--pontos-geo malha.geo] [--zoom x0 x1 y0 y1] [--titulo T] [--dpi 200]
"""
import argparse
import csv
import re
import sys

import numpy as np

# geometria do talude (Secao 6): pe em (15, 5), crista em (10, 10)
GEOM = np.array([(0, 0), (25, 0), (25, 5), (15, 5), (10, 10), (0, 10)], float)
PE = (15.0, 5.0)
CRISTA = (10.0, 10.0)
ZOOM_TALUDE = (4.0, 21.0, 2.0, 10.5)
CORES = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#6250d6', '#e34948']
COR_ARESTA = '#3c3c3c'
COR_BANDA = '#e34948'    # banda da rodada mostrada no painel
COR_BANDA1 = '#2a78d6'   # banda da etapa 1 (pontos que guiaram o refino)


def fmt(v, nd=3, sinal=False):
    """Numero com virgula decimal (pt-BR)."""
    s = f'{v:+.{nd}f}' if sinal else f'{v:.{nd}f}'
    return s.replace('.', ',')


def ler_msh(caminho):
    """Le um .msh ASCII v4.x: devolve as coordenadas (x, y) dos nos e os triangulos lineares (tipo 2)."""
    with open(caminho) as f:
        linhas = f.read().split('\n')
    i, n = 0, len(linhas)
    tags, xy, tri = [], [], []
    while i < n:
        lin = linhas[i].strip()
        if lin == '$MeshFormat':
            ver = linhas[i + 1].split()
            if not ver[0].startswith('4') or ver[1] != '0':
                raise ValueError(f'{caminho}: formato {ver} nao suportado (use 4.1 ASCII)')
            i += 2
        elif lin == '$Nodes':
            nblocos = int(linhas[i + 1].split()[0])
            i += 2
            for _ in range(nblocos):
                m = int(linhas[i].split()[3])
                i += 1
                tags.extend(int(linhas[i + k]) for k in range(m))
                i += m
                for k in range(m):
                    c = linhas[i + k].split()
                    xy.append((float(c[0]), float(c[1])))
                i += m
        elif lin == '$Elements':
            nblocos = int(linhas[i + 1].split()[0])
            i += 2
            for _ in range(nblocos):
                _, _, tipo, m = map(int, linhas[i].split())
                i += 1
                if tipo == 2:
                    for k in range(m):
                        c = linhas[i + k].split()
                        tri.append((int(c[1]), int(c[2]), int(c[3])))
                i += m
        else:
            i += 1
    tags = np.array(tags, int)
    xy = np.array(xy, float)
    idx = np.full(tags.max() + 1, -1, int)
    idx[tags] = np.arange(len(tags))
    t = idx[np.array(tri, int)]
    if (t < 0).any():
        raise ValueError(f'{caminho}: triangulo com no inexistente')
    return xy, t


def arestas(tri):
    e = np.sort(np.concatenate([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]]), axis=1)
    return np.unique(e, axis=0)


def estatisticas(xy, tri):
    p = xy[tri]
    area = 0.5 * np.abs((p[:, 1, 0] - p[:, 0, 0]) * (p[:, 2, 1] - p[:, 0, 1])
                        - (p[:, 2, 0] - p[:, 0, 0]) * (p[:, 1, 1] - p[:, 0, 1]))
    e = arestas(tri)
    comp = np.linalg.norm(xy[e[:, 1]] - xy[e[:, 0]], axis=1)
    return dict(nos=len(xy), triangulos=len(tri), arestas=len(e), area_m2=float(area.sum()),
                aresta_min_m=float(comp.min()), aresta_media_m=float(comp.mean()),
                aresta_max_m=float(comp.max()), area_tri_min_m2=float(area.min()),
                area_tri_max_m2=float(area.max()))


def pontos_banda_geo(caminho):
    """Pontos Point(...) acrescentados pelo refina_banda.py (depois do comentario '// refina_banda.py')."""
    pts, dentro = [], False
    rx = re.compile(r'\s*Point\(\s*\d+\s*\)\s*=\s*\{\s*([-+0-9.eE]+)\s*,\s*([-+0-9.eE]+)')
    with open(caminho) as f:
        for lin in f:
            if lin.startswith('// refina_banda.py'):
                dentro = True
                continue
            if dentro:
                m = rx.match(lin)
                if m:
                    pts.append((float(m.group(1)), float(m.group(2))))
    return np.array(pts, float).reshape(-1, 2)


def ler_csv(caminho):
    with open(caminho) as f:
        return list(csv.DictReader(f))


def instante_fator(fs_csv):
    """(t_h, FS) da linha de menor FS_SRM do fs.csv."""
    lin = min(ler_csv(fs_csv), key=lambda r: float(r['FS_SRM']))
    return float(lin['t_h']), float(lin['FS_SRM'])


def banda(mecanismo_csv, fs_csv=None, limiar=0.1, tol=1e-6, t_h=None):
    """Pontos (x, y) com sqrt(J2) >= limiar * maximo no instante do fator (t_h do menor FS_SRM do fs.csv,
    comparado com tolerancia relativa tol). Sem fs.csv nem t_h, usa o ultimo instante do mecanismo.csv."""
    d = np.genfromtxt(mecanismo_csv, delimiter=',', names=True)
    d = np.atleast_1d(d)
    t = d['t_h']
    if t_h is None:
        t_h = instante_fator(fs_csv)[0] if fs_csv else float(t.max())
    sel = np.abs(t - t_h) <= tol * max(abs(t_h), 1e-300)
    if not sel.any():
        raise ValueError(f'{mecanismo_csv}: nenhum ponto no instante t_h = {t_h}')
    s, x, y = d['epsp_raizJ2'][sel], d['x'][sel], d['y'][sel]
    m = s >= limiar * s.max()
    return x[m], y[m]


def gerar_msh(geo, msh):
    """Gera a malha 2D de um .geo com o modulo python do gmsh, em 1 thread; devolve o caminho do .msh."""
    import gmsh
    gmsh.initialize()
    try:
        gmsh.option.setNumber('General.Terminal', 0)
        gmsh.option.setNumber('General.NumThreads', 1)
        gmsh.open(geo)
        gmsh.option.setNumber('General.NumThreads', 1)
        gmsh.option.setNumber('Mesh.MaxNumThreads2D', 1)
        gmsh.option.setNumber('Mesh.MshFileVersion', 4.1)
        gmsh.option.setNumber('Mesh.Binary', 0)
        gmsh.model.mesh.generate(2)
        gmsh.write(msh)
    finally:
        gmsh.finalize()
    return msh


def _fmt_eixo(v, _pos=None):
    s = f'{v:.10g}'
    return s.replace('.', ',')


def desenhar(ax, xy, tri, banda_xy=None, zoom=None, cor_banda=COR_BANDA, s_banda=3.0, lw=0.22,
             titulo=None, contorno=True, rotulo_banda=None, eixos=True, fs_titulo=8.5, marcador='o'):
    """Desenha a malha (arestas finas cinza-escuro), o contorno do talude e, opcionalmente, a banda."""
    from matplotlib.collections import LineCollection
    from matplotlib.ticker import FuncFormatter
    e = arestas(tri)
    ax.add_collection(LineCollection(xy[e], colors=COR_ARESTA, linewidths=lw, zorder=1,
                                     rasterized=False))
    if contorno:
        g = np.vstack([GEOM, GEOM[:1]])
        ax.plot(g[:, 0], g[:, 1], color='k', lw=0.8, zorder=2)
    if banda_xy is not None and len(banda_xy[0]):
        ax.scatter(banda_xy[0], banda_xy[1], s=s_banda, c=cor_banda, linewidths=0, zorder=3,
                   label=rotulo_banda, marker=marcador)
    ax.set_aspect('equal')
    if zoom:
        ax.set_xlim(zoom[0], zoom[1])
        ax.set_ylim(zoom[2], zoom[3])
    else:
        ax.set_xlim(-0.3, 25.3)
        ax.set_ylim(-0.3, 10.3)
    ax.xaxis.set_major_formatter(FuncFormatter(_fmt_eixo))
    ax.yaxis.set_major_formatter(FuncFormatter(_fmt_eixo))
    ax.tick_params(labelsize=7, length=2.5, width=0.6, pad=1.5)
    for s in ax.spines.values():
        s.set_linewidth(0.6)
    if not eixos:
        ax.set_xticks([])
        ax.set_yticks([])
    if titulo:
        ax.set_title(titulo, fontsize=fs_titulo, pad=3)
    return ax


def _main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    a = sub.add_parser('estat', help='estatisticas de malhas .msh')
    a.add_argument('msh', nargs='+')
    g = sub.add_parser('gerar', help='gera .msh de um .geo (gmsh python, 1 thread)')
    g.add_argument('geo')
    g.add_argument('msh')
    f = sub.add_parser('figura', help='figura de uma malha, com banda e zoom opcionais')
    f.add_argument('msh')
    f.add_argument('-o', '--saida', required=True)
    f.add_argument('--mecanismo')
    f.add_argument('--fs')
    f.add_argument('--limiar', type=float, default=0.1)
    f.add_argument('--pontos-geo')
    f.add_argument('--zoom', type=float, nargs=4, metavar=('X0', 'X1', 'Y0', 'Y1'))
    f.add_argument('--titulo')
    f.add_argument('--dpi', type=int, default=200)
    f.add_argument('--largura', type=float, default=7.0, help='largura da figura (pol.)')
    args = ap.parse_args(argv)
    if args.cmd == 'estat':
        for c in args.msh:
            est = estatisticas(*ler_msh(c))
            print(c + ': ' + ', '.join(f'{k}={v:.6g}' if isinstance(v, float) else f'{k}={v}'
                                       for k, v in est.items()))
    elif args.cmd == 'gerar':
        gerar_msh(args.geo, args.msh)
        est = estatisticas(*ler_msh(args.msh))
        print(f"{args.msh}: {est['triangulos']} triangulos, {est['nos']} nos")
    else:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        xy, tri = ler_msh(args.msh)
        z = args.zoom
        lx, ly = (z[1] - z[0], z[3] - z[2]) if z else (25.6, 10.6)
        fig, ax = plt.subplots(figsize=(args.largura, args.largura * ly / lx + 0.5))
        b = banda(args.mecanismo, args.fs, args.limiar) if args.mecanismo else None
        desenhar(ax, xy, tri, b, zoom=z, titulo=args.titulo or f'{len(tri)} triângulos')
        if args.pontos_geo:  # pontos da banda da etapa 1 gravados no .geo, por baixo da banda da rodada
            p = pontos_banda_geo(args.pontos_geo)
            ax.scatter(p[:, 0], p[:, 1], s=4, c=COR_BANDA1, marker='s', linewidths=0, zorder=2.5)
        ax.set_xlabel('x (m)', fontsize=8)
        ax.set_ylabel('y (m)', fontsize=8)
        fig.tight_layout()
        fig.savefig(args.saida, dpi=args.dpi)
        print(args.saida)


if __name__ == '__main__':
    sys.exit(_main())
