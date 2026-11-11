#!/usr/bin/env python3
"""figuras_malhas.py: figuras de malhas da parte de refinamento do relatorio (usa malhas.py).

  python3 figuras_malhas.py comparacao   # ref_malhas_fixas, ref_metodo_det, ref_refinadas
  python3 figuras_malhas.py niveis       # niv_malhas_iniciais, niv_niveis_<caso>
  python3 figuras_malhas.py todas

Dados (somente leitura): COMP = comparacao anterior (<metodo>/<caso>/ e tabela.csv), NIV = campanha de
niveis (<inicial>/<config>/<caso>/), MALHAS = malhas gmsh fixas. As malhas refinadas da comparacao
(cujo .msh foi apagado) sao regeneradas a partir de malha_refinada.geo em GERADAS. Os caminhos podem ser
trocados pelas variaveis de ambiente REL_COMP, REL_NIV, REL_MALHAS e REL_GERADAS.
"""
import csv
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
import malhas as mm  # noqa: E402

SCR = '/tmp/claude-0/-home-user-neopz-master/674eb771-9e25-5098-bcbd-eda87ce4d4e9/scratchpad'
COMP = os.environ.get('REL_COMP', f'{SCR}/refina/comparacao')
NIV = os.environ.get('REL_NIV', f'{SCR}/niveis')
MALHAS = os.environ.get('REL_MALHAS', os.path.normpath(
    os.path.join(AQUI, '..', '..', 'LearningMohrCoulomb2', 'cluster', 'malhas')))
GERADAS = os.environ.get('REL_GERADAS', f'{SCR}/rel_malhas')
FIG = os.path.join(AQUI, 'figuras')
DADOS = os.path.join(AQUI, 'dados')
DPI = 200
CASOS = ['det'] + [f'r{i:05d}' for i in range(8)]
H_FIXAS = {'foco375': 0.375, 'foco25': 0.25, 'foco1875': 0.1875, 'foco125': 0.125}

plt.rcParams.update({'font.size': 8, 'axes.titlesize': 8.5, 'font.family': 'DejaVu Sans'})


def nome_caso(c):
    return 'det' if c == 'det' else f'r{int(c[1:])}'


def kv(caminho):
    """Le 'chave=valor chave=valor ...' de um .txt (resultado.txt, tempo.txt, malha.txt)."""
    d = {}
    with open(caminho) as f:
        for tok in f.read().split():
            if '=' in tok:
                k, v = tok.split('=', 1)
                d[k] = v
    return d


def tabela_comp():
    t = {}
    with open(os.path.join(COMP, 'tabela.csv')) as f:
        for r in csv.DictReader(f):
            c = 'det' if r['caso'] == 'det' else f"r{int(r['caso']):05d}"
            t[(r['metodo'], c)] = r
    return t


def gamma_foco125(caso):
    return float(kv(os.path.join(COMP, 'foco125', caso, 'resultado.txt'))['FS_min'])


def malha_refinada_comp(caso):
    """Malha refinada da comparacao (h = 0,125 m), regenerada do .geo se ainda nao existir."""
    msh = os.path.join(GERADAS, f'refinada_{caso}.msh')
    if not os.path.exists(msh):
        os.makedirs(GERADAS, exist_ok=True)
        mm.gerar_msh(os.path.join(COMP, 'refinada', caso, 'malha_refinada.geo'), msh)
    return mm.ler_msh(msh)


def banda_dir(d):
    return mm.banda(os.path.join(d, 'mecanismo.csv'), os.path.join(d, 'fs.csv'))


def salvar(fig, nome):
    caminho = os.path.join(FIG, nome)
    fig.savefig(caminho, dpi=DPI)
    plt.close(fig)
    try:  # paleta indexada (256 cores, sem pontilhado): arquivos bem menores, sem perda visivel
        from PIL import Image
        im = Image.open(caminho).convert('RGB')
        im.quantize(colors=256, method=Image.Quantize.MEDIANCUT, dither=Image.Dither.NONE).save(
            caminho, optimize=True, dpi=(DPI, DPI))
    except ImportError:
        pass
    print(f'{caminho}: {os.path.getsize(caminho) / 1e6:.2f} MB')


def legenda_banda(ax, itens, loc='lower left', fs=7):
    h = [Line2D([], [], ls='', marker=m, ms=ms, color=c, label=r) for (r, c, m, ms) in itens]
    ax.legend(handles=h, loc=loc, fontsize=fs, frameon=True, framealpha=0.9, borderpad=0.4,
              handletextpad=0.3, borderaxespad=0.4)


# ---------------------------------------------------------------- comparacao
def fig_malhas_fixas():
    fig, axs = plt.subplots(2, 2, figsize=(7.2, 3.55), sharex=True, sharey=True)
    for ax, nome in zip(axs.flat, H_FIXAS):
        xy, tri = mm.ler_msh(os.path.join(MALHAS, f'ceron_h5_a45_{nome}.msh'))
        mm.desenhar(ax, xy, tri, lw=0.18,
                    titulo=f'{nome} (h = {mm.fmt(H_FIXAS[nome], 4).rstrip("0")} m): {len(tri)} triângulos')
    for ax in axs[1]:
        ax.set_xlabel('x (m)')
    for ax in axs[:, 0]:
        ax.set_ylabel('y (m)')
    fig.tight_layout(pad=0.4, h_pad=0.6, w_pad=0.8)
    salvar(fig, 'ref_malhas_fixas.png')


def fig_metodo_det():
    z = mm.ZOOM_TALUDE
    fig, axs = plt.subplots(2, 1, figsize=(7.0, 7.3), sharex=True)
    xy, tri = mm.ler_msh(os.path.join(MALHAS, 'ceron_h5_a45_foco375.msh'))
    p = mm.pontos_banda_geo(os.path.join(COMP, 'refinada', 'det', 'malha_refinada.geo'))
    g1 = float(kv(os.path.join(COMP, 'refinada', 'det', 'grossa', 'resultado.txt'))['FS_min']) \
        if os.path.exists(os.path.join(COMP, 'refinada', 'det', 'grossa', 'resultado.txt')) else None
    tit = f'(a) etapa 1: malha foco375, {len(tri)} triângulos'
    if g1:
        tit += f', Γ = {mm.fmt(g1)}'
    mm.desenhar(axs[0], xy, tri, (p[:, 0], p[:, 1]), zoom=z, cor_banda=mm.COR_BANDA1, s_banda=9,
                marcador='s', lw=0.35, titulo=tit + f'; {len(p)} pontos da banda', fs_titulo=9)
    legenda_banda(axs[0], [('pontos da banda da etapa 1 (grade de 0,125 m)', mm.COR_BANDA1, 's', 3.5)])
    xy2, tri2 = malha_refinada_comp('det')
    b = banda_dir(os.path.join(COMP, 'refinada', 'det'))
    g2 = float(kv(os.path.join(COMP, 'refinada', 'det', 'resultado.txt'))['FS_min'])
    mm.desenhar(axs[1], xy2, tri2, b, zoom=z, s_banda=4, lw=0.3, fs_titulo=9,
                titulo=f'(b) etapa 2: malha refinada (h = 0,125 m, folga 0,375 m), {len(tri2)} triângulos, '
                       f'Γ = {mm.fmt(g2)}')
    legenda_banda(axs[1], [('banda da rodada refinada', mm.COR_BANDA, 'o', 3.0)])
    for ax in axs:
        ax.set_ylabel('y (m)')
    axs[1].set_xlabel('x (m)')
    fig.tight_layout(pad=0.4, h_pad=0.8)
    salvar(fig, 'ref_metodo_det.png')
    print(f'  det: etapa 1 {len(tri)} tri, {len(p)} pontos, Gamma_grossa={g1}; refinada {len(tri2)} tri, '
          f'banda {len(b[0])} pontos, Gamma={g2}')


def fig_refinadas():
    z = (5.0, 19.5, 3.0, 10.5)
    fig, axs = plt.subplots(4, 2, figsize=(7.2, 8.3), sharex=True, sharey=True)
    for ax, caso in zip(axs.flat, CASOS[1:]):
        xy, tri = malha_refinada_comp(caso)
        b = banda_dir(os.path.join(COMP, 'refinada', caso))
        g = float(kv(os.path.join(COMP, 'refinada', caso, 'resultado.txt'))['FS_min'])
        mm.desenhar(ax, xy, tri, b, zoom=z, s_banda=2.2, lw=0.2,
                    titulo=f'{nome_caso(caso)}: {len(tri)} triângulos, Γ = {mm.fmt(g)}')
        print(f'  {caso}: {len(tri)} tri, banda {len(b[0])} pontos, Gamma={g}')
    for ax in axs[-1]:
        ax.set_xlabel('x (m)')
    for ax in axs[:, 0]:
        ax.set_ylabel('y (m)')
    fig.tight_layout(pad=0.4, h_pad=0.5, w_pad=0.6)
    salvar(fig, 'ref_refinadas.png')


# ---------------------------------------------------------------- niveis
NIVEIS = [0.25, 0.1875, 0.125]
INICIAIS = {'foco375': 'ceron_h5_a45_foco375.msh', 'g15': 'ceron_h5_a45_g15.msh'}


def tag_h(h):
    return f'{h:g}'


def painel_nivel(inicial, h, caso):
    """Dados de um painel: malha, banda, Gamma, erro relativo a foco125 e CPU da cadeia (s)."""
    ini = os.path.join(NIV, inicial, 'inicial', caso)
    if inicial == 'foco375' and h == 0.125:
        d = os.path.join(COMP, 'refinada', caso)
        xy, tri = malha_refinada_comp(caso)
        cpu = float(kv(os.path.join(d, 'tempo.txt'))['cpu_s'])
        cpu_partes = (None, None, cpu)
        fonte = f'comparacao/refinada/{caso}'
    else:
        d = os.path.join(NIV, inicial, f'um_passo_h{tag_h(h)}', caso)
        if not os.path.exists(os.path.join(d, 'tempo.txt')) or not os.path.exists(os.path.join(ini, 'tempo.txt')):
            return None
        xy, tri = mm.ler_msh(os.path.join(d, 'malha.msh'))
        t1 = float(kv(os.path.join(ini, 'tempo.txt'))['cpu_s'])
        tk = kv(os.path.join(d, 'tempo.txt'))
        cpu_partes = (t1, float(tk['cpu_malha_s']), float(tk['cpu_s']))
        cpu = sum(cpu_partes)
        fonte = f'niveis/{inicial}/um_passo_h{tag_h(h)}/{caso}'
        mt = kv(os.path.join(d, 'malha.txt'))
        if int(mt['triangulos']) != len(tri):
            print(f'  AVISO {fonte}: malha.txt diz {mt["triangulos"]} triangulos, .msh tem {len(tri)}')
    if not os.path.exists(os.path.join(d, 'mecanismo.csv')):
        return None
    g = float(kv(os.path.join(d, 'resultado.txt'))['FS_min'])
    g125 = gamma_foco125(caso)
    b = banda_dir(d)
    return dict(xy=xy, tri=tri, banda=b, gamma=g, g125=g125, e=100 * (g / g125 - 1), cpu=cpu,
                cpu_partes=cpu_partes, fonte=fonte)


def fig_malhas_iniciais():
    z = (0.0, 25.0, 0.0, 10.5)
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 1.95), sharey=True)
    for ax, (ini, arq) in zip(axs, INICIAIS.items()):
        xy, tri = mm.ler_msh(os.path.join(MALHAS, arq))
        d = os.path.join(NIV, ini, 'inicial', 'det')
        b = banda_dir(d)
        g = float(kv(os.path.join(d, 'resultado.txt'))['FS_min'])
        mm.desenhar(ax, xy, tri, b, zoom=z, cor_banda=mm.COR_BANDA1, s_banda=2.5, lw=0.3,
                    titulo=f'{ini}: {len(tri)} triângulos, Γ(det) = {mm.fmt(g)}')
        ax.set_xlabel('x (m)')
        print(f'  inicial {ini}: {len(tri)} tri, banda det {len(b[0])} pontos, Gamma={g}')
    axs[0].set_ylabel('y (m)')
    legenda_banda(axs[1], [('banda do caso det', mm.COR_BANDA1, 'o', 3.0)], loc='upper right')
    fig.tight_layout(pad=0.4, w_pad=0.8)
    salvar(fig, 'niv_malhas_iniciais.png')


def fig_niveis(caso, linhas_csv):
    z = (5.0, 19.0, 3.5, 10.5)
    fig, axs = plt.subplots(2, 3, figsize=(7.4, 3.5), sharex=True, sharey=True)
    for i, ini in enumerate(INICIAIS):
        for j, h in enumerate(NIVEIS):
            ax = axs[i, j]
            p = painel_nivel(ini, h, caso)
            hs = mm.fmt(h, 4).rstrip('0')
            if p is None:
                ax.text(0.5, 0.5, 'em execução', ha='center', va='center', transform=ax.transAxes)
                ax.set_title(f'h = {hs} m', fontsize=7)
                ax.set_aspect('equal')
                ax.set_xlim(z[0], z[1])
                ax.set_ylim(z[2], z[3])
                continue
            tit = (f'h = {hs} m: {len(p["tri"])} tri., Γ = {mm.fmt(p["gamma"])}\n'
                   f'({mm.fmt(p["e"], 2, sinal=True)}% da foco125), {mm.fmt(p["cpu"], 1)} s')
            mm.desenhar(ax, p['xy'], p['tri'], p['banda'], zoom=z, s_banda=1.3, lw=0.16, titulo=tit,
                        fs_titulo=7)
            t1, tm, tr = p['cpu_partes']
            print(f'  {caso} {ini} h={h}: tri={len(p["tri"])} Gamma={p["gamma"]:.6f} '
                  f'Gamma_foco125={p["g125"]:.6f} e={p["e"]:+.3f}% cpu_cadeia={p["cpu"]:.3f} s '
                  f'(etapa1={t1} malha={tm} rodada={tr}) banda={len(p["banda"][0])} pts [{p["fonte"]}]')
            linhas_csv.append(dict(caso=caso, inicial=ini, h=h, triangulos=len(p['tri']),
                                   gamma=f'{p["gamma"]:.9g}', gamma_foco125=f'{p["g125"]:.9g}',
                                   dif_pct=f'{p["e"]:.4f}', cpu_cadeia_s=f'{p["cpu"]:.3f}',
                                   cpu_etapa1_s='' if t1 is None else t1,
                                   cpu_malha_s='' if tm is None else tm, cpu_rodada_s=tr,
                                   banda_pts=len(p['banda'][0]), fonte=p['fonte']))
        axs[i, 0].set_ylabel(f'{ini}\ny (m)', fontsize=8)
    legenda_banda(axs[0, -1], [('banda da rodada', mm.COR_BANDA, 'o', 3.0)], loc='upper right', fs=6.5)
    for ax in axs[-1]:
        ax.set_xlabel('x (m)')
    fig.tight_layout(pad=0.3, h_pad=0.4, w_pad=0.5)
    salvar(fig, f'niv_niveis_{caso}.png')


def niveis():
    fig_malhas_iniciais()
    linhas = []
    for caso in ['det', 'r00005', 'r00003', 'r00007']:
        fig_niveis(caso, linhas)
    if linhas:
        caminho = os.path.join(DADOS, 'malhas_niveis_paineis.csv')
        with open(caminho, 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=list(linhas[0]))
            w.writeheader()
            w.writerows(linhas)
        print(caminho)


def tabela_malhas_comparacao():
    """dados/malhas_comparacao.csv: estatisticas das malhas fixas e das refinadas (regeneradas) da comparacao."""
    linhas = []
    for nome in H_FIXAS:
        est = mm.estatisticas(*mm.ler_msh(os.path.join(MALHAS, f'ceron_h5_a45_{nome}.msh')))
        linhas.append(dict(malha=nome, caso='', **est, banda_etapa1_pts='', gamma=''))
    for caso in CASOS:
        d = os.path.join(COMP, 'refinada', caso)
        est = mm.estatisticas(*malha_refinada_comp(caso))
        linhas.append(dict(malha='refinada', caso=caso, **est,
                           banda_etapa1_pts=len(mm.pontos_banda_geo(os.path.join(d, 'malha_refinada.geo'))),
                           gamma=kv(os.path.join(d, 'resultado.txt'))['FS_min']))
    caminho = os.path.join(DADOS, 'malhas_comparacao.csv')
    with open(caminho, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(linhas[0]))
        w.writeheader()
        for r in linhas:
            w.writerow({k: (f'{v:.6g}' if isinstance(v, float) else v) for k, v in r.items()})
    print(caminho)


def comparacao():
    fig_malhas_fixas()
    fig_metodo_det()
    fig_refinadas()
    tabela_malhas_comparacao()


if __name__ == '__main__':
    qual = sys.argv[1] if len(sys.argv) > 1 else 'todas'
    if qual in ('comparacao', 'todas'):
        comparacao()
    if qual in ('niveis', 'todas'):
        niveis()
