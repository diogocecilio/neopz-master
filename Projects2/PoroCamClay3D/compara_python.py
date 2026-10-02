"""
compara_python.py
=================

Compara os resultados de AterroCamClay e TriaxialAbaqusCamClay (NeoPZ, C++) com o FE u-p de referência em
Python (python/fe3d_up.py, aterro_elementos.py, triaxial_abaqus.py, que usam
../PlasticityTestsCamClay/python/modified_cam_clay.py) e com os dados digitalizados do FLAC3D e do Abaqus.

Uso (no diretório onde os executáveis gravaram os CSV):
    ./AterroCamClay && ./TriaxialAbaqusCamClay
    python3 <fonte>/Projects2/PoroCamClay3D/compara_python.py [diretório_dos_csv]

Roda o código Python para os mesmos casos (alguns minutos) e, se o matplotlib estiver disponível, grava
comparacao_aterro.png e comparacao_triaxial.png.
"""
import csv
import json
import os
import sys

import numpy as np

AQUI = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.path.join(AQUI, 'python'))
sys.path.insert(0, os.path.join(AQUI, '..', 'PlasticityTestsCamClay', 'python'))
import aterro_elementos as ae  # noqa: E402
import triaxial_abaqus as ta  # noqa: E402

PASTA = sys.argv[1] if len(sys.argv) > 1 else '.'


def ler(nome):
    caminho = os.path.join(PASTA, nome)
    if not os.path.exists(caminho):
        return None
    linhas = list(csv.DictReader(open(caminho)))
    return {k: np.array([float(r[k]) for r in linhas]) for k in linhas[0]}


# ---------------------------------------------------------------------------------------- aterro
print('Aterro sobre fundação Cam-Clay (FLAC3D): NeoPZ x Python (todos os instantes)')
ATERRO = {}
for et in ('hex8', 'hex20'):
    c = ler(f'aterro_{et}_historico.csv')
    if c is None:
        continue
    r = ae.run(et, init='gp', verbose=False)
    h = r['hist']
    ATERRO[et] = (c, h, r['n_load'])
    dif = {k: np.abs(c[k] - h[k]).max() for k in ('uz0', 'uz2', 'uz4', 'uz6', 'pp1', 'pp2')}
    print(f'   {et:6s} ' + '  '.join(f'max|Δ{k}| = {v:.1e}' for k, v in dif.items()))

# ---------------------------------------------------------------------------------------- triaxial
print('Abaqus 1.15.2 (malha 2x2x4, 150 passos): NeoPZ x Python (todos os passos)')
TRI = {}
for placa, platen in (('lisa', 'smooth'), ('rugosa', 'rough')):
    for et in ('hex20', 'hex20r'):
        c = ler(f'triaxial_{placa}_{et}_224.csv')
        if c is None:
            continue
        h = ta.run(platen, 'porous', 2, 2, 4, 150, et, verbose=False)['hist']
        TRI[(placa, et)] = (c, h)
        print(f'   {placa:6s} {et:6s} max|Δq_A| = {np.abs(c["q_A"] - h["q"]).max():.1e}  '
              f'max|Δp_A| = {np.abs(c["p_A"] - h["p"]).max():.1e}  '
              f'max|Δσ_a| = {np.abs(c["sigma_a_placa"] - h["sa"]).max():.1e} kPa;  '
              f'iterações de Newton iguais: {bool((c["iteracoes"] == h["it"]).all())}')

# ---------------------------------------------------------------------------------------- figuras
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except ImportError:
    print('matplotlib não disponível: figuras não geradas')
    sys.exit(0)

C_PZ, C_PY, C_REF = '#eb6834', '#2a78d6', '#0b0b0b'
if ATERRO:
    flac = json.load(open(os.path.join(AQUI, 'python', 'flac_historicos_digitalizados.json')))
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.6))
    for et, (c, h, nl) in ATERRO.items():
        ls = '-' if et == 'hex20' else ':'
        t = c['t'][nl:].copy(); t[0] = 10.0
        for k, cor in zip(('uz0', 'uz4'), ('#0b0b0b', '#1baf7a')):
            axs[0].semilogx(t, -h[k][nl:], ls, color=C_PY, lw=1.4)
            axs[0].semilogx(t[::2], -c[k][nl::2], 'o', mfc='none', mec=C_PZ, ms=4)
        for k in ('pp1', 'pp2'):
            axs[1].semilogx(t, h[k][nl:], ls, color=C_PY, lw=1.4)
            axs[1].semilogx(t[::2], c[k][nl::2], 'o', mfc='none', mec=C_PZ, ms=4)
    for k, ax in (('uz_x0', axs[0]), ('uz_x4', axs[0]), ('pp1', axs[1]), ('pp2', axs[1])):
        tt, vv = np.array(flac[k][0]), np.array(flac[k][1]); o = np.argsort(tt)
        ax.semilogx(tt[o], vv[o], '--', color=C_REF, lw=1.0)
    axs[0].set_ylabel('recalque (m), x = 0 e x = 4'); axs[1].set_ylabel('poropressão de zona (kPa), pp1 e pp2')
    for ax in axs:
        ax.set_xlabel('tempo (s)  [t = 10 s: fim da etapa não drenada]'); ax.grid(alpha=.3)
    axs[0].plot([], [], '-', color=C_PY, label='Python hex20 (traço: hex8)')
    axs[0].plot([], [], 'o', mfc='none', mec=C_PZ, label='NeoPZ')
    axs[0].plot([], [], '--', color=C_REF, label='FLAC3D (digitalizado)')
    axs[0].legend(fontsize=8.5)
    fig.suptitle('Aterro sobre fundação Cam-Clay (exemplo do FLAC3D)')
    fig.tight_layout(); fig.savefig(os.path.join(PASTA, 'comparacao_aterro.png'), dpi=120); plt.close(fig)

if TRI:
    abq = ta.load_abaqus()
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.6))
    for ax, placa, chave in ((axs[0], 'lisa', 'lisa'), (axs[1], 'rugosa', 'rugosa')):
        a = np.array(abq[chave]['qd'])
        ax.plot(a[:, 0], a[:, 1], 'o', color=C_REF, ms=6, label='Abaqus (digitalizado)')
        for et, ls in (('hex20', '-'), ('hex20r', '--')):
            if (placa, et) not in TRI:
                continue
            c, h = TRI[(placa, et)]
            ax.plot(h['dh'], h['q'], ls, color=C_PY, lw=1.4, label=f'Python {et}')
            ax.plot(c['dH'][::6], c['q_A'][::6], 'o', mfc='none', mec=C_PZ, ms=4, label=f'NeoPZ {et}')
        ax.set_title(f'placa {placa}'); ax.set_xlabel('δ/H'); ax.set_ylabel('q no ponto A (kPa)')
        ax.grid(alpha=.3); ax.legend(fontsize=8.5)
    fig.suptitle('Abaqus 1.15.2: adensamento de um corpo de prova triaxial (hex20, malha 2x2x4)')
    fig.tight_layout(); fig.savefig(os.path.join(PASTA, 'comparacao_triaxial.png'), dpi=120); plt.close(fig)
print('figuras: comparacao_aterro.png, comparacao_triaxial.png')
