"""
compara_python.py
=================

Compara os resultados do NeoPZ (C++) com o FE u-p de referência em Python (python/fe3d_up.py,
aterro_elementos.py, triaxial_abaqus.py, que usam ../PlasticityTestsCamClay/python/modified_cam_clay.py) e com
os dados digitalizados do FLAC3D e do Abaqus:

  * versão nativa (AterroCamClay, TriaxialAbaqusCamClay: TPZMultiphysicsCompMesh + TPZMatPoroElastoPlastic3DMem):
    hex8 é o mesmo elemento do Python (diferença ~0); o hex20 nativo tem u Q2 (27 funções), o do Python é o
    serendipity de 20 nós, então a diferença é a do elemento;
  * versão serendipity (AterroCamClaySerendipity, TriaxialAbaqusCamClaySerendipity): os mesmos elementos do
    Python, inclusive o hex20 serendipity e o hex20r (diferença ~0).

Uso (no diretório onde os executáveis gravaram os CSV):
    ./AterroCamClay && ./TriaxialAbaqusCamClay
    ./AterroCamClaySerendipity && ./TriaxialAbaqusCamClaySerendipity     # opcional
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
VERSOES = (('', 'nativo'), ('serendipity_', 'serendipity'))
print('Aterro sobre fundação Cam-Clay (FLAC3D): NeoPZ x Python (todos os instantes)')
ATERRO, PY_ATERRO = {}, {}
for pref, versao in VERSOES:
    for et in ('hex8', 'hex20'):
        c = ler(f'{pref}aterro_{et}_historico.csv')
        if c is None:
            continue
        if et not in PY_ATERRO:
            r = ae.run(et, init='gp', verbose=False)
            PY_ATERRO[et] = (r['hist'], r['n_load'])
        h, nl = PY_ATERRO[et]
        ATERRO[(versao, et)] = (c, h, nl)
        dif = {k: np.abs(c[k] - h[k]).max() for k in ('uz0', 'uz2', 'uz4', 'uz6', 'pp1', 'pp2')}
        obs = '  (u Q2 x serendipity do Python)' if (versao, et) == ('nativo', 'hex20') else ''
        print(f'   {versao:11s} {et:6s} ' + '  '.join(f'max|Δ{k}| = {v:.1e}' for k, v in dif.items()) + obs)

# ---------------------------------------------------------------------------------------- triaxial
print('Abaqus 1.15.2 (malha 2x2x4, 150 passos): NeoPZ x Python (todos os passos)')
TRI, PY_TRI = {}, {}
for pref, versao in VERSOES:
    for placa, platen in (('lisa', 'smooth'), ('rugosa', 'rough')):
        for et in ('hex20', 'hex20r'):
            c = ler(f'{pref}triaxial_{placa}_{et}_224.csv')
            if c is None:
                continue
            if (placa, et) not in PY_TRI:
                PY_TRI[(placa, et)] = ta.run(platen, 'porous', 2, 2, 4, 150, et, verbose=False)['hist']
            h = PY_TRI[(placa, et)]
            TRI[(versao, placa, et)] = (c, h)
            obs = '  (u Q2 x serendipity do Python)' if versao == 'nativo' else ''
            print(f'   {versao:11s} {placa:6s} {et:6s} max|Δq_A| = {np.abs(c["q_A"] - h["q"]).max():.1e}  '
                  f'max|Δp_A| = {np.abs(c["p_A"] - h["p"]).max():.1e}  '
                  f'max|Δσ_a| = {np.abs(c["sigma_a_placa"] - h["sa"]).max():.1e} kPa;  '
                  f'iterações de Newton iguais: {bool((c["iteracoes"] == h["it"]).all())}' + obs)

# ---------------------------------------------------------------------------------------- figuras
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except ImportError:
    print('matplotlib não disponível: figuras não geradas')
    sys.exit(0)

C_PZ, C_SER, C_PY, C_REF = '#eb6834', '#7a3fb0', '#2a78d6', '#0b0b0b'
MARCA = {'nativo': dict(marker='o', mfc='none', mec=C_PZ, ms=4, ls='none'),
         'serendipity': dict(marker='x', color=C_SER, ms=4, ls='none')}
if ATERRO:
    flac = json.load(open(os.path.join(AQUI, 'python', 'flac_historicos_digitalizados.json')))
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.6))
    for et, (h, nl) in PY_ATERRO.items():
        ls = '-' if et == 'hex20' else ':'
        t = np.array(h['t'][nl:]); t[0] = 10.0
        for k in ('uz0', 'uz4'):
            axs[0].semilogx(t, -np.array(h[k][nl:]), ls, color=C_PY, lw=1.4)
        for k in ('pp1', 'pp2'):
            axs[1].semilogx(t, np.array(h[k][nl:]), ls, color=C_PY, lw=1.4)
    for (versao, et), (c, h, nl) in ATERRO.items():
        t = c['t'][nl:].copy(); t[0] = 10.0
        for k in ('uz0', 'uz4'):
            axs[0].semilogx(t[::2], -c[k][nl::2], **MARCA[versao])
        for k in ('pp1', 'pp2'):
            axs[1].semilogx(t[::2], c[k][nl::2], **MARCA[versao])
    for k, ax in (('uz_x0', axs[0]), ('uz_x4', axs[0]), ('pp1', axs[1]), ('pp2', axs[1])):
        tt, vv = np.array(flac[k][0]), np.array(flac[k][1]); o = np.argsort(tt)
        ax.semilogx(tt[o], vv[o], '--', color=C_REF, lw=1.0)
    axs[0].set_ylabel('recalque (m), x = 0 e x = 4'); axs[1].set_ylabel('poropressão de zona (kPa), pp1 e pp2')
    for ax in axs:
        ax.set_xlabel('tempo (s)  [t = 10 s: fim da etapa não drenada]'); ax.grid(alpha=.3)
    axs[0].plot([], [], '-', color=C_PY, label='Python hex20 serendipity (pontilhado: hex8)')
    axs[0].plot([], [], **MARCA['nativo'], label='NeoPZ nativo (hex20 Q2 e hex8)')
    axs[0].plot([], [], **MARCA['serendipity'], label='NeoPZ serendipity')
    axs[0].plot([], [], '--', color=C_REF, label='FLAC3D (digitalizado)')
    axs[0].legend(fontsize=8.5)
    fig.suptitle('Aterro sobre fundação Cam-Clay (exemplo do FLAC3D)')
    fig.tight_layout(); fig.savefig(os.path.join(PASTA, 'comparacao_aterro.png'), dpi=120); plt.close(fig)

if TRI:
    abq = ta.load_abaqus()
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.6))
    for ax, placa in ((axs[0], 'lisa'), (axs[1], 'rugosa')):
        a = np.array(abq[placa]['qd'])
        ax.plot(a[:, 0], a[:, 1], 'o', color=C_REF, ms=6, label='Abaqus (digitalizado)')
        for et, ls in (('hex20', '-'), ('hex20r', '--')):
            if (placa, et) in PY_TRI:
                h = PY_TRI[(placa, et)]
                ax.plot(h['dh'], h['q'], ls, color=C_PY, lw=1.4, label=f'Python {et} serendipity')
        for (versao, pl, et), (c, h) in TRI.items():
            if pl == placa:
                rot = f'NeoPZ nativo {et} (Q2)' if versao == 'nativo' else f'NeoPZ serendipity {et}'
                ax.plot(c['dH'][::6], c['q_A'][::6], **MARCA[versao], label=rot)
        ax.set_title(f'placa {placa}'); ax.set_xlabel('δ/H'); ax.set_ylabel('q no ponto A (kPa)')
        ax.grid(alpha=.3); ax.legend(fontsize=8.5)
    fig.suptitle('Abaqus 1.15.2: adensamento de um corpo de prova triaxial (malha 2x2x4)')
    fig.tight_layout(); fig.savefig(os.path.join(PASTA, 'comparacao_triaxial.png'), dpi=120); plt.close(fig)
print('figuras: comparacao_aterro.png, comparacao_triaxial.png')
