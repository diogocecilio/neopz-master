"""
compara_python.py
=================

Compara os resultados do executável PlasticityTestsCamClay (TPZModifiedCamClay, C++) com o código
Python de referência (python/modified_cam_clay.py e python/poro_camclay_fem.py).

Uso (no diretório onde o executável gravou os CSV):
    ./PlasticityTestsCamClay
    python3 <fonte>/Projects2/PlasticityTestsCamClay/compara_python.py [diretório_dos_csv]

Se o matplotlib estiver disponível, grava comparacao_triaxial.png e comparacao_itasca.png.
"""
import os
import sys

import numpy as np

AQUI = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.path.join(AQUI, 'python'))
import modified_cam_clay as mcc  # noqa: E402
import poro_camclay_fem as fe  # noqa: E402

PASTA = sys.argv[1] if len(sys.argv) > 1 else '.'
XX, XY, XZ, YY, YZ, ZZ = range(6)
TAB81 = dict(M=1.2, lam=0.077, kap=0.0066, N=1.788, nu=0.3, G=20000.)


def ler(nome):
    return np.genfromtxt(os.path.join(PASTA, nome), delimiter=',', names=True, dtype=None, encoding='utf-8')


def relerr(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    return np.abs(a - b).max() / max(np.abs(b).max(), 1e-300)


# ------------------------------------------------------------------------------ 1) estados pontuais
print('1) Estados de ponto material (return mapping e tangente consistente)')
E = ler('estados_pontuais.csv')
PX = mcc.mcc_parameters(**TAB81, v0=1.7, p0=150., pc0=200., pt=20., beta=0.6, elasticity='pressure_dependent',
                        shear='constant_nu', sigma0=[-150., 10., -5., -140., 3., -160.])
PH = mcc.mcc_parameters(M=1.0, lam=0.174, kap=0.026, v0=2.08, pc0=2 * 58.3, p0=100., elasticity='pressure_dependent',
                        shear='hypo_nu', nu=0.3, sigma0=[-100., 0., 0., -100., 0., -100.])
sig_n = np.array([-95., 4., -2., -105., 1., -120.])
eps_n = np.array([0.001, 0.0005, 0., 0.0008, -0.0003, -0.002])
epsp_n = np.array([0.0001, 0., 0.0002, -0.0001, 0., 0.0003])
deps = [np.array([0.0002, 0.0001, 0., 0.0002, 0., -0.0005]), np.array([0.003, 0.002, -0.001, 0.003, 0.0005, -0.009]),
        np.array([-0.002, 0., 0., -0.002, 0., -0.002])]
worst = dict(stress=0., Dep=0., alpha=0., epsp=0.)
for row in E:
    k = int(row['k'])
    if row['caso'] == 'px':
        e = np.array([np.sin(1.3 * k + i) for i in range(6)]) * 0.005
        ep = np.array([np.cos(0.7 * k + 2 * i) for i in range(6)]) * 0.001
        r = mcc.return_mapping(PX, e, ep, 0.0025 * (1 + np.sin(k)))
    else:
        r = mcc.return_mapping(PH, eps_n + deps[k], epsp_n, 0.001, sig_n=sig_n, eps_n=eps_n)
    s = np.array([row[f's{i}'] for i in range(6)])
    D = np.array([row[f'D{i}'] for i in range(36)]).reshape(6, 6)
    epsp = np.array([row[f'ep{i}'] for i in range(6)])
    assert bool(row['plastic']) == bool(r['plastic'])
    worst['stress'] = max(worst['stress'], relerr(s, r['stress']))
    worst['Dep'] = max(worst['Dep'], relerr(D, r['Dep']))
    worst['alpha'] = max(worst['alpha'], abs(row['alpha'] - r['alpha']) / max(abs(r['alpha']), 1e-12))
    worst['epsp'] = max(worst['epsp'], relerr(epsp, r['plastic_strain']))
print(f'   {len(E)} estados: maior diferença relativa C++ x Python: ' +
      ', '.join(f'{k} = {v:.2e}' for k, v in worst.items()))

# ------------------------------------------------------------------------------ 2) triaxiais drenados
print('2) Ensaios triaxiais drenados (ponto material)')
CASOS = {
    'Fig8.5': mcc.mcc_parameters(**TAB81, v0=1.7, p0=200., pc0=200., shear='constant_nu'),
    'Fig8.6': mcc.mcc_parameters(**TAB81, v0=1.7, p0=200., pc0=200., shear='constant_G'),
    'Fig8.7': mcc.mcc_parameters(**TAB81, v0=1.7, p0=100., pc0=200., shear='constant_nu'),
    'Fig8.8': mcc.mcc_parameters(**TAB81, v0=1.7, p0=100., pc0=500., shear='constant_nu'),
    'OCR5_K_p': mcc.mcc_parameters(**TAB81, v0=1.7, p0=100., pc0=500., elasticity='pressure_dependent',
                                   shear='constant_nu'),
}
TRI = {}
for nome, P in CASOS.items():
    c = ler(f'triaxial_{nome}.csv')
    num = mcc.triaxial_drained(P, 0.2, 400)
    cpp = np.c_[c['ea'], c['p'], c['q'], c['ev'], c['eq'], c['sa'], c['pc']]
    TRI[nome] = (cpp, num, ler(f'triaxial_{nome}_analitico.csv'))
    print(f'   {nome:9s} max|Δq| = {np.abs(cpp[:, 2] - num[:, 2]).max():.2e} kPa   '
          f'max|Δεv| = {np.abs(cpp[:, 3] - num[:, 3]).max():.2e}')

# ------------------------------------------------------------------------------ 3) Itasca (FE)
print('3) Benchmark Itasca (FE com 1 hexaedro): C++ (TPZMatElastoPlastic) x Python (poro_camclay_fem)')
ITA = {}
for nome, R, drenado in [('drenado_R1.6', 1.6, True), ('drenado_R8', 8.0, True),
                         ('nao_drenado_R1.6', 1.6, False), ('nao_drenado_R8', 8.0, False)]:
    P = mcc.mcc_parameters(M=1.02, lam=0.2, kap=0.05, N=3.32, pc0=R * 5., p0=5., elasticity='pressure_dependent',
                           shear='constant_G', G=250.)
    if drenado:
        h, _ = fe.triaxial_test(P, drained=True, ea_max=0.5, nsteps=500)
    else:
        n0 = (P['v0'] - 1) / P['v0']
        h, _ = fe.triaxial_test(P, drained=False, ea_max=0.1, nsteps=400, biot_modulus=2e4 / n0)
    c = ler(f'itasca_{nome}.csv')
    ITA[nome] = (c, h)
    print(f'   {nome:17s} max|Δp\'| = {np.abs(c["p"] - h["p"]).max():.2e}  max|Δq| = {np.abs(c["q"] - h["q"]).max():.2e}'
          f'  max|Δv| = {np.abs(c["v"] - h["v"]).max():.2e}  max|Δu| = {np.abs(c["u"] - h["u"]).max():.2e}')

# ------------------------------------------------------------------------------ 4) hipoelástico
print('4) Hipoelástico (Abaqus 1.15.2): FE C++ x ponto material Python')
eps = np.zeros(6); epsp = np.zeros(6); al = 0.0; sig = PH['sigma0'].copy(); ratio = 0.3
dea = 0.1 / 200; sc = PH['sigma0'][XX]; hyp = [(0.0, 100.0, 0.0)]
for k in range(200):
    e = eps.copy(); e[ZZ] -= dea; e[XX] -= ratio * dea; e[YY] -= ratio * dea
    for it in range(30):
        r = mcc.return_mapping(PH, e, epsp, al, sig_n=sig, eps_n=eps)
        res = np.array([r['stress'][XX] - sc, r['stress'][YY] - sc])
        if np.abs(res).max() <= 1e-10 * abs(sc):
            break
        d = np.linalg.solve(r['Dep'][np.ix_([XX, YY], [XX, YY])], -res)
        e[XX] += d[0]; e[YY] += d[1]
    ratio = (e[XX] - eps[XX]) / (-dea)
    eps, epsp, al, sig = e, r['plastic_strain'], r['alpha'], r['stress']
    p, q = mcc.invariants(sig)
    hyp.append((-eps[ZZ], -p, q))
hyp = np.array(hyp)
c = ler('hypo_fe.csv')
print(f'   max|Δp\'| = {np.abs(c["p"] - hyp[:, 1]).max():.2e} kPa   max|Δq| = {np.abs(c["q"] - hyp[:, 2]).max():.2e} kPa')

# ------------------------------------------------------------------------------ figuras
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except ImportError:
    print('matplotlib não disponível: figuras não geradas')
    sys.exit(0)

C_CPP, C_PY, C_ANA = '#eb6834', '#2a78d6', '#0b0b0b'
fig, ax = plt.subplots(2, len(TRI), figsize=(4 * len(TRI), 7))
for j, (nome, (cpp, num, ana)) in enumerate(TRI.items()):
    for i, (iy, lab) in enumerate([(2, 'q (kPa)'), (3, 'εv')]):
        a = ax[i, j]
        a.plot(ana['ea'], ana['q'] if iy == 2 else ana['ev'], '-', color=C_ANA, lw=1.4, label='analítica')
        a.plot(num[:, 0], num[:, iy], '--', color=C_PY, lw=1.4, label='Python')
        a.plot(cpp[::10, 0], cpp[::10, iy], 'o', mfc='none', mec=C_CPP, ms=5, label='NeoPZ (C++)')
        a.set_xlim(0, 0.2); a.set_xlabel('εa'); a.set_ylabel(lab); a.grid(alpha=.3)
        if i == 0:
            a.set_title(nome)
ax[0, 0].legend(fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(PASTA, 'comparacao_triaxial.png'), dpi=120); plt.close(fig)

fig, ax = plt.subplots(1, 3, figsize=(13, 4))
for nome, (c, h) in ITA.items():
    ax[0].plot(h['p'], h['q'], '-', color=C_PY, lw=1.2)
    ax[0].plot(c['p'][::10], c['q'][::10], 'o', mfc='none', mec=C_CPP, ms=4)
    ax[1].plot(h['ea'], h['q'], '-', color=C_PY, lw=1.2)
    ax[1].plot(c['ea'][::10], c['q'][::10], 'o', mfc='none', mec=C_CPP, ms=4)
    ax[2].plot(h['ea'], h['v'], '-', color=C_PY, lw=1.2)
    ax[2].plot(c['ea'][::10], c['v'][::10], 'o', mfc='none', mec=C_CPP, ms=4)
for a, (xl, yl) in zip(ax, [("p' (kPa)", 'q (kPa)'), ('εa', 'q (kPa)'), ('εa', 'v')]):
    a.set_xlabel(xl); a.set_ylabel(yl); a.grid(alpha=.3)
ax[0].plot([], [], '-', color=C_PY, label='Python (poro_camclay_fem)')
ax[0].plot([], [], 'o', mfc='none', mec=C_CPP, label='NeoPZ (TPZMatElastoPlastic)')
ax[0].legend(fontsize=8); fig.suptitle('Benchmark Itasca: triaxial drenado e não drenado (R = 1.6 e 8)')
fig.tight_layout(); fig.savefig(os.path.join(PASTA, 'comparacao_itasca.png'), dpi=120); plt.close(fig)
print('figuras: comparacao_triaxial.png, comparacao_itasca.png')
