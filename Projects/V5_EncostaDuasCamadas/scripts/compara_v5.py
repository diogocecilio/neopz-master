#!/usr/bin/env python3
# V5_EncostaDuasCamadas: comparacao com Abbasov et al. (2024) de duas execucoes do main.cpp: CASO=tabela1 (mapas das
# Figs. 4 a 6) e CASO=tabela4 (Fig. 14 e, com NU_LISTA, C_LISTA e PHI_LISTA, a Fig. 13). Imprime a tabela
# "artigo x nosso" (markdown) e grava as figuras: perfis nas verticais x_zona = 2, 10, 17 e 20,5 m contra os perfis
# digitalizados (perfis_v5.png), barras da Fig. 14 (fig14_v5.png), efeitos marginais da Fig. 13 (fig13_v5.png) e os
# mapas (mapas_tabela1.png, com scripts/mapas_v5.py).
#
# Uso: python3 compara_v5.py --tab1 <pasta tabela1> --tab4 <pasta tabela4> [--figuras ../figuras] [--csv tabela.csv]
import argparse
import os
import subprocess
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

AQUI = os.path.dirname(os.path.abspath(__file__))
DADOS = os.path.join(AQUI, '..', 'dados')
SECOES = {2.: 17., 10.: 25., 17.: 32., 20.5: 35.5}  # x_zona -> x


def le(arq):
    with open(arq) as f:
        pula = 1 if f.readline().startswith('#') else 0
    return np.genfromtxt(arq, delimiter=',', names=True, skip_header=pula, dtype=None, encoding=None)


def campo(r, nome):
    """valor de uma coluna (nan se o arquivo nao a tem: execucoes de versoes anteriores)"""
    return float(r[nome]) if nome in r.dtype.names else np.nan


def cruza_desc(prof, v, alvo, pmax=2.):
    """frente de umedecimento na camada superior (ate pmax): a primeira profundidade, descendo, em que v cai abaixo de
    alvo partindo de v >= alvo; pmax se a camada toda passa de alvo; nan se ja comeca abaixo"""
    m = (prof <= pmax) & np.isfinite(v)
    p, w = prof[m], v[m]
    for k in range(1, len(w)):
        if w[k - 1] >= alvo > w[k]:
            return p[k - 1] + (p[k] - p[k - 1]) * (w[k - 1] - alvo) / (w[k - 1] - w[k])
    return pmax if len(w) and w.min() >= alvo else np.nan


def perfil_artigo(dig, var, t, xz):
    s = dig[(dig['variavel'] == var) & (dig['t_h'] == t) & (dig['x_zona'] == xz)]
    o = np.argsort(s['prof'])
    return s['prof'][o], s['valor'][o]


def perfil_nosso(per, col, t, x):
    s = per[(np.abs(per['t_h'] - t) < 1e-9) & (np.abs(per['x'] - x) < 1e-9)]
    o = np.argsort(s['prof'])
    return s['prof'][o], s[col][o]


def rms_mapa(pasta, arqdig, var, t, margem=0.15):
    """diferenca media e RMS entre o nosso mapa (campos_t<t>.csv: theta e tensao de succao lineares no elemento, LFS
    do elemento limitado a 3) e o mapa digitalizado, nos pontos a mais de 'margem' da interface"""
    import matplotlib.tri as mtri
    arq = os.path.join(pasta, f'campos_t{t:g}.csv')
    if not os.path.exists(arq):
        return np.nan, np.nan, 0
    nos = le(os.path.join(pasta, 'nos.csv'))
    els = {int(e['el']): e for e in le(os.path.join(pasta, 'elementos.csv')) if e['zona'] == 1}
    c = [r for r in le(arq) if int(r['el']) in els]
    tri = np.array([[int(els[int(r['el'])][f'n{j}']) for j in range(3)] for r in c])
    T = mtri.Triangulation(nos['x'] - 15., nos['z'] - 17., tri)
    busca = T.get_trifinder()
    d = le(os.path.join(DADOS, arqdig))
    nome = d.dtype.names[3]
    d = d[np.abs(d['t_h'] - t) < 1e-9]
    dif = []
    for p in d:
        x, z = p['x_zona'], p['z_zona']
        zi = 6. if x <= 4. else (0. if x >= 19. else 6. - 6. * (x - 4.) / 15.)
        if abs(z - zi) < margem:
            continue
        k = int(busca(x, z))
        if k < 0:
            continue
        r = c[k]
        if var == 'lfs':
            v = min(r['lfs'], 3.)
        else:
            # coordenadas baricentricas no triangulo k
            xs, zs = T.x[tri[k]], T.y[tri[k]]
            det = (xs[1] - xs[0]) * (zs[2] - zs[0]) - (xs[2] - xs[0]) * (zs[1] - zs[0])
            l1 = ((x - xs[0]) * (zs[2] - zs[0]) - (xs[2] - xs[0]) * (z - zs[0])) / det
            l2 = ((xs[1] - xs[0]) * (z - zs[0]) - (x - xs[0]) * (zs[1] - zs[0])) / det
            w = (1. - l1 - l2, l1, l2)
            v = sum(w[j] * r[f'{var}{j}'] for j in range(3))
        dif.append(v - p[nome])
    dif = np.array(dif)
    return (float(dif.mean()), float(np.sqrt((dif ** 2).mean())), len(dif)) if len(dif) else (np.nan, np.nan, 0)


def nivel_dagua(prof, sig, z_sup=8.):
    """profundidade do nivel d'agua: abaixo da interface (2 m), a primeira em que a tensao de succao cai abaixo de 0,3 kPa"""
    m = prof > 2.05
    p, s = prof[m], sig[m]
    for k in range(len(s)):
        if np.isfinite(s[k]) and s[k] < 0.3:
            return p[k]
    return np.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--tab1', required=True)
    ap.add_argument('--tab4', required=True)
    ap.add_argument('--figuras', default=os.path.join(AQUI, '..', 'figuras'))
    ap.add_argument('--csv', default='')
    ap.add_argument('--fig9', default='', help='pasta da execucao CASO=fig9 (Fig. 9, 15 h)')
    args = ap.parse_args()
    os.makedirs(args.figuras, exist_ok=True)
    dig = le(os.path.join(DADOS, 'perfis_fig4_5_6.csv'))
    per1 = le(os.path.join(args.tab1, 'perfis.csv'))
    met1 = le(os.path.join(args.tab1, 'metricas.csv'))
    met4 = le(os.path.join(args.tab4, 'metricas.csv'))
    f14 = le(os.path.join(DADOS, 'fig14_lfs.csv'))
    linhas = []  # grandeza, artigo, nosso, diferenca, fonte

    def add(g, a, n, fonte, rel=False, fmt='{:.3f}'):
        if np.isfinite(a) and np.isfinite(n):
            d = (n - a) / a * 100. if rel else n - a
            ds = f'{d:+.1f}%' if rel else ('{:+.3f}'.format(d))
        else:
            ds = '-'
        linhas.append((g, fmt.format(a) if np.isfinite(a) else '-', fmt.format(n) if np.isfinite(n) else '-', ds, fonte))

    # ---- fluxo (tabela1, Figs. 4 e 5), vertical da crista x_zona = 2 m
    pa, ta = perfil_artigo(dig, 'theta', 0., 2.)
    pn, tn = perfil_nosso(per1, 'theta', 0., 17.)
    add('theta camada superior, crista, 1 m, 0 h', np.interp(1., pa, ta), np.interp(1., pn, tn), 'Fig. 4 (digit.)')
    add('theta topo da camada inferior, crista, 2,1 m, 0 h', np.interp(2.1, pa, ta), np.interp(2.1, pn, tn), 'Fig. 4 (digit.)')
    add('theta camada inferior, crista, 3 m, 0 h', np.interp(3., pa, ta), np.interp(3., pn, tn), 'Fig. 4 (digit.)')
    pa, sa = perfil_artigo(dig, 'sigs_kPa', 0., 2.)
    pn, sn = perfil_nosso(per1, 'sigs_kPa', 0., 17.)
    add('tensao de succao camada superior, crista, 1 m, 0 h (kPa)', np.interp(1., pa, sa), np.interp(1., pn, sn), 'Fig. 5 (digit.)')
    add('tensao de succao topo da camada inferior, crista, 2,1 m, 0 h (kPa)', np.interp(2.1, pa, sa), np.interp(2.1, pn, sn), 'Fig. 5 (digit.)')
    add('profundidade do nivel d\'agua na crista, 0 h (m)', nivel_dagua(pa, sa), nivel_dagua(pn, sn), 'Fig. 5 (digit.)')
    for t in (7.5, 15.):
        pa, ta = perfil_artigo(dig, 'theta', t, 2.)
        pn, tn = perfil_nosso(per1, 'theta', t, 17.)
        add(f'frente de umedecimento na crista (theta = 0,32), {t:g} h (m)', cruza_desc(pa, ta, 0.32),
            cruza_desc(pn, tn, 0.32), 'Fig. 4 (digit.)')
        pa, ta = perfil_artigo(dig, 'theta', t, 10.)
        pn, tn = perfil_nosso(per1, 'theta', t, 25.)
        add(f'frente de umedecimento no meio do talude, {t:g} h (m)', cruza_desc(pa, ta, 0.32),
            cruza_desc(pn, tn, 0.32), 'Fig. 4 (digit.)')
    pa, sa = perfil_artigo(dig, 'sigs_kPa', 15., 2.)
    pn, sn = perfil_nosso(per1, 'sigs_kPa', 15., 17.)
    add('tensao de succao na zona umedecida, crista, 0,5 m, 15 h (kPa)', np.interp(0.5, pa, sa), np.interp(0.5, pn, sn), 'Fig. 5 (digit.)')
    # LFS (tabela1, Fig. 6): na crista, a 7,6 m (fundo da zona), e no meio do talude a 1 m
    for t in (0., 15.):
        pa, la = perfil_artigo(dig, 'lfs', t, 2.)
        pn, ln = perfil_nosso(per1, 'lfs', t, 17.)
        add(f'LFS crista, 7,6 m de profundidade, {t:g} h', np.interp(7.6, pa, la), np.interp(7.6, pn, ln), 'Fig. 6 (digit.)')
        add(f'LFS crista, 4 m de profundidade, {t:g} h', np.interp(4., pa, la), np.interp(4., pn, ln), 'Fig. 6 (digit.)')
        pa, la = perfil_artigo(dig, 'lfs', t, 10.)
        pn, ln = perfil_nosso(per1, 'lfs', t, 25.)
        add(f'LFS meio do talude, 1 m, {t:g} h', np.interp(1., pa, la), min(np.interp(1., pn, ln), 3.), 'Fig. 6 (digit., <= 3)')
        add(f'LFS meio do talude, 4 m, {t:g} h', np.interp(4., pa, la), min(np.interp(4., pn, ln), 3.), 'Fig. 6 (digit., <= 3)')
    # mapas inteiros (tabela1): diferenca RMS nos pontos digitalizados
    for t in (0., 7.5, 15.):
        for var, arq, rot in (('theta', 'fig4_theta.csv', 'theta'), ('sigs', 'fig5_sigs.csv', 'tensao de succao (kPa)'),
                              ('lfs', 'fig6_lfs.csv', 'LFS limitado a 3')):
            med, rms, n = rms_mapa(args.tab1, arq, var, t)
            linhas.append((f'mapa de {rot}, {t:g} h: diferenca media / RMS ({n} pontos)', '-',
                           f'{med:+.3f} / {rms:.3f}', '-', arq.split('_')[0].replace('fig', 'Fig. ') + ' (digit.)'))
    # Fig. 9 (caso fig9, 15 h)
    if args.fig9:
        dig9 = dig[dig['variavel'] == 'lfs_fig9']
        per9 = le(os.path.join(args.fig9, 'perfis.csv'))
        for xz, x, rot in ((2., 17., 'crista'), (10., 25., 'meio do talude')):
            pa, la = perfil_artigo(dig, 'lfs_fig9', 15., xz)
            pn, ln = perfil_nosso(per9, 'lfs', 15., x)
            for pr in (1., 4., 7.5 if xz == 2. else 5.):
                if pr <= pn.max():
                    add(f'Fig. 9: LFS {rot}, {pr:g} m, 15 h', np.interp(pr, pa, la), min(np.interp(pr, pn, ln), 3.),
                        'Fig. 9 (digit., <= 3)')
        med, rms, n = rms_mapa(args.fig9, 'fig9_lfs.csv', 'lfs', 15.)
        linhas.append((f'Fig. 9: mapa do LFS limitado a 3, 15 h: diferenca media / RMS ({n} pontos)', '-',
                       f'{med:+.3f} / {rms:.3f}', '-', 'Fig. 9 (digit.)'))
        m9 = le(os.path.join(args.fig9, 'metricas.csv'))
        m9 = m9[np.abs(m9['t_h'] - 15.) < 1e-9]
        if len(m9):
            linhas.append(('Fig. 9: A(LFS < 1) e LFS medio limitado a 3, 15 h', '(mapa)',
                           f"{float(m9['A_LFS_m2'][0]):.2f} m2 / {float(m9['LFS_medio_lim3'][0]):.3f}", '-', 'caso fig9'))
    # ---- metricas (tabela4, Fig. 14)
    for r in f14:
        t = float(r['t_h'])
        m = met4[np.abs(met4['t_h'] - t) < 1e-9]
        if len(m) == 0:
            continue
        m = m[0]
        add(f'LFS medio (eq. 11), {t:g} h', r['LFS_medio_det'], m['LFS_medio'], 'Fig. 14', rel=True)
        add(f'LFS medio limitado a 3, {t:g} h', r['LFS_medio_det'], campo(m, 'LFS_medio_lim3'), 'Fig. 14', rel=True)
        add(f'A(LFS < 1), {t:g} h (m2)', r['A_LFS_det_m2'], m['A_LFS_m2'], 'Fig. 14' + (' e sec. 4.6' if t == 15 else ''),
            rel=True, fmt='{:.2f}')
    add('area da zona de investigacao (m2)', 113., float(met4['area_zona_m2'][0]), 'sec. 4.6', rel=True, fmt='{:.2f}')
    # ---- sensibilidade (tabela4, Fig. 13) a 15 h, se houver
    sen_arq = os.path.join(args.tab4, 'sensibilidade.csv')
    f13 = le(os.path.join(DADOS, 'fig13_efeitos_marginais.csv'))
    sen = le(sen_arq) if os.path.exists(sen_arq) else None
    if sen is not None:
        for par in ('nu', 'phi', 'c'):
            ref = f13[f13['parametro'] == par]
            s = sen[(sen['parametro'] == par) & (np.abs(sen['t_h'] - 15.) < 1e-9)]
            for v in ref['valor']:
                k = np.where(np.abs(s['valor'] - v) < 1e-9)[0]
                if len(k) == 0:
                    continue
                rr = ref[ref['valor'] == v][0]
                add(f'Fig. 13, {par} = {v:g}: LFS medio limitado a 3 (15 h)', rr['LFS_medio'], s['LFS_medio_lim3'][k[0]],
                    'Fig. 13 (digit.)', rel=True)
                add(f'Fig. 13, {par} = {v:g}: A(LFS < 1) (15 h, m2)', rr['A_LFS_m2'], s['A_LFS_m2'][k[0]],
                    'Fig. 13 (digit.)', fmt='{:.1f}')
    # ---- tabela
    print('| grandeza | artigo | nosso | diferenca | fonte |')
    print('|---|---|---|---|---|')
    for l in linhas:
        print('| ' + ' | '.join(l) + ' |')
    if args.csv:
        with open(args.csv, 'w') as f:
            f.write('grandeza,artigo,nosso,diferenca,fonte\n')
            for l in linhas:
                f.write(','.join('"' + x + '"' for x in l) + '\n')

    # ---- figura dos perfis (tabela1)
    fig, ax = plt.subplots(3, 3, figsize=(13, 11), sharey=True)
    for i, t in enumerate((0., 7.5, 15.)):
        for j, (var, col, rot, lim) in enumerate((('theta', 'theta', 'theta', (0.2, 0.42)),
                                                  ('sigs_kPa', 'sigs_kPa', 'tensao de succao (kPa)', (0, 9)),
                                                  ('lfs', 'lfs', 'LFS (artigo satura em 3)', (0, 4)))):
            a = ax[i, j]
            for k, (xz, x) in enumerate(SECOES.items()):
                pa, va = perfil_artigo(dig, var, t, xz)
                pn, vn = perfil_nosso(per1, col, t, x)
                a.plot(va, pa, '.', ms=2.5, color=f'C{k}', alpha=0.6)
                a.plot(vn, pn, '-', lw=1.2, color=f'C{k}', label=f'x_zona = {xz:g} m')
            a.set_xlim(*lim)
            a.set_title(f'{rot}, t = {t:g} h', fontsize=9)
            a.grid(alpha=0.3)
            if j == 0:
                a.set_ylabel('profundidade (m)')
    ax[0, 0].invert_yaxis()
    ax[0, 0].legend(fontsize=7, title='linha: nosso; pontos: artigo', title_fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(args.figuras, 'perfis_v5.png'), dpi=110)

    # ---- Fig. 14
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    ts = [float(r['t_h']) for r in f14]
    w = 0.25
    for k, t in enumerate(ts):
        m = met4[np.abs(met4['t_h'] - t) < 1e-9]
        r = f14[k]
        ax[0].bar(k - w, r['LFS_medio_det'], w, color='C0', label='artigo' if k == 0 else None)
        if len(m):
            ax[0].bar(k, m['LFS_medio'][0], w, color='C1', label='nosso, eq. (11)' if k == 0 else None)
            ax[0].bar(k + w, campo(m[0], 'LFS_medio_lim3'), w, color='C2', label='nosso, LFS limitado a 3' if k == 0 else None)
            ax[1].bar(k + w / 2, m['A_LFS_m2'][0], w, color='C1', label='nosso' if k == 0 else None)
        ax[1].bar(k - w / 2, r['A_LFS_det_m2'], w, color='C0', label='artigo' if k == 0 else None)
    for a, rot in zip(ax, ('LFS medio', 'A(LFS < 1) (m2)')):
        a.set_xticks(range(len(ts)))
        a.set_xticklabels([f'{t:g} h' for t in ts])
        a.set_title(rot + ' - caso tabela4', fontsize=9)
        a.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(args.figuras, 'fig14_v5.png'), dpi=110)

    # ---- Fig. 13
    if sen is not None:
        fig, ax = plt.subplots(2, 3, figsize=(12, 6))
        for j, (par, rot) in enumerate((('nu', 'nu'), ('phi', "phi' (graus)"), ('c', "c' (kPa)"))):
            ref = f13[f13['parametro'] == par]
            for i, t in enumerate((0., 15.)):
                s = sen[(sen['parametro'] == par) & (np.abs(sen['t_h'] - t) < 1e-9)]
                if len(s) == 0:
                    continue
                ax[0, j].plot(s['valor'], s['LFS_medio_lim3'], 'o-', color=f'C{i + 1}', ms=3, label=f'nosso, limitado a 3, {t:g} h')
                ax[0, j].plot(s['valor'], s['LFS_medio'], 'x--', color=f'C{i + 1}', ms=3, label=f'nosso, eq. (11), {t:g} h')
                ax[1, j].plot(s['valor'], s['A_LFS_m2'], 'o-', color=f'C{i + 1}', ms=3, label=f'nosso, {t:g} h')
            ax[0, j].plot(ref['valor'], ref['LFS_medio'], 'k-', lw=2, alpha=0.5, label='artigo (Fig. 13)')
            ax[1, j].plot(ref['valor'], ref['A_LFS_m2'], 'k-', lw=2, alpha=0.5, label='artigo (Fig. 13)')
            ax[1, j].set_xlabel(rot)
        ax[0, 0].set_ylabel('LFS medio')
        ax[1, 0].set_ylabel('A(LFS < 1) (m2)')
        ax[0, 0].legend(fontsize=6)
        ax[1, 0].legend(fontsize=6)
        fig.tight_layout()
        fig.savefig(os.path.join(args.figuras, 'fig13_v5.png'), dpi=110)

    # ---- mapas
    subprocess.run([sys.executable, os.path.join(AQUI, 'mapas_v5.py'), args.tab1, '--saida',
                    os.path.join(args.figuras, 'mapas_tabela1.png')], check=False)


if __name__ == '__main__':
    main()
