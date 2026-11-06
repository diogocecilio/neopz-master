#!/usr/bin/env python3
# V5_EncostaDuasCamadas: ajuste da curva de retencao de van Genuchten (a unica proxima de Brooks-Corey na biblioteca) as
# curvas de Brooks-Corey de Abbasov et al. (2024), Tabela 1, das duas camadas.
#
# Brooks-Corey (eq. 3 do artigo): S_e = (alpha_BC h)^(-n_BC) se alpha_BC h > 1, senao 1 (h = succao em coluna d'agua).
# van Genuchten da biblioteca (TPZUnsaturatedHydraulics::EVanGenuchten): S_e = [1 + (a_VG p_c)^b]^(-m), m = 1 - 1/b,
# com p_c em kPa (p_c = gamma_w h, gamma_w = 9.81 kN/m3).
# Ajuste: minimos quadrados de S_e em 0 < h <= H_MAX (padrao 1,5 m: a faixa de succao do problema depois dos 10 anos,
# no maximo ~1,1 m na camada superior e ~0,9 m na inferior com Brooks-Corey), pontos igualmente espacados, a_VG e b
# livres; busca em grade seguida de refino (sem scipy).
# A condutividade nao precisa de ajuste: K = K_s S_e^(l + 2 + 2/n_BC) (eq. 4) e a lei de potencia k_rw = S_e^p da
# biblioteca com p = 9,211 (superior) e 18,267 (inferior), com o S_e = (S_w - S_wr)/(1 - S_wr) de van Genuchten.
#
# Uso: python3 ajuste_vg.py [--hmax 1.5] [--figura ../figuras/ajuste_vg.png]
# Imprime a_VG (1/kPa) e b de cada camada (os valores do problema.h) e os erros nos estados do problema.
import argparse
import numpy as np

GW = 9.81
CAMADAS = {  # Tabela 1 (theta_s, theta_r, K_s em cm/h, alpha_BC em 1/cm, n_BC); l_BC = 1
    'superior': dict(ts=0.412, tr=0.041, ks=2.59, alpha=0.068, n=0.322),
    'inferior': dict(ts=0.385, tr=0.09, ks=0.06, alpha=0.027, n=0.131),
}
Q600 = 0.6 / (365.25 * 86400.)  # 600 mm/ano em m/s


def se_bc(h, L):
    x = L['alpha'] * 100. * h
    return np.where(x > 1., np.maximum(x, 1e-300) ** (-L['n']), 1.)


def se_vg(h, a, b):  # a em 1/m
    return (1. + (a * h) ** b) ** (-(1. - 1. / b))


def h_vg(se, a, b):
    m = 1. - 1. / b
    return (se ** (-1. / m) - 1.) ** (1. / b) / a


def h_bc(se, L):
    return se ** (-1. / L['n']) / (L['alpha'] * 100.)


def ajusta(L, hmax, n=300):
    h = np.linspace(0., hmax, n + 1)[1:]
    y = se_bc(h, L)

    def erro(a, b):
        e = se_vg(h, a, b) - y
        return np.mean(e * e)
    melhor = (np.inf, 1., 1.5)
    for la in np.linspace(np.log(0.05), np.log(50.), 250):
        for b in np.linspace(1.02, 4., 250):
            s = erro(np.exp(la), b)
            if s < melhor[0]:
                melhor = (s, np.exp(la), b)
    s, a, b = melhor
    da, db = 0.05, 0.02
    for _ in range(6):
        for la in np.linspace(np.log(a) - da, np.log(a) + da, 41):
            for bb in np.linspace(max(1.001, b - db), b + db, 41):
                s2 = erro(np.exp(la), bb)
                if s2 < melhor[0]:
                    melhor = (s2, np.exp(la), bb)
        s, a, b = melhor
        da, db = da / 4., db / 4.
    return a, b, np.sqrt(s)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--hmax', type=float, default=1.5)
    ap.add_argument('--figura', default='')
    args = ap.parse_args()
    res = {}
    for nome, L in CAMADAS.items():
        p = 3. + 2. / L['n']
        ks = L['ks'] / 100. / 3600.
        a, b, rmse = ajusta(L, args.hmax)
        res[nome] = (a, b)
        se_perm = (Q600 / ks) ** (1. / p)       # gradiente unitario com 600 mm/ano
        se_chuva = min(1., (20. / (L['ks'] * 10.)) ** (1. / p)) if L['ks'] * 10. > 20. else 1.
        print(f"{nome}: a_VG = {a / GW:.5f} 1/kPa ({a:.4f} 1/m), b_VG = {b:.4f}, m = {1 - 1 / b:.4f}, "
              f"S_wr = {L['tr'] / L['ts']:.4f}, p = {p:.3f}, RMSE(S_e) = {rmse:.4f} em 0 < h <= {args.hmax} m")
        for rot, se in (('600 mm/ano', se_perm), ('20 mm/h', se_chuva)):
            if se >= 1.:
                print(f"   {rot}: S_e = 1 (K_s < chuva)")
                continue
            hb, hv = h_bc(se, L), h_vg(se, a, b)
            print(f"   {rot}: S_e = {se:.4f}, theta = {L['tr'] + (L['ts'] - L['tr']) * se:.4f}, h_BC = {hb:.3f} m, "
                  f"h_VG = {hv:.3f} m; tensao de succao S_e gamma_w h: BC {se * GW * hb:.2f} kPa, VG {se * GW * hv:.2f} kPa")
    if args.figura:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
        h = np.logspace(-3, 1.2, 400)
        for k, (nome, L) in enumerate(CAMADAS.items()):
            a, b = res[nome]
            p = 3. + 2. / L['n']
            th = lambda se: L['tr'] + (L['ts'] - L['tr']) * se
            ax[0].semilogx(h, th(se_bc(h, L)), f'C{k}-', label=f'{nome}: Brooks-Corey (artigo)')
            ax[0].semilogx(h, th(se_vg(h, a, b)), f'C{k}--', label=f'{nome}: van Genuchten ajustado')
            ax[1].loglog(h, L['ks'] * 10. * se_bc(h, L) ** p, f'C{k}-')
            ax[1].loglog(h, L['ks'] * 10. * se_vg(h, a, b) ** p, f'C{k}--')
        for x in ax:
            x.axvspan(1e-3, args.hmax, color='0.92', zorder=0)
            x.set_xlabel('succao h (m de agua)')
        ax[0].set_ylabel('theta')
        ax[1].set_ylabel('K (mm/h)')
        ax[1].axhline(600. / 365.25 / 24., color='k', lw=0.6, ls=':')
        ax[1].axhline(20., color='k', lw=0.6, ls=':')
        ax[0].legend(fontsize=7)
        ax[0].set_title('retencao (faixa do ajuste em cinza)')
        ax[1].set_title('K = K_s S_e^p (linhas: 600 mm/ano e 20 mm/h)')
        fig.tight_layout()
        fig.savefig(args.figura, dpi=130)


if __name__ == '__main__':
    main()
