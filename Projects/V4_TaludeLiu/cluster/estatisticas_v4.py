#!/usr/bin/env python3
# estatisticas_v4.py <pasta> [<pasta> ...]: estatistica das realizacoes concluidas de uma combinacao do Monte Carlo do
# V4 (as pastas r*/ com resumo.csv, de rodar_v4.sh): media, desvio e COV do FS (SRM), o beta e o P_f da lognormal
# ajustada (eqs. 12 e 13 de Liu et al. 2024: beta = ln[mu/sqrt(1+COV^2)]/sqrt(ln(1+COV^2)), P_f = Phi(-beta)), o P_f
# contado (FS < 1) com o intervalo de Wilson de 95 %, a area da massa que desliza (media e desvio) e os modos de ruptura
# (D_f: talude 0-5 m, pe 5-6 m, base 6-10 m, Fig. 9). Compara com dados/referencias.csv quando houver o valor do artigo.
# Grava <pasta>/mc_resumo.csv (todas as realizacoes juntas). So usa a biblioteca padrao do Python.
import csv, math, os, sys, glob

def wilson(k, n, z=1.96):
    if n == 0: return (float('nan'),) * 2
    p = k / n; d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d; h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(0., c - h), min(1., c + h)

def refs():
    arq = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'dados', 'referencias.csv')
    r = {}
    try:
        with open(arq) as f:
            for row in csv.DictReader(l for l in f if not l.startswith('#')):
                r[(row['grandeza'], row['cenario'], row['alfa'])] = row['valor']
    except OSError:
        pass
    return r

R = refs()
for pasta in sys.argv[1:]:
    linhas = []
    for a in sorted(glob.glob(os.path.join(pasta, 'r*', 'resumo.csv'))):
        if not os.path.exists(os.path.join(os.path.dirname(a), 'resultado.txt')): continue
        with open(a) as f:
            linhas += list(csv.DictReader(f))
    if not linhas:
        print(f"{pasta}: nenhuma realizacao concluida"); continue
    with open(os.path.join(pasta, 'mc_resumo.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(linhas[0].keys())); w.writeheader(); w.writerows(linhas)
    ok = [l for l in linhas if l['situacao'] == 'ok' and l['FS_min'] not in ('', '-1')]
    fs = [float(l['FS_min']) for l in ok]
    n = len(fs)
    mu = sum(fs) / n
    sd = math.sqrt(sum((x - mu) ** 2 for x in fs) / (n - 1)) if n > 1 else float('nan')
    cov = sd / mu
    beta = math.log(mu / math.sqrt(1 + cov * cov)) / math.sqrt(math.log(1 + cov * cov)) if n > 1 else float('nan')
    pf_ln = 0.5 * math.erfc(beta / math.sqrt(2))
    k = sum(1 for x in fs if x < 1.)
    lo, hi = wilson(k, n)
    areas = [float(l['massa_area_m2']) for l in ok if l['massa_area_m2']]
    am = sum(areas) / len(areas) if areas else float('nan')
    asd = math.sqrt(sum((x - am) ** 2 for x in areas) / (len(areas) - 1)) if len(areas) > 1 else float('nan')
    modos = {}
    for l in ok: modos[l['modo']] = modos.get(l['modo'], 0) + 1
    caso, alfa = ok[0]['caso'], ok[0]['angulo']
    cen = caso[1:] if caso.startswith('s') else caso
    print(f"{pasta}: caso {caso}, alfa {alfa}, rho {ok[0]['rho']}, {n} realizacoes ok de {len(linhas)}")
    print(f"  FS: media {mu:.4f}, desvio {sd:.4f}, COV {cov:.4f}, min {min(fs):.4f}, max {max(fs):.4f}")
    print(f"  lognormal (eqs. 12-13): beta = {beta:.3f}, P_f = {pf_ln:.3e}")
    print(f"  contado: P_f = {k}/{n} = {k / n:.3e}, Wilson 95 %: [{lo:.2e}, {hi:.2e}]"
          f"{' (P_f do artigo ~1e-4 ou menor: nao resolvivel com este N)' if n < 10000 else ''}")
    print(f"  area da massa que desliza: media {am:.1f} m2, desvio {asd:.1f} m2; modos: {modos}")
    for g, v in (('media_FS', mu), ('dp_FS', sd), ('beta', beta), ('Pf', pf_ln), ('media_area_m2', am), ('dp_area_m2', asd)):
        ref = R.get((g, cen, alfa))
        if ref is not None:
            print(f"    {g:14s} nosso {v:10.4g}   artigo {float(ref):10.4g}")
