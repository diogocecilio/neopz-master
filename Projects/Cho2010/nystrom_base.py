#!/usr/bin/env python3
"""nystrom_base.py: a base de Karhunen-Loeve de Nystrom de um talude para o RandomField (opcao nystrom= do Cho2010 e do
SlopeRainMonteCarlo).

O dominio e o poligono do talude, os 6 pontos do kTalude dos drivers (base esquerda, base direita, terreno do pe, pe,
crista, crista a esquerda); padrao: o talude do Cho (2010). A quadratura e a regra do ponto medio numa grade de celulas
dx x dy sobre a caixa do talude: os nos sao os centros das celulas que caem dentro do talude, e o peso de cada no e a
area da celula. Os autopares saem de W^(1/2) C W^(1/2) v = lambda v, com f_k = v_k / sqrt(w) nos nos
(sum_j w_j f_k^2 = 1). O RandomField interpola f_k em qualquer ponto pela formula de Nystrom:
f_k(x) = (1/lambda_k) sum_j w_j C(x, x_j) f_k(x_j), com o mesmo nucleo, que vai no cabecalho do arquivo.

Nucleos (sem o fator 2 e sem o pi):
  exp    C = exp(-|dx|/lh - |dy|/lv)                                    (Cho 2010, Hu et al. 2024; o padrao)
  gauss  C = exp(-(dx cos a + dy sin a)^2/lh^2 - (-dx sin a + dy cos a)^2/lv^2), a = angulo em graus
         (Liu et al. 2024, Eq. 1; o exp{-pi[(dx/30)^2 + (dy/3)^2]} de Li et al. 2023 e lh = 30/sqrt(pi), lv = 3/sqrt(pi))

Uso: python3 -I nystrom_base.py <saida.txt | -> [M=200] [talude=x0,y0,...,x5,y5] [nucleo=exp] [lh=20] [lv=2]
                                [angulo=0] [dx=0.5] [dy=0.25] [copias=] [nos=] [autopares=] [compara=] [mercer=0]
  <saida.txt>      a base; "-": nao grava (so o resumo, a comparacao e o teste de Mercer)
  talude=          os 6 pontos (12 numeros, x e y alternados), ex. o caso A (Hu et al. 2024):
                   talude=0,0,70,0,70,10,40,10,30,20,0,20
  dx, dy           a celula; se a caixa do talude nao tem um numero inteiro de celulas, dx e dy sao ajustados (avisa)
  copias=20,100    tambem grava <saida>_M20.txt e <saida>_M100.txt com os 20 e os 100 primeiros modos (os mesmos
                   autovetores, com os mesmos sinais: os M primeiros termos de uma valem para a outra). O FORM com poucos
                   termos le o arquivo pequeno; a amostragem com todos os termos, o grande
  nos=<arquivo>    so grava os nos (x y por linha) e sai: para calcular os autopares fora (base grande, ver abaixo)
  autopares=<pref> le os autopares de <pref>.lam (M autovalores) e <pref>.f (n x M valores f_k(x_j), linha por no), em
                   Real64, no lugar do calculo daqui (os sinais ficam os do arquivo)
  compara=<arq>    autovalores de referencia (<pref>.lam do Mathematica, ou outra base .txt): o erro relativo de
                   lambda_1..lambda_M desta grade contra eles (a escolha da grade)
  mercer=<pares>   teste de Mercer em <pares> pares de pontos sorteados no talude (metade perto, com |dx| < lh e
                   |dy| < lv, metade em qualquer lugar): C_M(x, y) = sum_{k <= M} lambda_k f_k(x) f_k(y), com f_k
                   interpolado por Nystrom, contra o nucleo C(x, y), para M = 100, 150, 200 e o M da base
Os autopares daqui: nucleo separavel (exp, ou gauss com angulo 0) -> eigsh (ARPACK, scipy.sparse.linalg) com o
produto matriz-vetor sem a matriz, C = C_y (x) C_x na grade cheia (memoria O(n M), qualquer n); gauss girado -> a matriz
densa (n^2 x 8 bytes) e eigsh. Vetor inicial fixo e o sinal de cada autovetor escolhido por uma projecao fixa: a mesma
grade da sempre a mesma base, com os mesmos sinais.

Base grande (milhares de nos e de modos, M > ~300): o eigsh fica caro, e o caminho e o Mathematica:
  python3 -I nystrom_base.py base.txt M=2000 dx=0.2 dy=0.1 nos=base.nos
  wolframscript -file nystrom_autopares.wl base.nos 2000 20 2 0.02 base      # o Eigensystem do Mathematica (MKL)
  python3 -I nystrom_base.py base.txt M=2000 dx=0.2 dy=0.1 autopares=base copias=100
A grade de 0,2 x 0,1 m (15 000 nos) retem 99,7% da variancia com 2000 modos; a de 0,5 x 0,25 m (2400 nos) erra os
autovalores acima de ~50 modos (lambda_100 23% e lambda_1000 57% contra a de 0,2 x 0,1).
Arquivo: 'n M nucleo=... lh=... lv=... [angulo=...] dx=... dy=...' (o cabecalho antigo, so 'n M', e o nucleo exp com
as distancias do driver) e, depois, n linhas 'x y w', 1 linha com os M autovalores e n linhas com os M valores
f_k(x_j).

Coluna 1D (Cho 2014, o SlopeRainColumn): python3 -I nystrom_base.py <saida.txt | -> coluna=H [dy=H/(8M)] [nucleo=exp]
                                [lv=2] [M=200] [copias=] [compara=] [mercer=0]
  a base do segmento [0, H] em x = 0: nos (0, (j + 1/2) dy) com peso dy (n = H/dy nos; sem dy=, n = 8M, ex. dy = 0.0025
  com H = 2 e M = 100) e os M maiores autopares da matriz densa w C (eigh, so o fator 1D em y do nucleo); o sinal de
  cada modo e v_k . cos(2.345 y + 0.5) > 0 (a projecao do talude em x = 0). O cabecalho leva lh=1e+30 e dx=1.0: o nucleo
  do RandomField fica exatamente 1D e a grade rapida roda com nx = 1. Imprime lambda_1..5, a variancia retida e, com o
  nucleo exp, o erro contra a KL analitica do KLSeparavel.h (omega tan(omega a) = c e omega + c tan(omega a) = 0 em
  [-a, a], lambda = 2c/(omega^2 + c^2), c = 1/lv); mercer= sorteia os pares em [0, H]. Recusa lambda_M <= 1e-12 lambda_1
  (M menor) e as opcoes do talude (talude=, lh=, dx=, angulo=, nos=, autopares=)."""
import math
import sys
import time

import numpy as np

CHO = [(0, 0), (30, 0), (30, 5), (20, 5), (10, 15), (0, 15)]  # o cho2010.cpp (kTalude)


def dentro(talude, x, y):
    """ponto dentro do poligono do talude (raio)"""
    c = False
    n = len(talude)
    for i in range(n):
        x1, y1 = talude[i]
        x2, y2 = talude[(i + 1) % n]
        if (y1 > y) != (y2 > y) and x < x1 + (y - y1) * (x2 - x1) / (y2 - y1):
            c = not c
    return c


def area_poligono(t):
    return 0.5 * abs(sum(t[i][0] * t[(i + 1) % len(t)][1] - t[(i + 1) % len(t)][0] * t[i][1] for i in range(len(t))))


class Nucleo:
    """o nucleo de correlacao: C(dx, dy) e, se separavel, os fatores 1D"""

    def __init__(self, tipo, lh, lv, angulo):
        if tipo not in ("exp", "gauss"):
            sys.exit(f"nucleo={tipo}: exp ou gauss")
        self.tipo, self.lh, self.lv, self.angulo = tipo, lh, lv, angulo
        self.separavel = tipo == "exp" or angulo == 0.

    def __call__(self, dx, dy):
        if self.tipo == "exp":
            return np.exp(-np.abs(dx) / self.lh - np.abs(dy) / self.lv)
        a = math.radians(self.angulo)
        u, v = dx * math.cos(a) + dy * math.sin(a), -dx * math.sin(a) + dy * math.cos(a)
        return np.exp(-(u / self.lh) ** 2 - (v / self.lv) ** 2)

    def fator(self, d, L):
        """o fator 1D de um nucleo separavel"""
        return np.exp(-np.abs(d) / L) if self.tipo == "exp" else np.exp(-(d / L) ** 2)

    def cabecalho(self):
        s = f"nucleo={self.tipo} lh={self.lh!r} lv={self.lv!r}"
        return s + (f" angulo={self.angulo!r}" if self.tipo == "gauss" else "")


def grade(talude, dx, dy):
    """os centros das celulas dentro do talude (y por fora, x por dentro: a ordem do arquivo), a mascara na grade cheia e
    os dx, dy efetivos"""
    xmin, xmax = min(p[0] for p in talude), max(p[0] for p in talude)
    ymin, ymax = min(p[1] for p in talude), max(p[1] for p in talude)
    nx, ny = max(1, int(round((xmax - xmin) / dx))), max(1, int(round((ymax - ymin) / dy)))
    dxe, dye = (xmax - xmin) / nx, (ymax - ymin) / ny
    if abs(dxe - dx) > 1e-12 * dx or abs(dye - dy) > 1e-12 * dy:
        print(f"  aviso: a caixa nao tem um numero inteiro de celulas de {dx} x {dy}: uso {dxe:.6g} x {dye:.6g}")
        dx, dy = dxe, dye
    # o mesmo arredondamento do np.arange(dx/2, L, dx) da versao antiga (x0 = 0): os nos da base do Cho nao mudam
    xs = (xmin + dx / 2) + dx * np.arange(nx)
    ys = (ymin + dy / 2) + dy * np.arange(ny)
    mask = np.array([[dentro(talude, x, y) for x in xs] for y in ys])
    pts = np.array([(x, y) for j, y in enumerate(ys) for i, x in enumerate(xs) if mask[j, i]])
    return pts, xs, ys, mask, dx, dy


def autopares(pts, xs, ys, mask, w, nuc, M):
    """os M maiores autopares de w C (a grade e uniforme): eigsh com o produto pela matriz sem forma-la (separavel) ou
    com a matriz densa (gauss girado)"""
    from scipy.sparse.linalg import LinearOperator, eigsh
    n = len(pts)
    if nuc.separavel:
        Cx = nuc.fator(xs[:, None] - xs[None, :], nuc.lh)
        Cy = nuc.fator(ys[:, None] - ys[None, :], nuc.lv)
        V = np.zeros(mask.shape)

        def mv(v):
            V[...] = 0.
            V[mask] = np.ravel(v)
            return w * (Cy @ V @ Cx.T)[mask]

        A = LinearOperator((n, n), matvec=mv, dtype=float)
    else:
        mem = n * n * 8 / 2 ** 30
        print(f"  nucleo girado (nao separavel): matriz densa de {n} x {n} ({mem:.1f} GB)")
        A = w * nuc(pts[:, 0, None] - pts[None, :, 0], pts[:, 1, None] - pts[None, :, 1])
    v0 = np.cos(0.37 * pts[:, 0] + 0.73 * pts[:, 1]) + 1.5  # vetor inicial fixo: a mesma base em toda rodada
    t0 = time.time()
    lam, v = eigsh(A, k=M, which="LA", v0=v0, ncv=min(n, max(2 * M + 1, 20)), tol=0.)
    ordem = np.argsort(lam)[::-1]
    lam, v = lam[ordem], v[:, ordem]
    g = np.cos(1.234 * pts[:, 0] + 2.345 * pts[:, 1] + 0.5)  # sinal: v_k . g > 0 (projecao fixa, generica)
    v *= np.where(v.T @ g < 0., -1., 1.)[None, :]
    res = max(np.linalg.norm((A @ v[:, k]) - lam[k] * v[:, k]) for k in (0, M // 2, M - 1))
    print(f"  eigsh: {M} autopares em {time.time() - t0:.0f} s (residuo |A v - lambda v| <= {res:.1e})")
    return lam, v


def grava(arq, pts, w, lam, f, cab):
    """a base no formato que o RandomField::UseNystrom le"""
    n, M = f.shape
    with open(arq, "w") as a:
        a.write(f"{n} {M} {cab}\n")
        np.savetxt(a, np.column_stack([pts, w]), fmt=["%.10f", "%.10f", "%.12e"])
        a.write(" ".join(f"{x:.12e}" for x in lam) + "\n")
        np.savetxt(a, f, fmt="%.10e")


def le_lambdas(arq):
    """os autovalores de uma base .txt (a linha depois dos nos) ou de um .lam (Real64)"""
    if not arq.endswith(".txt"):
        return np.fromfile(arq, dtype="<f8")
    with open(arq) as a:
        n = int(a.readline().split()[0])
        for _ in range(n):
            a.readline()
        return np.array([float(x) for x in a.readline().split()])


def interpola(nuc, pts, w, lam, f, x):
    """f_k(x) nos pontos x (m x 2) pela formula de Nystrom: (1/lambda_k) sum_j w C(x, x_j) f_k(x_j)"""
    out = np.empty((len(x), f.shape[1]))
    for i0 in range(0, len(x), 200):
        c = nuc(x[i0:i0 + 200, 0, None] - pts[None, :, 0], x[i0:i0 + 200, 1, None] - pts[None, :, 1])
        out[i0:i0 + 200] = (c * w) @ f / lam[None, :]
    return out


def mercer(talude, nuc, pts, w, lam, f, npares, Ms):
    """C_M(x, y) contra C(x, y) em pares sorteados (semente fixa)"""
    rng = np.random.default_rng(2024)
    xmin, xmax = min(p[0] for p in talude), max(p[0] for p in talude)
    ymin, ymax = min(p[1] for p in talude), max(p[1] for p in talude)

    def sorteia(m, perto=None):
        out = []
        while len(out) < m:
            if perto is None:
                x, y = rng.uniform(xmin, xmax), rng.uniform(ymin, ymax)
            else:
                x0, y0 = perto[len(out)]
                x, y = x0 + rng.uniform(-nuc.lh, nuc.lh), y0 + rng.uniform(-nuc.lv, nuc.lv)
            if dentro(talude, x, y):
                out.append((x, y))
        return np.array(out)

    a = sorteia(npares)
    b = np.vstack([sorteia(npares // 2, perto=a[:npares // 2]), sorteia(npares - npares // 2)])
    fa, fb = interpola(nuc, pts, w, lam, f, a), interpola(nuc, pts, w, lam, f, b)
    C = nuc(a[:, 0] - b[:, 0], a[:, 1] - b[:, 1])
    print(f"  teste de Mercer ({npares} pares; metade perto, |dx| < lh e |dy| < lv; C medio {C.mean():.3f}):")
    print("      M   variancia C_M(x,x) (media; min-max)   erro |C_M - C|: medio    rms     max   | (C_M - C) medio")
    for M in Ms:
        if M > len(lam):
            continue
        var = (fa[:, :M] ** 2 * lam[None, :M]).sum(axis=1)
        CM = (fa[:, :M] * fb[:, :M] * lam[None, :M]).sum(axis=1)
        e = CM - C
        print(f"   {M:4d}   {var.mean():.4f} ({var.min():.4f}-{var.max():.4f})              "
              f"{np.abs(e).mean():.4f}  {np.sqrt((e ** 2).mean()):.4f}  {np.abs(e).max():.4f}  |  {e.mean():+.4f}")


def kl_exp_1d(H, L, n):
    """os n maiores autovalores e as autofuncoes (em y de [0, H]) do nucleo exp(-|s - t|/L) em forma fechada, como o
    CalculaModos1D do KLSeparavel.h (Ghanem e Spanos 1991, sec. 2.3.3): com a = H/2 e c = 1/L, os pares de
    omega tan(omega a) = c (cos) e os impares de omega + c tan(omega a) = 0 (sin), lambda = 2c/(omega^2 + c^2)"""
    from scipy.optimize import brentq
    a, c, pi, e = 0.5 * H, 1. / L, math.pi, 1e-12
    v = []
    for k in range(1, n + 1):
        w = brentq(lambda x: x * math.tan(x * a) - c, ((k - 1) * pi + e) / a, ((k - 0.5) * pi - e) / a, xtol=1e-15,
                   maxiter=500)
        v.append((2. * c / (w * w + c * c), w, 1))
        w = brentq(lambda x: x + c * math.tan(x * a), ((k - 0.5) * pi + e) / a, (k * pi - e) / a, xtol=1e-15, maxiter=500)
        v.append((2. * c / (w * w + c * c), w, 0))
    v = sorted(v, key=lambda t: -t[0])[:n]
    lam = np.array([t[0] for t in v])

    def f(k, y):
        _, w, par = v[k]
        s2 = math.sin(2. * w * a) / (2. * w)
        return (np.cos(w * (y - a)) / math.sqrt(a + s2)) if par else (np.sin(w * (y - a)) / math.sqrt(a - s2))
    return lam, f


def mercer_coluna(H, nuc, pts, w, lam, f, npares, Ms):
    """o teste de Mercer da coluna: C_M(y, y') contra C(y, y') em pares sorteados em [0, H] (semente fixa; metade perto,
    com |dy| < lv)"""
    rng = np.random.default_rng(2024)
    a = rng.uniform(0., H, npares)
    b = np.empty(npares)
    for i in range(npares):
        while True:
            b[i] = a[i] + rng.uniform(-nuc.lv, nuc.lv) if i < npares // 2 else rng.uniform(0., H)
            if 0. <= b[i] <= H:
                break
    pa, pb = np.column_stack([np.zeros(npares), a]), np.column_stack([np.zeros(npares), b])
    fa, fb = interpola(nuc, pts, w, lam, f, pa), interpola(nuc, pts, w, lam, f, pb)
    C = nuc(0. * a, a - b)
    print(f"  teste de Mercer ({npares} pares em [0, {H:g}]; metade perto, |dy| < lv; C medio {C.mean():.3f}):")
    print("      M   variancia C_M(y,y) (media; min-max)   erro |C_M - C|: medio    rms     max   | (C_M - C) medio")
    for M in Ms:
        if M > len(lam):
            continue
        var = (fa[:, :M] ** 2 * lam[None, :M]).sum(axis=1)
        CM = (fa[:, :M] * fb[:, :M] * lam[None, :M]).sum(axis=1)
        e = CM - C
        print(f"   {M:4d}   {var.mean():.4f} ({var.min():.4f}-{var.max():.4f})              "
              f"{np.abs(e).mean():.4f}  {np.sqrt((e ** 2).mean()):.4f}  {np.abs(e).max():.4f}  |  {e.mean():+.4f}")


def coluna(p, dados):
    """a base 1D da coluna vertical [0, H] (opcao coluna=H): nos (0, (j + 1/2) dy) com peso dy, os M maiores autopares da
    matriz densa w C e o cabecalho com lh = 1e30 e dx = 1 (o nucleo do RandomField fica 1D; a grade rapida, nx = 1)"""
    if not p["coluna"] > 0.:
        sys.exit("coluna=H: H > 0 (a altura da coluna, m)")
    proibidas = sorted(dados & {"talude", "lh", "dx", "nos", "autopares"})
    if proibidas or p["angulo"] != 0.:
        sys.exit(f"coluna=: sem {', '.join(proibidas + (['angulo'] if p['angulo'] != 0. else []))} (a coluna fica em"
                 f" x = 0; o cabecalho leva lh=1e+30 e dx=1.0)")
    H, M = p["coluna"], p["M"]
    n = max(1, int(round(H / p["dy"]))) if "dy" in dados else 8 * M
    dy = H / n
    if "dy" in dados and abs(dy - p["dy"]) > 1e-12 * p["dy"]:
        print(f"  aviso: H = {H:g} nao tem um numero inteiro de celulas de {p['dy']}: uso dy = {dy:.6g}")
    M = min(M, n)
    if n < 8 * M:
        print(f"  aviso: n = {n} nos < 8M = {8 * M}: os ultimos modos ficam menos precisos (o padrao e n = 8M)")
    nuc = Nucleo(p["nucleo"], 1e30, p["lv"], 0.)
    ys = (0. + dy / 2) + dy * np.arange(n)  # a mesma formula dos centros da grade do talude
    pts = np.column_stack([np.zeros(n), ys])
    w = dy
    A = w * nuc.fator(ys[:, None] - ys[None, :], nuc.lv)
    t0 = time.time()
    lam, v = np.linalg.eigh(A)
    ordem = np.argsort(lam)[::-1][:M]
    lam, v = lam[ordem], v[:, ordem]
    if not lam[-1] > 1e-12 * lam[0]:
        k = int(np.argmax(~(lam > 1e-12 * lam[0])))
        sys.exit(f"coluna=: lambda_{k + 1} = {lam[k]:.3e} <= 1e-12 lambda_1 ({lam[0]:.3e}): use M <= {k}")
    g = np.cos(2.345 * ys + 0.5)  # sinal: v_k . g > 0 (a projecao do talude, cos(1.234 x + 2.345 y + 0.5), em x = 0)
    proj = v.T @ g
    v *= np.where(proj < 0., -1., 1.)[None, :]
    res = max(np.linalg.norm((A @ v[:, k]) - lam[k] * v[:, k]) for k in (0, M // 2, M - 1))
    f = v / math.sqrt(w)
    cab = f"{nuc.cabecalho()} dx={1.0!r} dy={dy!r}"
    print(f"{n} nos de dy = {dy:.6g} m na coluna de H = {H:g} m (x = 0), {cab}, M = {M}")
    print(f"  eigh: {M} autopares em {time.time() - t0:.1f} s (residuo |A v - lambda v| <= {res:.1e}; "
          f"sinal: min |v_k . g| = {np.abs(proj).min():.2e})")
    print(f"  normalizacao: max |sum_j w_j f_k^2 - 1| = {np.abs((w * f * f).sum(axis=0) - 1).max():.2e}")
    print("  lambda_1..5: " + " ".join(f"{x:.6f}" for x in lam[:5]))
    for m in (10, 20, 50, 100, 150, 200, 500, 1000, 2000):
        if m <= M:
            print(f"  variancia retida com M = {m:4d}: {lam[:m].sum() / H:.4f}")
    if M not in (10, 20, 50, 100, 150, 200, 500, 1000, 2000):
        print(f"  variancia retida com M = {M:4d}: {lam.sum() / H:.4f}")
    if nuc.tipo == "exp":
        an, fan = kl_exp_1d(H, nuc.lv, M)
        e = lam / an - 1.
        K = min(20, M)
        # as autofuncoes analiticas com o sinal da regra acima
        ef = []
        for k in range(min(5, M)):
            fa = fan(k, ys)
            fa *= -1. if fa @ g < 0. else 1.
            ef.append(np.abs(f[:, k] - fa).max())
        print(f"  KL analitica (KLSeparavel.h): erro relativo max em lambda_1..{K} = {np.abs(e[:K]).max():.2e}" +
              (f", em lambda_1..{M} = {np.abs(e).max():.2e} (no modo {int(np.abs(e).argmax()) + 1})" if M > K else "") +
              f"; variancia retida analitica com M = {M}: {an.sum() / H:.4f}; max |f_k - f_k analitica| "
              f"(k = 1..{len(ef)}) = {max(ef):.2e}")
    if p["compara"]:
        ref = le_lambdas(p["compara"])
        K = min(M, len(ref))
        e = lam[:K] / ref[:K] - 1.
        print(f"  contra {p['compara']} (lambda_1 = {ref[0]:.6g}): max |erro| em lambda_1..{K} = {100 * np.abs(e).max():.2f}%"
              f" (no modo {int(np.abs(e).argmax()) + 1}); " +
              ", ".join(f"lambda_{k} {100 * e[k - 1]:+.2f}%" for k in (1, 10, 50, 100, 150, 200) if k <= K))
    if p["mercer"] > 0:
        mercer_coluna(H, nuc, pts, np.full(n, w), lam, f, p["mercer"], sorted({20, 50, 100, M}))
    if sys.argv[1] == "-":
        return
    grava(sys.argv[1], pts, np.full(n, w), lam, f, cab)
    print(f"  gravada: {sys.argv[1]}")
    base = sys.argv[1][:-4] if sys.argv[1].endswith(".txt") else sys.argv[1]
    for m in (int(c) for c in p["copias"].split(",") if c):
        if m < M:
            grava(f"{base}_M{m}.txt", pts, np.full(n, w), lam[:m], f[:, :m], cab)
            print(f"  copia com os {m} primeiros modos: {base}_M{m}.txt")


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    p = {"M": 200, "talude": "", "nucleo": "exp", "lh": 20., "lv": 2., "angulo": 0., "dx": 0.5, "dy": 0.25,
         "copias": "", "nos": "", "autopares": "", "compara": "", "mercer": 0, "coluna": 0.}
    dados = set()  # as opcoes dadas (a coluna tem outro padrao de dy)
    for a in sys.argv[2:]:
        k, v = a.split("=", 1)
        if k not in p:
            sys.exit(f"opcao desconhecida: {a}")
        p[k] = type(p[k])(v)
        dados.add(k)
    if "coluna" in dados:
        coluna(p, dados)
        return
    talude = CHO
    if p["talude"]:
        v = [float(x) for x in p["talude"].split(",")]
        if len(v) != 12:
            sys.exit("talude: 12 numeros (os 6 pontos x,y do kTalude)")
        talude = list(zip(v[0::2], v[1::2]))
    nuc = Nucleo(p["nucleo"], p["lh"], p["lv"], p["angulo"])
    pts, xs, ys, mask, dx, dy = grade(talude, p["dx"], p["dy"])
    n = len(pts)
    w = dx * dy
    if p["nos"]:
        np.savetxt(p["nos"], pts, fmt="%.10f")
        print(f"{n} nos em {p['nos']}")
        return
    M = min(p["M"], n)
    if p["autopares"]:
        lam = np.fromfile(p["autopares"] + ".lam", dtype="<f8")
        f = np.fromfile(p["autopares"] + ".f", dtype="<f8")
        if len(lam) < M or f.size != n * len(lam):
            sys.exit(f"autopares de {p['autopares']}: {len(lam)} autovalores e {f.size} valores, esperados {M} e n x M"
                     f" com n = {n}")
        f = f.reshape(n, len(lam))[:, :M]
        lam = lam[:M]
        if np.any(np.diff(lam) > 0):
            sys.exit("os autovalores nao estao em ordem decrescente")
    else:
        lam, v = autopares(pts, xs, ys, mask, w, nuc, M)
        f = v / math.sqrt(w)
    area = w * n
    print(f"{n} nos de {dx:.4g} x {dy:.4g} m (area {area:.1f} m2, a do talude e {area_poligono(talude):.1f}), "
          f"{nuc.cabecalho()}, M = {M}; matriz densa equivalente {n * n * 8 / 2 ** 30:.1f} GB, base {n * M * 8 / 2 ** 20:.0f} MB"
          f" na memoria")
    print(f"  normalizacao: max |sum_j w_j f_k^2 - 1| = {np.abs((w * f * f).sum(axis=0) - 1).max():.2e}")
    for m in (10, 20, 50, 100, 150, 200, 500, 1000, 2000):
        if m <= M:
            print(f"  variancia retida com M = {m:4d}: {lam[:m].sum() / area:.4f}")
    if p["compara"]:
        ref = le_lambdas(p["compara"])
        K = min(M, len(ref))
        e = lam[:K] / ref[:K] - 1.
        print(f"  contra {p['compara']} (lambda_1 = {ref[0]:.6g}): max |erro| em lambda_1..{K} = {100 * np.abs(e).max():.2f}%"
              f" (no modo {int(np.abs(e).argmax()) + 1}); " +
              ", ".join(f"lambda_{k} {100 * e[k - 1]:+.2f}%" for k in (1, 10, 50, 100, 150, 200) if k <= K))
    if p["mercer"] > 0:
        mercer(talude, nuc, pts, np.full(n, w), lam, f, p["mercer"], sorted({100, 150, 200, M}))
    if sys.argv[1] == "-":
        return
    cab = f"{nuc.cabecalho()} dx={dx!r} dy={dy!r}"
    grava(sys.argv[1], pts, np.full(n, w), lam, f, cab)
    print(f"  gravada: {sys.argv[1]}")
    base = sys.argv[1][:-4] if sys.argv[1].endswith(".txt") else sys.argv[1]
    for m in (int(c) for c in p["copias"].split(",") if c):
        if m < M:
            grava(f"{base}_M{m}.txt", pts, np.full(n, w), lam[:m], f[:, :m], cab)
            print(f"  copia com os {m} primeiros modos: {base}_M{m}.txt")


if __name__ == "__main__":
    main()
