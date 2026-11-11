"""compara_gstools.py: o campo da KL em forma fechada do retangulo (a do SlopeRandomField, M = 2000) contra o metodo
de randomizacao do GSTools (Hesse et al. 2014: u(x) = sqrt(1/N) sum Z1 cos(k.x) + Z2 sin(k.x), k sorteado da
densidade espectral), nos pontos de uma grade de 0,25 m no talude de Vargas Ceron et al. (2025).

Tres geradores, todos com media 0 e variancia 1:
  KL      a forma fechada do retangulo [0,25] x [0,10] com M = 2000 termos (klfechada.KL2D), restrita ao talude;
  GS-exp  o gs.Exponential do GSTools com len_scale = [20, 2]: exp(-r) com r = sqrt((dx/20)^2 + (dy/2)^2), a distancia
          euclidiana escalada (anisotropia elipsoidal) -- outro nucleo que o do Ceron;
  RM-sep  o mesmo metodo de randomizacao com o nucleo do Ceron, exp(-|dx|/20 - |dy|/2): como o nucleo e o produto de
          dois exponenciais 1D, a densidade espectral e o produto de duas Cauchy, k_x ~ Cauchy(0, 1/20) e
          k_y ~ Cauchy(0, 1/2) (a funcao caracteristica da Cauchy e exp(-|h|/L)); o GSTools nao tem esse modelo.
Medidas: a covariancia-alvo C(p, q) = exp(-|dx|/20 - |dy|/2); a covariancia que cada gerador produz -- exata para a KL
(sum lambda_k f_k(p) f_k(q)) e para o GS-exp (o nucleo dele), empirica (R realizacoes) para os de randomizacao; a
variancia pontual; a curtose (normalidade) num ponto; e o tempo por realizacao.

Uso: python3 compara_gstools.py [R=2000]
"""
import sys
import time

import numpy as np

import klfechada as kf

R = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
rng = np.random.default_rng(1)

# pontos: grade de 0,25 m no talude
xg, yg = np.meshgrid(np.arange(0.125, 25, 0.25), np.arange(0.125, 10, 0.25))
dentro = kf.no_talude(xg, yg)
X, Y = xg[dentro], yg[dentro]
n = len(X)
# pares para a covariancia: 400 aleatorios e 400 proximos (ate 3 m em x e 1 m em y)
ia = rng.integers(0, n, 400)
ib = rng.integers(0, n, 400)
c = rng.integers(0, n, 4000)
d = np.array([np.argmin((X - X[i] - rng.uniform(-3, 3)) ** 2 + (Y - Y[i] - rng.uniform(-1, 1)) ** 2) for i in c[:400]])
pa, pb = np.concatenate([ia, c[:400]]), np.concatenate([ib, d])
dx, dy = np.abs(X[pa] - X[pb]), np.abs(Y[pa] - Y[pb])
alvo = np.exp(-dx / 20 - dy / 2)
gs_nucleo = np.exp(-np.sqrt((dx / 20) ** 2 + (dy / 2) ** 2))


def resumo(nome, cov, var, t, kurt=None):
    e = np.abs(cov - alvo)
    print("%-8s cov: erro medio %.4f, max %.4f (proximos: medio %.4f) | variancia %.3f a %.3f | %s%.4f s/realizacao"
          % (nome, e.mean(), e.max(), e[400:].mean(), var.min(), var.max(),
             "" if kurt is None else "curtose %.2f | " % kurt, t))


# ---- KL (exata para a covariancia truncada) ----
kl = kf.KL2D(2000)
F, G = kl.fx(X), kl.gy(Y)
Phi = F[kl.ix] * G[kl.iy] * np.sqrt(kl.lam)[:, None]          # (M, n)
cov_kl = np.einsum('kp,kp->p', Phi[:, pa], Phi[:, pb])
var_kl = (Phi ** 2).sum(0)
t0 = time.time()
for _ in range(20):
    z = rng.standard_normal(2000) @ Phi
t_kl = (time.time() - t0) / 20

# ---- randomizacao (Hesse 2014), N modos ----
N = 1000


def rm(kx, ky):
    ph = np.outer(kx, X) + np.outer(ky, Y)                      # (N, n)
    z1, z2 = rng.standard_normal(N), rng.standard_normal(N)
    return np.sqrt(1. / N) * (z1 @ np.cos(ph) + z2 @ np.sin(ph))


def ensemble(sorteia):
    S = np.zeros(len(pa)); S2 = np.zeros(n); p0 = []
    t0 = time.time()
    for _ in range(R):
        u = rm(*sorteia())
        S += u[pa] * u[pb]; S2 += u ** 2; p0.append(u[n // 2])
    t = (time.time() - t0) / R
    p0 = np.array(p0)
    return S / R, S2 / R, t, np.mean(p0 ** 4) / np.mean(p0 ** 2) ** 2


def sep():
    return rng.standard_cauchy(N) / 20., rng.standard_cauchy(N) / 2.


cov_rm, var_rm, t_rm, k_rm = ensemble(sep)

try:
    import gstools as gs
    mod = gs.Exponential(dim=2, var=1., len_scale=[20., 2.])
    gs_ok = True
except ImportError:
    gs_ok = False

print("talude: %d pontos (grade de 0,25 m), %d pares; R = %d realizacoes nos de randomizacao (ruido amostral da "
      "covariancia ~ %.3f)" % (n, len(pa), R, 1 / np.sqrt(R)))
resumo("KL", cov_kl, var_kl, t_kl)
resumo("RM-sep", cov_rm, var_rm, t_rm, k_rm)
if gs_ok:
    srf = gs.SRF(mod, mode_no=N)
    S = np.zeros(len(pa)); S2 = np.zeros(n); p0 = []
    t0 = time.time()
    for s in range(R):
        u = srf((X, Y), seed=s)
        S += u[pa] * u[pb]; S2 += u ** 2; p0.append(u[n // 2])
    t_gs = (time.time() - t0) / R
    p0 = np.array(p0)
    resumo("GS-exp", S / R, S2 / R, t_gs, np.mean(p0 ** 4) / np.mean(p0 ** 2) ** 2)
    e = np.abs(gs_nucleo - alvo)
    print("nucleo do GS-exp contra o do Ceron (sem amostragem): diferenca media %.4f, maxima %.4f" % (e.mean(), e.max()))
