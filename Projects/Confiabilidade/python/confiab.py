"""confiab.py: FORM (HL-RF com busca de Armijo), Monte Carlo e amostragem por importancia no espaco normal padrao, com a
transformacao de Nataf (marginais com correlacao). E o mesmo algoritmo do confiabilidade.wl (Mathematica), em Python, para
gerar os numeros e as figuras do documento (doc/confiabilidade.tex) e para conferir o codigo do Mathematica.

Convencoes (as do documento):
  X  variaveis fisicas, com marginais F_i e correlacao RX;
  Z  normais padrao correlacionadas, z_i = Phi^-1(F_i(x_i)), com correlacao RZ (a de Nataf);
  U  normais padrao independentes, z = L u, RZ = L L^T (Cholesky, L triangular inferior); u = L^-1 z.
A funcao de estado limite e g(x) <= 0 na falha; no espaco U, G(u) = g(x(u)).
"""
import numpy as np
from scipy import optimize, stats


# ---------------------------------------------------------------- marginais
class Normal:
    def __init__(self, media, dp):
        self.media, self.dp = float(media), float(dp)

    def x_de_z(self, z):
        return self.media + self.dp * np.asarray(z, float)

    def z_de_x(self, x):
        return (np.asarray(x, float) - self.media) / self.dp

    def __repr__(self):
        return "N(%g; %g)" % (self.media, self.dp)


class LogNormal:
    """lognormal pela media e pelo coeficiente de variacao: ln X ~ N(lam, zeta^2)"""

    def __init__(self, media, cov):
        self.media, self.cov = float(media), float(cov)
        self.dp = self.media * self.cov
        self.zeta = np.sqrt(np.log1p(self.cov ** 2))
        self.lam = np.log(self.media) - 0.5 * self.zeta ** 2

    def x_de_z(self, z):
        return np.exp(self.lam + self.zeta * np.asarray(z, float))

    def z_de_x(self, x):
        return (np.log(np.asarray(x, float)) - self.lam) / self.zeta

    def __repr__(self):
        return "LN(%g; CoV %g)" % (self.media, self.cov)


class Generica:
    """qualquer distribuicao do scipy.stats (congelada), pela transformacao exata x = F^-1(Phi(z))"""

    def __init__(self, dist):
        self.d = dist
        self.media, self.dp = dist.mean(), dist.std()

    def x_de_z(self, z):
        z = np.asarray(z, float)
        # cauda superior pela funcao de sobrevivencia, para nao perder precisao com Phi(z) perto de 1
        return np.where(z <= 0, self.d.ppf(stats.norm.cdf(z)), self.d.isf(stats.norm.sf(z)))

    def z_de_x(self, x):
        x = np.asarray(x, float)
        c = self.d.cdf(x)
        return np.where(c <= 0.5, stats.norm.ppf(c), -stats.norm.ppf(self.d.sf(x)))


# ---------------------------------------------------------------- correlacao de Nataf
def rho_x_de_rho_z(di, dj, r, n=96):
    """correlacao de X_i, X_j quando a de Z_i, Z_j e r (Gauss-Hermite 2D)"""
    t, w = np.polynomial.hermite_e.hermegauss(n)
    w = w / np.sqrt(2 * np.pi)
    a = (di.x_de_z(t) - di.media) / di.dp                                     # (n,)
    zj = r * t[:, None] + np.sqrt(1 - r * r) * t[None, :]                     # (n, n)
    b = (dj.x_de_z(zj) - dj.media) / dj.dp
    return float(np.einsum("i,j,i,ij->", w, w, a, b))


def rho_nataf(di, dj, rho):
    """correlacao no espaco Z (Liu e Der Kiureghian 1986): formas fechadas para normal e lognormal; nas outras, a
    equacao integral resolvida numericamente"""
    if rho == 0:
        return 0.
    if isinstance(di, Normal) and isinstance(dj, Normal):
        return rho
    if isinstance(di, LogNormal) and isinstance(dj, LogNormal):
        return np.log1p(rho * di.cov * dj.cov) / (di.zeta * dj.zeta)
    if isinstance(di, Normal) and isinstance(dj, LogNormal):
        return rho * dj.cov / dj.zeta
    if isinstance(di, LogNormal) and isinstance(dj, Normal):
        return rho * di.cov / di.zeta
    return optimize.brentq(lambda r: rho_x_de_rho_z(di, dj, r) - rho, -0.9999, 0.9999)


class Modelo:
    def __init__(self, dists, RX=None):
        self.dists, self.n = list(dists), len(dists)
        self.RX = np.eye(self.n) if RX is None else np.asarray(RX, float)
        self.RZ = np.array([[1. if i == j else rho_nataf(self.dists[i], self.dists[j], self.RX[i, j])
                             for j in range(self.n)] for i in range(self.n)])
        self.L = np.linalg.cholesky(self.RZ)                                    # RZ = L L^T, L triangular inferior
        assert np.abs(self.L @ self.L.T - self.RZ).max() < 1e-12

    def x_de_u(self, u):
        """u: (n,) ou (N, n) -> x com a mesma forma"""
        z = np.asarray(u, float) @ self.L.T
        return np.stack([d.x_de_z(z[..., i]) for i, d in enumerate(self.dists)], axis=-1)

    def u_de_x(self, x):
        x = np.asarray(x, float)
        z = np.stack([d.z_de_x(x[..., i]) for i, d in enumerate(self.dists)], axis=-1)
        return np.linalg.solve(self.L, z.T).T


# ---------------------------------------------------------------- FORM
def gradiente(G, u, h):
    """diferencas centrais (2n avaliacoes de G)"""
    n = len(u)
    d = np.empty(n)
    for i in range(n):
        e = np.zeros(n)
        e[i] = h
        d[i] = (G(u + e) - G(u - e)) / (2 * h)
    return d


def form(G, n, passo=1e-5, tol_g=1e-6, tol_u=1e-4, max_iter=100, armijo=True):
    """HL-RF com busca de Armijo na funcao de merito m(u) = |u|^2/2 + c |G(u)| (Zhang e Der Kiureghian 1995).
    G: funcao de um vetor u (n,). Devolve beta, u*, alfa, pf = Phi(-beta), iteracoes e numero de avaliacoes de G."""
    u = np.zeros(n)
    g = G(u)
    nav = 1
    g0 = abs(g) if g != 0 else 1.
    hist = []
    conv = False
    for it in range(1, max_iter + 1):
        dg = gradiente(G, u, passo)
        nav += 2 * n
        ndg = np.linalg.norm(dg)
        alfa = -dg / ndg
        hist.append((it, u.copy(), g, np.linalg.norm(u)))
        if abs(g) / g0 < tol_g and np.linalg.norm(u - (alfa @ u) * alfa) < tol_u:
            conv = True
            break
        uhl = (dg @ u - g) / ndg ** 2 * dg                                     # ponto HL-RF
        d = uhl - u
        lam = 1.
        if armijo:
            c = (2 * np.linalg.norm(u) + 1.) / ndg
            m0 = 0.5 * u @ u + c * abs(g)
            dm = (u + c * np.sign(g) * dg) @ d                                  # derivada direcional do merito
            while True:
                un = u + lam * d
                gn = G(un)
                nav += 1
                if 0.5 * un @ un + c * abs(gn) <= m0 + 1e-4 * lam * dm or lam < 2. ** -10:
                    break
                lam *= 0.5
        else:
            un = u + d
            gn = G(un)
            nav += 1
        u, g = un, gn
    beta = float(alfa @ u)
    return dict(beta=beta, u=u, alfa=alfa, pf=stats.norm.cdf(-beta), iteracoes=it, avaliacoes=nav,
                convergiu=conv, historico=hist)


# ---------------------------------------------------------------- Monte Carlo e amostragem por importancia
def monte_carlo(Gv, n, N, rng, bloco=2_000_000):
    """Gv: funcao vetorizada de (m, n) -> (m,). Devolve pf, CoV e os indicadores"""
    I = np.empty(N, bool)
    for a in range(0, N, bloco):
        m = min(bloco, N - a)
        I[a:a + m] = Gv(rng.standard_normal((m, n))) <= 0
    pf = I.mean()
    cov = np.sqrt((1 - pf) / (N * pf)) if pf > 0 else np.inf
    return dict(pf=pf, cov=cov, q=I.astype(float), N=N)


def importancia(Gv, ustar, N, rng):
    """amostragem por importancia com h = N(u*, I): q = I(G <= 0) phi(v)/phi(v - u*)"""
    ustar = np.asarray(ustar, float)
    v = ustar + rng.standard_normal((N, len(ustar)))
    lw = -v @ ustar + 0.5 * ustar @ ustar
    q = (Gv(v) <= 0) * np.exp(lw)
    pf = q.mean()
    cov = q.std(ddof=1) / np.sqrt(N) / pf if pf > 0 else np.inf
    return dict(pf=pf, cov=cov, q=q, N=N, falhas=int((q > 0).sum()))


def media_acumulada(q):
    k = np.arange(1, len(q) + 1)
    m = np.cumsum(q) / k
    s2 = np.maximum(np.cumsum(q * q) / k - m * m, 0.)
    return m, np.sqrt(s2 / k)


def beta_de_pf(pf):
    return -stats.norm.ppf(pf) if 0 < pf < 1 else np.nan
