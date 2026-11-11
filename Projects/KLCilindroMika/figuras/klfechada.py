"""klfechada.py: a expansao de Karhunen-Loeve (KL) em forma fechada do nucleo exponencial separavel
    C(x, x') = exp(-|x - x'| / Lx - |y - y'| / Ly)
no retangulo, em numpy: a mesma conta de ../../SlopeRandomField/KLSeparavel.h (klsep::KL2D), usada pelo
SlopeRandomField no retangulo que contem a malha, e as integrais exatas no talude de Vargas Ceron et al. (2025).

1D: nucleo exp(-c |s - t|), c = 1/L, em [-a, a]
    pares:   omega tan(omega a) = c,      f(s) = cos(omega s) / sqrt(a + sin(2 omega a) / (2 omega))
    impares: omega + c tan(omega a) = 0,  f(s) = sin(omega s) / sqrt(a - sin(2 omega a) / (2 omega))
    lambda = 2 c / (omega^2 + c^2)
2D: os autopares sao os produtos lambda_i^x lambda_j^y, f_i(x) g_j(y); a KL truncada usa os M maiores produtos.

Tambem: o gerador de normais do SlopeRandomField (xoshiro256** semeado por splitmix64, Box-Muller) e a transformacao
lognormal, para refazer em Python o campo que o programa poe nos pontos de integracao.
"""
import numpy as np

# talude de Vargas Ceron et al. (2025, Secao 6) e o retangulo que o contem (a caixa dos nos das malhas)
TALUDE = np.array([[0., 0.], [25., 0.], [25., 5.], [15., 5.], [10., 10.], [0., 10.]])
X0, Y0, DX, DY = 0., 0., 25., 10.
AREA_TALUDE, AREA_RET = 187.5, 250.
LX, LY = 20., 2.


def no_talude(x, y):
    """pontos do talude: 0 <= x <= 25, 0 <= y <= 10 e, acima de y = 5, a esquerda da face x + y = 20"""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    return (x >= 0) & (x <= 25) & (y >= 0) & (y <= 10) & ((y <= 5) | (x + y <= 20))


def _bissecao(g, lo, hi, it=200):
    """raizes de g (crescente) em (lo, hi), vetorizado, como detalhe::Bissecao de KLSeparavel.h"""
    lo = np.array(lo, float)
    hi = np.array(hi, float)
    glo = g(lo)
    for _ in range(it):
        m = 0.5 * (lo + hi)
        gm = g(m)
        mesmo = (gm > 0) == (glo > 0)
        lo = np.where(mesmo, m, lo)
        glo = np.where(mesmo, gm, glo)
        hi = np.where(mesmo, hi, m)
        if np.all(hi - lo < 1e-15 * np.maximum(1., hi)):
            break
    return 0.5 * (lo + hi)


class Modos1D:
    """os n primeiros autopares do nucleo exp(-|s - t| / L) em [-a, a], a = comprimento / 2, em ordem decrescente de
    lambda (o modo i e par, cosseno, se i e par: as raizes alternam par, impar, par, ...)"""

    def __init__(self, comprimento, L, n):
        self.a = a = 0.5 * comprimento
        self.L = L
        self.c = c = 1. / L
        eps = 1e-12
        npar, nimp = (n + 1) // 2, n // 2
        k = np.arange(1, npar + 1)
        wp = _bissecao(lambda x: x * np.tan(x * a) - c, ((k - 1) * np.pi + eps) / a, ((k - 0.5) * np.pi - eps) / a)
        k = np.arange(1, nimp + 1)
        wi = _bissecao(lambda x: x + c * np.tan(x * a), ((k - 0.5) * np.pi + eps) / a, (k * np.pi - eps) / a)
        w = np.empty(n)
        w[0::2], w[1::2] = wp, wi
        self.omega = w
        self.par = (np.arange(n) % 2 == 0)
        self.lam = 2. * c / (w * w + c * c)
        s2 = np.sin(2. * w * a) / (2. * w)
        self.norma = 1. / np.sqrt(np.where(self.par, a + s2, a - s2))

    def avalia(self, i, s):
        """f_i(s), s medido a partir do centro do intervalo; i escalar ou vetor (linhas), s vetor (colunas)"""
        i = np.atleast_1d(i)
        arg = self.omega[i][:, None] * np.asarray(s, float)[None, :]
        return self.norma[i][:, None] * np.where(self.par[i][:, None], np.cos(arg), np.sin(arg))


class KL2D:
    """os M maiores produtos lambda_i^x lambda_j^y no retangulo [x0, x0 + Dx] x [y0, y0 + Dy], como klsep::KL2D
    (candidatos (i + 1)(j + 1) <= M, ordem decrescente de lambda)"""

    def __init__(self, M, x0=X0, y0=Y0, Dx=DX, Dy=DY, Lx=LX, Ly=LY):
        self.x0, self.y0, self.Dx, self.Dy = x0, y0, Dx, Dy
        self.mx = Modos1D(Dx, Lx, M)
        self.my = Modos1D(Dy, Ly, M)
        ii, jj = [], []
        for i in range(M):
            nj = min(M, M // (i + 1))
            ii.append(np.full(nj, i))
            jj.append(np.arange(nj))
        ii = np.concatenate(ii)
        jj = np.concatenate(jj)
        lam = self.mx.lam[ii] * self.my.lam[jj]
        o = np.argsort(-lam, kind='stable')[:M]
        self.ix, self.iy, self.lam = ii[o], jj[o], lam[o]
        self.M = M

    def fx(self, x):
        """tabela f_i(x_p) dos modos 1D em x usados (linhas: indice 1D 0..max)"""
        return self.mx.avalia(np.arange(self.ix.max() + 1), np.asarray(x) - (self.x0 + 0.5 * self.Dx))

    def gy(self, y):
        return self.my.avalia(np.arange(self.iy.max() + 1), np.asarray(y) - (self.y0 + 0.5 * self.Dy))

    def modo(self, k, x, y):
        """modo 2D k (0..M-1) nos pontos (x, y)"""
        return (self.mx.avalia(self.ix[k], np.asarray(x) - (self.x0 + 0.5 * self.Dx))[0] *
                self.my.avalia(self.iy[k], np.asarray(y) - (self.y0 + 0.5 * self.Dy))[0])

    def campo(self, xi, x, y, M=None):
        """Z(x) = sum_k sqrt(lambda_k) xi_k f_k(x) nos pontos (x, y), com os M primeiros termos"""
        M = self.M if M is None else M
        F, G = self.fx(x), self.gy(y)
        z = np.zeros(len(np.atleast_1d(x)))
        for k in range(M):
            z += np.sqrt(self.lam[k]) * xi[k] * F[self.ix[k]] * G[self.iy[k]]
        return z

    def variancia_grade(self, xg, yg, Ms):
        """sigma_M^2 = sum_{k<=M} lambda_k f_k^2 na grade tensorial xg x yg, para cada M de Ms: lista de (ny, nx)"""
        F2, G2 = self.fx(xg) ** 2, self.gy(yg) ** 2
        res = []
        for M in Ms:
            Lam = np.zeros((F2.shape[0], G2.shape[0]))
            np.add.at(Lam, (self.ix[:M], self.iy[:M]), self.lam[:M])
            res.append((F2.T @ Lam @ G2).T)
        return res

    def integral_talude(self):
        """integral no talude de f_i(x)^2 g_j(y)^2 de cada termo, exata, para o retangulo [0,25] x [0,10]:
        1/2 na faixa 0 <= y <= 5 (g^2 e par em t = y - 5) mais a faixa de cima, 5 <= y <= 10, 0 <= x <= 20 - y, em
        forma fechada (integrais de polinomio de grau 1 vezes funcoes trigonometricas)"""
        assert (self.x0, self.y0, self.Dx, self.Dy) == (X0, Y0, DX, DY)
        w, nu = self.mx.omega[self.ix], self.my.omega[self.iy]
        sx = np.where(self.mx.par[self.ix], 1., -1.)
        sy = np.where(self.my.par[self.iy], 1., -1.)
        n2, m2 = self.mx.norma[self.ix] ** 2, self.my.norma[self.iy] ** 2

        def I(alfa, beta):  # int_0^5 sin(alfa - beta t) dt
            return 5. * np.sin(alfa - 2.5 * beta) * np.sinc(2.5 * beta / np.pi)

        k = 2. * nu
        J = 10. * np.sin(5. * k) / k + 2. * np.sin(2.5 * k) ** 2 / k ** 2  # int_0^5 (15 - t) cos(k t) dt
        s25 = np.sin(25. * w) / (2. * w)
        U = (62.5 + 5. * sx * s25 + sx / (2. * w) * I(5. * w, 2. * w) + sy * J
             + sx * sy * s25 * np.sin(10. * nu) / (2. * nu)
             + sx * sy / (2. * w) * 0.5 * (I(5. * w, 2. * (w - nu)) + I(5. * w, 2. * (w + nu))))
        return 0.5 + 0.25 * n2 * m2 * U

    def erro_talude(self):
        """epsilon_M = 1 - media de sigma_M^2 no talude, M = 1..self.M"""
        return 1. - np.cumsum(self.lam * self.integral_talude()) / AREA_TALUDE


# ---- gerador do SlopeRandomField (KLSeparavel.h, klsep::Normal) ----
_M64 = (1 << 64) - 1


class Normal:
    """xoshiro256** semeado por splitmix64; normal por Box-Muller (devolve r cos t e depois r sin t)"""

    def __init__(self, semente):
        x = semente & _M64
        self.s = []
        for _ in range(4):
            x = (x + 0x9e3779b97f4a7c15) & _M64
            z = x
            z = ((z ^ (z >> 30)) * 0xbf58476d1ce4e5b9) & _M64
            z = ((z ^ (z >> 27)) * 0x94d049bb133111eb) & _M64
            self.s.append(z ^ (z >> 31))
        self.tem, self.guarda = False, 0.

    @staticmethod
    def _rot(x, k):
        return ((x << k) | (x >> (64 - k))) & _M64

    def _proximo(self):
        s = self.s
        r = (self._rot((s[1] * 5) & _M64, 7) * 9) & _M64
        t = (s[1] << 17) & _M64
        s[2] ^= s[0]
        s[3] ^= s[1]
        s[1] ^= s[2]
        s[0] ^= s[3]
        s[2] ^= t
        s[3] = self._rot(s[3], 45)
        return r

    def _uniforme(self):
        return ((self._proximo() >> 11) + 0.5) * (1.0 / 9007199254740992.0)

    def __call__(self):
        if self.tem:
            self.tem = False
            return self.guarda
        u1, u2 = self._uniforme(), self._uniforme()
        r, t = np.sqrt(-2. * np.log(u1)), 2. * np.pi * u2
        self.guarda, self.tem = r * np.sin(t), True
        return r * np.cos(t)

    def vetor(self, n):
        return np.array([self() for _ in range(n)])


def xi_slope(campo, real=0, semente=20250416, M=2000):
    """as M normais do campo (0: c', 1: phi', 2: k) da realizacao real, como no main.cpp do SlopeRandomField"""
    s = (semente * 0x9E3779B97F4A7C15 + 3 * real + campo) & _M64
    return Normal(s).vetor(M)


def lognormal(media, cov):
    """(mu_ln, sigma_ln) da normal de ln H (Vargas Ceron et al. 2025, eqs. 65-66)"""
    s = np.sqrt(np.log(1. + cov * cov))
    return np.log(media / np.sqrt(1. + cov * cov)), s


if __name__ == '__main__':
    kl = KL2D(2000)
    print('lambda 1, 10, 20:', kl.lam[0], kl.lam[9], kl.lam[19])
    print('fracao retida no retangulo, M = 2000:', kl.lam.sum() / AREA_RET)
    e = kl.erro_talude()
    print('epsilon_M no talude, M = 50, 200, 1000, 2000:', e[[49, 199, 999, 1999]])
