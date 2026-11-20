#!/usr/bin/env python3
"""campos.py: campos aleatorios lognormais (c', phi', k) do Monte Carlo do Ceron (2025, Secao 6) pela expansao de
Karhunen-Loeve ANALITICA do nucleo exponencial separavel, o mesmo esquema do C++
Projects/SlopeRandomFieldRefinement/KLSeparavel.h + RandomField.cpp:

    C(x1, x2) = exp(-|x1 - x2|/Lx - |y1 - y2|/Ly)            (artigo, Eq. 64; RandomField.h: a mesma convencao,
                                                               sem o fator 2)
  * retangulo [x0, x0 + Dx] x [y0, y0 + Dy] que contem o dominio; em cada direcao os autopares do nucleo
    exp(-|s - t|/L) em [-a, a] (raizes transcendentais: pares w tan(w a) = 1/L, impares w + tan(w a)/L = 0),
    lambda = 2 (1/L) / (w^2 + 1/L^2);
  * termos 2D = produtos lambda_i^x lambda_j^y; ficam os M maiores (M = 2000 por campo), Z = sum_k sqrt(lambda_k)
    xi_k f_k(x) g_k(y); fracao de variancia retida = sum lambda_k / (Dx Dy);
  * lognormal (Eqs. 65-66): H = exp(mu_ln + sigma_ln Z), sigma_ln = sqrt(ln(1 + CoV^2)), mu_ln = ln(media) -
    sigma_ln^2 / 2. Os tres campos sao independentes e usam os mesmos modos, cada um com os seus xi.

Avaliacao: f e g so dos modos 1D usados. Num conjunto de pontos soltos guarda Fx (n x nx) e Gy (n x ny) e
Z = sum_a Fx[:, a] * (Gy @ A^T)[:, a], com A[a, b] = soma de sqrt(lambda_k) xi_k dos termos k = (a, b). Numa grade
tensorial (xs, ys): Z = Fx(xs) @ A @ Gy(ys)^T. As duas formas custam milissegundos por campo.

Sorteio: realizacao i -> numpy.random.default_rng(semente + i), 3 M normais na ordem (c, phi, k)."""
from __future__ import annotations

import math
import time

import numpy as np

# ---------------------------------------------------------------------------------------------------------------------
# Dados do caso de referencia (Tabela 2 do artigo)
# ---------------------------------------------------------------------------------------------------------------------
LX, LY = 20.0, 2.0
STATS = {"c": (10.0, 0.30), "phi": (30.0, 0.10), "k": (1.0, 0.60)}   # media, CoV (phi em graus; k relativo)
M_PADRAO = 2000


class Lognormal:
    def __init__(self, media, cov):
        self.sigma_ln = math.sqrt(math.log(1.0 + cov * cov))
        self.mu_ln = math.log(media / math.sqrt(1.0 + cov * cov))

    def __call__(self, z):
        return np.exp(self.mu_ln + self.sigma_ln * z)


def modos_1d(comprimento, L, n):
    """Os n primeiros autopares de exp(-|s - t|/L) em [-a, a] (a = comprimento/2), em ordem decrescente de lambda.
    Devolve (lambda, omega, par, norma). Bissecao vetorizada, como a do KLSeparavel.h."""
    a = 0.5 * comprimento
    c = 1.0 / L
    k = np.arange(1, n // 2 + 2)
    eps = 1e-12

    def bissecao(g, lo, hi):
        glo = g(lo)
        for _ in range(200):
            m = 0.5 * (lo + hi)
            gm = g(m)
            mesmo = (gm > 0) == (glo > 0)
            lo = np.where(mesmo, m, lo)
            glo = np.where(mesmo, gm, glo)
            hi = np.where(mesmo, hi, m)
            if np.all(hi - lo < 1e-15 * np.maximum(1.0, hi)):
                break
        return 0.5 * (lo + hi)

    wp = bissecao(lambda x: x * np.tan(x * a) - c, ((k - 1) * np.pi + eps) / a, ((k - 0.5) * np.pi - eps) / a)
    wi = bissecao(lambda x: x + c * np.tan(x * a), ((k - 0.5) * np.pi + eps) / a, (k * np.pi - eps) / a)
    w = np.empty(2 * len(k))
    w[0::2], w[1::2] = wp, wi
    par = np.tile([1, 0], len(k))
    w, par = w[:n], par[:n]
    lam = 2.0 * c / (w * w + c * c)
    o = np.argsort(-lam, kind="stable")
    w, par, lam = w[o], par[o], lam[o]
    s2 = np.sin(2.0 * w * a) / (2.0 * w)
    norma = 1.0 / np.sqrt(np.where(par == 1, a + s2, a - s2))
    return lam, w, par, norma


class KL2D:
    """Expansao 2D truncada nos M maiores produtos (o KL2D do KLSeparavel.h)."""

    def __init__(self, x0, y0, Dx, Dy, Lx=LX, Ly=LY, M=M_PADRAO):
        self.x0, self.y0, self.Dx, self.Dy = x0, y0, Dx, Dy
        self.Lx, self.Ly = Lx, Ly
        n1 = M
        lx, self.wx, self.px, self.nx = modos_1d(Dx, Lx, n1)
        ly, self.wy, self.py, self.ny = modos_1d(Dy, Ly, n1)
        # candidatos (i, j) com (i+1)(j+1) <= M (os outros nao podem estar entre os M maiores)
        ii, jj = [], []
        for i in range(n1):
            jmax = min(n1, M // (i + 1))
            if jmax <= 0:
                break
            ii.append(np.full(jmax, i))
            jj.append(np.arange(jmax))
        ii, jj = np.concatenate(ii), np.concatenate(jj)
        lam = lx[ii] * ly[jj]
        o = np.argsort(-lam, kind="stable")[:M]
        self.ix, self.iy, self.lam = ii[o], jj[o], lam[o]
        self.M = len(self.lam)
        self.raiz = np.sqrt(self.lam)
        self.fracao = float(self.lam.sum() / (Dx * Dy))
        # modos 1D usados e a coluna de cada termo
        self.ux, self.cx = np.unique(self.ix, return_inverse=True)
        self.uy, self.cy = np.unique(self.iy, return_inverse=True)

    def Fx(self, x):
        """f_a(x) dos modos x usados: (len(x), nux)."""
        s = np.asarray(x, float)[:, None] - (self.x0 + 0.5 * self.Dx)
        w, p, nm = self.wx[self.ux], self.px[self.ux], self.nx[self.ux]
        return nm * np.where(p == 1, np.cos(w * s), np.sin(w * s))

    def Gy(self, y):
        s = np.asarray(y, float)[:, None] - (self.y0 + 0.5 * self.Dy)
        w, p, nm = self.wy[self.uy], self.py[self.uy], self.ny[self.uy]
        return nm * np.where(p == 1, np.cos(w * s), np.sin(w * s))

    def coef(self, xi):
        """A[a, b] = sum dos sqrt(lambda_k) xi_k dos termos com modos (a, b)."""
        A = np.zeros((len(self.ux), len(self.uy)))
        np.add.at(A, (self.cx, self.cy), self.raiz * xi)
        return A

    def variancia(self, Fx, Gy):
        """Variancia truncada sum lambda_k f_k^2 g_k^2 nos pontos (Fx, Gy da mesma lista de pontos)."""
        return ((Fx[:, self.cx] * Gy[:, self.cy]) ** 2 * self.lam).sum(axis=1)


class PontosKL:
    """Base guardada para pontos soltos (ex.: centroides da malha do Darcy)."""

    def __init__(self, kl: KL2D, x, y):
        self.kl = kl
        self.Fx = kl.Fx(x)
        self.Gy = kl.Gy(y)

    def Z(self, xi):
        A = self.kl.coef(xi)
        return np.einsum("pa,pa->p", self.Fx, self.Gy @ A.T)


class GradeKL:
    """Base guardada para uma grade tensorial xs x ys (devolve Z[ix, iy])."""

    def __init__(self, kl: KL2D, xs, ys):
        self.kl = kl
        self.Fx = kl.Fx(xs)
        self.Gy = kl.Gy(ys)

    def Z(self, xi):
        return self.Fx @ (self.kl.coef(xi) @ self.Gy.T)


def sortear_xi(semente, i, M, ncampos=3):
    """xi da realizacao i: (ncampos, M), gerador default_rng(semente + i); ordem dos campos (c, phi, k)."""
    rng = np.random.default_rng(int(semente) + int(i))
    return rng.standard_normal((ncampos, M))


# ---------------------------------------------------------------------------------------------------------------------
# Malhas, Darcy permanente e a realizacao completa
# ---------------------------------------------------------------------------------------------------------------------
GW = 9.81


def _nos_graduados(a, b, h_fino, cap, focos, taxa=0.25):
    """Nos de a a b: espacamento h_fino + taxa * (distancia ao foco mais proximo), limitado por cap(x)."""
    xs = [a]
    while xs[-1] < b - 1e-9:
        d = min(abs(xs[-1] - f) for f in focos)
        h = min(cap(xs[-1]), h_fino + taxa * d)
        xs.append(min(b, xs[-1] + h))
    return np.array(xs)


def grade_tensorial(xs, ys, tal):
    """Triangulos P1 da grade tensorial xs x ys (y para baixo) so no solo: cada celula vira dois triangulos pela
    diagonal (i, j)-(i+1, j+1) (paralela a face com beta = 45 graus); nas celulas da face fica so o triangulo de
    baixo (y > x). Exige xs e ys com o mesmo espacamento uniforme em [0, xT] e [0, H].
    Devolve (xy, tri, no_solo, celula_do_triangulo)."""
    nx, ny = len(xs), len(ys)
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    xy = np.c_[X.ravel(), Y.ravel()]
    idx = np.arange(nx * ny).reshape(nx, ny)
    tol = 1e-9
    xa, xb = xs[:-1, None], xs[1:, None]
    ya = ys[None, :-1]
    ysb = tal.y_sup(xb)                       # superficie no lado direito da celula (y_sup nao decresce)
    cheia = ya >= ysb - tol
    face = (np.abs(ya - xa) < tol) & (xa >= -tol) & (xb <= tal.xT + tol)
    a, b, c, d = idx[:-1, :-1], idx[1:, :-1], idx[1:, 1:], idx[:-1, 1:]
    cel = np.arange((nx - 1) * (ny - 1)).reshape(nx - 1, ny - 1)
    tA = np.c_[a[cheia], b[cheia], c[cheia]]          # s >= t
    tB = np.c_[a[cheia | face], c[cheia | face], d[cheia | face]]
    tri = np.r_[tA, tB]
    cel_tri = np.r_[cel[cheia], cel[cheia | face]]
    usados = np.unique(tri)
    no_solo = np.zeros(nx * ny, bool)
    no_solo[usados] = True
    return xy, tri, no_solo, cel_tri


class Problema:
    """Geometria, malhas e bases do Monte Carlo do Ceron (Secao 6, configuracao de referencia).

    * dominio do Darcy: 14H atras da crista, 2,8H a frente do pe, 8,4H abaixo do pe (Le = 70, Ld = 14, D = 42 m com
      H = 5), o contorno do Ceron (crista u = 0, face u = -gw y ate hw e -gw hw abaixo, terreno do pe -gw hw, lateral
      esquerda u = 0, base u = 0, lateral direita sem fluxo);
    * malha do Darcy: grade tensorial graduada (hD_fino no talude, ate hD_janela na janela dos mecanismos, ate 2 m
      fora), k constante por celula (o valor do campo no centro da celula, os dois triangulos);
    * grade dos campos (c, phi, u): grade uniforme de passo h na janela x em [x_jan, xT + Ld], y em [0, y_jan];
      c e phi = o KL nos nos (P1), u = o P1 do Darcy interpolado nos nos. E nela que o cinematico trabalha (o talude
      do mecanismo tem x_esq = x_jan e y_base = y_jan).
    * KL: retangulo [-Le, xT + Ld] x [0, H + D], M termos."""

    def __init__(self, h=0.05, x_jan=-25.0, y_jan=15.0, hD_fino=0.0625, hD_janela=0.25, M=M_PADRAO,
                 H=5.0, beta_graus=45.0, Le=70.0, Ld=14.0, D=42.0, hw=None, gamma=20.0, stats=None):
        import cinematico as kin
        self.H, self.beta = H, math.radians(beta_graus)
        self.hw = H if hw is None else hw
        self.gamma_ef = gamma - GW
        self.stats = dict(STATS if stats is None else stats)
        self.dom = kin.Talude(H=H, beta=self.beta, x_esq=-Le, x_dir=H / math.tan(self.beta) + Ld, y_base=H + D)
        xT = self.dom.xT
        assert abs(xT - H) < 1e-12, "a grade tensorial exige beta = 45 graus"
        # ---- KL ----
        self.kl = KL2D(-Le, 0.0, Le + xT + Ld, H + D, M=M)
        self.lognormal = {nome: Lognormal(*self.stats[nome]) for nome in ("c", "phi", "k")}
        # ---- malha do Darcy ----
        capx = lambda x: hD_janela if x >= x_jan - 1e-9 else 2.0
        xs = np.unique(np.r_[_nos_graduados(-Le, 0.0, hD_fino, capx, [0.0])[:-1],
                             np.linspace(0.0, xT, int(round(xT / hD_fino)) + 1),
                             _nos_graduados(xT, xT + Ld, hD_fino, capx, [xT])[1:]])
        capy = lambda y: hD_janela if y <= y_jan + 1e-9 else 2.0
        ys = np.unique(np.r_[np.linspace(0.0, H, int(round(H / hD_fino)) + 1),
                             _nos_graduados(H, H + D, hD_fino, capy, [H])[1:]])
        self.xsD, self.ysD = xs, ys
        xyD, triD, soloD, celD = grade_tensorial(xs, ys, self.dom)
        # renumera so os nos do solo
        nova = -np.ones(len(xyD), int)
        nova[soloD] = np.arange(soloD.sum())
        self.xyD, self.triD = xyD[soloD], nova[triD]
        ic, jc = np.unravel_index(celD, (len(xs) - 1, len(ys) - 1))
        self.celD_i, self.celD_j = ic, jc                    # celula (i, j) de cada triangulo
        self._prepara_darcy()
        self.gradeKL_D = GradeKL(self.kl, 0.5 * (xs[1:] + xs[:-1]), 0.5 * (ys[1:] + ys[:-1]))
        # ---- grade dos campos ----
        self.h = h
        nxg = int(round((xT + Ld - x_jan) / h)) + 1
        nyg = int(round(y_jan / h)) + 1
        self.xg = x_jan + h * np.arange(nxg)
        self.yg = h * np.arange(nyg)
        assert abs(self.xg[-1] - (xT + Ld)) < 1e-9 and abs(round(-x_jan / h) * h + x_jan) < 1e-9
        self.gradeKL_G = GradeKL(self.kl, self.xg, self.yg)
        self.tal = kin.Talude(H=H, beta=self.beta, x_esq=x_jan, x_dir=xT + Ld, y_base=y_jan)
        self._prepara_interp()

    # ---- Darcy ----
    def _prepara_darcy(self):
        import scipy.sparse as sp
        xy, tri = self.xyD, self.triD
        n = len(xy)
        P = xy[tri]
        B = np.stack([P[:, 1] - P[:, 0], P[:, 2] - P[:, 0]], axis=2)
        A = 0.5 * np.abs(np.linalg.det(B))
        G = np.einsum("kji,jl->kil", np.linalg.inv(B), np.array([[-1.0, 1, 0], [-1.0, 0, 1]]))
        self._Ke = np.einsum("k,kia,kib->kab", A, G, G).reshape(len(tri), 9)
        x, y = xy[:, 0], xy[:, 1]
        tal, t, gw, hw = self.dom, 1e-7, GW, self.hw
        ud = np.full(n, np.nan)
        ud[np.abs(y - tal.y_base) < t] = 0.0                  # base u = 0
        ud[np.abs(x - tal.x_esq) < t] = 0.0                   # lateral esquerda u = 0 (direita: sem fluxo)
        sup = np.abs(y - tal.y_sup(x)) < t
        ud[sup & (x <= t)] = 0.0                              # crista
        face = sup & (x > t) & (x < tal.xT - t)
        ud[face] = np.where(y[face] <= hw, -gw * y[face], -gw * hw)
        ud[sup & (x >= tal.xT - t)] = -gw * hw                # terreno do pe
        D = ~np.isnan(ud)
        self._uD = np.where(D, ud, 0.0)
        self._livre = np.nonzero(~D)[0]
        self._fixo = np.nonzero(D)[0]
        # padrao esparso: K = S @ (k_e * Ke) com S fixo (de cada entrada elementar para a posicao CSR)
        rows = np.repeat(tri, 3, 1).ravel()
        cols = np.tile(tri, (1, 3)).ravel()
        # mapa entrada elementar -> posicao no csr (linhas e colunas em ordem lexicografica, como o csr ordenado)
        Kc = sp.coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, n)).tocsr()
        Kc.sum_duplicates()
        Kc.sort_indices()
        pos = np.empty(len(rows), int)
        ordem = np.lexsort((cols, rows))
        rs, cs = rows[ordem], cols[ordem]
        novo = np.r_[True, (rs[1:] != rs[:-1]) | (cs[1:] != cs[:-1])]
        pos[ordem] = np.cumsum(novo) - 1
        self._Kc, self._pos = Kc, pos
        mapa = -np.ones(n, int)
        mapa[self._livre] = np.arange(len(self._livre))
        self._mapa_livre = mapa

    def darcy(self, k_tri):
        """u (excesso de poropressao) nos nos do Darcy para k por triangulo; fatoracao esparsa nova (SuperLU)."""
        import scipy.sparse.linalg as spl
        K = self._Kc.copy()
        K.data = np.bincount(self._pos, weights=(k_tri[:, None] * self._Ke).ravel(), minlength=K.nnz)
        F, Dn = self._livre, self._fixo
        KFF = K[F][:, F].tocsc()
        rhs = -(K[F][:, Dn] @ self._uD[Dn])
        u = self._uD.copy()
        u[F] = spl.splu(KFF, permc_spec="MMD_AT_PLUS_A").solve(rhs)
        return u

    # ---- interpolacao Darcy -> grade ----
    def _prepara_interp(self):
        import scipy.sparse as sp
        from matplotlib.tri import Triangulation
        X, Y = np.meshgrid(self.xg, self.yg, indexing="ij")
        x, y = X.ravel(), Y.ravel()
        T = Triangulation(self.xyD[:, 0], self.xyD[:, 1], self.triD)
        k = np.asarray(T.get_trifinder()(x, y), int)
        solo = y >= self.tal.y_sup(x) - 1e-9
        perdidos = solo & (k < 0)
        if np.any(perdidos):     # na fronteira: o triangulo de centroide mais proximo
            from scipy.spatial import cKDTree
            cen = self.xyD[self.triD].mean(axis=1)
            k[perdidos] = cKDTree(cen).query(np.c_[x[perdidos], y[perdidos]])[1]
        ok = k >= 0
        P = self.xyD[self.triD[k[ok]]]
        M = np.concatenate([np.ones((ok.sum(), 3, 1)), P], axis=2)
        lam = np.einsum("pi,pij->pj", np.c_[np.ones(ok.sum()), x[ok], y[ok]], np.linalg.inv(M))
        lin = np.repeat(np.nonzero(ok)[0], 3)
        self._W = sp.csr_matrix((lam.ravel(), (lin, self.triD[k[ok]].ravel())), shape=(len(x), len(self.xyD)))
        self._ar = ~ok                                       # nos da grade no ar (fora do Darcy)
        self._u_ar = -GW * np.minimum(y, self.hw)            # extensao hidrostatica (so por arredondamento)
        self.forma_grade = X.shape

    # ---- realizacoes ----
    def campos_medios(self):
        """c, phi (rad) e u na grade com os valores medios (k uniforme)."""
        shp = self.forma_grade
        t0 = time.perf_counter()
        uD = self.darcy(np.ones(len(self.triD)))
        td = time.perf_counter() - t0
        return dict(c=np.full(shp, self.stats["c"][0]), phi=np.full(shp, math.radians(self.stats["phi"][0])),
                    u=self._u_grade(uD), uD=uD, k=np.ones(len(self.triD)), t_campos=0.0, t_darcy=td)

    def _u_grade(self, uD):
        u = self._W @ uD
        u[self._ar] = self._u_ar[self._ar]
        return u.reshape(self.forma_grade)

    def realizacao(self, i, semente):
        """Campos da realizacao i: c (kPa), phi (rad) nos nos da grade, k por triangulo do Darcy, u nos nos da grade."""
        t0 = time.perf_counter()
        xi = sortear_xi(semente, i, self.kl.M)
        c = self.lognormal["c"](self.gradeKL_G.Z(xi[0]))
        phi = np.radians(self.lognormal["phi"](self.gradeKL_G.Z(xi[1])))
        kc = self.lognormal["k"](self.gradeKL_D.Z(xi[2]))       # por celula do Darcy
        k_tri = kc[self.celD_i, self.celD_j]
        t1 = time.perf_counter()
        uD = self.darcy(k_tri)
        t2 = time.perf_counter()
        u = self._u_grade(uD)
        t3 = time.perf_counter()
        return dict(c=c, phi=phi, u=u, uD=uD, k=k_tri, xi=xi, t_campos=t1 - t0, t_darcy=t3 - t1)


if __name__ == "__main__":
    import time
    t0 = time.perf_counter()
    kl = KL2D(-70.0, 0.0, 89.0, 47.0)
    print(f"M = {kl.M}, fracao de variancia = {kl.fracao:.4f}, modos x usados {len(kl.ux)} (max {kl.ix.max()}),"
          f" modos y usados {len(kl.uy)} (max {kl.iy.max()}), {time.perf_counter() - t0:.2f} s")
