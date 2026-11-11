#!/usr/bin/env python3
"""nystrom.py: a KL "exata" do talude (referencia independente do NeoPZ) pelo metodo de Nystrom.

Grade uniforme de passo h = 5/n no talude (nos (ih, jh)); as dobras do nucleo (x = x', y = y') ficam nas arestas das
celulas; trapezio tensorial nos quadrados e regra dos vertices (h^2/6 por vertice) nos meios-quadrados cortados pela
face x + y = 20. O problema discreto e W^1/2 K W^1/2 y = lambda y (K_ij = C(x_i, x_j), W = diag(pesos)); os modos nos
nos sao phi_k(x_j) = y_kj / sqrt(w_j) e, num ponto qualquer, a formula de Nystrom phi_k(x) = sum_j w_j C(x, x_j)
phi_k(x_j) / lambda_k. O produto por K nao monta a matriz: o nucleo separavel na grade do retangulo envolvente e
Ky U Kx^T (mesmo esquema de nystrom_talude_mf.py da verificacao independente).

Uso: nystrom.py [n] [k]   (padrao n = 160, k = 210; com n = 160 sao 193 mil nos, cerca de 1 GB e 1 a 3 min)
Grava em ../../SlopeRandomFieldRefina/relatorio/dados:
  kl_nystrom_modos.csv   os 9 primeiros modos nos nos de uma grade de 0,25 m no talude
  kl_nystrom_s2.csv      sigma_M^2 (M = 50 e 200) nos nos de uma grade de 0,125 m
  kl_nystrom_resumo.csv  por M: minimo e maximo de sigma_M^2 em todos os nos, 1 - sum lambda / area e os erros da
                         covariancia truncada nos 1000 pares de pares.py (os do KLTaludeCeron)
  kl_nystrom_custo.csv   tempo e memoria
"""
import os
import resource
import sys
import time

import numpy as np
from scipy.sparse.linalg import LinearOperator, eigsh

import estilo as E
from pares import pontos

LX, LY, AREA = 20., 2., 187.5


def grade(n):
    """pesos W[j, i] dos nos (i h, j h) da grade do retangulo [0,25] x [0,10] (zero fora do talude)"""
    h = 5. / n
    nx, ny = 5 * n + 1, 2 * n + 1
    W = np.zeros((ny, nx))
    j, i = np.meshgrid(np.arange(2 * n), np.arange(5 * n), indexing="ij")
    cheia = (j < n) | (i + j <= 4 * n - 2)
    tri = (j >= n) & (i + j == 4 * n - 1)
    for dj, di in ((0, 0), (0, 1), (1, 0), (1, 1)):
        np.add.at(W, (j[cheia] + dj, i[cheia] + di), h * h / 4)
    for dj, di in ((0, 0), (0, 1), (1, 0)):
        np.add.at(W, (j[tri] + dj, i[tri] + di), h * h / 6)
    return h, W


class Nystrom:
    def __init__(self, n):
        self.n = n
        self.h, self.W = grade(n)
        self.jj, self.ii = np.nonzero(self.W)
        self.w = self.W[self.jj, self.ii]
        self.sw = np.sqrt(self.w)
        ny, nx = self.W.shape
        self.gx = np.arange(nx) * self.h
        self.gy = np.arange(ny) * self.h
        self.Kx = np.exp(-np.abs(self.gx[:, None] - self.gx[None, :]) / LX)
        self.Ky = np.exp(-np.abs(self.gy[:, None] - self.gy[None, :]) / LY)
        self.N = len(self.w)

    def mv(self, v):
        U = np.zeros(self.W.shape)
        U[self.jj, self.ii] = self.sw * np.ravel(v)
        return self.sw * (self.Ky @ U @ self.Kx.T)[self.jj, self.ii]

    def resolve(self, k, tol=1e-10):
        op = LinearOperator((self.N, self.N), matvec=self.mv, dtype=float)
        lam, Y = eigsh(op, k=k, which="LA", tol=tol, ncv=min(self.N - 1, 2 * k + 40))
        o = np.argsort(-lam)
        self.lam, self.Y = lam[o], Y[:, o]
        # sinal: o valor de maior modulo positivo
        for m in range(k):
            a = np.argmax(np.abs(self.Y[:, m]))
            if self.Y[a, m] < 0:
                self.Y[:, m] *= -1
        return self.lam

    def modos_nos(self, m):
        """phi_m nos nos (normalizados: sum w phi^2 = 1)"""
        return self.Y[:, m] / self.sw

    def modos_pontos(self, P, M):
        """phi_k(P), k < M, pela formula de Nystrom: (M, npontos)"""
        GX = np.exp(-np.abs(P[:, 0][:, None] - self.gx[None, :]) / LX)
        GY = np.exp(-np.abs(P[:, 1][:, None] - self.gy[None, :]) / LY)
        PHI = np.empty((M, len(P)))
        U = np.zeros(self.W.shape)
        for m in range(M):
            U[self.jj, self.ii] = self.sw * self.Y[:, m]
            PHI[m] = np.einsum("pj,pj->p", GY @ U, GX) / self.lam[m]
        return PHI

    def sub(self, passo_m):
        """nos da grade de passo passo_m (multiplo de h) dentro do talude: indices em jj/ii"""
        q = max(1, int(round(passo_m / self.h)))
        s = (self.jj % q == 0) & (self.ii % q == 0)
        return np.nonzero(s)[0]


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 160
    k = int(sys.argv[2]) if len(sys.argv) > 2 else 210
    t0 = time.time()
    ny = Nystrom(n)
    print("n = %d, h = %g m, %d nos, area %.10f" % (n, ny.h, ny.N, ny.w.sum()), flush=True)
    lam = ny.resolve(k)
    t_auto = time.time() - t0
    print("autopares: %.1f s; 10 maiores: %s" % (t_auto, " ".join("%.5f" % v for v in lam[:10])), flush=True)
    x, y = ny.ii * ny.h, ny.jj * ny.h
    # modos 1 a 9 na grade de 0,25 m
    s = ny.sub(0.25)
    with open(os.path.join(E.DADOS, "kl_nystrom_modos.csv"), "w") as f:
        f.write("x,y," + ",".join("phi%d" % (m + 1) for m in range(9)) + "\n")
        for p in s:
            f.write("%.4f,%.4f," % (x[p], y[p]) + ",".join("%.5f" % (ny.Y[p, m] / ny.sw[p]) for m in range(9)) + "\n")
    # sigma_M^2 nos nos
    Ms = [M for M in (10, 20, 50, 100, 200) if M <= k]
    s2 = {}
    acc = np.zeros(ny.N)
    for m in range(max(Ms)):
        acc += lam[m] * ny.Y[:, m] ** 2 / ny.w
        if m + 1 in Ms:
            s2[m + 1] = acc.copy()
    s = ny.sub(0.125)
    with open(os.path.join(E.DADOS, "kl_nystrom_s2.csv"), "w") as f:
        f.write("x,y,s2_M50,s2_M200\n")
        for p in s:
            f.write("%.4f,%.4f,%.5f,%.5f\n" % (x[p], y[p], s2[50][p], s2[200][p]))
    # covariancia truncada nos pares do KLTaludeCeron
    P = np.array(pontos())
    npar = len(P) // 4
    t1 = time.time()
    PHI = ny.modos_pontos(P, max(Ms))
    t_interp = time.time() - t1
    G = np.exp(-np.abs(P[0::2, 0] - P[1::2, 0]) / LX - np.abs(P[0::2, 1] - P[1::2, 1]) / LY)
    with open(os.path.join(E.DADOS, "kl_nystrom_resumo.csv"), "w") as f:
        f.write("M,n,erro_mq,s2_min_nos,s2_max_nos,s2_min_pares,s2_max_pares,cov_unif_max,cov_unif_media,cov_prox_max,"
                "cov_prox_media\n")
        for M in Ms:
            C = np.sum(lam[:M, None] * PHI[:M, 0::2] * PHI[:M, 1::2], axis=0)
            e = np.abs(C - G)
            s2p = np.sum(lam[:M, None] * PHI[:M] ** 2, axis=0)
            lin = (M, n, 1 - lam[:M].sum() / AREA, s2[M].min(), s2[M].max(), s2p.min(), s2p.max(), e[:npar].max(),
                   e[:npar].mean(), e[npar:].max(), e[npar:].mean())
            f.write("%d,%d,%.8f,%.6f,%.6f,%.6f,%.6f,%.4e,%.4e,%.4e,%.4e\n" % lin)
            print("M = %d: 1 - sum/area %.5f; sigma2 nos nos %.4f a %.4f; cov. uniformes max %.2e media %.2e, proximos "
                  "max %.2e media %.2e" % (M, lin[2], lin[3], lin[4], lin[7], lin[8], lin[9], lin[10]), flush=True)
    mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.
    with open(os.path.join(E.DADOS, "kl_nystrom_custo.csv"), "w") as f:
        f.write("n,h,nos,k,t_autopares_s,t_interpolacao_s,t_total_s,memoria_MB\n")
        f.write("%d,%g,%d,%d,%.1f,%.1f,%.1f,%.0f\n" % (n, ny.h, ny.N, k, t_auto, t_interp, time.time() - t0, mem))
    np.save(os.path.join(E.SCRATCH, "kl_rel_tal", "nystrom_n%d_lam.npy" % n), lam)
    print("total %.1f s, pico de memoria %.0f MB" % (time.time() - t0, mem))


if __name__ == "__main__":
    main()
