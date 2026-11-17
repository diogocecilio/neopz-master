#!/usr/bin/env python3
"""diagnostico_ai.py: a AI pelo ponto de projeto do FORM tem vies? O teste com o proprio direto, nas bases analitica (KL do
retangulo) e de Nystrom (KL do talude), sem nenhuma analise nova de elementos finitos.

1. Os xi de cada realizacao saem de uma copia do gerador do RandomField (xoshiro256** semeado por splitmix64, Box-Muller;
   uma semente por realizacao e campo). A copia confere com o log_peso de cada AI (ate 1e-6) e, junto com a KL analitica
   abaixo, com o campo c', phi' que o Cho2010 grava (5e-6, a precisao do CSV).
2. Modo de ruptura fora da AI: as falhas do direto projetadas na direcao do ponto de projeto, t = alfa . xi, e a variancia
   verdadeira da AI com deslocamento mu, pela identidade E_h[I w^2] = E_phi[I w], com as falhas do direto. Um modo que a
   AI nao amostra aparece como falhas com t pequeno e peso w grande, e o CoV verdadeiro passa muito do estimado pela AI.
3. Variavel de controle (o sorteio de cada amostra): h(xi) = Phi((-g(xi) - c)/s), com g = a + b . xi a regressao de ln FS
   nos 2K primeiros xi (K = 200 por campo) e (c, s) o probit das falhas, tudo ajustado no direto. A media exata e
   E_phi[h] = Phi((-a - c)/sqrt(s^2 + |b|^2)). Na AI, P_cv = media(I w) - beta (media(h w) - E[h]); no proprio direto,
   com validacao cruzada em 5 dobras (h ajustado sem a dobra em que e aplicado). O controle nao muda o valor esperado:
   tira a parte do acaso que h explica.
4. Base de Nystrom: o campo de cada realizacao do direto (analitica, 2000 termos) nos 2400 nos do talude, projetado nos
   modos de Nystrom: xi_N = Lambda^(-1/2) F^T W Z_a (T T^T ~ I: a mesma lei N(0, I)). Os itens 2 e 3 valem entao para as
   AI de Nystrom (o FS das falhas e o do campo analitico, uma aproximacao: o campo de Nystrom e a projecao dele).
Uso: python3 -I diagnostico_ai.py   (grava <runs>/cho2010/diagnostico/diagnostico.json; caches .npy na mesma pasta)
     No relatorio: runpy.run_path(...)["roda"]() devolve o dicionario, com os vetores das figuras.
"""
import csv
import glob
import json
import math
import os

import numpy as np
from scipy.optimize import minimize
from scipy.stats import norm

RUNS = os.path.expanduser("~/projects/runs")
E1, E3 = f"{RUNS}/fermi_etapa1/etapa1", f"{RUNS}/fermi_etapa3/etapa3"
G = f"{RUNS}/cho2010/grossa/campanha"
BASE_NYS = f"{RUNS}/cho2010/nystrom/base_nystrom.txt"
CACHE = f"{RUNS}/cho2010/diagnostico"
K = 200  # termos por campo na regressao e nas projecoes
SEMENTES_NOVAS = range(11, 19)  # fase 4 da etapa3


# ---------------------------------------------------------------- gerador (copia do klsep::Normal e do RandomField::Sample)
def _u(x):
    return np.uint64(x)


def _rot(x, k):
    return (x << _u(k)) | (x >> _u(64 - k))


def normais(sementes, n):
    """os n primeiros normais de cada fluxo (sementes: uint64)"""
    with np.errstate(over="ignore"):
        x = sementes.astype(np.uint64).copy()
        s = []
        for _ in range(4):
            x = x + _u(0x9E3779B97F4A7C15)
            z = (x ^ (x >> _u(30))) * _u(0xBF58476D1CE4E5B9)
            z = (z ^ (z >> _u(27))) * _u(0x94D049BB133111EB)
            s.append(z ^ (z >> _u(31)))
        s0, s1, s2, s3 = s
        out = np.empty((len(sementes), n))
        for j in range(0, n, 2):
            u = []
            for _ in range(2):
                r = _rot(s1 * _u(5), 7) * _u(9)
                t = s1 << _u(17)
                s2 = s2 ^ s0
                s3 = s3 ^ s1
                s1 = s1 ^ s2
                s0 = s0 ^ s3
                s2 = s2 ^ t
                s3 = _rot(s3, 45)
                u.append(((r >> _u(11)).astype(np.float64) + 0.5) * (1.0 / 9007199254740992.0))
            raio, ang = np.sqrt(-2.0 * np.log(u[0])), 2.0 * np.pi * u[1]
            out[:, j] = raio * np.cos(ang)
            if j + 1 < n:
                out[:, j + 1] = raio * np.sin(ang)
    return out


def xi(semente, reais, n, ncampos=2):
    """os n primeiros xi de cada campo, sem deslocamento: (len(reais), ncampos, n)"""
    base = (semente * 0x9E3779B97F4A7C15) & ((1 << 64) - 1)
    res = np.empty((len(reais), ncampos, n))
    for c in range(ncampos):
        sem = np.array([(base + ncampos * int(r) + c) & ((1 << 64) - 1) for r in reais], dtype=np.uint64)
        res[:, c, :] = normais(sem, n)
    return res


# ---------------------------------------------------------------- KL analitica (copia do klsep::KL2D) e base de Nystrom
def _bissecao(g, lo, hi):
    glo = g(lo)
    for _ in range(200):
        m = 0.5 * (lo + hi)
        gm = g(m)
        if (gm > 0.) == (glo > 0.):
            lo, glo = m, gm
        else:
            hi = m
        if hi - lo < 1.e-15 * max(1., hi):
            break
    return 0.5 * (lo + hi)


def _modos1d(comp, L, n):
    a, c, eps = 0.5 * comp, 1.0 / L, 1e-12
    v, k = [], 1
    while len(v) < n:
        w = _bissecao(lambda x: x * math.tan(x * a) - c, ((k - 1) * math.pi + eps) / a, ((k - 0.5) * math.pi - eps) / a)
        v.append((2 * c / (w * w + c * c), w, 1))
        if len(v) >= n:
            break
        w = _bissecao(lambda x: x + c * math.tan(x * a), ((k - 0.5) * math.pi + eps) / a, (k * math.pi - eps) / a)
        v.append((2 * c / (w * w + c * c), w, 0))
        k += 1
    v.sort(key=lambda p: -p[0])
    om, par = np.array([p[1] for p in v]), np.array([p[2] for p in v])
    s2 = np.sin(2 * om * a) / (2 * om)
    return dict(lam=np.array([p[0] for p in v]), om=om, par=par, norma=1 / np.sqrt(np.where(par == 1, a + s2, a - s2)))


class KL2D:
    """o retangulo [0, 30] x [0, 15] do Cho2010, l_h = 20 m, l_v = 2 m, M = 2000"""

    def __init__(self, Dx=30., Dy=15., Lx=20., Ly=2., M=2000):
        self.Dx, self.Dy = Dx, Dy
        self.mx, self.my = _modos1d(Dx, Lx, M), _modos1d(Dy, Ly, M)
        cand = sorted(((self.mx["lam"][i] * self.my["lam"][j], i, j) for i in range(M) for j in range(M)
                       if (i + 1) * (j + 1) <= M), key=lambda t: -t[0])[:M]
        self.ix, self.iy = np.array([t[1] for t in cand]), np.array([t[2] for t in cand])
        self.raiz = np.sqrt(np.array([t[0] for t in cand]))
        self.fracao = (self.raiz ** 2).sum() / (Dx * Dy)

    def psi(self, x, y):
        """sqrt(lambda_k) f_k nos pontos: (npontos, M)"""
        def f1(m, i, s):
            return m["norma"][i] * (np.cos(m["om"][i] * s) if m["par"][i] else np.sin(m["om"][i] * s))
        sx, sy = np.asarray(x) - 0.5 * self.Dx, np.asarray(y) - 0.5 * self.Dy
        fx = {i: f1(self.mx, i, sx) for i in np.unique(self.ix)}
        gy = {j: f1(self.my, j, sy) for j in np.unique(self.iy)}
        return np.stack([self.raiz[k] * fx[self.ix[k]] * gy[self.iy[k]] for k in range(len(self.raiz))], axis=1)


def le_nystrom(arq):
    with open(arq) as f:
        n, M = map(int, f.readline().split())
        nos = np.array([list(map(float, f.readline().split())) for _ in range(n)])
        lam = np.array(list(map(float, f.readline().split())))
        F = np.loadtxt(f, max_rows=n).reshape(n, M)
    return dict(x=nos[:, 0], y=nos[:, 1], w=nos[:, 2], lam=lam, F=F)


# ---------------------------------------------------------------- leitura das rodadas
def realizacoes(pasta):
    out = {}
    for arq in glob.glob(os.path.join(pasta, "mc_[0-9]*_[0-9]*.csv")):
        with open(arq) as f:
            for l in csv.DictReader(f):
                if int(l["real"]) >= 0:
                    out[int(l["real"])] = l
    return out


def parametros(pasta):
    p = {}
    for l in open(os.path.join(pasta, "mc_parametros.txt")):
        if "=" in l:
            k, v = l.strip().split("=", 1)
            p[k] = v
    return p


def desloc(pasta):
    """mu (2, K) do mc_desloc_*.txt de uma AI com o deslocamento dado (FORM)"""
    arqs = sorted(glob.glob(os.path.join(pasta, "mc_desloc_*.txt")))
    c, p, dentro = [], [], False
    for l in open(arqs[0]):
        if l.startswith("termo,"):
            dentro = True
        elif dentro and l.strip():
            _, a, b = l.split(",")
            c.append(float(a)), p.append(float(b))
    return np.array([c, p])


# ---------------------------------------------------------------- o substituto e o controle
def ajusta(X, fs):
    """ln FS = a + b . X (minimos quadrados) e o probit P(falha | X) = Phi((-g - c)/s); devolve (a, b, c, s, R2)"""
    A = np.hstack([np.ones((len(X), 1)), X])
    coef, *_ = np.linalg.lstsq(A, np.log(fs), rcond=None)
    g = A @ coef
    I = fs < 1

    def nll(p):
        q = np.clip(norm.cdf((-g - p[0]) / np.exp(p[1])), 1e-300, 1 - 1e-16)
        return -(I * np.log(q) + (~I) * np.log(1 - q)).sum()

    r = minimize(nll, [0., math.log(0.03)], method="Nelder-Mead", options=dict(xatol=1e-8, fatol=1e-8, maxiter=4000))
    res = np.log(fs) - g
    return coef[0], coef[1:], r.x[0], math.exp(r.x[1]), 1 - res.var() / np.log(fs).var()


def h_de(sub, X):
    a, b, c, s, _ = sub
    return norm.cdf((-(a + X @ b) - c) / s)


def media_exata(sub):
    a, b, c, s, _ = sub
    return norm.cdf((-a - c) / math.sqrt(s * s + b @ b))


def com_controle(y, x, Eh):
    """P_cv e o seu desvio, rho e beta (y = I w, x = h w)"""
    beta = np.cov(y, x)[0, 1] / x.var(ddof=1)
    p = y.mean() - beta * (x.mean() - Eh)
    return p, (y - beta * x).std(ddof=1) / math.sqrt(len(y)), float(np.corrcoef(y, x)[0, 1]), beta


def direto_com_controle(X, fs, dobras=5, semente=12345):
    I = (fs < 1).astype(float)
    dob = np.random.default_rng(semente).permutation(len(I)) % dobras
    h, Eh = np.empty(len(I)), np.empty(len(I))
    for j in range(dobras):
        sub = ajusta(X[dob != j], fs[dob != j])
        h[dob == j], Eh[dob == j] = h_de(sub, X[dob == j]), media_exata(sub)
    beta = np.cov(I, h)[0, 1] / h.var(ddof=1)
    p = I.mean() - beta * (h - Eh).mean()
    sd = (I - beta * (h - Eh)).std(ddof=1) / math.sqrt(len(I))
    return dict(p=float(p), sd=float(sd), rho=float(np.corrcoef(I, h)[0, 1]))


# ---------------------------------------------------------------- o diagnostico
def casos():
    """(r, base, rotulo, pasta) das AI pelo ponto de projeto"""
    v = [("r05", "kl", "FORM $M = 10$, dif. centrais", f"{E1}/r05_ai_form"),
         ("r05", "kl", "FORM $M = 20$, dif. centrais", f"{E1}/r05_ai_form20"),
         ("r05", "kl", "FORM $M = 20$, dif. centrais", f"{E3}/r05_ai_form20_s4"),
         ("r05", "kl", "FORM $M = 20$, dif. centrais", f"{E3}/r05_ai_form20_s5"),
         ("r05", "kl", "FORM modal $M = 100$", f"{G}/ai_form100_modal"),
         ("r05", "kl", "FORM modal $M = 200$", f"{G}/ai_form200_modal"),
         ("r05", "kl", "FORM modal $M = 200$, campo de 200 termos", f"{G}/ai_form200_trunc"),
         ("r05", "kl", "FORM modal $M = 100$", f"{E3}/ai_form100_modal_s6"),
         ("r05", "kl", "FORM modal $M = 200$", f"{E3}/ai_form200_modal_s7"),
         ("r05", "kl", "direção da Tri(0), $M = 100$", f"{G}/ai_form100_tri0")]
    v += [("r05", "kl", "FORM modal $M = 100$", f"{E3}/ai_form100_modal_s{s}") for s in SEMENTES_NOVAS]
    v += [("r05", "nys", "FORM modal $M = 20$", f"{G}/ai_nys_form20"),
          ("r05", "nys", "FORM modal $M = 200$", f"{G}/ai_nys_form200"),
          ("r05", "nys", "FORM modal $M = 20$", f"{E3}/ai_nys_form20_s8"),
          ("r05", "nys", "FORM modal $M = 200$", f"{E3}/ai_nys_form200_s9")]
    v += [("r05", "nys", "FORM modal $M = 200$", f"{E3}/ai_nys_form200_s{s}") for s in SEMENTES_NOVAS]
    v += [("r07", "kl", "FORM $M = 10$, dif. centrais", f"{E3}/r07_ai_form10"),
          ("r07", "kl", "FORM $M = 20$, dif. centrais", f"{E3}/r07_ai_form20"),
          ("r07", "kl", "FORM $M = 20$, dif. centrais", f"{E3}/r07_ai_form20_s4"),
          ("r07", "kl", "FORM $M = 20$, dif. centrais", f"{E3}/r07_ai_form20_s5")]
    return [c for c in v if glob.glob(f"{c[3]}/mc_[0-9]*_[0-9]*.csv")]


def roda():
    os.makedirs(CACHE, exist_ok=True)
    diretos = {"r05": f"{E1}/r05_direto", "r07": f"{E3}/r07_direto"}
    D = {}
    for r, pasta in diretos.items():
        R = realizacoes(pasta)
        reais = np.array(sorted(k for k, l in R.items() if l["situacao"] == "ok"))
        D[r] = dict(reais=reais, fs=np.array([float(R[k]["FS_min"]) for k in reais]),
                    semente=int(parametros(pasta)["SEMENTE"]), n_total=len(R))
    # os xi do direto (as duas campanhas usam a semente 2010: os mesmos xi), analiticos (K) e de Nystrom
    todos = np.union1d(D["r05"]["reais"], D["r07"]["reais"])
    assert D["r05"]["semente"] == D["r07"]["semente"]
    arq = f"{CACHE}/xi_direto_{D['r05']['semente']}_{len(todos)}_{todos[-1]}.npz"
    if os.path.exists(arq):
        z = np.load(arq)
        XA, XN = z["XA"], z["XN"]
    else:
        kl, B = KL2D(), le_nystrom(BASE_NYS)
        T = (B["F"] * B["w"][:, None]).T @ kl.psi(B["x"], B["y"]) / np.sqrt(B["lam"])[:, None]  # xi_N = T xi_a
        XA, XN = np.empty((len(todos), 2, K)), np.empty((len(todos), 2, T.shape[0]))
        for i in range(0, len(todos), 2000):
            za = xi(D["r05"]["semente"], todos[i:i + 2000], 2000)
            XA[i:i + 2000], XN[i:i + 2000] = za[:, :, :K], za @ T.T
        TT = T @ T.T
        np.savez(arq, XA=XA, XN=XN, TTdiag=np.diag(TT), TTfora=np.abs(TT - np.diag(np.diag(TT)))[:50, :50].max())
    pos = {k: i for i, k in enumerate(todos)}
    z = np.load(arq)
    out = {"diretos": {}, "ai": [], "K": K,
           "TT": dict(d10=float(z["TTdiag"][:10].min()), d50=float(z["TTdiag"][49]), d100=float(z["TTdiag"][99]),
                      d200=float(z["TTdiag"][199]), fora=float(z["TTfora"]))}
    sub = {}
    for r in D:
        idx = np.array([pos[k] for k in D[r]["reais"]])
        D[r]["XA"], D[r]["XN"] = XA[idx].reshape(len(idx), -1), XN[idx].reshape(len(idx), -1)
        fs = D[r]["fs"]
        I = fs < 1
        P = I.mean()
        sub[(r, "kl")] = ajusta(D[r]["XA"], fs)
        sub[(r, "nys")] = ajusta(D[r]["XN"], fs)
        cv = direto_com_controle(D[r]["XA"], fs)
        out["diretos"][r] = dict(n=len(fs), falhas=int(I.sum()), p=float(P), sd=float(math.sqrt(P * (1 - P) / len(fs))),
                                 cv=cv, n_equiv=float(len(fs) * (P * (1 - P) / len(fs)) / cv["sd"] ** 2),
                                 r2_kl=float(sub[(r, "kl")][4]), r2_nys=float(sub[(r, "nys")][4]))
    for r, base, rot, pasta in casos():
        par = parametros(pasta)
        sem = int(par["SEMENTE"])
        R = realizacoes(pasta)
        reais = np.array(sorted(k for k, l in R.items() if l["situacao"] == "ok"))
        mu = desloc(pasta)
        k = mu.shape[1]
        z = xi(sem, reais, max(K, k))
        z[:, :, :k] += mu
        lw = np.array([float(R[j]["log_peso"]) for j in reais])
        erro_peso = float(np.abs(lw - (-(z[:, :, :k] * mu).sum(axis=(1, 2)) + 0.5 * (mu ** 2).sum())).max())
        w = np.exp(lw)
        y = (np.array([float(R[j]["FS_min"]) for j in reais]) < 1) * w
        n = len(y)
        s = sub[(r, base)]
        x = h_de(s, z[:, :, :K].reshape(n, -1)) * w
        Eh = media_exata(s)
        pcv, sdcv, rho, _ = com_controle(y, x, Eh)
        # pelo direto: as falhas projetadas na direcao de mu e a variancia verdadeira desta AI
        X = D[r]["XA" if base == "kl" else "XN"]
        fs = D[r]["fs"]
        I = fs < 1
        m = math.sqrt((mu ** 2).sum())
        zm = X.reshape(len(X), 2, -1)[:, :, :k].__mul__(mu).sum(axis=(1, 2))
        t = zm / m
        wd = np.exp(-zm + 0.5 * m * m)
        P = I.mean()
        m2 = float((I * wd).mean())
        bordas = np.arange(-4., 5.01, 0.25)
        hist = dict(bordas=bordas.tolist(), todas=np.histogram(t, bordas)[0].tolist(),
                    falhas=np.histogram(t[I], bordas)[0].tolist())
        out["ai"].append(dict(r=r, base=base, rotulo=rot, pasta=pasta, semente=sem, n=n, p=float(y.mean()),
                              sd=float(y.std(ddof=1) / math.sqrt(n)), p_cv=float(pcv), sd_cv=float(sdcv), rho=rho,
                              sorte=float(x.mean() / Eh), erro_peso=erro_peso, mu=m, M_desloc=k,
                              cov_verdadeiro=float(math.sqrt(max(m2 - P * P, 0) / n) / P),
                              t_min=float(t[I].min()), t_q=[float(q) for q in np.percentile(t[I], [5, 50, 95])],
                              hits_esperados=float(n * (I / wd).mean()), hist=hist,
                              abaixo_15=int((t[I] < 1.5).sum())))
    with open(f"{CACHE}/diagnostico.json", "w") as f:
        json.dump(out, f, indent=1)
    return out


def ponto_de_projeto(ustar_kl=f"{RUNS}/cho2010/grossa/ustar_M100_modal.txt",
                     ustar_nys=f"{RUNS}/cho2010/nystrom/ustar_nys_M200.txt", rho=-0.5, passo=0.25):
    """os campos c' e phi' no ponto de projeto das duas bases (grade no talude) e o angulo entre os dois pontos, com o da
    analitica levado para as coordenadas de Nystrom (xi_N = T xi_a)"""
    kl, B = KL2D(), le_nystrom(BASE_NYS)
    ua, un = np.loadtxt(ustar_kl), np.loadtxt(ustar_nys)
    ka, kn = len(ua) // 2, len(un) // 2
    xs, ys = np.arange(passo / 2, 30, passo), np.arange(passo / 2, 15, passo)
    X, Y = np.meshgrid(xs, ys)
    dentro = (Y <= np.where(X <= 10, 15, np.where(X >= 20, 5, 15 - (X - 10))))  # o talude espelhado (crista a esquerda)
    px, py = X[dentro], Y[dentro]
    Pa = kl.psi(px, py)[:, :ka]
    Cg = np.exp(-np.abs(px[:, None] - B["x"][None, :]) / 20. - np.abs(py[:, None] - B["y"][None, :]) / 2.)
    Pn = Cg @ (B["F"][:, :kn] * B["w"][:, None]) / np.sqrt(B["lam"][:kn])
    sc, sf = math.sqrt(math.log(1 + 0.3 ** 2)), math.sqrt(math.log(1 + 0.2 ** 2))
    mc, mf = math.log(10 / math.sqrt(1 + 0.3 ** 2)), math.log(30 / math.sqrt(1 + 0.2 ** 2))
    out = {}
    for nome, P, u, k in (("kl", Pa, ua, ka), ("nys", Pn, un, kn)):
        Zc, Zf = P @ u[:k], P @ u[k:]
        cc, ff = np.full(X.shape, np.nan), np.full(X.shape, np.nan)
        cc[dentro] = np.exp(mc + sc * Zc)
        ff[dentro] = np.exp(mf + sf * (rho * Zc + math.sqrt(1 - rho * rho) * Zf))
        out[nome] = dict(c=cc, phi=ff, beta=float(np.linalg.norm(u)))
    T = (B["F"] * B["w"][:, None]).T @ kl.psi(B["x"], B["y"])[:, :ka] / np.sqrt(B["lam"])[:, None]
    ua_n = np.concatenate([T[:kn] @ ua[:ka], T[:kn] @ ua[ka:]])
    out["angulo"] = float(np.degrees(np.arccos(ua_n @ un / (np.linalg.norm(ua_n) * np.linalg.norm(un)))))
    out["beta_kl_em_nys"] = float(np.linalg.norm(ua_n))
    out["x"], out["y"] = X, Y
    return out


if __name__ == "__main__":
    pp = ponto_de_projeto()
    print(f"ponto de projeto: beta analitica {pp['kl']['beta']:.3f}, Nystrom {pp['nys']['beta']:.3f}; o analitico nas coordenadas "
          f"de Nystrom tem |u| = {pp['beta_kl_em_nys']:.3f} e fica a {pp['angulo']:.1f} graus do de Nystrom; phi' minimo "
          f"{np.nanmin(pp['kl']['phi']):.1f} e {np.nanmin(pp['nys']['phi']):.1f} graus")
    o = roda()
    for r, d in o["diretos"].items():
        print(f"direto {r}: {d['n']} realizacoes, {d['falhas']} falhas, P_f = {100*d['p']:.3f}% +- {100*d['sd']:.3f}; com "
              f"controle {100*d['cv']['p']:.3f}% +- {100*d['cv']['sd']:.3f} (rho {d['cv']['rho']:.3f}, = {d['n_equiv']:.0f} "
              f"do direto puro); R2 ln FS: analitica {d['r2_kl']:.3f}, Nystrom {d['r2_nys']:.3f}")
    print(f"{'r':4s} {'base':4s} {'AI':42s} {'sem':>3s} {'N':>5s} {'P_f':>6s} {'CoV':>5s} {'CoV v.':>6s} {'t_min':>5s} | "
          f"{'sorte':>6s} {'rho':>5s} {'P_cv':>6s} {'CoV':>5s}")
    for a in o["ai"]:
        print(f"{a['r']:4s} {a['base']:4s} {a['rotulo'][:42]:42s} {a['semente']:3d} {a['n']:5d} {100*a['p']:6.3f} "
              f"{100*a['sd']/a['p']:5.1f} {100*a['cov_verdadeiro']:6.1f} {a['t_min']:5.2f} | {100*a['sorte']:5.1f}% "
              f"{a['rho']:5.2f} {100*a['p_cv']:6.3f} {100*a['sd_cv']/a['p_cv']:5.1f}   (peso: erro {a['erro_peso']:.1e})")
