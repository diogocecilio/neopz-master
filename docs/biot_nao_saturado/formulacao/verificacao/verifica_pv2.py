"""Confere TPZYCMohrCoulombPV2 (expressoes C++ do repositorio) contra:
 (a) as Tabelas 4-6 e Eqs. (57)-(58) de Lira Cecilio (2026, J Eng Math 157:10);
 (b) a projecao de ponto mais proximo na norma de energia complementar (Eq. 18), por enumeracao de
     conjuntos ativos dos seis planos (Eqs. 44-49) com KKT;
 (c) diferencas finitas das proprias tensoes projetadas (Grad3x3).
Uso: python3 verifica_pv2.py <caminho de TPZYCMohrCoulombPV2.cpp>
"""
import itertools
import math
import re
import sys

import numpy as np

src = open(sys.argv[1]).read()


def body(fname):
    i = src.index("TPZYCMohrCoulombPV2::" + fname + "(")
    j = src.index("{", i)
    depth, k = 0, j
    while True:
        if src[k] == "{":
            depth += 1
        elif src[k] == "}":
            depth -= 1
            if depth == 0:
                return src[j + 1:k]
        k += 1


def split_top(s):
    out, depth, cur = [], 0, ""
    for ch in s:
        if ch in "([{":
            depth += 1
        elif ch in ")]}":
            depth -= 1
        if ch == "," and depth == 0:
            out.append(cur)
            cur = ""
        else:
            cur += ch
    out.append(cur)
    return [o.strip() for o in out]


def code_exprs(fname):
    b = body(fname)
    m = re.search(r"sigpr\s*=\s*\{(.*?)\};", b, re.S)
    sig = split_top(m.group(1))
    grad = {}
    for mm in re.finditer(r"Grad3x3\((\d),(\d)\)\s*=\s*(.*?);", b, re.S):
        grad[(int(mm.group(1)), int(mm.group(2)))] = mm.group(3)
    dl = {}
    for mm in re.finditer(r"dlambda\[(\d)\]\s*=\s*(.*?);", b, re.S):
        dl[int(mm.group(1))] = mm.group(2)
    return sig, grad, dl


ENV = {"pow": math.pow, "sin": math.sin, "cos": math.cos}


def ev(expr, G, K, c, phi, s):
    env = dict(ENV, G=G, K=K, c=c, phi=phi, sinphi=math.sin(phi), cosphi=math.cos(phi),
               s1=s[0], s2=s[1], s3=s[2])
    return eval(" ".join(expr.split()), {"__builtins__": {}}, env)


CODE = {k: code_exprs(f) for k, f in [("main", "ComputeLambdaSigmaMainPlane"),
                                      ("right", "ComputeLambdaSigmaRigth"),
                                      ("left", "ComputeLambdaSigmaLeft")]}


def code_proj(kind, G, K, c, phi, s):
    sig, grad, dl = CODE[kind]
    sp = np.array([ev(e, G, K, c, phi, s) for e in sig])
    J = np.zeros((3, 3))
    for (i, j), e in grad.items():
        J[i, j] = ev(e, G, K, c, phi, s)
    lam = np.array([ev(dl[k], G, K, c, phi, s) for k in sorted(dl)])
    return sp, J, lam


# ---------------- article, Tables 4-6 (transcribed from the PDF), Eqs. (57)-(58)
def article_A(kind, G, K, c, phi):
    sn, cs = math.sin(phi), math.cos(phi)
    if kind == "main":
        D = 2. * (3. * G + (G + 3. * K) * sn ** 2)
        A = [[-3. * G * (-1. + sn), (G + 3. * K) * (-1. + sn) * sn, 0., 2. * c * cs * (3. * K * sn + G * (3. + sn))],
             [-6. * K * sn ** 2 + G * (-3. + sn ** 2), -3. * K * sn * (1. + sn) - G * (3. + (sn - 2.) * sn),
              6. * G + 2. * (G + 3. * K) * sn ** 2, (-4. * c * G + 6. * c * K) * cs * sn],
             [3. * G * (1. + sn), -(G + 3. * K) * sn * (1. + sn), 0., 2. * c * cs * (G * (sn - 3.) + 3. * K * sn)]]
    elif kind == "right":
        D = 9. * G + 6. * G * sn + G * sn ** 2 + 12. * K * sn ** 2
        a = -6. * K * sn + G * (3. + sn)
        b = 12. * K * sn + G * (3. + sn)
        d = 3. * K * sn + G * (3. + sn)
        r23 = [0.5 * (1. + sn) * a, -0.5 * (1. + sn) * b, (1. + sn) * d, -2. * c * cs * a]
        A = [[-0.5 * (-1. + sn) * a, 0.5 * (-1. + sn) * b, -(-1. + sn) * d, 4. * c * cs * d], r23, list(r23)]
    else:
        D = 9. * G - 6. * G * sn + G * sn ** 2 + 12. * K * sn ** 2
        a = G * (-3. + sn) - 6. * K * sn
        d = G * (-3. + sn) + 3. * K * sn
        r12 = [0.5 * (-1. + sn) * a, -0.5 * (-1. + sn) * a, (-1. + sn) * d, 2. * c * cs * (-G * (-3. + sn) + 6. * K * sn)]
        A = [r12, list(r12), [-0.5 * (1. + sn) * a, 0.5 * (1. + sn) * a, -(1. + sn) * d, 4. * c * cs * d]]
    return np.array(A), D


def article_proj(kind, G, K, c, phi, s):
    A, D = article_A(kind, G, K, c, phi)
    s1, s2, s3 = s
    sp = np.array([((s1 + s3) * A[i, 0] + (s1 - s3) * A[i, 1] + (s1 + s2) * A[i, 2] + A[i, 3]) / D for i in range(3)])
    J = np.array([[A[i, 0] + A[i, 1] + A[i, 2], A[i, 2], A[i, 0] - A[i, 1]] for i in range(3)]) / D
    return sp, J


# ---------------- closest-point projection in the energy norm, Eqs. (15)-(18), (44)-(49)
def planes(phi):
    sn = math.sin(phi)
    # Phi = n . sigma - 2 c cos(phi); pairs (i, j): (sigma_i - sigma_j) + (sigma_i + sigma_j) sin(phi)
    pairs = [(0, 2), (1, 2), (1, 0), (2, 0), (2, 1), (0, 1)]  # Phi1..Phi6
    N = []
    for i, j in pairs:
        n = np.zeros(3)
        n[i] += 1. + sn
        n[j] += -1. + sn
        N.append(n)
    return np.array(N)


def Dhw(G, K):
    a, b = K + 4. * G / 3., K - 2. * G / 3.
    return np.array([[a, b, b], [b, a, b], [b, b, a]])


def cpp(G, K, c, phi, s, tol=1e-9):
    """KKT by enumeration of active sets (convex polyhedron): returns projection and active set"""
    N = planes(phi)
    r = 2. * c * math.cos(phi)
    D = Dhw(G, K)
    s = np.asarray(s, float)
    if np.all(N @ s - r <= tol * (1 + abs(r))):
        return s, ()
    scale = 1. + np.abs(s).max()
    for k in (1, 2):
        for act in itertools.combinations(range(6), k):
            Na = N[list(act)]
            M = Na @ D @ Na.T
            if abs(np.linalg.det(M)) < 1e-12 * np.abs(M).max() ** k:
                continue
            lam = np.linalg.solve(M, Na @ s - r)
            sp = s - D @ Na.T @ lam
            if np.all(lam >= -tol * scale) and np.all(N @ sp - r <= tol * scale):
                return sp, act
    apex = np.full(3, c / math.tan(phi))
    return apex, ("apex",)


def code_select(G, K, c, phi, s):
    """ProjectSigma of the code (lines 265-316), with the code expressions"""
    sn = math.sin(phi)
    s1, s2, s3 = s
    if (s1 - s3) + (s1 + s3) * sn - 2. * c * math.cos(phi) <= 0:
        return np.array(s), "elastic"

    def ordered(sp):
        tol = 1e-10  # IsZero
        return (sp[0] > sp[1] or abs(sp[0] - sp[1]) < tol) and (sp[1] > sp[2] or abs(sp[1] - sp[2]) < tol)

    sp, _, _ = code_proj("main", G, K, c, phi, s)
    if ordered(sp):
        return sp, "main"
    if (1 - sn) * s1 - 2 * s2 + (1 + sn) * s3 > 0:
        sp, _, _ = code_proj("right", G, K, c, phi, s)
        if ordered(sp):
            return sp, "right"
    else:
        sp, _, _ = code_proj("left", G, K, c, phi, s)
        if ordered(sp):
            return sp, "left"
    return np.full(3, c / math.tan(phi)), "apex"


rng = np.random.default_rng(2026)
worst = {"art_sig": 0., "art_J": 0., "fd_J": 0., "cpp": 0., "lam": 0.}
count = {}
for trial in range(20000):
    E = 10 ** rng.uniform(3, 7)
    nu = rng.uniform(0.0, 0.49)
    G, K = E / (2 * (1 + nu)), E / (3 * (1 - 2 * nu))
    phi = math.radians(rng.uniform(5, 45))
    c = 10 ** rng.uniform(0, 3)
    s = np.sort(rng.normal(0., 1., 3) * c * rng.uniform(0.5, 20) + rng.normal(0, 1) * c * 5)[::-1]
    sref = max(np.abs(s).max(), c)
    # (a) article tables vs code for the three candidates (formulas, regardless of admissibility)
    for kind in ("main", "right", "left"):
        spc, Jc, lam = code_proj(kind, G, K, c, phi, s)
        spa, Ja = article_proj(kind, G, K, c, phi, s)
        worst["art_sig"] = max(worst["art_sig"], np.abs(spc - spa).max() / sref)
        worst["art_J"] = max(worst["art_J"], np.abs(Jc - Ja).max())
        # (c) Grad3x3 vs central differences of the code stresses
        h = 1e-6 * sref
        Jfd = np.zeros((3, 3))
        for j in range(3):
            sp_, sm_ = s.copy(), s.copy()
            sp_[j] += h
            sm_[j] -= h
            Jfd[:, j] = (code_proj(kind, G, K, c, phi, sp_)[0] - code_proj(kind, G, K, c, phi, sm_)[0]) / (2 * h)
        worst["fd_J"] = max(worst["fd_J"], np.abs(Jfd - Jc).max())
    # (b) selection of the code vs the closest point in the energy norm
    spc, kind = code_select(G, K, c, phi, s)
    spx, act = cpp(G, K, c, phi, s)
    count[kind] = count.get(kind, 0) + 1
    err = np.abs(spc - spx).max() / sref
    if err > worst["cpp"]:
        worst["cpp"] = err
        worst["cpp_case"] = (kind, act, s.tolist(), math.degrees(phi), c, G, K)
    # multipliers of the code (main plane / corners) vs the CPP multipliers
    if kind in ("main", "right", "left"):
        _, _, lam = code_proj(kind, G, K, c, phi, s)
        N = planes(phi)
        Dm = Dhw(G, K)
        # the code multiplier of the main plane multiplies the gradient of Phi1 (Delta sigma = -lam D n1)
        if kind == "main":
            rec = s - Dm @ N[0] * lam[0]
            worst["lam"] = max(worst["lam"], np.abs(rec - spc).max() / sref)

print("casos por regiao (codigo):", count)
for k, v in worst.items():
    print(k, v)
