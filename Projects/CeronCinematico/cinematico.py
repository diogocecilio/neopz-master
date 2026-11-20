"""Abordagem cinematica (limite superior) de Vargas Ceron et al. (2025),
IJNAMG 49(11):2459-2491, doi 10.1002/nag.3993.

Convencoes (as do artigo, Figs. 1, 2 e 7):
  * origem O na aresta da crista; x horizontal para fora do talude;
    y VERTICAL PARA BAIXO (profundidade abaixo da crista);
  * crista: y = 0, x <= 0; face: de O ate o pe T = (H/tan(beta), H);
    terreno a frente do pe: y = H, x >= H/tan(beta);
  * o bloco gira em torno de C = (Cx, -Cy) (Cy > 0 = altura de C acima da
    crista); theta e medido a partir da direcao -x, girando para baixo:
        P(theta) = (Cx - r cos(theta), -Cy + r sin(theta)),
        U = omega r e_theta,  e_theta = (sin(theta), cos(theta));
  * u = p - gamma_w y e o excesso de poropressao (p total, compressao > 0);
    as forcas de corpo sao gamma' e_y - grad u (gamma' = peso efetivo).

Fator de estabilidade (Eq. 54): Gamma = P_mr / (P_gamma + P_u), com
    P_mr    = dissipacao na descontinuidade,
    P_gamma = int_bloco gamma' e_y . U dA,
    P_u     = int_bloco (-grad u) . U dA.
Gamma multiplica TODAS as forcas de corpo (peso e percolacao).

Duas implementacoes:
  1. classica (homogenea): mecanismos (I) e (II) em forma fechada
     (Eqs. 42-53; o (II) com o triangulo C-T-B, que o artigo omite);
  2. geral (heterogenea): espiral por trechos com foco C unico, arcos de
     comprimento ~H/100, phi_i = max de phi(x) ao longo do proprio arco
     (o minimo que a normalidade pede), c e phi locais na dissipacao, e
     integracao numerica sobre o bloco (tambem com C atras da reta da face).
"""
from __future__ import annotations

import math
import time
from dataclasses import dataclass, field

import numpy as np
from scipy import optimize
from scipy.spatial import cKDTree
from matplotlib.tri import Triangulation

GAMMA_W = 9.81
INF = float("inf")


# ----------------------------------------------------------------------------
# Geometria
# ----------------------------------------------------------------------------
@dataclass
class Talude:
    """Talude simples. Coordenadas do artigo (y para baixo, O na crista).

    x_esq, x_dir: limites laterais do dominio; y_base: profundidade da base
    (rigida ou fundo do dominio), medida a partir da crista."""
    H: float
    beta: float                 # inclinacao da face (rad)
    x_esq: float = -INF
    x_dir: float = INF
    y_base: float = INF

    @property
    def xT(self) -> float:
        return self.H / math.tan(self.beta)

    def y_sup(self, x):
        """Profundidade da superficie do terreno em x (vetorizado)."""
        x = np.asarray(x, float)
        return np.where(x <= 0.0, 0.0,
                        np.where(x < self.xT, x * math.tan(self.beta), self.H))

    def y_sup_escalar(self, x: float) -> float:
        if x <= 0.0:
            return 0.0
        if x < self.xT:
            return x * math.tan(self.beta)
        return self.H


# ----------------------------------------------------------------------------
# Campos (c, phi, excesso de poropressao u)
# ----------------------------------------------------------------------------
class CampoUniforme:
    """c e phi constantes; u = 0 (talude seco ou sem percolacao).

    gamma_ef: peso especifico usado em P_gamma (gamma' com agua, gamma total
    no talude seco)."""
    uniforme = True
    tem_potencial = True

    def __init__(self, c: float, phi: float, gamma_ef: float):
        self.c0, self.phi0, self.gamma_ef = float(c), float(phi), float(gamma_ef)
        self.phi_ref = self.phi0
        self.tem_u = False

    def c(self, x, y):
        return np.full(np.shape(x), self.c0)

    def phi(self, x, y):
        return np.full(np.shape(x), self.phi0)

    def u(self, x, y):
        return np.zeros(np.shape(x))

    def grad_u(self, x, y):
        z = np.zeros(np.shape(x))
        return z, z

    def phi_max_arco(self, xs, ys, flecha=0.0) -> float:
        return self.phi0


class CampoMalha:
    """Campos nodais P1 numa malha de triangulos (a malha do transiente
    exportada pelo NeoPZ). Coordenadas ja no sistema do artigo.

    * c, phi, u: interpolacao linear por triangulo (os nos sao os vertices);
    * grad u: gradiente constante por triangulo (P1);
    * phi_max_arco: limite superior do maximo de phi(x) P1 ao longo de um
      arco amostrado (localizacao exata pelo trifinder).
    Pontos fora da malha (ar, por arredondamento) usam extrapolacao linear
    do triangulo mais proximo (pelo centroide)."""
    uniforme = False
    tem_potencial = True

    def __init__(self, xy, tri, c, phi, u, gamma_ef):
        xy = np.asarray(xy, float)
        self.tri = np.asarray(tri, int)
        self.T = Triangulation(xy[:, 0], xy[:, 1], self.tri)
        self.finder = self.T.get_trifinder()
        self.xy = xy
        self.cn = np.asarray(c, float)
        self.phin = np.asarray(phi, float)
        self.un = np.asarray(u, float)
        self.gamma_ef = float(gamma_ef)
        self.tem_u = bool(np.any(self.un != 0.0))
        self.phi_ref = float(np.mean(self.phin))
        # coeficientes da funcao linear f = a + bx x + by y por triangulo
        P = xy[self.tri]                                   # (nt,3,2)
        M = np.concatenate([np.ones((len(self.tri), 3, 1)), P], axis=2)
        self._Minv = np.linalg.inv(M)                      # (nt,3,3)
        self._coef = {}
        for nome, v in (("c", self.cn), ("phi", self.phin), ("u", self.un)):
            self._coef[nome] = np.einsum("kij,kj->ki", self._Minv, v[self.tri])
        self._kd = cKDTree(P.mean(axis=1))

    # -- utilitarios --
    def _triangulo(self, x, y):
        x = np.asarray(x, float).ravel()
        y = np.asarray(y, float).ravel()
        k = np.asarray(self.finder(x, y), int)
        fora = k < 0
        if np.any(fora):
            # fora da malha (ar, ou um ponto da superficie que o trifinder
            # perde por arredondamento): entre os 8 triangulos de centroide
            # mais proximo, o de maior coordenada baricentrica minima (o que
            # contem o ponto ou fica mais perto dele)
            xf, yf = x[fora], y[fora]
            cand = self._kd.query(np.c_[xf, yf], k=min(8, len(self.tri)))[1]
            p = np.stack([np.ones_like(xf), xf, yf], axis=1)          # (nf,3)
            lam = np.einsum("fr,fcrj->fcj", p, self._Minv[cand])      # (nf,8,3)
            k[fora] = cand[np.arange(len(xf)), np.argmax(lam.min(axis=2), axis=1)]
        return k

    def _avalia(self, nome, x, y):
        shp = np.shape(x)
        x = np.asarray(x, float).ravel()
        y = np.asarray(y, float).ravel()
        k = self._triangulo(x, y)
        a = self._coef[nome][k]
        return (a[:, 0] + a[:, 1] * x + a[:, 2] * y).reshape(shp)

    def c(self, x, y):
        return self._avalia("c", x, y)

    def phi(self, x, y):
        return self._avalia("phi", x, y)

    def u(self, x, y):
        return self._avalia("u", x, y)

    def grad_u(self, x, y):
        shp = np.shape(x)
        k = self._triangulo(x, y)
        a = self._coef["u"][k]
        return a[:, 1].reshape(shp), a[:, 2].reshape(shp)

    def phi_max_arco(self, xs, ys, flecha=0.0) -> float:
        """Limite superior de max phi(x) no arco amostrado em xs, ys (em ordem).

        phi e P1: dentro de um triangulo e linear. Entre duas amostras
        vizinhas j, j+1 (em triangulos k_j, k_j+1, cruzando no maximo uma
        aresta) o maximo e max(f_kj(x_j), f_kj(x_c), f_kj+1(x_j+1)), com x_c
        na aresta; como f_kj e linear na corda, f_kj(x_c) <= max(f_kj(x_j),
        f_kj(x_j+1)). Entao max phi no arco <= max sobre j de f_kj(x_j),
        f_kj(x_j+1) e f_kj+1(x_j) (as extrapolacoes cobrem o ponto de troca
        de triangulo). A curvatura do arco afasta a curva da corda de no
        maximo `flecha` (= h^2/(8 r), h o espacamento das amostras); soma-se
        |grad phi| * flecha dos triangulos tocados. potencias() confere de
        novo phi(x) <= phi_i nos pontos de Gauss."""
        x = np.asarray(xs, float).ravel()
        y = np.asarray(ys, float).ravel()
        k = self.finder(x, y)
        if k.min() < 0:
            k = self._triangulo(x, y)
        a = self._coef["phi"][k]                            # (n,3)
        a0, a1, a2 = a[:, 0], a[:, 1], a[:, 2]
        v = (a0 + a1 * x + a2 * y).max()
        if len(x) > 1:
            # f_kj(x_j+1) e f_kj+1(x_j); onde k nao muda e so o valor no ponto
            v = max(v, (a0[:-1] + a1[:-1] * x[1:] + a2[:-1] * y[1:]).max(),
                    (a0[1:] + a1[1:] * x[:-1] + a2[1:] * y[:-1]).max())
        if flecha > 0.0:
            v += flecha * np.sqrt((a1 * a1 + a2 * a2).max())
        return float(v)


# ----------------------------------------------------------------------------
# Leitura das exportacoes do NeoPZ (SlopeDrawdownMonteCarlo)
# ----------------------------------------------------------------------------
def ler_exportacao(pasta, origem=(35.0, 35.0), H=5.0, beta_graus=45.0,
                   largura=70.0, altura=35.0, snapshot=-1):
    """Le vertices.csv, triangles.csv, snaps_pw.csv e campos.csv de uma pasta
    <saida>_real<n>/ e devolve (talude, campo) no sistema do artigo.

    Exportacao: x para a direita, y para cima, aresta da crista em `origem`,
    p_w TOTAL (kPa, compressao > 0) nos vertices; c em kPa, phi em graus.
    Transformacao: x' = x - x0, y' = y0 - y; u = p_w - gamma_w y'."""
    import os
    x0, y0 = origem
    V = np.loadtxt(os.path.join(pasta, "vertices.csv"), delimiter=",", skiprows=1)
    tri = np.loadtxt(os.path.join(pasta, "triangles.csv"), delimiter=",",
                     skiprows=1).astype(int)
    linhas = np.loadtxt(os.path.join(pasta, "snaps_pw.csv"), delimiter=",", ndmin=2)
    pw = linhas[snapshot, 1:]
    C = np.loadtxt(os.path.join(pasta, "campos.csv"), delimiter=",", skiprows=1)
    no = V[:, 0].astype(int)
    assert np.all(no == np.arange(len(no))) and len(pw) == len(no)
    assert np.allclose(C[:, 1:3], V[:, 1:3])
    xy = np.c_[V[:, 1] - x0, y0 - V[:, 2]]
    u = pw - GAMMA_W * xy[:, 1]
    gamma = 20.0
    talude = Talude(H=H, beta=math.radians(beta_graus), x_esq=-x0,
                    x_dir=largura - x0, y_base=y0)
    campo = CampoMalha(xy, tri, C[:, 3], np.radians(C[:, 4]), u,
                       gamma_ef=gamma - GAMMA_W)
    campo.pw = pw
    campo.tempo = linhas[snapshot, 0]
    return talude, campo


# ----------------------------------------------------------------------------
# 1. Mecanismos classicos homogeneos (forma fechada)
# ----------------------------------------------------------------------------
def _dissip_fechada(c, r0, dth, tp):
    """P_mr/omega = c r0^2 (exp(2 dth tan phi) - 1)/(2 tan phi)  (Eq. 48)."""
    if tp < 1e-12:
        return c * r0 * r0 * dth
    return c * r0 * r0 * np.expm1(2.0 * dth * tp) / (2.0 * tp)


def _f1(r0, rh, t1, t2, tp):
    """Setor C-A-B: int (Cx - x) dA  (Eq. 51)."""
    return ((rh ** 3 * (3 * tp * math.cos(t2) + math.sin(t2))
             - r0 ** 3 * (3 * tp * math.cos(t1) + math.sin(t1)))
            / (3.0 * (9.0 * tp * tp + 1.0)))


def _tri_mom(C, P, Q):
    """int_{triangulo C P Q} (Cx - x) dA = A (2 Cx - xP - xQ)/3, com A a area
    COM SINAL (positiva quando o raio de C gira de P para Q no sentido de
    theta crescente). Com C atras da reta da face o triangulo C-O-T fica com
    sinal trocado e deve ser somado, nao subtraido."""
    A = -0.5 * ((P[0] - C[0]) * (Q[1] - C[1]) - (Q[0] - C[0]) * (P[1] - C[1]))
    return A * (2.0 * C[0] - P[0] - Q[0]) / 3.0


@dataclass
class MecanismoClassico:
    classe: str                     # 'I' (sai na face) ou 'II' (abaixo do pe)
    t1: float
    t2: float
    p3: float                       # eta (I) ou D (II)
    r0: float = 0.0
    rh: float = 0.0
    Cx: float = 0.0
    Cy: float = 0.0
    xA: float = 0.0
    xB: float = 0.0
    yB: float = 0.0
    Pmr: float = 0.0                # por omega
    Pg: float = 0.0                 # por omega
    Gamma: float = INF
    valido: bool = False
    f: tuple = ()


    def espiral(self, phi, n=400):
        th = np.linspace(self.t1, self.t2, n)
        r = self.r0 * np.exp((th - self.t1) * math.tan(phi))
        return self.Cx - r * np.cos(th), -self.Cy + r * np.sin(th)


def mecanismo_classico(classe, t1, t2, p3, talude: Talude, c, phi, gamma_ef,
                       checar=True) -> MecanismoClassico:
    """Mecanismo (I) (Eqs. 42-53) ou (II) (B a D do pe no terreno a frente).

    P_mr e P_gamma em forma fechada. Vinculos da Eq. 56 mais o dominio
    (x_esq, x_dir, y_base) e a espiral abaixo da superficie."""
    m = MecanismoClassico(classe, t1, t2, p3)
    H, b = talude.H, talude.beta
    tp = math.tan(phi)
    tmax = math.pi - b if classe == "I" else math.pi
    if not (0.0 < t1 < t2 < tmax):
        return m
    if classe == "I" and not (0.0 < p3 <= 1.0):
        return m
    if classe == "II" and p3 < 0.0:
        return m
    E = math.exp((t2 - t1) * tp)
    yB = p3 * H if classe == "I" else H
    den = E * math.sin(t2) - math.sin(t1)
    if den <= 0.0:
        return m
    r0 = yB / den                                          # Eq. 43
    rh = r0 * E                                            # Eq. 44
    xB = p3 * H / math.tan(b) if classe == "I" else talude.xT + p3
    Cx = rh * math.cos(t2) + xB                            # Eq. 45
    Cy = r0 * math.sin(t1)
    xA = Cx - r0 * math.cos(t1)                            # = -L (Eq. 46)
    m.r0, m.rh, m.Cx, m.Cy, m.xA, m.xB, m.yB = r0, rh, Cx, Cy, xA, xB, yB
    m.sB = yB - H if classe == "I" else p3
    if xA >= 0.0 or xA < talude.x_esq or xB > talude.x_dir:
        return m
    if Cy <= 0.0:
        return m
    if classe == "I" and Cy + math.tan(b) * Cx <= 0.0:
        # limite da forma fechada, nao do mecanismo: f3 (Eq. 53) nao foi
        # conferida com C atras da reta da face (em 120000 sorteios nenhum
        # mecanismo (I) valido caiu nesse caso). O codigo geral (potencias)
        # nao tem esse limite.
        return m
    if checar:
        X, Y = m.espiral(phi, 200)
        if Y.max() > talude.y_base + 1e-9:
            return m
        if np.any(Y[1:-1] < talude.y_sup(X[1:-1]) - 1e-9):
            return m
    f1 = _f1(r0, rh, t1, t2, tp)
    C = (Cx, -Cy)
    if classe == "I":
        tO = math.atan2(Cy, Cx)                            # Eq. 47
        f2 = r0 ** 3 * math.sin(t1) / 6 * (1 - math.sin(t1) ** 2 / math.sin(tO) ** 2)
        f3 = ((Cy + math.tan(b) * Cx) ** 3 / 6
              * (1 / (math.tan(tO) + math.tan(b)) ** 2
                 - 1 / (math.tan(t2) + math.tan(b)) ** 2))  # Eqs. 52-53
        m.f = (f1, f2, f3)
        Pg = f1 - f2 - f3
    else:
        A, O, T, B = (xA, 0.0), (0.0, 0.0), (talude.xT, H), (xB, H)
        tris = (_tri_mom(C, A, O), _tri_mom(C, O, T), _tri_mom(C, T, B))
        m.f = (f1,) + tris
        Pg = f1 - sum(tris)
    m.Pg = gamma_ef * Pg                                   # Eq. 50 (por omega)
    m.Pmr = _dissip_fechada(c, r0, t2 - t1, tp)            # Eq. 48 (por omega)
    if m.Pg <= 0.0:
        return m
    m.Gamma = m.Pmr / m.Pg
    m.valido = True
    return m


def minimizar_classico(talude, c, phi, gamma_ef, classes=("I", "II"), seeds=(0, 1),
                       Dmax=None):
    """Eq. 55 por evolucao diferencial (+ polimento). Devolve dict classe->mec."""
    H = talude.H
    if Dmax is None:
        Dmax = min(3.0 * H, talude.x_dir - talude.xT) if np.isfinite(talude.x_dir) else 3 * H
    res = {}
    for cl in classes:
        if cl == "I":
            lim = [(1e-3, math.pi / 2), (1e-2, math.pi - talude.beta), (1e-3, 1.0)]
        else:
            lim = [(1e-3, math.pi / 2), (1e-2, math.pi), (0.0, max(Dmax, 1e-6))]

        def fobj(p, cl=cl):
            m = mecanismo_classico(cl, p[0], p[1], p[2], talude, c, phi, gamma_ef)
            return m.Gamma if m.valido else 1e9
        best = None
        for s in seeds:
            r = optimize.differential_evolution(fobj, lim, seed=s, tol=1e-12,
                                                maxiter=2000, polish=True)
            if best is None or r.fun < best.fun:
                best = r
        res[cl] = mecanismo_classico(cl, *best.x, talude, c, phi, gamma_ef)
    return res


def fs_classico(talude, c, phi, gamma_ef, classes=("I", "II"), Fmin=0.3, Fmax=5.0):
    """FS pelas Eqs. 57-58: Gamma(c/FS, atan(tan(phi)/FS)) = 1."""
    def g(F):
        r = minimizar_classico(talude, c / F, math.atan(math.tan(phi) / F),
                               gamma_ef, classes, seeds=(0,))
        return min(m.Gamma for m in r.values()) - 1.0
    return optimize.brentq(g, Fmin, Fmax, xtol=1e-5)


def ponto_superficie(sB, talude: Talude):
    """Ponto da superficie com coordenada sB (origem no pe T):
    sB in [-H, 0]: face, profundidade H + sB;  sB > 0: terreno, x = xT + sB."""
    if sB <= 0.0:
        y = talude.H + sB
        return y / math.tan(talude.beta), y
    return talude.xT + sB, talude.H


def coordenada_superficie(xB, yB, talude: Talude):
    return yB - talude.H if xB <= talude.xT else xB - talude.xT


def classico_para_geral(m: MecanismoClassico):
    """Parametros (sB, theta2, rB) do mecanismo geral equivalente."""
    return m.sB, m.t2, m.rh


# ----------------------------------------------------------------------------
# 2. Mecanismo geral: espiral por trechos (Secao 5.2, Fig. 11)
# ----------------------------------------------------------------------------
@dataclass
class Mecanismo:
    sB: float                   # parametros: B na superficie, theta2, rB = |CB|
    t2: float
    rB: float
    xA: float = 0.0
    xB: float = 0.0
    yB: float = 0.0
    Cx: float = 0.0
    Cy: float = 0.0
    th: np.ndarray = field(default_factory=lambda: np.zeros(0))   # n+1 juncoes (A..B)
    r: np.ndarray = field(default_factory=lambda: np.zeros(0))    # n+1 raios
    tp: np.ndarray = field(default_factory=lambda: np.zeros(0))   # n tan(phi_i)
    classe: str = ""            # 'I' (B na face ou no pe), 'II' (abaixo do pe)
    motivo: str = ""            # por que e invalido
    Pmr: float = 0.0            # potencias por omega
    Pg: float = 0.0
    Pu: float = 0.0
    Pu_contorno: float = float("nan")
    viol_normalidade: float = 0.0   # maior phi(x) - phi_i nos pontos de Gauss (graus)
    Gamma: float = INF

    @property
    def valido(self):
        return self.classe != "" and np.isfinite(self.Gamma)

    @property
    def t1(self):
        return float(self.th[0])

    @property
    def r0(self):
        return float(self.r[0])

    def pontos(self, n_por_arco=4):
        """Pontos da curva AB (para desenho)."""
        X, Y = [], []
        for i in range(len(self.tp)):
            t = np.linspace(self.th[i], self.th[i + 1], n_por_arco + 1)[:-1]
            rr = self.r[i] * np.exp((t - self.th[i]) * self.tp[i])
            X.append(self.Cx - rr * np.cos(t))
            Y.append(-self.Cy + rr * np.sin(t))
        X.append([self.Cx - self.r[-1] * math.cos(self.th[-1])])
        Y.append([-self.Cy + self.r[-1] * math.sin(self.th[-1])])
        return np.concatenate(X), np.concatenate(Y)


def _dtheta_arco_tras(r, tp, l):
    """Incremento angular (para tras, theta decrescente) de um arco de
    comprimento l que termina com raio r:
    l = r (1 - exp(-dth tan phi))/sin phi."""
    if tp < 1e-12:
        return l / r
    sp = tp / math.sqrt(1.0 + tp * tp)
    q = l * sp / r
    if q >= 1.0:
        return INF
    return -math.log1p(-q) / tp


def construir_espiral(sB, t2, rB, campo, talude: Talude, l_arco=None,
                      n_amostras=6, max_iter=60, max_arcos=20000) -> Mecanismo:
    """Marcha a espiral por trechos PARA TRAS, de B ate A na crista.

    B = ponto_superficie(sB); C = B + rB (cos t2, -sin t2), isto e,
    P(theta) = (Cx - r cos theta, -Cy + r sin theta) com P(t2) = B.
    Arco i (de theta_i+1 para theta_i, comprimento l_arco = H/100): dado um
    phi provisorio, o arco e amostrado em n_amostras pontos igualmente
    espacados em theta (extremos incluidos, cerca de 1 cm entre eles) e
    phi_i = max de phi(x) ao longo do arco (CampoMalha.phi_max_arco, limite
    superior exato para phi P1). E o que a normalidade pede (Secao 5.2:
    [U].n >= |[U]| sin phi(x) em cada ponto do arco), e nada mais. Como o
    arco muda com phi_i, repete-se ate phi_i parar de crescer (monotono);
    sem convergencia em max_iter, o mecanismo e invalido. A
    marcha para quando a curva chega a crista (y = 0, x < 0); o ultimo arco
    e cortado na intersecao exata (brentq) e A = (xA, 0).
    Parametrizar por B (e nao por A) deixa a classe fixa pelo parametro:
    sB <= 0 -> mecanismo (I) (sB = 0: pelo pe, eta = 1), sB > 0 -> (II)."""
    H = talude.H
    if l_arco is None:
        l_arco = H / 100.0
    m = Mecanismo(sB, t2, rB)
    if not (-H < sB and 0.0 < t2 < math.pi and rB > 0.0):
        m.motivo = "parametros"
        return m
    if sB <= 0.0 and t2 >= math.pi - talude.beta:
        m.motivo = "theta2 >= pi - beta"
        return m
    xB, yB = ponto_superficie(sB, talude)
    if xB > talude.x_dir:
        m.motivo = "B fora do dominio"
        return m
    Cx = xB + rB * math.cos(t2)
    Cy = rB * math.sin(t2) - yB
    m.xB, m.yB, m.Cx, m.Cy = xB, yB, Cx, Cy
    if Cy <= 0.0:          # C abaixo da crista: a curva nao chega a crista
        m.motivo = "C abaixo da crista"
        return m
    ths, rs, tps = [t2], [rB], []
    th, r = t2, rB
    frac = np.linspace(0.0, 1.0, n_amostras)
    ysup = talude.y_sup_escalar
    uniforme = campo.uniforme
    tp_u = math.tan(campo.phi_ref) if uniforme else 0.0
    xe = ye = 0.0
    for _ in range(max_arcos):
        if uniforme:
            tp = tp_u
            dth = _dtheta_arco_tras(r, tp, l_arco)
        else:
            # partida: phi no ponto inicial (1o arco) ou o max de phi no arco
            # tracado com o phi_i do arco anterior (quase o mesmo arco: em
            # geral basta mais uma avaliacao para confirmar)
            ph = -INF
            if tps:
                tp = tps[-1]
                dth = _dtheta_arco_tras(r, tp, l_arco)
                if np.isfinite(dth):
                    t = th - frac * dth
                    rr = r * np.exp((t - th) * tp)
                    ph = campo.phi_max_arco(Cx - rr * np.cos(t), -Cy + rr * np.sin(t))
            if ph == -INF:
                ph = campo.phi_max_arco([Cx - r * math.cos(th)], [-Cy + r * math.sin(th)])
            conv = False
            for _it in range(max_iter):
                tp = math.tan(ph)
                dth = _dtheta_arco_tras(r, tp, l_arco)
                if not np.isfinite(dth):
                    conv = True
                    break
                t = th - frac * dth
                rr = r * np.exp((t - th) * tp)
                hs = l_arco / (n_amostras - 1)
                ph2 = campo.phi_max_arco(Cx - rr * np.cos(t), -Cy + rr * np.sin(t),
                                         hs * hs / (8.0 * rr.min()))
                if ph2 <= ph + 1e-12:
                    conv = True
                    break
                if ph2 - ph < 1e-4:
                    # subida pequena: trocar phi_i por ph2 move os pontos do
                    # arco de ~ (l^2/r) d(tan phi) < 1e-6 m, e phi(x) muda de
                    # |grad phi| vezes isso, muito menos que a folga
                    # ph2 - ph somada aqui. Evita a cauda da iteracao.
                    ph = ph2 + (ph2 - ph) + 1e-9
                    conv = True
                    break
                ph = ph2
            if not conv:
                m.motivo = "phi_i nao convergiu"
                return m
            tp = math.tan(ph)
            dth = _dtheta_arco_tras(r, tp, l_arco)
        t_ini = th - dth
        if not np.isfinite(dth) or t_ini <= 0.0:
            m.motivo = "nao chega a crista"
            return m
        r_ini = r * math.exp(-dth * tp)
        xe, ye = Cx - r_ini * math.cos(t_ini), -Cy + r_ini * math.sin(t_ini)
        if ye < ysup(xe):
            if len(tps) == 0:
                m.motivo = "nao entra no solo"
                return m

            def g(t, th=th, r=r, tp=tp):
                rr = r * math.exp((t - th) * tp)
                x, y = Cx - rr * math.cos(t), -Cy + rr * math.sin(t)
                return y - ysup(x)
            t_ini = optimize.brentq(g, t_ini, th, xtol=1e-13)
            r_ini = r * math.exp((t_ini - th) * tp)
            xe, ye = Cx - r_ini * math.cos(t_ini), 0.0
            if t_ini < th:
                ths.append(t_ini); rs.append(r_ini); tps.append(tp)
            break
        if ye > talude.y_base or xe < talude.x_esq or xe > talude.x_dir:
            m.motivo = "sai do dominio"
            return m
        ths.append(t_ini); rs.append(r_ini); tps.append(tp)
        th, r = t_ini, r_ini
    else:
        m.motivo = "arcos demais"
        return m
    if xe >= 0.0:
        m.motivo = "A fora da crista"
        return m
    if xe < talude.x_esq:
        m.motivo = "sai do dominio"
        return m
    m.th, m.r, m.tp = np.array(ths[::-1]), np.array(rs[::-1]), np.array(tps[::-1])
    m.xA = xe
    m.classe = "I" if sB <= 0.0 else "II"
    return m


_GL = {n: np.polynomial.legendre.leggauss(n) for n in range(1, 13)}


def _gauss_intervalos(a, b, n):
    """Pontos/pesos de Gauss-Legendre (n) em cada intervalo [a_k, b_k]."""
    xi, wi = _GL[n]
    h = 0.5 * (b - a)
    pts = (0.5 * (a + b))[:, None] + h[:, None] * xi[None, :]
    w = h[:, None] * wi[None, :]
    return pts, w


def _cruzamentos(tq, m, talude: Talude):
    """Intervalos de rho em que o raio theta (a partir de C) esta dentro do
    bloco, para cada theta em tq.

    O contorno do bloco e a espiral A-B mais a superficie de B ate A:
    crista (y = 0, xA <= x <= 0), face (de O ate B no (I), ou ate T no
    (II)) e terreno a frente do pe (y = H, xT <= x <= xB, so no (II)).
    C esta no ar (Cy > 0), entao, andando pelo raio a partir de C, cada
    cruzamento com o contorno alterna fora/dentro. Isso vale tambem quando
    C esta atras da reta da face (Cy + tan(beta) Cx < 0): ai um raio pode
    entrar pela crista, sair pela face, entrar pelo terreno e sair pela
    espiral (dois intervalos).
    Devolve (ray, a, b, da, db): indice do raio, extremos do intervalo e
    d rho/d theta de cada extremo (para a forma de contorno de P_u)."""
    H, tb = talude.H, math.tan(talude.beta)
    Cx, Cy = m.Cx, m.Cy
    s, co = np.sin(tq), np.cos(tq)
    n = len(tq)
    R = np.full((n, 4), np.inf)
    D = np.zeros((n, 4))
    eps = 1e-12
    # crista
    rc = Cy / s
    xc = Cx - rc * co
    ok = (xc >= m.xA - eps) & (xc <= eps)
    R[ok, 0] = rc[ok]
    D[ok, 0] = (-Cy * co / s ** 2)[ok]
    # face
    K = Cy + Cx * tb
    d = s + co * tb
    with np.errstate(divide="ignore", invalid="ignore"):
        rf = K / d
    yf = -Cy + rf * s
    yfim = m.yB if m.classe == "I" else H
    ok = np.isfinite(rf) & (rf > 0.0) & (yf >= -eps) & (yf <= yfim + eps)
    R[ok, 1] = rf[ok]
    D[ok, 1] = (-K * (co - s * tb) / d ** 2)[ok]
    # terreno a frente do pe
    if m.classe == "II":
        rt = (Cy + H) / s
        xt = Cx - rt * co
        ok = (xt >= talude.xT - eps) & (xt <= m.xB + eps)
        R[ok, 2] = rt[ok]
        D[ok, 2] = (-(Cy + H) * co / s ** 2)[ok]
    # espiral
    th, r, tp = m.th, m.r, m.tp
    ok = (tq >= th[0]) & (tq <= th[-1])
    k = np.clip(np.searchsorted(th, tq[ok], side="right") - 1, 0, len(tp) - 1)
    rq = r[k] * np.exp((tq[ok] - th[k]) * tp[k])
    R[ok, 3] = rq
    D[ok, 3] = rq * tp[k]
    o = np.argsort(R, axis=1)
    R = np.take_along_axis(R, o, axis=1)
    D = np.take_along_axis(D, o, axis=1)
    nv = np.isfinite(R).sum(axis=1)
    ray, a, b, da, db = [], [], [], [], []
    for j in (0, 2):
        sel = nv >= j + 2
        idx = np.nonzero(sel)[0]
        ray.append(idx)
        a.append(R[sel, j]); b.append(R[sel, j + 1])
        da.append(D[sel, j]); db.append(D[sel, j + 1])
    return (np.concatenate(ray), np.concatenate(a), np.concatenate(b),
            np.concatenate(da), np.concatenate(db))


TOL_NORMALIDADE = 1e-6      # graus


def potencias(m: Mecanismo, campo, talude: Talude, h_quad=None, ng_theta=3,
              ng_rho=2, dissipacao="exata", contorno=True):
    """P_mr, P_gamma e P_u (por omega) de um mecanismo construido.

    Dissipacao (arco a arco, 2 pontos de Gauss por arco para c e phi locais):
      'exata'  : dP_mr = c(x) tan(phi_i)/tan(phi(x)) r^2 dtheta
                 (pi = c cot(phi(x)) [U].n, [U].n = omega r sin(phi_i),
                  ds = r dtheta / cos(phi_i)) -> limite superior rigoroso
                 desde que phi(x) <= phi_i em todo o arco;
      'simples': dP_mr = c(x) r^2 dtheta (= 'exata' quando phi(x) = phi_i).
    Normalidade: confere phi(x) <= phi_i nos pontos de Gauss. Se falhar por
    mais de TOL_NORMALIDADE, o mecanismo fica invalido (motivo
    'normalidade violada'; m.viol_normalidade guarda a maior violacao em
    graus); nada e cortado em silencio.
    Bloco: coordenadas polares em torno de C. Em theta: Gauss com ng_theta
    pontos entre cortes (juncoes dos arcos e os angulos de A, O, T e B, onde
    o contorno tem quina), cobrindo todo o intervalo angular do contorno.
    Em cada raio, os intervalos de rho dentro do bloco vem de _cruzamentos
    (um ou dois). Em rho: P_gamma exato, P_u por Gauss composto com ng_rho
    pontos em segmentos <= h_quad (padrao H/50).
      P_gamma = gamma' int rho cos(theta) dA       (e_y.U = omega rho cos)
      P_u     = int -grad u . e_theta rho dA       (grad u da malha)
    Com contorno=True (e um campo com potencial u) calcula tambem
      P_u = soma nos extremos dos intervalos de +-u rho drho/dtheta dtheta
    (+ na saida, - na entrada), igual a de volume pelo teorema da
    divergencia (div U = 0)."""
    H = talude.H
    if h_quad is None:
        h_quad = H / 50.0
    Cx, Cy = m.Cx, m.Cy
    th, r, tp = m.th, m.r, m.tp
    # ---- dissipacao: Gauss por arco ----
    xi, wi = _GL[2]
    a, b = th[:-1], th[1:]
    h = 0.5 * (b - a)
    tg = (0.5 * (a + b))[:, None] + h[:, None] * xi[None, :]
    rg = r[:-1, None] * np.exp((tg - a[:, None]) * tp[:, None])
    xg, yg = Cx - rg * np.cos(tg), -Cy + rg * np.sin(tg)
    # int r^2 dtheta exato por arco; c(x) e a razao pela media de Gauss
    tpc = np.maximum(tp, 1e-300)
    I_r2 = np.where(tp < 1e-12, r[:-1] ** 2 * (b - a),
                    r[:-1] ** 2 * np.expm1(2 * (b - a) * tpc) / (2 * tpc))
    cg = campo.c(xg, yg)
    phig = campo.phi(xg, yg)
    m.viol_normalidade = float(max(0.0, np.degrees(np.max(phig - np.arctan(tp)[:, None]))))
    if m.viol_normalidade > TOL_NORMALIDADE:
        m.motivo = "normalidade violada"
        m.Gamma = INF
        return m
    if dissipacao == "exata":
        tl = np.tan(phig)
        razao = np.where(tl > 1e-12, tp[:, None] / np.maximum(tl, 1e-300), 1.0)
        dens = cg * razao
    elif dissipacao == "simples":
        dens = cg
    else:
        raise ValueError(dissipacao)
    m.Pmr = float(np.sum(I_r2 * (dens * wi[None, :]).sum(axis=1) / wi.sum()))
    # ---- bloco: quadratura polar ----
    thO = math.atan2(Cy, Cx)
    angs = [th[0], th[-1], thO]
    if m.classe == "II":
        angs.append(math.atan2(Cy + H, Cx - talude.xT))
    lo, hi = min(angs), max(angs)
    cortes = np.unique(np.concatenate([th, angs]))
    cortes = cortes[(cortes >= lo) & (cortes <= hi)]
    # fora de [theta_A, theta_B] (C atras da reta da face) nao ha juncoes de
    # arcos: corta os intervalos longos no passo angular tipico dos arcos
    dmax = float(np.max(np.diff(th)))
    nsub = np.maximum(1, np.ceil(np.diff(cortes) / dmax).astype(int))
    if np.any(nsub > 1):
        cortes = np.concatenate([np.linspace(cortes[i], cortes[i + 1], nsub[i] + 1)[:-1]
                                 for i in range(len(nsub))] + [cortes[-1:]])
    tq, wq = _gauss_intervalos(cortes[:-1], cortes[1:], ng_theta)
    tq, wq = tq.ravel(), wq.ravel()
    ray, ra, rb, da, db = _cruzamentos(tq, m, talude)
    ti, wt = tq[ray], wq[ray]
    m.Pg = float(campo.gamma_ef * np.sum(wt * np.cos(ti) * (rb ** 3 - ra ** 3) / 3.0))
    if campo.tem_u:
        L = np.maximum(rb - ra, 0.0)
        ns = np.maximum(1, np.ceil(L / h_quad).astype(int))
        iseg = np.repeat(np.arange(len(L)), ns)
        jseg = np.arange(ns.sum()) - np.repeat(np.cumsum(ns) - ns, ns)
        hs = (L / ns)[iseg]
        r0s = ra[iseg] + jseg * hs
        xr, wr = _GL[ng_rho]
        rho = (r0s[:, None] + 0.5 * hs[:, None] * (1 + xr[None, :])).ravel()
        w = ((0.5 * hs * wt[iseg])[:, None] * wr[None, :]).ravel()
        tt = np.repeat(ti[iseg], ng_rho)
        st, ct = np.sin(tt), np.cos(tt)
        x, y = Cx - rho * ct, -Cy + rho * st
        ux, uy = campo.grad_u(x, y)
        m.Pu = float(np.sum(w * rho * rho * (-(ux * st + uy * ct))))
        if contorno and campo.tem_potencial:
            ci, si = np.cos(ti), np.sin(ti)
            ub = campo.u(Cx - rb * ci, -Cy + rb * si)
            ua = campo.u(Cx - ra * ci, -Cy + ra * si)
            m.Pu_contorno = float(np.sum(wt * (ub * rb * db - ua * ra * da)))
    else:
        m.Pu = 0.0
        m.Pu_contorno = 0.0
    Pext = m.Pg + m.Pu
    m.Gamma = m.Pmr / Pext if Pext > 0.0 else INF
    return m


def avaliar(sB, t2, rB, campo, talude, l_arco=None, **kw) -> Mecanismo:
    """Constroi a espiral e calcula Gamma (Eq. 54) para um mecanismo."""
    m = construir_espiral(sB, t2, rB, campo, talude, l_arco)
    if m.classe == "":
        return m
    return potencias(m, campo, talude, **kw)


# ----------------------------------------------------------------------------
# Minimizacao (Eq. 55) para o mecanismo geral
# ----------------------------------------------------------------------------
def _dmax(talude):
    H = talude.H
    return min(3.0 * H, talude.x_dir - talude.xT) if np.isfinite(talude.x_dir) else 3.0 * H


def _candidatos(talude, phi_ref, classe, n1=9, n2=9, n3=6):
    """Pontos de partida: grade (theta1, theta2, eta|D) do mecanismo classico
    com phi_ref, convertida em (sB, theta2, rB)."""
    H, b = talude.H, talude.beta
    out = []
    t1s = np.linspace(0.15, 1.45, n1)
    if classe == "I":
        p3s = np.linspace(1.0, 0.3, n3)
        t2max = math.pi - b - 0.02
    else:
        p3s = np.linspace(0.0, max(min(2 * H, _dmax(talude) - 0.1), 0.0), n3)
        t2max = math.pi - 0.05
    for t1 in t1s:
        for t2 in np.linspace(t1 + 0.2, t2max, n2):
            for p3 in p3s:
                mc = mecanismo_classico(classe, t1, t2, p3, talude, 1.0, phi_ref, 1.0)
                if mc.valido:
                    out.append(classico_para_geral(mc))
    return out


@dataclass
class Resultado:
    Gamma: float
    mecanismo: Mecanismo
    por_classe: dict            # 'I' (face ou pe, eta <= 1) e 'II' (abaixo do pe)
    n_aval: int
    tempo_s: float
    n_normalidade: int = 0      # avaliacoes descartadas por phi(x) > phi_i (deve ser 0)
    n_phi_nao_conv: int = 0     # avaliacoes descartadas: phi_i sem convergencia


def _z_para_param(z, cl, talude):
    """z = (w, theta2, ln(rB/H)); w e cortado no intervalo da classe:
    I: sB = -H clip(w, 0, 0.98) (w <= 0 -> pelo pe); II: sB = clip(w, 0, Dmax/H) H."""
    H = talude.H
    if cl == "I":
        sB = -H * min(max(z[0], 0.0), 0.98)
    else:
        sB = H * min(max(z[0], 0.0), _dmax(talude) / H)
    return sB, z[1], H * math.exp(z[2])


def _param_para_z(sB, t2, rB, talude):
    H = talude.H
    w = -sB / H if sB <= 0.0 else sB / H
    return np.array([w, t2, math.log(rB / H)])


def minimizar(campo, talude, classes=("I", "II"), l_arco=None, l_grade=None,
              n_partidas=3, dissipacao="exata", xatol=1e-5, fatol=1e-6,
              grade=(9, 9, 6), verbose=False) -> Resultado:
    """Minimiza Gamma (Eq. 55) para cada classe de mecanismo.

    1) busca global: grade (theta1, theta2, eta|D) do mecanismo classico com
       phi medio, convertida para (sB, theta2, rB) e avaliada com o campo real
       e arcos grossos (l_grade = H/25);
    2) refino local por Nelder-Mead em z = (w, theta2, ln(rB/H)) com arcos
       l_arco (= H/100), a partir das n_partidas melhores de cada classe.
       w e cortado ao intervalo da classe, entao a fronteira (eta = 1, B no
       pe) e atingida exatamente: Nelder-Mead enxerga um patamar."""
    H = talude.H
    l_arco = H / 100.0 if l_arco is None else l_arco
    l_grade = H / 25.0 if l_grade is None else l_grade
    t0 = time.perf_counter()
    nav = [0]
    nviol = {"normalidade violada": 0, "phi_i nao convergiu": 0}
    PEN = 1e3

    def G(z, cl, l):
        nav[0] += 1
        m = avaliar(*_z_para_param(z, cl, talude), campo, talude, l,
                    dissipacao=dissipacao, contorno=False)
        if m.motivo in nviol:
            nviol[m.motivo] += 1
        return m.Gamma if m.valido else PEN

    por = {}
    for cl in classes:
        cands = _candidatos(talude, campo.phi_ref, cl, *grade)
        vals = []
        for p in cands:
            z = _param_para_z(*p, talude)
            vals.append((G(z, cl, l_grade), z))
        vals.sort(key=lambda v: v[0])
        best = None
        for g0, z0 in vals[:n_partidas]:
            if g0 >= PEN:
                break
            simplex = np.vstack([z0, z0 + np.diag([0.04, 0.04, 0.04])])
            r = optimize.minimize(G, z0, args=(cl, l_arco), method="Nelder-Mead",
                                  options=dict(xatol=xatol, fatol=fatol, maxiter=3000,
                                               initial_simplex=simplex))
            if best is None or r.fun < best.fun:
                best = r
        if best is not None and best.fun < PEN:
            por[cl] = avaliar(*_z_para_param(best.x, cl, talude), campo, talude, l_arco,
                              dissipacao=dissipacao, contorno=True)
    if not por:
        raise RuntimeError("nenhum mecanismo valido")
    if verbose:
        for cl, m in por.items():
            print(f"  {cl:>3}: Gamma = {m.Gamma:.4f}  xA = {m.xA:.3f}  xB = {m.xB:.3f}"
                  f"  yB = {m.yB:.3f}")
    mbest = min(por.values(), key=lambda mm: mm.Gamma)
    return Resultado(mbest.Gamma, mbest, por, nav[0], time.perf_counter() - t0,
                     nviol["normalidade violada"], nviol["phi_i nao convergiu"])


# ----------------------------------------------------------------------------
# Uma realizacao a partir dos arquivos exportados
# ----------------------------------------------------------------------------
def avaliar_realizacao(pasta, verbose=False, **kw):
    """Le a exportacao de uma realizacao (<saida>_real<n>/) e devolve o
    Resultado da minimizacao (Gamma cinematico e mecanismo critico).
    kw extras vao para ler_exportacao (geometria) ou minimizar."""
    chaves_ler = {"origem", "H", "beta_graus", "largura", "altura", "snapshot"}
    kl = {k: v for k, v in kw.items() if k in chaves_ler}
    km = {k: v for k, v in kw.items() if k not in chaves_ler}
    talude, campo = ler_exportacao(pasta, **kl)
    res = minimizar(campo, talude, verbose=verbose, **km)
    res.talude, res.campo = talude, campo
    return res
