"""Fator de estabilidade Gamma de Vargas Ceron et al. (2025) pela abordagem cinematica com mecanismos rotacionais em
espiral logaritmica (Eqs. 42-56), com as pressoes de agua de uma rodada do SlopeRandomField (pw.csv e triangulos.csv)
ou sem agua.

Mecanismo: bloco rigido girando (sentido anti-horario, omega = 1) em torno de O, limitado pela superficie do talude e pela
espiral r(th) = r_B exp((th - th_B) tan phi), que passa por B na superficie (na face ou no terreno do pe) e volta ate
cortar a superficie em A. Com p a pressao de agua (tensao efetiva de Terzaghi, superficie com p = 0):
  P_mr = c (r_B^2 - r_A^2) / (2 tan phi)                       (dissipacao na espiral, Eq. 48)
  P_g + P_u = gamma * int (-cos th) (r_sp^3 - r_s^3)/3 dth      (peso total, gamma e_y . U)
            + tan(phi) * int p(S(th)) r_sp^2 dth               (-int grad p . U = int_AB p |U| sin(phi) ds)
  Gamma = min P_mr / (P_g + P_u).
O peso efetivo mais as forcas de percolacao do artigo, gamma' e_y - grad u, somam gamma e_y - grad p; pelo teorema da
divergencia (div U = 0 no bloco rigido) o termo de p fica so na espiral, porque p = 0 na superficie.

Uso: python3 -I espiral_ceron.py seco <c> <phi_graus> <gamma>
     python3 -I espiral_ceron.py <pasta da rodada> <t_h> <c> <phi_graus> <gamma>"""
import csv
import sys

import numpy as np
from scipy.optimize import minimize

XC, YC, XT, YT, XMAX = 10., 10., 15., 5., 25.


def superficie_raio(Ox, Oy, th):
    """menor r > 0 em que o raio O + r (cos th, sin th) corta a superficie (crista, face, terreno do pe)"""
    dx, dy = np.cos(th), np.sin(th)
    best = np.full(np.shape(th), np.inf)
    for (Px, Py, Qx, Qy) in ((-1e3, YC, XC, YC), (XC, YC, XT, YT), (XT, YT, 1e3, YT)):
        ex, ey = Qx - Px, Qy - Py
        den = dx * (-ey) - dy * (-ex)
        with np.errstate(divide='ignore', invalid='ignore'):
            r = ((Px - Ox) * (-ey) - (Py - Oy) * (-ex)) / den
            s = (dx * (Py - Oy) - dy * (Px - Ox)) / den
        ok = (np.abs(den) > 1e-14) & (r > 1e-12) & (s >= -1e-12) & (s <= 1 + 1e-12)
        best = np.where(ok & (r < best), r, best)
    return best


def ysup(x):
    return np.where(x <= XC, YC, np.where(x >= XT, YT, YC - (x - XC) * (YC - YT) / (XT - XC)))


class Pressao:
    def __init__(self, pasta=None, t_h=None):
        self.f = None
        if pasta is None:
            return
        import matplotlib.tri as mtri
        nos = {}
        for r in csv.DictReader(open(pasta + '/pw.csv')):
            if abs(float(r['t_h']) - t_h) < 1e-9:
                nos[int(r['no'])] = (float(r['x']), float(r['y']), float(r['pw_kPa']))
        ids = sorted(nos)
        ind = {n: i for i, n in enumerate(ids)}
        X = np.array([nos[n][0] for n in ids]); Y = np.array([nos[n][1] for n in ids])
        P = np.array([nos[n][2] for n in ids])
        T = np.array([[ind[int(r['n0'])], ind[int(r['n1'])], ind[int(r['n2'])]]
                      for r in csv.DictReader(open(pasta + '/triangulos.csv'))])
        self.f = mtri.LinearTriInterpolator(mtri.Triangulation(X, Y, T), P)

    def __call__(self, x, y):
        if self.f is None:
            return np.zeros_like(x)
        v = self.f(x, y)
        return np.ma.filled(v, np.nan)


def gamma_mecanismo(Ox, Oy, Bx, By, c, phi, g, pres, n=600):
    """Gamma do mecanismo com centro O e espiral pela ponta B (na superficie); None se nao for admissivel"""
    tp = np.tan(phi)
    if Oy <= ysup(np.array(Ox)) + 1e-6:  # O acima da superficie
        return None
    thB = np.arctan2(By - Oy, Bx - Ox)
    if thB > 0:
        thB -= 2 * np.pi  # angulos para baixo: de -pi a 0 (cos < 0 a esquerda de O)
    rB = np.hypot(Bx - Ox, By - Oy)
    # volta de B ate a espiral cortar a superficie (A): r_sp(th) > r_s(th) dentro do solo
    ths = thB - np.linspace(0., np.pi, 4001)[1:]
    rsp = rB * np.exp((ths - thB) * tp)
    rs = superficie_raio(Ox, Oy, ths)
    fora = ~(rsp > rs)
    if not fora.any():
        return None
    k = int(np.argmax(fora))
    if k == 0:
        return None  # a espiral sai do solo logo ao deixar B: nao ha bloco
    # refina A por bissecao entre ths[k-1] (dentro) e ths[k] (fora)
    a, b = ths[k - 1], ths[k]
    for _ in range(60):
        m = 0.5 * (a + b)
        if rB * np.exp((m - thB) * tp) > superficie_raio(Ox, Oy, np.array([m]))[0]:
            a = m
        else:
            b = m
    thA = 0.5 * (a + b)
    th = np.linspace(thA, thB, n + 1)
    rsp = rB * np.exp((th - thB) * tp)
    rs = np.minimum(superficie_raio(Ox, Oy, th), rsp)
    Sx, Sy = Ox + rsp * np.cos(th), Oy + rsp * np.sin(th)
    if Sy.min() < -1e-9 or Sx.min() < -1e-9 or Sx.max() > XMAX + 1e-9:
        return None  # a espiral sai do dominio
    rA = rsp[0]
    Pmr = c * (rB ** 2 - rA ** 2) / (2. * tp)
    Pg = g * np.trapezoid(-np.cos(th) * (rsp ** 3 - rs ** 3) / 3., th)
    p = pres(Sx, Sy)
    if np.isnan(p).any():
        return None
    Pu = tp * np.trapezoid(p * rsp ** 2, th)
    if Pg + Pu <= 1e-9:
        return None
    return Pmr / (Pg + Pu)


def ponto_B(s):
    """s em [0, 1): na face, do pe (0) a crista; s em [-1, 0): no terreno do pe, x = 15 - 8 s"""
    if s >= 0:
        return XT + s * (XC - XT), YT + s * (YC - YT)
    return XT - 8. * s, YT


def busca(c, phi, g, pres):
    best = []
    for Ox in np.linspace(6., 30., 25):
        for Oy in np.linspace(6., 40., 35):
            for s in np.concatenate([np.linspace(-0.6, -0.02, 8), np.linspace(0., 0.9, 10)]):
                Bx, By = ponto_B(s)
                G = gamma_mecanismo(Ox, Oy, Bx, By, c, phi, g, pres, n=200)
                if G is not None:
                    best.append((G, Ox, Oy, s))
    best.sort()
    res = []
    for G0, Ox, Oy, s in best[:8]:
        def f(v):
            Bx, By = ponto_B(min(max(v[2], -1.), 0.99))
            G = gamma_mecanismo(v[0], v[1], Bx, By, c, phi, g, pres, n=600)
            return 1e3 if G is None else G
        r = minimize(f, [Ox, Oy, s], method='Nelder-Mead',
                     options={'xatol': 1e-5, 'fatol': 1e-7, 'maxiter': 4000, 'initial_simplex': None})
        res.append((r.fun, r.x))
    res.sort(key=lambda t: t[0])
    return res[0], len(best)


if __name__ == '__main__':
    if sys.argv[1] == 'seco':
        c, phi, g = float(sys.argv[2]), np.radians(float(sys.argv[3])), float(sys.argv[4])
        pres = Pressao()
        rotulo = 'seco'
    else:
        pasta, t_h = sys.argv[1], float(sys.argv[2])
        c, phi, g = float(sys.argv[3]), np.radians(float(sys.argv[4])), float(sys.argv[5])
        pres = Pressao(pasta, t_h)
        rotulo = pasta.rstrip('/').split('/')[-1] + ' t=%g h' % t_h
    (G, v), nval = busca(c, phi, g, pres)
    Bx, By = ponto_B(min(max(v[2], -1.), 0.99))
    print('%s: Gamma = %.4f (c = %g, phi = %g, gamma = %g; %d mecanismos validos na grade); O = (%.3f, %.3f), '
          'B = (%.3f, %.3f)' % (rotulo, G, c, np.degrees(phi), g, nval, v[0], v[1], Bx, By))
