#!/usr/bin/env python3
"""rapido.py: avaliador rapido do Gamma cinematico (o mesmo metodo do cinematico.py) para o Monte Carlo.

O que muda em relacao ao cinematico.py (o metodo e as regras ficam iguais):
  1. CampoGrade: c, phi e u sao P1 numa grade UNIFORME de passo h (cada celula dividida pela diagonal paralela a
     face de 45 graus; nas celulas da face so o triangulo de baixo e solo). Localizar um ponto e aritmetica
     (i = floor((x - x0)/h), j = floor(y/h), s > t ou nao), sem trifinder. Os valores sao exatamente os de um
     CampoMalha na mesma triangulacao (para_malha() devolve esse CampoMalha, para a conferencia).
  2. construir_espiral_rapido: a mesma marcha por arcos de comprimento l (H/100), o mesmo phi_i = limite superior
     de max phi(x) no arco (6 amostras, extrapolacoes nas trocas de triangulo, mais |grad phi| vezes a flecha), mas
     resolvida para TODOS os arcos de uma vez: dado phi_i de cada arco, r e theta saem de somas acumuladas
     (r_i = r_{i+1} - l sin phi_i; dtheta_i = ln(r_{i+1}/r_i)/tan phi_i), e repete-se phi_i <- limite(arco_i) ate
     parar de mudar (iteracao de Jacobi; o arco i so depende dos arcos entre ele e B, entao converge de B para A).
     No fim confere-se que phi_i >= limite do proprio arco em todos os arcos (normalidade); se nao convergir, cai
     na marcha sequencial do cinematico.py (com o mesmo campo) e conta.
  3. potencias: a mesma funcao do cinematico.py (dissipacao 'exata', quadratura polar), so que com o CampoGrade.
  4. minimizar_rapido: partida quente (o mecanismo critico deterministico de cada classe) + grade pequena de
     mecanismos classicos com arcos grossos + Nelder-Mead a partir das melhores partidas."""
from __future__ import annotations

import math
import time

import numpy as np
from scipy import optimize

import cinematico as kin

INF = float("inf")


# ---------------------------------------------------------------------------------------------------------------------
# Campo P1 numa grade uniforme
# ---------------------------------------------------------------------------------------------------------------------
class CampoGrade:
    """Campos nodais P1 na grade x0 + h i (i = 0..nx-1), h j (j = 0..ny-1) (y para baixo, crista em y = 0).

    Triangulos de cada celula (i, j), coordenadas locais s = (x - x_i)/h, t = (y - y_j)/h:
      A (s >= t): nos 00, 10, 11;  B (t > s): nos 00, 01, 11.
    Celulas acima da superficie sao ar; nas celulas da face (y_j = x_i, 0 <= x_i < xT) so B e solo. Pontos no ar (so
    por arredondamento) usam o triangulo de solo da coluna mais perto (crista e terreno: a celula de baixo; face: B)."""
    uniforme = False
    tem_potencial = True

    def __init__(self, x0, h, c, phi, u, gamma_ef, talude):
        self.x0, self.h = float(x0), float(h)
        self.nx, self.ny = c.shape
        self.cn = np.ascontiguousarray(c, float).ravel()
        self.phin = np.ascontiguousarray(phi, float).ravel()
        self.un = np.ascontiguousarray(u, float).ravel()
        self.gamma_ef = float(gamma_ef)
        self.tem_u = bool(np.any(self.un != 0.0))
        self.tal = talude
        H, xT = talude.H, talude.xT
        xi = self.x0 + self.h * np.arange(self.nx - 1)
        # primeira linha de solo de cada coluna de celulas e se a coluna e de face
        jH = int(round(H / h))
        self.jmin = np.where(xi + h <= 1e-9, 0, np.where(xi >= xT - 1e-9, jH, np.rint(xi / h).astype(int)))
        self.face = (xi > -1e-9) & (xi < xT - 1e-9)
        self._ih = 1.0 / self.h
        self._jminf = self.jmin.astype(float)
        self._jface = np.where(self.face, self.jmin, -1)      # linha da celula de face de cada coluna (-1: nenhuma)
        self._solo_nos = self._nos_solo()
        self.phi_ref = float(np.mean(self.phin[self._solo_nos]))

    def _nos_solo(self):
        X = self.x0 + self.h * np.arange(self.nx)
        Y = self.h * np.arange(self.ny)
        return (Y[None, :] >= self.tal.y_sup(X)[:, None] - 1e-9).ravel()

    # -- localizacao --
    def _loc(self, x, y):
        fx = (x - self.x0) * self._ih
        fy = y * self._ih
        i = np.minimum(np.maximum(np.floor(fx), 0.0), self.nx - 2).astype(np.intp)
        j = np.minimum(np.maximum(np.floor(fy), self._jminf[i]), self.ny - 2).astype(np.intp)
        s = fx - i
        t = fy - j
        A = (s >= t) & (j != self._jface[i])
        n00 = i * self.ny + j
        nm = np.where(A, n00 + self.ny, n00 + 1)
        n11 = n00 + (self.ny + 1)
        a = np.where(A, s, t)
        b = np.where(A, t, s)
        return n00, nm, n11, a, b, A

    def _valor(self, v, x, y):
        shp = np.shape(x)
        x = np.asarray(x, float).ravel()
        y = np.asarray(y, float).ravel()
        n00, nm, n11, a, b, _ = self._loc(x, y)
        v0, vm, v1 = v[n00], v[nm], v[n11]
        return (v0 + a * (vm - v0) + b * (v1 - vm)).reshape(shp)

    def _grad(self, v, x, y):
        shp = np.shape(x)
        x = np.asarray(x, float).ravel()
        y = np.asarray(y, float).ravel()
        n00, nm, n11, a, b, A = self._loc(x, y)
        g1 = (v[nm] - v[n00]) / self.h
        g2 = (v[n11] - v[nm]) / self.h
        return np.where(A, g1, g2).reshape(shp), np.where(A, g2, g1).reshape(shp)

    def c(self, x, y):
        return self._valor(self.cn, x, y)

    def phi(self, x, y):
        return self._valor(self.phin, x, y)

    def u(self, x, y):
        return self._valor(self.un, x, y)

    def grad_u(self, x, y):
        return self._grad(self.un, x, y)

    # -- limite superior de max phi num arco (o mesmo de CampoMalha.phi_max_arco), varios arcos de uma vez --
    def phi_max_arcos(self, X, Y, flecha=None):
        """X, Y: (n_arcos, n_amostras), as amostras de cada arco em ordem. Devolve (n_arcos,)."""
        na, ns = X.shape
        x = X.ravel()
        y = Y.ravel()
        n00, nm, n11, a, b, A = self._loc(x, y)
        v = self.phin
        v0, vm, v1 = v[n00], v[nm], v[n11]
        # f = c0 + cx x + cy y no triangulo de cada amostra
        g1 = (vm - v0) / self.h
        g2 = (v1 - vm) / self.h
        cx = np.where(A, g1, g2)
        cy = np.where(A, g2, g1)
        val = v0 + a * (vm - v0) + b * (v1 - vm)
        c0 = val - cx * x - cy * y
        cx, cy, c0 = cx.reshape(na, ns), cy.reshape(na, ns), c0.reshape(na, ns)
        out = val.reshape(na, ns).max(axis=1)
        if ns > 1:
            e1 = c0[:, :-1] + cx[:, :-1] * X[:, 1:] + cy[:, :-1] * Y[:, 1:]
            e2 = c0[:, 1:] + cx[:, 1:] * X[:, :-1] + cy[:, 1:] * Y[:, :-1]
            out = np.maximum(out, np.maximum(e1.max(axis=1), e2.max(axis=1)))
        if flecha is not None:
            out = out + flecha * np.sqrt((cx * cx + cy * cy).max(axis=1))
        return out

    def phi_max_arco(self, xs, ys, flecha=0.0) -> float:
        """Interface do cinematico.py (um arco): permite usar a marcha sequencial original com este campo."""
        X = np.asarray(xs, float).reshape(1, -1)
        Y = np.asarray(ys, float).reshape(1, -1)
        return float(self.phi_max_arcos(X, Y, np.array([flecha]) if flecha > 0.0 else None)[0])

    # -- a mesma triangulacao como CampoMalha (para a conferencia com o cinematico.py exato) --
    def triangulacao(self):
        nx, ny, h = self.nx, self.ny, self.h
        X, Y = np.meshgrid(self.x0 + h * np.arange(nx), h * np.arange(ny), indexing="ij")
        xy = np.c_[X.ravel(), Y.ravel()]
        idx = np.arange(nx * ny).reshape(nx, ny)
        jj = np.arange(ny - 1)[None, :]
        jm = self.jmin[:, None]
        cheia = jj >= jm + self.face[:, None]
        face = self.face[:, None] & (jj == jm)
        a, b, c, d = idx[:-1, :-1], idx[1:, :-1], idx[1:, 1:], idx[:-1, 1:]
        tri = np.r_[np.c_[a[cheia], b[cheia], c[cheia]], np.c_[a[cheia | face], c[cheia | face], d[cheia | face]]]
        usados = np.unique(tri)
        nova = -np.ones(nx * ny, int)
        nova[usados] = np.arange(len(usados))
        return xy[usados], nova[tri], usados

    def para_malha(self):
        xy, tri, usados = self.triangulacao()
        cm = kin.CampoMalha(xy, tri, self.cn[usados], self.phin[usados], self.un[usados], self.gamma_ef)
        return cm


# ---------------------------------------------------------------------------------------------------------------------
# Marcha da espiral por trechos, todos os arcos de uma vez
# ---------------------------------------------------------------------------------------------------------------------
class Estat:
    """Contadores (por processo)."""

    def __init__(self):
        self.zera()

    def zera(self):
        self.n = 0
        self.varreduras = 0
        self.fallback = 0
        self.max_varr = 0
        self.correcoes = 0
        self.fb_varr = 0


EST = Estat()


def _geometria(phis, t2, rB, l):
    """Arcos de B para A com os phi_i dados. Arco i: do lado de B (th_b, r_b) ao lado de A (th_a, r_a);
    r_a = r_b - l sin(phi_i) (o mesmo que r_b exp(-dtheta tan phi_i)) e dtheta = -ln(1 - l sin(phi_i)/r_b)/tan(phi_i)
    (cinematico._dtheta_arco_tras)."""
    sp, tp = np.sin(phis), np.tan(phis)
    r_a = rB - l * np.cumsum(sp)
    r_b = np.empty_like(r_a)
    r_b[0] = rB
    r_b[1:] = r_a[:-1]
    with np.errstate(divide="ignore", invalid="ignore"):
        q = l * sp / r_b
        dth = np.where(q < 1.0, -np.log1p(-np.minimum(q, 0.999999999999)) / tp, np.inf)
    th_a = t2 - np.cumsum(dth)
    th_b = np.empty_like(th_a)
    th_b[0] = t2
    th_b[1:] = th_a[:-1]
    return th_b, r_b, th_a, r_a, dth, tp


def _limites(campo, sl, th_b, r_b, dth, tp, r_a, Cx, Cy, frac, hs):
    """Limite superior de max phi nos arcos do intervalo sl (amostras em frac do arco; flecha h^2/(8 r_min))."""
    tb, rb = th_b[sl, None], r_b[sl, None]
    T = tb - frac[None, :] * dth[sl, None]
    R = rb * np.exp((T - tb) * tp[sl, None])
    with np.errstate(divide="ignore", invalid="ignore"):
        fl = hs * hs / (8.0 * r_a[sl])
    return campo.phi_max_arcos(Cx - R * np.cos(T), -Cy + R * np.sin(T), fl)


def construir_espiral_rapido(sB, t2, rB, campo: CampoGrade, talude, l_arco=None, n_amostras=6,
                             max_varr=80, tol=1e-11, quente=True, verificar=False) -> kin.Mecanismo:
    """A marcha do cinematico.construir_espiral, vetorizada (ver o cabecalho). Devolve o mesmo Mecanismo.

    Iteracao: phi_i <- limite(arco_i(phi)) so nos arcos a partir do primeiro que ainda mudou mais que tol (os
    anteriores e a sua geometria ja estao fixos: o arco i so depende de phi_0..phi_i). verificar=True refaz o limite
    em todos os arcos com a geometria final e sobe phi_i onde preciso (monotono)."""
    H = talude.H
    if l_arco is None:
        l_arco = H / 100.0
    m = kin.Mecanismo(sB, t2, rB)
    if not (-H < sB and 0.0 < t2 < math.pi and rB > 0.0):
        m.motivo = "parametros"
        return m
    if sB <= 0.0 and t2 >= math.pi - talude.beta:
        m.motivo = "theta2 >= pi - beta"
        return m
    xB, yB = kin.ponto_superficie(sB, talude)
    if xB > talude.x_dir:
        m.motivo = "B fora do dominio"
        return m
    Cx = xB + rB * math.cos(t2)
    Cy = rB * math.sin(t2) - yB
    m.xB, m.yB, m.Cx, m.Cy = xB, yB, Cx, Cy
    if Cy <= 0.0:
        m.motivo = "C abaixo da crista"
        return m
    EST.n += 1
    l = l_arco
    # numero de arcos: estimativa pela espiral com phi_ref, com folga; cresce se preciso
    nest = max(16, int(1.3 * _n_estimado(campo.phi_ref, t2, rB, Cy, l)) + 8)
    cache = getattr(campo, "_cache_phis", {}).get(l) if quente else None
    if cache is not None:
        # partida quente: os phi_i da ultima marcha com o mesmo campo e o mesmo l (Nelder-Mead anda pouco)
        phis = np.concatenate([cache, np.full(max(0, nest - len(cache)), cache[-1])])
    else:
        phis = np.full(nest, campo.phi_ref)
        phis[0] = campo.phi_max_arco([Cx - rB * math.cos(t2)], [-Cy + rB * math.sin(t2)])
    with np.errstate(invalid="ignore", over="ignore"):
        return _marcha(m, sB, t2, rB, campo, talude, l_arco, n_amostras, max_varr, tol, quente, verificar,
                       Cx, Cy, phis)


def _marcha(m, sB, t2, rB, campo, talude, l_arco, n_amostras, max_varr, tol, quente, verificar, Cx, Cy, phis,
            n_jacobi=15):
    l = l_arco
    frac = np.linspace(0.0, 1.0, n_amostras)
    hs = l / (n_amostras - 1)
    ysup = talude.y_sup
    kc_ant = -1
    i0 = 0
    for varr in range(max_varr):
        th_b, r_b, th_a, r_a, dth, tp = _geometria(phis, t2, rB, l)
        xa, ya = Cx - r_a * np.cos(th_a), -Cy + r_a * np.sin(th_a)
        idx = np.flatnonzero(~(ya >= ysup(xa)))   # lado de A do arco acima da superficie (ou indefinido)
        if len(idx) == 0:
            if len(phis) > 20000:
                m.motivo = "arcos demais"
                return m
            phis = np.concatenate([phis, np.full(len(phis), phis[-1])])
            continue
        kc = idx[0]
        n = kc + 1
        if i0 >= n:
            if kc == kc_ant:
                break
            kc_ant = kc
            i0 = min(i0, n)
            continue
        sl = slice(i0, n)
        ok = np.isfinite(dth[sl])
        b = phis[sl].copy()
        if np.all(ok):
            b = _limites(campo, sl, th_b, r_b, dth, tp, r_a, Cx, Cy, frac, hs)
        elif np.any(ok):
            j = i0 + np.flatnonzero(ok)
            b[ok] = _limites(campo, j, th_b, r_b, dth, tp, r_a, Cx, Cy, frac, hs)
        if varr >= n_jacobi:
            # ciclo de periodo 2 num arco (uma amostra troca de triangulo quando phi_i muda: limite(phi_alto) =
            # phi_baixo e limite(phi_baixo) = phi_alto): passa a subir so (phi_i <- max(phi_i, limite)), como a
            # iteracao do cinematico.construir_espiral, que para quando o limite nao passa de phi_i
            b = np.maximum(b, phis[sl])
        mud = np.flatnonzero(~(np.abs(b - phis[sl]) < tol))
        phis[sl] = b
        phis[n:] = b[-1]
        if len(mud) == 0 and kc == kc_ant:
            break
        i0 = i0 + mud[0] if len(mud) else n
        kc_ant = kc
    else:
        EST.fallback += 1
        EST.fb_varr += 1
        return kin.construir_espiral(sB, t2, rB, campo, talude, l_arco, n_amostras)
    EST.varreduras += varr + 1
    EST.max_varr = max(EST.max_varr, varr + 1)
    if verificar:
        # conferencia da normalidade com a geometria final em todos os arcos
        for _ in range(6):
            th_b, r_b, th_a, r_a, dth, tp = _geometria(phis, t2, rB, l)
            xa, ya = Cx - r_a * np.cos(th_a), -Cy + r_a * np.sin(th_a)
            n = np.flatnonzero(~(ya >= ysup(xa)))[0] + 1
            if not np.all(np.isfinite(dth[:n])):
                break
            bfin = _limites(campo, slice(0, n), th_b, r_b, dth, tp, r_a, Cx, Cy, frac, hs)
            if not np.any(bfin > phis[:n] + 1e-12):
                break
            EST.correcoes += 1
            phis[:n] = np.maximum(phis[:n], bfin + 1e-12)
            phis[n:] = phis[n - 1]
        else:
            EST.fallback += 1
            return kin.construir_espiral(sB, t2, rB, campo, talude, l_arco, n_amostras)
    else:
        th_b, r_b, th_a, r_a, dth, tp = _geometria(phis, t2, rB, l)
        xa, ya = Cx - r_a * np.cos(th_a), -Cy + r_a * np.sin(th_a)
    if quente:
        if not hasattr(campo, "_cache_phis"):
            campo._cache_phis = {}
        campo._cache_phis[l] = phis[:n].copy()
    # ---- as mesmas regras de parada do cinematico.construir_espiral ----
    for i in range(n - 1):   # arcos inteiros antes do corte
        if not np.isfinite(dth[i]) or th_a[i] <= 0.0:
            m.motivo = "nao chega a crista"
            return m
        if ya[i] > talude.y_base or xa[i] < talude.x_esq or xa[i] > talude.x_dir:
            m.motivo = "sai do dominio"
            return m
    i = n - 1
    if not np.isfinite(dth[i]) or th_a[i] <= 0.0:
        m.motivo = "nao chega a crista"
        return m
    if n == 1:
        m.motivo = "nao entra no solo"
        return m
    th, r, tpi = th_b[i], r_b[i], tp[i]
    ysc = talude.y_sup_escalar

    def g(t):
        rr = r * math.exp((t - th) * tpi)
        x, y = Cx - rr * math.cos(t), -Cy + rr * math.sin(t)
        return y - ysc(x)
    t_ini = optimize.brentq(g, th_a[i], th, xtol=1e-13)
    r_ini = r * math.exp((t_ini - th) * tpi)
    xe = Cx - r_ini * math.cos(t_ini)
    ths = np.r_[t2, th_a[:i]]
    rs = np.r_[rB, r_a[:i]]
    tps = tp[:i]
    if t_ini < th:
        ths, rs, tps = np.r_[ths, t_ini], np.r_[rs, r_ini], np.r_[tps, tpi]
    if xe >= 0.0:
        m.motivo = "A fora da crista"
        return m
    if xe < talude.x_esq:
        m.motivo = "sai do dominio"
        return m
    m.th, m.r, m.tp = ths[::-1].copy(), rs[::-1].copy(), tps[::-1].copy()
    m.xA = xe
    m.classe = "I" if sB <= 0.0 else "II"
    return m


def _n_estimado(phi, t2, rB, Cy, l):
    """Numero de arcos da espiral classica (phi uniforme) de B ate a altura da crista (y = 0)."""
    tp = math.tan(phi)
    # r(theta) = rB exp((theta - t2) tp); a crista: r sin(theta) = Cy; percorre theta de t2 para tras
    th = np.linspace(t2, 1e-3, 400)
    r = rB * np.exp((th - t2) * tp)
    k = np.flatnonzero(r * np.sin(th) <= Cy)
    th1 = th[k[0]] if len(k) else th[-1]
    comp = rB * (1.0 - math.exp((th1 - t2) * tp)) / math.sin(phi) if tp > 1e-12 else rB * (t2 - th1)
    return comp / l


def potencias_contorno(m, campo, talude, ng_theta=3, dissipacao="exata"):
    """O mesmo que cinematico.potencias (dissipacao arco a arco, normalidade nos pontos de Gauss, P_gamma exato em
    rho), mas com P_u SO pela forma de contorno (div U = 0: P_u = soma nos extremos dos intervalos de
    +-u rho drho/dtheta dtheta), que o cinematico.potencias calcula como verificacao. Para u P1 continuo as duas
    formas sao iguais; a diferenca e so de quadratura (cerca de 1e-4 em Gamma). Evita os ~50 mil pontos de grad u da
    forma de volume: e o que a busca usa. O Gamma final da busca e recalculado com cinematico.potencias (volume)."""
    H = talude.H
    Cx, Cy = m.Cx, m.Cy
    th, r, tp = m.th, m.r, m.tp
    xi, wi = kin._GL[2]
    a, b = th[:-1], th[1:]
    h = 0.5 * (b - a)
    tg = (0.5 * (a + b))[:, None] + h[:, None] * xi[None, :]
    rg = r[:-1, None] * np.exp((tg - a[:, None]) * tp[:, None])
    xg, yg = Cx - rg * np.cos(tg), -Cy + rg * np.sin(tg)
    tpc = np.maximum(tp, 1e-300)
    I_r2 = np.where(tp < 1e-12, r[:-1] ** 2 * (b - a), r[:-1] ** 2 * np.expm1(2 * (b - a) * tpc) / (2 * tpc))
    cg = campo.c(xg, yg)
    phig = campo.phi(xg, yg)
    m.viol_normalidade = float(max(0.0, np.degrees(np.max(phig - np.arctan(tp)[:, None]))))
    if m.viol_normalidade > kin.TOL_NORMALIDADE:
        m.motivo = "normalidade violada"
        m.Gamma = INF
        return m
    if dissipacao == "exata":
        tl = np.tan(phig)
        dens = cg * np.where(tl > 1e-12, tp[:, None] / np.maximum(tl, 1e-300), 1.0)
    else:
        dens = cg
    m.Pmr = float(np.sum(I_r2 * (dens * wi[None, :]).sum(axis=1) / wi.sum()))
    thO = math.atan2(Cy, Cx)
    angs = [th[0], th[-1], thO]
    if m.classe == "II":
        angs.append(math.atan2(Cy + H, Cx - talude.xT))
    lo, hi = min(angs), max(angs)
    cortes = np.unique(np.concatenate([th, angs]))
    cortes = cortes[(cortes >= lo) & (cortes <= hi)]
    dmax = float(np.max(np.diff(th)))
    nsub = np.maximum(1, np.ceil(np.diff(cortes) / dmax).astype(int))
    if np.any(nsub > 1):
        cortes = np.concatenate([np.linspace(cortes[i], cortes[i + 1], nsub[i] + 1)[:-1]
                                 for i in range(len(nsub))] + [cortes[-1:]])
    tq, wq = kin._gauss_intervalos(cortes[:-1], cortes[1:], ng_theta)
    tq, wq = tq.ravel(), wq.ravel()
    ray, ra, rb, da, db = kin._cruzamentos(tq, m, talude)
    ti, wt = tq[ray], wq[ray]
    ci, si = np.cos(ti), np.sin(ti)
    m.Pg = float(campo.gamma_ef * np.sum(wt * ci * (rb ** 3 - ra ** 3) / 3.0))
    if campo.tem_u:
        ub = campo.u(Cx - rb * ci, -Cy + rb * si)
        ua = campo.u(Cx - ra * ci, -Cy + ra * si)
        m.Pu = m.Pu_contorno = float(np.sum(wt * (ub * rb * db - ua * ra * da)))
    else:
        m.Pu = m.Pu_contorno = 0.0
    Pext = m.Pg + m.Pu
    m.Gamma = m.Pmr / Pext if Pext > 0.0 else INF
    return m


def avaliar_rapido(sB, t2, rB, campo, talude, l_arco=None, pu="volume", verificar=False, **kw) -> kin.Mecanismo:
    """Mecanismo rapido + potencias: pu = 'volume' (cinematico.potencias, o Gamma de referencia) ou 'contorno'
    (potencias_contorno, para a busca)."""
    m = construir_espiral_rapido(sB, t2, rB, campo, talude, l_arco, verificar=verificar)
    if m.classe == "":
        return m
    if pu == "contorno":
        return potencias_contorno(m, campo, talude, dissipacao=kw.get("dissipacao", "exata"))
    return kin.potencias(m, campo, talude, **kw)


# ---------------------------------------------------------------------------------------------------------------------
# Busca rapida
# ---------------------------------------------------------------------------------------------------------------------
def minimizar_rapido(campo, talude, partidas=None, classes=("I", "II"), grade=(5, 5, 4), l_arco=None,
                     l_grade=None, n_nm=1, folga_2a=0.03, xatol=1e-4, fatol=1e-5, maxfev=150,
                     dissipacao="exata", pu_busca="contorno") -> kin.Resultado:
    """Gamma minimo (Eq. 55) com poucas avaliacoes.

    partidas: dict classe -> (sB, theta2, rB) (ex.: o mecanismo critico deterministico de cada classe).
    Para cada classe: avalia a partida quente e a grade pequena de mecanismos classicos com arcos grossos
    (l_grade = H/25); Nelder-Mead (arcos H/100) a partir da melhor; tambem a partir da segunda melhor de outra bacia
    se o seu Gamma grosso estiver a menos de folga_2a da melhor."""
    H = talude.H
    l_arco = H / 100.0 if l_arco is None else l_arco
    l_grade = H / 25.0 if l_grade is None else l_grade
    t0 = time.perf_counter()
    nav = [0]
    nviol = {"normalidade violada": 0, "phi_i nao convergiu": 0}
    PEN = 1e3

    def G(z, cl, l):
        nav[0] += 1
        m = avaliar_rapido(*kin._z_para_param(z, cl, talude), campo, talude, l, pu=pu_busca,
                           dissipacao=dissipacao)
        if m.motivo in nviol:
            nviol[m.motivo] += 1
        return m.Gamma if m.valido else PEN

    por = {}
    for cl in classes:
        cands = []
        if partidas and cl in partidas:
            cands.append(kin._param_para_z(*partidas[cl], talude))
        cands += [kin._param_para_z(*p, talude) for p in kin._candidatos(talude, campo.phi_ref, cl, *grade)]
        vals = sorted(((G(z, cl, l_grade), z) for z in cands), key=lambda v: v[0])
        if not vals or vals[0][0] >= PEN:
            continue
        starts = [vals[0]]
        for g0, z0 in vals[1:]:
            if g0 >= PEN or g0 > vals[0][0] * (1.0 + folga_2a) or len(starts) > n_nm:
                break
            if all(np.max(np.abs(z0 - s[1])) > 0.1 for s in starts):
                starts.append((g0, z0))
        best = None
        for g0, z0 in starts:
            simplex = np.vstack([z0, z0 + np.diag([0.04, 0.04, 0.04])])
            r = optimize.minimize(G, z0, args=(cl, l_arco), method="Nelder-Mead",
                                  options=dict(xatol=xatol, fatol=fatol, maxfev=maxfev, initial_simplex=simplex))
            if best is None or r.fun < best.fun:
                best = r
        if best is not None and best.fun < PEN:
            por[cl] = avaliar_rapido(*kin._z_para_param(best.x, cl, talude), campo, talude, l_arco,
                                     verificar=True, dissipacao=dissipacao, contorno=True)
    if not por:
        raise RuntimeError("nenhum mecanismo valido")
    mbest = min(por.values(), key=lambda mm: mm.Gamma)
    return kin.Resultado(mbest.Gamma, mbest, por, nav[0], time.perf_counter() - t0,
                         nviol["normalidade violada"], nviol["phi_i nao convergiu"])


def partidas_deterministicas(problema, arquivo=None, verbose=False):
    """Mecanismo critico de cada classe com os campos medios (cinematico.minimizar, o metodo exato, no CampoGrade
    dos campos medios). Guarda em arquivo (json) para nao refazer. Devolve (dict classe -> (sB, theta2, rB), Gamma)."""
    import json
    import os
    if arquivo and os.path.exists(arquivo):
        with open(arquivo) as f:
            d = json.load(f)
        return {k: tuple(v) for k, v in d["partidas"].items()}, d["Gamma"]
    d = problema.campos_medios()
    cg = CampoGrade(problema.xg[0], problema.h, d["c"], d["phi"], d["u"], problema.gamma_ef, problema.tal)
    r = kin.minimizar(cg, problema.tal, verbose=verbose)
    part = {cl: (m.sB, m.t2, m.rB) for cl, m in r.por_classe.items()}
    if arquivo:
        with open(arquivo, "w") as f:
            json.dump({"partidas": part, "Gamma": r.Gamma,
                       "por_classe": {cl: m.Gamma for cl, m in r.por_classe.items()}}, f, indent=1)
    return part, r.Gamma
