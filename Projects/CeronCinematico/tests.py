"""Testes da abordagem cinematica.  Rodar:  python3 -I tests.py

(a) espiral por trechos (codigo geral) = forma fechada quando c, phi sao
    uniformes, mecanismo a mecanismo e no minimo; mecanismo com C atras da
    reta da face (forma fechada x codigo geral x poligono independente);
(b) taludes secos de Cho (Secao 5.3 do artigo): Gamma = 1,354 (coesivo) e
    Gamma = 1,777 / FS = 1,203 (c-phi);
(c) rebaixamento deterministico com as poropressoes do TRANSIENTE NeoPZ
    (problema diferente do artigo: sem referencia do artigo) x GIM da mesma
    exportacao extrapolado; e, a parte, a percolacao PERMANENTE (Eqs. 19-21)
    x o 1,336 do artigo (informativo: o contorno do MEF do artigo nao e dado);
(d) uma realizacao heterogenea exportada (avaliar_realizacao) x GIM
    extrapolado, e a normalidade phi(x) <= phi_i conferida arco a arco;
(e) P_u: equivalencia hidrostatica (u = -gamma_w y) e o campo semianalitico
    v'_opt do artigo (Secao 4.1) x Figs. 8 e 9.
Tambem gera figs/mecanismo_det.pdf e figs/mecanismo_real0.pdf."""
import math
import os
import sys
import time

import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
import cinematico as k              # noqa: E402
import permanente as pm             # noqa: E402
import percolacao_analitica as pa   # noqa: E402

DADOS = os.path.join(AQUI, "dados")
linhas, falhas = [], []


def reg(teste, descricao, obtido, referencia, criterio, ok):
    linhas.append((teste, descricao, obtido, referencia, criterio, "OK" if ok else "FALHA"))
    if not ok:
        falhas.append(f"{teste}: {descricao}")


def rel(a, b):
    return abs(a - b) / abs(b)


def gim_niveis(pasta):
    """Niveis de refinamento do GIM. Usa linha_mc_nref5.csv (mesma
    realizacao rodada com NREF=5) se existir, senao linha_mc.csv (NREF=3)."""
    for nome in ("linha_mc_nref5.csv", "linha_mc.csv"):
        arq = os.path.join(pasta, nome)
        if os.path.exists(arq):
            with open(arq) as f:
                d = dict(zip(f.readline().strip().split(","), f.readline().strip().split(",")))
            return np.array([float(v) for v in d["niveis"].split(";")])
    raise FileNotFoundError(pasta)


def richardson(v):
    """Limite pela razao dos dois ultimos incrementos (convergencia
    geometrica): v_n - d_n r/(1 - r), r = d_n/d_(n-1). Devolve tambem o
    mesmo limite com um nivel a menos (estimativa da incerteza)."""
    def lim(w):
        d = -np.diff(w)
        r = d[-1] / d[-2]
        return w[-1] - d[-1] * r / (1 - r), r
    g, r = lim(v)
    g0, _ = lim(v[:-1])
    return g, r, g0


def normalidade_densa(m, campo, n=200):
    """max (phi(x) - phi_i) em graus com n pontos por arco, e quantos arcos
    passam de 1e-6 graus (conferencia independente de phi_max_arco)."""
    pior, nviol = -1e9, 0
    for i in range(len(m.tp)):
        t = np.linspace(m.th[i], m.th[i + 1], n)
        rr = m.r[i] * np.exp((t - m.th[i]) * m.tp[i])
        ph = campo.phi(m.Cx - rr * np.cos(t), -m.Cy + rr * np.sin(t))
        d = float(np.degrees(ph.max() - math.atan(m.tp[i])))
        pior = max(pior, d)
        nviol += d > 1e-6
    return pior, nviol


def pg_poligono(m, tal, gamma_ef):
    """P_gamma por um poligono denso do bloco (formula do laco), independente
    da quadratura polar."""
    X, Y = m.pontos(200)
    pts = list(zip(X, Y))
    if m.classe == "II":
        pts.append((tal.xT, tal.H))
    pts.append((0.0, 0.0))
    P = np.array(pts)
    x, y = P[:, 0], P[:, 1]
    x1, y1 = np.roll(x, -1), np.roll(y, -1)
    cr = x * y1 - x1 * y
    A, Mx = cr.sum() / 2, ((x + x1) * cr).sum() / 6
    return gamma_ef * np.sign(A) * (m.Cx * A - Mx)


# ---------------------------------------------------------------- (a)
def teste_a():
    phi = math.radians(30.0)
    tal = k.Talude(5.0, math.radians(45.0), -35.0, 35.0, 35.0)
    unif = k.CampoUniforme(10.0, phi, 10.19)
    # CampoMalha com c, phi constantes e u = 0: exercita o caminho
    # heterogeneo (phi_max_arco, Gauss nos arcos)
    _, cdet = k.ler_exportacao(os.path.join(DADOS, "det"))
    n = len(cdet.xy)
    malha = k.CampoMalha(cdet.xy, cdet.tri, np.full(n, 10.0), np.full(n, phi), np.zeros(n), 10.19)
    casos = [("I", 0.60, 1.90, 0.80), ("I", 0.69, 1.96, 1.00), ("I", 1.10, 1.60, 0.50),
             ("I", 0.20, 2.10, 1.00), ("II", 0.30, 2.30, 2.00), ("II", 0.60, 2.20, 2.00)]
    pior = 0.0
    for cl, t1, t2, p3 in casos:
        mc = k.mecanismo_classico(cl, t1, t2, p3, tal, 10.0, phi, 10.19)
        assert mc.valido, (cl, t1, t2, p3)
        for nome, campo in (("uniforme", unif), ("malha", malha)):
            m = k.avaliar(*k.classico_para_geral(mc), campo, tal)
            d = max(rel(m.Pmr, mc.Pmr), rel(m.Pg, mc.Pg), rel(m.Gamma, mc.Gamma),
                    abs(m.xA - mc.xA) / tal.H)
            pior = max(pior, d)
            reg("a", f"mec {cl} ({t1},{t2},{p3}) campo {nome}: Gamma", m.Gamma, mc.Gamma,
                "dif rel < 1e-3", d < 1e-3)
    # C atras da reta da face (Cy + tan(beta) Cx < 0), C acima da crista:
    # mecanismo admissivel, antes recusado ('C dentro do solo')
    ph2 = math.radians(0.9)
    mc = k.mecanismo_classico("II", 0.075, 2.661, 1.019, tal, 10.0, ph2, 10.19)
    m = k.avaliar(*k.classico_para_geral(mc), k.CampoUniforme(10.0, ph2, 10.19), tal)
    pp = pg_poligono(m, tal, 10.19)
    reg("a", f"C atras da face (C=({m.Cx:.2f},{m.Cy:.2f} acima)), mec II: forma fechada x geral",
        mc.Gamma, m.Gamma, "valido; dif rel < 1e-4",
        mc.valido and m.valido and rel(mc.Gamma, m.Gamma) < 1e-4)
    reg("a", "  idem: P_gamma geral x poligono independente", m.Pg, pp, "dif rel < 1e-4",
        rel(m.Pg, pp) < 1e-4)
    # minimo: fechado (evolucao diferencial) x geral (grade + Nelder-Mead)
    for nome, talm, c, ph, g in (
            ("rebaix. seco H=5 b=45 c=10 phi=30 g'=10,19", tal, 10.0, 30.0, 10.19),
            ("Cho c-phi H=10 b=45 c=10 phi=30 g=20", k.Talude(10.0, math.radians(45), -10, 20, 15),
             10.0, 30.0, 20.0)):
        rc = k.minimizar_classico(talm, c, math.radians(ph), g)
        gc = min(m.Gamma for m in rc.values())
        rg = k.minimizar(k.CampoUniforme(c, math.radians(ph), g), talm)
        reg("a", f"minimo {nome}", rg.Gamma, gc, "dif rel < 1e-3", rel(rg.Gamma, gc) < 1e-3)
    return pior


# ---------------------------------------------------------------- (b)
def teste_b():
    # Secao 5.3.1: coesivo, H=5, 2:1, cu=23, gamma=20 total, base 10 m abaixo da crista
    tal = k.Talude(5.0, math.atan(0.5), -10.0, 20.0, 10.0)
    r = k.minimizar_classico(tal, 23.0, 0.0, 20.0)
    m = min(r.values(), key=lambda mm: mm.Gamma)
    reg("b", f"Cho coesivo (cu=23): Gamma, mec {m.classe} D={m.p3:.2f} m", m.Gamma, 1.354,
        "|dif| < 2e-3 (Cho LEM 1,356)", abs(m.Gamma - 1.354) < 2e-3)
    rg = k.minimizar(k.CampoUniforme(23.0, 0.0, 20.0), tal)
    reg("b", "Cho coesivo, codigo geral (phi = 0)", rg.Gamma, 1.354, "|dif| < 2e-3",
        abs(rg.Gamma - 1.354) < 2e-3)
    # Secao 5.3.2: c-phi, 1:1, H = 10 m (Fig. 19; o texto diz 5 m)
    tal = k.Talude(10.0, math.radians(45.0), -10.0, 20.0, 15.0)
    r = k.minimizar_classico(tal, 10.0, math.radians(30.0), 20.0)
    g = min(mm.Gamma for mm in r.values())
    reg("b", "Cho c-phi (H=10, 1:1): Gamma", g, 1.777, "|dif| < 2e-3", abs(g - 1.777) < 2e-3)
    fs = k.fs_classico(tal, 10.0, math.radians(30.0), 20.0, classes=("I",), Fmin=1.0, Fmax=1.6)
    reg("b", "Cho c-phi: FS (Eqs. 57-58)", fs, 1.203, "|dif| < 2e-3 (Cho 1,204)",
        abs(fs - 1.203) < 2e-3)


# ---------------------------------------------------------------- (c)
def teste_c():
    pasta = os.path.join(DADOS, "det")
    tal, campo = k.ler_exportacao(pasta)
    t0 = time.perf_counter()
    res = k.minimizar(campo, tal)
    t_min = time.perf_counter() - t0
    res.talude, res.campo = tal, campo
    m = res.mecanismo
    niv = gim_niveis(pasta)
    g_ext, r_ext, g_ext0 = richardson(niv)
    reg("c", f"TRANSIENTE (t=0,1 h, nao saturado; nao e o problema do artigo): Gamma cin. "
             f"(mec {m.classe}, xA={m.xA:.2f}, B={'pe' if abs(m.xB - tal.xT) < 1e-6 else f'{m.xB:.2f}'})",
        res.Gamma, g_ext, "ref = GIM extrapolado", True)
    reg("c", f"  GIM NeoPZ, mesma exportacao, NREF={len(niv) - 1} (niveis "
             f"{';'.join(f'{v:.4f}' for v in niv)})", niv[-1], g_ext, "informativo", True)
    reg("c", f"  GIM extrapolado (Richardson, r={r_ext:.2f}; com um nivel a menos {g_ext0:.4f})",
        g_ext, g_ext0, "|dif| < 1% (incerteza da extrapolacao)", rel(g_ext, g_ext0) < 0.01)
    reg("c", "  |cinematico - GIM extrap.|/GIM extrap.", rel(res.Gamma, g_ext), 0.0,
        "< 1% (campo homogeneo: o log-espiral deve ser quase o colapso)",
        rel(res.Gamma, g_ext) < 0.01)
    dP = rel(m.Pu, m.Pu_contorno)
    reg("c", "  P_u volume x P_u contorno (div U = 0)", m.Pu, m.Pu_contorno, "dif rel < 5e-3",
        dP < 5e-3)
    gs = [k.avaliar(m.sB, m.t2, m.rB, campo, tal, l).Gamma for l in (tal.H / 50, tal.H / 200)]
    reg("c", "  convergencia: Gamma com arcos H/50 e H/200", gs[0], gs[1], "dif rel < 2e-3",
        rel(gs[0], gs[1]) < 2e-3)
    # percolacao PERMANENTE saturada (Eqs. 19-21) na mesma malha: e o
    # problema do artigo, mas o dominio e o contorno do MEF do artigo nao
    # sao dados, entao a comparacao com 1,336 e so informativa
    gp = {}
    for lat in ("dirichlet", "neumann"):
        perm = pm.campo_permanente(campo, tal, lateral=lat)
        gp[lat] = k.minimizar(perm, tal).Gamma
        if lat == "dirichlet":
            du = np.abs(perm.un - campo.un).max()
    reg("c", "PERMANENTE (Eqs. 19-21), lateral esq. u=0: Gamma x artigo (Fig. 24, MEF)",
        gp["dirichlet"], 1.336, "informativo: contorno do MEF do artigo desconhecido", True)
    reg("c", "  PERMANENTE, lateral esq. sem fluxo", gp["neumann"], 1.336, "informativo", True)
    # tempo por avaliacao de Gamma
    t0 = time.perf_counter()
    for _ in range(50):
        k.avaliar(m.sB, m.t2, m.rB, campo, tal, contorno=False)
    t_av = (time.perf_counter() - t0) / 50
    return res, t_min, t_av, du


# ---------------------------------------------------------------- (d)
def teste_d():
    pasta = os.path.join(DADOS, "real0")
    t0 = time.perf_counter()
    res = k.avaliar_realizacao(pasta)
    t_min = time.perf_counter() - t0
    niv = gim_niveis(pasta)
    g_ext, r_ext, g_ext0 = richardson(niv)
    m = res.mecanismo
    phis = np.degrees(np.arctan(m.tp))
    reg("d", f"realizacao 0: Gamma cinematico ('exata'; phi_i de {phis.min():.1f} a "
             f"{phis.max():.1f} graus)", res.Gamma, g_ext, "informativo; ref = GIM extrapolado", True)
    rs = k.minimizar(res.campo, res.talude, dissipacao="simples")
    reg("d", "  idem com dP_mr = c r^2 dtheta ('simples')", rs.Gamma, res.Gamma,
        "dif rel < 1% (phi_i ~ phi(x))", rel(rs.Gamma, res.Gamma) < 0.01)
    reg("d", f"  GIM NeoPZ NREF={len(niv) - 1} (niveis {';'.join(f'{v:.4f}' for v in niv)})",
        niv[-1], g_ext, "informativo", True)
    reg("d", f"  GIM extrapolado (r={r_ext:.2f}; com um nivel a menos {g_ext0:.4f})", g_ext, g_ext0,
        "|dif| < 1%", rel(g_ext, g_ext0) < 0.01)
    reg("d", "  cinematico / GIM extrapolado", res.Gamma / g_ext, 1.0,
        ">= 0,99 (limite superior) e < 1,10", 0.99 <= res.Gamma / g_ext < 1.10)
    # normalidade: mecanismo critico e 400 mecanismos sorteados, 200 pontos por arco
    pior, nv = normalidade_densa(m, res.campo)
    reg("d", f"  normalidade no critico: max(phi(x) - phi_i), 200 pts/arco ({len(m.tp)} arcos)",
        pior, 0.0, "<= 1e-6 graus", pior <= 1e-6)
    rng = np.random.default_rng(7)
    nval = narcos = nviol = 0
    piorr = -1e9
    tal = res.talude
    while nval < 400:
        sB = rng.uniform(-0.9, 1.5) * tal.H
        mm = k.construir_espiral(sB, rng.uniform(0.3, math.pi - 0.05), tal.H * math.exp(rng.uniform(-0.5, 1.5)),
                                 res.campo, tal)
        if mm.classe == "":
            continue
        nval += 1
        p, v = normalidade_densa(mm, res.campo)
        piorr, nviol, narcos = max(piorr, p), nviol + v, narcos + len(mm.tp)
    reg("d", f"  normalidade em 400 mecanismos sorteados ({narcos} arcos): arcos violados", nviol, 0,
        f"0 (pior {piorr:.1e} graus)", nviol == 0)
    reg("d", "  avaliacoes descartadas na minimizacao (normalidade; phi_i sem conv.)",
        res.n_normalidade + res.n_phi_nao_conv, 0, "0", res.n_normalidade + res.n_phi_nao_conv == 0)
    t0 = time.perf_counter()
    for _ in range(50):
        k.avaliar(m.sB, m.t2, m.rB, res.campo, res.talude, contorno=False)
    t_av = (time.perf_counter() - t0) / 50
    return res, t_min, t_av


# ---------------------------------------------------------------- (e)
def teste_e():
    # 1) equivalencia hidrostatica: u = -gamma_w y (p_w = 0 em todo o dominio)
    #    da -grad u = gamma_w e_y, ou seja, peso gamma' + gamma_w = gamma:
    #    Gamma = Gamma_seco(gamma') gamma'/gamma. Trava sinal e escala de P_u.
    tal, cdet = k.ler_exportacao(os.path.join(DADOS, "det"))
    n = len(cdet.xy)
    phi = math.radians(30.0)
    hid = k.CampoMalha(cdet.xy, cdet.tri, np.full(n, 10.0), np.full(n, phi),
                       -k.GAMMA_W * cdet.xy[:, 1], 10.19)
    gh = k.minimizar(hid, tal).Gamma
    gseco = min(m.Gamma for m in k.minimizar_classico(tal, 10.0, phi, 20.0).values())
    reg("e", "P_u hidrostatico (u = -gamma_w y, gamma'=10,19) x seco com gamma = 20", gh, gseco,
        "dif rel < 1e-3", rel(gh, gseco) < 1e-3)
    # 2) campo semianalitico v'_opt (Secao 4.1), alfa = 1, h_w = H
    out = {}
    H = 5.0
    for bd, ref in ((45, 1.95), (60, 1.50), (75, 1.05)):
        b = math.radians(bd)
        f = pa.CampoVopt(10.0, math.radians(30.0), 10.19, H, b)
        g = k.minimizar(f, k.Talude(H, b, -60.0, 80.0, 40.0)).Gamma
        out[bd] = g
        reg("e", f"v'_opt, Fig. 9 tracejada (alfa=1, H=5, c=10, phi=30, g=20), beta={bd}: Gamma", g, ref,
            "|dif| < 6% (leitura do grafico)", rel(g, ref) < 0.06)
    # Fig. 8, h_w/H = 1 (gamma = 18, gamma' = 8,19): H_crit = H Gamma(H).
    # Os paineis da Fig. 8 estao trocados em relacao a Tabela 1: com
    # h_w = 0 (sem percolacao) os parametros 'London' (c=6, phi=32) dao
    # H_crit = 13,0 a 60 graus e 'Israeli' (c=11,7, phi=24,7) dao 18,0,
    # e a figura mostra 18 no painel London e 13 no Israeli.
    Hc = 10.0
    for nome, c, ph, bd, painel, ref in (("Londres", 6.0, 32.0, 60, "Israeli", 4.84),
                                         ("Londres", 6.0, 32.0, 35, "Israeli", 8.6),
                                         ("Israel", 11.7, 24.7, 60, "London", 8.05),
                                         ("Israel", 11.7, 24.7, 30, "London", 14.7)):
        b = math.radians(bd)
        f = pa.CampoVopt(c, math.radians(ph), 18.0 - k.GAMMA_W, Hc, b)
        hcrit = Hc * k.minimizar(f, k.Talude(Hc, b, -200.0, 300.0, 200.0)).Gamma
        reg("e", f"v'_opt, Fig. 8 (h_w/H=1): H_crit, param. {nome} (Tab. 1), beta={bd} "
                 f"(painel {painel})", hcrit, ref, "|dif| < 4% (leitura em escala log)",
            rel(hcrit, ref) < 0.04)
    f = pa.CampoVopt(10.0, math.radians(30.0), 10.19, 5.0, math.radians(45.0))
    xs = np.random.default_rng(3).uniform(-6, 3, 4000)
    ys = np.random.default_rng(4).uniform(0.05, 4.9, 4000)
    ok = (ys > np.maximum(0, xs)) & (np.hypot(xs, ys) > 0.1) & (np.abs(np.hypot(xs, ys) - f.Rw) > 1e-3)
    reg("e", "  v'_opt: max |div v/k| (diferencas finitas), beta=45", float(np.abs(f.div_v(xs[ok], ys[ok])).max()),
        0.0, "< 1e-5", float(np.abs(f.div_v(xs[ok], ys[ok])).max()) < 1e-5)
    reg("e", "  v'_opt: -J*/(k gamma_w^2 H^2), beta=45 (Fig. 5, alfa=1, curva cheia ~0,47)",
        -f.J_normalizado(), 0.47, "informativo", True)
    return out


if __name__ == "__main__":
    T0 = time.perf_counter()
    pior = teste_a()
    teste_b()
    rc, tmin_c, tav_c, du = teste_c()
    rd, tmin_d, tav_d = teste_d()
    teste_e()
    # tempo com campo uniforme
    tal = k.Talude(5.0, math.radians(45), -35, 35, 35)
    mc = k.mecanismo_classico("I", 0.69, 1.96, 1.0, tal, 10.0, math.radians(30), 10.19)
    t0 = time.perf_counter()
    for _ in range(200):
        k.avaliar(*k.classico_para_geral(mc), k.CampoUniforme(10.0, math.radians(30), 10.19), tal)
    tav_u = (time.perf_counter() - t0) / 200

    print()
    print(f"{'teste':5} | {'descricao':90} | {'obtido':>10} | {'referencia':>10} | criterio | ok")
    print("-" * 160)
    for t, d, o, r, c, ok in linhas:
        print(f"{t:5} | {d:90.90} | {o:10.4f} | {r:10.4f} | {c} | {ok}")
    print("-" * 160)
    print(f"(a) maior diferenca relativa mecanismo a mecanismo: {pior:.2e}")
    print(f"(c) max |u_permanente - u_transiente| na malha: {du:.1f} kPa")
    print("tempo por avaliacao de Gamma (arcos H/100, P_u por volume):")
    print(f"    campo uniforme, sem agua: {1e3 * tav_u:6.1f} ms")
    print(f"    det. (u da malha)       : {1e3 * tav_c:6.1f} ms   minimizacao: {tmin_c:5.1f} s "
          f"({rc.n_aval} avaliacoes)")
    print(f"    realizacao 0            : {1e3 * tav_d:6.1f} ms   minimizacao: {tmin_d:5.1f} s "
          f"({rd.n_aval} avaliacoes)")
    try:
        import figura
        os.makedirs(os.path.join(AQUI, "figs"), exist_ok=True)
        figura.figura(os.path.join(DADOS, "det"), os.path.join(AQUI, "figs", "mecanismo_det.pdf"),
                      "u", "Rebaixamento total (h_w = H), campos medios, p_w do transiente NeoPZ",
                      res=rc)
        figura.figura(os.path.join(DADOS, "real0"),
                      os.path.join(AQUI, "figs", "mecanismo_real0.pdf"), "phi",
                      "Realizacao 0 (campos lognormais de c, phi, k)", res=rd)
        print("figuras: figs/mecanismo_det.pdf, figs/mecanismo_real0.pdf")
    except Exception as e:  # figura nao e criterio de aprovacao
        print("aviso: figura nao gerada:", e)
    print(f"tempo total: {time.perf_counter() - T0:.0f} s")
    if falhas:
        print("FALHAS:", *falhas, sep="\n  ")
        sys.exit(1)
    print("TODOS OS TESTES PASSARAM")
