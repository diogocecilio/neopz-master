"""
aterro_itasca.py
================

Exemplo da Itasca (FLAC3D) "Embankment Loading on a Cam-Clay Foundation" resolvido com o FE u-p de
hexaedros + Cam-Clay modificado (poro_camclay_fem.py / modified_cam_clay.py). Unidades: kPa, m, s.

Fatia de 1 m (deformação plana) com meia simetria: 20 × 1 × 10 m, malha 20 × 1 × 10 (como no FLAC3D).
Fundação Cam-Clay modificado: ν = 0.3, M = 0.888, λ = 0.161, κ = 0.062, p_ref = 1 kPa, v_λ = 2.858,
p'c0 = 160 kPa; densidade seca 2000 kg/m³, n = 0.3 (ρ_sat = 2300 kg/m³);
fluido: K_f = 2·10⁵ kPa (biot off: α = 1, 1/M = n/K_f), mobilidade k = 10⁻¹² m²/(Pa·s) = 10⁻⁹ m²/(kPa·s).
Estado inicial: N.A. na superfície; 'zone initialize-stresses ratio 0.7' aplicado antes da poropressão
(razão sobre tensões TOTAIS) → σ'h/σ'v ≈ 6/13, como diz o texto do exemplo. Topo drenante (p = 0).
Etapa 1: sobrecarga de 50 kPa em 0 ≤ x ≤ 4 m sem fluxo (não drenada), em incrementos.
Etapa 2: adensamento acoplado até t = 10⁸ s (passos em escala logarítmica).

Uso: python aterro_itasca.py   -> figuras e resultados em ./figuras_aterro
"""
import os
import time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import modified_cam_clay as mcc
import poro_camclay_fem as fe

AQUI = os.path.dirname(os.path.abspath(__file__))
SAIDA = os.path.join(AQUI, 'figuras_aterro')

# ------------------------------------------------------------------ dados do problema
LX, LY, LZ = 20.0, 1.0, 10.0
NX, NY, NZ = 20, 1, 10
NU, M_CS, LAM, KAP, VLAM, PC0 = 0.3, 0.888, 0.161, 0.062, 2.858, 160.0
RHO_DRY, POR, RHO_W, GRAV = 2.0, 0.3, 1.0, 10.0          # t/m³ e m/s² -> kN/m³ = kPa/m
GAM_SAT = (RHO_DRY + POR * RHO_W) * GRAV                    # 23 kPa/m
GAM_W = RHO_W * GRAV                                        # 10 kPa/m
K0 = 0.7                                                    # razão σh/σv (totais)
KF = 2.0e5                                                  # kPa
MB = KF / POR                                               # módulo de Biot (α = 1, 1/M = n/Kf)
PERM = 1.0e-9                                               # m²/(kPa·s)
LOAD, XLOAD = 50.0, 4.0
TIMES = np.logspace(2, 8, 25)                               # instantes da etapa 2 (s)
FLAC = dict(undrained=0.14, drained=0.19)                   # recalques máximos citados no texto da Itasca


def initial_effective_stress(z):
    """σ' inicial (Voigt do Cam-Clay {xx,xy,xz,yy,yz,zz}, tração +): σh = 0.7 σv em tensões totais."""
    d = LZ - z
    sv_tot, pw = -GAM_SAT * d, GAM_W * d
    sv = sv_tot + pw                                        # σ'v
    sh = K0 * sv_tot + pw                                   # σ'h  (σ'h/σ'v = 6.1/13)
    return np.array([sh, 0.0, 0.0, sh, 0.0, sv])


def par_at(x, shear='hypo_nu', elastic=False):
    """parâmetros do Cam-Clay no ponto x (v0 da NCL + linha κ com o p'0 local).
    elastic=True: mesmo v0, mas p'c0 enorme (só a parte elástica, para estudos de sensibilidade)."""
    s0 = initial_effective_stress(x[2])
    p0 = -(s0[0] + s0[3] + s0[5]) / 3.0
    v0 = VLAM - LAM * np.log(PC0) + KAP * np.log(PC0 / p0)
    return mcc.mcc_parameters(M=M_CS, lam=LAM, kap=KAP, N=VLAM, pc0=1e9 if elastic else PC0, p0=p0, v0=v0,
                              elasticity='pressure_dependent', shear=shear, nu=NU, sigma0=s0)


class PoroField:
    """FE u-p com parâmetros/estado inicial por elemento (centróide, como as zonas do FLAC3D),
    peso próprio, poropressão hidrostática inicial e termo gravitacional no fluxo de Darcy.

    stab: False (padrão, formulação Q1-Q1 do código original), ('ppp', c) = projeção polinomial da pressão
    com τ = c/(2G0), ou True = laplaciano da pressão β = h²/(4(K+4G/3)).  As estabilizações eliminam o
    tabuleiro de xadrez da poropressão nodal no limite não drenado, mas aqui introduzem drenagem artificial
    nos elementos moles do topo (recalque não drenado +20%), por isso não são usadas nos resultados."""

    def __init__(self, mesh, par_of_x, alpha=1.0, biot_modulus=MB, perm=PERM, gam_w=GAM_W, gam_sat=GAM_SAT,
                 stab=False):
        self.mesh = mesh
        self.nn = len(mesh['coords']); self.nu = 3 * self.nn
        self.geo = [fe.element_geometry(mesh['coords'][e]) for e in mesh['elements']]
        self.par = [par_of_x(mesh['coords'][e].mean(axis=0)) for e in mesh['elements']]
        # estado por ponto de Gauss: ε^p, α, σ' e ε do último passo convergido
        self.state = [[dict(epsp=np.zeros(6), al=0.0, sig=P['sigma0'].copy(), eps=np.zeros(6))
                       for g in gps] for gps, P in zip(self.geo, self.par)]
        Qd = sp.lil_matrix((self.nu, self.nn)); Sd = sp.lil_matrix((self.nn, self.nn))
        Hd = sp.lil_matrix((self.nn, self.nn)); self.fg = np.zeros(self.nn); self.Fb = np.zeros(self.nu)
        for e, gps, P_e in zip(mesh['elements'], self.geo, self.par):
            Q, S, H = fe.contribute_porous(gps, alpha, 1.0 / biot_modulus, perm)
            du = np.array([[3 * a, 3 * a + 1, 3 * a + 2] for a in e]).ravel()
            for a in range(24):
                for b in range(8): Qd[du[a], e[b]] += Q[a, b]
            for a in range(8):
                for b in range(8):
                    Sd[e[a], e[b]] += S[a, b]; Hd[e[a], e[b]] += H[a, b]
            if isinstance(stab, tuple) and stab[0] == 'ppp':
                # projeção polinomial da pressão (Dohrmann & Bochev 2004): τ ∫ (N - N̄)ᵀ (N - N̄) dΩ, τ = c/(2 G0)
                tau = stab[1] / (2.0 * P_e['G0'])
                Me = sum(np.outer(g['N'], g['N']) * g['wdJ'] for g in gps)
                me = sum(g['N'] * g['wdJ'] for g in gps); Ve = sum(g['wdJ'] for g in gps)
                Sst = tau * (Me - np.outer(me, me) / Ve)
                for a in range(8):
                    for b in range(8): Sd[e[a], e[b]] += Sst[a, b]
            elif stab:
                # laplaciano da pressão (Aguilar et al. 2008): β ∫ ∇Nᵀ ∇N dΩ (P - Pn), β = h²/(4 (K + 4G/3))
                X = mesh['coords'][e]
                h = min(np.ptp(X[:, 0]), np.ptp(X[:, 1]), np.ptp(X[:, 2]))
                beta = h**2 / (4.0 * (P_e['K0'] + 4.0 * P_e['G0'] / 3.0))
                Sst = beta * sum(g['G'].T @ g['G'] * g['wdJ'] for g in gps)
                for a in range(8):
                    for b in range(8): Sd[e[a], e[b]] += Sst[a, b]
            for g in gps:
                # f_g = ∫ k ∇Nᵀ (ρw g⃗) dΩ  (g⃗ = -g e_z): fluxo nulo no estado hidrostático
                self.fg[e] += perm * (g['G'].T @ np.array([0.0, 0.0, -gam_w])) * g['wdJ']
                # F_b = ∫ Nᵀ b dΩ, b = (0, 0, -γ_sat)
                self.Fb[du[2::3]] += -gam_sat * g['N'] * g['wdJ']
        self.Q, self.S, self.H = Qd.tocsr(), Sd.tocsr(), Hd.tocsr()
        self.U = np.zeros(self.nu)
        self.P = gam_w * (LZ - mesh['coords'][:, 2])         # poropressão hidrostática
        self.trial = None

    def assemble(self, U):
        """forças internas efetivas e tangente consistente (esparsa)."""
        F = np.zeros(self.nu); rows, cols, vals = [], [], []; trial = []
        for e, gps, st, P in zip(self.mesh['elements'], self.geo, self.state, self.par):
            du = np.array([[3 * a, 3 * a + 1, 3 * a + 2] for a in e]).ravel()
            ue = U[du]; fe_ = np.zeros(24); ke = np.zeros((24, 24)); tr = []
            for g, s_n in zip(gps, st):
                eps = (g['B'] @ ue)[fe.FE2MCC]
                r = mcc.return_mapping(P, eps, s_n['epsp'], s_n['al'], sig_n=s_n['sig'], eps_n=s_n['eps'])
                r['eps'] = eps
                fe_ += g['B'].T @ r['stress'][fe.MCC2FE] * g['wdJ']
                ke += g['B'].T @ r['Dep'][np.ix_(fe.MCC2FE, fe.MCC2FE)] @ g['B'] * g['wdJ']
                tr.append(r)
            F[du] += fe_
            rows.append(np.repeat(du, 24)); cols.append(np.tile(du, 24)); vals.append(ke.ravel())
            trial.append(tr)
        K = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                          shape=(self.nu, self.nu))
        return F, K, trial

    def step(self, Fext, fixed_u, fixed_p, dt, flow, tol=1e-8, maxit=30, verbose=False):
        """Euler implícito + Newton.  Ru = Fint - Q P - Fext ;
        Rp = Qᵀ(U - Un) + S(P - Pn) + Δt (H P - f_g)  (H e f_g só com fluxo)."""
        nu, nn = self.nu, self.nn
        U = self.U.copy(); P = self.P.copy()
        for i, v in fixed_u.items(): U[i] = v
        for i, v in fixed_p.items(): P[i] = v
        fixed = np.array(sorted(list(fixed_u) + [nu + i for i in fixed_p]), int)
        free = np.setdiff1d(np.arange(nu + nn), fixed)
        Hs = self.H if flow else 0.0 * self.H
        fgs = self.fg if flow else 0.0 * self.fg
        ref = max(np.linalg.norm(Fext), 1.0)
        for it in range(1, maxit + 1):
            F, K, trial = self.assemble(U)
            Ru = F - self.Q @ P - Fext
            Rp = self.Q.T @ (U - self.U) + self.S @ (P - self.P) + dt * (Hs @ P - fgs)
            R = np.r_[Ru, Rp]
            nr = np.linalg.norm(R[free]) / ref
            if verbose:
                print(f'      it {it}: |R|/|F| = {nr:.3e}')
            if nr < tol:
                break
            J = sp.bmat([[K, -self.Q], [self.Q.T, self.S + dt * Hs]], format='csr')
            dx = np.zeros(nu + nn)
            dx[free] = spla.spsolve(J[free][:, free].tocsc(), -R[free])
            U += dx[:nu]; P += dx[nu:]
        else:
            raise RuntimeError(f'Newton não convergiu (|R|/|F| = {nr:.3e})')
        self.U, self.P = U, P
        self.state = [[dict(epsp=r['plastic_strain'], al=r['alpha'], sig=r['stress'], eps=r['eps'])
                       for r in tr] for tr in trial]
        self.trial = trial
        return it


def build(nx=NX, nz=NZ, par_of_x=par_at, stab=False):
    mesh = fe.box_mesh((LX, LY, LZ), (nx, NY, nz))
    model = PoroField(mesh, par_of_x, stab=stab)
    fixed_u = {}
    for nd in fe.nodes_on(mesh, 1) + fe.nodes_on(mesh, 2): fixed_u[3 * nd] = 0.0      # ux = 0 em x = 0 e x = 20
    for nd in fe.nodes_on(mesh, 3) + fe.nodes_on(mesh, 4): fixed_u[3 * nd + 1] = 0.0  # uy = 0 (deformação plana)
    for nd in fe.nodes_on(mesh, 5):                                                    # base fixa
        for c in range(3): fixed_u[3 * nd + c] = 0.0
    fixed_p = {nd: 0.0 for nd in fe.nodes_on(mesh, 6)}                                  # topo drenante
    # sobrecarga: faces do topo com centro em 0 <= x <= 4 (como 'range position-x 0 4' do FLAC3D)
    Fs = np.zeros(model.nu)
    for f, m in mesh['faces']:
        if m == 6 and mesh['coords'][f][:, 0].mean() <= XLOAD:
            Fs += fe.face_load(dict(coords=mesh['coords'], faces=[(f, 99)]), 99, [0.0, 0.0, -LOAD])
    return mesh, model, fixed_u, fixed_p, Fs


def node_at(mesh, x, y, z):
    return int(np.argmin(np.linalg.norm(mesh['coords'] - np.array([x, y, z]), axis=1)))


def zone_pp(mesh, model, x, z):
    """poropressão média (nós) do elemento que contém (x, z), como 'zone history pore-pressure'."""
    for e in mesh['elements']:
        X = mesh['coords'][e]
        if X[:, 0].min() <= x <= X[:, 0].max() and X[:, 2].min() <= z <= X[:, 2].max():
            return model.P[e].mean()


def run(n_load=10, times=TIMES, nx=NX, nz=NZ, par_of_x=par_at, verbose=True, drained_only=False, stab=False):
    """etapa 1 (não drenada, n_load incrementos) + etapa 2 (adensamento nos instantes 'times').
    drained_only=True: aplica a carga com poropressão hidrostática imposta em todos os nós (drenado)."""
    mesh, model, fixed_u, fixed_p, Fs = build(nx, nz, par_of_x, stab)
    F0, _, _ = model.assemble(model.U)
    R0 = F0 - model.Q @ model.P - model.Fb
    freeu = np.setdiff1d(np.arange(model.nu), list(fixed_u))
    if verbose:
        print(f'malha {nx}x1x{nz}: {model.nn} nós; |R_u(t=0)| nos gdl livres = {np.linalg.norm(R0[freeu]):.2e}')
    mon_u = {x: node_at(mesh, x, 0.0, LZ) for x in (0.0, 2.0, 4.0, 6.0)}
    hist = dict(t=[0.0], **{f'uz{int(x)}': [0.0] for x in mon_u},
                pp1=[zone_pp(mesh, model, 0.5, 9.5)], pp2=[zone_pp(mesh, model, 1.5, 7.5)])

    def rec(t):
        hist['t'].append(t)
        for x, nd in mon_u.items(): hist[f'uz{int(x)}'].append(model.U[3 * nd + 2])
        hist['pp1'].append(zone_pp(mesh, model, 0.5, 9.5)); hist['pp2'].append(zone_pp(mesh, model, 1.5, 7.5))

    t_start = time.time(); ncut = [0]
    fp = {i: model.P[i] for i in range(model.nn)} if drained_only else fixed_p

    def advance(lam0, lam1, t0, t1, flow, level=0):
        """avança de (λ0, t0) para (λ1, t1); se falhar, divide o incremento ao meio (até 10 níveis)."""
        saved = (model.U.copy(), model.P.copy(), model.state)
        try:
            return model.step(model.Fb + lam1 * Fs, fixed_u, fp, dt=max(t1 - t0, 1.0), flow=flow)
        except (mcc.ReturnMappingError, RuntimeError, np.linalg.LinAlgError):
            model.U, model.P, model.state = saved
            if level >= 10:
                raise
            ncut[0] += 1
            lm, tm = 0.5 * (lam0 + lam1), 0.5 * (t0 + t1)
            return advance(lam0, lm, t0, tm, flow, level + 1) + advance(lm, lam1, tm, t1, flow, level + 1)

    def snapshot():
        return dict(U=model.U.copy(), P=model.P.copy(),
                    plastic=np.array([[r['plastic'] for r in tr] for tr in model.trial]),
                    yielded=np.array([[s['al'] != 0.0 for s in st] for st in model.state]),
                    stress=np.array([[r['stress'] for r in tr] for tr in model.trial]))

    for k in range(1, n_load + 1):
        it = advance((k - 1) / n_load, k / n_load, 0.0, 0.0, False)
        rec(0.0)
        if verbose:
            print(f'  carga {k / n_load:5.2f}: {it} it.  uz(0) = {model.U[3 * mon_u[0.0] + 2]:.4f} m')
    snap_u = snapshot()
    t = 0.0
    if not drained_only:
        for tn in times:
            it = advance(1.0, 1.0, t, tn, True)
            t = tn
            rec(t)
            if verbose:
                print(f'  t = {t:9.3e} s: {it} it.  uz(0) = {model.U[3 * mon_u[0.0] + 2]:.4f} m   '
                      f'pp1 = {hist["pp1"][-1]:7.3f}  pp2 = {hist["pp2"][-1]:7.3f} kPa')
    if verbose:
        print(f'  tempo total {time.time() - t_start:.1f} s, subdivisões de passo: {ncut[0]}')
    snap_d = snapshot()
    return dict(mesh=mesh, model=model, hist={k: np.array(v) for k, v in hist.items()},
                undrained=snap_u, drained=snap_d, n_load=n_load, cuts=ncut[0])


# ---------------------------------------------------------------------------------------------- figuras
def figures(res, sens, txt):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    INK, INK2, MUTED, GRID = '#0b0b0b', '#52514e', '#898781', '#e1e0d9'
    C_FE, C_CSL, C_REF = '#eb6834', '#1baf7a', '#2a78d6'
    plt.rcParams.update({
        'figure.facecolor': '#fcfcfb', 'axes.facecolor': '#fcfcfb', 'savefig.facecolor': '#fcfcfb',
        'axes.edgecolor': '#c3c2b7', 'axes.labelcolor': INK, 'xtick.color': INK2, 'ytick.color': INK2,
        'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6, 'font.size': 10,
        'axes.titlesize': 10.5, 'axes.titleweight': 'bold', 'legend.frameon': False})
    os.makedirs(SAIDA, exist_ok=True)
    mesh, hist = res['mesh'], res['hist']
    X = mesh['coords']
    sel = np.where(np.abs(X[:, 1]) < 1e-9)[0]                  # nós do plano y = 0
    xs, zs = np.unique(X[sel, 0]), np.unique(X[sel, 2])
    grid = {(round(X[i, 0], 9), round(X[i, 2], 9)): i for i in sel}
    ig = np.array([[grid[(round(x, 9), round(z, 9))] for x in xs] for z in zs])
    XX, ZZ = np.meshgrid(xs, zs)
    stages = [('undrained', 'Fim da etapa não drenada'), ('drained', 'Fim do adensamento (t = 10⁸ s)')]

    def frame(ax):
        ax.set_xlim(-0.3, LX + 0.3); ax.set_ylim(-0.3, LZ + 1.2); ax.set_aspect('equal')
        ax.set_xlabel('x (m)'); ax.set_ylabel('z (m)'); ax.grid(False)
        ax.add_patch(plt.Rectangle((0, LZ), XLOAD, 0.35, color=MUTED, alpha=0.55, lw=0))
        ax.text(XLOAD / 2, LZ + 0.55, '50 kPa', ha='center', va='bottom', fontsize=8.5, color=INK2)

    def mesh_lines(ax, lw=0.4, c='#c3c2b7'):
        for x in xs: ax.plot([x, x], [0, LZ], color=c, lw=lw)
        for z in zs: ax.plot([0, LX], [z, z], color=c, lw=lw)

    # ---- Fig 1: geometria
    fig, ax = plt.subplots(figsize=(9, 5))
    mesh_lines(ax, 0.6, '#9a998f'); frame(ax)
    ax.plot([0, LX], [LZ, LZ], color=C_REF, lw=2.5)
    ax.text(LX - 0.2, LZ + 0.25, 'p = 0 (drenante)', ha='right', color=C_REF, fontsize=9)
    for (x, z), lab in [((0.5, 9.5), 'pp1'), ((1.5, 7.5), 'pp2')]:
        ax.plot(x, z, 's', color=C_CSL, ms=7); ax.text(x + 0.25, z - 0.1, lab, color=C_CSL, fontsize=9)
    for x in (0, 2, 4, 6):
        ax.plot(x, LZ, 'o', color=C_FE, ms=6, mec='white')
        ax.text(x, LZ - 0.55, f'z{x}', ha='center', color=C_FE, fontsize=8.5)
    ax.text(10, -0.25, 'base fixa, impermeável', ha='center', va='top', fontsize=9, color=INK2)
    ax.text(-0.2, 5, 'u_x = 0 (simetria)', rotation=90, ha='right', va='center', fontsize=9, color=INK2)
    ax.text(LX + 0.2, 5, 'u_x = 0', rotation=90, ha='left', va='center', fontsize=9, color=INK2)
    ax.set_ylim(-1.0, LZ + 1.2); ax.set_xlim(-1.0, LX + 1.0)
    ax.set_title('Geometria, malha 20 × 1 × 10, condições de contorno e pontos monitorados')
    fig.tight_layout(); fig.savefig(os.path.join(SAIDA, 'aterro_geometria.png'), dpi=130); plt.close(fig)

    # ---- Fig 2: vetores de deslocamento (Figs. 2 e 3 da Itasca)
    fig, axs = plt.subplots(2, 1, figsize=(9, 9.2))
    for ax, (key, ttl) in zip(axs, stages):
        U = res[key]['U'].reshape(-1, 3)
        mesh_lines(ax); frame(ax)
        ux, uz = U[ig, 0], U[ig, 2]
        mag = np.hypot(ux, uz).max()
        q = ax.quiver(XX, ZZ, ux, uz, color=C_FE, angles='xy', scale_units='xy', scale=mag / 1.4,
                      width=0.0035, headwidth=3.5)
        ax.set_title(f'{ttl}: vetores de deslocamento (máx. {mag:.3f} m)')
        # pontos de Gauss que plastificaram (α ≠ 0) e plásticos no último passo
        yi, pl = res[key]['yielded'], res[key]['plastic']
        gx = np.array([[g['N'] @ X[e] for g in gps] for e, gps in zip(mesh['elements'], res['model'].geo)])
        on = gx[:, :, 1] < 0.5
        ax.plot(gx[yi & on][:, 0], gx[yi & on][:, 2], 'x', color=C_REF, ms=4, mew=1.0)
        ax.plot(gx[pl & on][:, 0], gx[pl & on][:, 2], '+', color=INK, ms=6, mew=1.0)
    fig.legend(handles=[Line2D([], [], color=C_FE, lw=1.5, label='deslocamento (escala arbitrária)'),
                        Line2D([], [], ls='none', marker='x', color=C_REF, label='PG já plastificado (α ≠ 0)'),
                        Line2D([], [], ls='none', marker='+', color=INK, label='PG plástico no último passo')],
               loc='lower center', ncol=3, fontsize=9)
    fig.tight_layout(rect=(0, 0.04, 1, 1)); fig.savefig(os.path.join(SAIDA, 'aterro_vetores.png'), dpi=130)
    plt.close(fig)

    # ---- Fig 3 e 4: contornos de u_z e de poropressão (Figs. 4-7 da Itasca)
    # zonas (elementos) do plano y = 0: valores médios dos 8 nós (como a poropressão de zona do FLAC3D)
    zone = {}
    for k, e in enumerate(mesh['elements']):
        c = X[e].mean(axis=0)
        zone[(int(np.searchsorted(xs, c[0])) - 1, int(np.searchsorted(zs, c[2])) - 1)] = k
    iz = np.array([[zone[(i, j)] for i in range(len(xs) - 1)] for j in range(len(zs) - 1)])
    for fname, comp, lab, cmap in [('aterro_uz_contornos.png', 'uz', 'u_z (m)', 'viridis'),
                                   ('aterro_pp_contornos.png', 'pp', 'poropressão (kPa)', 'cividis')]:
        fig, axs = plt.subplots(2, 1, figsize=(9, 8.6))
        for ax, (key, ttl) in zip(axs, stages):
            if comp == 'uz':
                V = res[key]['U'].reshape(-1, 3)[ig, 2]
                cs = ax.contourf(XX, ZZ, V, levels=14, cmap=cmap)
                ax.contour(XX, ZZ, V, levels=14, colors='white', linewidths=0.3)
                ext = f'mín. {V.min():.3f} m'
            else:
                # fim não drenado: excesso de poropressão; fim do adensamento: poropressão total
                exc = key == 'undrained'
                Pn = res[key]['P'] - (GAM_W * (LZ - X[:, 2]) if exc else 0.0)
                Pz = np.array([Pn[e].mean() for e in mesh['elements']])[iz]
                cs = ax.pcolormesh(xs, zs, Pz, cmap=cmap, shading='flat', edgecolors='#fcfcfb', linewidth=0.2)
                lab = 'excesso de poropressão (kPa)' if exc else 'poropressão (kPa)'
                ext = f'máx. {Pz.max():.1f} kPa, valores de zona'
            frame(ax)
            cb = fig.colorbar(cs, ax=ax, shrink=0.85, pad=0.02); cb.set_label(lab)
            ax.set_title(f'{ttl}: {lab} ({ext})')
        fig.tight_layout(); fig.savefig(os.path.join(SAIDA, fname), dpi=130); plt.close(fig)

    # ---- Fig 5: históricos (Figs. 8 e 9 da Itasca)
    t = hist['t']; nl = res['n_load']; tc = t[nl:].copy(); tc[0] = 10.0     # fim da etapa 1 desenhado em t = 10 s
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.6))
    cols = [INK, C_FE, C_REF, C_CSL]
    for c, x in zip(cols, (0, 2, 4, 6)):
        axs[0].semilogx(tc, -hist[f'uz{x}'][nl:], '-o', ms=3, color=c, lw=1.4, label=f'z{x}: x = {x} m')
    for v, lab, dy in [(FLAC['undrained'], 'FLAC3D, fim não drenado (≈0.14 m)', -0.014),
                       (FLAC['drained'], 'FLAC3D, final (≈0.19 m)', 0.004)]:
        axs[0].axhline(v, color=MUTED, ls='--', lw=1.0)
        axs[0].text(12, v + dy, lab, fontsize=8, color=INK2)
    axs[0].set_xlabel('tempo (s)   [ponto em t = 10 s = fim da etapa não drenada]'); axs[0].set_ylabel('recalque −u_z (m)')
    axs[0].set_title('Recalques no topo'); axs[0].legend(fontsize=8.5, loc='upper left')
    for c, k, lab, hyd in [(INK, 'pp1', 'pp1 (0.5, 9.5)', 5.0), (C_FE, 'pp2', 'pp2 (1.5, 7.5)', 25.0)]:
        axs[1].semilogx(tc, hist[k][nl:], '-o', ms=3, color=c, lw=1.4, label=lab)
        axs[1].axhline(hyd, color=c, ls=':', lw=1.0)
    axs[1].set_xlabel('tempo (s)'); axs[1].set_ylabel('poropressão (kPa)')
    axs[1].set_title('Poropressão (pontilhado: hidrostática)'); axs[1].legend(fontsize=8.5)
    fig.tight_layout(); fig.savefig(os.path.join(SAIDA, 'aterro_historicos.png'), dpi=130); plt.close(fig)

    # ---- Fig 6: sensibilidade do recalque final à lei elástica
    fig, ax = plt.subplots(figsize=(7.6, 4.8))
    for (lab, h, c, ls) in sens:
        ax.semilogx(tc, -h['uz0'][nl:], ls, color=c, lw=1.5, ms=3, label=lab)
    for v in FLAC.values():
        ax.axhline(v, color=MUTED, ls='--', lw=1.0)
    ax.text(12, FLAC['drained'] + 0.004, 'FLAC3D final ≈ 0.19 m', fontsize=8, color=INK2)
    ax.text(12, FLAC['undrained'] + 0.004, 'FLAC3D não drenado ≈ 0.14 m', fontsize=8, color=INK2)
    ax.set_xlabel('tempo (s)'); ax.set_ylabel('recalque em x = 0 (m)')
    ax.set_title('Sensibilidade do recalque à forma de G(p)'); ax.legend(fontsize=8.5, loc='upper left')
    fig.tight_layout(); fig.savefig(os.path.join(SAIDA, 'aterro_sensibilidade.png'), dpi=130); plt.close(fig)
    open(os.path.join(SAIDA, 'resultados_aterro.txt'), 'w', encoding='utf-8').write(txt + '\n')


def main():
    print('Caso principal (Cam-Clay modificado, G hipoelástico com ν = 0.3, como o FLAC3D):')
    res = run()
    h = res['hist']; nl = res['n_load']
    lines = ['Aterro sobre fundação Cam-Clay (exemplo Itasca/FLAC3D), FE u-p hexa8, malha 20x1x10', '']
    lines.append(f"{'':28s}{'x = 0':>9s}{'x = 2':>9s}{'x = 4':>9s}{'x = 6':>9s}   (u_z, m)")
    for lab, i in [('fim da etapa não drenada', nl), ('t = 1e5 s', None), ('t = 1e6 s', None), ('t = 1e8 s (final)', -1)]:
        if i is None:
            i = int(np.argmin(np.abs(h['t'] - float(lab.split('=')[1].split()[0]))))
        lines.append(f"{lab:28s}" + ''.join(f"{h[f'uz{x}'][i]:9.4f}" for x in (0, 2, 4, 6)))
    lines.append('')
    lines.append(f"poropressão pp1 (0.5, 9.5): inicial {h['pp1'][0]:.2f}, não drenado {h['pp1'][nl]:.2f}, final {h['pp1'][-1]:.2f} kPa (hidrostática 5)")
    lines.append(f"poropressão pp2 (1.5, 7.5): inicial {h['pp2'][0]:.2f}, não drenado {h['pp2'][nl]:.2f}, "
                 f"máx {h['pp2'].max():.2f}, final {h['pp2'][-1]:.2f} kPa (hidrostática 25)")
    tt = h['t'][nl:]; uu = -h['uz0'][nl:]
    t19 = np.exp(np.interp(FLAC['drained'], uu, np.log(np.maximum(tt, 1.0))))
    lines.append(f"FLAC3D (texto): recalque máximo ≈ {FLAC['undrained']} m (não drenado) -> ≈ {FLAC['drained']} m (drenado)")
    lines.append(f"este código: {-h['uz0'][nl]:.3f} m -> {-h['uz0'][-1]:.3f} m; o recalque passa por 0.19 m em t ≈ {t19:.2e} s")
    lines.append(f"PGs plastificados (α ≠ 0): não drenado {res['undrained']['yielded'].sum()}, final {res['drained']['yielded'].sum()} de {res['drained']['yielded'].size}")
    lines.append('')
    print('Sensibilidade / convergência:')
    sens_runs = {}
    t0 = time.time()
    sens_runs['el_hypo'] = run(par_of_x=lambda x: par_at(x, 'hypo_nu', True), verbose=False)
    sens_runs['el_sec'] = run(par_of_x=lambda x: par_at(x, 'constant_nu', True), verbose=False)
    sens_runs['dren'] = run(drained_only=True, n_load=20, verbose=False)
    sens_runs['fina'] = run(nx=40, nz=20, verbose=False)
    sens_runs['dt'] = run(n_load=20, times=np.logspace(2, 8, 61), verbose=False)
    sens_runs['ppp'] = run(stab=('ppp', 1.0), verbose=False)
    sens_runs['fpl'] = run(stab=True, verbose=False)
    print(f'  ({time.time() - t0:.0f} s)')
    g = lambda r, i: -r['hist']['uz0'][i]
    lines.append('Sensibilidade e convergência (recalque em x = 0, m):')
    lines.append(f"  caso principal (hipoelástico, elastoplástico, 20x10, 10 + 25 passos): {g(res, nl):.4f} -> {g(res, -1):.4f}")
    lines.append(f"  malha 40x1x20:                                                       {g(sens_runs['fina'], nl):.4f} -> {g(sens_runs['fina'], -1):.4f}")
    lines.append(f"  20 incrementos de carga + 61 passos de tempo:                        {g(sens_runs['dt'], 20):.4f} -> {g(sens_runs['dt'], -1):.4f}")
    lines.append(f"  só elástico, G hipoelástico (ds = 2 G(p_n) de):                      {g(sens_runs['el_hypo'], nl):.4f} -> {g(sens_runs['el_hypo'], -1):.4f}")
    lines.append(f"  só elástico, G secante (s = 2 G(p) e^e):                             {g(sens_runs['el_sec'], nl):.4f} -> {g(sens_runs['el_sec'], -1):.4f}")
    lines.append(f"  carga aplicada drenada (sem etapa não drenada):                      final {g(sens_runs['dren'], -1):.4f}")
    lines.append('')
    lines.append('Estabilização da pressão (Q1-Q1 viola a condição LBB no limite não drenado):')
    for k, lab in [('ppp', 'projeção polinomial, τ = 1/(2G0)'), ('fpl', 'laplaciano, β = h²/(4(K+4G/3))')]:
        r = sens_runs[k]
        lines.append(f"  {lab:34s}: recalque {g(r, nl):.4f} -> {g(r, -1):.4f} m; pp1 não drenado {r['hist']['pp1'][nl]:.1f} kPa "
                     f"(sem estabilização {h['pp1'][nl]:.1f})")
    txt = '\n'.join(lines)
    print(txt)
    sens = [('caso principal (hipoelástico, elastoplástico)', res['hist'], '#eb6834', '-o'),
            ('só elástico, G hipoelástico', sens_runs['el_hypo']['hist'], '#0b0b0b', '--'),
            ('só elástico, G secante s = 2G(p)e^e', sens_runs['el_sec']['hist'], '#2a78d6', '-.'),
            ('malha 40 × 1 × 20', sens_runs['fina']['hist'], '#1baf7a', ':')]
    figures(res, sens, txt)
    np.savez(os.path.join(SAIDA, 'aterro_referencia.npz'), **{k: v for k, v in res['hist'].items()},
             U_undrained=res['undrained']['U'], P_undrained=res['undrained']['P'],
             U_drained=res['drained']['U'], P_drained=res['drained']['P'])
    print('Figuras em', SAIDA)


if __name__ == '__main__':
    main()
