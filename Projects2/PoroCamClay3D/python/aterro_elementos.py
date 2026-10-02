"""
aterro_elementos.py
===================

Aterro da Itasca (FLAC3D, "Embankment Loading on a Cam-Clay Foundation") com o FE u-p generalizado
(fe3d_up.py) e quatro tipos de elemento: hex8 (Q1-Q1), hex20 (Q2-Q1), tet4 (P1-P1) e tet10 (P2-P1).
Mesmos dados de aterro_itasca.py. Inclui a verificação de equilíbrio global (somas das forças).
"""
import time
import numpy as np
import fe3d_up as fe
import aterro_itasca as at


def build(etype='hex8', nx=at.NX, nz=at.NZ, par_of_x=at.par_at, init='gp'):
    mesh = fe.box_mesh((at.LX, at.LY, at.LZ), (nx, at.NY, nz), etype)
    model = fe.PoroField(mesh, par_of_x, alpha=1.0, biot_modulus=at.MB, perm=at.PERM, gam_w=at.GAM_W,
                         body=(0.0, 0.0, -at.GAM_SAT), p0_of_x=lambda x: at.GAM_W * (at.LZ - x[2]), init=init)
    fixed_u = {}
    for nd in fe.nodes_on(mesh, 1) + fe.nodes_on(mesh, 2): fixed_u[3 * nd] = 0.0      # u_x = 0 em x = 0 e 20
    for nd in fe.nodes_on(mesh, 3) + fe.nodes_on(mesh, 4): fixed_u[3 * nd + 1] = 0.0  # u_y = 0 (fatia)
    for nd in fe.nodes_on(mesh, 5):                                                    # base fixa
        for c in range(3): fixed_u[3 * nd + c] = 0.0
    fixed_p = {int(model.pmap[nd]): 0.0 for nd in fe.nodes_on(mesh, 6) if model.pmap[nd] >= 0}
    C = mesh['coords']
    # 'range position-x 0 4': faces do topo com centro em x <= 4
    Fs = fe.face_load(mesh, lambda f, m: m == 6 and C[f, 0].mean() <= at.XLOAD + 1e-9, [0.0, 0.0, -at.LOAD])
    return mesh, model, fixed_u, fixed_p, Fs


def node_at(mesh, x, y, z):
    return int(np.argmin(np.linalg.norm(mesh['coords'] - np.array([x, y, z]), axis=1)))


def equilibrium(model, Fext, fixed_u, label=''):
    """somas das forças: R = F_int - G P - F_ext (reações); ΣR + ΣF_ext = 0 em x, y, z;
    R ≈ 0 nos graus de liberdade livres. Reações por apoio, cada componente contada uma única vez:
    R_x nas faces x = 0 e x = Lx (sem os nós da base), R_y nas faces y = 0 e y = Ly (sem a base),
    e (R_x, R_y, R_z) na base."""
    R = model.Fint - model.Q @ model.P - Fext
    free = np.setdiff1d(np.arange(model.nu), list(fixed_u))
    out = {}
    for c, nm in enumerate('xyz'):
        out['F' + nm] = Fext[c::3].sum(); out['R' + nm] = R[c::3].sum()
        out['I' + nm] = (model.Fint - model.Q @ model.P)[c::3].sum()      # forças internas: soma nula
    out['Rfree'] = np.abs(R[free]).max()
    mesh = model.mesh
    base = set(fe.nodes_on(mesh, 5))
    side = {m: np.array(sorted(set(fe.nodes_on(mesh, m)) - base)) for m in (1, 2, 3, 4)}
    nb = np.array(sorted(base))
    out['R_x0'] = R[3 * side[1]].sum(); out['R_xL'] = R[3 * side[2]].sum()
    out['R_y0'] = R[3 * side[3] + 1].sum(); out['R_yL'] = R[3 * side[4] + 1].sum()
    out['R_base'] = [R[3 * nb + c].sum() for c in range(3)]
    return out


def run(etype='hex8', n_load=10, times=at.TIMES, nx=at.NX, nz=at.NZ, par_of_x=at.par_at, verbose=True, init='gp'):
    mesh, model, fixed_u, fixed_p, Fs = build(etype, nx, nz, par_of_x, init)
    F0, _, _ = model.assemble(model.U)
    model.Fint = F0
    eq0 = equilibrium(model, model.Fb, fixed_u)
    if verbose:
        print(f'{etype} {nx}x1x{nz}: {model.nn} nós, {len(mesh["elements"])} elementos, {model.npd} gdl de p; '
              f'|R| livre inicial = {eq0["Rfree"]:.1e}')
    mon_u = {x: node_at(mesh, x, 0.0, at.LZ) for x in (0.0, 2.0, 4.0, 6.0)}
    hist = dict(t=[0.0], **{f'uz{int(x)}': [0.0] for x in mon_u},
                pp1=[model.zone_pp(0.5, 9.5)], pp2=[model.zone_pp(1.5, 7.5)])

    def rec(t):
        hist['t'].append(t)
        for x, nd in mon_u.items(): hist[f'uz{int(x)}'].append(model.U[3 * nd + 2])
        hist['pp1'].append(model.zone_pp(0.5, 9.5)); hist['pp2'].append(model.zone_pp(1.5, 7.5))

    t0 = time.time(); ncut = [0]

    def advance(lam0, lam1, ta, tb, flow, level=0):
        try:
            return model.step(model.Fb + lam1 * Fs, fixed_u, fixed_p, dt=max(tb - ta, 1.0), flow=flow)
        except (fe.mcc.ReturnMappingError, RuntimeError, np.linalg.LinAlgError):
            if level >= 10:
                raise
            ncut[0] += 1
            lm, tm = 0.5 * (lam0 + lam1), 0.5 * (ta + tb)
            return advance(lam0, lm, ta, tm, flow, level + 1) + advance(lm, lam1, tm, tb, flow, level + 1)

    for k in range(1, n_load + 1):
        advance((k - 1) / n_load, k / n_load, 0.0, 0.0, False)
        rec(0.0)
    snap_u = dict(U=model.U.copy(), P=model.P.copy(), eq=equilibrium(model, model.Fb + Fs, fixed_u))
    if verbose:
        print(f'   fim não drenado: uz(0) = {model.U[3 * mon_u[0.0] + 2]:.4f}  pp1 = {hist["pp1"][-1]:.2f}  '
              f'pp2 = {hist["pp2"][-1]:.2f}')
    t = 0.0
    for tn in times:
        advance(1.0, 1.0, t, tn, True)
        t = tn
        rec(t)
    snap_d = dict(U=model.U.copy(), P=model.P.copy(), eq=equilibrium(model, model.Fb + Fs, fixed_u))
    if verbose:
        print(f'   t = 1e8: uz(0) = {model.U[3 * mon_u[0.0] + 2]:.4f}  uz(2) = {model.U[3 * mon_u[2.0] + 2]:.4f}  '
              f'uz(4) = {model.U[3 * mon_u[4.0] + 2]:.4f}  uz(6) = {model.U[3 * mon_u[6.0] + 2]:.4f}  '
              f'pp1 = {hist["pp1"][-1]:.2f}  pp2 = {hist["pp2"][-1]:.2f}  ({time.time() - t0:.0f} s, cortes {ncut[0]})')
    return dict(mesh=mesh, model=model, hist={k: np.array(v) for k, v in hist.items()}, undrained=snap_u,
                drained=snap_d, n_load=n_load, cuts=ncut[0], eq0=eq0)


if __name__ == '__main__':
    import sys
    init = 'centroid' if 'centroid' in sys.argv else 'gp'
    ets = [a for a in sys.argv[1:] if a != 'centroid']
    for et in (ets or ['hex8', 'hex20', 'tet10', 'tet4']):
        r = run(et, init=init)
        e = r['drained']['eq']
        print(f'   equilíbrio final: ΣFext = ({e["Fx"]:.3e}, {e["Fy"]:.3e}, {e["Fz"]:.3f}) kN;  '
              f'Σreações = ({e["Rx"]:.3e}, {e["Ry"]:.3e}, {e["Rz"]:.3f}) kN;  max|R| livre = {e["Rfree"]:.1e}')
        print(f'   reações em x: face x=0 {e["R_x0"]:.3f}, face x=20 {e["R_xL"]:.3f}, base {e["R_base"][0]:.3f} '
              f'-> soma {e["R_x0"] + e["R_xL"] + e["R_base"][0]:.2e} kN')
        print(f'   reações em y: face y=0 {e["R_y0"]:.3f}, face y=1 {e["R_yL"]:.3f}, base {e["R_base"][1]:.3f} '
              f'-> soma {e["R_y0"] + e["R_yL"] + e["R_base"][1]:.2e} kN')
        print(f'   reação em z (base) {e["R_base"][2]:.3f} kN;  Σ(F_int - G P) = ({e["Ix"]:.1e}, {e["Iy"]:.1e}, {e["Iz"]:.1e})')
