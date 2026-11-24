# mc2.py: python3 mc2.py cov lv nreal {ip001|ip01|tab} semente N
import numpy as np, sys
from col1d import *
cov, lv, nreal, var, seed, N = float(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]), sys.argv[4], int(sys.argv[5]), int(sys.argv[6])
H = 2.0; Nf = 400; dzf = H / Nf; yf = (np.arange(Nf) + 0.5) * dzf
lam, V = np.linalg.eigh(np.exp(-np.abs(yf[:, None] - yf[None, :]) / lv) * dzf)
idx = np.argsort(lam)[::-1][:100]; lam = lam[idx]; phi = V[:, idx] / np.sqrt(dzf)
sl = np.sqrt(np.log(1 + cov ** 2)); mu = np.log(3e-6) - 0.5 * sl ** 2
rng = np.random.default_rng(seed)
so = {'ip001': Solo(psi_e_kpa=0.01), 'ip01': Solo(psi_e_kpa=0.1), 'tab': SoloTab(eps=1e-3)}[var]
y = (np.arange(N) + 0.5) * H / N; res = []
for r in range(nreal):
    xi = rng.standard_normal(100)
    lnk = np.interp(y, yf, mu + sl * phi @ (np.sqrt(lam) * xi))
    col = Coluna(so, N=N, kmult=np.exp(lnk) / so.ks); linha = []
    for T in [12, 18, 24]:
        col.advance_to(T * 3600., dt=360., metodo='newton')
        z, FS, h = fs_perfil(col); k = np.argmin(FS); linha += [FS[k], z[k]]
    res.append(linha)
res = np.array(res)
for j, T in enumerate([12, 18, 24]):
    print(T, 'h: Pf =', np.mean(res[:, 2 * j] < 1), ' base =', np.mean(res[:, 2 * j + 1] > H - H / N))
