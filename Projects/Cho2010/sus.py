#!/usr/bin/env python3
"""sus.py: simulacao por subconjuntos (subset simulation, Au e Beck 2001) com o executavel do Cho2010 ou do rebaixamento.

    python3 sus.py <exe> <pasta> [prog=cho] [M=2000] [rcphi=-0.7] [fator=tri2] [nref=4] [metodo=srm] [srm_f0=0.8]
                   [transiente=tri3] [nystrom=] [N=500] [p0=0.1] [rho=0.8] [semente=1] [paralelos=24] [linear=0]

P_f = P(F_1) P(F_2 | F_1) ... P(F | F_m): F_j = {FS < b_j}, com os limiares b_1 > b_2 > ... escolhidos para que cada
probabilidade condicional seja p0 (o quantil p0 dos FS do nivel). Nivel 0: MC cru com N amostras. Em cada nivel, as
N p0 amostras de menor FS sao sementes de cadeias de Markov de 1/p0 estados que ficam em F_j. O passo e o pCN
(Papaioannou et al. 2015): xi' = rho xi + sqrt(1 - rho^2) z, que preserva a N(0, I) sem razao de aceitacao: o candidato
fica se FS(xi') < b_j, senao a cadeia repete o estado. Quando o quantil fica abaixo de 1, a ultima probabilidade e a
fracao de falhas (FS < 1) do nivel. Os xi tem todos os termos de cada campo (2M no Cho, 3M no rebaixamento): nao usa o
FORM nem nenhuma direcao.

O CoV de uma rodada sai pela formula de Au e Beck (2001): delta_j^2 = (1 - p_j)/(N p_j) (1 + gamma_j), com gamma_j da
correlacao das indicadoras ao longo das cadeias, e delta^2 ~ soma dos delta_j^2 (os niveis tratados como independentes).

linear=beta: teste sem elementos finitos, FS = exp(beta - alfa . xi) com alfa fixo nos 20 primeiros termos de cada campo:
P_f exata = Phi(-beta).

Saidas na pasta: analises.csv (nivel, cadeia, passo, FS), nivel_<j>.json (o limiar, as sementes do proximo nivel, as
indicadoras) para retomar e sus.json com o resultado.
"""
import concurrent.futures as cf
import csv
import json
import math
import os
import random
import subprocess
import sys
import threading

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from form_kl import comando  # noqa: E402


def phi_normal(z):
    return 0.5 * math.erfc(-z / math.sqrt(2.))


def parametros(argv):
    p = {"prog": "cho", "M": 2000, "rcphi": -0.7, "fator": "tri2", "nref": 4, "metodo": "srm", "srm_f0": 0.8,
         "transiente": "tri3", "nystrom": "", "N": 500, "p0": 0.1, "rho": 0.8, "semente": 1, "paralelos": 24,
         "linear": 0.}
    for a in argv:
        k, v = a.split("=", 1)
        if k not in p:
            sys.exit(f"opcao desconhecida: {a}")
        p[k] = type(p[k])(v)
    if p["prog"] == "rebaixamento":
        p["metodo"] = "gim"
    return p


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    exe, pasta = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2])
    p = parametros(sys.argv[3:])
    os.makedirs(pasta, exist_ok=True)
    nc = 3 if p["prog"] == "rebaixamento" else 2
    dim = nc * p["M"]
    N, p0, rho = p["N"], p["p0"], p["rho"]
    ns = int(round(N * p0))       # sementes (cadeias) por nivel
    lc = int(round(1. / p0))      # estados por cadeia
    assert ns * lc == N, "N p0 e 1/p0 precisam ser inteiros"
    trava = threading.Lock()
    reg = open(os.path.join(pasta, "analises.csv"), "a")
    cont = [0]

    alfa = None
    if p["linear"] > 0.:
        alfa = [0.] * dim
        g = random.Random(12345)
        for f in range(nc):
            for k in range(20):
                alfa[f * p["M"] + k] = g.gauss(0., 1.) / math.sqrt(k + 1.)
        n = math.sqrt(sum(a * a for a in alfa))
        alfa = [a / n for a in alfa]

    def fs(x, rotulo, info):
        if alfa is not None:
            v = math.exp(p["linear"] - sum(a * b for a, b in zip(alfa, x)))
        else:
            arq = os.path.join(pasta, rotulo + ".xi")
            with open(arq, "w") as f:
                f.write("\n".join(repr(t) for t in x) + "\n")
            env = dict(os.environ, PZ_THREADS="0", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
            r = subprocess.run(comando(exe, p, arq, rotulo), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                               env=env, cwd=pasta)
            v = None
            try:
                with open(os.path.join(pasta, rotulo + "_0_0.csv"), newline="") as f:
                    linha = next(csv.DictReader(f))
                if r.returncode == 0 and linha["situacao"] == "ok":
                    v = float(linha["FS_min"])
            except (FileNotFoundError, StopIteration):
                pass
            for suf in (".xi", "_0_0.csv", "_mapa_0_0.csv", "_parametros.txt"):
                try:
                    os.remove(os.path.join(pasta, rotulo + suf))
                except FileNotFoundError:
                    pass
        with trava:
            cont[0] += 1
            reg.write(f"{info},{v!r}\n")
            reg.flush()
        return v  # None: a analise falhou (no nivel 0 a amostra e refeita; na cadeia o candidato e recusado)

    rng = random.Random(p["semente"])
    gauss = lambda: [rng.gauss(0., 1.) for _ in range(dim)]
    resumo = []

    # nivel 0: MC cru (as que falham sao trocadas por outro sorteio, na mesma ordem para retomar igual)
    arq0 = os.path.join(pasta, "nivel_0.json")
    if os.path.exists(arq0):
        est = json.load(open(arq0))
        rng.setstate((est["rng"][0], tuple(est["rng"][1]), est["rng"][2]))
    else:
        xs = [gauss() for _ in range(N)]
        with cf.ThreadPoolExecutor(p["paralelos"]) as ex:
            vs = list(ex.map(lambda i: fs(xs[i], f"n0_{i}", f"0,{i},0"), range(N)))
        while any(v is None for v in vs):
            for i in [i for i, v in enumerate(vs) if v is None]:
                xs[i] = gauss()
                vs[i] = fs(xs[i], f"n0_{i}b", f"0,{i},0")
        est = {"xs": xs, "fs": vs, "cadeias": [[v] for v in vs], "pf": 1., "nivel": 0}
        est["rng"] = list(rng.getstate())
        json.dump(est, open(arq0, "w"))
        print(f"nivel 0: {N} analises", flush=True)

    j = 0
    pf = 1.
    deltas2 = []
    while True:
        xs, vs, cad = est["xs"], est["fs"], est["cadeias"]
        ordem = sorted(range(N), key=lambda i: vs[i])
        b = 0.5 * (vs[ordem[ns - 1]] + vs[ordem[ns]])
        final = b <= 1. or j >= 12
        lim = 1. if final else b
        ind = [[1 if v < lim else 0 for v in c] for c in cad]
        pj = sum(map(sum, ind)) / N
        # gamma de Au e Beck: correlacao das indicadoras ao longo das cadeias
        gam = 0.
        L = len(cad[0])
        if L > 1 and 0. < pj < 1.:
            r0 = pj * (1. - pj)
            for t in range(1, L):
                rt = sum(c[i] * c[i + t] for c in ind for i in range(L - t)) / (len(ind) * (L - t)) - pj * pj
                gam += 2. * (1. - t / L) * rt / r0
        d2 = (1. - pj) / (N * pj) * (1. + gam) if pj > 0. else float("inf")
        deltas2.append(d2)
        pf *= pj
        resumo.append({"nivel": j, "limiar": lim, "p": pj, "gamma": gam, "cov": math.sqrt(d2)})
        print(f"nivel {j}: limiar FS < {lim:.4f}  p = {pj:.4f}  gamma = {gam:.2f}  CoV do nivel = {100 * math.sqrt(d2):.1f}%  "
              f"P_f ate aqui = {pf:.3e}", flush=True)
        if final:
            break
        j += 1
        arqj = os.path.join(pasta, f"nivel_{j}.json")
        if os.path.exists(arqj):
            est = json.load(open(arqj))
            rng.setstate((est["rng"][0], tuple(est["rng"][1]), est["rng"][2]))
            continue
        sem = ordem[:ns]
        zs = [[gauss() for _ in range(lc - 1)] for _ in sem]   # os z no processo principal, na ordem
        a = math.sqrt(1. - rho * rho)

        def cadeia(k):
            x, v = xs[sem[k]], vs[sem[k]]
            out = [(x, v)]
            for t, z in enumerate(zs[k]):
                xc = [rho * xi + a * zi for xi, zi in zip(x, z)]
                vc = fs(xc, f"n{j}_{k}_{t}", f"{j},{k},{t + 1}")
                if vc is not None and vc < b:
                    x, v = xc, vc
                out.append((x, v))
            return out

        with cf.ThreadPoolExecutor(p["paralelos"]) as ex:
            res = list(ex.map(cadeia, range(ns)))
        est = {"xs": [s[0] for c in res for s in c], "fs": [s[1] for c in res for s in c],
               "cadeias": [[s[1] for s in c] for c in res], "nivel": j, "limiar_anterior": b}
        est["rng"] = list(rng.getstate())
        json.dump(est, open(arqj, "w"))

    cov = math.sqrt(sum(deltas2))
    nan = N + j * ns * (lc - 1)
    out = {"parametros": p, "P_f": pf, "CoV": cov, "analises": nan, "niveis": resumo,
           "analises_cov5": nan * (cov / 0.05) ** 2}
    if alfa is not None:
        out["P_f_exata"] = phi_normal(-p["linear"])
    json.dump(out, open(os.path.join(pasta, "sus.json"), "w"), indent=1)
    print(f"P_f = {pf:.4e}  CoV (Au-Beck) = {100 * cov:.1f}%  analises = {nan}  analises p/ CoV 5% = {out['analises_cov5']:.0f}"
          + (f"  (exata {out['P_f_exata']:.4e})" if alfa is not None else ""), flush=True)


if __name__ == "__main__":
    main()
