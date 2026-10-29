"""Roda os casos dos benchmarks com a literatura e do patamar inferior do FS (figuras e tabelas: benchmarks.py).

Grupos
  zhang    Zhang, Fredlund, Fredlund & Wilson (2014): casos estáticos A e B, FS inicial com as leis de chi (8 casos)
  patamar  superfície alagada (I -> infinito) até o regime permanente, modelo do artigo de Rahardjo (6 casos)
  tang     Tang, Huang, Sheng & Sloan (2018): chuvas uniformes, adiantada e atrasada (12 casos)
  pradhan  Pradhan et al. (2022): corte de estrada e encosta sem o corte, chuva de 1 a 12/07/2018 (2 casos)

Uso, na pasta dos resultados:
  python3 <neopz>/Projects/SlopeRainfall/scripts/rodar_benchmarks.py rodar --exe <SlopeRainfall> [--jobs 4] [grupos]
  python3 <neopz>/Projects/SlopeRainfall/scripts/rodar_benchmarks.py status
  python3 <neopz>/Projects/SlopeRainfall/scripts/rodar_benchmarks.py lista [grupos]

Os resultados ficam em benchmarks/<grupo>/<caso> (log em benchmarks/<grupo>/<caso>.log). Casos com summary.txt são
pulados: depois de uma interrupção basta repetir o comando. Cada caso usa --threads núcleos (padrão 1) e --jobs casos
rodam ao mesmo tempo, os mais longos primeiro. Para continuar com o terminal fechado: nohup python3 ... &
Figuras: python3 <neopz>/Projects/SlopeRainfall/scripts/benchmarks.py --runs benchmarks --figuras figuras
"""
import argparse
import os
import re
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

P = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
GRUPOS = ("pradhan", "patamar", "tang", "zhang")


def zhang():
    a = (f"mesh={P}/meshes/zhang_A.msh na=0:14,10:12.4,20:10.2,30:7.6,40:4.6,50:3.4,60:1.7,70:0 "
         "gamma=18 c=10 phi=34").split()
    b = f"mesh={P}/meshes/zhang_B.msh na=0:22,15:21,20:19.5,30:17.3,40:14.5,55:9,70:6.3 gamma=18 c=10 phi=24".split()
    com = "I=0 smax=1e6 acoplado=0 passivo=1 gi=0 ttotal=0.01 tfs=1".split()
    return [("A_chi0", com + a + ["chi=0"]), ("A_phib15", com + a + ["chi=0.3973"]), ("A_chi1", com + a + ["chi=1"]),
            ("A_sw_swcc1", com + a + ["chi=sw", "swcc=fx:1,2,1"]), ("B_chi0", com + b + ["chi=0"]),
            ("B_phib15", com + b + ["chi=0.6018"]), ("B_chi1", com + b + ["chi=1"]),
            ("A_chi1_smax75", com + a + ["chi=1", "smax=75"])]


def patamar():
    com = "mecanica=0 alpha=45 hs=10 hw=5 I=0 alagado=1 acoplado=0 passivo=1 gi=0".split()
    t768 = "train=768 ttotal=768 dtmax=24 tfs=3,6,12,24,48,96,192,384,768".split()
    return [("f100_iso", com + "solo=f100 train=6144 ttotal=6144 dtmax=96 tfs=6,24,48,96,192,384,768,1536,3072,6144"
             .split()),
            ("f50_iso", com + "solo=f50 train=1536 ttotal=1536 dtmax=24 tfs=3,6,12,24,48,96,192,384,768,1536".split()),
            ("f50_th135", com + "solo=f50 kx=1e-4 ky=1e-5 theta=135".split() + t768),
            ("f50_th0", com + "solo=f50 kx=1e-4 ky=1e-5 theta=0".split() + t768),
            ("f50_th45", com + "solo=f50 kx=1e-4 ky=1e-5 theta=45".split() + t768),
            ("f10_iso", com + "solo=f10 train=384 ttotal=384 dtmax=12 tfs=1,3,6,12,24,48,96,192,384".split())]


def tang():
    """talude de Rahardjo com o solo f100, k_s = 2e-6 m/s, c' = 0, phi^b = 17,3 graus (chi = 0,6387 na sucção),
    sucção inicial limitada a 50 kPa; percolação desacoplada e p_a = 0 (SEEP/W); FS nos instantes da grade do
    artigo (0,25 h na chuva, 1 h nas 24 h seguintes, 4 h depois)"""
    def grade(train, ttotal):
        t, ts = 0.0, []
        while t < ttotal - 1e-9:
            dt = 0.25 if t < train - 1e-9 else (1.0 if t < train + 24 - 1e-9 else 4.0)
            dt = min(dt, ttotal - t)
            if t < train - 1e-9:
                dt = min(dt, train - t)
            t = round(t + dt, 6)
            ts.append(t)
        return ts

    def tfs(train, ttotal):
        g = grade(train, ttotal)
        want = [x for x in (4, 8, 12, 16, 20, 24) if x < train] + [train] + \
               [train + d for d in (2, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40)] + [ttotal]
        out = []
        for w in want:
            best = min(g, key=lambda v: abs(v - w))
            if best not in out and best <= ttotal:
                out.append(best)
        return ",".join(f"{v:g}" for v in sorted(out))

    base = "mecanica=0 solo=f100 ks=2e-6 c=0 chi=0.6387 smax=50 acoplado=0 passivo=1 gi=0".split()
    casos = [("controle", 0.0, 25.6, "uniforme")]
    casos += [(f"u25.6_I{I:g}", I, 25.6, "uniforme") for I in (2.0, 2.4, 2.7, 3.0, 3.5)]
    casos += [(f"u6.4_I{I:g}", I, 6.4, "uniforme") for I in (8.0, 10.0, 12.0)]
    casos += [(f"{p}25.6_I2.64", 2.64, 25.6, p) for p in ("uniforme", "adiantado", "atrasado")]
    out = []
    for nome, I, dr, pad in casos:
        ttotal = round(2.5 * dr, 6)
        out.append((nome, base + [f"I={I:g}", f"train={dr:g}", f"ttotal={ttotal:g}", f"padrao={pad}",
                                  f"tfs={tfs(dr, ttotal)}"]))
    return out


def pradhan():
    com = (f"mecanica=0 dtchuva=0.5 camadas={P}/pradhan/camadas.txt na=0:13,34:27,69:27 lateral=0 basena=1 smax=300 "
           f"acoplado=0 passivo=1 gi=0 regime=0.1188 chuva={P}/pradhan/chuva_julho2018.csv evap=4 psimin=294.3 "
           "ttotal=271 tfs=34,42,96,101,167,173,184,200,216,222,261,264,267,271").split()
    return [("corte", com + [f"mesh={P}/meshes/pradhan_corte.msh"]),
            ("semcorte", com + [f"mesh={P}/meshes/pradhan_semcorte.msh"])]


def casos(grupos, raiz):
    fun = {"zhang": zhang, "patamar": patamar, "tang": tang, "pradhan": pradhan}
    return [(g, nome, args, os.path.join(raiz, g)) for g in GRUPOS if g in grupos for nome, args in fun[g]()]


def roda(c, exe, env):
    g, nome, args, pasta = c
    if os.path.exists(os.path.join(pasta, nome, "summary.txt")):
        return g, nome, "já pronto"
    os.makedirs(pasta, exist_ok=True)
    t0 = time.time()
    with open(os.path.join(pasta, nome + ".log"), "w") as log:
        r = subprocess.run([exe] + args + [f"out={nome}"], cwd=pasta, stdout=log, stderr=subprocess.STDOUT, env=env)
    ok = os.path.exists(os.path.join(pasta, nome, "summary.txt"))
    dur = f"{(time.time() - t0) / 60:.1f} min"
    return g, nome, f"pronto em {dur}" if ok and r.returncode == 0 else \
        f"falhou em {dur} (código {r.returncode}, veja {nome}.log)"


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("acao", choices=("rodar", "status", "lista"))
    p.add_argument("grupos", nargs="*", help=f"{', '.join(GRUPOS)} (padrão: todos)")
    p.add_argument("--exe", help="executável SlopeRainfall")
    p.add_argument("--jobs", type=int, default=4, help="casos simultâneos")
    p.add_argument("--threads", type=int, default=1, help="núcleos por caso (MKL/OpenMP)")
    p.add_argument("--runs", default="benchmarks", help="pasta dos resultados")
    a = p.parse_args()
    grupos = a.grupos or GRUPOS
    if set(grupos) - set(GRUPOS):
        p.error(f"grupos válidos: {', '.join(GRUPOS)}")
    cs = casos(grupos, os.path.abspath(a.runs))
    if a.acao == "lista":
        for g, nome, args, pasta in cs:
            print(f"cd {pasta} && {a.exe or 'SlopeRainfall'} {' '.join(args)} out={nome}")
    elif a.acao == "status":
        em = subprocess.run(["pgrep", "-fa", "SlopeRainfall "], capture_output=True, text=True).stdout
        rodando = set(re.findall(r"out=(\S+)", em))
        for g in grupos:
            gc = [c for c in cs if c[0] == g]
            prontos = [c[1] for c in gc if os.path.exists(os.path.join(c[3], c[1], "summary.txt"))]
            emx = [c[1] for c in gc if c[1] in rodando and c[1] not in prontos]
            print(f"{g:8s} {len(prontos):2d} de {len(gc):2d} prontos" + (f"; rodando: {', '.join(emx)}" if emx else ""))
    else:
        if not a.exe or not os.access(a.exe, os.X_OK):
            sys.exit("indique o executável com --exe <build>/Projects/SlopeRainfall/SlopeRainfall")
        n = str(a.threads)
        env = dict(os.environ, MKL_NUM_THREADS=n, OMP_NUM_THREADS=n, MKL_DYNAMIC="FALSE")
        exe = os.path.abspath(a.exe)
        print(f"{len(cs)} casos, {a.jobs} simultâneos, {a.threads} núcleo(s) por caso; resultados em "
              f"{os.path.abspath(a.runs)}", flush=True)
        # submitted in the order of the list (longest first) and reported as each case ends
        with ThreadPoolExecutor(max_workers=a.jobs) as ex:
            for f in as_completed([ex.submit(roda, c, exe, env) for c in cs]):
                g, nome, st = f.result()
                print(f"{g}/{nome}: {st}", flush=True)


if __name__ == "__main__":
    main()
