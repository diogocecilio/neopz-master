"""Roda as séries paramétricas de Rahardjo et al. (2007) com o executável SlopeRainfall (NeoPZ), em paralelo.

Séries (Tabela 2 do artigo):
  A  solos f10, f50, f100 x alpha 26,6/33,7/45/63,4 graus x I_r 9, 80 mm/h e 1 k_s; Hs = 10 m, Hw = 5 m; até 240 h
  B  alpha = 45 graus; intensidades em múltiplos de k_s para cada solo; até 48 h
  C  f50, alpha = 45 graus, Hs = 10 m; Hw = 2,5/5/7,5/10/15 m x I_r 9, 36, 80 mm/h; até 48 h
  D  f50, alpha = 45 graus, Hw = 5 m; Hs = 5/10/20/40 m x I_r 9, 36, 80 mm/h; até 48 h
Cada caso grava resultados/<série>/<caso>/ (fs.csv, hist.csv, snaps_*.csv, summary.txt...); casos com summary.txt
são pulados (--refazer recalcula).

Exemplos:
  python rodar_series.py --serie B --listar
  python rodar_series.py --serie B --jobs 4 --exe ../../../build/Projects/SlopeRainfall/SlopeRainfall
  python rodar_series.py --serie todas --jobs 4 --sem-gi        # só redução de resistência (mais rápido)
  python rodar_series.py --serie C --chuvas 36 --rapido          # poucos instantes de FS
Depois: python figuras.py --serie todas ; python animacao.py resultados/B/f50_a45_H10_Hw5_I36
"""
import argparse
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

ALPHAS = (26.6, 33.7, 45.0, 63.4)
KS_MMH = {"f10": 360.0, "f50": 36.0, "f100": 3.6}
SERIE_B = {"f10": (9, 18, 36, 90, 180, 360, 900), "f50": (3.6, 9, 18, 36, 54, 90, 180),
           "f100": (0.9, 1.8, 3.6, 5.4, 9, 18, 36)}
T_LONGO = (6, 12, 18, 24, 30, 36, 48, 72, 96, 144, 240)
T_CURTO = (3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 36, 48)
T_RAPIDO = (6, 12, 18, 24, 36, 48)


def caso(solo="f50", alpha=45.0, hs=10.0, hw=5.0, I=36.0, ttotal=48.0, tfs=T_CURTO):
    return dict(solo=solo, alpha=alpha, hs=hs, hw=hw, I=float(I), ttotal=ttotal, tfs=tfs)


def rotulo(c):
    return f"{c['solo']}_a{c['alpha']:g}_H{c['hs']:g}_Hw{c['hw']:g}_I{c['I']:g}"


def casos_da_serie(serie):
    out = []
    if serie == "A":
        for solo in ("f10", "f50", "f100"):
            for a in ALPHAS:
                for I in sorted({9.0, 80.0, KS_MMH[solo]}):
                    out.append(caso(solo=solo, alpha=a, I=I, ttotal=240.0, tfs=T_LONGO))
    elif serie == "B":
        for solo, Is in SERIE_B.items():
            for I in Is:
                out.append(caso(solo=solo, I=I))
    elif serie == "C":
        for hw in (2.5, 5.0, 7.5, 10.0, 15.0):
            for I in (9.0, 36.0, 80.0):
                out.append(caso(hw=hw, I=I))
    elif serie == "D":
        for hs in (5.0, 10.0, 20.0, 40.0):
            for I in (9.0, 36.0, 80.0):
                out.append(caso(hs=hs, I=I))
    else:
        raise ValueError(serie)
    return out


def main():
    aqui = os.path.dirname(os.path.abspath(__file__))
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--serie", default="B", help="A, B, C, D ou todas")
    p.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2)))
    p.add_argument("--exe", default=os.path.join(aqui, "..", "..", "..", "..", "build", "Projects", "SlopeRainfall",
                                                  "SlopeRainfall"))
    p.add_argument("--saida", default="resultados")
    p.add_argument("--chuvas", nargs="*", type=float, help="filtra as intensidades (mm/h)")
    p.add_argument("--solos", nargs="*", help="filtra os solos (f10 f50 f100)")
    p.add_argument("--passivo", action="store_true", help="p_a = 0 em todo o domínio (como o SEEP/W)")
    p.add_argument("--desacoplado", action="store_true",
                   help="percolação sem o acoplamento hidromecânico nas equações de massa (como o SEEP/W)")
    p.add_argument("--sem-gi", action="store_true", help="não calcula o FS por aumento da gravidade")
    p.add_argument("--rapido", action="store_true", help="poucos instantes de FS")
    p.add_argument("--vtk", action="store_true", help="grava VTK nos instantes do FS")
    p.add_argument("--refazer", action="store_true")
    p.add_argument("--listar", action="store_true")
    a = p.parse_args()
    series = ["A", "B", "C", "D"] if a.serie == "todas" else [a.serie]
    tarefas = []
    for s in series:
        pasta = os.path.join(a.saida, s + ("_passivo" if a.passivo else "") + ("_desacoplado" if a.desacoplado else ""))
        for c in casos_da_serie(s):
            if a.chuvas and c["I"] not in a.chuvas:
                continue
            if a.solos and c["solo"] not in a.solos:
                continue
            if a.rapido:
                c["tfs"] = tuple(t for t in T_RAPIDO if t <= c["ttotal"]) + ((c["ttotal"],) if c["ttotal"] > 48 else ())
            out = os.path.join(pasta, rotulo(c))
            if not a.refazer and os.path.exists(os.path.join(out, "summary.txt")):
                continue
            args = [a.exe] + [f"{k}={c[k]:g}" if k != "solo" else f"solo={c[k]}" for k in
                              ("solo", "alpha", "hs", "hw", "I", "ttotal")]
            args += ["tfs=" + ",".join(f"{t:g}" for t in c["tfs"]), f"out={out}", f"passivo={int(a.passivo)}", f"acoplado={0 if a.desacoplado else 1}",
                     f"gi={0 if a.sem_gi else 1}", f"vtk={int(a.vtk)}"]
            tarefas.append((out, args))
    print(f"{len(tarefas)} casos com {a.jobs} processo(s)")
    for out, args in tarefas:
        print("  ", " ".join(args[1:]))
    if a.listar or not tarefas:
        return
    for out, _ in tarefas:
        os.makedirs(out, exist_ok=True)

    def roda(t):
        out, args = t
        t0 = time.time()
        with open(os.path.join(out, "log.txt"), "w") as log:
            r = subprocess.run(args, stdout=log, stderr=subprocess.STDOUT)
        return out, r.returncode, (time.time() - t0) / 60

    falhas = []
    with ThreadPoolExecutor(max_workers=a.jobs) as ex:
        for f in as_completed([ex.submit(roda, t) for t in tarefas]):
            out, rc, mins = f.result()
            ok = rc == 0 and os.path.exists(os.path.join(out, "summary.txt"))
            print(f"{'ok' if ok else 'FALHOU'}  {out}  ({mins:.1f} min)", flush=True)
            if not ok:
                falhas.append(out)
    print(f"{len(falhas)} falha(s)")
    sys.exit(1 if falhas else 0)


if __name__ == "__main__":
    main()
