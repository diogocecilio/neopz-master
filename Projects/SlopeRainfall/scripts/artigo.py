"""Cálculos do artigo sobre a intensidade crítica (docs/biot_nao_saturado/artigo): os casos que faltavam para as
seções de decomposição do viés, generalização, malha, confinamento do ar, robustez do instante da ruptura e validação.

Talude de Rahardjo et al. (2007), H_s = 10 m, alpha = 45 graus, H_w = 5 m, chuva de 24 h, salvo indicação. Modelos:
  i    completo: acoplado, ar ativo (acoplado=1 passivo=0)
  ii   acoplado, ar passivo (acoplado=1 passivo=1)
  iii  o do artigo: desacoplado, ar passivo (acoplado=0 passivo=1)
  iv   desacoplado, ar ativo (acoplado=0 passivo=0; fluxo bifásico sem deformação)
Grupos:
  D   decomposição (f50): (ii) com 54, 66 e 78 mm/h; (ii) com psi = 0 e 42, 54 mm/h; (iv) com 20, 25, 30, 36, 42 e 48 mm/h
  G   generalização: f10 em (i) com 100, 120, 140, 165 e em (iii) com 120, 140 mm/h; f50 com chi = S_e em (i) com 24, 30, 36 e
      em (iii) com 36, 42, 48, 54 mm/h
  M   malha (f50, h_sup = 0,35 m): (i) com 30 mm/h e (iii) com 48 mm/h, FS de 18 a 26 h
  T   passo de tempo (f50, (i) com psi = 0, 36 mm/h): passo de 0,125 h na chuva, instante da perda de equilíbrio
  Ar  ar também pelas laterais e pela base: G1 do item A2 (f100, alpha = 63,4 graus, 36 mm/h) e f50 (i) com 30 mm/h
  L   validação: drenagem da coluna de Liakopoulos (1965), ar ativo e passivo (acoplado; elástico)
FS de hora em hora até 30 h, tolerância do SRM 0,002 e parada no primeiro FS <= 1,01 (parar) nos casos de intensidade
crítica. A intensidade crítica I_c é a menor intensidade com FS mínimo <= 1,01 (interpolação linear entre os casos que
a cercam, inclusive os já calculados no item A1c, na Série B (iii) e no item A2).

Uso:
  python3 artigo.py rodar --exe <dir>/SlopeRainfall [--jobs 12] [--grupos D G M T Ar L]
  python3 artigo.py tabela                       # resumo dos casos
O executável SlopeRainfallPV (psi = 0) é procurado na mesma pasta; se ela tiver libpz.so, ela entra no
LD_LIBRARY_PATH (binários congelados, independentes de recompilações).
"""
import argparse
import os
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

AQUI = os.path.dirname(os.path.abspath(__file__))
MALHAS = os.path.normpath(os.path.join(AQUI, "..", "meshes"))
MODELOS = {"i": dict(acoplado=1, passivo=0), "ii": dict(acoplado=1, passivo=1), "iii": dict(acoplado=0, passivo=1),
           "iv": dict(acoplado=0, passivo=0)}
GRUPOS = ["D", "G", "M", "T", "Ar", "L"]
HORARIO = list(range(1, 31)) + [33, 36, 42, 48]
JANELA = [12, 18, 20, 21, 22, 23, 24, 25, 26]
TFS_A2 = [3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 36, 48]
# Liakopoulos (1965), propriedades de Lewis e Schrefler (1998) como em Khoei e Mohammadnejad (2011):
# k = 4,5e-13 m2 -> k_s = k rho_w g / mu_w; S_w = 1 - 0,10152 (p_c / rho_w g)^2,4279; k_rw = 1 - 2,207 (1 - S_w)^1,0121;
# k_ra de Brooks-Corey com lambda = 3 e S_r = 0,2; E = 1,3 MPa, nu = 0,4, n = 0,2975, rho_s = 2000 kg/m3
LIAK = ("E=1300 nu=0.4 c=100000 phi=30 psi=30 gamma=16.7 ks=4.4127e-6 poro=0.2975 swcc=liak:3.9748e-4,2.4279,0.2 "
        "krw=liak:2.207,1.0121 chi=sw").split()
LIAK_T = [k / 120 for k in (10, 20, 40, 60, 120, 240)]       # 5, 10, 20, 30, 60 e 120 min


def casos(grupos):
    out = []
    base = dict(solo="f50", alpha=45.0, hs=10.0, hw=5.0, train=24.0, ttotal=48.0, tfs=HORARIO, tolfs=0.002,
                parar=1.01)
    if "D" in grupos:
        out += [dict(base, grupo="D", modelo="ii", I=I) for I in (54.0, 66.0, 78.0)]
        out += [dict(base, grupo="D", modelo="ii", I=I, psi=0.0) for I in (42.0, 54.0)]
        out += [dict(base, grupo="D", modelo="iv", I=I) for I in (20.0, 25.0, 30.0, 36.0, 42.0, 48.0)]
    if "G" in grupos:
        out += [dict(base, grupo="G", modelo="i", solo="f10", I=I) for I in (100.0, 120.0, 140.0, 165.0)]
        out += [dict(base, grupo="G", modelo="iii", solo="f10", I=I) for I in (120.0, 140.0)]
        out += [dict(base, grupo="G", modelo="i", chi="se", I=I) for I in (24.0, 30.0, 36.0)]
        out += [dict(base, grupo="G", modelo="iii", chi="se", I=I) for I in (36.0, 42.0, 48.0, 54.0)]
    if "M" in grupos:
        for m, I in (("i", 30.0), ("iii", 48.0)):
            out.append(dict(base, grupo="M", modelo=m, I=I, malha=0.35, tfs=JANELA, parar=0.0))
    if "T" in grupos:
        out.append(dict(base, grupo="T", modelo="i", I=36.0, psi=0.0, dtchuva=0.125, ttotal=24.0,
                        tfs=list(range(1, 25)), parar=0.0))
    if "Ar" in grupos:
        out.append(dict(base, grupo="Ar", modelo="i", solo="f100", alpha=63.4, I=36.0, tfs=TFS_A2, parar=0.0,
                        arbase=1, arlateral=1))
        out.append(dict(base, grupo="Ar", modelo="i", I=30.0, arbase=1, arlateral=1))
    if "L" in grupos:
        for m in ("i", "ii"):
            out.append(dict(grupo="L", modelo=m, liak=True))
    return out


def rotulo(c):
    if c.get("liak"):
        return f"L_liakopoulos_{c['modelo']}"
    r = f"{c['grupo']}_{c['modelo']}_{c['solo']}_a{c['alpha']:g}_I{c['I']:g}"
    for k, f in (("psi", "psi{:g}"), ("chi", "chi{}"), ("malha", "m{:g}"), ("dtchuva", "dt{:g}")):
        if k in c:
            r += "_" + f.format(c[k])
    if c.get("arbase"):
        r += "_araberto"
    return r


def argumentos(c, exe, out):
    m = MODELOS[c["modelo"]]
    if c.get("liak"):
        arq = os.path.join(os.path.dirname(out), "coluna_liakopoulos.msh")
        if not os.path.exists(arq):
            subprocess.run([sys.executable, os.path.join(AQUI, "malha_coluna.py"), "1", "0.1", "2", "80", arq],
                           check=True, capture_output=True)
        return [exe, f"mesh={arq}", "na=0:0,0.1:0", "lateral=0", "basena=1", "arbase=1", "smax=0", "I=0", "train=2",
                "ttotal=2", "dtchuva=0.008333333333", "fs=0", "gi=0", f"passivo={m['passivo']}",
                f"acoplado={m['acoplado']}", "tsnap=" + ",".join(f"{t:.10g}" for t in LIAK_T), f"out={out}"] + LIAK
    args = [exe, f"solo={c['solo']}", f"alpha={c['alpha']:g}", f"hs={c['hs']:g}", f"hw={c['hw']:g}", f"I={c['I']:g}",
            f"train={c['train']:g}", f"ttotal={c['ttotal']:g}", f"passivo={m['passivo']}", f"acoplado={m['acoplado']}",
            "gi=0", "tfs=" + ",".join(f"{t:g}" for t in c["tfs"]), f"tolfs={c['tolfs']:g}", f"out={out}"]
    if c.get("parar"):
        args.append(f"parar={c['parar']:g}")
    for k in ("psi", "dtchuva", "arbase", "arlateral"):
        if k in c:
            args.append(f"{k}={c[k]:g}")
    if "chi" in c:
        args.append(f"chi={c['chi']}")
    if "malha" in c:
        args.append(f"mesh={os.path.join(MALHAS, f'talude_H{c['hs']:g}_a{c['alpha']:g}_hsup{c['malha']:g}.msh')}")
    return args


def rodando(out):
    return subprocess.run(["pgrep", "-f", f"out={re.escape(out)}( |$)"], capture_output=True).returncode == 0


def rodar(a):
    pasta = os.path.dirname(os.path.abspath(a.exe))
    env = dict(os.environ, MKL_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_DYNAMIC="FALSE", TBB_NUM_THREADS="1")
    if os.path.exists(os.path.join(pasta, "libpz.so")):
        env["LD_LIBRARY_PATH"] = pasta + os.pathsep + env.get("LD_LIBRARY_PATH", "")
    tarefas = []
    for c in casos(a.grupos):
        out = os.path.join(os.path.abspath(a.saida), rotulo(c))
        if (os.path.exists(os.path.join(out, "summary.txt")) and not a.refazer) or rodando(out):
            continue
        os.makedirs(out, exist_ok=True)
        exe = os.path.abspath(a.exe) + ("PV" if "psi" in c else "")
        # os mais longos primeiro: malha fina, modelo acoplado, solo f10 (chuvas fortes)
        custo = (5 if "malha" in c else 1) * (2 if MODELOS[c["modelo"]]["acoplado"] else 1) * (0 if c.get("liak") else 1)
        tarefas.append((custo, out, argumentos(c, exe, out)))
    tarefas.sort(key=lambda t: -t[0])
    print(f"{len(tarefas)} casos, {a.jobs} simultâneos", flush=True)

    def roda(t):
        _, out, args = t
        with open(os.path.join(out, "log.txt"), "w") as log:
            return out, subprocess.run(args, stdout=log, stderr=subprocess.STDOUT, env=env).returncode

    with ThreadPoolExecutor(max_workers=a.jobs) as ex:
        for f in as_completed([ex.submit(roda, t) for t in tarefas]):
            print(*f.result(), flush=True)


def ler(pasta):
    s = {}
    p = os.path.join(pasta, "summary.txt")
    if not os.path.exists(p):
        return None
    for linha in open(p):
        k, v = linha.split(None, 1)
        try:
            s[k] = float(v)
        except ValueError:
            s[k] = v.strip()
    fs = []
    for linha in list(open(os.path.join(pasta, "fs.csv")))[1:]:
        t, f = linha.split(",")[:2]
        fs.append((float(t), float(f)))
    s["fs_t"] = fs
    s["t_cruza"] = next((t for t, f in fs if f <= 1.01), -1.)
    return s


def tabela(a):
    for c in casos(GRUPOS):
        if c.get("liak"):
            continue
        s = ler(os.path.join(a.saida, rotulo(c)))
        if s is None:
            print(f"{rotulo(c):45s} --")
            continue
        print(f"{rotulo(c):45s} FS ini {s['fs_ini_srm']:.4f}  FS min {s['fs_min_srm']:.4f} ({s['t_min_srm']:g} h)  "
              f"FS<=1,01 em {s['t_cruza']:g} h  t_f {s['t_fail']:g} h  {s['minutes']:.0f} min")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("acao", choices=("rodar", "tabela"))
    p.add_argument("--exe", help="executável SlopeRainfall")
    p.add_argument("--jobs", type=int, default=12)
    p.add_argument("--grupos", nargs="+", default=GRUPOS, choices=GRUPOS)
    p.add_argument("--saida", default="resultados/artigo")
    p.add_argument("--refazer", action="store_true")
    a = p.parse_args()
    if a.acao == "rodar":
        if not a.exe or not os.access(a.exe, os.X_OK):
            raise SystemExit("indique o executável com --exe")
        rodar(a)
    else:
        tabela(a)


if __name__ == "__main__":
    main()
