#!/usr/bin/env python3
# pf_v2.py <pasta> [<pasta> ...]: para cada pasta de combinacao do Monte Carlo do V2 (<SAIDA>/<COMB>), junta os
# r*/resumo.csv em mc_resumo.csv, calcula P_f = P(FS < 1) com o intervalo de Wilson de 95%, a media, o desvio e quantis
# do FS, grava a FDA empirica em fs_fda.csv e compara com a P_f de Liao et al. (2025). Realizacoes sem equilibrio
# geostatico (FS < 0,1) contam como ruptura com FS = 0; outras situacoes ficam fora da P_f. Com FS_MAX > 0 (busca
# censurada) so a P_f e valida, nao a media do FS.
import csv
import glob
import math
import os
import sys

REF = {"cho1": (0.1012, "Liao et al. 2025, Tabela 2 (MCS, 10 000 realizacoes)"),
       "cho2": (0.0012, "Liao et al. 2025, Tabela 4 (MCS, 100 000 realizacoes)")}


def wilson(k, n, z=1.96):
    if n == 0:
        return 0.0, 1.0
    p, z2 = k / n, z * z
    d = 1 + z2 / n
    c, m = (p + z2 / (2 * n)) / d, z * math.sqrt(p * (1 - p) / n + z2 / (4 * n * n)) / d
    return max(0.0, c - m), min(1.0, c + m)


def processa(pasta):
    linhas = []
    for arq in sorted(glob.glob(os.path.join(pasta, "r*", "resumo.csv"))):
        with open(arq, newline="") as f:
            linhas += list(csv.DictReader(f))
    if not linhas:
        print(f"{pasta}: nenhuma realizacao concluida")
        return
    cols = list(linhas[0].keys())
    for l in linhas:
        cols += [c for c in l if c not in cols]
    with open(os.path.join(pasta, "mc_resumo.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(linhas)
    ok = [l for l in linhas if l.get("situacao") == "ok"]
    # sem equilibrio nem com c/F0 e tan(phi)/F0 (F0 ate 0,1): FS < 0,1, uma ruptura (entra na P_f com FS = 0)
    sem = [l for l in linhas if l.get("situacao") == "sem_equilibrio_geostatico"]
    fs = sorted([float(l["FS_min"]) for l in ok] + [0.0] * len(sem))
    if not fs:
        print(f"{pasta}: nenhuma realizacao com FS")
        return
    n, k = len(fs), sum(1 for v in fs if v < 1.0)
    lo, hi = wilson(k, n)
    caso = linhas[0].get("caso", "?")
    nuc = linhas[0].get("nucleo", "?")
    censurado = any(float(l.get("FS_max") or 0) > 0 for l in ok)
    media = sum(fs) / n
    dp = math.sqrt(sum((v - media) ** 2 for v in fs) / max(1, n - 1))
    q = lambda p: fs[min(n - 1, int(p * n))]
    print(f"{pasta}: caso {caso}, nucleo {nuc}, KL_M {linhas[0].get('KL_M')}, malha {linhas[0].get('h_malha')} m")
    print(f"  realizacoes: {len(linhas)} ({len(ok)} ok, {len(sem)} sem equilibrio geostatico, contadas como ruptura; "
          f"{len(linhas) - len(ok) - len(sem)} com outra situacao, fora da P_f)")
    print(f"  P_f = {k}/{n} = {k / n:.4f}, Wilson 95% [{lo:.4f}, {hi:.4f}]")
    if caso in REF:
        r, fonte = REF[caso]
        dentro = "dentro" if lo <= r <= hi else "FORA"
        print(f"  referencia: {r:.4f} ({fonte}): {dentro} do intervalo")
    if censurado:
        print("  busca censurada (FS_MAX > 0): media e quantis do FS nao valem")
    else:
        print(f"  FS: media {media:.4f}, desvio {dp:.4f}, CoV {dp / media:.4f}, minimo {fs[0]:.4f}, "
              f"quantis 5/50/95% {q(0.05):.4f} {q(0.5):.4f} {q(0.95):.4f}, maximo {fs[-1]:.4f}")
    with open(os.path.join(pasta, "fs_fda.csv"), "w") as f:
        f.write("FS,F\n")
        for i, v in enumerate(fs):
            f.write(f"{v},{(i + 1) / n}\n")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__ or "uso: pf_v2.py <pasta> [...]")
        sys.exit(2)
    for p in sys.argv[1:]:
        processa(p)
