#!/usr/bin/env python3
"""mc_la.py: Monte Carlo direto e ISMC do Cho (2010) com a analise limite do Ceron (la.py), nos MESMOS xi do Cho2010.

    python3 -I mc_la.py saida=<arquivo.csv> r=-0.5 [campo=analitico] [semente=2010] [primeira=0] [n=200000]
                        [cov=0.05] [nmin=0] [procs=4] [bloco=10] [h=0.025] [fs_ate=0] [fs_gamma=0]
                        [centro=<ustar.txt ou form.json>] [base=<base_nys2000.txt>] [cache=<psi.npy>]
                        [robusto=1.5]
    python3 -I mc_la.py resumo=1 saida=<arquivo.csv> [<outros.csv> ...]

  * realizacao n da semente s: os xi do RandomField::Sample do Cho2010 (campos_cho.xi_cpp); com centro= (ISMC) o
    deslocamento mu do FORM (os K primeiros termos de cada campo) e somado, e log_peso = sum(-xi mu + mu^2/2), o mesmo
    do C++ (RandomField::LogWeight). Entao a realizacao n daqui e a realizacao n do MEF com a mesma semente, o mesmo
    campo e o mesmo centro: o MC direto do MEF (analitico, semente 2010) e as ISMC da etapa 3b (nystrom, sementes
    21-25, centros form_r0x_<M20|M100|gim20>/ustar.txt) ficam PAREADOS realizacao a realizacao (pareado.py juntar).
  * busca do mecanismo: la.Avaliador.minimizar (rapida; a forte de novo se Gamma < robusto; robusto=99: forte em todas).
  * falha: Gamma < 1 (o mesmo que FS < 1). FS (la.Avaliador.fs, ~3 minimizacoes a mais) so nas realizacoes
    n < fs_ate ou com Gamma < fs_gamma.
  * parada: os blocos sao submetidos em ordem (no maximo 3 procs blocos em voo) e gravados na ordem em que terminam
    (o CSV nao fica ordenado por real); depois de nmin realizacoes, quando o CoV(P_f) estimado fica <= cov para de
    submeter (as em voo terminam). n e o maximo. Retoma: as realizacoes ja no CSV nao sao refeitas; uma ultima linha
    cortada (job morto no meio da gravacao) e apagada antes de continuar. Um arquivo PARAR ao lado do CSV tambem para.
  * CSV: real,semente,Gamma,FS,modo,classe,sB,xA,xB,log_peso,n_aval,descartes,t_s,fs_conv (uma linha por realizacao;
    fs_conv = 1 se a secante do FS convergiu, 0 se esgotou as iteracoes, vazio sem FS).
  * resumo: P_f, CoV, IC 95 %, Gamma medio, modos (pe/face/base), FS medio, tempo medio; varios CSV sao somados."""
from __future__ import annotations

import csv
import math
import multiprocessing as mp
import os
import sys
import time

import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
import campos_cho as cc  # noqa: E402
import la  # noqa: E402

COLS = ["real", "semente", "Gamma", "FS", "modo", "classe", "sB", "xA", "xB", "log_peso", "n_aval", "descartes", "t_s", "fs_conv"]
_G = {}


def prepara(a):
    """o avaliador (no processo pai, antes do fork: os modos e a base ficam compartilhados)"""
    campo = a.get("campo", "analitico")
    h = float(a.get("h", la.H_PADRAO))
    campos = cc.CamposCho(float(a["r"]), campo, h, base=a.get("base"), cache=a.get("cache"))
    av = la.Avaliador(campos, robusto=float(a.get("robusto", la.G_ROBUSTO)))
    mu = cc.le_desloc(a["centro"]) if a.get("centro") else None
    return av, mu


def _tarefa(args):
    semente, reais, fs_ate, fs_gamma = args
    av, mu = _G["av"], _G["mu"]
    X = cc.xi_cpp(semente, reais, av.cc.M)
    out = []
    for n, xi in zip(reais, X):
        lw = 0.0
        if mu is not None:
            xi = xi + mu
            lw = float(np.sum(-xi * mu + 0.5 * mu * mu))
        try:
            r = av.avalia(xi, fs=n < fs_ate, fs_gamma=fs_gamma)
        except Exception as err:          # noqa: BLE001  uma realizacao com erro nao derruba o job: fica fora do CSV
            import traceback
            print(f"ERRO na realizacao {n} (semente {semente}): {err!r}\n{traceback.format_exc()}", file=sys.stderr,
                  flush=True)
            continue
        out.append(dict(real=n, semente=semente, Gamma=f"{r['Gamma']:.6f}",
                        FS=f"{r['FS']:.6f}" if "FS" in r else "", modo=r["modo"], classe=r["classe"],
                        sB=f"{r['sB']:.4f}", xA=f"{r['xA']:.4f}", xB=f"{r['xB']:.4f}", log_peso=repr(lw),
                        n_aval=r["n_aval"], descartes=r["descartes"], t_s=f"{r['t']:.3f}",
                        fs_conv=("1" if r["fs_conv"] else "0") if "FS" in r else ""))
    return out


def corta_parcial(arq):
    """apaga uma ultima linha sem '\\n' (gravacao interrompida), para a retomada nao colar a proxima linha nela"""
    if not os.path.exists(arq) or os.path.getsize(arq) == 0:
        return
    with open(arq, "rb+") as f:
        f.seek(0, 2)
        n = f.tell()
        k = min(n, 1 << 16)
        f.seek(n - k)
        cauda = f.read(k)
        if cauda.endswith(b"\n"):
            return
        i = cauda.rfind(b"\n")
        f.truncate(n - k + i + 1 if i >= 0 else 0)
        print(f"  {arq}: ultima linha cortada apagada ({k - i - 1 if i >= 0 else n} bytes)", flush=True)


def ler(arqs):
    rs = {}
    for arq in arqs:
        if os.path.exists(arq):
            with open(arq, newline="") as f:
                for l in csv.DictReader(f):
                    try:
                        float(l["Gamma"]), float(l["t_s"]), float(l["log_peso"] or 0.0)
                        rs[(int(l["semente"]), int(l["real"]))] = l
                    except (ValueError, KeyError, TypeError):
                        pass                                   # linha cortada ou com campos faltando
    return list(rs.values())


def estat(linhas):
    """P_f (com peso), CoV, IC 95 %"""
    if not linhas:
        return dict(N=0, falhas=0, Pf=float("nan"), CoV=float("nan"), ic=(float("nan"),) * 2)
    # falha: Gamma < 1; com Gamma = inf (a busca nao achou mecanismo admissivel com F = 1) o FS decide
    I = np.array([float(l["Gamma"]) < 1.0 or (not math.isfinite(float(l["Gamma"])) and bool(l.get("FS"))
                  and float(l["FS"]) < 1.0) for l in linhas], float)
    w = np.exp(np.array([float(l["log_peso"] or 0.0) for l in linhas]))
    y = I * w
    N = len(y)
    p = y.mean()
    sd = y.std(ddof=1) if N > 1 else float("nan")
    cov = sd / (math.sqrt(N) * p) if p > 0 else float("inf")
    return dict(N=N, falhas=int(I.sum()), Pf=p, CoV=cov, ic=(p - 1.96 * sd / math.sqrt(N), p + 1.96 * sd / math.sqrt(N)))


def resumo(linhas, rotulo=""):
    e = estat(linhas)
    if not e["N"]:
        print(f"{rotulo}: nenhuma realizacao")
        return e
    G = np.array([float(l["Gamma"]) for l in linhas])
    ninf = int(np.sum(~np.isfinite(G)))
    G = G[np.isfinite(G)]
    t = np.array([float(l["t_s"]) for l in linhas])
    fs = np.array([float(l["FS"]) for l in linhas if l.get("FS")])
    modos = {m: sum(l["modo"] == m for l in linhas) for m in ("pe", "face", "base", "nenhum")}
    falhas = [l for l in linhas if float(l["Gamma"]) < 1.0]
    mf = {m: sum(l["modo"] == m for l in falhas) for m in ("pe", "face", "base")}
    pond = any(float(l["log_peso"] or 0) != 0 for l in linhas)
    nec = e["N"] * (e["CoV"] / 0.05) ** 2 if e["Pf"] > 0 else float("inf")
    print(f"{rotulo}: N = {e['N']}, falhas = {e['falhas']}{' (com peso: ISMC)' if pond else ''}\n"
          f"  P_f = {100 * e['Pf']:.4f} %  CoV = {100 * e['CoV']:.2f} %  IC 95 % = [{100 * e['ic'][0]:.4f}, "
          f"{100 * e['ic'][1]:.4f}] %  -> p/ CoV 5 %: {nec:.0f} realizacoes\n"
          f"  Gamma: media {G.mean():.4f}, dp {G.std(ddof=1) if len(G) > 1 else 0:.4f}, min {G.min():.4f}, "
          f"mediana {np.median(G):.4f} (sem as {ninf} com Gamma = inf: nenhum mecanismo admissivel)"
          + (f"\n  FS ({len(fs)} realizacoes): media {fs.mean():.4f}, min {fs.min():.4f}; secante do FS nao "
             f"convergida em {sum(l.get('fs_conv') == '0' for l in linhas)}" if len(fs) else "")
          + f"\n  modo (todas): " + ", ".join(f"{m} {100 * v / e['N']:.1f} %" for m, v in modos.items())
          + (f"\n  modo (falhas, sem peso): " + ", ".join(f"{m} {100 * v / len(falhas):.1f} %" for m, v in mf.items())
             if falhas else "")
          + f"\n  tempo medio {t.mean():.2f} s por realizacao (1 processo)")
    return e


def main(argv):
    a = la.le_kv(argv, campo="analitico", semente="2010", primeira="0", n="200000", cov="0.05", nmin="0",
                 procs="4", bloco="10", h=str(la.H_PADRAO), fs_ate="0", fs_gamma="0")
    saida = a.get("saida")
    if not saida:
        sys.exit(__doc__)
    if a.get("resumo"):
        extra = [x for x in argv if "=" not in x]
        resumo(ler([saida] + extra), saida)
        return
    if "r" not in a:
        sys.exit("falta r=")
    semente, primeira, nmax = int(a["semente"]), int(a["primeira"]), int(a["n"])
    alvo, nmin, procs, bloco = float(a["cov"]), int(a["nmin"]), int(a["procs"]), int(a["bloco"])
    fs_ate, fs_gamma = int(a["fs_ate"]), float(a["fs_gamma"])
    os.makedirs(os.path.dirname(os.path.abspath(saida)), exist_ok=True)
    corta_parcial(saida)
    feitas = ler([saida])
    ja = {int(l["real"]) for l in feitas if int(l["semente"]) == semente}
    novo = not os.path.exists(saida) or os.path.getsize(saida) == 0
    t0 = time.perf_counter()
    av, mu = prepara(a)
    print(f"mc_la: r = {a['r']}, campo {a['campo']} (variancia retida {av.cc.fracao:.4f}), semente {semente}, "
          f"realizacoes {primeira}..{primeira + nmax - 1}, {'ISMC centro ' + a['centro'] if mu is not None else 'direto'}, "
          f"cov alvo {alvo} (nmin {nmin}), {procs} processos, ja feitas {len(ja)}; preparo {time.perf_counter() - t0:.1f} s",
          flush=True)
    _G.update(av=av, mu=mu)
    pend = [n for n in range(primeira, primeira + nmax) if n not in ja]
    blocos = [pend[i:i + bloco] for i in range(0, len(pend), bloco)]
    parar_arq = os.path.join(os.path.dirname(os.path.abspath(saida)), "PARAR")
    f = open(saida, "a", newline="")
    w = csv.DictWriter(f, fieldnames=COLS)
    if novo:
        w.writeheader()
    linhas = [l for l in feitas if int(l["semente"]) == semente]
    ultimo = time.time()

    def chega():
        if os.path.exists(parar_arq):
            return "PARAR"
        if alvo > 0 and len(linhas) >= max(nmin, 100):
            e = estat(linhas)
            if e["falhas"] >= 10 and e["CoV"] <= alvo:
                return f"CoV {100 * e['CoV']:.2f} % <= {100 * alvo:g} %"
        return ""

    motivo = chega() or ""
    with mp.get_context("fork").Pool(procs) as pool:
        voo = []
        i = 0
        while (i < len(blocos) and not motivo) or voo:
            while not motivo and i < len(blocos) and len(voo) < 3 * procs:
                voo.append(pool.apply_async(_tarefa, ((semente, blocos[i], fs_ate, fs_gamma),)))
                i += 1
            pronto = [v for v in voo if v.ready()]
            if not pronto:
                time.sleep(0.05)
                continue
            for v in pronto:
                voo.remove(v)
                for l in v.get():
                    w.writerow(l)
                    linhas.append({k: str(x) for k, x in l.items()})
            f.flush()
            if not motivo:
                motivo = chega()
            if time.time() - ultimo > 60 or (motivo and not voo):
                e = estat(linhas)
                print(f"  {time.strftime('%H:%M:%S')} N = {e['N']}, falhas {e['falhas']}, P_f = {100 * e['Pf']:.4f} %, "
                      f"CoV = {100 * e['CoV']:.2f} %{', parando: ' + motivo if motivo else ''}", flush=True)
                ultimo = time.time()
    f.close()
    print(f"fim ({motivo or 'n maximo'}), {time.perf_counter() - t0:.0f} s", flush=True)
    resumo(linhas, saida)


if __name__ == "__main__":
    main(sys.argv[1:])
