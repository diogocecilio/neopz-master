#!/usr/bin/env python3
"""resumo_cho2014.py: o resumo da campanha cho2014_casos.sbatch (a coluna 1D de Cho 2014, SlopeRainColumn) por caso:
  - o deterministico contra a Fig. 3 (FS de 6/12/18/24 h e as frentes de 6/12/18 h, as faixas do V3.1);
  - o Monte Carlo direto: realizacoes, as que pararam (o transiente nao convergiu) e as perdidas (pararam sem travessia
    antes: a falha delas nao se sabe), a P_f no instante (FS_min < 1, a do artigo) e a de primeira passagem
    (t_passagem <= t) em cada instante, com o erro-padrao, o IC 95% de Wilson e os limites da censura (as perdidas como
    seguras e como falhas; a P_f e a das conhecidas). Uma parada conta como falha se t_passagem <= t (a travessia veio
    antes da parada; o FS so cai com chuva constante e base fechada, Secao 3 do cho2014.md);
  - z em 18 e 24 h contra a Fig. 10 (o desvio combinado com o do artigo, N = 5000) e contra a referencia em Python
    (--python: ~/projects/runs/cho2014_python, a mesma estrutura <caso>/*.csv);
  - os alvos das Figs. 5 a 9: a fracao de superficies criticas na base em 12/18/24 h (Fig. 9; 0,29 no caso base), as
    medianas de z_crit (1,08/1,38/1,77 m, Fig. 7) e a fracao acima de 0,75 m (Fig. 5: nenhuma), e o FS de 24 h (Fig. 6:
    media 0,8376, desvio 0,3706, min 0,3626, max 1,4806, pico de 0,29 na classe 0,36-0,39);
  - os FORM + line sampling (casos_cho2014_ls.txt): beta, Phi(-beta), a P_f do LS (CoV, IC 95%, linhas, analises), a
    bissecao, contra o MC do mesmo caso com a mesma falha (primeira passagem em T): z = (LS - MC) / desvio combinado
    (passa com |z| <= 2); cos(alfa_FORM, media dos xi das falhas do MC) (os xi regerados como o RandomField::Sample,
    o gerador do coluna_ref.py gerador=pz; >= 0,9 esperado), alfa^2 nos 10/20/50 primeiros termos e o custo do LS ate
    CoV 5% contra o do MC, (1 - P) / (0,05^2 P).

    python3 -I resumo_cho2014.py <SAIDA> [--casos casos_cho2014.txt] [--ls casos_cho2014_ls.txt] [--python <pasta>]
                                 [--json <arquivo>] [--curto]

--curto: so a tabela por caso e a do LS (o status_cho2014.sh). Sem numpy (roda no fermi); do resumo_hu2024.py so o
wilson."""
import csv
import glob
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from resumo_hu2024 import wilson  # noqa: E402

AQUI = os.path.dirname(os.path.abspath(__file__))
N_ARTIGO = 5000
# Fig. 3 (as faixas do V3.1): FS por instante e a frente (o pe da frente, frente_m) por instante
FIG3_FS = {6: (1.45, 1.49), 12: (1.22, 1.27), 18: (1.04, 1.08), 24: (0.35, 0.39)}
FIG3_FRENTE = {6: (0.60, 0.65), 12: (1.15, 1.25), 18: (1.75, 1.85)}
# Fig. 9: a fracao na base em 12/18/24 h
FIG9 = {"c05_l0.5": (0.29, 0.10, 0.44), "c10_l0.5": (0.555, 0.36, 0.43), "c15_l0.5": (0.665, 0.51, 0.52),
        "c05_l2": (0.255, 0.13, 0.51), "c05_l8": (0.235, 0.15, 0.54)}
FIG7_MEDIANA = {12: 1.08, 18: 1.38, 24: 1.77}  # caso base
FIG6 = {"media": 0.8376, "dp": 0.3706, "min": 0.3626, "max": 1.4806, "pico": 0.29}  # caso base, FS de 24 h
# o gerador do RandomField (KLSeparavel.h, klsep::Normal), como o coluna_ref.py
G_SEMENTE = 0x9E3779B97F4A7C15
MASK = (1 << 64) - 1


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def le_casos(arq):
    """id, Pf18, Pf24 ('-': None), N_mc e as opcoes (dict) do casos_cho2014.txt"""
    out = []
    for ln in open(arq):
        t = ln.split()
        if not t or t[0].startswith("#") or len(t) < 4:
            continue
        op = dict(x.split("=", 1) for x in t[4:] if "=" in x)
        out.append({"id": t[0], "pf18": num(t[1]), "pf24": num(t[2]), "N_mc": int(t[3]), "opcoes": op})
    return out


def le_ls_casos(arq):
    out = []
    if not os.path.exists(arq):
        return out
    for ln in open(arq):
        t = ln.split()
        if not t or t[0].startswith("#") or len(t) < 5:
            continue
        out.append({"id": t[0], "caso": t[1], "T": float(t[2]), "falha": t[3], "N_ls": int(t[4])})
    return out


def le_realizacoes(arquivos):
    """{real: {t: linha}} dos csv (do SlopeRainColumn ou do coluna_ref.py); linhas incompletas ficam de fora"""
    R = {}
    for a in arquivos:
        try:
            with open(a, newline="") as f:
                for r in csv.DictReader(f):
                    t = num(r.get("t_h"))
                    if t is None or r.get("real") is None or not r.get("KL_variancia"):
                        continue
                    R.setdefault(int(r["real"]), {})[t] = r
        except OSError:
            pass
    return R


def classifica(r, t, tipo):
    """1 falha, 0 segura, None censurada (o transiente parou antes, sem travessia) no instante t"""
    tp = num(r["t_passagem"])
    if r["situacao"] == "transiente":
        return 1 if tp is not None and 0. <= tp <= t else None
    if tipo == "passagem":
        return 1 if tp is not None and 0. <= tp <= t else 0
    if r["situacao"] == "falha":
        return 1
    return 1 if num(r["FS_min"]) < 1. else 0


def pf(R, t, tipo):
    """P_f das conhecidas, erro-padrao, Wilson, limites da censura, n, falhas e censuradas no instante t"""
    k = u = n = 0
    for linhas in R.values():
        r = linhas.get(t)
        if r is None:
            continue
        n += 1
        c = classifica(r, t, tipo)
        if c is None:
            u += 1
        else:
            k += c
    nc = n - u
    if nc == 0:
        return {"n": n, "falhas": k, "censuradas": u, "P_f": None}
    p = k / nc
    return {"n": n, "falhas": k, "censuradas": u, "P_f": p, "SE": math.sqrt(p * (1 - p) / nc), "IC95": wilson(k, nc),
            "limites": [k / n, (k + u) / n]}


def z(p, se, ref, se_ref):
    if p is None or ref is None:
        return None
    s = math.sqrt(se * se + se_ref * se_ref)
    return (p - ref) / s if s > 0 else (0. if p == ref else math.inf)


def mediana(v):
    v = sorted(v)
    if not v:
        return None
    m = len(v) // 2
    return v[m] if len(v) % 2 else 0.5 * (v[m - 1] + v[m])


def media_dp(v):
    if len(v) < 2:
        return (v[0] if v else None), None
    m = sum(v) / len(v)
    return m, math.sqrt(sum((x - m) ** 2 for x in v) / (len(v) - 1))


def alvos(R):
    """fracao na base, z_crit e o FS de 24 h (so as linhas ok)"""
    out = {}
    for t in (12., 18., 24.):
        ok = [linhas[t] for linhas in R.values() if t in linhas and linhas[t]["situacao"] == "ok"]
        if not ok:
            continue
        zc = [num(r["z_crit"]) for r in ok]
        out[f"{t:g}"] = {"base": sum(int(r["base"]) for r in ok) / len(ok), "z_crit_mediana": mediana(zc),
                         "z_crit_min": min(zc), "frac_z_menor_0.75": sum(x < 0.75 for x in zc) / len(zc), "n_ok": len(ok)}
    fs = [num(linhas[24.]["FS_min"]) for linhas in R.values() if 24. in linhas and linhas[24.]["situacao"] == "ok"]
    if fs:
        m, d = media_dp(fs)
        out["FS24"] = {"media": m, "dp": d, "min": min(fs), "max": max(fs),
                       "pico": sum(0.36 <= x < 0.39 for x in fs) / len(fs)}
    return out


def le_det(d):
    f = os.path.join(d, "det", "det_0_-1.csv")
    R = le_realizacoes([f])
    if -1 not in R:
        return None
    L = R[-1]
    out = {"FS": {f"{t:g}": num(r["FS_min"]) for t, r in sorted(L.items())},
           "frente": {f"{t:g}": num(r["frente_m"]) for t, r in sorted(L.items())},
           "t_passagem": num(L[max(L)]["t_passagem"]), "situacao": [r["situacao"] for r in L.values()]}
    ruins = [f"FS {t}" for t, (a, b) in FIG3_FS.items() if not a - 1e-9 <= (out["FS"].get(f"{t:g}") or -1) <= b + 1e-9]
    ruins += [f"frente {t}" for t, (a, b) in FIG3_FRENTE.items()
              if not a - 1e-9 <= (out["frente"].get(f"{t:g}") or -1) <= b + 1e-9]
    out["fig3"] = "ok" if not ruins else "fora: " + ", ".join(ruins)
    return out


class NormalPZ:
    """klsep::Normal (xoshiro256** semeado por splitmix64, Box-Muller), como o coluna_ref.py"""

    def __init__(s, semente):
        x = semente & MASK
        s.s = []
        for _ in range(4):
            x = (x + 0x9E3779B97F4A7C15) & MASK
            zz = x
            zz = ((zz ^ (zz >> 30)) * 0xBF58476D1CE4E5B9) & MASK
            zz = ((zz ^ (zz >> 27)) * 0x94D049BB133111EB) & MASK
            s.s.append(zz ^ (zz >> 31))
        s.tem, s.guarda = False, 0.

    def proximo(s):
        s0, s1, s2, s3 = s.s
        r1 = (s1 * 5) & MASK
        r = ((((r1 << 7) | (r1 >> 57)) & MASK) * 9) & MASK
        t = (s1 << 17) & MASK
        s2 ^= s0
        s3 ^= s1
        s1 ^= s2
        s0 ^= s3
        s2 ^= t
        s3 = ((s3 << 45) | (s3 >> 19)) & MASK
        s.s = [s0, s1, s2, s3]
        return r

    def __call__(s):
        if s.tem:
            s.tem = False
            return s.guarda
        u1 = ((s.proximo() >> 11) + 0.5) * (1.0 / 9007199254740992.0)
        u2 = ((s.proximo() >> 11) + 0.5) * (1.0 / 9007199254740992.0)
        r, t = math.sqrt(-2. * math.log(u1)), 2. * math.pi * u2
        s.guarda, s.tem = r * math.sin(t), True
        return r * math.cos(t)


def xi_realizacao(semente, n, M):
    g = NormalPZ((semente * G_SEMENTE + n) & MASK)
    return [g() for _ in range(M)]


def le_ls(d):
    """o FORM, o line sampling (pronto ou parcial pelo linhas.csv) e a bissecao de $SAIDA/ls/<id>"""
    out = {"estado": "falta"}
    tr = os.path.join(d, ".trava")
    if os.path.exists(tr) and "ERRO" in open(tr).read():
        out["estado"] = "ERRO"
    elif os.path.exists(tr):
        out["estado"] = "rodando"
    if os.path.exists(os.path.join(d, "FEITO")):
        out["estado"] = "pronto"
    fj = os.path.join(d, "form", "form.json")
    if os.path.exists(fj):
        f = json.load(open(fj))
        a2 = f["alfa2"]
        out["form"] = {"beta": f["beta"], "P_f": f["P_f"], "convergiu": f["convergiu"], "iteracoes": f["iteracoes"],
                       "avaliacoes": f["avaliacoes"], "u": f["u"],
                       "alfa2_10_20_50": [sum(a2[:m]) for m in (10, 20, 50)]}
    else:
        fp = os.path.join(d, "form", "form_parcial.json")
        if os.path.exists(fp):
            h = json.load(open(fp))["historico"]
            out["form"] = {"parcial": True, "it": h[-1]["it"], "g": h[-1]["g"], "beta": h[-1]["beta"]}
    lj = os.path.join(d, "ls", "ls.json")
    ps, na = [], 0
    fl = os.path.join(d, "ls", "linhas.csv")
    if os.path.exists(fl):
        for r in csv.DictReader(open(fl)):
            if num(r.get("p")) is not None:
                ps.append(num(r["p"]))
                na += int(num(r.get("analises")) or 0)
    if os.path.exists(lj) and os.path.exists(os.path.join(d, "FEITO_LS")):
        j = json.load(open(lj))
        out["ls"] = {k: j.get(k) for k in ("P_f", "CoV", "IC95", "linhas", "perdidas", "analises", "analises_cov5",
                                           "motivo")}
    elif len(ps) >= 2:
        m = sum(ps) / len(ps)
        sd = math.sqrt(sum((x - m) ** 2 for x in ps) / (len(ps) - 1))
        out["ls"] = {"P_f": m, "CoV": sd / math.sqrt(len(ps)) / m if m > 0 else None, "linhas": len(ps),
                     "analises": na, "parcial": True,
                     "IC95": [m - 1.96 * sd / math.sqrt(len(ps)), m + 1.96 * sd / math.sqrt(len(ps))]}
    bj = os.path.join(d, "bis", "bissecao.json")
    if os.path.exists(bj):
        b = json.load(open(bj))
        out["bis"] = {"max_abs_dc": b["max_abs_dc"], "P_f_com_bissecao": b["P_f_com_bissecao"], "P_f_ls": b["P_f_ls"],
                      "nao_cercadas": b["nao_cercadas"]}
    return out


def fmt(x, f="%.4f"):
    return "-" if x is None else f % x


def main():
    args = list(sys.argv[1:])
    if not args or args[0].startswith("--"):
        sys.exit(__doc__)
    saida = args.pop(0)
    op = {"--casos": os.path.join(AQUI, "casos_cho2014.txt"), "--ls": os.path.join(AQUI, "casos_cho2014_ls.txt"),
          "--python": None, "--json": None}
    curto = "--curto" in args
    if curto:
        args.remove("--curto")
    while args:
        k = args.pop(0)
        if k not in op or not args:
            sys.exit(f"opcao desconhecida ou sem valor: {k}")
        op[k] = args.pop(0)
    casos = le_casos(op["--casos"])
    tudo = {"casos": [], "ls": []}
    MCs = {}
    print("\n%-17s | %-7s | %5s %3s %3s | %-25s | %-25s | %6s %6s | %6s %6s | %5s"
          % ("caso", "det", "MC n", "par", "cen", "P_f 18 h (Wilson)  Cho", "P_f 24 h (Wilson)  Cho", "z18", "z24",
             "zpy18", "zpy24", "base12"))
    for c in casos:
        d = os.path.join(saida, c["id"])
        det = le_det(d)
        R = le_realizacoes(sorted(glob.glob(os.path.join(d, "mc", "mc_*.csv"))))
        MCs[c["id"]] = R
        tempos = sorted({t for L in R.values() for t in L})
        parou = sum(any(r["situacao"] == "transiente" for r in L.values()) for L in R.values())
        res = {"id": c["id"], "opcoes": c["opcoes"], "det": det, "n": len(R), "pararam": parou, "instante": {},
               "passagem": {}, "alvos": alvos(R), "python": {}}
        for t in tempos:
            res["instante"][f"{t:g}"] = pf(R, t, "instante")
            res["passagem"][f"{t:g}"] = pf(R, t, "passagem")
        res["censuradas_24"] = res["instante"].get("24", {}).get("censuradas")
        Rp = None
        if op["--python"]:
            Rp = le_realizacoes([a for a in sorted(glob.glob(os.path.join(op["--python"], c["id"], "*.csv")))
                                 if "_parametros" not in a and "_perfil_" not in a])
            for t in sorted({t for L in Rp.values() for t in L}):
                res["python"][f"{t:g}"] = {"instante": pf(Rp, t, "instante"), "passagem": pf(Rp, t, "passagem")}
        for t, ref in (("18", c["pf18"]), ("24", c["pf24"])):
            q = res["instante"].get(t)
            if q and q.get("P_f") is not None:
                q["z_cho"] = z(q["P_f"], q["SE"], ref, math.sqrt(ref * (1 - ref) / N_ARTIGO)) if ref is not None else None
                qp = res["python"].get(t, {}).get("instante")
                q["z_python"] = z(q["P_f"], q["SE"], qp["P_f"], qp["SE"]) if qp and qp.get("P_f") is not None else None
        tudo["casos"].append(res)

        def col(t, ref):
            q = res["instante"].get(t)
            if not q or q.get("P_f") is None:
                return "-"
            return "%.4f [%.3f;%.3f] %s" % (q["P_f"], q["IC95"][0], q["IC95"][1], fmt(ref, "%.3f"))
        zz = lambda t, k: fmt((res["instante"].get(t) or {}).get(k), "%+.2f")
        b12 = res["alvos"].get("12", {}).get("base")
        print("%-17s | %-7s | %5d %3d %3s | %-25s | %-25s | %6s %6s | %6s %6s | %5s"
              % (c["id"], fmt(det["FS"].get("24") if det else None), len(R), parou,
                 fmt(res["censuradas_24"], "%d"), col("18", c["pf18"]), col("24", c["pf24"]), zz("18", "z_cho"),
                 zz("24", "z_cho"), zz("18", "z_python"), zz("24", "z_python"), fmt(b12, "%.3f")))
    print("\ndet: FS de 24 h do deterministico; MC n: realizacoes; par: as que pararam (transiente); cen: censuradas em"
          "\n24 h (pararam sem travessia: fora da P_f, nos limites); P_f no instante (FS_min < 1, a do artigo) com o IC 95%"
          "\nde Wilson e a Fig. 10 (Cho); z contra o Cho (N = 5000) e contra o Python (--python); base12: fracao das"
          "\nsuperficies criticas na base em 12 h (Fig. 9: 0,29 no caso base).")

    if not curto:
        print("\nDeterministicos contra a Fig. 3 (FS 1,47/[1,22;1,27]/1,06/0,37 em 6/12/18/24 h; frentes 0,60-0,65/"
              "1,15-1,25/1,75-1,85 m):")
        for res in tudo["casos"]:
            det = res["det"]
            if det:
                print("  %-17s FS %s; frente %s; t_passagem %s h; %s"
                      % (res["id"], " ".join(f"{t}:{fmt(v)}" for t, v in det["FS"].items()),
                         " ".join(f"{t}:{fmt(v, '%.3f')}" for t, v in det["frente"].items() if float(t) <= 18),
                         fmt(det["t_passagem"]), det["fig3"]))
        print("\nP_f no instante e de primeira passagem por instante (NeoPZ; Python entre parenteses; [limites da censura]):")
        for res in tudo["casos"]:
            if not res["n"]:
                continue
            partes = []
            for t, q in res["instante"].items():
                p2 = res["passagem"][t]
                py = res["python"].get(t, {}).get("instante", {}).get("P_f")
                s = f"{t} h {fmt(q['P_f'])}" + (f" ({py:.4f})" if py is not None else "")
                if q.get("censuradas"):
                    s += " [%.4f;%.4f]" % tuple(q["limites"])
                if p2.get("P_f") is not None and q.get("P_f") is not None and abs(p2["P_f"] - q["P_f"]) > 1e-12:
                    s += f" pp {p2['P_f']:.4f}"
                partes.append(s)
            print(f"  {res['id']:17s} " + "; ".join(partes))
        print("\nAlvos das Figs. 5 a 9 (fracao na base, mediana de z_crit e fracao com z_crit < 0,75 m em 12/18/24 h;"
              " FS de 24 h):")
        for res in tudo["casos"]:
            a = res["alvos"]
            if not a:
                continue
            f9 = FIG9.get(res["id"])
            partes = []
            for i, t in enumerate(("12", "18", "24")):
                if t in a:
                    partes.append(f"{t} h base {a[t]['base']:.3f}" + (f" (Fig. 9 {f9[i]})" if f9 else "")
                                  + f" med {a[t]['z_crit_mediana']:.3f}"
                                  + (f" (Fig. 7 {FIG7_MEDIANA[int(t)]})" if res["id"] == "c05_l0.5" else "")
                                  + f" <0,75 {a[t]['frac_z_menor_0.75']:.3f}")
            print(f"  {res['id']:17s} " + "; ".join(partes))
            if "FS24" in a:
                f = a["FS24"]
                print(f"  {'':17s} FS 24 h: media {f['media']:.4f} dp {fmt(f['dp'])} min {f['min']:.4f} max "
                      f"{f['max']:.4f} pico 0,36-0,39 {f['pico']:.3f}"
                      + ("  (Fig. 6: 0,8376 0,3706 0,3626 1,4806 0,29)" if res["id"] == "c05_l0.5" else ""))

    lss = le_ls_casos(op["--ls"])
    if lss:
        print("\n%-4s %-9s %4s %-6s | %-7s %-8s %10s %5s | %-6s %10s %6s %21s %5s %6s | %10s %5s %6s | %10s %6s %6s"
              % ("id", "caso", "T", "falha", "estado", "beta", "P_FORM", "aval", "LS", "P_f", "CoV%", "IC95", "lin",
                 "anal", "MC P_f", "n", "z", "Python", "cos", "bis%"))
    for L in lss:
        d = os.path.join(saida, "ls", L["id"])
        r = le_ls(d)
        r.update(L)
        # o MC do caso com a mesma falha: primeira passagem em T (t_passagem <= T: FS_passagem(T) < 1, as duas falhas)
        R = MCs.get(L["caso"], {})
        q = pf(R, L["T"], "passagem") if R else None
        r["MC"] = q
        if op["--python"]:
            Rp = le_realizacoes([a for a in sorted(glob.glob(os.path.join(op["--python"], L["caso"], "*.csv")))
                                 if "_parametros" not in a and "_perfil_" not in a])
            r["python"] = pf(Rp, L["T"], "passagem") if Rp else None
        ls, fm = r.get("ls"), r.get("form", {})
        if ls and ls.get("P_f") is not None and q and q.get("P_f") is not None:
            se_ls = (ls["CoV"] or 0.) * ls["P_f"]
            r["z_MC"] = z(ls["P_f"], se_ls, q["P_f"], q["SE"])
            r["passa"] = abs(r["z_MC"]) <= 2.
            if ls.get("analises_cov5") and q["P_f"] > 0:
                r["custo"] = {"LS_cov5": ls["analises_cov5"] + (fm.get("avaliacoes") or 0),
                              "MC_cov5": (1 - q["P_f"]) / (0.05 ** 2 * q["P_f"])}
        # cos(alfa_FORM, media dos xi das falhas do MC)
        if fm.get("u") and R:
            u = fm["u"]
            M = len(u)
            sem = int(next((c["opcoes"].get("semente") for c in casos if c["id"] == L["caso"]), None) or 2025)
            falhas = [n for n, Lr in R.items() if L["T"] in Lr and classifica(Lr[L["T"]], L["T"], "passagem") == 1]
            if falhas:
                m = [0.] * M
                for n in falhas:
                    x = xi_realizacao(sem, n, M)
                    m = [a + b / len(falhas) for a, b in zip(m, x)]
                nu, nm = math.sqrt(sum(a * a for a in u)), math.sqrt(sum(a * a for a in m))
                r["cos_alfa_mc"] = sum(a * b for a, b in zip(u, m)) / (nu * nm) if nu > 0 and nm > 0 else None
                r["falhas_mc"] = len(falhas)
        tudo["ls"].append(dict({k: v for k, v in r.items() if k != "form"},
                               form={k: v for k, v in fm.items() if k != "u"}))
        ic = ls.get("IC95") if ls else None
        bis = r.get("bis", {})
        print("%-4s %-9s %4g %-6s | %-7s %-8s %10s %5s | %-6s %10s %6s %21s %5s %6s | %10s %5s %6s | %10s %6s %6s"
              % (L["id"], L["caso"], L["T"], L["falha"], r["estado"],
                 (fmt(fm.get("beta")) + ("" if fm.get("convergiu", True) else "*")) if fm else "-",
                 fmt(fm.get("P_f"), "%.3e"), fmt(fm.get("avaliacoes"), "%d"),
                 ("parc" if ls and ls.get("parcial") else "ok") if ls else "-", fmt(ls.get("P_f") if ls else None, "%.3e"),
                 fmt(100 * ls["CoV"] if ls and ls.get("CoV") else None, "%.1f"),
                 ("[%.2e;%.2e]" % tuple(ic)) if ic else "-", fmt(ls.get("linhas") if ls else None, "%d"),
                 fmt(ls.get("analises") if ls else None, "%d"), fmt(q.get("P_f") if q else None, "%.3e"),
                 fmt(q.get("n") if q else None, "%d"), fmt(r.get("z_MC"), "%+.2f"),
                 fmt((r.get("python") or {}).get("P_f"), "%.3e"), fmt(r.get("cos_alfa_mc"), "%.3f"),
                 fmt(100 * (bis["P_f_com_bissecao"] / bis["P_f_ls"] - 1) if bis.get("P_f_ls") else None, "%+.2f")))
        if not curto and fm.get("alfa2_10_20_50"):
            print("      alfa^2 nos 10/20/50 primeiros termos: %s; FORM %s em %d iteracoes%s"
                  % ("/".join("%.3f" % x for x in fm["alfa2_10_20_50"]), "convergiu" if fm["convergiu"] else "NAO convergiu",
                     fm["iteracoes"], "; custo ate CoV 5%%: LS %.0f analises (com o FORM), MC %.0f"
                     % (r["custo"]["LS_cov5"], r["custo"]["MC_cov5"]) if r.get("custo") else ""))
    if lss:
        print("\nbeta (* = o FORM nao convergiu) e Phi(-beta); LS: P_f, CoV, IC 95%, linhas e analises; MC: a primeira"
              "\npassagem em T no MC do caso (z = (LS - MC) / desvio combinado: passa com |z| <= 2); Python: o mesmo na"
              "\nreferencia; cos: cos(alfa_FORM, media dos xi das falhas do MC); bis%: a P_f com as raizes da bissecao"
              "\ncontra a do LS.")
    if op["--json"]:
        json.dump(tudo, open(op["--json"], "w"), indent=1)


if __name__ == "__main__":
    main()
