#!/usr/bin/env python3
"""resumo_ref.py: o resumo do Monte Carlo de referencia em Python (lote_ref.py) da coluna de Cho (2014): por caso, as
realizacoes (completas, paradas), a P_f no instante (FS_min < 1) e a de primeira passagem (FS_passagem < 1) em cada
instante gravado, com o IC de 95% de Wilson, e o valor do artigo (Fig. 10 em 18 e 24 h; 0 em 6 e 12 h, texto p. 17)
com z = (P - P_Cho) / sqrt(P(1 - P)/n + P_Cho(1 - P_Cho)/5000); mais a fracao de superficies criticas na base e a
mediana de z_crit em 12, 18 e 24 h (Figs. 5, 7 e 9) e as estatisticas do FS em 24 h (Fig. 6). Funciona com o lote
pela metade (so os blocos prontos).

    python3 -I resumo_ref.py <pasta> [casos.txt] [--json <arquivo>]

Censura: uma realizacao cujo transiente parou (situacao=transiente) conta como falha se t_passagem <= t_parou; senao,
entra nos dois limites (sem falha e com falha), impressos se houver alguma."""
import csv
import glob
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

FIG6 = {"media": 0.8376, "dp": 0.3706, "min": 0.3626, "max": 1.4806}  # Fig. 6 (caso base, 24 h)
FIG9 = {"c05_l0.5": (0.29, 0.10, 0.44), "c10_l0.5": (0.555, 0.36, 0.43), "c15_l0.5": (0.665, 0.51, 0.52),
        "c05_l2": (0.255, 0.13, 0.51), "c05_l8": (0.235, 0.15, 0.54)}  # fracao na base em 12, 18 e 24 h
FIG7 = (1.08, 1.38, 1.77)  # medianas de z_crit (caso base) em 12, 18 e 24 h


def wilson(k, n, z=1.96):
    if n == 0:
        return None
    p = k / n
    c = (p + z * z / (2 * n)) / (1 + z * z / n)
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return [max(0., c - h), min(1., c + h)]


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return math.nan


def le_casos(arq):
    casos = []
    for linha in open(arq):
        t = linha.split()
        if t and not t[0].startswith("#"):
            casos.append({"id": t[0], "Pf18": num(t[1]), "Pf24": num(t[2]), "N": int(t[3]), "opcoes": " ".join(t[4:])})
    return casos


def le_caso(pasta, cid):
    """as linhas de cada realizacao: {real: {t: linha}}"""
    reais = {}
    for arq in glob.glob(os.path.join(pasta, cid, f"{cid}_*_*.csv")):
        if "_perfil_" in arq:
            continue
        with open(arq) as a:
            for d in csv.DictReader(a):
                reais.setdefault(int(d["real"]), {})[num(d["t_h"])] = d
    return reais


def resumo_caso(c, reais):
    ts = sorted({t for r in reais.values() for t in r})
    out = {"id": c["id"], "n": len(reais), "N": c["N"], "opcoes": c["opcoes"], "instantes": {}}
    paradas = sum(1 for r in reais.values() if any(d["situacao"] != "ok" for d in r.values()))
    out["paradas"] = paradas
    for t in ts:
        ok = [r[t] for r in reais.values() if t in r and r[t]["situacao"] == "ok"]
        par = [r[t] for r in reais.values() if t in r and r[t]["situacao"] != "ok"]
        ki = sum(1 for d in ok if num(d["FS_min"]) < 1.)
        kp = sum(1 for d in ok if num(d["FS_passagem"]) < 1.)
        # paradas: falha certa se ja passou antes de parar; senao, os dois limites
        kpar_certa = sum(1 for d in par if 0. <= num(d["t_passagem"]) <= num(d["t_parou"]) + 1e-9)
        kpar_duvida = len(par) - kpar_certa
        n = len(ok) + len(par)
        base = sum(1 for d in ok if d["base"] == "1")
        zc = sorted(num(d["z_crit"]) for d in ok)
        fs = [num(d["FS_min"]) for d in ok]
        e = {"n_ok": len(ok), "n_paradas": len(par),
             "Pf_inst": ki / len(ok) if ok else None, "IC_inst": wilson(ki, len(ok)),
             "Pf_pass": (kp + kpar_certa) / n if n else None, "IC_pass": wilson(kp + kpar_certa, n),
             "Pf_pass_sup": (kp + kpar_certa + kpar_duvida) / n if n else None,
             "base": base / len(ok) if ok else None, "z_crit_mediana": zc[len(zc) // 2] if zc else None}
        if fs:
            m = sum(fs) / len(fs)
            e.update({"FS_media": m, "FS_dp": math.sqrt(sum((x - m) ** 2 for x in fs) / max(len(fs) - 1, 1)),
                      "FS_min": min(fs), "FS_max": max(fs)})
        cho = {6.: 0., 12.: 0., 18.: c["Pf18"], 24.: c["Pf24"]}.get(t)
        if cho is not None and not math.isnan(cho) and e["Pf_inst"] is not None:
            p, nn = e["Pf_inst"], len(ok)
            se = math.sqrt(p * (1 - p) / nn + cho * (1 - cho) / 5000.)
            e["Cho"], e["z"] = cho, ((p - cho) / se if se > 0 else (0. if p == cho else math.inf))
        out["instantes"][t] = e
    return out


def pct(x):
    return "   -   " if x is None else f"{100 * x:6.2f}%"


def ic(x):
    return "       -       " if x is None else f"[{100 * x[0]:5.2f}, {100 * x[1]:5.2f}]"


def main():
    args = sys.argv[1:]
    js = None
    if "--json" in args:
        i = args.index("--json")
        js = args[i + 1] if i + 1 < len(args) else sys.exit("--json <arquivo>")
        del args[i:i + 2]
    if not args:
        sys.exit(__doc__)
    pasta = args[0]
    arq_casos = args[1] if len(args) > 1 else os.path.join(pasta, "casos.txt")
    casos = le_casos(arq_casos)
    todos = []
    print(f"MC de referencia (coluna_ref.py) em {pasta}: P_f no instante e de primeira passagem, IC 95% de Wilson; Cho: "
          f"Fig. 10 (18 e 24 h; 0 em 6 e 12 h); z = (P - Cho)/EP combinado (n e 5000)")
    for c in casos:
        reais = le_caso(pasta, c["id"])
        if not reais:
            print(f"\n{c['id']}: nada ainda")
            continue
        r = resumo_caso(c, reais)
        todos.append(r)
        print(f"\n{c['id']} ({c['opcoes']}): {r['n']} de {c['N']} realizacoes" +
              (f", {r['paradas']} com o transiente parado" if r["paradas"] else ""))
        print("     t (h)   n_ok   P_f instante  IC 95%            1a passagem   IC 95%            Cho      z     "
              "base  z_crit med")
        for t, e in sorted(r["instantes"].items()):
            cho = f"{100 * e['Cho']:5.1f}% {e['z']:+6.2f}" if "Cho" in e else "      -       "
            sup = f" (ate {pct(e['Pf_pass_sup'])})" if e["n_paradas"] else ""
            print(f"   {t:6g}   {e['n_ok']:5d}   {pct(e['Pf_inst'])}    {ic(e['IC_inst'])}   {pct(e['Pf_pass'])}    "
                  f"{ic(e['IC_pass'])}{sup}  {cho}  {pct(e['base'])}  {e['z_crit_mediana'] if e['z_crit_mediana'] is not None else '-':>6}")
        e24 = r["instantes"].get(24.)
        if e24 and "FS_media" in e24:
            txt = (f"  FS(24 h): media {e24['FS_media']:.4f}, dp {e24['FS_dp']:.4f}, min {e24['FS_min']:.4f}, max "
                   f"{e24['FS_max']:.4f}")
            if c["id"] == "c05_l0.5":
                txt += f" (Fig. 6: {FIG6['media']}, {FIG6['dp']}, {FIG6['min']}, {FIG6['max']})"
            print(txt)
        if c["id"] in FIG9:
            f9 = FIG9[c["id"]]
            print("  fracao na base em 12/18/24 h: " + "/".join(pct(r["instantes"].get(t, {}).get("base")).strip()
                                                         for t in (12., 18., 24.)) +
                  f" (Fig. 9: {f9[0]}/{f9[1]}/{f9[2]})" +
                  (f"; mediana de z_crit {'/'.join(str(r['instantes'].get(t, {}).get('z_crit_mediana')) for t in (12., 18., 24.))}"
                   f" (Fig. 7: {FIG7[0]}/{FIG7[1]}/{FIG7[2]})" if c["id"] == "c05_l0.5" else ""))
    if todos:
        print("\nquadro (P_f no instante; 1a passagem entre parenteses quando difere):")
        print("  caso                       n     18 h               Cho     z      24 h               Cho     z")
        for r in todos:
            lin = f"  {r['id']:22s} {r['n']:5d}"
            for t in (18., 24.):
                e = r["instantes"].get(t)
                if not e:
                    lin += "     -" * 4
                    continue
                pp = f" ({100 * e['Pf_pass']:.2f})" if e["Pf_pass"] is not None and abs(e["Pf_pass"] - e["Pf_inst"]) > 1e-12 else ""
                cho = f"{100 * e['Cho']:5.1f}% {e['z']:+6.2f}" if "Cho" in e else "      -       "
                lin += f"   {pct(e['Pf_inst'])}{pp:9s} {cho}"
            print(lin)
    if js:
        with open(js, "w") as f:
            json.dump(todos, f, indent=1)
        print(f"\njson: {js}")


if __name__ == "__main__":
    main()
