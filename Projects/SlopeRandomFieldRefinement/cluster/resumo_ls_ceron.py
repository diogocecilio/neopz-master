#!/usr/bin/env python3
"""resumo_ls_ceron.py: o resumo da campanha ls_ceron_casos.sbatch (Ceron 2025, Secoes 6.1 a 6.5, agua do artigo) por
caso, contra as Tabelas 3 a 6 do artigo: o Gamma deterministico, o FORM modal (beta e P_f = Phi(-beta)), o line
sampling (P_f, CoV, IC 95%, linhas e analises; parcial pelo linhas.csv enquanto roda) e o Monte Carlo direto pareado
(realizacoes, P_f, media e CoV de Gamma, ajuste log-normal). z: (P_f do LS - P_f do artigo) / desvio combinado.

    python3 -I resumo_ls_ceron.py <SAIDA> [casos_ls_ceron.txt] [--json <arquivo>] [--csv <arquivo>]

Sem numpy (roda no fermi). O --json grava tudo por caso (com os Gamma do MC e o p de cada linha do LS), para o
figuras_ls_ceron.py; o --csv, a tabela."""
import csv
import glob
import json
import math
import os
import sys


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def le_casos(arq):
    casos = []
    for linha in open(arq):
        t = linha.split()
        if not t or t[0].startswith("#"):
            continue
        # "-" nas colunas do artigo: um caso sem numero no artigo (os testes de casos_ls_ceron_testes.txt)
        pf, cv = num(t[2]), num(t[3])
        casos.append({"id": t[0], "secao": t[1], "pf_artigo": None if pf is None else pf / 100.,
                      "cov_pf_artigo": None if cv is None else cv / 100., "det_artigo": num(t[4]),
                      "opcoes": "" if t[5] == "-" else " ".join(t[5:])})
    return casos


def familia(c):
    if c["id"] == "ref":
        return "ref", None
    if "=" not in c["id"]:  # um teste fora das familias (ex.: alfa4_covk0)
        return c["id"], None
    k, v = c["id"].split("=", 1)
    return k, num(v)


def le_det(d):
    f = os.path.join(d, "det", "det_0_-1.csv")
    if not os.path.exists(f):
        return None
    with open(f, newline="") as fh:
        r = next(csv.DictReader(fh), None)
    if not r:
        return None
    return {"Gamma": num(r["FS_min"]), "situacao": r["situacao"]}


def le_form(d):
    f = os.path.join(d, "form", "form.json")
    if os.path.exists(f):
        j = json.load(open(f))
        return {"estado": "pronto", "beta": j["beta"], "P_f": j["P_f"], "convergiu": j["convergiu"],
                "iteracoes": j["iteracoes"], "analises": j["avaliacoes"], "alfa2_c": j["alfa2_c"],
                "alfa2_phi": j["alfa2_phi"], "alfa2_k": j["alfa2_k"]}
    trava = os.path.join(d, ".form.trava")
    if os.path.exists(trava) and "ERRO" in open(trava).read():
        return {"estado": "ERRO"}
    log = os.path.join(d, "form.log")
    if os.path.exists(log):
        its = [l.split() for l in open(log) if l[:6].strip().isdigit() and len(l.split()) >= 3]
        if its:
            return {"estado": "rodando", "it": int(its[-1][0]), "g": num(its[-1][1]), "beta": num(its[-1][2])}
        return {"estado": "rodando", "it": None}
    return {"estado": "falta"} if not os.path.exists(trava) else {"estado": "rodando", "it": None}


def le_ls(d):
    pasta = os.path.join(d, "ls")
    ps = []
    f = os.path.join(pasta, "linhas.csv")
    if os.path.exists(f):
        with open(f, newline="") as fh:
            for r in csv.DictReader(fh):
                p = num(r.get("p"))
                if p is not None:
                    ps.append((int(r["linha"]), p, int(num(r.get("analises")) or 0)))
    ps.sort()
    out = {"p": [p for _, p, _ in ps]}
    fj = os.path.join(pasta, "ls.json")
    trava = os.path.join(d, ".ls.trava")
    if os.path.exists(fj) and os.path.exists(os.path.join(d, "FEITO_LS")):
        j = json.load(open(fj))
        out.update({"estado": "pronto", "P_f": j["P_f"], "CoV": j["CoV"], "IC95": j["IC95"], "linhas": j["linhas"],
                    "perdidas": j["perdidas"], "analises": j["analises"], "motivo": j["motivo"]})
        return out
    if os.path.exists(trava) and "ERRO" in open(trava).read():
        out["estado"] = "ERRO"
    else:
        out["estado"] = "rodando" if os.path.exists(trava) else ("parcial" if ps else "falta")
    n = len(ps)
    if n >= 2:
        m = sum(p for _, p, _ in ps) / n
        v = sum((p - m) ** 2 for _, p, _ in ps) / (n - 1)
        cov = math.sqrt(v / n) / m if m > 0 else None
        out.update({"P_f": m, "CoV": cov, "linhas": n, "analises": sum(a for _, _, a in ps),
                    "IC95": [m - 1.96 * math.sqrt(v / n), m + 1.96 * math.sqrt(v / n)]})
    return out


def le_mc(d):
    gs = []
    for f in glob.glob(os.path.join(d, "mc", "mc_*_*.csv")):
        if os.path.basename(f).startswith("mc_mapa_"):
            continue
        with open(f, newline="") as fh:
            for r in csv.DictReader(fh):
                g, s = num(r.get("FS_min")), r.get("situacao")
                if s == "ok" and g is not None and math.isfinite(g):
                    gs.append((int(r["real"]), g))
                elif s == "fator" and g == 0.:  # GIM sem equilibrio nem com 1% da gravidade: falha
                    gs.append((int(r["real"]), 0.))
    gs.sort()
    g = [x for _, x in gs]
    n = len(g)
    if n < 2:
        return {"S": n, "gamma": g}
    nf = sum(1 for x in g if x < 1.)
    pf = nf / n
    media = sum(g) / n
    dp = math.sqrt(sum((x - media) ** 2 for x in g) / (n - 1))
    lg = [math.log(x) for x in g if x > 0]
    ml = sum(lg) / len(lg)
    sl = math.sqrt(sum((x - ml) ** 2 for x in lg) / (len(lg) - 1))
    return {"S": n, "n_f": nf, "P_f": pf, "CoV_Pf": math.sqrt((1 - pf) / (n * pf)) if nf else None, "media": media,
            "desvio": dp, "CoV": dp / media, "mu_ln": ml, "sigma_ln": sl, "gamma": g, "reais": [r for r, _ in gs]}


def pct(x, f="%.2f"):
    return (f % (100 * x)) if x is not None else "-"


def main():
    args = [a for a in sys.argv[1:]]
    if not args:
        sys.exit(__doc__)
    saida = args.pop(0)
    arq_json = arq_csv = None
    if "--json" in args:
        i = args.index("--json")
        arq_json = args[i + 1]
        del args[i:i + 2]
    if "--csv" in args:
        i = args.index("--csv")
        arq_csv = args[i + 1]
        del args[i:i + 2]
    arq_casos = args[0] if args else os.path.join(os.path.dirname(os.path.abspath(__file__)), "casos_ls_ceron.txt")
    casos = le_casos(arq_casos)
    tudo = []
    print("\n%-12s %4s | %7s %6s | %-12s %7s | %-8s %7s %6s %17s %5s %6s | %5s %7s %6s %6s | %7s %6s %5s"
          % ("caso", "sec", "Gdet", "artigo", "FORM", "Pf%", "LS", "Pf%", "CoV%", "IC95%", "lin", "anal",
             "MC n", "Pf%", "media", "CoV%", "artigo%", "CoV%", "z"))
    tot_an = 0
    for c in casos:
        d = os.path.join(saida, c["id"])
        fam, val = familia(c)
        det, form, ls, mc = le_det(d), le_form(d), le_ls(d), le_mc(d)
        z = None
        if ls.get("P_f") is not None and ls.get("CoV") is not None and c["pf_artigo"] is not None:
            s2 = (ls["CoV"] * ls["P_f"]) ** 2 + (c["cov_pf_artigo"] * c["pf_artigo"]) ** 2
            z = (ls["P_f"] - c["pf_artigo"]) / math.sqrt(s2) if s2 > 0 else None
        tot_an += (form.get("analises") or 0) + (ls.get("analises") or 0) + mc["S"]
        fs = form["estado"]
        if fs == "pronto":
            fs = "%.3f%s" % (form["beta"], "" if form["convergiu"] else "*")
        elif fs == "rodando" and form.get("it") is not None:
            fs = "it%d b%.2f" % (form["it"], form["beta"] or 0)
        ic = ("[%.2f;%.2f]" % (100 * ls["IC95"][0], 100 * ls["IC95"][1])) if ls.get("IC95") else "-"
        print("%-12s %4s | %7s %6s | %-12s %7s | %-8s %7s %6s %17s %5s %6s | %5d %7s %6s %6s | %7s %6s %5s"
              % (c["id"], c["secao"], ("%.4f" % det["Gamma"]) if det and det["Gamma"] is not None else "-",
                 ("%.3f" % c["det_artigo"]) if c["det_artigo"] else "-", fs,
                 pct(form.get("P_f")) if form["estado"] == "pronto" else "-", ls["estado"], pct(ls.get("P_f")),
                 pct(ls.get("CoV"), "%.1f"), ic, ls.get("linhas", "-"), ls.get("analises", "-"), mc["S"],
                 pct(mc.get("P_f")), ("%.3f" % mc["media"]) if mc.get("media") else "-", pct(mc.get("CoV"), "%.1f"),
                 pct(c["pf_artigo"]), pct(c["cov_pf_artigo"], "%.2f"), ("%+.1f" % z) if z is not None else "-"))
        tudo.append({**c, "familia": fam, "valor": val, "det": det, "form": form, "ls": ls, "mc": mc, "z": z})
    print("\nGdet: Gamma com as medias (det=1); FORM: beta (* = nao convergiu) e P_f = Phi(-beta); LS: P_f, CoV e IC 95%%"
          "\ndo line sampling (parcial pelo linhas.csv enquanto roda); MC: Monte Carlo direto pareado (as mesmas"
          "\nrealizacoes em todos os casos); artigo: P_f e CoV(P_f) das Tabelas 3 a 6; z: LS contra o artigo."
          "\nanalises do MEF ate agora (FORM + LS + MC): %d" % tot_an)
    if arq_json:
        json.dump(tudo, open(arq_json, "w"))
    if arq_csv:
        with open(arq_csv, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["caso", "secao", "Gamma_det", "Gamma_det_artigo", "beta_FORM", "Pf_FORM", "FORM_iteracoes",
                        "Pf_LS", "CoV_LS", "IC95_lo", "IC95_hi", "linhas", "analises_LS", "MC_n", "Pf_MC", "media_MC",
                        "CoV_Gamma_MC", "Pf_artigo", "CoVPf_artigo", "z_LS_artigo"])
            for t in tudo:
                f, l, m = t["form"], t["ls"], t["mc"]
                w.writerow([t["id"], t["secao"], (t["det"] or {}).get("Gamma"), t["det_artigo"],
                            f.get("beta") if f["estado"] == "pronto" else None,
                            f.get("P_f"), f.get("iteracoes"), l.get("P_f"), l.get("CoV"),
                            (l.get("IC95") or [None, None])[0], (l.get("IC95") or [None, None])[1], l.get("linhas"),
                            l.get("analises"), m["S"], m.get("P_f"), m.get("media"), m.get("CoV"), t["pf_artigo"],
                            t["cov_pf_artigo"], t["z"]])


if __name__ == "__main__":
    main()
