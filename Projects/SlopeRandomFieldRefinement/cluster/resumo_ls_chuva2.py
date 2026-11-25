#!/usr/bin/env python3
"""resumo_ls_chuva2.py: o resumo do FORM modal + line sampling da chuva2 (ls_chuva2.sbatch) por LS, contra o Monte Carlo
direto da chuva2 (job.sbatch CAMPANHA=chuva2) no MESMO evento: o menor FS nos instantes T do LS abaixo de 1, nas
realizacoes do MC com todos esses instantes ok (o MC tem tambem os outros instantes do caso; a P_f com todos eles, a da
primeira passagem do MC, sai ao lado, com a fracao das realizacoes em que o menor FS cai em T).
Por LS: o deterministico (FS das medias em cada instante do caso), o FORM (beta, P_f = Phi(-beta), iteracoes, analises,
alfa^2 de c', phi' e k), o line sampling (P_f, CoV, IC 95%, linhas, perdidas, analises, analises por linha, linhas e
analises para CoV 5%; parcial pelo linhas.csv enquanto roda; com o custos.csv, o tempo medio das analises que calculam o
transiente e das que o leem) e o MC (realizacoes, falhas, P_f, CoV). z = (P_f LS - P_f MC) / desvio combinado.

    python3 -I resumo_ls_chuva2.py <SAIDA do LS> [--mc <SAIDA do MC>] [--casos casos_chuva2.txt] [--ls casos_chuva2_ls.txt]
                                   [--json <arquivo>]

Sem numpy (roda no fermi)."""
import csv
import glob
import json
import math
import os
import sys

AQUI = os.path.dirname(os.path.abspath(__file__))


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def le_ls_txt(arq):
    itens = []
    for linha in open(arq):
        t = linha.split()
        if not t or t[0].startswith("#"):
            continue
        itens.append({"id": t[0], "caso": t[1], "T": [float(x) for x in t[2].split(",")], "n_max": int(t[3])})
    return itens


def le_det(d):
    f = os.path.join(d, "det", "det_0_-1.csv")
    if not os.path.exists(f):
        return None
    with open(f, newline="") as fh:
        rs = list(csv.DictReader(fh))
    if not rs:
        return None
    return {"FS": {float(r["t_h"]): num(r["FS_min"]) for r in rs}, "situacao": {float(r["t_h"]): r["situacao"] for r in rs}}


def le_form(d):
    f = os.path.join(d, "form", "form.json")
    trava = os.path.join(d, ".form.trava")
    if os.path.exists(f):
        j = json.load(open(f))
        h = j["historico"][-1] if j.get("historico") else {}
        return {"estado": "pronto", "beta": j["beta"], "P_f": j["P_f"], "convergiu": j["convergiu"],
                "iteracoes": j["iteracoes"], "analises": j["avaliacoes"], "alfa2_c": j["alfa2_c"],
                "alfa2_phi": j["alfa2_phi"], "alfa2_k": j["alfa2_k"], "fs_ultimo": (h.get("fs_pontos") or [None])[0],
                "reuso": j["parametros"].get("reuso", 0)}
    if os.path.exists(trava) and "ERRO" in open(trava).read():
        return {"estado": "ERRO"}
    log = os.path.join(d, "form.log")
    if os.path.exists(log):
        its = [l.split() for l in open(log) if l[:6].strip().isdigit() and len(l.split()) >= 3]
        if its:
            return {"estado": "rodando", "it": int(its[-1][0]), "g": num(its[-1][1]), "beta": num(its[-1][2])}
    return {"estado": "rodando" if os.path.exists(trava) else "falta"}


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
    fc = os.path.join(pasta, "custos.csv")
    if os.path.exists(fc):  # o tempo medio por modo (calcula/grava: com o transiente; le: so o fator)
        modos = {}
        with open(fc, newline="") as fh:
            for r in csv.DictReader(fh):
                modos.setdefault(r["modo"], []).append((num(r["s"]) or 0., num(r["s_transiente"]) or 0.,
                                                        num(r["s_fator"]) or 0.))
        out["custos"] = {m: {"n": len(v), "s": sum(x[0] for x in v) / len(v), "s_transiente": sum(x[1] for x in v) / len(v),
                             "s_fator": sum(x[2] for x in v) / len(v)} for m, v in modos.items()}
    fj = os.path.join(pasta, "ls.json")
    trava = os.path.join(d, ".ls.trava")
    if os.path.exists(fj) and os.path.exists(os.path.join(d, "FEITO_LS")):
        j = json.load(open(fj))
        out.update({"estado": "pronto", "P_f": j["P_f"], "CoV": j["CoV"], "IC95": j["IC95"], "linhas": j["linhas"],
                    "perdidas": j["perdidas"], "analises": j["analises"], "analises_linha": j["analises_linha"],
                    "analises_cov5": j["analises_cov5"], "motivo": j["motivo"], "reuso": j.get("reuso", 0)})
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
        na = sum(a for _, _, a in ps)
        out.update({"P_f": m, "CoV": cov, "linhas": n, "analises": na, "analises_linha": na / n,
                    "IC95": [m - 1.96 * math.sqrt(v / n), m + 1.96 * math.sqrt(v / n)],
                    "analises_cov5": na / n * n * (cov / 0.05) ** 2 if cov else None})
    return out


def le_mc(d, T):
    """o MC do caso em d: por realizacao, o FS de cada instante (so as linhas ok; Gamma/FS = 0 com situacao fator e
    falha); a P_f no menor FS dos instantes T (as realizacoes com todos os T ok) e com todos os instantes do caso"""
    fs = {}
    for f in glob.glob(os.path.join(d, "mc_[0-9]*_[0-9]*.csv")):
        with open(f, newline="") as fh:
            for r in csv.DictReader(fh):
                real, x, s = int(r["real"]), num(r.get("FS_min")), r.get("situacao")
                if real < 0:
                    continue
                t = float(r["t_h"])
                if s == "ok" and x is not None and math.isfinite(x):
                    fs.setdefault(real, {})[t] = x
                elif s == "fator" and x == 0.:  # sem equilibrio: falha
                    fs.setdefault(real, {})[t] = 0.
                else:
                    fs.setdefault(real, {})[t] = None
    out = {"S_total": len(fs)}
    if not fs:
        return out
    todos = sorted({t for v in fs.values() for t in v})

    def pf(ts):
        ms = [min(v[t] for t in ts) for v in fs.values() if all(v.get(t) is not None for t in ts)]
        n, nf = len(ms), sum(1 for x in ms if x < 1.)
        p = nf / n if n else None
        return {"S": n, "n_f": nf, "P_f": p,
                "CoV": math.sqrt((1 - p) / (n * p)) if n and nf else None}

    out["T"] = pf(T)
    out["todos"] = dict(pf(todos), instantes=todos)
    # em que instante cai o menor FS (as realizacoes com todos os instantes ok)
    onde = {}
    for v in fs.values():
        if all(v.get(t) is not None for t in todos):
            tm = min(todos, key=lambda t: v[t])
            onde[tm] = onde.get(tm, 0) + 1
    out["instante_do_menor"] = {str(t): c for t, c in sorted(onde.items())}
    return out


def pct(x, f="%.2f"):
    return (f % (100 * x)) if x is not None else "-"


def main():
    args = sys.argv[1:]
    if not args or args[0].startswith("--"):
        sys.exit(__doc__)
    saida = args.pop(0)
    op = {"--mc": None, "--casos": os.path.join(AQUI, "casos_chuva2.txt"),
          "--ls": os.path.join(AQUI, "casos_chuva2_ls.txt"), "--json": None}
    while args:
        k = args.pop(0)
        if k not in op or not args:
            sys.exit(__doc__)
        op[k] = args.pop(0)
    itens = le_ls_txt(op["--ls"])
    tudo = []
    print("\n%-8s %-6s %-8s | %-34s | %-12s %7s %3s %6s | %-8s %7s %5s %15s %4s %5s %5s %6s | %5s %7s %5s %7s %5s"
          % ("LS", "caso", "T (h)", "FS det (instantes)", "FORM", "Pf%", "it", "a2_k", "LS", "Pf%", "CoV%",
             "IC95%", "lin", "anal", "a/l", "a_5%", "MC n", "Pf%", "CoV%", "Pf_tod%", "z"))
    for it in itens:
        d = os.path.join(saida, it["id"])
        det, form, ls = le_det(d), le_form(d), le_ls(d)
        mc = le_mc(os.path.join(op["--mc"], it["caso"]), it["T"]) if op["--mc"] else {"S_total": 0}
        z = None
        mt = mc.get("T") or {}
        if ls.get("P_f") is not None and ls.get("CoV") and mt.get("P_f") is not None and mt.get("CoV"):
            s2 = (ls["CoV"] * ls["P_f"]) ** 2 + (mt["CoV"] * mt["P_f"]) ** 2
            z = (ls["P_f"] - mt["P_f"]) / math.sqrt(s2)
        fs = form["estado"]
        if fs == "pronto":
            fs = "b %.4f%s" % (form["beta"], "" if form["convergiu"] else "*")
        elif fs == "rodando" and form.get("it") is not None:
            fs = "it%d b%.2f" % (form["it"], form["beta"] or 0)
        dets = " ".join("%g:%.3f" % (t, x) for t, x in sorted(det["FS"].items()) if x is not None) if det else "-"
        ic = ("[%.2f;%.2f]" % (100 * ls["IC95"][0], 100 * ls["IC95"][1])) if ls.get("IC95") else "-"
        print("%-8s %-6s %-8s | %-34s | %-12s %7s %3s %6s | %-8s %7s %5s %15s %4s %5s %5s %6s | %5s %7s %5s %7s %5s"
              % (it["id"], it["caso"], ",".join("%g" % t for t in it["T"]), dets[:34], fs,
                 pct(form.get("P_f")) if form["estado"] == "pronto" else "-",
                 form.get("iteracoes", "-"), ("%.3f" % form["alfa2_k"]) if form["estado"] == "pronto" else "-",
                 ls["estado"], pct(ls.get("P_f")), pct(ls.get("CoV"), "%.1f"), ic, ls.get("linhas", "-"),
                 ls.get("analises", "-"), ("%.2f" % ls["analises_linha"]) if ls.get("analises_linha") else "-",
                 ("%.0f" % ls["analises_cov5"]) if ls.get("analises_cov5") else "-", mt.get("S", 0),
                 pct(mt.get("P_f")), pct(mt.get("CoV"), "%.1f"), pct((mc.get("todos") or {}).get("P_f")),
                 ("%+.1f" % z) if z is not None else "-"))
        if ls.get("custos"):
            print("         custo por analise (custos.csv): " + "; ".join(
                "%s %d: %.1f s (transiente %.1f s, fator %.1f s)" % (m, c["n"], c["s"], c["s_transiente"], c["s_fator"])
                for m, c in sorted(ls["custos"].items())))
        if mc.get("instante_do_menor"):
            print("         MC: o menor FS caiu em " + ", ".join("%s h: %d" % (t, c) for t, c in
                                                             mc["instante_do_menor"].items()) +
                  " (das %d realizacoes com todos os instantes %s)" % (mc["todos"]["S"],
                                                                       ",".join("%g" % t for t in mc["todos"]["instantes"])))
        tudo.append({**it, "det": det, "form": form, "ls": ls, "mc": mc, "z": z})
    print("\nFS det: o das medias (det=1) em cada instante do caso; FORM: beta (* = nao convergiu), P_f = Phi(-beta),"
          "\niteracoes (it) e alfa^2 de k (0 no modal); LS: P_f, CoV, IC 95%, linhas, analises, analises por linha (a/l)"
          "\ne analises para CoV 5% (a_5%); MC: o Monte Carlo direto da chuva2 no mesmo evento (o menor FS nos instantes T"
          "\nabaixo de 1), e Pf_tod: com todos os instantes do caso; z: LS contra o MC do mesmo evento.")
    if op["--json"]:
        json.dump(tudo, open(op["--json"], "w"), indent=1)


if __name__ == "__main__":
    main()
