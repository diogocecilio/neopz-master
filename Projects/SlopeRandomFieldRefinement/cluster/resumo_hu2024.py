#!/usr/bin/env python3
"""resumo_hu2024.py: o resumo da campanha hu2024_casos.sbatch (caso A, Hu et al. 2024) por cenario de chuva, contra a
P(FS < 1) do artigo (primeira passagem ate 72 h; Figs. 16 e 19, 0,7781 impresso): o deterministico (FS de 6 em 6 h e o
menor), o FORM modal (beta e P_f = Phi(-beta); beta < 0 quando as medias ja falham), o line sampling (P_f, CoV, IC 95%,
linhas e analises; parcial pelo linhas.csv enquanto roda) e o Monte Carlo direto (realizacoes, falhas, P_f, CoV e IC 95%
de Wilson; media e desvio do FS no fim da chuva), mais a conferencia da primeira passagem (mc13/: as primeiras realizacoes
com o FS em todos os instantes: em quantas o menor FS cai no fim da chuva, e a maior diferenca entre o FS no fim da chuva
e o menor da historia). z: (P_f do LS - P_f do MC) / desvio combinado.

    python3 -I resumo_hu2024.py <SAIDA> [casos_hu2024.txt] [--json <arquivo>] [--csv <arquivo>]

Sem numpy (roda no fermi). Um instante sem fator (o SRM sem equilibrio nem com F = 0,2: FS_min = 0, situacao fator) e
falha; outro problema (o transiente que parou) tira a realizacao da conta (e aparece em 'perd.')."""
import csv
import glob
import json
import math
import os
import re
import sys

SEM_EQUILIBRIO = 0.005  # o FS de um instante sem equilibrio (falha), como no form_kl.py


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
        opc = " ".join(t[3:])
        d = re.search(r"chuva=([0-9.]+)", opc)
        i = re.search(r"\bI=([0-9.]+)", opc)
        casos.append({"id": t[0], "pf_artigo": num(t[1]), "n_mc": int(t[2]), "opcoes": opc,
                      "D": float(d.group(1)) if d else None, "I": float(i.group(1)) if i else None})
    return casos


def le_linhas(arquivos):
    """as linhas dos csv do SlopeRainMonteCarlo: realizacao -> {t: FS} (sem equilibrio = SEM_EQUILIBRIO) e as
    realizacoes com algum instante perdido"""
    reais, perdidas = {}, set()
    for f in arquivos:
        with open(f, newline="") as fh:
            for r in csv.DictReader(fh):
                n, t, fs, s = int(r["real"]), num(r["t_h"]), num(r["FS_min"]), r["situacao"]
                if s == "ok" and fs is not None and math.isfinite(fs):
                    reais.setdefault(n, {})[t] = fs
                elif s == "fator" and fs == 0.:
                    reais.setdefault(n, {})[t] = SEM_EQUILIBRIO
                else:
                    perdidas.add(n)
    for n in perdidas:
        reais.pop(n, None)
    return reais, perdidas


def le_det(d):
    f = os.path.join(d, "det", "det_0_-1.csv")
    if not os.path.exists(f):
        return None
    reais, _ = le_linhas([f])
    fs = reais.get(-1)
    if not fs:
        return None
    t = min(fs, key=fs.get)
    return {"FS": fs, "FS_min": fs[t], "t_min": t}


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
        out.update({"P_f": m, "CoV": math.sqrt(v / n) / m if m > 0 else None, "linhas": n,
                    "analises": sum(a for _, _, a in ps), "IC95": [m - 1.96 * math.sqrt(v / n), m + 1.96 * math.sqrt(v / n)]})
    return out


def wilson(k, n, z=1.96):
    if n == 0:
        return None
    p = k / n
    c = (p + z * z / (2 * n)) / (1 + z * z / n)
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return [max(0., c - h), min(1., c + h)]


def media_dp(v):
    if len(v) < 2:
        return (v[0] if v else None), None
    m = sum(v) / len(v)
    return m, math.sqrt(sum((x - m) ** 2 for x in v) / (len(v) - 1))


def le_mc(d, sub, D):
    arqs = [f for f in glob.glob(os.path.join(d, sub, "mc_*_*.csv")) if not os.path.basename(f).startswith("mc_mapa_")]
    reais, perdidas = le_linhas(arqs)
    nts = max((len(v) for v in reais.values()), default=0)
    reais = {n: v for n, v in reais.items() if len(v) == nts}  # so as realizacoes com todos os instantes do bloco
    mins = {n: min(v.values()) for n, v in reais.items()}
    nr, nf = len(mins), sum(1 for x in mins.values() if x < 1.)
    out = {"S": nr, "n_f": nf, "perdidas": len(perdidas), "instantes": nts,
           "P_f": nf / nr if nr else None, "CoV_Pf": math.sqrt((1 - nf / nr) / nf) if nf else None,
           "IC95": wilson(nf, nr), "fs_min": [mins[n] for n in sorted(mins)]}
    if D is not None and nr:
        fd = [v[D] for v in reais.values() if D in v]
        out["media_D"], out["dp_D"] = media_dp(fd)
        # primeira passagem: em quantas o menor cai em t = D, e o quanto o FS em D fica acima do menor
        no_fim = sum(1 for v in reais.values() if D in v and v[D] <= min(v.values()) + 1e-12)
        out["argmin_D"] = no_fim
        out["dif_max"] = max((v[D] - min(v.values()) for v in reais.values() if D in v), default=None)
        # media e desvio por instante (para a Fig. 9)
        ts = sorted({t for v in reais.values() for t in v})
        out["por_instante"] = {t: media_dp([v[t] for v in reais.values() if t in v]) for t in ts}
    return out


def pct(x, f="%.2f"):
    return (f % (100 * x)) if x is not None else "-"


def main():
    args = list(sys.argv[1:])
    if not args:
        sys.exit(__doc__)
    saida = args.pop(0)
    arq_json = arq_csv = None
    for opt in ("--json", "--csv"):
        if opt in args:
            i = args.index(opt)
            if opt == "--json":
                arq_json = args[i + 1]
            else:
                arq_csv = args[i + 1]
            del args[i:i + 2]
    arq_casos = args[0] if args else os.path.join(os.path.dirname(os.path.abspath(__file__)), "casos_hu2024.txt")
    casos = le_casos(arq_casos)
    tudo = []
    print("\n%-8s | %6s %4s | %-12s %7s | %-8s %8s %6s %19s %5s %6s | %5s %4s %8s %6s %19s %6s | %8s %5s | %5s %6s"
          % ("caso", "FSdet", "t", "FORM beta", "Pf%", "LS", "Pf%", "CoV%", "IC95%", "lin", "anal", "MC n", "nf",
             "Pf%", "CoV%", "IC95% (Wilson)", "FS(D)", "artigo%", "LS/a", "z", "D/13"))
    tot = 0
    for c in casos:
        d = os.path.join(saida, c["id"])
        det, form, ls, mc, mc13 = le_det(d), le_form(d), le_ls(d), le_mc(d, "mc", c["D"]), le_mc(d, "mc13", c["D"])
        z = None
        if ls.get("P_f") is not None and ls.get("CoV") is not None and mc.get("P_f") and mc.get("CoV_Pf"):
            s2 = (ls["CoV"] * ls["P_f"]) ** 2 + (mc["CoV_Pf"] * mc["P_f"]) ** 2
            z = (ls["P_f"] - mc["P_f"]) / math.sqrt(s2) if s2 > 0 else None
        tot += (form.get("analises") or 0) + (ls.get("analises") or 0) + mc["S"] + mc13["S"]
        fs = form["estado"]
        if fs == "pronto":
            fs = "%.3f%s" % (form["beta"], "" if form["convergiu"] else "*")
        elif fs == "rodando" and form.get("it") is not None:
            fs = "it%d b%.2f" % (form["it"], form["beta"] or 0)
        ic = ("[%.3f;%.3f]" % (100 * ls["IC95"][0], 100 * ls["IC95"][1])) if ls.get("IC95") else "-"
        icm = ("[%.3f;%.3f]" % (100 * mc["IC95"][0], 100 * mc["IC95"][1])) if mc.get("IC95") else "-"
        razao = ls["P_f"] / c["pf_artigo"] if ls.get("P_f") is not None and c["pf_artigo"] else None
        conf = ("%d/%d" % (mc13["argmin_D"], mc13["S"])) if mc13.get("argmin_D") is not None else "-"
        print("%-8s | %6s %4s | %-12s %7s | %-8s %8s %6s %19s %5s %6s | %5d %4d %8s %6s %19s %6s | %8s %5s | %5s %6s"
              % (c["id"], ("%.4f" % det["FS_min"]) if det else "-", ("%g" % det["t_min"]) if det else "-", fs,
                 pct(form.get("P_f"), "%.3f") if form["estado"] == "pronto" else "-", ls["estado"],
                 pct(ls.get("P_f"), "%.4f"), pct(ls.get("CoV"), "%.1f"), ic, ls.get("linhas", "-"),
                 ls.get("analises", "-"), mc["S"], mc["n_f"], pct(mc.get("P_f"), "%.3f"), pct(mc.get("CoV_Pf"), "%.1f"),
                 icm, ("%.3f" % mc["media_D"]) if mc.get("media_D") is not None else "-", pct(c["pf_artigo"], "%.2f"),
                 ("%.2f" % razao) if razao is not None else "-", ("%+.1f" % z) if z is not None else "-", conf))
        tudo.append({**c, "det": det, "form": form, "ls": ls, "mc": mc, "mc13": mc13, "z": z, "razao_ls_artigo": razao})
    print("\nFSdet, t: o menor FS do deterministico (as medias) e o instante; FORM: beta (* = nao convergiu) e P_f = Phi(-beta);"
          "\nLS: line sampling (parcial pelo linhas.csv enquanto roda); MC: Monte Carlo direto (primeira passagem nos instantes"
          "\nda campanha; FS(D): media do FS no fim da chuva); artigo: P(FS < 1) de Hu et al. (2024); LS/a: LS / artigo; z: LS"
          "\ncontra o MC; D/13: conferencia mc13 (em quantas realizacoes o menor FS de 0 a 72 h cai no fim da chuva)."
          "\nanalises do MEF ate agora (FORM + LS + MC + conferencia): %d" % tot)
    for t in tudo:  # a media e o desvio do FS no tempo na conferencia (o I50_D24 e a Fig. 9 do artigo)
        pi = t["mc13"].get("por_instante")
        if pi:
            print("  %s, %d realizacoes (mc13): FS medio (desvio) " % (t["id"], t["mc13"]["S"]) + ", ".join(
                "%g h %.3f (%.3f)" % (k, m, s if s is not None else 0.) for k, (m, s) in sorted(pi.items())) +
                  "; maior FS(D) - FS_min = %.4f" % (t["mc13"]["dif_max"] or 0.))
    if arq_json:
        json.dump(tudo, open(arq_json, "w"))
    if arq_csv:
        with open(arq_csv, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["caso", "I", "D", "FS_det_min", "t_det_min", "beta_FORM", "Pf_FORM", "FORM_iteracoes", "Pf_LS",
                        "CoV_LS", "IC95_lo", "IC95_hi", "linhas", "analises_LS", "MC_n", "MC_falhas", "Pf_MC", "CoV_MC",
                        "MC_IC95_lo", "MC_IC95_hi", "media_FS_D", "dp_FS_D", "Pf_artigo", "LS_sobre_artigo", "z_LS_MC",
                        "mc13_n", "mc13_argmin_D", "mc13_dif_max"])
            for t in tudo:
                f, l, m, k = t["form"], t["ls"], t["mc"], t["mc13"]
                w.writerow([t["id"], t["I"], t["D"], (t["det"] or {}).get("FS_min"), (t["det"] or {}).get("t_min"),
                            f.get("beta") if f["estado"] == "pronto" else None, f.get("P_f"), f.get("iteracoes"),
                            l.get("P_f"), l.get("CoV"), (l.get("IC95") or [None, None])[0],
                            (l.get("IC95") or [None, None])[1], l.get("linhas"), l.get("analises"), m["S"], m["n_f"],
                            m.get("P_f"), m.get("CoV_Pf"), (m.get("IC95") or [None, None])[0],
                            (m.get("IC95") or [None, None])[1], m.get("media_D"), m.get("dp_D"), t["pf_artigo"],
                            t["razao_ls_artigo"], t["z"], k["S"], k.get("argmin_D"), k.get("dif_max")])


if __name__ == "__main__":
    main()
