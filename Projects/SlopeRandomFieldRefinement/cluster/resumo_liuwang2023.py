#!/usr/bin/env python3
"""resumo_liuwang2023.py: o resumo da campanha liuwang2023_casos.sbatch (caso C, Liu e Wang 2023) por caso, contra a P_F1
do artigo: CRPC (o metodo semi-analitico: N efetivas, descartadas com FS0 < 1, P_F1 = media de P_1, CoV da Eq. 14, IC
95%, SRM feitos; e a mesma conta com a interpolacao da curva no lugar do ajuste a I^b + c) e Monte Carlo direto (N
efetivas, descartadas, falhas, P_F1, CoV e IC 95%). As contas sao as do crpc_liuwang.py (gpd= troca a chuva).

    python3 -I resumo_liuwang2023.py <SAIDA> [casos_liuwang2023.txt] [gpd=...] [--json <arquivo>]

Sem numpy (roda no fermi)."""
import glob
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import crpc_liuwang as C  # noqa: E402


def main():
    args = list(sys.argv[1:])
    if not args:
        sys.exit(__doc__)
    saida = args.pop(0)
    js = None
    if "--json" in args:
        i = args.index("--json")
        js = args[i + 1]
        del args[i:i + 2]
    gpd = C.GPD_PADRAO
    for a in list(args):
        if a.startswith("gpd="):
            gpd = tuple(float(x) for x in a[4:].split(","))
            args.remove(a)
    arq = args[0] if args else os.path.join(os.path.dirname(os.path.abspath(__file__)), "casos_liuwang2023.txt")
    print("chuva: xi_I %g, mu_I %g, sigma_I %g; xi_D %g, mu_D %g, sigma_D %g; omega %g" % gpd)
    tudo = []
    for l in open(arq):
        t = l.split()
        if not t or t[0].startswith("#"):
            continue
        cid, tipo, pf, n = t[0], t[1], float(t[2]), int(t[3])
        d = os.path.join(saida, cid)
        if tipo == "crpc":
            arqs = sorted(glob.glob(os.path.join(d, "crpc", "*_crpc_*.csv")))
            if not arqs:
                print("%-9s crpc: nada ainda (artigo %.2e)" % (cid, pf))
                continue
            r = C.resumo_crpc(arqs, gpd, 200.0, "potencia")
            ri = C.resumo_crpc(arqs, gpd, 200.0, "interp")
            e = lambda x: "%.3e" % x if x is not None else "-"
            print("%-9s crpc: N %d efetivas (de %d pedidas), %d descartadas, %d erros, %d SRM | P_F1 %s (CoV %s, IC95 %s) | "
                  "interp %s | artigo %.2e" % (cid, r["N"], max(n, 1), r["descartadas"], r["erros"], r["srm"],
                                               e(r.get("P_F1")), "%.3f" % r["CoV"] if "CoV" in r else "-",
                                               "[%.2e; %.2e]" % tuple(r["IC95"]) if "IC95" in r else "-",
                                               e(ri.get("P_F1")), pf))
            if n == 0:
                for x in r["realizacoes"]:
                    fit = x["ajuste"]
                    print("    solo fixo: FS0 %.4f, pontos %s, I_nf %s%s, I_lb %s" % (
                        x["FS0"], " ".join("(%g; %.1f)" % p for p in sorted(x["pontos"])), x["I_nf"],
                        ", D_C = %.4g I^%.3f + %.2f" % tuple(fit) if fit else "",
                        "%.2f" % x["I_lb"] if x["I_lb"] else "-"))
            tudo.append({"id": cid, "tipo": tipo, "pf_artigo": pf, "crpc": r, "interp": {k: v for k, v in ri.items()
                                                                                          if k != "realizacoes"}})
        else:
            r = C.resumo_mc(os.path.join(d, "mc"))
            print("%-9s mc:   %d realizacoes, %d efetivas, %d descartadas, %d perdidas, %d falhas | P_F1 %s (CoV %s, IC95 %s) "
                  "| artigo %.2e" % (cid, r["NT"], r["N"], r["descartadas"], r["perdidas"], r["falhas"],
                                     "%.3e" % r["P_F1"] if "P_F1" in r else "-",
                                     "%.3f" % r["CoV"] if "CoV" in r else "-",
                                     "[%.2e; %.2e]" % tuple(r["IC95"]) if "IC95" in r else "-", pf))
            tudo.append({"id": cid, "tipo": tipo, "pf_artigo": pf, "mc": r})
    if js:
        json.dump(tudo, open(js, "w"))


if __name__ == "__main__":
    main()
