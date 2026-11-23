#!/usr/bin/env python3
"""resumo_liu2024.py: o resumo da campanha liu2024_casos.sbatch (caso B, Liu et al. 2024) por cenario: o deterministico
(o FS das medias; o artigo da 1,525 pelo Bishop), o FORM modal (beta e Phi(-beta)), o line sampling (P_f, CoV, IC 95%;
parcial pelo linhas.csv enquanto roda) e o Monte Carlo direto (realizacoes, falhas, P_f, CoV, IC 95% de Wilson, media e
desvio do FS), contra a Tabela 4 do artigo (alfa = 0, q/k_s = 0,01; lognormais ajustadas) e a verificacao da Tabela 2
(4,0e-3 por hipercubo latino com 2500 amostras). z: (P_f do LS - P_f do MC) / desvio combinado. As leituras sao as do
resumo_hu2024.py.

    python3 -I resumo_liu2024.py <SAIDA> [casos_liu2024.txt] [--json <arquivo>]

Sem numpy (roda no fermi)."""
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import resumo_hu2024 as R  # noqa: E402

PF_TABELA2 = 4.0e-3  # LHS, 2500 amostras (este estudo e Jiang et al. 2014)
FS_ARTIGO = 1.525    # Bishop, "Entry and Exit"


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
    arq = args[0] if args else os.path.join(os.path.dirname(os.path.abspath(__file__)), "casos_liu2024.txt")
    tudo = []
    print("\n%-4s | %7s | %-12s %9s | %-8s %10s %6s %23s %5s %6s | %6s %4s %10s %6s %23s %6s %6s | %9s"
          % ("cen", "FSdet", "FORM beta", "Pf", "LS", "Pf", "CoV%", "IC95", "lin", "anal", "MC n", "nf", "Pf", "CoV%",
             "IC95 (Wilson)", "FSmed", "FSdp", "Tab.4"))
    for c in R.le_casos(arq):
        d = os.path.join(saida, c["id"])
        det, form, ls, mc = R.le_det(d), R.le_form(d), R.le_ls(d), R.le_mc(d, "mc", None)
        fm, fd = R.media_dp(mc["fs_min"])
        fs = form["estado"]
        if fs == "pronto":
            fs = "%.3f%s" % (form["beta"], "" if form["convergiu"] else "*")
        elif fs == "rodando" and form.get("it") is not None:
            fs = "it%d b%.2f" % (form["it"], form["beta"] or 0)
        e = lambda x: ("%.3e" % x) if x is not None else "-"
        ic = lambda v: ("[%.2e;%.2e]" % tuple(v)) if v else "-"
        print("%-4s | %7s | %-12s %9s | %-8s %10s %6s %23s %5s %6s | %6d %4d %10s %6s %23s %6s %6s | %9s"
              % (c["id"], ("%.4f" % det["FS_min"]) if det else "-", fs,
                 e(form.get("P_f")) if form["estado"] == "pronto" else "-", ls["estado"], e(ls.get("P_f")),
                 R.pct(ls.get("CoV"), "%.1f"), ic(ls.get("IC95")), ls.get("linhas", "-"), ls.get("analises", "-"),
                 mc["S"], mc["n_f"], e(mc.get("P_f")), R.pct(mc.get("CoV_Pf"), "%.1f"), ic(mc.get("IC95")),
                 ("%.3f" % fm) if fm is not None else "-", ("%.3f" % fd) if fd is not None else "-", e(c["pf_artigo"])))
        tudo.append({**c, "det": det, "form": form, "ls": ls, "mc": {k: v for k, v in mc.items() if k != "fs_min"},
                     "FS_media": fm, "FS_dp": fd})
    print("\nFSdet: o FS das medias (SRM; o artigo da %.3f pelo Bishop); FORM: beta (* = nao convergiu) e Phi(-beta); LS:"
          "\nline sampling; MC: Monte Carlo direto (FS < 1; FSmed e FSdp: media e desvio do FS); Tab.4: a P_f do artigo"
          "\n(lognormal ajustada a 2000 realizacoes). A verificacao da Tabela 2 (o cenario 3) da %.1e por hipercubo latino;"
          "\nas Tabelas 2 e 4 do artigo nao fecham entre si (~40 vezes)." % (FS_ARTIGO, PF_TABELA2))
    if js:
        json.dump(tudo, open(js, "w"))


if __name__ == "__main__":
    main()
