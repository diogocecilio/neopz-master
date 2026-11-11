#!/usr/bin/env python3
"""tabelas_comparacao.py [csv]: as tabelas da comparacao refinada x malhas fixas do relatorio.

Le dados/ref_comparacao.csv (copia do tabela.csv de 'compara.py tabela <pasta> --csv'; 9 casos: det e realizacoes
0 a 7; 5 metodos) e escreve, ao lado deste script:
  tabelas/ref_comparacao.tex        por metodo: triangulos, |G/G_foco125 - 1| medio e maximo, a diferenca media com
                                    sinal, os tempos medios de parede e de CPU por caso e a razao deles para a foco125;
  tabelas/ref_comparacao_casos.tex  G de cada caso e metodo, com a diferenca para a foco125.
So o tabular (com booktabs); a legenda fica na secao. Virgula decimal.
"""
import csv
import os
import sys

AQUI = os.path.dirname(os.path.abspath(__file__))
METODOS = ["foco375", "foco25", "foco1875", "refinada", "foco125"]
NOMES = {"foco375": r"\texttt{foco375}", "foco25": r"\texttt{foco25}", "foco1875": r"\texttt{foco1875}",
         "refinada": r"\texttt{refinada}", "foco125": r"\texttt{foco125}"}
CASOS = ["det"] + [str(i) for i in range(8)]


def br(x, nd):
    """numero com virgula decimal; o sinal de menos em modo matematico"""
    s = "%.*f" % (nd, x)
    if s.startswith("-"):
        s = "$-$" + s[1:]
    return s.replace(".", "{,}")


def pct(x, nd=2, sinal=False):
    s = "%+.*f" % (nd, 100.0 * x) if sinal else "%.*f" % (nd, 100.0 * x)
    if s.startswith("-"):
        s = "$-$" + s[1:]
    elif s.startswith("+"):
        s = "$+$" + s[1:]
    return s.replace(".", "{,}") + r"\,\%"


def media(v):
    return sum(v) / len(v)


def main():
    arq = sys.argv[1] if len(sys.argv) > 1 else os.path.join(AQUI, "dados", "ref_comparacao.csv")
    with open(arq) as f:
        linhas = list(csv.DictReader(f))
    d = {(r["metodo"], r["caso"]): r for r in linhas}
    falta = [(m, c) for m in METODOS for c in CASOS if (m, c) not in d]
    if falta:
        raise SystemExit("faltam casos: %s" % falta)
    ref = {c: float(d[("foco125", c)]["fator"]) for c in CASOS}
    dif = {(m, c): float(d[(m, c)]["fator"]) / ref[c] - 1.0 for m in METODOS for c in CASOS}
    parede125 = media([float(d[("foco125", c)]["parede_s"]) for c in CASOS])
    cpu125 = media([float(d[("foco125", c)]["cpu_s"]) for c in CASOS])

    # por metodo
    out = [r"\begin{tabular}{@{}lrrrrrrrrr@{}}", r"\toprule",
           r" & & & \multicolumn{3}{c}{$\Gamma/\Gamma_{\mathrm{foco125}}-1$}"
           r" & \multicolumn{2}{c}{tempo por caso (s)} & \multicolumn{2}{c}{razão p/ \texttt{foco125}} \\",
           r"\cmidrule(lr){4-6}\cmidrule(lr){7-8}\cmidrule(l){9-10}",
           r"método & triâng. & equações & $|\cdot|$ médio & $|\cdot|$ máx. & média c/ sinal"
           r" & parede & CPU & parede & CPU \\", r"\midrule"]
    resumo = []
    for m in METODOS:
        tri = media([float(d[(m, c)]["triangulos"]) for c in CASOS])
        neq = media([float(d[(m, c)]["neq"]) for c in CASOS])
        ad = [abs(dif[(m, c)]) for c in CASOS]
        sd = media([dif[(m, c)] for c in CASOS])
        par = media([float(d[(m, c)]["parede_s"]) for c in CASOS])
        cpu = media([float(d[(m, c)]["cpu_s"]) for c in CASOS])
        difs = ["---"] * 3 if m == "foco125" else [pct(media(ad)), pct(max(ad)), pct(sd, sinal=True)]
        out.append("%s & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\" % (
            NOMES[m], "%.0f" % tri, "%.0f" % neq, difs[0], difs[1], difs[2],
            br(par, 1), br(cpu, 1), br(par / parede125, 2), br(cpu / cpu125, 2)))
        resumo.append((m, tri, neq, media(ad), max(ad), sd, par, cpu, par / parede125, cpu / cpu125))
    out += [r"\bottomrule", r"\end{tabular}"]
    with open(os.path.join(AQUI, "tabelas", "ref_comparacao.tex"), "w") as f:
        f.write("% gerado por tabelas_comparacao.py a partir de dados/ref_comparacao.csv\n")
        f.write("\n".join(out) + "\n")

    # por caso
    out = [r"\begin{tabular}{@{}lccccc@{}}", r"\toprule",
           "caso & " + " & ".join(NOMES[m] for m in METODOS) + r" \\", r"\midrule"]
    for c in CASOS:
        cel = []
        for m in METODOS:
            g = br(float(d[(m, c)]["fator"]), 4)
            cel.append(g if m == "foco125" else "%s (%s)" % (g, pct(dif[(m, c)], sinal=True)))
        out.append("%s & %s \\\\" % ("det" if c == "det" else "r" + c, " & ".join(cel)))
    out += [r"\bottomrule", r"\end{tabular}"]
    with open(os.path.join(AQUI, "tabelas", "ref_comparacao_casos.tex"), "w") as f:
        f.write("% gerado por tabelas_comparacao.py a partir de dados/ref_comparacao.csv\n")
        f.write("\n".join(out) + "\n")

    # as etapas da refinada (para o texto)
    ref_l = [d[("refinada", c)] for c in CASOS]
    eg = media([float(r["seg_grossa"]) for r in ref_l])
    em = media([float(r["seg_malha"]) for r in ref_l])
    ef = media([float(r["seg_fina"]) for r in ref_l])
    par = media([float(r["parede_s"]) for r in ref_l])
    print("metodo tri neq dif_media dif_max dif_sinal parede cpu parede/125 cpu/125")
    for r in resumo:
        print("%-9s %7.1f %8.1f %7.4f %7.4f %+8.4f %7.1f %7.1f %5.3f %5.3f" % r)
    print("refinada: etapas medias grossa %.2f s, gmsh %.2f s, refinada %.2f s; parede %.2f s (%.1f%%, %.1f%%, %.1f%%)"
          % (eg, em, ef, par, 100 * eg / par, 100 * em / par, 100 * ef / par))
    print("refinada: triangulos de %d a %d; celulas de banda em media %.0f; fator grossa = foco375: %s" % (
        min(int(r["triangulos"]) for r in ref_l), max(int(r["triangulos"]) for r in ref_l),
        media([float(r["celulas_banda"]) for r in ref_l]),
        all(r["fator_grossa"] == d[("foco375", c)]["fator"] for r, c in zip(ref_l, CASOS))))


if __name__ == "__main__":
    main()
