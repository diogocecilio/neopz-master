#!/bin/bash
# gera.sh: refaz o relatorio do Kang2020 (relatorio_kang2020.tex): copia as figuras das rodadas (~/projects/runs/kang2020,
# passos 3 a 7), refaz as tabelas (o FS nosso contra o lido da Fig. 9) e compila o LaTeX em duas passadas
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
K=~/projects/runs/kang2020/passo4
cp ~/projects/runs/kang2020/malha_coluna.pdf figs/kang_malha.pdf
# a malha do artigo (a faixa do FLAC, Passo 6) e a comparacao faixa x coluna x artigo
P6=~/projects/runs/kang2020/passo6
cp "$P6/faixa_malha.pdf" figs/kang_faixa_malha.pdf
cp "$P6/faixa_coluna_kang.png" figs/kang_faixa_coluna.png
# Passo 7: a coluna normal ao talude contra a faixa (chuva por metro horizontal)
P7=~/projects/runs/kang2020/passo7
cp "$P7/passo7.png" figs/kang_normal.png
python3 - "$P7" "$K" ~/projects/neopz-master/Projects/Kang2020/data/fig9_lido.csv > tab/kang_normal.tex <<'PY'
import csv, sys
P7, P4, lidoArq = sys.argv[1:4]
lido = [r for r in csv.reader(open(lidoArq)) if r and not r[0].startswith("#")][1:]
L = {(int(z), m, int(t)): float(f) for z, m, t, f in lido}
def serie(a, c): return [float(r[c]) for r in csv.DictReader(open(a))]
v = lambda x: f"{x:.2f}".replace(".", "{,}")
print(r"\begin{tabular}{cc|cccc|c}")
print(r"\toprule")
print(r"zona & modelo & coluna normal & faixa & coluna vertical & Kang & $\max|\Delta|$ normal--faixa\\")
print(r"\midrule")
for z in range(1, 6):
    for m in ("sp", "tp"):
        n = serie(f"{P7}/normal_zona{z}_{m}_fs.csv", "FS_beta"); f = serie(f"{P7}/faixah_zona{z}_{m}_fs.csv", "FS_beta")
        c = serie(f"{P4}/kang_zona{z}_{m}_fs.csv", "FS_medio")
        d = max(abs(a - b) for a, b in zip(n, f))
        print(f"{z} & {m.upper()} & {v(n[-1])} & {v(f[-1])} & {v(c[-1])} & {v(L[(z, m, 17)])} & {v(d)}" + r" \\")
print(r"\bottomrule")
print(r"\end{tabular}")
PY
cp "$K/fig9_nosso.png" figs/kang_fig9_nosso.png
# os perfis com ar ativo: a coluna (cheia) e a faixa do FLAC na vertical do meio (tracejada)
python3 ~/projects/neopz-master/Projects/Kang2020/perfis_faixa_coluna.py "$K" ~/projects/runs/kang2020/passo6 tp \
    figs/kang_perfis_tp.png > /dev/null
# a tabela do FS no angulo medio (nosso / Kang lido da Fig. 9), TP e SP
python3 - "$K" ~/projects/neopz-master/Projects/Kang2020/data/fig9_lido.csv > tab/kang_fs.tex <<'PY'
import csv, sys
P, lidoArq = sys.argv[1], sys.argv[2]
lido = [r for r in csv.reader(open(lidoArq)) if r and not r[0].startswith("#")][1:]
L = {(int(z), m, int(t)): float(f) for z, m, t, f in lido}
N = {}
for z in range(1, 6):
    for m in ("sp", "tp"):
        for r in csv.DictReader(open(f"{P}/kang_zona{z}_{m}_fs.csv")):
            N[(z, m, int(r["t_h"]))] = float(r["FS_medio"])
v = lambda x: f"{x:.2f}".replace(".", "{,}")
print(r"\begin{tabular}{c|cccc|cccc}")
print(r"\toprule")
print(r" & \multicolumn{4}{c|}{TP (ar ativo)} & \multicolumn{4}{c}{SP (ar passivo)}\\")
print(r"zona & 3~h & 9~h & 12~h & 17~h & 3~h & 9~h & 12~h & 17~h\\")
print(r"\midrule")
for z in range(1, 6):
    cel = [f"{v(N[(z, m, t)])} / {v(L[(z, m, t)]) if (z, m, t) in L else '--'}" for m in ("tp", "sp") for t in (3, 9, 12, 17)]
    print(f"{z} & " + " & ".join(cel) + r" \\")
print(r"\bottomrule")
print(r"\end{tabular}")
PY
# o FS final (17 h) no angulo medio: faixa / coluna / artigo, SP e TP, e a conferencia analitica do FS saturado pela
# Eq. 14 (S_w = 1, W pela Eq. 15): p_w = gamma_w z cos^2(beta) (o escoamento paralelo ao talude, o da faixa drenando) e
# p_w = gamma_w z (o hidrostatico, o da coluna)
python3 - "$P6" "$K" ~/projects/neopz-master/Projects/Kang2020/data/fig9_lido.csv > tab/kang_faixa.tex <<'PY'
import csv, math, sys
P6, P4, lidoArq = sys.argv[1], sys.argv[2], sys.argv[3]
lido = [r for r in csv.reader(open(lidoArq)) if r and not r[0].startswith("#")][1:]
L = {(int(z), m, int(t)): float(f) for z, m, t, f in lido}
# as zonas (Tabela 1 do artigo; n = theta_s da Fig. 4a; o angulo medio da Fig. 9)
ZONAS = {1: (25.3, 9.6, 15.0, 0.50, 30.7), 2: (28.5, 5.8, 14.9, 0.50, 31.0), 3: (37.6, 7.7, 14.0, 0.51, 31.2),
         4: (30.9, 7.6, 14.5, 0.51, 31.8), 5: (28.2, 6.3, 15.0, 0.51, 31.1)}
gw, H = 9.81, 2.
def fs_sat(z_, fator):
    """o minimo da Eq. 14 na profundidade, saturado, com p_w = fator gamma_w z"""
    phi, c, gd, n, b = ZONAS[z_]
    b, phi = math.radians(b), math.radians(phi)
    m = float("inf")
    for i in range(1, 201):
        z = H * i / 200
        W = (gd + n * gw) * z
        m = min(m, ((W * math.cos(b) ** 2 - fator(b) * gw * z) * math.tan(phi) + c) / (W * math.sin(b) * math.cos(b)))
    return m
def fim(arq, col):
    R = list(csv.DictReader(open(arq)))
    return float(R[-1][col])
v = lambda x: f"{x:.2f}".replace(".", "{,}")
print(r"\begin{tabular}{c|ccc|ccc|cc}")
print(r"\toprule")
print(r" & \multicolumn{3}{c|}{SP (ar passivo)} & \multicolumn{3}{c|}{TP (ar ativo)} & "
      r"\multicolumn{2}{c}{Eq.~14 saturado, $S_w = 1$}\\")
print(r"zona & faixa & coluna & Kang & faixa & coluna & Kang & $p_w = \gamma_w z\cos^2\beta$ & $p_w = \gamma_w z$\\")
print(r"\midrule")
for z in range(1, 6):
    cel = []
    for m in ("sp", "tp"):
        cel += [v(fim(f"{P6}/faixa_zona{z}_{m}_fs.csv", "FS_beta")), v(fim(f"{P4}/kang_zona{z}_{m}_fs.csv", "FS_medio")),
                v(L[(z, m, 17)]) if (z, m, 17) in L else "--"]
    cel += [v(fs_sat(z, lambda b: math.cos(b) ** 2)), v(fs_sat(z, lambda b: 1.))]
    print(f"{z} & " + " & ".join(cel) + r" \\")
print(r"\bottomrule")
print(r"\end{tabular}")
PY

cp "$K/fig8_nosso.png" figs/kang_fig8_nosso.png
cp "$K/fig7_nosso.png" figs/kang_fig7_nosso.png
python3 ~/projects/neopz-master/Projects/Kang2020/perfis_faixa_coluna.py "$K" "$P6" sp figs/kang_perfis_sp.png > /dev/null
v="$(date +%d/%m/%Y)"
for i in 1 2; do pdflatex -interaction=nonstopmode "\\def\\versao{$v}\\input{relatorio_kang2020.tex}" > /dev/null || true; done
grep -E "^!|undefined" relatorio_kang2020.log | head -20 || true
rm -f relatorio_kang2020.{aux,log,out,toc}
echo "pronto: relatorio_kang2020.pdf ($(pdfinfo relatorio_kang2020.pdf | awk '/Pages/ {print $2}') p.)"
