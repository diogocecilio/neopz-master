#!/bin/bash
# gera.sh: refaz o relatorio da proposta da curva de fragilidade P_f(I) da chuva (relatorio_fragilidade.tex):
#   1. figuras, tabelas e numeros (gera_dados.py), so com os dados que ja existem em ~/projects/runs;
#   2. o LaTeX em duas passadas. A capa leva a data e o commit (+ se ha mudancas nesta pasta).
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
export SOURCE_DATE_EPOCH=1767225600
python3 -I gera_dados.py > /dev/null
v="$(date +%d/%m/%Y), commit $(git rev-parse --short HEAD)$(git diff --quiet HEAD -- . ':!relatorio_fragilidade.pdf' 2>/dev/null || echo +)"
for i in 1 2; do pdflatex -interaction=nonstopmode "\\def\\versao{$v}\\input{relatorio_fragilidade.tex}" > /dev/null || true; done
grep -E "^!|undefined|Overfull" relatorio_fragilidade.log | head -20 || true
rm -f relatorio_fragilidade.{aux,log,out,toc}
echo "pronto: relatorio_fragilidade.pdf ($(pdfinfo relatorio_fragilidade.pdf | awk '/Pages/ {print $2}') p., versao $v)"
