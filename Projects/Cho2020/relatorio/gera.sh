#!/bin/bash
# gera.sh: refaz o relatorio do Cho2020 (relatorio_cho2020.tex): figuras e tabelas (gera_dados.py, com os dados que ja
# existem em ~/projects/runs/cho2020) e o LaTeX em duas passadas
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
python3 -I gera_dados.py > /dev/null
v="$(date +%d/%m/%Y)"
for i in 1 2; do pdflatex -interaction=nonstopmode "\\def\\versao{$v}\\input{relatorio_cho2020.tex}" > /dev/null || true; done
grep -E "^!|undefined" relatorio_cho2020.log | head -20 || true
rm -f relatorio_cho2020.{aux,log,out,toc}
echo "pronto: relatorio_cho2020.pdf ($(pdfinfo relatorio_cho2020.pdf | awk '/Pages/ {print $2}') p.)"
