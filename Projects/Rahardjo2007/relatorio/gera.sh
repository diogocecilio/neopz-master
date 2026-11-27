#!/bin/bash
# gera.sh: refaz o relatorio do Rahardjo2007 (relatorio_rahardjo2007.tex): figuras e tabelas (gera_dados.py, com os dados
# que ja existem em ~/projects/runs/rahardjo2007) e o LaTeX em duas passadas. A capa leva a data e o commit
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
python3 -I gera_dados.py
v="$(date +%d/%m/%Y), commit $(git rev-parse --short HEAD)$(git diff --quiet HEAD -- .. ':!relatorio_rahardjo2007.pdf' 2>/dev/null || echo +)"
for i in 1 2; do pdflatex -interaction=nonstopmode "\\def\\versao{$v}\\input{relatorio_rahardjo2007.tex}" > /dev/null || true; done
grep -E "^!|undefined" relatorio_rahardjo2007.log | head -20 || true
rm -f relatorio_rahardjo2007.{aux,log,out,toc}
echo "pronto: relatorio_rahardjo2007.pdf ($(pdfinfo relatorio_rahardjo2007.pdf | awk '/Pages/ {print $2}') p.)"
