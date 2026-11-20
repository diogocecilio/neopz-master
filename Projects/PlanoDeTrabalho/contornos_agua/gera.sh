#!/bin/bash
# gera.sh: refaz as figuras e as tabelas (figuras_contornos.py, que so le o Projects/CeronCinematico e as rodadas em
# ~/projects/runs) e compila relatorio_contornos.tex em duas passadas. A capa leva a data e o commit (+ se ha mudancas).
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
python3 -I figuras_contornos.py
v="$(date +%d/%m/%Y), commit $(git rev-parse --short HEAD)$(git diff --quiet -- . 2>/dev/null || echo +)"
for i in 1 2; do pdflatex -interaction=nonstopmode "\\def\\versao{$v}\\input{relatorio_contornos.tex}" > /dev/null || true; done
grep -E "^!|undefined|Overfull" relatorio_contornos.log | head -20 || true
echo "pronto: relatorio_contornos.pdf ($(pdfinfo relatorio_contornos.pdf | awk '/Pages/ {print $2}') p., versao $v)"
