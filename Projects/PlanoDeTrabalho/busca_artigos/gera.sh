#!/bin/bash
# gera.sh: compila relatorio_busca.tex (pdflatex duas vezes, para o sumario e as citacoes).
# A capa leva o commit curto (+ se o .tex ou este script diferem do commit ou ainda nao foram commitados).
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
v="$(git rev-parse --short HEAD)$([ -z "$(git status --porcelain -- relatorio_busca.tex gera.sh)" ] || echo +)"
for i in 1 2; do nice -n 10 pdflatex -interaction=nonstopmode "\\def\\versao{$v}\\input{relatorio_busca.tex}" > /dev/null || true; done
grep -E "^!|undefined" relatorio_busca.log | head || true
echo "pronto: relatorio_busca.pdf ($(pdfinfo relatorio_busca.pdf | awk '/Pages/ {print $2}') p., versao $v)"
