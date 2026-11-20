#!/bin/bash
# gera_relatorio.sh: refaz as figuras e as tabelas com os dados que houver (gera_dados.py) e compila
# relatorio_comparacao.tex. A capa leva o commit (+ se ha mudancas).
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
python3 -I gera_dados.py
v="$(git rev-parse --short HEAD)$(git diff --quiet -- . || echo +)"
for i in 1 2; do pdflatex -interaction=nonstopmode "\\def\\versao{$v}\\input{relatorio_comparacao.tex}" > /dev/null || true; done
grep -E "^!|undefined" relatorio_comparacao.log | head || true
echo "pronto: relatorio_comparacao.pdf ($(pdfinfo relatorio_comparacao.pdf | awk '/Pages/ {print $2}') p., versao $v)"
