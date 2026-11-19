#!/bin/bash
# gera_relatorio.sh: compila relatorio_viga.tex (a viga com campo aleatorio em Mathematica: MC, FORM, SORM e ISMC). As
# figuras em figs/ saem do proprio codigo (codigo/viga_FORM5_AI200.txt grava os PDFs em figuras_AI_M5_M200/ na pasta do
# notebook; os PNG do estudo em M saem de codigo/viga_estudo_M_ISMC.txt). A capa leva o commit (+ se ha mudancas).
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
v="$(git rev-parse --short HEAD)$(git diff --quiet -- . || echo +)"
for i in 1 2; do pdflatex -interaction=nonstopmode "\\def\\versao{$v}\\input{relatorio_viga.tex}" > /dev/null || true; done
grep -E "^!|undefined" relatorio_viga.log | head || true
echo "pronto: relatorio_viga.pdf ($(pdfinfo relatorio_viga.pdf | awk '/Pages/ {print $2}') p., versao $v)"
