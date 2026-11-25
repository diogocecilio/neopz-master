#!/bin/bash
# gera.sh: refaz o relatorio det_3solos (os 3 solos do Rahardjo, 3 chuvas de 24 h cada, variantes iii, ii e i, FS de
# hora em hora ate 48 h; rodadas do PC em ~/projects/runs/det_pc_3solos):
#   1. figuras, tabelas e numeros (gera_dados.py);
#   2. o LaTeX em duas passadas. A capa leva a data e o commit (+ se ha mudancas nesta pasta).
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
export SOURCE_DATE_EPOCH=1767225600
VARIANTES="${VARIANTES:-iii,ii,i}" python3 -I gera_dados.py "${1:-$HOME/projects/runs/det_pc_3solos}"
v="$(date +%d/%m/%Y), commit $(git rev-parse --short HEAD)$(git diff --quiet HEAD -- . ':!relatorio_det_3solos.pdf' 2>/dev/null || echo +)"
for i in 1 2; do pdflatex -interaction=nonstopmode "\\def\\versao{$v}\\input{relatorio_det_3solos.tex}" > /dev/null || true; done
grep -E "^!|undefined|Overfull" relatorio_det_3solos.log | head -20 || true
rm -f relatorio_det_3solos.{aux,log,out,toc}
echo "pronto: relatorio_det_3solos.pdf ($(pdfinfo relatorio_det_3solos.pdf | awk '/Pages/ {print $2}') p., versao $v)"
