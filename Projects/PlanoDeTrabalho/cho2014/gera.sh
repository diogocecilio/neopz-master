#!/bin/bash
# gera.sh: refaz o relatorio da coluna 1D de Cho (2014) sob chuva:
#   1. o resumo da campanha cho2014 (resumo_cho2014.py, com a referencia em Python) sobre o que ja foi baixado (para
#      baixar de novo do fermi, rode antes o baixar_cho2014.sh no seu terminal); a saida vai para dados/ (a pasta da
#      campanha nao e tocada);
#   2. as tabelas, os numeros e as figuras (gera_dados.py): so le dados que ja existem (a campanha, o MC em Python, os
#      estudos da coluna e o line sampling); nenhuma simulacao nova. As figuras de Cho sao recortadas do PDF;
#   3. o LaTeX em duas passadas. A capa leva a data e o commit (+ se ha mudancas nesta pasta).
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
export SOURCE_DATE_EPOCH=1767225600
CAMPANHA="${CAMPANHA:-$HOME/projects/runs/fermi_cho2014}"
R=../../SlopeRandomFieldRefinement/cluster/resumo_cho2014.py
mkdir -p dados
python3 -I "$R" "$CAMPANHA" --python "$HOME/projects/runs/cho2014_python" --json dados/resumo.json > dados/resultado.txt
python3 -I gera_dados.py "$CAMPANHA"
v="$(date +%d/%m/%Y), commit $(git rev-parse --short HEAD)$(git diff --quiet HEAD -- . ':!relatorio_cho2014.pdf' 2>/dev/null || echo +)"
for i in 1 2; do pdflatex -interaction=nonstopmode "\\def\\versao{$v}\\input{relatorio_cho2014.tex}" > /dev/null || true; done
grep -E "^!|undefined|Overfull" relatorio_cho2014.log | head -20 || true
rm -f relatorio_cho2014.{aux,log,out,toc}
echo "pronto: relatorio_cho2014.pdf ($(pdfinfo relatorio_cho2014.pdf | awk '/Pages/ {print $2}') p., versao $v)"
