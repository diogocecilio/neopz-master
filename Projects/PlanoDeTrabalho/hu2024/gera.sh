#!/bin/bash
# gera.sh: refaz o relatorio do caso A (Hu et al. 2024) da chuva:
#   1. o resumo da campanha hu2024 (resumo_hu2024.py) sobre o que ja foi baixado (para baixar de novo do fermi, rode
#      antes o baixar_hu2024.sh no seu terminal, que pede a senha);
#   2. as tabelas, os numeros e as figuras (gera_dados.py), com as conferencias do PC em ~/projects/runs/hu2024 e as
#      curvas de Hu lidas dos vetores do PDF;
#   3. o LaTeX em duas passadas. A capa leva a data e o commit (+ se ha mudancas nesta pasta).
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
# datas fixas nos PDFs do matplotlib: as figuras so mudam quando o conteudo muda (e a capa nao ganha "+" a toa)
export SOURCE_DATE_EPOCH=1767225600
CAMPANHA="${CAMPANHA:-$HOME/projects/runs/fermi_hu2024_ks2e-5}"
R=../../SlopeRandomFieldRefinement/cluster/resumo_hu2024.py
python3 -I "$R" "$CAMPANHA" --json "$CAMPANHA/resumo.json" --csv "$CAMPANHA/resultado.csv" > "$CAMPANHA/resultado.txt"
python3 -I gera_dados.py "$CAMPANHA"
v="$(date +%d/%m/%Y), commit $(git rev-parse --short HEAD)$(git diff --quiet HEAD -- . ':!relatorio_hu2024.pdf' 2>/dev/null || echo +)"
for i in 1 2; do pdflatex -interaction=nonstopmode "\\def\\versao{$v}\\input{relatorio_hu2024.tex}" > /dev/null || true; done
grep -E "^!|undefined|Overfull" relatorio_hu2024.log | head -20 || true
rm -f relatorio_hu2024.{aux,log,out,toc}
echo "pronto: relatorio_hu2024.pdf ($(pdfinfo relatorio_hu2024.pdf | awk '/Pages/ {print $2}') p., versao $v)"
