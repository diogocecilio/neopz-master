#!/bin/bash
# gera.sh: refaz o relatorio das Secoes 6.1 a 6.5 do Ceron com a agua do artigo (campanha lsCeronC):
#   1. o resumo e as figuras do artigo sobre o que ja foi baixado (SEM_BAIXAR=1 baixar_ls_ceron.sh; para baixar de
#      novo do fermi, rode antes o baixar_ls_ceron.sh no seu terminal, que pede a senha);
#   2. as tabelas, os numeros e as figuras novas (gera_dados.py);
#   3. o LaTeX em duas passadas. A capa leva a data e o commit (+ se ha mudancas nesta pasta).
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
CAMPANHA="${CAMPANHA:-$HOME/projects/runs/fermi_ceron_ls_casos}"
SEM_BAIXAR=1 DEST="$CAMPANHA" ../../SlopeRandomFieldRefinement/cluster/baixar_ls_ceron.sh > /dev/null
python3 -I gera_dados.py "$CAMPANHA"
v="$(date +%d/%m/%Y), commit $(git rev-parse --short HEAD)$(git diff --quiet -- . 2>/dev/null || echo +)"
for i in 1 2; do pdflatex -interaction=nonstopmode "\\def\\versao{$v}\\input{relatorio_ceron_secao6.tex}" > /dev/null || true; done
grep -E "^!|undefined|Overfull" relatorio_ceron_secao6.log | head -20 || true
rm -f relatorio_ceron_secao6.{aux,log,out,toc}
echo "pronto: relatorio_ceron_secao6.pdf ($(pdfinfo relatorio_ceron_secao6.pdf | awk '/Pages/ {print $2}') p., versao $v)"
