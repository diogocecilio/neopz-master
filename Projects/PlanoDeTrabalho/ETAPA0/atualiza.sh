#!/bin/bash
# atualiza.sh: refaz o relatorio da etapa 0 com os dados de agora. Roda o figuras_etapa0.py (le as rodadas, refaz
# figs/fig_e0_*.pdf, numeros_etapa0.tex e as copias em resultados/; imprime cada numero com o arquivo de onde sai,
# na tela) e compila o relatorio_etapa0.tex duas vezes numa pasta temporaria, copiando so o PDF para ca.
# So le as pastas das rodadas. Roda com nice, para nao atrapalhar as simulacoes.
set -euo pipefail
cd "$(dirname "$0")"
nice -n 19 python3 -I figuras_etapa0.py "$@"
T=$(mktemp -d)
trap 'rm -rf "$T"' EXIT
cp relatorio_etapa0.tex numeros_etapa0.tex "$T"/
cp -r figs "$T"/
(cd "$T" && for i in 1 2; do nice -n 19 pdflatex -interaction=nonstopmode -halt-on-error relatorio_etapa0.tex \
    > /dev/null || { tail -40 relatorio_etapa0.log; exit 1; }; done)
grep -E "Overfull|Underfull|undefined|Warning" "$T/relatorio_etapa0.log" || true
cp "$T/relatorio_etapa0.pdf" .
echo "relatorio: $PWD/relatorio_etapa0.pdf ($(pdfinfo relatorio_etapa0.pdf 2>/dev/null | awk '/^Pages/ {print $2}') paginas)"
