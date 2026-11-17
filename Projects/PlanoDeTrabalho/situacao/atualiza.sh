#!/bin/bash
# atualiza.sh: regenera os relatorios com os dados baixados (gera_situacao.py roda o ../ETAPAS1a3/relatorio_etapas.py, que
# refaz os numeros, as figuras e o relatorio das etapas com todos os resultados), o indice do README (../documentos.py) e
# compila tres PDFs, nesta ordem, porque cada um cita os anteriores pelo .aux: a formulacao
# (docs/confiabilidade_campos/formulacao), o relatorio das etapas e a situacao. Na capa dos tres: a data e o commit, com +
# se as fontes mudaram desde ele (fontes: o formulacao.tex, os dois geradores, o documentos.py; os .tex dos relatorios, os
# PDFs e as figuras sao saidas).
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
F=../../../docs/confiabilidade_campos/formulacao
mudou="$(git status --porcelain -- . "$F" ../documentos.py ../ETAPAS1a3 |
    grep -v -E 'situacao\.(tex|pdf)$|relatorio_etapas\.(tex|pdf)$|formulacao\.pdf$|/figs/' || true)"
V="$(date +%d/%m/%Y), commit $(git log -1 --format=%h)${mudou:++}"

python3 -I gera_situacao.py
python3 -I ../documentos.py

# compila <pasta> <nome> [versao]: duas passadas (tres na primeira vez, para o sumario); mostra erros e referencias soltas
compila() {
    local d="$1" n="$2" v="${3:-}" cmd i
    cmd="\\input{$n.tex}"
    [ -n "$v" ] && cmd="\\def\\versao{$v}$cmd"
    for i in 1 2 3; do
        (cd "$d" && pdflatex -interaction=nonstopmode -halt-on-error -jobname="$n" "$cmd" > /dev/null) || {
            tail -25 "$d/$n.log"
            echo "ERRO ao compilar $d/$n.tex" >&2
            exit 1
        }
        grep -q "Rerun to get" "$d/$n.log" || [ "$i" -lt 2 ] || break
    done
    grep -E "undefined|multiply defined" "$d/$n.log" | sort -u || true
    echo "$d/$n.pdf: $(pdfinfo "$d/$n.pdf" | awk '/^Pages/ {print $2}') paginas"
}
compila "$F" formulacao "$V"
compila ../ETAPAS1a3 relatorio_etapas "$V"
compila . situacao "$V"
python3 -I ../documentos.py # de novo: as paginas dos PDFs novos
