#!/bin/bash
# atualiza.sh: regenera o relatorio de situacao com os dados baixados (gera_situacao.py, que tambem refaz os numeros e as
# figuras das etapas 1 a 3), o indice do README (../documentos.py) e compila tres PDFs: a formulacao
# (docs/confiabilidade_campos/formulacao, antes, porque a situacao cita as secoes dela pelo .aux), a situacao e o relatorio
# das etapas 1 a 3. Na capa da formulacao e da situacao: a data e o commit, com + se as fontes delas mudaram desde ele.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
F=../../../docs/confiabilidade_campos/formulacao
mudou="$(git status --porcelain -- . "$F" ../documentos.py | grep -v -E '(situacao|formulacao)\.(tex|pdf)$|/figs/' || true)"
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
compila . situacao "$V"
# o das etapas so quando o texto dele mudou (os numeros vem dos mesmos dados): evita um PDF novo igual ao anterior
if ! git diff --quiet -- ../ETAPAS1a3/relatorio_etapas.tex || [ ! -f ../ETAPAS1a3/relatorio_etapas.pdf ]; then
    compila ../ETAPAS1a3 relatorio_etapas
fi
python3 -I ../documentos.py # de novo: as paginas dos PDFs novos
