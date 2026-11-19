#!/bin/bash
# form.sh: os FORM do Cho2010 (Projects/Cho2010/form_kl.py) que os jobs do submeter.sh rodam antes da fila (job.sbatch),
# um depois do outro, cada um com PARALELO avaliacoes ao mesmo tempo e o executavel EXE:
#   FORM_M="10 20"   (o primeiro job) um FORM por M, em $SAIDA/form_M<M>, com as FORM_OPCOES da campanha (campanha.sh);
#   FORM_LISTA       (todos os jobs; campanha.sh) FORM com nome e opcoes proprias, "nome|opcoes;nome|opcoes;...", em
#                    $SAIDA/form_<nome>, com FORM_OPCOES antes das opcoes do nome. A trava de cada FORM divide a lista
#                    entre os jobs: ponha os mais caros primeiro. @SAIDA@ nas opcoes vira $SAIDA.
# Cada FORM grava form.log, form.json e ustar.txt (o ponto de projeto, um numero por linha: o desloc_xi= de um caso com
# espera=); se o FORM nao termina ou nao converge (convergiu falso no form.json), ustar.txt.ERRO (os casos que esperam
# por ele ficam de fora). A trava <pasta>/trava (noclobber) impede dois jobs no mesmo FORM; o submeter.sh a libera (e o
# ERRO) quando nao ha jobs da campanha na fila. Um FORM pronto e pulado.
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
source ./campanha.sh
: "${SAIDA:?}" "${EXE:?}" "${PARALELO:?}"
FORM_PY="$(cd ../../Cho2010 && pwd)/form_kl.py"
# a lista de FORM desta chamada: "pasta|rotulo|opcoes"
lista=()
for M in ${FORM_M:-}; do
    lista+=("form_M$M|M = $M|M=$M $FORM_OPCOES")
done
IFS=';' read -ra itens <<< "${FORM_LISTA:-}"
for it in "${itens[@]}"; do
    [ -n "${it// /}" ] || continue
    nome="${it%%|*}" opc="${it#*|}"
    lista+=("form_$nome|$nome|$FORM_OPCOES $opc")
done
for item in "${lista[@]}"; do
    IFS='|' read -r pasta rot opc <<< "$item"
    opc="${opc//@SAIDA@/$SAIDA}"
    F="$SAIDA/$pasta"
    [ -s "$F/ustar.txt" ] && continue
    mkdir -p "$F"
    (set -o noclobber; echo "$(hostname -s) job ${SLURM_JOB_ID:-local} $(date '+%F %T')" > "$F/trava") 2> /dev/null || {
        echo "$(date '+%F %T') FORM $rot: outro job ja o roda ($F/trava)"
        continue
    }
    rm -f "$F/ustar.txt.ERRO"
    echo "$(date '+%F %T') FORM $rot em $F: $PARALELO avaliacoes ao mesmo tempo, $opc"
    # shellcheck disable=SC2086  # as opcoes se separam nos espacos
    if python3 -I "$FORM_PY" "$EXE" "$F" $opc "paralelos=$PARALELO" > "$F/form.log" 2>&1 &&
        python3 -I -c 'import json, sys; d = json.load(open(sys.argv[1])); assert d.get("convergiu"), "nao convergiu"
print("\n".join(repr(x) for x in d["u"]))' \
            "$F/form.json" > "$F/ustar.tmp"; then
        mv -f "$F/ustar.tmp" "$F/ustar.txt"
        echo "$(date '+%F %T') FORM $rot:$(grep -m 1 '  beta = ' "$F/form.log")"
    else
        rm -f "$F/ustar.tmp"
        echo "ERRO $(date '+%F %T'): veja $F/form.log" > "$F/ustar.txt.ERRO"
        echo "$(date '+%F %T') FORM $rot: ERRO, veja $F/form.log (ultimas linhas abaixo)"
        tail -n 5 "$F/form.log"
    fi
done
