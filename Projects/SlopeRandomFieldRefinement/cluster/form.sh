#!/bin/bash
# form.sh: os FORM do Cho2010 (Projects/Cho2010/form_kl.py) com M = FORM_M termos por campo (ex.: "10 20"), que o
# primeiro job do submeter.sh roda antes da fila (job.sbatch): um depois do outro, cada um com PARALELO avaliacoes ao
# mesmo tempo, com as FORM_OPCOES da campanha (campanha.sh) e o executavel EXE. Cada um em $SAIDA/form_M<M>: form.log,
# form.json e ustar.txt (o ponto de projeto, um numero por linha: o desloc_xi= de um caso com espera=); se o FORM nao
# termina, ustar.txt.ERRO (os casos que esperam por ele ficam de fora). A trava form_M<M>/trava (noclobber) impede dois
# jobs no mesmo FORM; o submeter.sh a libera (e o ERRO) quando nao ha jobs da campanha na fila. Um FORM pronto e pulado.
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
source ./campanha.sh
: "${SAIDA:?}" "${EXE:?}" "${PARALELO:?}" "${FORM_M:?}"
FORM_PY="$(cd ../../Cho2010 && pwd)/form_kl.py"
for M in $FORM_M; do
    F="$SAIDA/form_M$M"
    [ -s "$F/ustar.txt" ] && continue
    mkdir -p "$F"
    (set -o noclobber; echo "$(hostname -s) job ${SLURM_JOB_ID:-local} $(date '+%F %T')" > "$F/trava") 2> /dev/null || {
        echo "$(date '+%F %T') FORM M = $M: outro job ja o roda ($F/trava)"
        continue
    }
    rm -f "$F/ustar.txt.ERRO"
    echo "$(date '+%F %T') FORM M = $M em $F: $PARALELO avaliacoes ao mesmo tempo, $FORM_OPCOES"
    # shellcheck disable=SC2086  # as opcoes se separam nos espacos
    if python3 -I "$FORM_PY" "$EXE" "$F" "M=$M" $FORM_OPCOES "paralelos=$PARALELO" > "$F/form.log" 2>&1 &&
        python3 -I -c 'import json, sys; print("\n".join(repr(x) for x in json.load(open(sys.argv[1]))["u"]))' \
            "$F/form.json" > "$F/ustar.tmp"; then
        mv -f "$F/ustar.tmp" "$F/ustar.txt"
        echo "$(date '+%F %T') FORM M = $M:$(grep -m 1 '  beta = ' "$F/form.log")"
    else
        echo "ERRO $(date '+%F %T'): veja $F/form.log" > "$F/ustar.txt.ERRO"
        echo "$(date '+%F %T') FORM M = $M: ERRO, veja $F/form.log (ultimas linhas abaixo)"
        tail -n 5 "$F/form.log"
    fi
done
