#!/bin/bash
# submeter.sh [nos = 4]: submete <nos> jobs rah.sbatch (um no inteiro cada, PARALELO processos), que dividem os casos de
# casos.txt. Volta em segundos: pode desconectar. Pode rodar de novo a qualquer momento: os jobs novos pegam o que
# estiver livre. Sem jobs rah2007 na fila, libera antes as travas dos casos que nao terminaram (jobs que morreram).
# Ambiente: SAIDA ($DATA/rahardjo2007), TEMPO (1-00:00:00), PARTICAO, PARALELO (7), THREADS (3), REFAZER_ERROS=1 (refaz os casos que
# terminaram com codigo != 0)
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
NOS="${1:-4}"
case "$NOS" in '' | *[!0-9]* | 0) echo "uso: $0 [numero de nos, padrao 4]" >&2; exit 2 ;; esac
export SAIDA="${SAIDA:-${DATA:-/mnt/data/${USER:-$(id -un)}}/rahardjo2007}"
[ -x bin/Rahardjo2007 ] || { echo "compile antes, no no de login: ./compilar.sh" >&2; exit 1; }
[ -e "$SAIDA/PARAR" ] && { echo "existe $SAIDA/PARAR: apague antes (rm $SAIDA/PARAR)" >&2; exit 1; }
mkdir -p logs "$SAIDA"
fila="$(squeue -h -u "${USER:-$(id -un)}" -n rah2007 -o %i)" || { echo "o squeue nao respondeu: tente de novo" >&2; exit 1; }
if [ -z "$fila" ]; then
    n=0
    for d in "$SAIDA"/*/; do
        [ -d "$d" ] || continue
        d="${d%/}"
        [ -f "$d/.trava" ] || continue
        if [ ! -f "$d/fim" ]; then
            rm -f "$d/.trava"; n=$((n + 1))
        elif [ -n "${REFAZER_ERROS:-}" ] && [ "$(cat "$d/codigo")" != 0 ]; then
            rm -f "$d/.trava" "$d/fim" "$d/codigo"; n=$((n + 1))
        fi
    done
    echo "sem jobs rah2007 na fila: $n travas liberadas"
fi
for ((j = 1; j <= NOS; j++)); do
    # shellcheck disable=SC2086
    sbatch --export=ALL ${TEMPO:+--time=$TEMPO} ${PARTICAO:+--partition=$PARTICAO} rah.sbatch
done
echo "saida: $SAIDA  (acompanhar: ./rah.sbatch resumo; tail $SAIDA/eventos.log)"
