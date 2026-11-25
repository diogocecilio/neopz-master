#!/bin/bash
# submeter_det_chuva.sh [nos = 4]: submete <nos> jobs det_chuva.sbatch (um no inteiro cada), que dividem os casos de
# casos_det_chuva.txt (a validacao deterministica da chuva). Volta em segundos: pode desconectar. Pode rodar de novo a
# qualquer momento: os jobs novos pegam o que estiver livre. Sem jobs detChuva na fila, libera antes as travas dos casos
# que nao terminaram (jobs que morreram). Ambiente: SAIDA ($DATA/det_chuva), TEMPO (2-00:00:00), PARTICAO, PARALELO
# (24), REFAZER_ERROS=1 (apaga tambem os casos que terminaram com codigo != 0, para rodarem de novo)
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
NOS="${1:-4}"
case "$NOS" in '' | *[!0-9]* | 0) echo "uso: $0 [numero de nos, padrao 4]" >&2; exit 2 ;; esac
export SAIDA="${SAIDA:-${DATA:-/mnt/data/${USER:-$(id -un)}}/det_chuva}"
EXE=bin_det/SlopeRainMonteCarlo
[ -x "$EXE" ] || { echo "compile antes, no no de login: ./compilar_det_chuva.sh" >&2; exit 1; }
grep -q "k_theta" "$EXE" || { echo "$EXE sem as opcoes k_theta=/diag=: recompile (./compilar_det_chuva.sh)" >&2; exit 1; }
[ -e "$SAIDA/PARAR" ] && { echo "existe $SAIDA/PARAR: apague antes (rm $SAIDA/PARAR)" >&2; exit 1; }
mkdir -p logs "$SAIDA"
fila="$(squeue -h -u "${USER:-$(id -un)}" -n detChuva -o %i)" || { echo "o squeue nao respondeu: tente de novo" >&2; exit 1; }
if [ -z "$fila" ]; then
    n=0
    for d in "$SAIDA"/*/; do
        d="${d%/}"
        [ -f "$d/.trava" ] || continue
        if [ ! -f "$d/fim" ]; then
            rm -f "$d/.trava"; n=$((n + 1))
        elif [ -n "${REFAZER_ERROS:-}" ] && [ "$(cat "$d/codigo")" != 0 ]; then
            rm -f "$d/.trava" "$d/fim" "$d/codigo"; n=$((n + 1))
        fi
    done
    echo "sem jobs detChuva na fila: $n travas liberadas"
fi
for ((j = 1; j <= NOS; j++)); do
    # shellcheck disable=SC2086
    sbatch --export=ALL ${TEMPO:+--time=$TEMPO} ${PARTICAO:+--partition=$PARTICAO} det_chuva.sbatch
done
echo "saida: $SAIDA  (acompanhar: ./det_chuva.sbatch resumo; tail $SAIDA/eventos.log)"
