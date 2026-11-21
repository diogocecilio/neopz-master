#!/bin/bash
# submeter_ls_ceron.sh [nos = 4]: submete <nos> jobs ls_ceron_casos.sbatch (um no inteiro cada), que dividem os casos
# de casos_ls_ceron.txt (Ceron 2025, Secoes 6.1 a 6.5, FORM modal + line sampling + Monte Carlo pareado). Volta em
# segundos: pode desconectar. Pode rodar de novo a qualquer momento: os jobs novos pegam o que estiver livre.
# Sem jobs lsCeronC na fila, libera antes as travas de tarefas que nao terminaram (de jobs que morreram: fim do tempo,
# queda do no); as com ERRO ficam, a menos que RETENTAR=1.
# Ambiente: SAIDA ($DATA/ceron_ls_casos), TEMPO (2-00:00:00), PARTICAO, DEPENDE (id de job: --dependency=afterany),
# e as do sbatch (PARALELO, COV, NMIN, N, FORM_IT, NMC, BLOCO_MC, CASOS), que vao junto (--export=ALL).
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
NOS="${1:-4}"
case "$NOS" in '' | *[!0-9]* | 0) echo "uso: $0 [numero de nos, padrao 4]" >&2; exit 2 ;; esac
export SAIDA="${SAIDA:-${DATA:-/mnt/data/${USER:-$(id -un)}}/ceron_ls_casos}"
EXE=bin_ceron/SlopeDrawdownMonteCarlo
[ -x "$EXE" ] || { echo "compile antes, no no de login: ./compilar_ceron.sh" >&2; exit 1; }
# o executavel tem de ter o contorno=ceron (commit 4842ae2) e o resto da agua do artigo
v="$(cut -d_ -f1 bin_ceron/versao.txt)"
v="${v%+mod}"
if ! git merge-base --is-ancestor 4842ae2 "$v" 2> /dev/null; then
    echo "bin_ceron/ e de $v, anterior ao contorno=ceron (4842ae2): rode ./compilar_ceron.sh" >&2
    exit 1
fi
[ -e "$SAIDA/PARAR" ] && { echo "existe $SAIDA/PARAR: apague antes (rm $SAIDA/PARAR)" >&2; exit 1; }
mkdir -p logs "$SAIDA"

fila="$(squeue -h -u "${USER:-$(id -un)}" -n lsCeronC -o %i)" || { echo "o squeue nao respondeu: tente de novo" >&2; exit 1; }
if [ -z "$fila" ]; then
    n=0
    libera() { # <trava> <feito 0|1>
        [ -e "$1" ] && [ "$2" = 0 ] || return 0
        if grep -q '^ERRO' "$1" && [ -z "${RETENTAR:-}" ]; then
            echo "  com ERRO (fica; RETENTAR=1 refaz): $1"
            return 0
        fi
        rm -f "$1"
        n=$((n + 1))
    }
    for d in "$SAIDA"/*/; do
        d="${d%/}"
        [ -d "$d" ] || continue
        f=0; [ -s "$d/form/form.json" ] && f=1; libera "$d/.form.trava" $f
        f=0; [ -s "$d/det/det_0_-1.csv" ] && f=1; libera "$d/.det.trava" $f
        f=0; [ -e "$d/FEITO_LS" ] && f=1; libera "$d/.ls.trava" $f
        for t in "$d"/mc/.tarefas/*; do
            [ -e "$t" ] || continue
            ab="$(basename "$t")" a="${ab%_*}" b="${ab#*_}" c="$d/mc/mc_${ab}.csv" f=0
            [ -s "$c" ] && [ $(($(wc -l < "$c") - 1)) -ge $((b - a + 1)) ] && f=1
            libera "$t" $f
        done
    done
    echo "sem jobs lsCeronC na fila: $n travas de tarefas interrompidas liberadas"
fi

for ((j = 1; j <= NOS; j++)); do
    # shellcheck disable=SC2086
    sbatch --export=ALL ${TEMPO:+--time=$TEMPO} ${PARTICAO:+--partition=$PARTICAO} \
        ${DEPENDE:+--dependency=afterany:$DEPENDE} ls_ceron_casos.sbatch
done
echo "saida: $SAIDA  (acompanhar: ./status_ls_ceron.sh)"
