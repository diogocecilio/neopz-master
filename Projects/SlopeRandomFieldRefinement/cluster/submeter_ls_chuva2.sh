#!/bin/bash
# submeter_ls_chuva2.sh [nos = 1]: submete <nos> jobs ls_chuva2.sbatch (um no inteiro cada), que dividem os LS de
# casos_chuva2_ls.txt (FORM modal + line sampling da campanha chuva2). Cada LS roda inteiro num no so (24 linhas em voo):
# um no por LS ao mesmo tempo basta (com so a referencia, 1). Volta em segundos: pode desconectar. Pode rodar de novo a
# qualquer momento: os jobs novos pegam o que estiver livre.
# Sem jobs lsChuva2 na fila, libera antes as travas de tarefas que nao terminaram (de jobs que morreram: fim do tempo,
# queda do no); as com ERRO ficam, a menos que RETENTAR=1.
# Ambiente: SAIDA ($DATA/chuva2_ls), TEMPO (2-00:00:00), PARTICAO, DEPENDE (id de job: --dependency=afterany), REUSO (0;
# 1 confere o patch do campo de k no executavel) e as do sbatch (SAIDA_MC, IREF, PARALELO, COV, NMIN, N, TOL, FORM_IT,
# FORM_THREADS, FORM_EXTRA, LS_CASOS, EXE), que vao junto (--export=ALL). O fermi roda 4 jobs por usuario e aceita 8 na
# fila: com o MC da chuva2 (jobs chuva2) nos 4 nos, use DEPENDE ou submeta depois.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
NOS="${1:-1}"
case "$NOS" in '' | *[!0-9]* | 0) echo "uso: $0 [numero de nos, padrao 1]" >&2; exit 2 ;; esac
export SAIDA="${SAIDA:-${DATA:-/mnt/data/${USER:-$(id -un)}}/chuva2_ls}"
export REUSO="${REUSO:-0}"
EXE="${EXE:-bin_chuva2/SlopeRainMonteCarlo}"
[ -x "$EXE" ] || { echo "compile antes, no no de login: ./compilar_ls_chuva2.sh (ou EXE=<outro>)" >&2; exit 1; }
export EXE="$(readlink -f "$EXE")"
# o executavel tem de ter xi=, gradiente= e tempos= (montecarlo.cpp desde 2bc901a); com REUSO=1, o patch do campo de k
v="$(cut -d_ -f1 "$(dirname "$EXE")/versao.txt" 2> /dev/null || echo sem_versao)"
v="${v%+mod}"
if ! git merge-base --is-ancestor 2bc901a "$v" 2> /dev/null; then
    echo "$EXE e de $v, anterior as opcoes xi=/gradiente= do montecarlo.cpp (2bc901a): recompile" >&2
    exit 1
fi
if [ "$REUSO" = 1 ] && ! grep -q "campo_k=" "$EXE"; then
    echo "REUSO=1: $EXE nao tem o patch do campo de k (grava_pressoes= com k aleatorio): rode sem REUSO" >&2
    exit 1
fi
[ -e "$SAIDA/PARAR" ] && { echo "existe $SAIDA/PARAR: apague antes (rm $SAIDA/PARAR)" >&2; exit 1; }
mkdir -p logs "$SAIDA"

fila="$(squeue -h -u "${USER:-$(id -un)}" -n lsChuva2 -o %i)" || { echo "o squeue nao respondeu: tente de novo" >&2; exit 1; }
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
        f=0; [ -s "$d/det/det_0_-1.csv" ] && [ "$(wc -l < "$d/det/det_0_-1.csv")" -gt 1 ] && f=1; libera "$d/.det.trava" $f
        f=0; [ -e "$d/FEITO_LS" ] && f=1; libera "$d/.ls.trava" $f
    done
    echo "sem jobs lsChuva2 na fila: $n travas de tarefas interrompidas liberadas"
fi

for ((j = 1; j <= NOS; j++)); do
    # shellcheck disable=SC2086
    sbatch --export=ALL ${TEMPO:+--time=$TEMPO} ${PARTICAO:+--partition=$PARTICAO} \
        ${DEPENDE:+--dependency=afterany:$DEPENDE} ls_chuva2.sbatch
done
echo "saida: $SAIDA  (acompanhar: ./status_ls_chuva2.sh)"
