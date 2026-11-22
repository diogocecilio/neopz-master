#!/bin/bash
# submeter_hu2024.sh [nos = 4]: submete <nos> jobs hu2024_casos.sbatch (um no inteiro cada), que dividem os casos de
# casos_hu2024.txt (caso A, Hu et al. 2024: transiente, deterministico, FORM modal, line sampling e Monte Carlo direto). Volta
# em segundos: pode desconectar. Pode rodar de novo a qualquer momento: os jobs novos pegam o que estiver livre.
# Sem jobs huA na fila, libera antes as travas de tarefas que nao terminaram (de jobs que morreram: fim do tempo, queda do
# no); as com ERRO ficam, a menos que RETENTAR=1.
# Ambiente: KS (2e-5), SAIDA ($DATA/hu2024_ks$KS), TEMPO (2-00:00:00), PARTICAO, DEPENDE (id de job: --dependency=afterany),
# e as do sbatch (M, BASE_NYS, TRANS, TEMPOS, PARALELO, COV, NMIN, N, FORM_IT, FORM_THREADS, NMC_MAX, BLOCO_MC, NCHECK,
# CASOS), que vao junto (--export=ALL). Duas campanhas (os dois k_s) dividem os 4 jobs que rodam: KS=5.2e-6 com DEPENDE=<o
# ultimo job da primeira> espera ela.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
NOS="${1:-4}"
case "$NOS" in '' | *[!0-9]* | 0) echo "uso: $0 [numero de nos, padrao 4]" >&2; exit 2 ;; esac
export KS="${KS:-2e-5}"
export SAIDA="${SAIDA:-${DATA:-/mnt/data/${USER:-$(id -un)}}/hu2024_ks$KS}"
EXE=bin_hu/SlopeRainMonteCarlo
[ -x "$EXE" ] || { echo "compile antes, no no de login: ./compilar_hu2024.sh" >&2; exit 1; }
# o executavel tem de ter as opcoes do caso A (geometria=, vg=, pressoes=, nystrom=, gradiente=; commit COMMIT_HU)
COMMIT_HU=2bc901a
v="$(cut -d_ -f1 bin_hu/versao.txt)"
v="${v%+mod}"
if ! git merge-base --is-ancestor "$COMMIT_HU" "$v" 2> /dev/null; then
    echo "bin_hu/ e de $v, anterior as opcoes do caso A ($COMMIT_HU): rode ./compilar_hu2024.sh" >&2
    exit 1
fi
[ -e "$SAIDA/PARAR" ] && { echo "existe $SAIDA/PARAR: apague antes (rm $SAIDA/PARAR)" >&2; exit 1; }
[ -x "${PYV:-$HOME/venv_la/bin/python}" ] || { echo "falta o python com numpy e scipy (${PYV:-$HOME/venv_la/bin/python}): a base de Nystrom e feita com ele" >&2; exit 1; }
mkdir -p logs "$SAIDA"

fila="$(squeue -h -u "${USER:-$(id -un)}" -n huA -o %i)" || { echo "o squeue nao respondeu: tente de novo" >&2; exit 1; }
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
    BASE_NYS="${BASE_NYS:-${DATA:-/mnt/data/${USER:-$(id -un)}}/hu2024/base_hu_M${M:-150}.txt}"
    f=0; [ -s "$BASE_NYS" ] && f=1; libera "$BASE_NYS.trava" $f
    for d in "$SAIDA"/*/; do
        d="${d%/}"
        [ -d "$d" ] || continue
        f=0; [ -s "$d/pressoes.txt" ] && f=1; libera "$d/.tr.trava" $f
        f=0; [ -s "$d/form/form.json" ] && f=1; libera "$d/.form.trava" $f
        f=0; [ -s "$d/det/det_0_-1.csv" ] && [ "$(wc -l < "$d/det/det_0_-1.csv")" -gt 13 ] && f=1; libera "$d/.det.trava" $f
        f=0; [ -e "$d/FEITO_LS" ] && f=1; libera "$d/.ls.trava" $f
        for p in mc mc13; do
            for t in "$d"/$p/.tarefas/*; do
                [ -e "$t" ] || continue
                # o bloco esta feito se o csv tem as realizacoes todas (qualquer numero de instantes: o log diz o fim)
                ab="$(basename "$t")" c="$d/$p/mc_${ab}.csv" g="$d/$p/log_${ab}.txt" f=0
                [ -s "$c" ] && [ -s "$g" ] && grep -q 'arquivos:' "$g" && f=1
                libera "$t" $f
            done
        done
    done
    echo "sem jobs huA na fila: $n travas de tarefas interrompidas liberadas"
fi

for ((j = 1; j <= NOS; j++)); do
    # shellcheck disable=SC2086
    sbatch --export=ALL ${TEMPO:+--time=$TEMPO} ${PARTICAO:+--partition=$PARTICAO} \
        ${DEPENDE:+--dependency=afterany:$DEPENDE} hu2024_casos.sbatch
done
echo "saida: $SAIDA  (acompanhar: KS=$KS ./status_hu2024.sh)"
