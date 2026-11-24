#!/bin/bash
# submeter_cho2014.sh [nos = 4]: submete <nos> jobs cho2014_casos.sbatch (um no inteiro cada), que dividem os casos de
# casos_cho2014.txt (coluna de Cho 2014: deterministico e Monte Carlo direto) e os FORM + line sampling de
# casos_cho2014_ls.txt. Volta em segundos: pode desconectar. Pode rodar de novo a qualquer momento: os jobs novos pegam
# o que estiver livre.
# Antes, no no de login (segundos): as bases de Nystrom 1D de cada lv dos casos (nystrom_base.py coluna=2 dy=0.0025,
# com o python do venv), em $BASES. Sem jobs cho14 na fila, libera as travas de tarefas que nao terminaram (de jobs que
# morreram: fim do tempo, queda do no); as com ERRO ficam, a menos que RETENTAR=1.
# Ambiente: SAIDA ($DATA/cho2014), BASES ($SAIDA/bases), M (100), TEMPO (2-00:00:00), PARTICAO, DEPENDE (id de job:
# --dependency=afterany), PYV ($HOME/venv_la/bin/python) e as do sbatch (NY, EXTRA, PARALELO, COV, NMIN, N, FORM_IT,
# NMC_MAX, BLOCO_MC, CASOS, LS_CASOS), que vao junto (--export=ALL). O fermi roda 4 jobs por usuario e aceita 8 na fila:
# com outras campanhas rodando (ex. jobs liuB), submeta menos nos ou com DEPENDE.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
NOS="${1:-4}"
case "$NOS" in '' | *[!0-9]* | 0) echo "uso: $0 [numero de nos, padrao 4]" >&2; exit 2 ;; esac
export SAIDA="${SAIDA:-${DATA:-/mnt/data/${USER:-$(id -un)}}/cho2014}"
export BASES="${BASES:-$SAIDA/bases}" M="${M:-100}"
PYV="${PYV:-$HOME/venv_la/bin/python}"
EXE=bin_coluna/SlopeRainColumn
[ -x "$EXE" ] || { echo "compile antes, no no de login: ./compilar_coluna.sh" >&2; exit 1; }
# o executavel tem de ter o t_passagem localizado (coluna.cpp; o grep confere)
grep -q "t_passagem" "$EXE" || { echo "bin_coluna/ sem o SlopeRainColumn atual: rode ./compilar_coluna.sh" >&2; exit 1; }
[ -e "$SAIDA/PARAR" ] && { echo "existe $SAIDA/PARAR: apague antes (rm $SAIDA/PARAR)" >&2; exit 1; }
[ -x "$PYV" ] || { echo "falta o python com numpy e scipy ($PYV): as bases de Nystrom sao feitas com ele" >&2; exit 1; }
mkdir -p logs "$SAIDA" "$BASES"

# as bases: uma por lv dos casos (os nao comentados)
for lv in $(awk '!/^#/ && NF >= 5 { for (i = 5; i <= NF; i++) if ($i ~ /^lv=/) print substr($i, 4) }' casos_cho2014.txt | sort -u); do
    b="$BASES/base_lv${lv}_M$M.txt"
    [ -s "$b" ] && continue
    if "$PYV" -I ../../Cho2010/nystrom_base.py "$b.tmp" coluna=2 lv="$lv" dy=0.0025 M="$M" > "$BASES/base_lv${lv}_M$M.log" 2>&1; then
        mv -f "$b.tmp" "$b"
        echo "base lv = $lv: $(grep -m 1 'variancia retida com M = *'"$M"':' "$BASES/base_lv${lv}_M$M.log" | tr -s ' ')"
    else
        echo "a base lv = $lv FALHOU: veja $BASES/base_lv${lv}_M$M.log" >&2
        exit 1
    fi
done

fila="$(squeue -h -u "${USER:-$(id -un)}" -n cho14 -o %i)" || { echo "o squeue nao respondeu: tente de novo" >&2; exit 1; }
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
        case "$(basename "$d")" in bases | ls) continue ;; esac
        f=0; [ -s "$d/det/det_0_-1.csv" ] && [ "$(wc -l < "$d/det/det_0_-1.csv")" -gt 1 ] && f=1; libera "$d/.det.trava" $f
        for t in "$d"/mc/.tarefas/*; do
            [ -e "$t" ] || continue
            # o bloco esta feito se o log tem a linha 'arquivos:' (o executavel terminou)
            ab="$(basename "$t")" g="$d/mc/log_${ab}.txt" f=0
            [ -s "$d/mc/mc_${ab}.csv" ] && [ -s "$g" ] && grep -q 'arquivos:' "$g" && f=1
            libera "$t" $f
        done
    done
    for d in "$SAIDA"/ls/*/; do
        d="${d%/}"
        [ -d "$d" ] || continue
        f=0; [ -e "$d/FEITO" ] && f=1; libera "$d/.trava" $f
    done
    echo "sem jobs cho14 na fila: $n travas de tarefas interrompidas liberadas"
fi

for ((j = 1; j <= NOS; j++)); do
    # shellcheck disable=SC2086
    sbatch --export=ALL ${TEMPO:+--time=$TEMPO} ${PARTICAO:+--partition=$PARTICAO} \
        ${DEPENDE:+--dependency=afterany:$DEPENDE} cho2014_casos.sbatch
done
echo "saida: $SAIDA  (acompanhar: ./status_cho2014.sh)"
