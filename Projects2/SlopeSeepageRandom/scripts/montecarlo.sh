#!/bin/bash
# Monte Carlo em paralelo: divide as amostras entre N processos (cada um grava o seu CSV) e junta no fim.
# Uso: montecarlo.sh <executável> <diretório> <nome> <n amostras> <processos> [argumentos do comando mc ...]
set -u
EXE=$(readlink -f "$1"); OUT=$2; NOME=$3; N=$4; P=$5; shift 5
mkdir -p "$OUT"; cd "$OUT"
# autoproblema KL calculado uma vez (cache lido pelos processos)
"$EXE" mc "$@" n=0 saida="${NOME}_preparo.csv" > "${NOME}_preparo.log" 2>&1
per=$(( (N + P - 1) / P ))
for ((k = 0; k < P; k++)); do
    ini=$((k * per)); n=$per; [ $((ini + n)) -gt $N ] && n=$((N - ini)); [ $n -le 0 ] && continue
    csv="${NOME}_parte${k}.csv"
    # retomada: amostras já gravadas
    feitas=0; [ -s "$csv" ] && feitas=$(( $(wc -l < "$csv") - 1 ))
    [ $feitas -ge $n ] && continue
    "$EXE" mc "$@" inicio=$((ini + feitas)) n=$((n - feitas)) saida="$csv" > "${NOME}_parte${k}.log" 2>&1 &
done
wait
head -1 "${NOME}_parte0.csv" > "${NOME}.csv"
for f in ${NOME}_parte*.csv; do tail -n +2 "$f"; done | sort -t, -k1,1n >> "${NOME}.csv"
echo "$(($(wc -l < "${NOME}.csv") - 1)) amostras em ${NOME}.csv"
