#!/bin/bash
# Monte Carlo em paralelo: divide as realizacoes 0 .. N-1 em P blocos continuos, um processo de 1 thread por bloco, na
# pasta corrente (os CSV e os logs ficam nela). Retomavel e extensivel: cada processo pula as realizacoes que ja estao em
# algum CSV da mesma configuracao na pasta (rodar de novo com N maior acrescenta realizacoes). Espera os processos e
# imprime o resumo de tudo.
# Uso: rodar_mc.sh [N = 10000] [P = nucleos fisicos] [nunif = 0] [nref = 4] [pmax = 2]
#      EXE=<executavel> rodar_mc.sh ...   (padrao: o da pasta de build release ao lado do repositorio)
N=${1:-10000}
P=${2:-$(lscpu -p=core,socket | grep -v '^#' | sort -u | wc -l)}
NUNIF=${3:-0}; NREF=${4:-4}; PMAX=${5:-2}
DIR=$(cd "$(dirname "$0")" && pwd)
EXE=${EXE:-$DIR/../../../neopz-master-build-release/Projects/SlopeMonteCarlo/SlopeMonteCarlo}
[ -x "$EXE" ] || { echo "executavel nao encontrado: $EXE (defina EXE)"; exit 1; }
BLOCO=$(( (N + P - 1) / P ))
echo "$N realizacoes em blocos de $BLOCO, ate $P processos | configuracao u${NUNIF}h${NREF}p${PMAX} | $(pwd)"
for ((p = 0; p < P; p++)); do
    PRIMEIRA=$((p * BLOCO))
    [ "$PRIMEIRA" -ge "$N" ] && break
    NB=$(( BLOCO < N - PRIMEIRA ? BLOCO : N - PRIMEIRA ))
    OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 "$EXE" "$NB" "$PRIMEIRA" "$NUNIF" "$NREF" "$PMAX" > "log_u${NUNIF}h${NREF}p${PMAX}_${PRIMEIRA}.txt" 2>&1 &
done
echo "acompanhe com: tail -qn 1 log_u${NUNIF}h${NREF}p${PMAX}_*.txt"
wait
"$EXE" resumo "$NUNIF" "$NREF" "$PMAX"
