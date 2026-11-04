#!/bin/bash
# rodar_caso.sh <n>: roda a linha n de casos.txt num unico processo de 1 CPU: transiente ate 48 h, depois os FS.
# Retomavel: rodar de novo pula o transiente ja concluido e os instantes ja gravados em fs.csv.
# Ambiente opcional: EXE (executavel), SAIDA (pasta de resultados), TF_H (horas, padrao 48).
# Tempos: cada fase executada acrescenta uma linha a tempos.csv (inicio, fim, segundos, maquina, job) e plano.txt guarda
# quantos passos e instantes de FS o caso deve ter. ./status.sh le esses arquivos, balanco.csv e fs.csv.
set -u
PKG="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXE="${EXE:-$PKG/bin/rahardjo}"
SAIDA="${SAIDA:-$PKG/resultados}"
TF="${TF_H:-48}"
linha="$(sed -n "${1:?informe o numero da linha de casos.txt}p" "$PKG/casos.txt")"
[ -n "$linha" ] || { echo "casos.txt nao tem a linha $1" >&2; exit 2; }
caso="${linha%% *}"; resto="${linha#* }"
# prioridade: t = 24 h, inicio, fim, depois preenche ate chegar a FS de 1 em 1 h
TODAS="$(seq -s, 0 "$TF")"
HORAS="24,0,12,48,6,18,30,36,42,3,9,15,21,27,33,39,45,$TODAS"
HORAS_USADAS=""; MESH="$PKG/malhas/rahardjo_a45.msh"; DT=0.1
for t in $resto; do
  case "$t" in
    HORAS=*) HORAS_USADAS="${t#HORAS=}" ;;
    MALHA=*) MESH="$PKG/malhas/${t#MALHA=}" ;;
    DT_H=*) DT="${t#DT_H=}" ;;
  esac
done
[ -n "$HORAS_USADAS" ] && HORAS="$HORAS_USADAS"
# remove instantes alem de TF e repetidos, mantendo a ordem
HORAS="$(echo "$HORAS" | tr ',' '\n' | awk -v tf="$TF" '$1<=tf && !v[$1]++' | paste -sd, -)"
# os estados gravados acompanham os instantes pedidos
DUMP="$(echo "$HORAS" | tr ',' '\n' | sort -n | paste -sd, -)"
export MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
d="$SAIDA/$caso"; mkdir -p "$d"; cd "$d" || exit 3
echo "passos=$(awk -v a="$TF" -v b="$DT" 'BEGIN { printf "%d", a / b + 0.5 }') fs=$(echo "$HORAS" | tr ',' '\n' | grep -c .) horas=$HORAS" > plano.txt
HOST="$(hostname -s 2>/dev/null || hostname)"
CPU="$(awk -F': *' '/^model name/ { print $2; exit }' /proc/cpuinfo 2>/dev/null | tr ',' ' ')"
JOB="${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-local}}${SLURM_ARRAY_TASK_ID:+_$SLURM_ARRAY_TASK_ID}"
[ -s tempos.csv ] || echo "fase,inicio,fim,seg,rc,host,job,cpu" > tempos.csv
marca() { echo "$1,$3,$(date '+%F %T'),$(($(date +%s) - $2)),$4,$HOST,$JOB,$CPU" >> tempos.csv; } # fase t0 inicio rc
ENVS=(MESH="$MESH" TF_H="$TF" CADA_FS=0)
for t in $resto; do case "$t" in HORAS=*|MALHA=*) ;; *) ENVS+=("$t") ;; esac; done
if ! grep -q "exit=0" exit.txt 2>/dev/null; then
  rm -f exit.txt
  t0=$(date +%s); i0="$(date '+%F %T')"
  env "${ENVS[@]}" DUMP_H="$DUMP" "$EXE" > run.log 2> run.err
  rc=$?
  echo "exit=$rc" > exit.txt
  marca transiente "$t0" "$i0" "$rc"
fi
# FS: valem os padroes abaixo, a menos que a linha de casos.txt os redefina
t0=$(date +%s); i0="$(date '+%F %T')"
env GIM_H=0,12,24,36,48 FS_PROJETO_H=0,6,12,18,24,30,36,42,48 VTK_H=24 "${ENVS[@]}" FS_ESTADOS="$HORAS" FS_CSV=fs.csv "$EXE" > fs.log 2> fs.err
rc=$?
marca fs "$t0" "$i0" "$rc"
echo "caso $caso: transiente $(cat exit.txt), FS rc=$rc, $(($(wc -l < fs.csv)-1)) instantes em fs.csv"
exit $rc
