#!/bin/bash
# rah.sbatch: as rodadas do Rahardjo2007 (casos.txt), cada uma serial na pasta $SAIDA/<id>/. PARALELO processos por no;
# cada um percorre a lista na ordem e pega o proximo caso livre criando $SAIDA/<id>/.trava com o noclobber do bash
# (atomico, tambem entre nos e entre jobs). Feito: existe $SAIDA/<id>/fim (o codigo de saida em codigo). Sem caso livre,
# ou com o arquivo $SAIDA/PARAR, o processo termina. Resumo: ./rah.sbatch resumo. Uso: ./submeter.sh [nos].
# Ambiente: SAIDA ($DATA/rahardjo2007), EXE (bin/Rahardjo2007), CASOS_ARQ (casos.txt), PARALELO (7 processos por no:
# 4 nos x 7 = os 28 casos de uma vez), THREADS (3 threads de montagem por processo: 21 dos 24 do no)
#SBATCH --job-name=rah2007
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=60G
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/%x_%j.out
[ -r /etc/profile.d/scratch.sh ] && source /etc/profile.d/scratch.sh
cd "${SLURM_SUBMIT_DIR:-$(dirname "$(readlink -f "$0")")}" || exit 1
[ -n "${SLURM_JOB_ID:-}" ] && source ../../LearningMohrCoulomb2/cluster/modulos.sh
DADOS="${DATA:-/mnt/data/${USER:-$(id -un)}}"
SAIDA="${SAIDA:-$DADOS/rahardjo2007}"
EXE="${EXE:-$PWD/bin/Rahardjo2007}"
CASOS_ARQ="${CASOS_ARQ:-$PWD/casos.txt}"
PARALELO="${PARALELO:-7}"
THREADS="${THREADS:-3}"

resumo() {
    local tot=0 feito=0 rod=0 erro=0
    while read -r id _; do
        [[ -z "$id" || "$id" == \#* ]] && continue
        tot=$((tot + 1))
        if [ -f "$SAIDA/$id/fim" ]; then
            feito=$((feito + 1))
            [ "$(cat "$SAIDA/$id/codigo" 2> /dev/null)" = 0 ] || { erro=$((erro + 1)); echo "  erro: $id (codigo $(cat "$SAIDA/$id/codigo"))"; }
        elif [ -f "$SAIDA/$id/.trava" ] && [ ! -s "$SAIDA/$id/.trava" ] && [ ! -f "$SAIDA/$id/inicio" ]; then
            # trava vazia e sem inicio: a criacao falhou no meio e ninguem roda o caso (o submeter.sh, sem job na fila, libera)
            echo "  trava vazia (ninguem rodando): $id"
        elif [ -f "$SAIDA/$id/.trava" ]; then
            rod=$((rod + 1))
            echo "  rodando: $id ($(tail -n 1 "$SAIDA/$id"/rah_*_fs.csv 2> /dev/null | cut -d, -f1) h)"
        fi
    done < "$CASOS_ARQ"
    echo "casos $tot: feitos $feito (com erro $erro), rodando $rod, livres $((tot - feito - rod))"
}
if [ "${1:-}" = resumo ]; then resumo; exit 0; fi

[ -x "$EXE" ] || { echo "falta $EXE (rode ./compilar.sh)"; exit 1; }
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PZ_THREADS="$THREADS" LD_LIBRARY_PATH="$(dirname "$EXE"):${LD_LIBRARY_PATH:-}"
mkdir -p "$SAIDA"
echo "$(date '+%F %T') [${SLURM_JOB_ID:-local}] inicio em $(hostname): $PARALELO processos de $THREADS threads; $(cat "$(dirname "$EXE")/versao.txt" 2> /dev/null)" >> "$SAIDA/eventos.log"
trabalhador() {
    local id opc d
    while read -r id opc; do
        [[ -z "$id" || "$id" == \#* ]] && continue
        [ -e "$SAIDA/PARAR" ] && return
        d="$SAIDA/$id"
        [ -f "$d/fim" ] && continue
        mkdir -p "$d" 2> /dev/null
        (set -o noclobber; echo "${SLURM_JOB_ID:-local} $(hostname) $(date '+%F %T')" > "$d/.trava") 2> /dev/null || continue
        (
            cd "$d" || exit 1
            echo "$EXE $opc saida=rah" > comando.txt
            date +%s > inicio
            # shellcheck disable=SC2086
            "$EXE" $opc saida=rah > log.txt 2>&1
            echo $? > codigo
            date +%s > fim
        )
        echo "$(date '+%F %T') [${SLURM_JOB_ID:-local}] $id codigo $(cat "$d/codigo") ($(( $(cat "$d/fim") - $(cat "$d/inicio") )) s)" >> "$SAIDA/eventos.log"
    done < "$CASOS_ARQ"
}
for ((p = 0; p < PARALELO; p++)); do
    trabalhador &
    sleep 0.3
done
wait
echo "$(date '+%F %T') [${SLURM_JOB_ID:-local}] fim em $(hostname)" >> "$SAIDA/eventos.log"
resumo
