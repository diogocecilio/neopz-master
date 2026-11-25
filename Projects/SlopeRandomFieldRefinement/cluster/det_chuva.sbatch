#!/bin/bash
# det_chuva.sbatch: a validacao deterministica da chuva (Projects/PlanoDeTrabalho/chuva2/plano_deterministico.md): os
# casos de casos_det_chuva.txt (gera_casos_det_chuva.py), cada um uma rodada serial do SlopeRainMonteCarlo (det=1 n=0,
# FS de hora em hora, diag=diag.csv, grava_pressoes=pw.txt) na pasta $SAIDA/<id>/. PARALELO processos por no; cada um
# percorre a lista na ordem (a prioridade do arquivo) e pega o proximo caso livre criando $SAIDA/<id>/.trava com o
# noclobber do bash (atomico, tambem entre nos e entre jobs). Feito: existe $SAIDA/<id>/fim (codigo em codigo). Sem caso
# livre, ou com o arquivo $SAIDA/PARAR, o processo termina.
# Retomada: o submeter_det_chuva.sh, sem jobs detChuva na fila, apaga as travas dos casos sem fim (jobs que morreram).
# Resumo: ./det_chuva.sbatch resumo (feitos, rodando, livres, com erro). Uso: ./submeter_det_chuva.sh [nos].
# Ambiente: SAIDA ($DATA/det_chuva), EXE (bin_det/SlopeRainMonteCarlo), CASOS_ARQ (casos_det_chuva.txt), PARALELO (24)
#SBATCH --job-name=detChuva
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=60G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/%x_%j.out
[ -r /etc/profile.d/scratch.sh ] && source /etc/profile.d/scratch.sh
cd "${SLURM_SUBMIT_DIR:-$(dirname "$(readlink -f "$0")")}" || exit 1
[ -n "${SLURM_JOB_ID:-}" ] && source ../../LearningMohrCoulomb2/cluster/modulos.sh
DADOS="${DATA:-/mnt/data/${USER:-$(id -un)}}"
SAIDA="${SAIDA:-$DADOS/det_chuva}"
EXE="${EXE:-$PWD/bin_det/SlopeRainMonteCarlo}"
CASOS_ARQ="${CASOS_ARQ:-$PWD/casos_det_chuva.txt}"
PARALELO="${PARALELO:-24}"
COMUM="det=1 n=0 covc=0 covphi=0 covk=0 saida=det diag=diag.csv grava_pressoes=pw.txt"

resumo() {
    local tot=0 feito=0 rod=0 erro=0
    while read -r id _; do
        [[ -z "$id" || "$id" == \#* ]] && continue
        tot=$((tot + 1))
        if [ -f "$SAIDA/$id/fim" ]; then
            feito=$((feito + 1))
            [ "$(cat "$SAIDA/$id/codigo" 2> /dev/null)" = 0 ] || { erro=$((erro + 1)); echo "  erro: $id (codigo $(cat "$SAIDA/$id/codigo"))"; }
        elif [ -f "$SAIDA/$id/.trava" ]; then rod=$((rod + 1)); fi
    done < "$CASOS_ARQ"
    echo "casos $tot: feitos $feito (com erro $erro), rodando $rod, livres $((tot - feito - rod))"
}
if [ "${1:-}" = resumo ]; then resumo; exit 0; fi

[ -x "$EXE" ] || { echo "falta $EXE (rode ./compilar_det_chuva.sh)"; exit 1; }
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PZ_THREADS=1 LD_LIBRARY_PATH="$(dirname "$EXE"):${LD_LIBRARY_PATH:-}"
mkdir -p "$SAIDA"
echo "$(date '+%F %T') [${SLURM_JOB_ID:-local}] inicio em $(hostname): $PARALELO processos; $(cat "$(dirname "$EXE")/versao.txt" 2> /dev/null)" >> "$SAIDA/eventos.log"
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
            echo "$EXE $COMUM $opc" > comando.txt
            date +%s > inicio
            # shellcheck disable=SC2086
            "$EXE" $COMUM $opc > log.txt 2>&1
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
