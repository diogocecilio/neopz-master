#!/bin/bash
# bissecao_ls.sbatch: confere por bissecao as raizes das linhas de menor c de um line sampling pronto
# (../../Cho2010/bissecao_ls.py): o z de cada linha e regerado da semente e do indice da linha (os .xi nao ficam), a
# 1a analise refaz a 1a analise da linha no analises.csv (conferencia da reconstrucao), a cerca e c_ls +- MEIA (mais os
# pontos EXTRAS, que devem dar seguro em c = 0 e falha em c = 4) e a bissecao vai ate |dc| < TOL. As linhas rodam em
# paralelo (uma analise por vez cada, com THREADS threads do PZ). Padrao: as 10 linhas de menor raiz do alfa=4 da
# campanha ls_ceron_casos.sbatch (Ceron 2025, Secao 6.4: o LS da 1,42% contra 1,00% do artigo; o teste descarta raizes
# espurias do passo de Newton da secante).
# Uso (no fermi, na pasta cluster): sbatch bissecao_ls.sbatch. Ambiente: PASTA_LS ($DATA/ceron_ls_casos/alfa=4/ls),
# SAIDA ($DATA/ceron_ls_testes/bissecao_alfa4), LINHAS (10), LISTA (vazio: as LINHAS de menor c; ou indices separados
# por virgulas), TOL (0.001), MEIA (0.25), EXTRAS (0,4), THREADS (2), EXE (bin_ceron/SlopeDrawdownMonteCarlo).
# Saidas em $SAIDA: bissecao.csv, bissecao.json, analises.csv e bissecao.log.
#SBATCH --job-name=bisLS
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=24G
#SBATCH --time=0-06:00:00
#SBATCH --output=logs/%x_%j.out
[ -r /etc/profile.d/scratch.sh ] && source /etc/profile.d/scratch.sh
cd "${SLURM_SUBMIT_DIR:-$(dirname "$(readlink -f "$0")")}" || exit 1
[ -n "${SLURM_JOB_ID:-}" ] && source ../../LearningMohrCoulomb2/cluster/modulos.sh
DADOS="${DATA:-/mnt/data/${USER:-$(id -un)}}"
PASTA_LS="${PASTA_LS:-$DADOS/ceron_ls_casos/alfa=4/ls}"
SAIDA="${SAIDA:-$DADOS/ceron_ls_testes/bissecao_alfa4}"
LINHAS="${LINHAS:-10}" LISTA="${LISTA:-}" TOL="${TOL:-0.001}" MEIA="${MEIA:-0.25}" EXTRAS="${EXTRAS:-0,4}"
THREADS="${THREADS:-2}" EXE="${EXE:-$PWD/bin_ceron/SlopeDrawdownMonteCarlo}" PY="$(cd ../../Cho2010 && pwd)"
[ -x "$EXE" ] || { echo "falta $EXE: rode ./compilar_ceron.sh no no de login" >&2; exit 1; }
[ -s "$PASTA_LS/ls.json" ] || { echo "falta $PASTA_LS/ls.json (o line sampling pronto)" >&2; exit 1; }
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p "$SAIDA"
echo "$(date '+%F %T') bissecao de $PASTA_LS em $(hostname -s), executavel $(cat "$(dirname "$EXE")/versao.txt" 2> /dev/null)"
python3 -I "$PY/bissecao_ls.py" "$EXE" "$PASTA_LS" "$SAIDA" linhas="$LINHAS" ${LISTA:+lista="$LISTA"} tol="$TOL" \
    meia="$MEIA" extras="$EXTRAS" paralelos="${SLURM_CPUS_PER_TASK:-24}" threads="$THREADS" 2>&1 | tee "$SAIDA/bissecao.log"
echo "$(date '+%F %T') fim"
