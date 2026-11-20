#!/bin/bash
# sus.sbatch: simulacao por subconjuntos (Projects/Cho2010/sus.py) do Cho2010 com r = -0,7 (o mesmo modelo do r07_direto:
# KL analitica com M = 2000, SRM, fator=tri2 nref=4), uma rodada independente por semente de SEMENTES, uma depois da
# outra, com PARALELO analises ao mesmo tempo. N = 500 por nivel e p0 = 0,1: ~1400 analises por rodada (~1,5 h no no).
# Uso (no fermi, nesta pasta):  SEMENTES="1 2 3 4" sbatch sus.sbatch      (resultados em $DATA/sus_r07/s<semente>)
#SBATCH --job-name=susCho
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=32G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/%x_%j.out
[ -r /etc/profile.d/scratch.sh ] && source /etc/profile.d/scratch.sh
cd "$SLURM_SUBMIT_DIR" || exit 1
source ../../LearningMohrCoulomb2/cluster/modulos.sh
: "${SEMENTES:?SEMENTES=\"1 2 3 4\" sbatch sus.sbatch}"
SAIDA="${SAIDA:-${DATA:-$PWD}/sus_r07}" PARALELO="${PARALELO:-24}" N="${N:-500}"
EXE="$PWD/bin/Cho2010"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PZ_THREADS=0
for s in $SEMENTES; do
    echo "$(date '+%F %T') semente $s em $(hostname -s)"
    python3 -I ../../Cho2010/sus.py "$EXE" "$SAIDA/s$s" prog=cho M=2000 rcphi=-0.7 fator=tri2 nref=4 metodo=srm \
        N="$N" p0=0.1 rho=0.8 semente="$s" paralelos="$PARALELO"
done
echo "$(date '+%F %T') fim"
