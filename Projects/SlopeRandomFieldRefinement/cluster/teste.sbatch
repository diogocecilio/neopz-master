#!/bin/bash
# teste.sbatch: teste rapido do ambiente antes da campanha: o deterministico e 2 realizacoes da referencia, em ./teste
# (chuva, uns 15 min: com 0 e 24 h) ou ./teste_ceron (CAMPANHA=ceron, uns 5 min). Submeter com  sbatch teste.sbatch  (ou
# CAMPANHA=ceron sbatch teste.sbatch) e ver  cat logs/chuvaTeste_<job>.out
#SBATCH --job-name=chuvaTeste
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x_%j.out
[ -r /etc/profile.d/scratch.sh ] && source /etc/profile.d/scratch.sh
cd "$SLURM_SUBMIT_DIR" || exit 1
source ../../LearningMohrCoulomb2/cluster/modulos.sh
source ./campanha.sh
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PZ_THREADS=0
if [ "$CAMPANHA" = ceron ]; then
    mkdir -p teste_ceron
    cd teste_ceron || exit 1
    echo "$(date '+%F %T') $(hostname -s), versao $(cat ../bin/versao.txt 2> /dev/null || echo -)"
    ../bin/SlopeDrawdownMonteCarlo n=2 det=1 saida=mc
    echo "$(date '+%F %T') fim: o deterministico deve dar Gamma ~1,325 (1,3253 no PC; o artigo da 1,336)"
else
    mkdir -p teste
    cd teste || exit 1
    echo "$(date '+%F %T') $(hostname -s), versao $(cat ../bin/versao.txt 2> /dev/null || echo -)"
    ../bin/SlopeRainMonteCarlo n=2 det=1 I="${IREF:-31.05}" tempos=0,24 saida=mc
    echo "$(date '+%F %T') fim: o deterministico deve dar FS(0 h) ~2,094 e FS(24 h) ~1,126 com I = 31,05 mm/h"
fi
