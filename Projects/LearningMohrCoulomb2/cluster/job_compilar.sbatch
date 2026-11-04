#!/bin/bash
# Compila num no de computacao. Submetido por submeter.sh.
#SBATCH --job-name=rahardjo45_build
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%j.out
cd "$SLURM_SUBMIT_DIR" || exit 1
source ./modulos.sh
JOBS="${SLURM_CPUS_PER_TASK:-4}" ./compilar.sh
