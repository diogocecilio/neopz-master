#!/bin/bash
# Job array: uma tarefa de 1 CPU por linha de casos.txt. Submeter de dentro desta pasta: sbatch job_slurm.sbatch
#SBATCH --job-name=rahardjo45
#SBATCH --array=1-21
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x_%A_%a.out
# particao e modulos: em modulos.sh (submeter.sh passa a particao)
cd "$SLURM_SUBMIT_DIR" || exit 1
source ./modulos.sh
# o executavel foi compilado no no de login: confere que este no acha todas as bibliotecas (MKL, libpz, libstdc++)
if ! [ -x bin/rahardjo ] || ldd bin/rahardjo | grep -q "not found"; then
  echo "$(hostname -s): bin/rahardjo ausente ou com bibliotecas faltando neste no:" >&2
  [ -x bin/rahardjo ] && ldd bin/rahardjo | grep "not found" >&2
  exit 1
fi
./rodar_caso.sh "$SLURM_ARRAY_TASK_ID"
