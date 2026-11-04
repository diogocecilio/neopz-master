#!/bin/bash
# Todos os casos num unico job, que conta como 1 job na cota do SLURM: roda ate SLURM_CPUS_PER_TASK casos ao mesmo tempo,
# cada um num processo de 1 CPU, e comeca o proximo da lista quando um termina. Submetido por submeter.sh, que define
# --cpus-per-task, --mem e a lista de casos (CASOS, por exemplo 1-21 ou 3,7).
#SBATCH --job-name=rahardjo45
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=6G
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/%x_%j.out
cd "$SLURM_SUBMIT_DIR" || exit 1
source ./modulos.sh
# o executavel foi compilado no no de login: confere que este no acha todas as bibliotecas (MKL, libpz, libstdc++)
if ! [ -x bin/rahardjo ] || ldd bin/rahardjo | grep -q "not found"; then
  echo "$(hostname -s): bin/rahardjo ausente ou com bibliotecas faltando neste no:" >&2
  [ -x bin/rahardjo ] && ldd bin/rahardjo | grep "not found" >&2
  exit 1
fi
K="${SLURM_CPUS_PER_TASK:-1}"
# "1-3,7" -> 1 2 3 7
casos="$(echo "${CASOS:-1-21}" | tr ',' '\n' | awk -F- 'NF == 2 { for (i = $1; i <= $2; i++) print i; next } NF { print $1 }')"
mkdir -p logs/casos
echo "$(date '+%F %T') job $SLURM_JOB_ID em $(hostname -s): casos $(echo $casos), ate $K ao mesmo tempo"
for n in $casos; do
  while [ "$(jobs -rp | wc -l)" -ge "$K" ]; do wait -n; done
  (
    echo "$SLURM_JOB_ID $(date '+%F %T')" > "logs/casos/$n.ativo"
    ./rodar_caso.sh "$n" > "logs/casos/$n.out" 2>&1
    echo "$(date '+%F %T') $(tail -n 1 "logs/casos/$n.out")"
    rm -f "logs/casos/$n.ativo"
  ) &
done
wait
echo "$(date '+%F %T') fim"
