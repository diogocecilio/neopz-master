#!/bin/bash
# v3.sbatch: os casos completos do V3 (Sharma et al. 2021: barragem retangular, face inclinada e dois blocos) e o
# estudo de sensibilidade a curva de retencao e a malha, no SLURM. Cada execucao roda num processo, na sua pasta
# ($SAIDA/<nome>), ao mesmo tempo (um nucleo cada; de 30 s a ~10 min cada num nucleo de 2,1 GHz); no fim junta os
# resumo.csv em $SAIDA/resumo.csv, roda os testes [.lento] se TESTES=1 e a figura de comparacao se houver python3 com
# numpy e matplotlib.
# Uso, de dentro desta pasta, depois de compilar (compilar_v3.sh):
#   sbatch v3.sbatch                                  # os tres casos e a sensibilidade (11 execucoes)
#   EXECUCOES="retangular inclinada" TESTES=1 sbatch v3.sbatch
# Ambiente: BUILD (o mesmo de compilar_v3.sh), SAIDA (padrao $DATA/v3 ou ./resultados_v3), EXECUCOES (nomes da lista
# abaixo), TESTES. Cada nome define as variaveis do executavel (CASO, NX, NY, AG, ...; ver main.cpp e o README).
#SBATCH --job-name=v3_freatica
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=6G
#SBATCH --time=03:00:00
#SBATCH --output=v3_%j.out
[ -r /etc/profile.d/scratch.sh ] && source /etc/profile.d/scratch.sh
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}" || exit 1
AQUI="$PWD"
NEOPZ="$(cd "$AQUI/../../.." && pwd)"
source "$NEOPZ/Projects/Cluster/modulos.sh"
export MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
BUILD="${BUILD:-$NEOPZ/build_verif}"
EXE="$BUILD/Projects/V3_SuperficieLivre/V3_SuperficieLivre"
SAIDA="${SAIDA:-${DATA:+$DATA/v3}}"
SAIDA="${SAIDA:-$AQUI/resultados_v3}"
# nome -> variaveis do executavel
declare -A CONF=(
  [retangular]="CASO=retangular"
  [inclinada]="CASO=inclinada"
  [dois_blocos]="CASO=dois_blocos"
  [retangular_a05]="CASO=retangular AG=0.5"
  [retangular_a2]="CASO=retangular AG=2"
  [retangular_a4]="CASO=retangular AG=4"
  [retangular_fina]="CASO=retangular NX=64 NY=96"
  [retangular_linear]="CASO=retangular CI=linear"
  [inclinada_a2]="CASO=inclinada AG=2"
  [inclinada_a4]="CASO=inclinada AG=4"
  [inclinada_fina]="CASO=inclinada NX=56 NY=80"
)
EXECUCOES="${EXECUCOES:-retangular inclinada dois_blocos retangular_a05 retangular_a2 retangular_a4 retangular_fina retangular_linear inclinada_a2 inclinada_a4 inclinada_fina}"
if ! [ -x "$EXE" ] || ldd "$EXE" | grep -q "not found"; then
    echo "$EXE ausente ou com bibliotecas faltando: rode compilar_v3.sh antes" >&2
    exit 1
fi
mkdir -p "$SAIDA"
echo "$(date '+%F %T') $(hostname -s): execucoes $EXECUCOES, saida $SAIDA"
# no maximo SLURM_CPUS_PER_TASK processos ao mesmo tempo
MAXP="${SLURM_CPUS_PER_TASK:-4}"
for e in $EXECUCOES; do
    [ -n "${CONF[$e]:-}" ] || { echo "execucao desconhecida: $e" >&2; continue; }
    while [ "$(jobs -rp | wc -l)" -ge "$MAXP" ]; do sleep 5; done
    mkdir -p "$SAIDA/$e"
    (cd "$SAIDA/$e" && env ${CONF[$e]} VERBOSE=1 "$EXE" > log.txt 2>&1; echo "$e: codigo $?" >> log.txt) &
done
wait
# junta os resumo.csv (o cabecalho da primeira execucao) com o nome de cada uma
primeiro=1
for e in $EXECUCOES; do
    [ -f "$SAIDA/$e/resumo.csv" ] || continue
    if [ "$primeiro" = 1 ]; then echo "execucao,$(head -n 1 "$SAIDA/$e/resumo.csv")" > "$SAIDA/resumo.csv"; primeiro=0; fi
    echo "$e,$(tail -n 1 "$SAIDA/$e/resumo.csv")" >> "$SAIDA/resumo.csv"
    grep -h "codigo" "$SAIDA/$e/log.txt"
done
if [ "${TESTES:-0}" = 1 ]; then
    (cd "$SAIDA" && "$BUILD/Projects/V3_SuperficieLivre/V3_testes" "[V3]" > testes.txt 2>&1;
     "$BUILD/Projects/V3_SuperficieLivre/V3_testes" "[.lento]" > testes_lento.txt 2>&1)
    tail -n 3 "$SAIDA/testes.txt" "$SAIDA/testes_lento.txt"
fi
if command -v python3 > /dev/null && python3 -c "import numpy, matplotlib" 2> /dev/null; then
    python3 "$AQUI/../scripts/plot_v3.py" "$SAIDA" --saida "$SAIDA/v3_comparacao.png" > "$SAIDA/comparacao.txt"
fi
echo "$(date '+%F %T') fim: $SAIDA"
