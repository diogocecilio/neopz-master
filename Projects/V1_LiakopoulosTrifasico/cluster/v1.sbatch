#!/bin/bash
# v1.sbatch: o caso completo do V1 (coluna de Liakopoulos ate 600 min) no SLURM. Cada caso roda num processo, na sua
# pasta, ao mesmo tempo (um nucleo cada; cerca de 3 min por caso num nucleo de 2,1 GHz); no fim junta os CSVs em
# $SAIDA, roda os testes [.lento] se TESTES=1 e a figura de comparacao se houver python3 com numpy e matplotlib.
# Uso, de dentro desta pasta, depois de compilar (compilar_v1.sh):
#   sbatch v1.sbatch                       # os 4 casos, 600 min
#   CASOS="bifasico trifasico_alert" T_FIM=7200 TESTES=1 sbatch v1.sbatch
# Ambiente: BUILD (o mesmo de compilar_v1.sh), SAIDA (padrao $DATA/v1 ou ./resultados_v1), CASOS, T_FIM, NY,
# DT_FATOR, TESTES. As variaveis CASO, T_FIM, NY e DT_FATOR sao as do executavel (ver main.cpp e o README).
#SBATCH --job-name=v1_liako
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=v1_%j.out
[ -r /etc/profile.d/scratch.sh ] && source /etc/profile.d/scratch.sh
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}" || exit 1
AQUI="$PWD"
NEOPZ="$(cd "$AQUI/../../.." && pwd)"
source "$NEOPZ/Projects/Cluster/modulos.sh"
export MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
BUILD="${BUILD:-$NEOPZ/build_verif}"
EXE="$BUILD/Projects/V1_LiakopoulosTrifasico/V1_LiakopoulosTrifasico"
SAIDA="${SAIDA:-${DATA:+$DATA/v1}}"
SAIDA="${SAIDA:-$AQUI/resultados_v1}"
CASOS="${CASOS:-bifasico trifasico trifasico_alert trifasico_laloui}"
if ! [ -x "$EXE" ] || ldd "$EXE" | grep -q "not found"; then
    echo "$EXE ausente ou com bibliotecas faltando: rode compilar_v1.sh antes" >&2
    exit 1
fi
mkdir -p "$SAIDA"
echo "$(date '+%F %T') $(hostname -s): casos $CASOS, T_FIM=${T_FIM:-36000}, saida $SAIDA"
for c in $CASOS; do
    mkdir -p "$SAIDA/$c"
    (cd "$SAIDA/$c" && CASO="$c" "$EXE" > log.txt 2>&1; echo "$c: codigo $?" >> log.txt) &
done
wait
# junta os CSVs (o cabecalho do primeiro caso)
primeiro=1
for c in $CASOS; do
    for f in historico perfis resumo; do
        if [ "$primeiro" = 1 ]; then cp "$SAIDA/$c/$f.csv" "$SAIDA/$f.csv"; else tail -n +2 "$SAIDA/$c/$f.csv" >> "$SAIDA/$f.csv"; fi
    done
    primeiro=0
    grep -h "codigo\|passos" "$SAIDA/$c/log.txt"
done
if [ "${TESTES:-0}" = 1 ]; then
    (cd "$SAIDA" && "$BUILD/Projects/V1_LiakopoulosTrifasico/V1_testes" "[V1]" > testes.txt 2>&1;
     "$BUILD/Projects/V1_LiakopoulosTrifasico/V1_testes" "[.lento]" > testes_lento.txt 2>&1)
    tail -n 3 "$SAIDA/testes.txt" "$SAIDA/testes_lento.txt"
fi
if command -v python3 > /dev/null && python3 -c "import numpy, matplotlib" 2> /dev/null; then
    python3 "$AQUI/../scripts/plot_v1.py" "$SAIDA" --saida "$SAIDA/v1_comparacao.png" > "$SAIDA/comparacao.txt"
fi
echo "$(date '+%F %T') fim: $SAIDA"
