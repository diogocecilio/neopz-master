#!/bin/bash
# v6.sbatch: as execucoes completas do V6 (Yuan et al. 2021: media do campo KL num arco) no SLURM. Cada execucao roda
# num processo, na sua pasta ($SAIDA/<nome>), ao mesmo tempo (um nucleo cada); no fim junta os resumo.csv em
# $SAIDA/resumo.csv, roda os testes (padrao e [.lento]) se TESTES=1 e as figuras se houver python3 com numpy e
# matplotlib. Tempos num nucleo de ~2 GHz: H = 10 m ~1 min; H = 50 m com 90% da variancia e 10^5 realizacoes ~20 min;
# com 95% (M ~ 7,5e5) e 2 10^4 realizacoes ~10 min; memoria de ~180 MB a ~360 MB por execucao (95%).
# Uso, de dentro desta pasta, depois de compilar (compilar_v6.sh):
#   sbatch v6.sbatch                               # todas as execucoes da lista (17)
#   EXECUCOES="tabelas fig15 b40_h50" TESTES=1 sbatch v6.sbatch
# Ambiente: BUILD (o mesmo de compilar_v6.sh), SAIDA (padrao $DATA/v6 ou ./resultados_v6), EXECUCOES (nomes da lista
# abaixo), N (realizacoes das execucoes com 90%; padrao 100000), TESTES. Cada nome define as variaveis do executavel
# (CASO, CONVENCAO, FRACAO, MKL, N, SEMENTE; ver main.cpp e o README).
#SBATCH --job-name=v6_klarco
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=12G
#SBATCH --time=04:00:00
#SBATCH --output=v6_%j.out
[ -r /etc/profile.d/scratch.sh ] && source /etc/profile.d/scratch.sh
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}" || exit 1
AQUI="$PWD"
NEOPZ="$(cd "$AQUI/../../.." && pwd)"
source "$NEOPZ/Projects/Cluster/modulos.sh"
export MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
BUILD="${BUILD:-$NEOPZ/build_verif}"
EXE="$BUILD/Projects/V6_MediaLocalKL/V6_MediaLocalKL"
SAIDA="${SAIDA:-${DATA:+$DATA/v6}}"
SAIDA="${SAIDA:-$AQUI/resultados_v6}"
NR="${N:-100000}"
# nome -> variaveis do executavel (sementes diferentes em cada execucao)
declare -A CONF=(
  [tabelas]="CASO=tabelas"
  [tabelas_delta]="CASO=tabelas CONVENCAO=delta"
  [b40_h10]="CASO=b40_h10 N=$NR SEMENTE=101"
  [b40_h20]="CASO=b40_h20 N=$NR SEMENTE=102"
  [b40_h30]="CASO=b40_h30 N=$NR SEMENTE=103"
  [b40_h40]="CASO=b40_h40 N=$NR SEMENTE=104"
  [b40_h50]="CASO=b40_h50 N=$NR SEMENTE=105"
  [h50_b30]="CASO=h50_b30 N=$NR SEMENTE=106"
  [h50_b35]="CASO=h50_b35 N=$NR SEMENTE=107"
  [h50_b45]="CASO=h50_b45 N=$NR SEMENTE=108"
  [h50_b50]="CASO=h50_b50 N=$NR SEMENTE=109"
  [fig15]="CASO=fig15 N=$NR SEMENTE=110"
  [fig15_delta]="CASO=fig15 CONVENCAO=delta N=$NR SEMENTE=111"
  [fig15_f95]="CASO=fig15 FRACAO=0.95 MKL=900000 N=20000 SEMENTE=112"
  [b40_h50_f95]="CASO=b40_h50 FRACAO=0.95 MKL=900000 N=20000 SEMENTE=113"
  [b40_h10_f99]="CASO=b40_h10 FRACAO=0.99 MKL=900000 N=20000 SEMENTE=114"
  [b40_h10_n1e6]="CASO=b40_h10 N=1000000 SEMENTE=115"
)
EXECUCOES="${EXECUCOES:-tabelas tabelas_delta b40_h50_f95 fig15_f95 b40_h10_f99 b40_h50 h50_b30 h50_b35 h50_b45 h50_b50 fig15 fig15_delta b40_h40 b40_h30 b40_h20 b40_h10 b40_h10_n1e6}"
if ! [ -x "$EXE" ] || ldd "$EXE" | grep -q "not found"; then
    echo "$EXE ausente ou com bibliotecas faltando: rode compilar_v6.sh antes" >&2
    exit 1
fi
mkdir -p "$SAIDA"
echo "$(date '+%F %T') $(hostname -s): execucoes $EXECUCOES, saida $SAIDA"
MAXP="${SLURM_CPUS_PER_TASK:-4}" # no maximo um processo por CPU do job
TEMPO=""; [ -x /usr/bin/time ] && TEMPO="/usr/bin/time -v" # tempo e memoria de cada execucao (tempo.txt), se houver
for e in $EXECUCOES; do
    [ -n "${CONF[$e]:-}" ] || { echo "execucao desconhecida: $e" >&2; continue; }
    [ -f "$SAIDA/$e/resultado.txt" ] && { echo "$e: ja feita (resultado.txt), pulada"; continue; }
    while [ "$(jobs -rp | wc -l)" -ge "$MAXP" ]; do sleep 5; done
    mkdir -p "$SAIDA/$e"
    (cd "$SAIDA/$e" && env ${CONF[$e]} VERBOSE=1 $TEMPO "$EXE" > log.txt 2> tempo.txt; echo "$e: codigo $?" >> log.txt) &
done
wait
# junta os resumo.csv das execucoes com KL (o cabecalho da primeira) com o nome de cada uma
primeiro=1
for e in $EXECUCOES; do
    [ -f "$SAIDA/$e/resumo.csv" ] || continue
    case "${CONF[$e]}" in *CASO=tabelas*) continue ;; esac
    if [ "$primeiro" = 1 ]; then echo "execucao,$(head -n 1 "$SAIDA/$e/resumo.csv")" > "$SAIDA/resumo.csv"; primeiro=0; fi
    echo "$e,$(tail -n 1 "$SAIDA/$e/resumo.csv")" >> "$SAIDA/resumo.csv"
done
for e in $EXECUCOES; do grep -h "codigo" "$SAIDA/$e/log.txt" 2>/dev/null; done
if [ "${TESTES:-0}" = 1 ]; then
    # sem filtro: os testes padrao ("[V6]" selecionaria tambem os ocultos [.lento])
    (cd "$SAIDA" && "$BUILD/Projects/V6_MediaLocalKL/V6_testes" --durations yes > testes.txt 2>&1;
     "$BUILD/Projects/V6_MediaLocalKL/V6_testes" "[.lento]" --durations yes > testes_lento.txt 2>&1)
    tail -n 3 "$SAIDA/testes.txt" "$SAIDA/testes_lento.txt"
fi
if command -v python3 > /dev/null && python3 -c "import numpy, matplotlib" 2> /dev/null; then
    python3 "$AQUI/../scripts/plot_v6.py" "$SAIDA" --saida "$SAIDA" > "$SAIDA/comparacao.txt" 2>&1
fi
echo "$(date '+%F %T') fim: $SAIDA"
