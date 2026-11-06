#!/bin/bash
# v5.sbatch: as execucoes do V5 (Abbasov et al. 2024: encosta em duas camadas, 10 anos a 600 mm/ano e 15 h a 20 mm/h,
# LFS) no SLURM. Cada execucao roda num processo, na sua pasta ($SAIDA/<nome>), ao mesmo tempo (um nucleo cada; ~25 a
# 40 min cada na malha malhas/encosta.msh num nucleo de 2,1 GHz, ~2 h na encosta_muito_fina.msh); no fim junta os
# resumo.csv em $SAIDA/resumo.csv, roda a comparacao com o artigo (scripts/compara_v5.py, se houver python3 com numpy e
# matplotlib) e os testes [.lento] se TESTES=1.
# Uso, de dentro desta pasta, depois de compilar (compilar_v5.sh):
#   sbatch v5.sbatch                                  # todas as execucoes da lista abaixo
#   EXECUCOES="tabela1 tabela4" TESTES=1 sbatch v5.sbatch
# Ambiente: BUILD (o mesmo de compilar_v5.sh), SAIDA (padrao $DATA/v5 ou ./resultados_v5), EXECUCOES (nomes da lista
# abaixo), TESTES. Cada nome define as variaveis do executavel (CASO, MALHA, KS_INF, NU, ...; ver main.cpp e o README).
#SBATCH --job-name=v5_encosta
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=v5_%j.out
[ -r /etc/profile.d/scratch.sh ] && source /etc/profile.d/scratch.sh
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}" || exit 1
AQUI="$PWD"
NEOPZ="$(cd "$AQUI/../../.." && pwd)"
source "$NEOPZ/Projects/Cluster/modulos.sh"
export MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
BUILD="${BUILD:-$NEOPZ/build_verif}"
EXE="$BUILD/Projects/V5_EncostaDuasCamadas/V5_EncostaDuasCamadas"
SAIDA="${SAIDA:-${DATA:+$DATA/v5}}"
SAIDA="${SAIDA:-$AQUI/resultados_v5}"
SENS="NU_LISTA=0.2,0.25,0.3,0.35,0.4 C_LISTA=3,3.5,4,4.5,5,5.5,6 PHI_LISTA=25,30,35,40,45"
# nome -> variaveis do executavel
declare -A CONF=(
  [tabela1]="CASO=tabela1 VTK=1 NU_LISTA=0.2,0.25,0.3,0.35,0.4"
  [tabela4]="CASO=tabela4 $SENS"
  [fig9]="CASO=fig9"
  [tabela4_ks001]="CASO=tabela4 KS_INF=0.01"
  [tabela4_ks003]="CASO=tabela4 KS_INF=0.03"
  [tabela4_ks007]="CASO=tabela4 KS_INF=0.07"
  [tabela4_ks01]="CASO=tabela4 KS_INF=0.1"
  [tabela4_muito_fina]="CASO=tabela4 MALHA=encosta_muito_fina.msh"
  [tabela4_dt]="CASO=tabela4 DT2_H=0.05 DT_MAX_D=10"
  [tabela4_plano_tensao]="CASO=tabela4 NU=0.2307692"
)
EXECUCOES="${EXECUCOES:-tabela1 tabela4 fig9 tabela4_ks001 tabela4_ks003 tabela4_ks007 tabela4_ks01 tabela4_muito_fina tabela4_dt tabela4_plano_tensao}"
if ! [ -x "$EXE" ] || ldd "$EXE" | grep -q "not found"; then
    echo "$EXE ausente ou com bibliotecas faltando: rode compilar_v5.sh antes" >&2
    exit 1
fi
mkdir -p "$SAIDA"
echo "$(date '+%F %T') $(hostname -s): execucoes $EXECUCOES, saida $SAIDA"
# no maximo SLURM_CPUS_PER_TASK processos ao mesmo tempo; o log sem os '*' da montagem serial da biblioteca
MAXP="${SLURM_CPUS_PER_TASK:-4}"
for e in $EXECUCOES; do
    [ -n "${CONF[$e]:-}" ] || { echo "execucao desconhecida: $e" >&2; continue; }
    while [ "$(jobs -rp | wc -l)" -ge "$MAXP" ]; do sleep 5; done
    mkdir -p "$SAIDA/$e"
    (cd "$SAIDA/$e" && env ${CONF[$e]} VERBOSO=1 "$EXE" 2>&1 | tr -d '*' | grep -v '^$' > log.txt;
     echo "$e: codigo ${PIPESTATUS[0]}" >> log.txt) &
done
wait
# junta os resumo.csv (o cabecalho da primeira execucao com as mesmas colunas) com o nome de cada uma
primeiro=1
for e in $EXECUCOES; do
    [ -f "$SAIDA/$e/resumo.csv" ] || { echo "$e: sem resumo.csv" >&2; continue; }
    if [ "$primeiro" = 1 ]; then echo "execucao,$(head -n 1 "$SAIDA/$e/resumo.csv")" > "$SAIDA/resumo.csv"; primeiro=0; fi
    echo "$e,$(tail -n 1 "$SAIDA/$e/resumo.csv")" >> "$SAIDA/resumo.csv"
    tail -n 1 "$SAIDA/$e/log.txt"
done
if command -v python3 > /dev/null && python3 -c "import numpy, matplotlib" 2> /dev/null &&
   [ -f "$SAIDA/tabela1/perfis.csv" ] && [ -f "$SAIDA/tabela4/metricas.csv" ]; then
    python3 "$AQUI/../scripts/compara_v5.py" --tab1 "$SAIDA/tabela1" --tab4 "$SAIDA/tabela4" --figuras "$SAIDA" \
        --csv "$SAIDA/comparacao.csv" > "$SAIDA/comparacao.txt"
fi
if [ "${TESTES:-0}" = 1 ]; then
    # sem filtro: os testes padrao (um filtro "[V5]" incluiria tambem o oculto [.lento])
    (cd "$SAIDA" && "$BUILD/Projects/V5_EncostaDuasCamadas/V5_testes" -d yes 2>&1 | tr -d '*' > testes.txt;
     "$BUILD/Projects/V5_EncostaDuasCamadas/V5_testes" "[.lento]" -d yes 2>&1 | tr -d '*' > testes_lento.txt)
    tail -n 3 "$SAIDA/testes.txt" "$SAIDA/testes_lento.txt"
fi
echo "$(date '+%F %T') fim: $SAIDA"
