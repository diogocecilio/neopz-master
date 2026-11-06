#!/bin/bash
# caso.sbatch: os casos completos de um projeto de verificacao (o executavel do main.cpp, com a malha e os tempos do
# artigo) num no de computacao. PROJETO=<V1..V6, ou o nome completo> e CASOS="<caso>[:VAR=valor,...] ...": cada entrada
# roda num processo de 1 CPU, na sua pasta $SAIDA_CASOS/<projeto>/<caso>/ (padrao $DATA/verificacao; sem $DATA,
# resultados_casos nesta pasta), com CASO=<caso> (a menos que as trocas definam CASO) e as trocas VAR=valor; PARALELO ao
# mesmo tempo (padrao: os nucleos fisicos do no, 12 no fermi). O nome do caso e o da pasta: para rodar o mesmo CASO com
# outras variaveis, de outro nome e CASO explicito (retangular_a2:CASO=retangular,AG=2). Um valor pode ser uma lista com
# virgulas: um pedaco sem '=' continua o valor anterior (tabela4:NU_LISTA=0.2,0.25,0.3,KS_INF=0.01).
# Na pasta: log.txt (a saida do programa), caso.txt (projeto, entrada, variaveis, versao do executavel, job, inicio e o
# codigo de saida de cada execucao), variaveis.txt (as variaveis da entrada) e os CSV do programa (resultado.txt por
# ultimo). Retomavel: uma pasta com resultado.txt e pulada; uma pasta rodada com outras variaveis e recusada (de outro
# nome ao caso); REFAZER=1 roda de novo (a pasta antiga vira <caso>.ant_<data>). Uma entrada que outro job esta
# rodando e pulada. Variaveis do ambiente tambem chegam ao programa (por exemplo T_FIM=7200), mas so as da entrada ficam
# em variaveis.txt: prefira pô-las na entrada.
# Submetido por submeter_caso.sh; a mao, de dentro de Projects/Cluster:
#   PROJETO=V3 CASOS="retangular inclinada" sbatch --job-name=verif-caso-V3 caso.sbatch
#SBATCH --job-name=verif-caso
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=24G
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/%x_%j.out
[ -r /etc/profile.d/scratch.sh ] && source /etc/profile.d/scratch.sh
PKG=""
for p in "${VERIF_CLUSTER:-}" "${SLURM_SUBMIT_DIR:-}" "${SLURM_SUBMIT_DIR:-}/.."; do
    if [ -n "$p" ] && [ -f "$p/comum.sh" ] && [ -f "$p/caso.sbatch" ]; then PKG="$(cd "$p" && pwd)"; break; fi
done
[ -n "$PKG" ] || { echo "nao achei Projects/Cluster: submeta de la (ou com submeter_caso.sh)" >&2; exit 1; }
cd "$PKG" || exit 1
source "$PKG/modulos.sh"
source "$PKG/comum.sh"
export MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 LC_ALL=C
PJ="$(projeto_completo "${PROJETO:-}")" || exit 2
[ -n "${CASOS:-}" ] || { echo "defina CASOS=\"<caso>[:VAR=valor,...] ...\"" >&2; exit 2; }
vistos=" "
for e in $CASOS; do
    m="$(valida_caso "$e")"
    [ -z "$m" ] || { echo "CASOS: $m" >&2; exit 2; }
    [[ "$vistos" == *" ${e%%:*} "* ]] && { echo "CASOS: o caso ${e%%:*} aparece duas vezes (cada caso e uma pasta)" >&2; exit 2; }
    vistos="$vistos${e%%:*} "
done
EXE="$PKG/bin/$(alvo_exe "$PJ")"
JOB="${SLURM_JOB_ID:-local$$}"
# o nome verif-caso-<Vn>, tambem num job submetido a mao com o nome padrao
[ -z "${SLURM_JOB_ID:-}" ] || [ "${SLURM_JOB_NAME:-}" = "verif-caso-${PJ%%_*}" ] ||
    scontrol update JobId="$SLURM_JOB_ID" JobName="verif-caso-${PJ%%_*}" > /dev/null 2>&1 || true
SAIDA_CASOS="${SAIDA_CASOS:-${DATA:+$DATA/verificacao}}"
SAIDA_CASOS="${SAIDA_CASOS:-$PKG/resultados_casos}"
D="$SAIDA_CASOS/$PJ"
mkdir -p "$D"
falta="$(aviso_ldd "$EXE")"
if [ -n "$falta" ]; then
    echo "$(hostname -s): $falta" >&2
    echo "$(date '+%F %T') $(hostname -s) job $JOB ($PJ): saiu no inicio: $falta" >> "$PKG/logs/avisos.txt"
    exit 1
fi
P="${PARALELO:-$(paralelo_job)}"
V="$(versao_de "$EXE")"
echo "$(date '+%F %T') job $JOB em $(hostname -s): $PJ ($V), casos: $CASOS; ate $P ao mesmo tempo; saida $D"

# dono_morto <pasta/.ativo>: o job (ou o processo, fora do SLURM) que a trava registra acabou
dono_morto() {
    local j h p out
    read -r j h p < "$1/dono" 2>/dev/null || return 0
    [ "$j" = "$JOB" ] && return 0 # este job, numa vez anterior (o SLURM o recolocou na fila)
    if [[ "$j" =~ ^[0-9]+$ ]]; then
        out="$(squeue -h -j "$j" -o %T 2>&1)" || { case "$out" in *[Ii]"nvalid job id"*) return 0 ;; esac; return 1; }
        case "$out" in '' | COMPLETED* | FAILED* | CANCELLED* | TIMEOUT* | NODE_FAIL* | PREEMPTED* | OUT_OF_MEMORY*) return 0 ;; esac
        return 1
    fi
    [ "$h" = "$(hostname -s)" ] && ! kill -0 "$p" 2>/dev/null
}
# roda <entrada>: uma entrada, na sua pasta
roda() {
    local e="$1" c="${1%%:*}" d="$D/${1%%:*}" rc vars=() t s0 tem_caso=0
    while IFS= read -r t; do
        vars+=("$t")
        [ "${t%%=*}" = CASO ] && tem_caso=1
    done < <(trocas_caso "$e")
    [ $tem_caso = 1 ] || vars=("CASO=$c" ${vars[@]+"${vars[@]}"})
    local novo
    novo="$(printf '%s\n' "${vars[@]}" | sort)"
    if [ -d "$d" ] && [ "${REFAZER:-0}" = 1 ] && [ ! -d "$d/.ativo" ]; then
        mv "$d" "$d.ant_$(date +%Y%m%d_%H%M%S)" && echo "$c: a pasta antiga foi para $(basename "$d").ant_*"
    fi
    if [ -f "$d/variaveis.txt" ] && [ "$(cat "$d/variaveis.txt")" != "$novo" ]; then
        echo "$c: RECUSADO: $d ja foi rodada com outras variaveis ($(tr '\n' ' ' < "$d/variaveis.txt")); de outro nome ao" \
            "caso (ou REFAZER=1)" >&2
        echo "$(date '+%F %T') job $JOB: $PJ/$c recusado: a pasta ja foi rodada com outras variaveis" >> "$PKG/logs/avisos.txt"
        return 0
    fi
    if [ -s "$d/resultado.txt" ]; then echo "$c: ja feito (resultado.txt), pulado"; return 0; fi
    mkdir -p "$d" || return 1
    if ! mkdir "$d/.ativo" 2>/dev/null; then
        if dono_morto "$d/.ativo"; then rm -rf "$d/.ativo"; mkdir "$d/.ativo" 2>/dev/null || return 0
        else echo "$c: outro job esta rodando ($(cat "$d/.ativo/dono" 2>/dev/null)), pulado"; return 0
        fi
    fi
    echo "$JOB $(hostname -s) $BASHPID" > "$d/.ativo/dono"
    echo "$novo" > "$d/variaveis.txt"
    {
        echo "projeto=$PJ entrada=$e"
        echo "variaveis: ${vars[*]}"
        echo "executavel: $EXE ($V)"
        echo "job $JOB em $(hostname -s), inicio $(date '+%F %T')"
    } >> "$d/caso.txt"
    s0=$(date +%s)
    (cd "$d" && env "${vars[@]}" "$EXE" > log.txt 2>&1)
    rc=$?
    echo "$(date '+%F %T') codigo de saida $rc em $(($(date +%s) - s0)) s" >> "$d/caso.txt"
    rm -rf "$d/.ativo"
    if [ -s "$d/resultado.txt" ]; then
        echo "$(date '+%F %T') $c: codigo $rc, $(tr ' ' '\n' < "$d/resultado.txt" | grep -E '^(situacao|seg)=' | tr '\n' ' ')"
    else
        echo "$(date '+%F %T') $c: codigo $rc SEM resultado.txt (veja $d/log.txt)" >&2
    fi
}
for e in $CASOS; do
    while [ "$(jobs -rp | wc -l)" -ge "$P" ]; do wait -n 2>/dev/null || sleep 5; done
    roda "$e" &
done
wait
# o resumo: o estado de cada entrada
ok=0 n=0
for e in $CASOS; do
    d="$D/${e%%:*}"
    n=$((n + 1))
    if [ -s "$d/resultado.txt" ]; then
        ok=$((ok + 1))
        echo "  ${e%%:*}: $(tr ' ' '\n' < "$d/resultado.txt" | grep -E '^(situacao|seg)=' | tr '\n' ' ')"
    else
        echo "  ${e%%:*}: sem resultado.txt"
    fi
done
echo "$(date '+%F %T') job $JOB: fim; $ok de $n casos com resultado.txt em $D"
[ "$ok" -eq "$n" ]
