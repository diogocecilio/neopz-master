#!/bin/bash
# submeter_mc.sh [nos]: o Monte Carlo num cluster SLURM (fermi). Compila, se preciso, no no de login em segundo plano, e
# submete <nos> jobs job_mc.sbatch (padrao 4), um no inteiro cada, que dividem a mesma fila: o caso deterministico de
# cada combinacao (DET=1, o padrao) e depois as realizacoes PRIMEIRA..ULTIMA (padrao 0..999) de todas as combinacoes de
# COMBS (padrao: os 8 casos de Rahardjo alem de f50:36), intercaladas. Volta em segundos: pode desconectar.
# Rodar de novo submete mais nos para a mesma fila (o que esta concluido ou rodando nao se repete), e e tambem o jeito
# de retomar uma campanha parada (apaga $SAIDA/PARAR). Um job que chega perto do limite de tempo (TEMPO, padrao
# 4-00:00:00) submete a sua continuacao. Se o SLURM recusar um job (cota de jobs por usuario), os ja submetidos seguem.
# Os avisos (cota, parametros, bibliotecas, falhas) ficam em $SAIDA/AVISOS, que status_mc.sh mostra.
# Ambiente: COMBS, PRIMEIRA, ULTIMA, DET, SAIDA (padrao $DATA/monte_carlo), TEMPO, PARTICAO, RECOMPILAR=1, FORCAR=1
# (recompila mesmo com jobs na fila), FORCAR_SEM_COMPILAR=1 (submete com o executavel atual, sem recompilar) e as
# variaveis de rodar_mc.sh e do programa (VARIANTE, MALHA, FS_H, COV_C, ...), que vao para os jobs.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
source ../LearningMohrCoulomb2/cluster/modulos.sh
NOS="${1:-4}"
case "$NOS" in '' | *[!0-9]* | 0) echo "uso: $0 [numero de nos, padrao 4]" >&2; exit 2 ;; esac
COMBS="$(echo ${COMBS:-f10:9 f10:80 f10:360 f50:9 f50:80 f100:9 f100:80 f100:3.6} | sed -E 's/,(f10|f50|f100):/ \1:/g')"
export COMBS
# a mesma validacao de rodar_mc.sh (o modo --pastas so le COMBS), antes de compilar e submeter
./rodar_mc.sh --pastas > /dev/null || exit 2
# resultados: $DATA/monte_carlo, a particao de dados de onde o manual do fermi manda rodar os jobs (temporaria, sem
# backup: copie para o PC ou para o $HOME ao terminar); sem $DATA, resultados_mc nesta pasta
SAIDA="${SAIDA:-${DATA:+$DATA/monte_carlo}}"
export PRIMEIRA="${PRIMEIRA:-0}" ULTIMA="${ULTIMA:-999}" DET="${DET:-1}" SAIDA="${SAIDA:-$PWD/resultados_mc}"
[ -n "${TEMPO:-}" ] && export TEMPO
mkdir -p logs "$SAIDA"
EU="${USER:-$(id -un)}"
MAQ="$(hostname -s 2>/dev/null || uname -n)"
MAQ="${MAQ%%.*}"
BOOT="$(cat /proc/sys/kernel/random/boot_id 2>/dev/null || echo -)"
aviso() { echo "$(date '+%F %T') $*" | tee -a "$SAIDA/AVISOS" >&2; }

# a campanha (para status_mc.sh), ja antes da compilacao
{
    echo "COMBS=\"$COMBS\""
    echo "PRIMEIRA=$PRIMEIRA"
    echo "ULTIMA=$ULTIMA"
    echo "DET=$DET"
    echo "VARIANTE=${VARIANTE:-iii}"
} > "$SAIDA/campanha.txt"

# compilacao no no de login (os nos de computacao do fermi nao tem cmake): logs/compilar_mc.pid guarda "pid no boot
# inicio" e logs/compilar_mc.pedido o pedido que a cadeia submete ao fim
compilando_aqui() {
    local pid no boot
    [ -f logs/compilar_mc.pid ] || return 1
    read -r pid no boot _ < logs/compilar_mc.pid || return 1
    [ "$no" = "$MAQ" ] && [ "$boot" = "$BOOT" ] && kill -0 "$pid" 2>/dev/null &&
        tr '\0' ' ' 2>/dev/null < "/proc/$pid/cmdline" | grep -q 'compilar.sh'
}
# algum fonte do executavel mais novo que bin/SlopeRandomField (por exemplo depois de um git pull): a biblioteca (tudo
# fora de Projects), este projeto e LearningMohrCoulomb2, cujo main.cpp ele inclui
fonte_nova() {
    find ../.. \( -name .git -o -name build -o -name 'build_*' -o -name deps -o -name bin -o -name 'resultados*' \
        -o -name .claude -o -name logs \) -prune -o \
        \( -path '../../Projects/*' ! -path '../../Projects/SlopeRandomField*' ! -path '../../Projects/LearningMohrCoulomb2*' \
        ! -path '../../Projects/CMakeLists.txt' \) -prune -o -newer bin/SlopeRandomField \
        \( -name '*.cpp' -o -name '*.h' -o -name '*.hpp' -o -name CMakeLists.txt -o -name '*.cmake' \) -print -quit 2>/dev/null || true
}
motivo=""
if ! [ -x bin/SlopeRandomField ]; then motivo="bin/SlopeRandomField ainda nao existe"
elif [ -n "${RECOMPILAR:-}" ]; then motivo="RECOMPILAR=1"
elif [ -z "${FORCAR_SEM_COMPILAR:-}" ]; then
    nova="$(fonte_nova)"
    [ -n "$nova" ] && motivo="fontes mais novos que bin/SlopeRandomField (por exemplo $nova)"
fi
if [ -n "$motivo" ]; then
    if [ -z "${CADEIA_MC:-}" ] && compilando_aqui; then
        echo "a compilacao ja esta rodando neste no (pid $(cut -d' ' -f1 logs/compilar_mc.pid)). Ao fim ela submete so o" >&2
        echo "pedido original ($(cat logs/compilar_mc.pedido 2>/dev/null)); repita este pedido depois dela." >&2
        exit 1
    fi
    # recompilar troca build/libpz.so, que os jobs rodando usam (o executavel em bin/ carrega a de build/): so sem jobs
    if ! fila="$(squeue -h -u "$EU" -n slopeRF -o %i 2>&1)"; then
        if [ -z "${FORCAR:-}" ]; then
            echo "$motivo, mas o squeue nao respondeu ($fila): sem saber se ha jobs usando a libpz.so, nao recompila." >&2
            echo "Tente de novo daqui a pouco (ou FORCAR=1)." >&2
            exit 1
        fi
        fila=""
    fi
    if [ -n "$fila" ] && [ -z "${FORCAR:-}" ]; then
        echo "$motivo, mas ha $(echo "$fila" | wc -l) job(s) slopeRF na fila: recompilar agora trocaria a libpz.so que eles usam." >&2
        echo "Espere terminarem (ou scancel -n slopeRF; as realizacoes continuam dos checkpoints) e rode de novo." >&2
        echo "Para so acrescentar nos com o executavel atual: FORCAR_SEM_COMPILAR=1 ./submeter_mc.sh <nos>" >&2
        exit 1
    fi
    # Catch2 (o CMake o baixaria ao configurar): o do pacote cluster, baixado aqui se faltar
    C2=../LearningMohrCoulomb2/cluster/deps/Catch2
    if ! [ -f "$C2/CMakeLists.txt" ]; then
        rm -rf "$C2"
        git -c advice.detachedHead=false clone -q --depth 1 --branch v2.13.10 https://github.com/catchorg/Catch2.git "$C2"
    fi
    # nohup + setsid: continua depois que voce desconectar. Ao fim, o executavel fica com a data do inicio da compilacao
    # (um fonte trocado durante ela, por um git pull, continua mais novo e e recompilado) e os jobs sao submetidos
    echo "$NOS no(s), COMBS=\"$COMBS\", realizacoes $PRIMEIRA..$ULTIMA, DET=$DET, SAIDA=$SAIDA" > logs/compilar_mc.pedido
    nohup setsid bash -c 'nice -n 10 ./compilar.sh && touch -r logs/compilar_mc.pid bin/SlopeRandomField &&
        CADEIA_MC=1 RECOMPILAR= ./submeter_mc.sh "$1"' _ "$NOS" > logs/compilar_mc.log 2>&1 < /dev/null &
    echo "$! $MAQ $BOOT $(date '+%F %T')" > logs/compilar_mc.pid
    echo "$motivo: compilando neste no de login em segundo plano (pid $!, log em logs/compilar_mc.log; a biblioteca"
    echo "inteira leva uns 15 a 30 min, so o programa 1 a 2 min). Ao terminar, o(s) $NOS job(s) e submetido(s) sozinho(s)."
    echo "Pode desconectar; acompanhe com ./status_mc.sh. Para desistir: kill -- -$!"
    exit 0
fi

# PARAR: submeter e um pedido de recomeco (apaga PARAR); so a cadeia da compilacao o respeita (o usuario parou depois)
if [ -e "$SAIDA/PARAR" ]; then
    if [ -n "${CADEIA_MC:-}" ]; then
        aviso "compilacao terminada, mas $SAIDA/PARAR existe: nenhum job submetido (apague-o e rode ./submeter_mc.sh $NOS)"
        exit 1
    fi
    rm -f "$SAIDA/PARAR"
    echo "apagado $SAIDA/PARAR: a campanha volta a rodar"
fi
# os parametros de cada combinacao conferem com os de parametros.txt (senao os jobs sairiam logo no inicio)
if ! verif="$(./rodar_mc.sh --verifica 0 0 2>&1 | tr '\n' ' '; exit "${PIPESTATUS[0]}")"; then
    aviso "nao submetido: $verif"
    exit 2
fi
# travas de parada de jobs que ja acabaram
for f in "$SAIDA"/.parar.*; do
    [ -e "$f" ] || continue
    j="${f##*.}"
    if out="$(squeue -h -j "$j" -o %i 2>&1)"; then [ -n "$out" ] || rm -f "$f"
    else case "$out" in *[Ii]"nvalid job id"*) rm -f "$f" ;; esac
    fi
done
ja="$(squeue -h -u "$EU" -n slopeRF -o %i 2>/dev/null | wc -l)" || ja=0
[ "$ja" -gt 0 ] && echo "ja ha $ja job(s) slopeRF na fila"
P=()
[ -n "${PARTICAO:-}" ] && P=(--partition="$PARTICAO")
ok=0
err="$(mktemp)"
for ((i = 1; i <= NOS; i++)); do
    if out="$(sbatch --parsable ${P[@]+"${P[@]}"} ${TEMPO:+--time="$TEMPO"} job_mc.sbatch 2> "$err")"; then
        j="${out##*$'\n'}"
        j="${j%%;*}"
        ok=$((ok + 1))
        echo "job $j submetido"
        echo "$(date '+%F %T') $j COMBS=\"$COMBS\" realizacoes $PRIMEIRA..$ULTIMA DET=$DET SAIDA=$SAIDA" >> logs/jobs_mc.txt
    else
        aviso "o SLURM recusou o job $i de $NOS: $(cat "$err") $out"
        break
    fi
done
rm -f "$err"
[ "$ok" -gt 0 ] && rm -f "$SAIDA/CONTINUACAO_PENDENTE"
echo "$ok de $NOS jobs submetidos: $(echo $COMBS | wc -w) combinacoes, realizacoes $PRIMEIRA..$ULTIMA$([ "$DET" = 1 ] &&
    echo ' e os deterministicos'), resultados em $SAIDA"
echo "acompanhe com: ./status_mc.sh"
[ "$ok" -gt 0 ]
