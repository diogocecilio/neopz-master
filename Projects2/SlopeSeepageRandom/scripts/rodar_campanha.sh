#!/bin/bash
# Compila o SlopeSeepageRandom e inicia (ou retoma) a campanha de Monte Carlo do artigo (Vargas Ceron et al. 2025).
# Os caminhos são descobertos a partir da localização deste script; nada precisa ser editado.
#
# Uso:  bash Projects2/SlopeSeepageRandom/scripts/rodar_campanha.sh [alvo] [processos] [diretório da campanha]
#       bash Projects2/SlopeSeepageRandom/scripts/rodar_campanha.sh parar [diretório da campanha]
#   alvo:      amostras por caso (padrão 1000), "artigo" (S do artigo) ou "cov5" (CoV(Pf) < 5 %)
#   processos: padrão = número de núcleos físicos
#   diretório: padrão ~/campanha_artigo
# Variável opcional: BUILD=<diretório de compilação> (padrão <neopz>/build-campanha)
#
# Rodar de novo com outro alvo continua a mesma campanha (nada do que já foi calculado é refeito).
set -eu
SCRIPTS=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
if [ "${1:-}" = parar ]; then
    CAMP=${2:-$HOME/campanha_artigo}
    if [ -s "$CAMP/fila.pid" ] && kill -0 "$(cat "$CAMP/fila.pid")" 2>/dev/null; then
        kill -- -"$(cat "$CAMP/fila.pid")" && echo "campanha em $CAMP parada (para continuar, rode o script de novo)"
    else
        echo "nenhuma campanha rodando em $CAMP"
    fi
    exit 0
fi
NEOPZ=$(cd "$SCRIPTS/../../.." && pwd)
FISICOS=$(lscpu -p=Core,Socket 2>/dev/null | grep -v '^#' | sort -u | wc -l)
[ "${FISICOS:-0}" -gt 0 ] || FISICOS=$(nproc)
ALVO=${1:-1000}
P=${2:-$FISICOS}
CAMP=${3:-$HOME/campanha_artigo}
BUILD=${BUILD:-$NEOPZ/build-campanha}
EXE=$BUILD/Projects2/SlopeSeepageRandom/SlopeSeepageRandom

echo "NeoPZ:     $NEOPZ ($(git -C "$NEOPZ" rev-parse --abbrev-ref HEAD 2>/dev/null) $(git -C "$NEOPZ" rev-parse --short HEAD 2>/dev/null))"
echo "build:     $BUILD"
echo "campanha:  $CAMP (alvo $ALVO, $P processos)"

# dependências
falta=""
for c in cmake g++ make python3 lscpu; do command -v $c > /dev/null || falta="$falta $c"; done
ldconfig -p | grep -q 'liblapack\.so ' || falta="$falta liblapack-dev"
ldconfig -p | grep -q 'libblas\.so ' || falta="$falta libblas-dev"
if [ -n "$falta" ]; then
    echo "faltam:$falta"
    echo "instale com: sudo apt install cmake g++ make python3 util-linux liblapack-dev libblas-dev"
    exit 1
fi
if [ -s "$CAMP/fila.pid" ] && kill -0 "$(cat "$CAMP/fila.pid")" 2>/dev/null; then
    echo "a campanha em $CAMP já está rodando; acompanhe com:"
    echo "  PROCESSOS=$P python3 $SCRIPTS/campanha_artigo.py status $CAMP"
    exit 0
fi

# compilação (incremental): Release, plasticidade e LAPACK (autoproblema da KL)
GEN="Unix Makefiles"
command -v ninja > /dev/null && GEN=Ninja
if [ ! -f "$BUILD/CMakeCache.txt" ]; then
    cmake -S "$NEOPZ" -B "$BUILD" -G "$GEN" -DCMAKE_BUILD_TYPE=Release -DBUILD_PLASTICITY_MATERIALS=ON -DUSING_LAPACK=ON
fi
cmake --build "$BUILD" --target SlopeSeepageRandom -j "$(nproc)"
if ! grep -q y_min_banda "$EXE"; then
    echo "o executável não é da versão do branch claude/great-clarke-xist30 (falta a saída .modo); atualize o NeoPZ"
    exit 1
fi

# fila e execução em segundo plano
mkdir -p "$CAMP"
python3 "$SCRIPTS/campanha_artigo.py" jobs "$EXE" "$CAMP" --alvo "$ALVO" > "$CAMP/jobs.txt"
echo "$(grep -c . "$CAMP/jobs.txt") comandos na fila"
# sessão própria: sobrevive ao fechamento do terminal e é parada inteira com "parar"
setsid bash -c 'echo $$ > "$1/fila.pid"; exec bash "$2/fila.sh" "$1/jobs.txt" "$3"' _ "$CAMP" "$SCRIPTS" "$P" \
    > "$CAMP/fila.log" 2>&1 < /dev/null &
echo
echo "rodando em segundo plano (pode fechar o terminal). Comandos úteis:"
echo "  andamento e previsão:  PROCESSOS=$P python3 $SCRIPTS/campanha_artigo.py status $CAMP"
echo "  parar:                 bash $0 parar $CAMP"
echo "  continuar/aumentar:    bash $0 <alvo> $P $CAMP"
echo "  pacote para enviar:    python3 $SCRIPTS/campanha_artigo.py pacote $CAMP ~/resultados_campanha.tar.gz"
