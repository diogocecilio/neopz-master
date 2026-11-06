#!/bin/bash
# compilar.sh [projeto ...]: configura e compila, no no de login (o unico do fermi com cmake), a biblioteca NeoPZ e os
# alvos dos projetos de verificacao, o executavel e os testes de cada um (padrao: todos os de Projects/CMakeLists.txt,
# V1 a V6), em ./build, e copia os executaveis para ./bin. A troca e atomica: uma realizacao que ja esta rodando continua
# com o executavel antigo. bin/versao.txt tem uma linha por executavel, "<executavel> <commit>[+mod] <data>" (+mod: havia
# fontes alterados ou novos fora do git), que o Monte Carlo grava no resumo.csv de cada realizacao. Os testes ficam em
# build/Projects/<projeto>/ (o ctest os acha la; ver testes.sbatch).
# Um projeto que nao compila nao impede os outros: ele fica de fora, a mensagem do compilador fica no log, e o codigo de
# saida e 1. Se a biblioteca mudou nesta compilacao, o executavel antigo dele sai de bin/ (vai para bin/antigos/): ele
# carregaria a libpz.so nova, de outra versao. Pelo mesmo motivo, com so alguns projetos pedidos e a biblioteca nova, os
# outros que estao em bin/ sao compilados tambem.
# Os executaveis em bin/ carregam a libpz.so de ./build: recompilar com jobs destes projetos na fila (verif-*) ou com
# processos rodando um executavel de bin/ ou de build/ nesta maquina os deixaria com a biblioteca nova e o executavel
# velho, entao compilar.sh recusa (FORCAR=1 compila mesmo assim). Duas compilacoes ao mesmo tempo tambem sao recusadas.
# Ambiente: USAR_MKL=ON|OFF (padrao ON: PARDISO com ordenacao METIS; OFF: skyline, mais lento), MKL_THREAD_MODEL (padrao
# seq), JOBS (processos do make, padrao 4), BUILD (padrao ./build), FORCAR=1, BAIXAR_CATCH2=0 (nao baixa o Catch2 para
# deps/; o CMake entao o baixa sozinho ao configurar). Com os modulos carregados (source modulos.sh; no fermi, o MKL).
set -euo pipefail
PKG="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=comum.sh
source "$PKG/comum.sh"
NEOPZ="$VERIF_NEOPZ"
BUILD="${BUILD:-$PKG/build}"
JOBS="${JOBS:-4}"
USAR_MKL="${USAR_MKL:-ON}"
mkdir -p "$PKG/logs" "$PKG/bin"

PROJS=()
if [ $# -gt 0 ]; then
    for a in "$@"; do
        p="$(projeto_completo "$a")" || exit 2
        PROJS+=("$p")
    done
else
    mapfile -t PROJS < <(projetos_verif)
fi
[ ${#PROJS[@]} -gt 0 ] || { echo "nenhum projeto em $NEOPZ/Projects/CMakeLists.txt" >&2; exit 2; }
mkdir -p "$BUILD"
BUILD="$(cd "$BUILD" && pwd)"

# uma compilacao por vez (a em segundo plano dos submeter_*.sh tambem passa por aqui)
if [ -z "${VERIF_EM_SEGUNDO_PLANO:-}" ] && compilando_aqui; then
    echo "a compilacao em segundo plano esta rodando neste no (pid $(cut -d' ' -f1 "$PKG/logs/compilar.pid"), log em" \
        "logs/compilar.log): espere ela terminar" >&2
    exit 1
fi
exec 8>> "$PKG/logs/compilar.trava"
if command -v flock > /dev/null 2>&1; then
    flock -n -E 99 8 || {
        rc=$?
        if [ $rc -eq 99 ]; then echo "outra compilacao esta rodando (logs/compilar.trava): espere ela terminar" >&2; exit 1; fi
        echo "aviso: flock falhou (codigo $rc); seguindo sem a trava" >&2
    }
fi

# o executavel em bin/ carrega a libpz.so de $BUILD: nada destes projetos rodando (aqui ou em jobs)
if [ -z "${FORCAR:-}" ]; then
    procs="$(processos_verif "$BUILD")"
    if [ -n "$procs" ]; then
        echo "ha processos rodando executaveis destes projetos nesta maquina (eles usam a libpz.so de $BUILD):" >&2
        echo "$procs" | head -n 5 | sed 's/^/  /' >&2
        echo "espere terminarem (ou FORCAR=1)" >&2
        exit 1
    fi
    if ! fila="$(jobs_verif)"; then
        echo "o squeue nao respondeu: sem saber se ha jobs usando a libpz.so de $BUILD, nao compila (ou FORCAR=1)" >&2
        exit 1
    fi
    if [ -n "$fila" ]; then
        echo "ha $(echo "$fila" | wc -l) job(s) destes projetos no SLURM usando a libpz.so de $BUILD: espere terminarem" \
            "(ou scancel; ou FORCAR=1):" >&2
        echo "$fila" | sed 's/^/  /' >&2
        exit 1
    fi
fi

# Catch2 (os testes): uma copia em deps/, para nao depender do GitHub a cada configuracao do zero
C2="$PKG/deps/Catch2"
if ! [ -f "$C2/CMakeLists.txt" ] && [ "${BAIXAR_CATCH2:-1}" = 1 ] && command -v git > /dev/null 2>&1; then
    rm -rf "$C2"
    mkdir -p "$PKG/deps"
    if ! git -c advice.detachedHead=false clone -q --depth 1 --branch v2.13.10 https://github.com/catchorg/Catch2.git \
        "$C2" 2>&1; then
        echo "aviso: nao consegui baixar o Catch2 para deps/ (o CMake tenta baixar sozinho)" >&2
        rm -rf "$C2"
    fi
fi

ARGS=(-S "$NEOPZ" -B "$BUILD" -DCMAKE_BUILD_TYPE=Release -DBUILD_PLASTICITY_MATERIALS=ON -DBUILD_PROJECTS=ON
      -DUSING_MKL="$USAR_MKL")
if [ "$USAR_MKL" = ON ]; then ARGS+=(-DMKL_THREAD_MODEL="${MKL_THREAD_MODEL:-seq}"); else ARGS+=(-DUSING_LAPACK=OFF); fi
[ -f "$C2/CMakeLists.txt" ] && ARGS+=(-DFETCHCONTENT_SOURCE_DIR_CATCH2="$C2")
[ -f "$BUILD/CMakeCache.txt" ] || ARGS+=(-G "Unix Makefiles")
echo "$(date '+%F %T') configurando $BUILD (projetos: ${PROJS[*]})"
cmake "${ARGS[@]}"

lib_antes="$(stat -c '%Y %i' "$BUILD/libpz.so" 2>/dev/null || echo -)"
# todos os alvos de uma vez (o make paraleliza entre eles), seguindo depois de um erro (-k), para que um projeto com
# erro nao impeca os outros; CMake < 3.15 nao aceita varios alvos: um por vez
case "$(sed -n 's/^CMAKE_GENERATOR:INTERNAL=//p' "$BUILD/CMakeCache.txt")" in
Ninja*) KG=(-k 0) ;;
*Makefiles*) KG=(-k) ;;
*) KG=() ;;
esac
read -r cmaj cmin < <(cmake --version | sed -n '1s/[^0-9]*\([0-9][0-9]*\)\.\([0-9][0-9]*\).*/\1 \2/p') || true
OK=() FALHA=()
# compila <projeto...>: compila os alvos dos projetos e poe cada um em OK ou FALHA
compila() {
    local p e t a rc=0 alvos=()
    for p in "$@"; do alvos+=("$(alvo_exe "$p")" "$(alvo_teste "$p")"); done
    echo "$(date '+%F %T') compilando ${alvos[*]} com $JOBS processos"
    if [ "${cmaj:-0}" -gt 3 ] || { [ "${cmaj:-0}" -eq 3 ] && [ "${cmin:-0}" -ge 15 ]; }; then
        cmake --build "$BUILD" --target "${alvos[@]}" -j "$JOBS" -- ${KG[@]+"${KG[@]}"} || rc=$?
    else
        for a in "${alvos[@]}"; do cmake --build "$BUILD" --target "$a" -j "$JOBS" || rc=$?; done
    fi
    # quais projetos compilaram (com erro em algum, cada um de novo: rapido se ja esta em dia)
    for p in "$@"; do
        e="$(alvo_exe "$p")" t="$(alvo_teste "$p")"
        if [ $rc -ne 0 ]; then
            echo "$(date '+%F %T') conferindo $p"
            if ! cmake --build "$BUILD" --target "$e" -j "$JOBS" || ! cmake --build "$BUILD" --target "$t" -j "$JOBS"; then
                FALHA+=("$p")
                continue
            fi
        fi
        if [ -x "$BUILD/Projects/$p/$e" ] && [ -x "$BUILD/Projects/$p/$t" ]; then OK+=("$p"); else FALHA+=("$p"); fi
    done
}
compila "${PROJS[@]}"
lib_depois="$(stat -c '%Y %i' "$BUILD/libpz.so" 2>/dev/null || echo -)"
# so alguns projetos pedidos e a biblioteca mudou: os outros que estao em bin/ carregariam a libpz.so nova com o
# executavel velho; compilados tambem
if [ "$lib_antes" != "$lib_depois" ]; then
    RESTO=()
    for p in $(projetos_verif); do
        [[ " ${PROJS[*]} " == *" $p "* ]] && continue
        [ -e "$PKG/bin/$(alvo_exe "$p")" ] && RESTO+=("$p")
    done
    if [ ${#RESTO[@]} -gt 0 ]; then
        echo "$(date '+%F %T') a biblioteca mudou: compilando tambem ${RESTO[*]} (os executaveis em bin/ a carregam)"
        compila "${RESTO[@]}"
    fi
fi

# a versao: o commit do git, com "+mod" se havia fontes alterados ou novos (fora do git), e a data
v="$(git -C "$NEOPZ" rev-parse --short HEAD 2>/dev/null || echo sem_git)"
if [ "$v" != sem_git ] && [ -n "$(git -C "$NEOPZ" status --porcelain --untracked-files=all -- '*.cpp' '*.h' '*.hpp' \
    '*CMakeLists.txt' '*.cmake' ':(exclude)*_deps/*' ':(exclude)Projects/Cluster/*' 2>/dev/null)" ]; then
    v="$v+mod"
fi
versao="$v $(date '+%F_%T')"
# versao_txt <executavel> <versao ou vazio>: troca (ou tira) a linha do executavel em bin/versao.txt
versao_txt() {
    local f="$PKG/bin/versao.txt"
    {
        if [ -f "$f" ]; then awk -v e="$1" '$1 != e' "$f"; fi
        if [ -n "$2" ]; then echo "$1 $2"; fi
    } | sort > "$f.tmp.$$"
    mv -f "$f.tmp.$$" "$f"
}
COPIADOS=()
for p in ${OK[@]+"${OK[@]}"}; do
    e="$(alvo_exe "$p")"
    # troca atomica: realizacoes rodando continuam com o executavel antigo
    cp "$BUILD/Projects/$p/$e" "$PKG/bin/.$e.novo"
    mv -f "$PKG/bin/.$e.novo" "$PKG/bin/$e"
    versao_txt "$e" "$versao"
    COPIADOS+=("$e")
done
for p in ${FALHA[@]+"${FALHA[@]}"}; do
    e="$(alvo_exe "$p")"
    if [ -e "$PKG/bin/$e" ] && [ "$lib_antes" != "$lib_depois" ]; then
        mkdir -p "$PKG/bin/antigos"
        mv -f "$PKG/bin/$e" "$PKG/bin/antigos/$e.$(date +%Y%m%d_%H%M%S)"
        versao_txt "$e" ""
        echo "aviso: $p nao compilou e a libpz.so mudou: bin/$e foi para bin/antigos/ (rodaria com outra biblioteca)" >&2
    elif [ -e "$PKG/bin/$e" ]; then
        echo "aviso: $p nao compilou; bin/$e continua o anterior ($(versao_de "$PKG/bin/$e"); a biblioteca nao mudou)" >&2
    fi
done
# onde esta o ctest (os nos de computacao do fermi nao tem cmake; testes.sbatch tenta este caminho)
command -v ctest > "$BUILD/ctest_caminho" 2>/dev/null || rm -f "$BUILD/ctest_caminho"
[ ${#COPIADOS[@]} -gt 0 ] && echo "Compilado: ${COPIADOS[*]} ($versao; em $PKG/bin; testes em $BUILD/Projects/<projeto>/)"
if [ ${#FALHA[@]} -gt 0 ]; then
    echo "Falhou: ${FALHA[*]} (veja os erros do compilador acima)"
    exit 1
fi
