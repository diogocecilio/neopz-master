#!/bin/bash
# Contraexemplos A1 (dilatância, executável SlopeRainfallPV) e A2 (fase de ar) em paralelo, um caso por núcleo.
#
# Rodar na pasta dos resultados (p. ex. ~/projects/runs):
#   <neopz>/Projects/SlopeRainfall/scripts/rodar_A1_A2.sh rodar <pasta do build> [núcleos A1 = 4] [núcleos A2 = 4]
#   <neopz>/Projects/SlopeRainfall/scripts/rodar_A1_A2.sh status
#   <neopz>/Projects/SlopeRainfall/scripts/rodar_A1_A2.sh zip
#
# rodar   compila SlopeRainfall e SlopeRainfallPV e lança os dois lotes em segundo plano (nohup: continuam com o
#         terminal fechado; logs em lote_A1.log e lote_A2.log). Um lote que já está rodando não é lançado de novo e
#         os casos com summary.txt são pulados: depois de uma interrupção basta repetir o comando.
# status  casos prontos (A1: 15, A2: 22) e casos em execução.
# zip     a1_a2.zip com os resultados para enviar (sem snaps_pa.csv; snaps_pw.csv é usado nas figuras do A2).

SCR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

case "${1:-status}" in
rodar)
    [ -n "$2" ] || { echo "uso: $0 rodar <pasta do build> [núcleos A1] [núcleos A2]"; exit 1; }
    BUILD=$(cd "$2" && pwd) || exit 1
    J1=${3:-4}
    J2=${4:-4}
    EXE="$BUILD/Projects/SlopeRainfall/SlopeRainfall"
    echo "== compilando SlopeRainfall e SlopeRainfallPV em $BUILD"
    cmake --build "$BUILD" --target SlopeRainfall -j 8 || { echo "falha ao compilar SlopeRainfall"; exit 1; }
    cmake --build "$BUILD" --target SlopeRainfallPV -j 8 || { echo "falha ao compilar SlopeRainfallPV"; exit 1; }
    if pgrep -f "[c]ontraexemplos.py rodar --itens A1" > /dev/null; then
        echo "A1 já está rodando: não foi lançado de novo"
    else
        nohup python3 "$SCR/contraexemplos.py" rodar --itens A1 --jobs "$J1" --exe "$EXE" > lote_A1.log 2>&1 &
        echo "A1 lançado com $J1 núcleos; log em $(pwd)/lote_A1.log"
    fi
    if pgrep -f "[f]ase_ar.py rodar" > /dev/null; then
        echo "A2 já está rodando: não foi lançado de novo"
    else
        nohup python3 "$SCR/fase_ar.py" rodar --jobs "$J2" --exe "$EXE" > lote_A2.log 2>&1 &
        echo "A2 lançado com $J2 núcleos; log em $(pwd)/lote_A2.log"
    fi
    echo "acompanhe com: $0 status"
    ;;
status)
    n1=$(ls resultados/contraexemplos/A1_*/summary.txt 2> /dev/null | wc -l)
    n2=$(ls resultados/fase_ar/*/summary.txt 2> /dev/null | wc -l)
    echo "prontos: A1 $n1 de 15, A2 $n2 de 22"
    em=$(pgrep -fa "[S]lopeRainfall(PV)? solo" | grep -o "out=[^ ]*" | sed 's/^out=/  /')
    echo "em execução:"
    echo "${em:-  nenhum caso}"
    ;;
zip)
    rm -f a1_a2.zip
    zip -qr a1_a2.zip resultados/contraexemplos/A1_* resultados/fase_ar -x "*/snaps_pa.csv"
    ls -lh a1_a2.zip
    ;;
*)
    echo "uso: $0 rodar <pasta do build> [núcleos A1] [núcleos A2] | status | zip"
    exit 1
    ;;
esac
