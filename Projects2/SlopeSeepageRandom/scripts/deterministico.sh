#!/bin/bash
# Casos determinísticos do artigo (Vargas Ceron et al. 2025) com o SlopeSeepageRandom.
# Uso: deterministico.sh <executável> <diretório de saída> [parte]   (parte: base | beta | alfa | hw | todas)
set -u
EXE=$(readlink -f "$1"); OUT=$2; PARTE=${3:-todas}
mkdir -p "$OUT"; cd "$OUT"
run() {  # nome, argumentos...
    local nome=$1; shift
    [ -s "$nome.log" ] && grep -q "tempo total" "$nome.log" && return
    "$EXE" det "$@" > "$nome.log" 2>&1
}
if [ "$PARTE" = base ] || [ "$PARTE" = todas ]; then
    # Cho (2010): coesivo 2:1 (FS = 1.356) e c-phi 1:1 com H = 10 m (FS = 1.204); artigo: Gamma = 1.354 e 1.777
    run cho_coesivo_mc  caso=cho_coesivo modelo=mc h=1 adapt=3
    run cho_coesivo_mc_h2 caso=cho_coesivo modelo=mc h=2 adapt=3   # malha do Monte Carlo (nível 2)
    run cho_cphi_mc     caso=cho_cphi    modelo=mc h=2 adapt=3
    run cho_cphi_mcc    caso=cho_cphi    modelo=mcc h=2 adapt=3
    # talude com rebaixamento rápido (Tabela 2): Gamma = 1.336
    run percolacao_mc   caso=percolacao  modelo=mc h=1 adapt=3
    run percolacao_mcc  caso=percolacao  modelo=mcc h=1 adapt=3
fi
if [ "$PARTE" = beta ] || [ "$PARTE" = todas ]; then
    # Fig. 9: Gamma x beta para alpha = 1, 5, 10 (H = hw = 5 m, c = 10, phi = 30)
    for a in 1 5 10; do for b in 15 30 45 60 75 90; do
        run beta${b}_alfa${a} caso=percolacao modelo=mc h=1 adapt=2 fs=0 beta=$b alpha=$a
    done; done
fi
if [ "$PARTE" = alfa ] || [ "$PARTE" = todas ]; then
    # Tabela 5 (coluna determinística): alpha = 1..5
    for a in 1 2 3 4 5; do run alfa${a} caso=percolacao modelo=mc h=1 adapt=2 fs=0 alpha=$a; done
fi
if [ "$PARTE" = hw ] || [ "$PARTE" = todas ]; then
    # Tabela 6 (coluna determinística): hw/H = 0.5..1
    for hw in 2.5 3 3.5 4 4.5 5; do run hw${hw} caso=percolacao modelo=mc h=1 adapt=2 fs=0 hw=$hw; done
fi
