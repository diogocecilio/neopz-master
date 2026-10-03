#!/bin/bash
# Rebaixamento acoplado (u-p, Biot) do talude de referência para várias durações adimensionais T_d = c_v t_d / H².
# Uso: rebaixamento.sh <executável> <diretório> <modelo mc|mcc> <processos> [argumentos extras do comando rebaixamento]
# Cada T_d grava rebaixamento_<modelo>_Td<T_d>.csv (FS e Γ com a poropressão congelada em cada tempo).
set -u
EXE=$(readlink -f "$1"); OUT=$2; MODELO=$3; P=$4; shift 4
mkdir -p "$OUT"; cd "$OUT"
for td in 0.001 0.01 0.1 1 10; do
    nome="rebaixamento_${MODELO}_Td${td}"
    [ -s "$nome.log" ] && grep -q "tempo total" "$nome.log" && continue
    echo "$EXE rebaixamento caso=percolacao modelo=$MODELO Td=$td saida=$nome.csv $* > $nome.log 2>&1"
done | xargs -P "$P" -I{} bash -c '{}'
