#!/bin/bash
# rebaixamento instantaneo (trebaix = 1e-4 h), deterministico, agua saturada (u-p) e nao saturada (u-pw-pa)
B=/home/diogo/projects/neopz-master-build-release/Projects/SlopeRandomFieldRefinement
export PZ_THREADS=2 OMP_NUM_THREADS=1
for agua in sat unsat; do
  for par in "0.0001 0.0001" "0.1 0.01" "1 0.1" "10 0.5"; do
    set -- $par; t=$1; dt=$2
    s=${agua}_t$t
    [ -f ${s}_-1_-1.csv ] || nice $B/SlopeDrawdownMonteCarlo det=1 n=0 figuras=-1 trebaix=0.0001 tfim=$t dt=$dt agua=$agua saida=$s > $s.log 2>&1
    suc=total; [ $agua = unsat ] && suc=sw
    [ -f up_$s.csv ] || nice $B/SlopeDrawdownUP trebaix=0.0001 tfim=$t dt=$dt agua=$agua succao=$suc saida=up_$s > up_$s.log 2>&1
    echo "$s feito $(date +%T)"
  done
done
echo FIM
