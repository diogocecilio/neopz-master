#!/bin/sh
# FEM gravity-increase stability factor (FEMStability.h, command 'fembatch') of the cases of the paper's Figs. 9 and 8,
# as an independent check of the limit-analysis curves (Ceron et al., IJNAMG 2025, Sect. 4.3), with the production
# settings of the convergence study (defaults of 'fembatch', README section 4.4):
#   Fig. 9: alpha = 1, 5, 10 x beta = 30, 45, 60, 75, 90 x {FE, vopt}            -> results/cpp/fem_fig9.csv (30 cases)
#   Fig. 8: London 30/60, Israeli 35/60 ((c, phi) swapped) x h_w / H = 0, 0.2, 0.5, 1 x {FE, vopt}, h_w = 0 once
#                                                                                -> results/cpp/fem_fig8.csv (28 runs)
# Each row also holds the limit analysis of the same case. Resumable: 'fembatch' skips the cases already in its csv,
# so this script can simply be re-run after an interruption (the case in progress is redone from the start).
#   sh Projects/SlopeSeepageForces/scripts/run_fem_batch.sh [<SlopeSeepageForces binary>] [9|8|9,8] [extra options]
# The environment variable CPUS (e.g. CPUS=2,3) pins the run to these CPUs with taskset (the elastoplastic assembly
# of SlopeAnalysis.h uses one thread per CPU it sees). Wall times of each step: results/cpp/fem_runtimes.txt.
DIR=$(cd "$(dirname "$0")/.." && pwd)
B=${1:-$DIR/../../build/Projects/SlopeSeepageForces/SlopeSeepageForces}
FIGS=${2:-9,8}
[ $# -ge 1 ] && shift
[ $# -ge 1 ] && shift
OUT=$DIR/results/cpp
mkdir -p "$OUT"
cd "$DIR" || exit 1
for fig in $(echo "$FIGS" | tr ',' ' '); do
    t0=$(date +%s)
    # stdbuf -oL: line-buffered stdout (the continuation trials of SlopeAnalysis are printed without flush)
    if [ -n "$CPUS" ]; then taskset -c "$CPUS" stdbuf -oL "$B" fembatch fig="$fig" "$@" >> "$OUT/fem_fig$fig.stdout" 2>&1
    else stdbuf -oL "$B" fembatch fig="$fig" "$@" >> "$OUT/fem_fig$fig.stdout" 2>&1; fi
    rc=$?
    echo "$(date -u '+%Y-%m-%d %H:%M:%S')  fembatch fig=$fig: $(($(date +%s) - t0)) s, exit code $rc  ($*)" | tee -a "$OUT/fem_runtimes.txt"
done
