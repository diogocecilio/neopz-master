#!/bin/sh
# Production runs of the C++ reproduction of Figs. 5, 8 and 9 (Ceron et al., IJNAMG 2025, Sect. 4.3) into
# results/cpp/ (resumable: the commands skip the cases already in their csv, so this script can simply be re-run after
# an interruption), then the comparison tables and figures of scripts/plot_results.py.
#   sh Projects/SlopeSeepageForces/scripts/run_cpp_figures.sh [<SlopeSeepageForces binary>] [threads]
# Wall times of each step are appended to results/cpp/runtimes.txt.
set -e
DIR=$(cd "$(dirname "$0")/.." && pwd)
B=${1:-$DIR/../../build/Projects/SlopeSeepageForces/SlopeSeepageForces}
T=${2:-2}
OUT=$DIR/results/cpp
mkdir -p "$OUT"
cd "$DIR"

run() { # label, command...
    label=$1
    shift
    t0=$(date +%s)
    "$@" > "$OUT/$label.stdout" 2>&1
    t1=$(date +%s)
    echo "$(date -u '+%Y-%m-%d %H:%M:%S')  $label: $((t1 - t0)) s  ($*)" | tee -a "$OUT/runtimes.txt"
}

# Fig. 5: normalized functionals, h_w = H, alpha = 1, 2, 4, 10, beta = 15..90 every 7.5 deg, href = 1
run fig5 "$B" fig5 out=results/cpp/fig5.csv
# Fig. 9: H = 5 m, h_w = H, gamma_w = 9.81, both fields, FE box 50 / 10 / 30 m
run fig9 "$B" fig9 threads="$T"
# Fig. 8: H_ref = 1 m, (c, phi) of Table 1 swapped between the panels, gamma_w = 9.8, both fields, FE box 50 / 10 / 30 m
run fig8 "$B" fig8 threads="$T"
# Fig. 8 with Table 1 as printed
run fig8_table1 "$B" fig8 soil=table1 threads="$T"
# gamma_w evidence: the h_w = 0 ends with 9.8 and 9.81 for both soil sets
for gw in 9.8 9.81; do
    for soil in swapped table1; do
        run "fig8_hw0_${soil}_gw$gw" "$B" fig8 hws=0 soil=$soil gammaw=$gw threads="$T" out=results/cpp/fig8_hw0_gammaw.csv
    done
done
# Fig. 9 with gamma_w = 9.8 (sensitivity)
run fig9_gw9.8 "$B" fig9 gammaw=9.8 threads="$T" out=results/cpp/fig9_gammaw9.8.csv
# comparison with Python and the paper, figures results/fig5.png, fig8.png, fig9.png
run plot python3 scripts/plot_results.py
