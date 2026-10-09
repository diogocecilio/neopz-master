#!/bin/sh
# Diagnosis of the FEM batch outliers (results/cpp/fem_fig9.csv, fem_fig8.csv; scripts/fem_batch_table.py): the cases
# where the order-1 extrapolation Gamma_FEM = 2 lambda_3 - lambda_2 is more than ~3 % below the limit-analysis upper
# bound Gamma_LA of the same case. Each outlier is re-run with the command 'fs' and the production settings of
# 'fembatch' (mark 5 %, initial mesh 0.25 H, sa = 2, Newton cap 100, tolfs 0.002, href 1), with more refinement cycles
# (nref = 4 or 5: cycles 0-3 must reproduce the batch row), and, for Fig. 9, with a larger hydraulic box (right side
# 30 m instead of 10 m, so that the stability box is no longer capped at 2 H right of the toe; the FE field changes
# slightly: the limit analysis of the same box is run too).
#   sh Projects/SlopeSeepageForces/scripts/run_fem_outliers.sh [<SlopeSeepageForces binary>] [<case-name filter>]
# One log per case in results/fem/outliers/<name>.log; a case whose log ends with 'total time' (fs) or holds a
# 'Gamma' line (la) is skipped, so the script can simply be re-run after an interruption. CPUS=0,1 pins the runs with
# taskset. Summary: scripts/fem_outliers_table.py -> results/fem/outliers/summary.csv.
DIR=$(cd "$(dirname "$0")/.." && pwd)
B=${1:-$DIR/../../build/Projects/SlopeSeepageForces/SlopeSeepageForces}
FILTER=${2:-}
OUT=$DIR/results/fem/outliers
mkdir -p "$OUT"
# production FEM settings of fembatch (FEMCommands.cpp FEMBatchSettings)
FEM="mark=0.05 sa=2 sh0=0.25 shs=0.25 sgrade=0.25 shmax=1 maxnewton=100 tolfs=0.002 href=1"
# Fig. 9 data (fembatch fig=9): hydraulic box 50 / 10 / 30 m, stability box capped at its right side (2 H from T)
FIG9="H=5 beta=30 c=10 phi=30 gamma=20 gammaw=9.81 hw=1"
# Fig. 8 panels at the reference height H_ref of the batch rows (H_crit of the limit analysis, 3 digits; hydraulic box
# 50 / 10 / 30 H_ref = the defaults hleft, hright, hdepth), (c, phi) of Table 1 swapped between the panels
ISR35="beta=35 c=6 phi=32 gamma=18 gammaw=9.8 alpha=1"
LON30="beta=30 c=11.7 phi=24.7 gamma=18 gammaw=9.8 alpha=1"
run() { # name, command (fs | la), options...
    name=$1
    cmd=$2
    shift
    shift
    case "$name" in *"$FILTER"*) ;; *) return ;; esac
    log=$OUT/$name.log
    if [ -f "$log" ] && { grep -q '^total time' "$log" || grep -q '^Gamma' "$log"; }; then return; fi
    t0=$(date +%s)
    echo "$(date -u '+%Y-%m-%d %H:%M:%S')  $name: start  ($cmd $*)" >> "$OUT/runtimes.txt"
    if [ -n "$CPUS" ]; then taskset -c "$CPUS" stdbuf -oL "$B" "$cmd" "$@" > "$log" 2>&1
    else stdbuf -oL "$B" "$cmd" "$@" > "$log" 2>&1; fi
    echo "$(date -u '+%Y-%m-%d %H:%M:%S')  $name: $(($(date +%s) - t0)) s  ($cmd $*)" | tee -a "$OUT/runtimes.txt"
}
# --- group A (Fig. 9, FE field, beta = 30, alpha = 10 and 5) ---------------------------------------------------------
# limit analysis of the same cases with the production box and with the right side at 30 m (threads=2 as fig9)
run A_la_a10_b30_box10 la $FIG9 alpha=10 water=fe hboxm=50,10,30 href=1 threads=2
run A_la_a10_b30_box30 la $FIG9 alpha=10 water=fe hboxm=50,30,30 href=1 threads=2
run A_la_a5_b30_box10 la $FIG9 alpha=5 water=fe hboxm=50,10,30 href=1 threads=2
run A_la_a5_b30_box30 la $FIG9 alpha=5 water=fe hboxm=50,30,30 href=1 threads=2
# alpha = 10: one more cycle (0-4) with the production box; cycles 0-3 with the right side at 30 m (stability box
# 3.73 H right of T, uncapped), the plastic strain of every cycle in vtk
run A_fs_a10_b30_n4 fs $FIG9 alpha=10 water=fe hboxm=50,10,30 sright=2 $FEM nref=4 vtk=$OUT/A_fs_a10_b30_n4
run A_fs_a10_b30_box30 fs $FIG9 alpha=10 water=fe hboxm=50,30,30 $FEM nref=3
# alpha = 5: one more cycle (0-4) with the production box
run A_fs_a5_b30_n4 fs $FIG9 alpha=5 water=fe hboxm=50,10,30 sright=2 $FEM nref=4
# --- group B (Fig. 8, beta = 30 / 35 panels) --------------------------------------------------------------------------
# Israeli 35, h_w = 0 (f = 0, gamma' = gamma - gamma_w; H_ref 229 m): cycles 0-5
run B_fs_isr35_hw0_n5 fs H=229 hw=0 $ISR35 water=none $FEM nref=5
# London 30, h_w = 0 (H_ref 157 m): cycles 0-4
run B_fs_lon30_hw0_n4 fs H=157 hw=0 $LON30 water=none $FEM nref=4
# London 30, K^-1 v'_opt, h_w / H = 0.2 (H_ref 58.7 m): cycles 0-4
run B_fs_lon30_vopt02_n4 fs H=58.7 hw=0.2 $LON30 water=analytical lm=10 $FEM nref=4
# --- group A, mesh family: alpha = 10 with a finer initial mesh (0.125 H at O, W, T and along the face), cycles 0-2:
#     its cycle 2 has the element size of cycle 3 of the production mesh in the plastic zone
run A_fs_a10_b30_s0125_n2 fs $FIG9 alpha=10 water=fe hboxm=50,10,30 sright=2 mark=0.05 sa=2 sh0=0.125 shs=0.125 sgrade=0.25 shmax=1 maxnewton=100 tolfs=0.002 href=1 nref=2
