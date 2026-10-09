#!/bin/sh
# Convergence study of the FEM gravity-increase factor (FEMStability.h, command 'fs') for the data of the paper's
# Fig. 9 (H = 5 m, c = 10 kPa, phi = 30 deg, gamma = 20, gamma_w = 9.81, h_w = H, alpha = 1): refinement cycles of the
# plastic zone, marking fraction, initial mesh size, uniform refinement, stability-domain extents and driver settings,
# with the FE field -grad u'_FE (hydraulic box fixed in metres: 50 / 10 / 30 m = 10 / 2 / 6 H, u = 0 on the left
# side and the base, href = 1) and one case with the analytical field K^-1 v'_opt.
#   sh Projects/SlopeSeepageForces/scripts/run_fem_convergence.sh [<SlopeSeepageForces binary>] [<case-name filter>]
# One log per case in results/fem/convergence/<name>.log; a case whose log ends with 'total time' is skipped, so the
# script can simply be re-run after an interruption (an interrupted case is redone from the start). The environment
# variable CPUS (e.g. CPUS=2,3) pins the runs to these CPUs with taskset. Summary: scripts/fem_convergence_table.py.
DIR=$(cd "$(dirname "$0")/.." && pwd)
B=${1:-$DIR/../../build/Projects/SlopeSeepageForces/SlopeSeepageForces}
FILTER=${2:-}
OUT=$DIR/results/fem/convergence
mkdir -p "$OUT"
# Fig. 9 data; FE field on the paper's box in metres; stability box: sa H + H / tan(beta) to the left of O and below
# T, 2 H (= the right side of the hydraulic box) to the right of T
FIG9="H=5 c=10 phi=30 gamma=20 gammaw=9.81 hw=1 alpha=1 hleft=10 hright=2 hdepth=6 href=1 sright=2"
run() { # name, fs options...
    name=$1
    shift
    case "$name" in *"$FILTER"*) ;; *) return ;; esac
    log=$OUT/$name.log
    if [ -f "$log" ] && grep -q '^total time' "$log"; then return; fi
    t0=$(date +%s)
    # stdbuf -oL: line-buffered output (the continuation trials of SlopeAnalysis are printed without flush)
    if [ -n "$CPUS" ]; then taskset -c "$CPUS" stdbuf -oL "$B" fs $FIG9 "$@" > "$log" 2>&1
    else stdbuf -oL "$B" fs $FIG9 "$@" > "$log" 2>&1; fi
    echo "$(date -u '+%Y-%m-%d %H:%M:%S')  $name: $(($(date +%s) - t0)) s  (fs $*)" | tee -a "$OUT/runtimes.txt"
}
# 1. refinement cycles of the plastic zone (2 % marking), beta = 60, FE field: cycles 0-5 (this run was stopped during
#    cycle 5, more than 22571 equations: its log holds cycles 0-4; b60_m02_n5 at the end repeats it to cycle 5)
run b60_m02 beta=60 nref=5 mark=0.02
# 2. driver settings (beta = 60, cycles 0-2): Newton cap 30 (SlopeMohrCoulomb / SlopeDrawdown) and 200, continuation
#    tolerance 1e-3 and 5e-3 (default 100 and 2e-3: b60_m02)
run b60_mn30_n2 beta=60 nref=2 mark=0.02 maxnewton=30
run b60_mn200_n2 beta=60 nref=2 mark=0.02 maxnewton=200
run b60_tol1e3_n2 beta=60 nref=2 mark=0.02 tolfs=0.001
run b60_tol5e3_n2 beta=60 nref=2 mark=0.02 tolfs=0.005
# 3. beta = 90 and 30, FE field, and beta = 60 with K^-1 v'_opt: cycles 0-4
run b90_m02 beta=90 nref=4 mark=0.02
run b30_m02 beta=30 nref=4 mark=0.02
run b60_vopt_m02 beta=60 nref=4 mark=0.02 water=analytical
# 4. marking fraction (beta = 60): 5 % and 10 % (cycles 0-4), 1 % (cycles 0-3)
run b60_m05 beta=60 nref=4 mark=0.05
run b60_m10 beta=60 nref=4 mark=0.1
run b60_m01 beta=60 nref=3 mark=0.01
# 5. initial mesh (beta = 60): sizes 0.125 H at O, W, T and along the face (default 0.25 H) with cycles; uniform
#    refinements of the default mesh without cycles (lambda versus h)
run b60_s0125 beta=60 nref=3 mark=0.02 sh0=0.125 shs=0.125
run b60_sref1 beta=60 nref=0 sref=1
run b60_sref2 beta=60 nref=0 sref=2
# 6. stability-domain extents (cycles 0-2; default sa = 2, right side at the hydraulic box, 2 H from T): sa = 1.5 and 3
#    (left of O and below T), right side 1.5 H from T
run b30_sa15_n2 beta=30 nref=2 mark=0.02 sa=1.5
run b30_sa3_n2 beta=30 nref=2 mark=0.02 sa=3
run b30_sr15_n2 beta=30 nref=2 mark=0.02 sright=1.5
run b90_sa15_n2 beta=90 nref=2 mark=0.02 sa=1.5
run b90_sa3_n2 beta=90 nref=2 mark=0.02 sa=3
run b90_sr15_n2 beta=90 nref=2 mark=0.02 sright=1.5
# 7. 5 % marking for beta = 30 and 90 (cycles 0-4) and the analytical field (cycles 0-4)
run b30_m05 beta=30 nref=4 mark=0.05
run b90_m05 beta=90 nref=4 mark=0.05
run b60_vopt_m05 beta=60 nref=4 mark=0.05 water=analytical
# 8. cycle 5 with 5 % marking (beta = 60). With 2 % marking cycle 5 has more than 22571 equations and takes about
#    1 h alone (serial skyline LDLt; b60_m02 was stopped there): not run, uncomment to add it
run b60_m05_n5 beta=60 nref=5 mark=0.05
# run b60_m02_n5 beta=60 nref=5 mark=0.02
