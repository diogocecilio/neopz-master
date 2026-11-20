#!/bin/bash
# ls_etapa3b.sbatch: line sampling (Projects/Cho2010/line_sampling.py) nos casos da etapa 3b: Cho2010 com r = -0,5 e -0,7
# no campo de Nystrom de 2000 modos (o mesmo das ISMC), com as analises pelo SRM (como o direto e as ISMC).
#   1) os FORM pelo GIM com o gradiente modal e os 2000 termos por campo (form_r05_gim2000, form_r07_gim2000; 1 analise
#      por iteracao), em paralelo com o passo 2;
#   2) o line sampling nas direcoes dos 6 FORM prontos (M20, M100 e gim20, r05 e r07);
#   3) o line sampling nas direcoes dos dois FORM de 2000 termos.
# N linhas por direcao (padrao 100), mesma semente em todas (numeros aleatorios comuns). Resultados em
# $SAIDA/ls/<r>_<direcao>/ (ls.json, linhas.csv, analises.csv) e o log no logs/ desta pasta. Retoma de onde parou.
# Uso (no fermi, nesta pasta):  sbatch ls_etapa3b.sbatch
#SBATCH --job-name=lsEtapa3b
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/%x_%j.out
[ -r /etc/profile.d/scratch.sh ] && source /etc/profile.d/scratch.sh
cd "$SLURM_SUBMIT_DIR" || exit 1
source ../../LearningMohrCoulomb2/cluster/modulos.sh
SAIDA="${SAIDA:-${DATA:?}/etapa3b}" PARALELO="${PARALELO:-24}" N="${N:-100}"
BASE="$SAIDA/pontos/base_nys2000.txt" EXE="$PWD/bin/Cho2010" PY="$(cd ../../Cho2010 && pwd)"
[ -r "$BASE" ] || { echo "falta $BASE" >&2; exit 1; }
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PZ_THREADS=0
mkdir -p "$SAIDA/ls"
rc() { [ "$1" = r05 ] && echo -0.5 || echo -0.7; }
ls_dir() {  # ls_dir <r> <direcao> <processos>
    local F="$SAIDA/form_$1_$2/form.json"
    [ -r "$F" ] || { echo "$(date '+%F %T') sem $F"; return; }
    [ -f "$SAIDA/ls/$1_$2/ls.json" ] && { echo "$(date '+%F %T') $1_$2 ja pronto"; return; }
    echo "$(date '+%F %T') line sampling $1_$2"
    python3 -I "$PY/line_sampling.py" "$EXE" "$F" "$SAIDA/ls/$1_$2" n="$N" semente=7 paralelos="$3" K=2000 \
        nystrom="$BASE" metodo=srm
}
echo "$(date '+%F %T') job $SLURM_JOB_ID em $(hostname -s): $PARALELO processos, $N linhas por direcao, SAIDA=$SAIDA"
for r in r05 r07; do
    F="$SAIDA/form_${r}_gim2000"
    [ -f "$F/form.json" ] && continue
    python3 -I "$PY/form_kl.py" "$EXE" "$F" M=2000 rcphi="$(rc $r)" metodo=gim gradiente=modal fator=tri2 nref=4 \
        nystrom="$BASE" iteracoes=15 paralelos=1 > "$SAIDA/form_${r}_gim2000.log" 2>&1 &
done
for r in r05 r07; do
    for d in gim20 M20 M100; do ls_dir $r $d $((PARALELO - 2)); done
done
wait
for r in r05 r07; do ls_dir $r gim2000 "$PARALELO"; done
echo "$(date '+%F %T') resumo:"
python3 -I - "$SAIDA/ls" <<'EOF'
import glob, json, os, sys
for f in sorted(glob.glob(os.path.join(sys.argv[1], "*", "ls.json"))):
    d = json.load(open(f))
    print("%-14s beta %.3f  P_f %.3f%%  CoV %.1f%% (%d linhas)  CoV/linha %.0f%%  analises/linha %.2f  FORM %s  p/ CoV 5%%: %.0f"
          % (os.path.basename(os.path.dirname(f)), d["beta"], 100 * d["P_f"], 100 * d["CoV"], d["linhas"],
             100 * d["CoV_linha"], d["analises_linha"], d["analises_form"], d["analises_cov5"]))
EOF
echo "$(date '+%F %T') fim"
