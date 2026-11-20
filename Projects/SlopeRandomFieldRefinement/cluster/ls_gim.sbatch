#!/bin/bash
# ls_gim.sbatch: line sampling (Projects/Cho2010/line_sampling.py) do Cho2010 com o fator pelo GIM (metodo=gim, Gamma;
# falha: Gamma < 1), pareado com o line sampling SRM da etapa 3b (rodado no PC: ~/projects/runs/etapa3b_ls/ls): mesmas
# direcoes, mesmo campo de Nystrom de 2000 modos, mesma semente 7 (a linha i usa o mesmo Z), so metodo=gim a mais.
# Direcoes (r05 e r07; FORM ja prontos, nenhum FORM novo):
#   M20      FORM com M = 20 e diferencas centrais do SRM ($DATA/etapa3b/form_<r>_M20/form.json, ja no fermi)
#   gim2000  FORM pelo GIM com o gradiente modal e M = 2000 (rodado no PC; suba antes os dois form.json para
#            $SAIDA/form_<r>_gim2000/form.json: o comando esta no README, secao etapa3b_gim)
# Parada pelo CoV(P_f) = 5% (cov=), com no minimo NMIN linhas (padrao 100: as 100 linhas do SRM e da analise limite,
# todas pareadas) e no maximo N (600). Resultados em $SAIDA/ls/<r>_<direcao>/ (ls.json, linhas.csv, analises.csv); o log
# no logs/ desta pasta. Retoma de onde parou.
# Uso (no fermi, nesta pasta):  sbatch ls_gim.sbatch      (ambiente: SAIDA, PARALELO (24), N (600), NMIN (100), COV)
#SBATCH --job-name=lsGim
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/%x_%j.out
[ -r /etc/profile.d/scratch.sh ] && source /etc/profile.d/scratch.sh
cd "$SLURM_SUBMIT_DIR" || exit 1
source ../../LearningMohrCoulomb2/cluster/modulos.sh
SAIDA="${SAIDA:-${DATA:?}/etapa3b_gim}" PARALELO="${PARALELO:-24}" N="${N:-600}" NMIN="${NMIN:-100}" COV="${COV:-0.05}"
E3B="${E3B:-$DATA/etapa3b}"
BASE="$E3B/pontos/base_nys2000.txt" EXE="$PWD/bin/Cho2010" PY="$(cd ../../Cho2010 && pwd)"
[ -r "$BASE" ] || { echo "falta $BASE" >&2; exit 1; }
# cada analise do LS e um processo novo que le a base (525 MB): copia para o disco local do no, se couber, para nao
# ler do /mnt/data a cada analise (24 em paralelo); se nao der, usa a do /mnt/data
LOCAL="${TMPDIR:-/tmp}/lsGim_${SLURM_JOB_ID:-$$}"
if mkdir -p "$LOCAL" && cp "$BASE" "$LOCAL/" 2> /dev/null; then
    BASE="$LOCAL/$(basename "$BASE")"; trap 'rm -rf "$LOCAL"' EXIT
else
    rm -rf "$LOCAL"; echo "$(date '+%F %T') sem espaco local: base lida de $BASE"
fi
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PZ_THREADS=0
mkdir -p "$SAIDA/ls"
echo "$(date '+%F %T') job $SLURM_JOB_ID em $(hostname -s): $PARALELO processos, CoV $COV, linhas $NMIN a $N," \
    "SAIDA=$SAIDA, base $BASE, versao $(cat bin/versao.txt 2> /dev/null || echo -)"
for c in r05_gim2000 r07_gim2000 r05_M20 r07_M20; do
    d=${c#*_}
    if [ "$d" = M20 ]; then F="$E3B/form_${c}/form.json"; else F="$SAIDA/form_${c}/form.json"; fi
    [ -r "$F" ] || { echo "$(date '+%F %T') sem $F: $c fica de fora"; continue; }
    [ -f "$SAIDA/ls/$c/ls.json" ] && { echo "$(date '+%F %T') $c ja pronto"; continue; }  # o ls.json sai no fim
    echo "$(date '+%F %T') line sampling GIM $c ($F)"
    python3 -I "$PY/line_sampling.py" "$EXE" "$F" "$SAIDA/ls/$c" n="$N" semente=7 paralelos="$PARALELO" K=2000 \
        nystrom="$BASE" metodo=gim cov="$COV" nmin="$NMIN" || echo "$(date '+%F %T') ERRO em $c"
done
echo "$(date '+%F %T') resumo:"
python3 -I - "$SAIDA/ls" <<'EOF'
import glob, json, os, sys
for f in sorted(glob.glob(os.path.join(sys.argv[1], "*", "ls.json"))):
    d = json.load(open(f))
    if d.get("P_f") is None:
        print(os.path.basename(os.path.dirname(f)), "sem resultado"); continue
    print("%-14s beta %.3f  P_f %.4f%%  CoV %.1f%% (%d linhas)  analises %d  analises/linha %.2f  motivo: %s"
          % (os.path.basename(os.path.dirname(f)), d["beta"], 100 * d["P_f"], 100 * (d["CoV"] or float("nan")), d["linhas"],
             d["analises"], d["analises_linha"], d.get("motivo", "-")))
EOF
echo "$(date '+%F %T') fim"
