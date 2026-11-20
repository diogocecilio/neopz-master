#!/bin/bash
# ls_ceron.sbatch: P_f do rebaixamento rapido do Ceron (2025) com a agua do artigo (contorno=ceron: percolacao
# saturada depois do rebaixamento instantaneo, Eqs. 19-21) por line sampling na direcao do FORM com o gradiente modal.
#   1) FORM (Projects/Cho2010/form_kl.py prog=rebaixamento gradiente=modal) com M = 2000 termos por campo: 1 analise por
#      iteracao, ate 15 iteracoes (pula se $SAIDA/form/form.json ja existe; um FORM interrompido recomeca do zero);
#   2) line sampling (line_sampling.py) com o campo de K = 2000 termos, parando quando CoV(P_f) <= COV (5%), com no
#      minimo NMIN (30) e no maximo N (600) linhas, PARALELO linhas em voo (24: os 12 nucleos com hyperthreading;
#      ~0,22 GB por processo, ~5 GB no total). Retoma de onde parou (linhas.csv); rodar de novo com o CoV ja atingido so
#      refaz o ls.json.
# Caso (o mesmo no FORM e nas linhas): contorno=ceron agua=sat trebaix=0.0001 tfim=20 dt=5 (4 passos: o permanente com
# |dGamma| <= 2e-4 mesmo com k = 0,05 k_v no solo todo), transiente tri3, fator tri1 + 5 niveis (NREF 3 tem vies de
# +4,4% no Gamma). Usa o executavel de bin_ceron/ (./compilar_ceron.sh), nao o bin/ das campanhas.
# Saida em $SAIDA (padrao $DATA/ceron_ls): form/ (form.json), form.log, ls/ (ls.json, linhas.csv, analises.csv,
# perdidas.csv) e o log no logs/ desta pasta. Parar sem perder linhas: touch $SAIDA/ls/PARAR (as linhas em voo terminam).
# Uso (no fermi, nesta pasta):  sbatch ls_ceron.sbatch        (ou COV=0.03 N=1500 sbatch ls_ceron.sbatch)
#SBATCH --job-name=lsCeron
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/%x_%j.out
[ -r /etc/profile.d/scratch.sh ] && source /etc/profile.d/scratch.sh
cd "$SLURM_SUBMIT_DIR" || exit 1
source ../../LearningMohrCoulomb2/cluster/modulos.sh
SAIDA="${SAIDA:-${DATA:?}/ceron_ls}" PARALELO="${PARALELO:-${SLURM_CPUS_PER_TASK:-24}}"
COV="${COV:-0.05}" NMIN="${NMIN:-30}" N="${N:-600}"
EXE="$PWD/bin_ceron/SlopeDrawdownMonteCarlo" PY="$(cd ../../Cho2010 && pwd)"
CASO="contorno=ceron agua=sat trebaix=0.0001 tfim=20 dt=5 transiente=tri3 fator=tri1 nref=5"
[ -x "$EXE" ] || { echo "falta $EXE: rode ./compilar_ceron.sh no no de login" >&2; exit 1; }
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PZ_THREADS=0
mkdir -p "$SAIDA"
echo "$(date '+%F %T') job $SLURM_JOB_ID em $(hostname -s): $PARALELO processos, CoV alvo $COV, $NMIN a $N linhas," \
    "SAIDA=$SAIDA, executavel $(cat bin_ceron/versao.txt 2> /dev/null)"
echo "  caso: $CASO"
if [ ! -f "$SAIDA/form/form.json" ]; then
    echo "$(date '+%F %T') FORM modal M = 2000 (log em $SAIDA/form.log)"
    # shellcheck disable=SC2086
    python3 -I "$PY/form_kl.py" "$EXE" "$SAIDA/form" prog=rebaixamento gradiente=modal M=2000 iteracoes=15 \
        paralelos=1 $CASO > "$SAIDA/form.log" 2>&1 || { echo "o FORM falhou:"; tail -n 20 "$SAIDA/form.log"; exit 1; }
fi
tail -n 5 "$SAIDA/form.log"
echo "$(date '+%F %T') line sampling"
python3 -I "$PY/line_sampling.py" "$EXE" "$SAIDA/form/form.json" "$SAIDA/ls" K=2000 semente=7 cov="$COV" \
    nmin="$NMIN" n="$N" paralelos="$PARALELO"
echo "$(date '+%F %T') resumo:"
python3 -I - "$SAIDA" <<'EOF'
import json, os, sys
f = json.load(open(os.path.join(sys.argv[1], "form", "form.json")))
print("FORM: beta %.4f  P_f %.3f%%  (%s, %d iteracoes, %d analises)  alfa^2: c %.3f phi %.3f k %.3f"
      % (f["beta"], 100 * f["P_f"], "convergiu" if f["convergiu"] else "NAO convergiu", f["iteracoes"],
         f["avaliacoes"], f["alfa2_c"], f["alfa2_phi"], f["alfa2_k"]))
arq = os.path.join(sys.argv[1], "ls", "ls.json")
if not os.path.exists(arq):
    sys.exit("LS:   sem ls.json (o line sampling falhou antes de terminar: veja o log acima)")
d = json.load(open(arq))
if d["P_f"] is None:
    print("LS:   sem P_f (%d linhas, %d perdidas); parada: %s" % (d["linhas"], d["perdidas"], d["motivo"]))
    sys.exit(0)
cov = d["CoV"] if d["CoV"] is not None else float("nan")
print("LS:   P_f %.3f%%  IC 95%% [%.3f%%, %.3f%%]  CoV %.1f%%  %d linhas (%d perdidas, %.1f%%)  %d analises (%.2f por linha)"
      % (100 * d["P_f"], 100 * d["IC95"][0], 100 * d["IC95"][1], 100 * cov, d["linhas"], d["perdidas"],
         100 * d.get("fracao_perdida", 0.), d["analises"], d["analises_linha"]))
if d.get("fracao_perdida", 0.) > 0.01:
    print("      AVISO: mais de 1% de linhas perdidas (perdidas.csv): a P_f pode sair baixa")
print("      parada: %s   (artigo: P_f = 11,50%%)" % d["motivo"])
EOF
echo "$(date '+%F %T') fim"
