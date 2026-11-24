#!/bin/bash
# casoC_diag.sbatch: diagnostico do Newton do caso C (Liu e Wang 2023) deterministico (as medias, CoV 0), chuva de
# 200 h a I = 36 e 4 mm/h, so o transiente (so_transiente=1), com vg_tab=1e-3 e o Newton da coluna (newton=1e-6,80,1,
# dt_min=0.00625). No PC, as 14 variantes do T3 e as 2 com esse Newton param em 33-37 h (ciclo do Newton num bolsao
# saturado interior, com a malha tri3 de 1,25 m). Aqui, 7 variantes x 2 I, todas em paralelo num no (PZ_THREADS=1):
#   ref       tri3 (a de antes)
#   kw4, kw3  tri3 com K_w = 2,2e4 e 2,2e3 kPa (armazenamento da zona saturada 100 e 1000 vezes maior: imita o esqueleto
#             deformavel do modelo acoplado)
#   tri4      malha de 0,625 m;  tri4kw4: tri4 com K_w = 2,2e4
#   cam       camadas de 0,25 m ate 5 m de profundidade (cam0.25p5);  camkw4: cam com K_w = 2,2e4
#   kw2       tri3 com K_w = 2,2e2;  tri5 (0,31 m) e tri5kw4
# Os 4 jobs de submeter_casoC.sh dividem as variantes (VARS) e usam PZT threads por rodada (PZ_THREADS).
# Pergunta: quem chega a 200 h (o armazenamento, a malha ou os dois).
# Saida em $SAIDA/<id>_I<I>/ (log.txt, nt_*.txt: newton_trace, tr_*.txt: pressoes se chegar ao fim); resumo em
# $SAIDA/resumo.txt (tambem: ./casoC_diag.sbatch resumo). Uso: sbatch casoC_diag.sbatch (depois do compilar_casoC.sh).
# Ambiente: SAIDA ($DATA/casoC_diag), EXE (bin_casoC/SlopeRainMonteCarlo), VARS (as 7 primeiras, separadas por virgula),
# INTS (36,4), PZT (1)
#SBATCH --job-name=casoCd
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=48G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/%x_%j.out
[ -r /etc/profile.d/scratch.sh ] && source /etc/profile.d/scratch.sh
cd "${SLURM_SUBMIT_DIR:-$(dirname "$(readlink -f "$0")")}" || exit 1
[ -n "${SLURM_JOB_ID:-}" ] && source ../../LearningMohrCoulomb2/cluster/modulos.sh
DADOS="${DATA:-/mnt/data/${USER:-$(id -un)}}"
SAIDA="${SAIDA:-$DADOS/casoC_diag}"
EXE="${EXE:-$PWD/bin_casoC/SlopeRainMonteCarlo}"
VARS="${VARS:-ref,kw4,kw3,tri4,tri4kw4,cam,camkw4}" INTS="${INTS:-36,4}"
TALUDE=0,0,74,0,74,20,43,20,31,30,0,30
BASE="geometria=$TALUDE na=14.18918919 na_beta=7.696051722 carga_esq=20 carga_dir=10 vg=0.358,0,0.149011,1.289"
BASE+=" ks=3e-6 q0=0.1369863 chi=ponto phib=15 c=5 phi=32 gamma=20 covc=0 covphi=0 covk=0 rcphi=0 fator=tri1 nref=4"
BASE+=" det=1 n=0 chuva=200 tempos=24,48,72,100,150,200 so_transiente=1 vg_tab=1e-3 newton=1e-6,80,1 dt_min=0.00625"
declare -A V=([ref]="transiente=tri3" [kw4]="transiente=tri3 kw=2.2e4" [kw3]="transiente=tri3 kw=2.2e3"
              [tri4]="transiente=tri4" [tri4kw4]="transiente=tri4 kw=2.2e4"
              [cam]="transiente=cam0.25p5" [camkw4]="transiente=cam0.25p5 kw=2.2e4"
              [kw2]="transiente=tri3 kw=2.2e2" [tri5]="transiente=tri5" [tri5kw4]="transiente=tri5 kw=2.2e4")

resumo() {
    printf "%-14s %-9s %-9s %-8s %s\n" id estado ultimo_t minutos detalhe
    for d in "$SAIDA"/*_I*/; do
        [ -f "$d/log.txt" ] || continue
        local id; id="$(basename "$d")"
        local ini fim est ult det=""
        ini=$(cat "$d/inicio" 2>/dev/null || echo 0); fim=$(cat "$d/fim" 2>/dev/null || date +%s)
        ult=$(grep -o 'transiente em [0-9.]* h' "$d/log.txt" | tail -1 | awk '{print $3}')
        if grep -q 'o transiente parou em' "$d/log.txt"; then
            est=parou; det="$(grep -m1 '  o transiente parou em' "$d/log.txt" | sed 's/^  o transiente //')"
        elif [ -f "$d/fim" ]; then
            if [ "$ult" = 200 ]; then est=completo; else est="saiu($(cat "$d/codigo"))"; fi
        else est=rodando; fi
        printf "%-14s %-9s %-9s %-8s %s\n" "$id" "$est" "${ult:--}" "$(( (fim - ini) / 60 ))" "$det"
    done
}
if [ "${1:-}" = resumo ]; then resumo; exit 0; fi

[ -x "$EXE" ] || { echo "falta $EXE (rode ./compilar_casoC.sh)"; exit 1; }
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PZ_THREADS="${PZT:-1}"
mkdir -p "$SAIDA"
echo "$(date '+%F %T') [${SLURM_JOB_ID:-local}] inicio: $VARS x I $INTS; $(cat "$(dirname "$EXE")/versao.txt" 2>/dev/null)" >> "$SAIDA/eventos.log"
IFS=, read -r -a vs <<< "$VARS"; IFS=, read -r -a is <<< "$INTS"
for v in "${vs[@]}"; do
    [ -n "${V[$v]:-}" ] || { echo "variante desconhecida: $v"; continue; }
    for I in "${is[@]}"; do
        (
            id="${v}_I$I"; d="$SAIDA/$id"; mkdir -p "$d" && cd "$d" || exit 1
            date +%s > inicio
            # shellcheck disable=SC2086
            echo "$EXE $BASE ${V[$v]} I=$I grava_pressoes=tr_$id.txt newton_trace=nt_$id.txt saida=$id" > comando.txt
            # shellcheck disable=SC2086
            "$EXE" $BASE ${V[$v]} I="$I" grava_pressoes="tr_$id.txt" newton_trace="nt_$id.txt" saida="$id" > log.txt 2>&1
            echo $? > codigo; date +%s > fim
            echo "$(date '+%F %T') [${SLURM_JOB_ID:-local}] $id terminou (codigo $(cat codigo))" >> "$SAIDA/eventos.log"
        ) &
    done
done
wait
resumo > "$SAIDA/resumo.txt"
echo "$(date '+%F %T') [${SLURM_JOB_ID:-local}] fim" >> "$SAIDA/eventos.log"
cat "$SAIDA/resumo.txt"
