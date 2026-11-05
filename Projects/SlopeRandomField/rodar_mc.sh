#!/bin/bash
# rodar_mc.sh <primeira> <ultima>: roda as realizacoes primeira..ultima do Monte Carlo de uma combinacao (solo, chuva e
# estatistica nas variaveis de ambiente), cada uma num processo de 1 CPU, PARALELO ao mesmo tempo (padrao: os nucleos
# fisicos da maquina; com hyper-threading, dois processos por nucleo deixam cada um ~2x mais lento).
# Cada realizacao fica em $SAIDA/<combinacao>/r<nnnnn>/ (fs.csv, resultado.txt, realizacao.txt, run.log); o resumo de
# todas vai para $SAIDA/<combinacao>/mc_fs.csv (FS por instante) e mc_resumo.csv. Retomavel: realizacoes com
# resultado.txt sao puladas. ./resumo_mc.sh <pasta da combinacao> calcula P_f e a estatistica do FS.
# Ambiente: SOLO (f10|f50|f100), I_MM_H (9, 80 ou o k_s do solo em mm/h: 360, 36, 3.6), VARIANTE (iii, i), MALHA
# (rahardjo_a45.msh ou rahardjo_a45_fina.msh), FS_H (instantes do FS, horas), COV_C, COV_PHI, COV_K, LX, LY, ALFA, KL_M,
# SEMENTE, CAMPOS=0 (deterministico), EXE, SAIDA, PARALELO, COMB (nome da combinacao). Pos-processamento: CAMPOS_VTK=1
# grava campos.vtk de cada realizacao, VTK_H=24 grava o mecanismo de ruptura em t = 24 h, e SO_CAMPOS=1 so gera os campos
# (campos.vtk e campos.csv) das realizacoes pedidas, sem simular: a mesma SEMENTE e o mesmo numero refazem o mesmo campo.
set -u
PKG="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
[ $# -eq 2 ] || { sed -n '2,14p' "$0"; exit 2; }
export EXE="${EXE:-$PKG/bin/SlopeRandomField}"
[ -x "$EXE" ] || { echo "executavel $EXE nao encontrado (compile o alvo SlopeRandomField; ver README.md)" >&2; exit 2; }
export MESH="${MESH:-$PKG/../LearningMohrCoulomb2/cluster/malhas/${MALHA:-rahardjo_a45.msh}}"
export SOLO="${SOLO:-f50}" I_MM_H="${I_MM_H:-36}" VARIANTE="${VARIANTE:-iii}"
export MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
if [ "${CAMPOS:-1}" = 0 ]; then COMB="${COMB:-${SOLO}_I${I_MM_H}_${VARIANTE}_det}"; else COMB="${COMB:-${SOLO}_I${I_MM_H}_${VARIANTE}}"; fi
export D="${SAIDA:-$PKG/resultados_mc}/$COMB"
mkdir -p "$D"
NUC="$(lscpu -p=CORE,SOCKET 2>/dev/null | grep -v '^#' | sort -u | wc -l)"
P="${PARALELO:-${NUC:-1}}"
[ "$P" -ge 1 ] 2>/dev/null || P=1
[ -s "$D/mc_fs.csv" ] || echo "real,t_h,FS_SRM,seg" > "$D/mc_fs.csv"
[ -s "$D/mc_resumo.csv" ] || echo "real,FS_min,t_min_h,situacao,passos,cortes,newton,seg" > "$D/mc_resumo.csv"
echo "$(date '+%F %T') $COMB: realizacoes $1..$2, $P ao mesmo tempo, $(hostname -s)" | tee -a "$D/execucoes.txt"

# uma realizacao; ao terminar, acrescenta as linhas ao resumo (flock: varios processos escrevem nos mesmos arquivos)
uma() {
    local n="$1" r
    r="$D/$(printf 'r%05d' "$n")"
    mkdir -p "$r" && cd "$r" || return 1
    if [ "${SO_CAMPOS:-0}" = 1 ]; then # so os campos desta realizacao, em VTK e CSV, sem simulacao
        REAL="$n" CAMPOS_VTK=1 CAMPOS_CSV=1 "$EXE" > campos.log 2>&1 && echo "real $n: campos em $r/campos.vtk"
        return 0
    fi
    [ -s resultado.txt ] && return 0
    REAL="$n" "$EXE" > run.log 2> run.err
    [ -s resultado.txt ] || { echo "real $n: sem resultado.txt (veja $r/run.err)" >&2; return 1; }
    (
        flock 9
        tail -n +2 fs.csv | sed "s/^/$n,/" >> "$D/mc_fs.csv"
        tr ' ' '\n' < resultado.txt | awk -F= -v n="$n" '{ v[$1] = $2 } END {
            print n "," v["FS_min"] "," v["t_min_h"] "," v["situacao"] "," v["passos"] "," v["cortes"] "," v["newton"] "," v["seg"] }' \
            >> "$D/mc_resumo.csv"
    ) 9> "$D/.trava"
    echo "$(date '+%T') real $n: $(cat resultado.txt)"
}
export -f uma
seq "$1" "$2" | xargs -P "$P" -I{} bash -c 'uma {}'
echo "$(date '+%F %T') $COMB: fim" | tee -a "$D/execucoes.txt"
