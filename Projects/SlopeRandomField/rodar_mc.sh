#!/bin/bash
# rodar_mc.sh <primeira> <ultima>: roda as realizacoes primeira..ultima do Monte Carlo de uma combinacao (solo, chuva e
# estatistica nas variaveis de ambiente), cada uma num processo de 1 CPU, PARALELO ao mesmo tempo (padrao: os nucleos
# fisicos da maquina; com hyper-threading, dois processos por nucleo deixam cada um ~2x mais lento).
# Cada realizacao fica em $SAIDA/<combinacao>/r<nnnnn>/ (resumo.csv, fs.csv, passos.csv, fs_tentativas.csv, pw.csv,
# mecanismo.csv, campos.csv, resultado.txt, realizacao.txt, run.log); no fim, junta_mc.sh junta as concluidas em
# $SAIDA/<combinacao>/mc_resumo.csv e mc_fs.csv. Retomavel: realizacoes com resultado.txt sao puladas, e uma realizacao
# interrompida (Ctrl-C, queda de energia) continua do seu checkpoint.bin. ./resumo_mc.sh <pasta> calcula P_f.
# Ambiente: SOLO (f10|f50|f100), I_MM_H (9, 80 ou o k_s do solo em mm/h: 360, 36, 3.6), VARIANTE (iii, i), MALHA
# (rahardjo_a45.msh ou rahardjo_a45_fina.msh), FS_H (instantes do FS, horas), COV_C, COV_PHI, COV_K, LX, LY, ALFA, KL_M,
# SEMENTE, CAMPOS=0 (deterministico), EXE, SAIDA, PARALELO, COMB (nome da combinacao). CSV opcionais: PW_CSV, MEC_CSV,
# CAMPOS_CSV (1 = grava, o padrao). VTK: CAMPOS_VTK=1 grava campos.vtk, VTK_H=24 grava o mecanismo de ruptura em
# t = 24 h, e SO_CAMPOS=1 so gera os campos (campos.vtk e campos.csv) das realizacoes pedidas, sem simular.
set -u
PKG="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
[ $# -eq 2 ] || { sed -n '2,13p' "$0"; exit 2; }
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
echo "$(date '+%F %T') $COMB: realizacoes $1..$2, $P ao mesmo tempo, $(hostname -s)" | tee -a "$D/execucoes.txt"

# uma realizacao, na sua pasta; a saida do programa vai para run.log/run.err (acrescentando, para as retomadas)
uma() {
    local n="$1" r
    r="$D/$(printf 'r%05d' "$n")"
    mkdir -p "$r" && cd "$r" || return 1
    if [ "${SO_CAMPOS:-0}" = 1 ]; then # so os campos desta realizacao, em VTK e CSV, sem simulacao
        REAL="$n" CAMPOS_VTK=1 CAMPOS_CSV=1 "$EXE" > campos.log 2>&1 && echo "real $n: campos em $r/campos.vtk"
        return 0
    fi
    [ -s resultado.txt ] && return 0
    [ -s checkpoint.bin ] && echo "$(date '+%T') real $n: retomando do checkpoint"
    REAL="$n" "$EXE" >> run.log 2>> run.err
    [ -s resultado.txt ] || { echo "real $n: sem resultado.txt (veja $r/run.err)" >&2; return 1; }
    echo "$(date '+%T') real $n: $(cat resultado.txt)"
}
export -f uma
seq "$1" "$2" | xargs -P "$P" -I{} bash -c 'uma {}'
[ "${SO_CAMPOS:-0}" = 1 ] || "$PKG/junta_mc.sh" "$D"
echo "$(date '+%F %T') $COMB: fim" | tee -a "$D/execucoes.txt"
