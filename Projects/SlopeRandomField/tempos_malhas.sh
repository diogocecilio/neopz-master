#!/bin/bash
# tempos_malhas.sh <saida> <malha> [<malha> ...]: o caso deterministico do exemplo do Ceron (README.md, "Exemplo do
# Ceron": rebaixamento do nivel 10 m ate o pe, 5 m, variante iii, fator em 0,01 h) em cada malha (o nome do arquivo em
# ../LearningMohrCoulomb2/cluster/malhas, sem .msh, ou um caminho para um .msh), com o GIM por bissecao (gim) e com
# arc-length (gim_arco). Uma rodada por vez, para os tempos serem comparaveis. No fim, a tabela (tambem em
# <saida>/tempos.md): malha, numero de equacoes, fator, valor, tempo do calculo do fator, iteracoes de Newton do fator,
# solucoes tentadas (passos da busca ou incrementos do arco) e tempo total da rodada (com o transiente e a montagem).
# Ambiente: FATORES (padrao "gim gim_arco"), FS_H (padrao 0.01), DT_H (padrao 0.01), SOLVER, EXE. No fermi, num job:
#   sbatch -J tempos -N 1 -n 1 -c 2 --time=04:00:00 -o tempos_%j.out --wrap "source ../LearningMohrCoulomb2/cluster/modulos.sh && ./tempos_malhas.sh \$DATA/tempos_malhas ceron_h5_a45 ceron_h5_a45_fina"
set -u
PKG="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
[ $# -ge 2 ] || { sed -n '2,/^set -u/{/^#/p}' "$0"; exit 2; }
SAIDA="$1"; shift
EXE="${EXE:-$PKG/bin/SlopeRandomField}"
[ -x "$EXE" ] || { echo "executavel $EXE nao encontrado (./compilar.sh)" >&2; exit 2; }
FATORES="${FATORES:-gim gim_arco}"
FS_H="${FS_H:-0.01}" DT_H="${DT_H:-0.01}"
mkdir -p "$SAIDA" || exit 1
SAIDA="$(cd "$SAIDA" && pwd)"
linhas=()
for m in "$@"; do
    case "$m" in
    */* | *.msh) msh="$m" ;;
    *) msh="$PKG/../LearningMohrCoulomb2/cluster/malhas/$m.msh" ;;
    esac
    [ -f "$msh" ] || { echo "malha $msh nao encontrada" >&2; continue; }
    msh="$(cd "$(dirname "$msh")" && pwd)/$(basename "$msh")"
    nome="$(basename "$msh" .msh)"
    for f in $FATORES; do
        d="$SAIDA/${nome}_$f"
        rm -rf "$d" && mkdir -p "$d" || exit 1
        echo "$(date '+%F %T') $nome $f" >&2
        (cd "$d" && env SOLO=f50 I_MM_H=0 VARIANTE=iii CAMPOS=0 MESH="$msh" PHI_DEG=30 HL=10 HR=10 XL=0 XR=25 XTOE=15 YTOE=5 \
            XCRISTA=10 YCRISTA=10 NIVEL_FINAL=5 T_REBAIX_H=0 FS_H="$FS_H" DT_H="$DT_H" FATOR="$f" SRF_VERSAO=tempos \
            MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 "$EXE" > run.log 2> run.err)
        if [ ! -s "$d/resultado.txt" ]; then
            linhas+=("| $nome | - | $f | falhou | | | | |")
            continue
        fi
        neq="$(grep -o 'neq=[0-9]*' "$d/realizacao.txt" | cut -d= -f2)"
        # o fator, o tempo do calculo do fator (ultima coluna de fs.csv) e as iteracoes e tentativas da busca
        read -r fs segf <<< "$(awk -F, 'NR > 1 { v = $2; s = $NF } END { print v, s }' "$d/fs.csv")"
        read -r it nt <<< "$(awk -F, 'NR > 1 { i += $7; n++ } END { print i + 0, n + 0 }' "$d/fs_tentativas.csv")"
        tot="$(grep -o 'seg=[0-9.]*' "$d/resultado.txt" | cut -d= -f2)"
        linhas+=("$(printf '| %s | %s | %s | %.4f | %.1f | %d | %d | %.1f |' "$nome" "$neq" "$f" "$fs" "$segf" "$it" "$nt" "$tot")")
    done
done
{
    echo "| malha | equações | fator | valor | tempo do fator (s) | iterações de Newton do fator | soluções tentadas | tempo total (s) |"
    echo "|---|---|---|---|---|---|---|---|"
    printf '%s\n' "${linhas[@]}"
} | tee "$SAIDA/tempos.md"
