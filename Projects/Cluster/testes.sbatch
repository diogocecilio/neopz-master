#!/bin/bash
# testes.sbatch: os testes dos projetos de verificacao (padrao: todos; ou PROJETOS="V1 V3 ...") num no de computacao.
# Primeiro os padrao, pelo ctest em build/Projects (a mesma coisa que rodar "ctest" la: cada executavel de testes sem
# argumentos, que exclui os [.lento]), todos ao mesmo tempo; depois os demorados de cada projeto ("<teste> [.lento]"), um
# processo de 1 CPU por projeto, todos ao mesmo tempo. Os nos de computacao do fermi nao tem cmake: sem o ctest (no
# PATH, ou no caminho que compilar.sh anota em build/ctest_caminho), os testes registrados no CTest (os add_test dos
# CTestTestfile.cmake) rodam direto, como o ctest os rodaria (cada um na sua pasta de build), e o log diz isso.
# Logs em logs/testes_<job>/: resumo.txt (passou ou nao, tempo e a ultima linha do Catch2 de cada grupo; tambem na saida
# do job, logs/verif-testes_<job>.out), ctest.log (ou ctest_<teste>.log) e <projeto>_lento.log. Os testes demorados
# rodam em pastas de trabalho proprias (alguns gravam arquivos na pasta corrente): $DATA/verificacao/testes/<job>/
# <projeto> (sem $DATA, logs/testes_<job>/<projeto>). Codigo de saida 0 so se tudo passou.
# Compile antes, no no de login (./compilar.sh; ./submeter_testes.sh compila se preciso e submete). A mao, de dentro de
# Projects/Cluster: sbatch testes.sbatch, ou PROJETOS="V2 V4" LENTOS=0 sbatch testes.sbatch.
# Ambiente: PROJETOS, LENTOS=0 (so os padrao), BUILD (padrao ./build), TIMEOUT_TESTE (s, limite de cada teste do
# ctest; padrao 14400) e as variaveis dos proprios testes (por exemplo V2_NMC, V4_N_MC: ver o README de cada projeto).
#SBATCH --job-name=verif-testes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=24G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%j.out
[ -r /etc/profile.d/scratch.sh ] && source /etc/profile.d/scratch.sh
PKG=""
for p in "${VERIF_CLUSTER:-}" "${SLURM_SUBMIT_DIR:-}" "${SLURM_SUBMIT_DIR:-}/.."; do
    if [ -n "$p" ] && [ -f "$p/comum.sh" ] && [ -f "$p/testes.sbatch" ]; then PKG="$(cd "$p" && pwd)"; break; fi
done
[ -n "$PKG" ] || { echo "nao achei Projects/Cluster: submeta de la (ou com submeter_testes.sh)" >&2; exit 1; }
cd "$PKG" || exit 1
source "$PKG/modulos.sh"
source "$PKG/comum.sh"
export MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 LC_ALL=C
BUILD="${BUILD:-$PKG/build}"
JOB="${SLURM_JOB_ID:-local$$}"
L="$PKG/logs/testes_$JOB"
mkdir -p "$L"
TRAB="${DATA:+$DATA/verificacao/testes/$JOB}"
TRAB="${TRAB:-$L}"
P="$(paralelo_job)"
PS=()
if [ -n "${PROJETOS:-}" ]; then
    for a in $PROJETOS; do p="$(projeto_completo "$a")" || exit 2; PS+=("$p"); done
else
    mapfile -t PS < <(projetos_verif)
fi
resumo() { echo "$*" | tee -a "$L/resumo.txt"; }
ok_tudo=1
resumo "$(date '+%F %T') testes em $(hostname -s) (job $JOB): ${PS[*]}; ate $P processos; build $BUILD"

# os executaveis de teste, com as bibliotecas; fontes mais novos que eles: um aviso (os testes sao os compilados)
OK=()
for p in "${PS[@]}"; do
    t="$BUILD/Projects/$p/$(alvo_teste "$p")"
    f="$(aviso_ldd "$t")"
    if [ -n "$f" ]; then
        resumo "$p: FALHOU antes de rodar: $f"
        ok_tudo=0
        continue
    fi
    n="$(fonte_nova "$p" "$t")"
    [ -n "$n" ] && resumo "aviso: $p: fontes mais novos que o teste compilado (${n#../../}): rode ./compilar.sh"
    OK+=("$p")
done

# ---- os testes padrao
if [ ${#OK[@]} -gt 0 ]; then
    RE="^($(for p in "${OK[@]}"; do alvo_teste "$p"; done | paste -sd'|' -))\$"
    CT=""
    for c in "$(command -v ctest 2>/dev/null)" "$(cat "$BUILD/ctest_caminho" 2>/dev/null)"; do
        if [ -n "$c" ] && "$c" --version > /dev/null 2>&1; then CT="$c"; break; fi
    done
    t0=$(date +%s)
    if [ -n "$CT" ]; then
        (cd "$BUILD/Projects" && "$CT" --output-on-failure -j "$P" --timeout "${TIMEOUT_TESTE:-14400}" -R "$RE") \
            > "$L/ctest.log" 2>&1
        rc=$?
        lin="$(grep -E 'tests passed|tests failed' "$L/ctest.log" | tail -n 1)"
        resumo "ctest ($CT): ${lin:-sem o resumo do ctest} ($(($(date +%s) - t0)) s; log $L/ctest.log)"
        sed -n '/The following tests FAILED/,$p' "$L/ctest.log" | sed 's/^/  /' | tee -a "$L/resumo.txt"
        [ $rc -eq 0 ] || ok_tudo=0
    else
        # sem ctest neste no: os testes registrados (add_test), como o ctest os rodaria
        resumo "sem ctest neste no (nem em build/ctest_caminho): os testes registrados no CTest rodam direto"
        nt=0 np=0
        while IFS=$'\t' read -r d nome exe; do
            [[ "$nome" =~ $RE ]] || continue
            while [ "$(jobs -rp | wc -l)" -ge "$P" ]; do wait -n 2>/dev/null || sleep 2; done
            (s0=$(date +%s); cd "$d" && timeout "${TIMEOUT_TESTE:-14400}" "$exe" > "$L/ctest_$nome.log" 2>&1
             echo "$? $(($(date +%s) - s0))" > "$L/.rc_$nome") &
        done < <(for f in "$BUILD"/Projects/*/CTestTestfile.cmake; do
            [ -f "$f" ] || continue
            sed -n 's/^add_test(\[*=*\[*\([A-Za-z0-9_.+-]*\)\]*=*\]* "\([^"]*\)".*/\1\t\2/p' "$f" |
                while IFS=$'\t' read -r nome exe; do printf '%s\t%s\t%s\n' "$(dirname "$f")" "$nome" "$exe"; done
        done)
        wait
        for f in "$L"/.rc_*; do
            [ -e "$f" ] || continue
            nome="${f##*/.rc_}"
            read -r rc s < "$f"
            rm -f "$f"
            nt=$((nt + 1))
            u="$(grep -E '^(All tests passed|test cases:)' "$L/ctest_$nome.log" | tail -n 1)"
            if [ "$rc" = 0 ]; then np=$((np + 1)); else ok_tudo=0; fi
            resumo "  $nome: $([ "$rc" = 0 ] && echo passou || echo "FALHOU (codigo $rc)") em $s s: ${u:-sem o resumo do Catch2}"
        done
        resumo "testes registrados: $np de $nt passaram ($(($(date +%s) - t0)) s; logs $L/ctest_<teste>.log)"
        [ "$nt" -gt 0 ] || ok_tudo=0
    fi
fi

# ---- os testes demorados ([.lento]) de cada projeto, todos ao mesmo tempo
if [ "${LENTOS:-1}" = 1 ]; then
    for p in ${OK[@]+"${OK[@]}"}; do
        t="$BUILD/Projects/$p/$(alvo_teste "$p")"
        while [ "$(jobs -rp | wc -l)" -ge "$P" ]; do wait -n 2>/dev/null || sleep 5; done
        mkdir -p "$TRAB/$p"
        (s0=$(date +%s); cd "$TRAB/$p" && "$t" "[.lento]" -d yes > "$L/${p}_lento.log" 2>&1
         echo "$? $(($(date +%s) - s0))" > "$L/.lento_$p") &
    done
    wait
    for p in ${OK[@]+"${OK[@]}"}; do
        f="$L/.lento_$p"
        [ -e "$f" ] || { resumo "$p [.lento]: FALHOU (o processo nao terminou)"; ok_tudo=0; continue; }
        read -r rc s < "$f"
        rm -f "$f"
        g="$L/${p}_lento.log"
        u="$(grep -E '^(All tests passed|test cases:)' "$g" | tail -n 1)"
        if grep -q '^No test cases matched' "$g"; then
            resumo "$p [.lento]: nenhum teste com essa marca"
        elif [ "$rc" = 0 ]; then
            resumo "$p [.lento]: passou em $s s: ${u:-sem o resumo do Catch2}"
        else
            resumo "$p [.lento]: FALHOU (codigo $rc) em $s s: ${u:-sem o resumo do Catch2} (log $g)"
            ok_tudo=0
        fi
    done
fi
resumo "$(date '+%F %T') fim: $([ $ok_tudo = 1 ] && echo 'tudo passou' || echo 'HA FALHAS (veja acima)')"
[ $ok_tudo = 1 ]
