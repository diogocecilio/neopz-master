#!/bin/bash
# la_cho.sbatch: os casos do Cho (2010, Exemplo 2; Tabelas 4 e 5 do relatorio_comparacao.pdf, r = -0,5 e -0,7) com a
# ANALISE LIMITE cinematica do Ceron (2025) no lugar do MEF, nos MESMOS campos (Projects/CeronCinematico/cho, Python).
# Uma etapa por job (ETAPA=...), um no cada; saidas em $SAIDA/h$H (padrao $DATA/cho_la/h0.025); retoma de onde parou.
# H = passo da grade dos campos na LA (padrao 0,025 m; com 0,1 o max phi no arco e subestimado e Gamma sai ~2,5 % baixo).
#   ETAPA=checa       python com numpy/scipy/matplotlib, o bin/Cho2010 com as bibliotecas (modulos.sh), o deterministico
#                     (Gamma 1,777, FS 1,203), o Psi de Nystrom na grade (cache, uma vez, ~10-20 min com 24 threads) e o
#                     tempo de 48 realizacoes com FS neste no. Rode primeiro (sozinho: os outros usam o cache).
#   ETAPA=direto R="r05 r00 r07"   (um ou mais r, em sequencia) MC direto LA ate CoV(P_f) = 5 %, no minimo as
#                     realizacoes do MEF (24 000 | 110 000: semente 2010, campo analitico = as mesmas realizacoes do MC
#                     direto do MEF, pareadas uma a uma), no maximo N_MAX (300 000); FS da LA (~13 s a mais por
#                     realizacao) so nas FS_ATE primeiras (padrao 2400, para a relacao FS_LA x FS_MEF) e onde
#                     Gamma < FS_GAMMA (1,3 no r05, 1,1 no r07: perto da falha). Depois junta com o MEF (SRM).
#                     R=r00: r = 0, ate CoV 5 % (~7000; sem MEF), so Gamma: a conferencia com o Ceron (P_f 5,50 %).
#   ETAPA=amostragem R="r05 r07"   (campo de Nystrom de 2000 modos, o das ISMC e do LS do MEF) Psi de Nystrom (cache);
#                     FORM LA M = 20 e 100 (diferencas centrais; FORM_2000=1 tambem M = 2000 a frente); line sampling LA
#                     (100 linhas, semente 7) nas direcoes do MEF (M20, gim2000) e da LA; ISMC LA nos centros do MEF
#                     (M20, M100, gim20) e da LA (M20), UMA semente (SEMENTES_ISMC, padrao 21: a mesma da ISMC do SRM e
#                     do GIM que entram na comparacao), N_ISMC (2500, cobre as k do SRM s21, 1281 a 2191) cada (as mesmas
#                     realizacoes das ISMC do MEF).
#   ETAPA=pareado     (LA e MEF no mesmo job) gera os xi de N_PAR realizacoes (padrao 240, semente SEMENTE_PAR=2010,
#                     r R_PAR=r05), Gamma e FS da LA (busca forte), o Cho2010 (bin/Cho2010) em cada xi com o fator das
#                     referencias, e junta (LA x MEF; com a semente 2010, tambem confere com o MC direto do MEF ja feito).
#   ETAPA=resumo      a tabela LA x MEF de tudo o que houver em $SAIDA/h$H.
# Uso (no fermi, nesta pasta):  J=$(ETAPA=checa sbatch --parsable la_cho.sbatch)   (o ETAPA vai pelo ambiente)
#   ETAPA=direto R="r05 r00 r07" sbatch --dependency=afterok:$J la_cho.sbatch
#   ETAPA=amostragem R="r05 r07" sbatch --dependency=afterok:$J la_cho.sbatch
#SBATCH --job-name=laCho
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=48G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/%x_%j.out
[ -r /etc/profile.d/scratch.sh ] && source /etc/profile.d/scratch.sh
cd "$SLURM_SUBMIT_DIR" || exit 1
# o MKL do bin/Cho2010 (libmkl_*.so.2): o mesmo modulo dos outros sbatch que rodam o Cho2010
source ../../LearningMohrCoulomb2/cluster/modulos.sh || echo "AVISO: modulos.sh falhou (o MEF do pareado pode nao rodar)"
ETAPA="${ETAPA:-checa}" R="${R:-r05}" H="${H:-0.025}"
SAIDA="${SAIDA:-${DATA:?}/cho_la}" P="${PARALELO:-24}"
RES="$SAIDA/h$H"
LA="$(cd ../../CeronCinematico/cho 2>/dev/null && pwd)"
[ -n "$LA" ] && [ -r "$LA/mc_la.py" ] || { echo "falta Projects/CeronCinematico/cho (mande a pasta para o fermi)" >&2; exit 1; }
BASE="${BASE:-$DATA/etapa3b/pontos/base_nys2000.txt}" EXE="${EXE:-$PWD/bin/Cho2010}"
CACHE="$SAIDA/cache/psi_nys2000_h$H.npy"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PZ_THREADS=0 PYTHONUNBUFFERED=1
mkdir -p "$RES"
rc() { case "$1" in r05) echo -0.5 ;; r07) echo -0.7 ;; r00) echo 0.0 ;; *) echo "R = r05, r07 ou r00" >&2; exit 1 ;; esac; }
data() { date '+%F %T'; }

# ---- python com numpy, scipy e matplotlib: $PY, ou o venv $HOME/venv_la, ou o python3 do sistema
if [ -z "$PY" ]; then
    for c in "$HOME/venv_la/bin/python3" python3; do
        if command -v "$c" >/dev/null 2>&1 && "$c" -c "import numpy, scipy, matplotlib" 2>/dev/null; then PY="$c"; break; fi
    done
fi
if [ -z "$PY" ]; then
    cat >&2 <<'EOF'
ERRO: nenhum python3 com numpy, scipy e matplotlib (o python3 do sistema nao tem; o venv $HOME/venv_la nao existe).
No no de login do fermi (precisa de internet; uma vez so):
    python3 -m venv $HOME/venv_la && $HOME/venv_la/bin/pip install numpy scipy matplotlib
(se o venv falhar: python3 -m pip install --user numpy scipy matplotlib). Depois submeta de novo.
Ou indique outro python: PY=/caminho/python3 ETAPA=... sbatch la_cho.sbatch
EOF
    exit 2
fi
echo "$(data) job ${SLURM_JOB_ID:-local} em $(hostname -s): ETAPA=$ETAPA R=$R h=$H, $P processos, saidas em $RES, python $PY"
"$PY" -c "import sys, numpy, scipy, matplotlib; print('  python', sys.version.split()[0], 'numpy', numpy.__version__, 'scipy', scipy.__version__)"

case "$ETAPA" in
checa)
    if [ -x "$EXE" ]; then
        if ldd "$EXE" | grep -q "not found"; then
            echo "ERRO: $EXE sem bibliotecas (o pareado vai falhar):" >&2; ldd "$EXE" | grep "not found" >&2; exit 1
        fi
        echo "  $EXE: bibliotecas ok"
    else
        echo "  AVISO: falta $EXE (so o pareado precisa; compilar.sh)"
    fi
    # os z do line sampling (random.Random(7).gauss) tem de ser os do PC e os do line_sampling.py do MEF
    z="$("$PY" -I -c "import random;r=random.Random(7);print(repr(r.gauss(0.,1.)),repr(r.gauss(0.,1.)))")"
    [ "$z" = "-0.2558802884476004 0.511431512516514" ] || { echo "ERRO: random.gauss diferente do PC: $z" >&2; exit 1; }
    echo "  random.Random(7).gauss: ok ($z)"
    "$PY" "$LA/testes.py" det || exit 1
    # os dados do MEF que a comparacao usa (so avisa: o pareado tambem pode ser feito no PC)
    for d in "$DATA/etapa1/r05_direto" "$DATA/etapa3/r07_direto" "$DATA"/etapa3b/r0?_nys_*_s21 "$DATA"/etapa3b_gim/r0?_nys_*_s21; do
        if [ -d "$d" ]; then
            n=$(cat "$d"/mc_[0-9]*_[0-9]*.csv 2>/dev/null | grep -vc '^real')
            echo "  MEF $d: $n realizacoes"
        else
            echo "  AVISO: falta $d (o pareado com ele fica para o PC)"
        fi
    done
    if [ -r "$BASE" ]; then echo "  base de Nystrom: $BASE ($(du -h "$BASE" | cut -f1))"; else echo "  AVISO: falta $BASE (a etapa amostragem precisa dela)"; fi
    for f in form_r05_M20 form_r05_gim2000 form_r07_M20 form_r07_gim2000; do
        [ -r "$LA/form_mef/$f/form.json" ] || echo "  AVISO: falta $LA/form_mef/$f/form.json"
    done
    if [ -r "$BASE" ] && [ ! -f "$CACHE" ]; then
        echo "$(data) Psi de Nystrom na grade h = $H (uma vez)"
        "$PY" "$LA/campos_cho.py" cache base="$BASE" cache="$CACHE" h="$H" procs="$P" || exit 1
    fi
    rm -f "$RES/checa/tempo.csv"
    "$PY" "$LA/mc_la.py" saida="$RES/checa/tempo.csv" r=-0.5 n=48 cov=0 procs="$P" bloco=1 fs_ate=48 h="$H" || exit 1
    ;;
direto)
  for R in $R; do
    echo "$(data) ===== direto $R"
    # r00 (r = 0, sem MEF): so para conferir a LA contra o Ceron (2025, Sec. 5.3.2: P_f 5,50 % com S = 50 000)
    case "$R" in
    r05) NMIN=24000 FSG=1.3 FSA=2400 MEFD0="$DATA/etapa1/r05_direto" ;;
    r07) NMIN=110000 FSG=1.1 FSA=2400 MEFD0="$DATA/etapa3/r07_direto" ;;
    r00) NMIN=2000 FSG=0 FSA=0 MEFD0="" ;;
    *) echo "R = r05, r07 ou r00" >&2; exit 1 ;;
    esac
    MEFD="$MEFD0"; [ "$R" != r00 ] && MEFD="${MEF_DIRETO:-$MEFD0}"   # r00 nao tem MEF
    "$PY" "$LA/mc_la.py" saida="$RES/direto/$R.csv" r="$(rc "$R")" campo=analitico semente=2010 n="${N_MAX:-300000}" \
        cov=0.05 nmin="$NMIN" procs="$P" fs_ate="${FS_ATE:-$FSA}" fs_gamma="${FS_GAMMA:-$FSG}" h="$H" || exit 1
    if [ -n "$MEFD" ] && [ -d "$MEFD" ]; then
        "$PY" "$LA/pareado.py" juntar la="$RES/direto/$R.csv" mef="$MEFD" saida="$RES/direto/${R}_pares.csv" || exit 1
    else
        echo "$(data) sem $MEFD (MC direto do MEF): pareie no PC com pareado.py juntar"
    fi
    # o direto do MEF pelo GIM (Gamma_GIM, campanha etapa3b_gim; 24 000 realizacoes), se ja existir: Gamma_LA x Gamma_GIM
    MEFG="${MEF_GIM:-$DATA/etapa3b_gim/${R}_direto_gim}"
    if [ "$R" != r00 ] && [ -d "$MEFG" ]; then
        echo "$(data) LA x MEF GIM ($MEFG): a coluna FS_MEF do pares e o Gamma do GIM"
        "$PY" "$LA/pareado.py" juntar la="$RES/direto/$R.csv" mef="$MEFG" saida="$RES/direto/${R}_pares_gim.csv" || true
    fi
  done
    ;;
amostragem)
  for R in $R; do
    echo "$(data) ===== amostragem $R"
    [ -r "$BASE" ] || { echo "falta $BASE" >&2; exit 1; }
    rcv="$(rc "$R")" NI="${N_ISMC:-2500}"
    mkdir -p "$RES/ismc"
    if [ ! -f "$CACHE" ]; then
        echo "$(data) Psi de Nystrom na grade h = $H (uma vez; melhor no checa, antes dos dois amostragem)"
        "$PY" "$LA/campos_cho.py" cache base="$BASE" cache="$CACHE" h="$H" procs="$P" || exit 1
    fi
    NY=(campo=nystrom base="$BASE" cache="$CACHE" h="$H")
    MS="20 100"; [ "${FORM_2000:-0}" = 1 ] && MS="$MS 2000"
    for M in $MS; do
        F="$RES/form/${R}_M$M"
        [ -f "$F/form.json" ] && continue
        dif=central; [ "$M" -gt 200 ] && dif=frente
        echo "$(data) FORM LA $R M = $M"
        "$PY" "$LA/form_la.py" pasta="$F" r="$rcv" M="$M" dif=$dif procs="$P" "${NY[@]}" || exit 1
    done
    for d in mef_M20 mef_gim2000 la_M20 la_M100 la_M2000; do
        case $d in mef_*) F="$LA/form_mef/form_${R}_${d#mef_}/form.json" ;; la_*) F="$RES/form/${R}_${d#la_}/form.json" ;; esac
        [ -r "$F" ] || continue
        [ -f "$RES/ls/${R}_$d/ls.json" ] && continue
        echo "$(data) line sampling LA $R direcao $d"
        "$PY" "$LA/ls_la.py" direcao="$F" pasta="$RES/ls/${R}_$d" r="$rcv" n="${N_LS:-100}" semente=7 K=2000 procs="$P" \
            "${NY[@]}" || exit 1
    done
    for c in mef_M20 mef_gim20 mef_M100 la_M20; do
        case $c in mef_*) F="$LA/form_mef/form_${R}_${c#mef_}/form.json" ;; la_*) F="$RES/form/${R}_${c#la_}/form.json" ;; esac
        [ -r "$F" ] || continue
        for s in ${SEMENTES_ISMC:-21}; do
            echo "$(data) ISMC LA $R centro $c semente $s ($NI)"
            fsg=0; [ "$s" = 21 ] && fsg="${FS_GAMMA:-1.3}"   # FS (Gamma < 1,3) so na semente 21: a pareada com o MEF
            "$PY" "$LA/mc_la.py" saida="$RES/ismc/${R}_${c}_s$s.csv" r="$rcv" semente=$s n="$NI" cov=0 procs="$P" \
                centro="$F" fs_gamma="$fsg" "${NY[@]}" > "$RES/ismc/${R}_${c}_s$s.log" 2>&1 \
                || { tail -n 30 "$RES/ismc/${R}_${c}_s$s.log" >&2; exit 1; }
            tail -n 9 "$RES/ismc/${R}_${c}_s$s.log"
        done
    done
  done
    ;;
pareado)
    RP="${R_PAR:-r05}" SP="${SEMENTE_PAR:-2010}" NP="${N_PAR:-240}" PP="${PRIMEIRA_PAR:-0}"
    D="$RES/pareado/${RP}_s${SP}_${PP}_$NP"
    [ -x "$EXE" ] || { echo "falta $EXE (compilar.sh)" >&2; exit 1; }
    if ldd "$EXE" | grep -q "not found"; then
        echo "ERRO: $EXE sem bibliotecas (modulos.sh / mkl):" >&2; ldd "$EXE" | grep "not found" >&2; exit 1
    fi
    if [ ! -f "$D/parametros.json" ]; then
        "$PY" "$LA/pareado.py" gerar pasta="$D" r="$(rc "$RP")" semente="$SP" primeira="$PP" n="$NP" h="$H" || exit 1
    fi
    mkdir -p "$D"
    "$PY" "$LA/pareado.py" la pasta="$D" procs="$P" > "$D/la.log" 2>&1 || { tail -n 30 "$D/la.log" >&2; exit 1; }
    tail -n 9 "$D/la.log"
    "$PY" "$LA/pareado.py" mef pasta="$D" exe="$EXE" procs="$P" > "$D/mef.log" 2>&1 || { tail -n 40 "$D/mef.log" >&2; exit 1; }
    tail -n 3 "$D/mef.log"
    "$PY" "$LA/pareado.py" juntar pasta="$D" || exit 1
    MEFD="$([ "$RP" = r05 ] && echo "$DATA/etapa1/r05_direto" || echo "$DATA/etapa3/r07_direto")"
    if [ "$SP" = 2010 ] && [ -d "$MEFD" ]; then
        echo "$(data) as mesmas realizacoes no MC direto do MEF ($MEFD): o FS pelo xi= deve ser identico"
        "$PY" - "$D/mef.csv" "$MEFD" <<'EOF'
import csv, glob, os, sys
a = {int(l["real"]): float(l["FS_MEF"]) for l in csv.DictReader(open(sys.argv[1])) if l["situacao"] == "ok"}
b = {}
for f in glob.glob(os.path.join(sys.argv[2], "mc_[0-9]*_[0-9]*.csv")):
    for l in csv.DictReader(open(f)):
        if int(l["real"]) in a:
            b[int(l["real"])] = float(l["FS_min"])
d = [abs(a[k] - b[k]) for k in b]
print(f"  {len(d)} realizacoes em comum; |FS(xi=) - FS(direto)| max = {max(d) if d else float('nan'):.2e}")
sys.exit(0 if d and max(d) < 1e-6 else 4)
EOF
        [ $? -eq 0 ] || { echo "ERRO: o FS do MEF pelo xi= nao confere com o do MC direto (ou nenhuma em comum)" >&2; exit 1; }
    fi
    ;;
resumo) ;;
*) echo "ETAPA desconhecida: $ETAPA (checa, direto, amostragem, pareado, resumo)" >&2; exit 1 ;;
esac
echo "$(data) resumo de $RES:"
"$PY" "$LA/resumo.py" "$RES"
echo "$(data) fim"
