#!/bin/bash
# submeter.sh [nos = 4]: submete <nos> jobs job.sbatch (um no inteiro cada), que dividem a fila de tarefas dos CASOS da
# campanha. Volta em segundos: pode desconectar. Pode rodar de novo a qualquer momento: os jobs novos pegam as tarefas
# livres. Sem jobs da campanha na fila, libera antes as tarefas interrompidas (as travas sem mapa), que sao refeitas.
# Ambiente:
#   CAMPANHA  chuva (padrao: casos.txt, jobs chuvaMC) ou ceron (casos_ceron.txt, jobs ceronMC); veja campanha.sh
#   SAIDA     pasta dos resultados (padrao $DATA/chuva_mc ou $DATA/ceron_ai, na particao de dados do fermi; sem
#             $DATA, ./resultados)
#   CASOS     todos (padrao), ids ou fases separados por virgula: CASOS=fase1, CASOS=ref,I=36
#   IREF      intensidade de referencia em mm/h (chuva; padrao 31.05); fica gravada em $SAIDA/campanha.txt e nao pode
#             mudar
#   BLOCO     realizacoes por tarefa (padrao 20, umas 2 h no fermi); tambem fixo por SAIDA
#   ULTIMA    ultima realizacao de todos os casos escolhidos (ex.: ULTIMA=9999 CASOS=ref acrescenta realizacoes)
#   PARALELO  processos por no (padrao 12, os nucleos fisicos; 24 usa o hyperthreading)
#   FORM_M    os FORM do Cho2010 que o primeiro job roda antes da fila (form.sh; padrao: FORM_PADRAO da campanha, ex.
#             "10 20" no etapa1; os ja prontos sao pulados; FORM_M= vazio nao roda nenhum)
#   TEMPO     limite de cada job (padrao 2-00:00:00), PARTICAO
#   DEPENDE   id de um job: os jobs novos so comecam quando ele acabar (--dependency=afterany; ficam na fila como
#             submetidos, sem ocupar vaga de job rodando)
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
source ../../LearningMohrCoulomb2/cluster/modulos.sh
source ./campanha.sh
NOS="${1:-4}"
case "$NOS" in '' | *[!0-9]* | 0) echo "uso: $0 [numero de nos, padrao 4]" >&2; exit 2 ;; esac
SAIDA="${SAIDA:-${DATA:+$DATA/$PASTA}}"
export SAIDA="${SAIDA:-$PWD/resultados_$CAMPANHA}" CASOS="${CASOS:-todos}" IREF="${IREF:-$IREF_PADRAO}"
export BLOCO="${BLOCO:-$BLOCO_PADRAO}"
export PARALELO="${PARALELO:-12}"
FORM_M="${FORM_M-$FORM_PADRAO}"
[ -n "${ULTIMA:-}" ] && export ULTIMA
export BIN="${BIN:-bin}" # BIN=: o bin dos jobs (padrao bin/; ex.: BIN=bin_frag, do compilar_frag.sh)
[ -x "$BIN/$PROG" ] || { echo "compile antes, no no de login: ./compilar.sh (ou o compilar do BIN=$BIN)" >&2; exit 1; }

# os casos pedidos existem?
if [ "$CASOS" != todos ]; then
    for c in ${CASOS//,/ }; do
        [[ "$c" =~ ^fase[0-9]+$ ]] && continue
        awk -v c="$c" '$1 == c { achou = 1 } END { exit !achou }' "$CASOS_TXT" || { echo "caso desconhecido: $c" >&2; exit 2; }
    done
fi

# IREF e BLOCO definem os arquivos de cada caso: fixos durante a campanha
mkdir -p logs "$SAIDA"
if [ -f "$SAIDA/campanha.txt" ]; then
    read -r i0 b0 < "$SAIDA/campanha.txt"
    if [ "$i0" != "$IREF" ] || [ "$b0" != "$BLOCO" ]; then
        echo "$SAIDA e da campanha IREF=$i0 BLOCO=$b0: use os mesmos valores (ou outra SAIDA)" >&2
        exit 1
    fi
else
    echo "$IREF $BLOCO" > "$SAIDA/campanha.txt"
fi

fila="$(squeue -h -u "${USER:-$(id -un)}" -n "$NOME" -o %i)" || { echo "o squeue nao respondeu: tente de novo" >&2; exit 1; }
if [ -z "$fila" ]; then
    # nenhum job na fila: uma pasta .tarefas que nao se lista (erro de E/S do sistema de arquivos: frag/f50_i em
    # 25/11/2026; criar e ler cada trava ainda funcionava) esconderia do find abaixo as travas das tarefas interrompidas,
    # que nunca seriam refeitas: vai para o lado (<pasta>_eio_<data>), e a fila cria uma nova; as tarefas sem mapa
    # serao refeitas (as com mapa a fila ja pula pelo mapa)
    for d in "$SAIDA"/*/.tarefas; do
        [ -d "$d" ] || continue
        ls "$d" > /dev/null 2>&1 && continue
        mv "$d" "${d}_eio_$(date +%Y%m%d%H%M%S)" &&
            echo "$d nao se lista (erro de E/S): trocada por uma nova; as tarefas sem mapa serao refeitas"
    done
    # uma trava sem o mapa e de uma tarefa interrompida (fim do tempo, queda do no): libera
    n=0
    while IFS= read -r -d '' trava; do
        t="$(basename "$trava")" caso="$(dirname "$(dirname "$trava")")"
        [ -f "$caso/mc_mapa_$t.csv" ] || { rm -rf "$trava"; n=$((n + 1)); }
    done < <(find "$SAIDA" -path '*/.tarefas/*' -mindepth 1 -print0 2> /dev/null)
    [ "$n" -eq 0 ] || echo "liberadas $n tarefas interrompidas (serao refeitas)"
    # um FORM sem o ponto de projeto, com a trava (job interrompido) ou com ERRO: refeito por um job novo
    for t in "$SAIDA"/form_*/trava; do
        [ -e "$t" ] || continue
        d="$(dirname "$t")"
        [ -s "$d/ustar.txt" ] || { rm -f "$t" "$d/ustar.txt.ERRO"; echo "o FORM $(basename "$d") sera refeito"; }
    done
else
    echo "ha $(echo "$fila" | wc -l) job(s) $NOME na fila: os novos entram na mesma fila de tarefas"
fi
rm -f "$SAIDA/PARAR"
for ((no = 0; no < NOS; no++)); do
    fm=""
    ((no == 0)) && fm="$FORM_M"
    j="$(FORM_M="$fm" sbatch --parsable --job-name="$NOME" ${PARTICAO:+--partition="$PARTICAO"} ${TEMPO:+--time="$TEMPO"} \
        ${DEPENDE:+--dependency=afterany:"$DEPENDE"} job.sbatch)" || {
        echo "o SLURM recusou o job $no (cota?): os ja submetidos seguem; submeta o resto depois" >&2
        break
    }
    echo "job ${j%%;*} ($((no + 1)) de $NOS)${fm:+: antes da fila, FORM com M = $fm}"
    echo "$(date '+%F %T') ${j%%;*} CAMPANHA=$CAMPANHA PARALELO=$PARALELO CASOS=$CASOS IREF=$IREF BLOCO=$BLOCO${ULTIMA:+ ULTIMA=$ULTIMA}${fm:+ FORM_M=\"$fm\"} SAIDA=$SAIDA" \
        >> logs/jobs.txt
done
echo "resultados em $SAIDA; acompanhe com CAMPANHA=$CAMPANHA ./status.sh"
