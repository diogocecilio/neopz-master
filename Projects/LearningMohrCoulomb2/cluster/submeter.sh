#!/bin/bash
# submeter.sh [casos]: submete a compilacao (se preciso) e, depois dela, os casos de casos.txt (padrao 1-21; ex.: 3,7 ou 1-9).
# Retorna em segundos; pode desligar a VPN logo em seguida. Os jobs seguem no cluster.
# Os casos vao num unico job (job_casos.sbatch), que conta como 1 job na cota do SLURM e roda ate NCPU casos ao mesmo tempo,
# 1 CPU cada (no fermi, um array de 21 tarefas excede a cota de jobs submetidos por usuario). MODO=array usa job_slurm.sbatch.
# Compilacao: COMPILAR=login (padrao no fermi, em modulos.sh) compila neste no em segundo plano e submete os casos ao fim;
# COMPILAR=job compila num job do SLURM (job_compilar.sbatch), do qual os casos dependem.
# Ambiente opcional: NCPU (padrao: numero de casos, limitado aos nucleos de um no); PARALELO (casos ao mesmo tempo, padrao
# NCPU; com hyper-threading, NCPU=2*PARALELO deixa um nucleo fisico por caso); SAIDA (pasta de resultados, padrao
# resultados); TF_H (horas simuladas, padrao 48); TEMPO (limite do job de casos, padrao
# 1-00:00:00; os casos sao retomaveis); RECOMPILAR=1; EMAIL=<endereco> (aviso do SLURM ao terminar ou falhar, se o cluster
# enviar e-mail); FORCAR=1 (submete mesmo com outro job de casos na fila). Os numeros dos jobs ficam em logs/jobs.txt.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
source ./modulos.sh
EU="${USER:-$(id -un)}"
P=(); [ -n "${PARTICAO:-}" ] && P=(--partition="$PARTICAO")
MB=(); MA=()
if [ -n "${EMAIL:-}" ]; then MB=(--mail-user="$EMAIL" --mail-type=FAIL); MA=(--mail-user="$EMAIL" --mail-type=END,FAIL); fi
LISTA="${1:-1-21}"
NCASOS="$(echo "$LISTA" | tr ',' '\n' | awk -F- 'NF == 2 { n += $2 - $1 + 1; next } NF { n++ } END { print n + 0 }')"
# dois jobs com os mesmos casos gravariam nas mesmas pastas
J="$(squeue -h -u "$EU" -n rahardjo45 -o %i 2>/dev/null | head -n 1)" || true
if [ -n "$J" ] && [ -z "${FORCAR:-}" ]; then
  echo "ja ha um job de casos na fila (job $J); veja ./status.sh. Para submeter mesmo assim: FORCAR=1 ./submeter.sh $LISTA" >&2
  exit 1
fi

# o CMake baixa o Catch2 ao configurar e os nos de computacao podem nao ter internet: baixa aqui, no no de login
baixar_catch2() {
  [ -f deps/Catch2/CMakeLists.txt ] && return 0
  rm -rf deps/Catch2
  git -c advice.detachedHead=false clone -q --depth 1 --branch v2.13.10 https://github.com/catchorg/Catch2.git deps/Catch2
  echo "Catch2 v2.13.10 baixado em deps/Catch2"
}
# compilacao em segundo plano neste no (COMPILAR=login): logs/compilar_login.pid guarda "pid no inicio"
compilando_aqui() {
  local pid no
  [ -f logs/compilar_login.pid ] || return 1
  read -r pid no _ < logs/compilar_login.pid || return 1
  [ "$no" = "$(hostname -s)" ] && kill -0 "$pid" 2>/dev/null
}

B=""
if [ -x bin/rahardjo ] && [ -z "${RECOMPILAR:-}" ]; then
  echo "bin/rahardjo ja existe; pulando a compilacao (RECOMPILAR=1 para refazer)"
elif [ "${COMPILAR:-job}" = login ]; then
  if compilando_aqui; then
    echo "a compilacao ja esta rodando neste no (pid $(cut -d' ' -f1 logs/compilar_login.pid)); os casos serao submetidos ao fim"
    exit 0
  fi
  baixar_catch2
  # nohup + setsid: continua depois que voce desconectar; ao fim, submete os casos (RECOMPILAR vazio: nao compila de novo)
  nohup setsid bash -c 'nice -n 10 ./compilar.sh && RECOMPILAR= ./submeter.sh "$1"' _ "$LISTA" \
    > logs/compilar_login.log 2>&1 < /dev/null &
  echo "$! $(hostname -s) $(date '+%F %T')" > logs/compilar_login.pid
  echo "compilando neste no de login em segundo plano (pid $!, log em logs/compilar_login.log, ~10 min com JOBS=${JOBS:-4})."
  echo "Ao terminar, os casos $LISTA sao submetidos sozinhos. Pode desconectar; acompanhe com ./status.sh"
  exit 0
elif B="$(squeue -h -u "$EU" -n rahardjo45_build -o %i 2>/dev/null | head -n 1)" && [ -n "$B" ]; then
  echo "compilacao ja na fila: job $B; os casos esperam por ela"
else
  baixar_catch2
  B=$(sbatch --parsable ${P[@]+"${P[@]}"} ${MB[@]+"${MB[@]}"} job_compilar.sbatch); B="${B%%;*}"
  echo "compilacao: job $B"
fi
DEP=(); [ -n "$B" ] && DEP=(--dependency=afterok:"$B" --kill-on-invalid-dep=yes)

if [ "${MODO:-pool}" = array ]; then
  A=$(sbatch --parsable ${P[@]+"${P[@]}"} ${MA[@]+"${MA[@]}"} ${DEP[@]+"${DEP[@]}"} --array="$LISTA" job_slurm.sbatch); A="${A%%;*}"
  echo "casos: job array $A ($LISTA)"
else
  # o job de casos fica num no so: limita as CPUs aos nucleos e a memoria a 90% da memoria do menor no da particao
  read -r CPN MEMN <<< "$(sinfo -h -N ${P[@]+"${P[@]}"} -o '%c %m' 2>/dev/null | tr -d '+' | sort -n | head -n 1)" || true
  case "${CPN:-}" in ''|*[!0-9]*) CPN="" ;; esac
  case "${MEMN:-}" in ''|*[!0-9]*) MEMN="" ;; esac
  K="${NCPU:-$NCASOS}"
  [ "$K" -gt "$NCASOS" ] && K="$NCASOS"
  [ -n "${CPN:-}" ] && [ "$K" -gt "$CPN" ] && K="$CPN"
  MEM=$((K * 1500))
  [ -n "${MEMN:-}" ] && [ "$MEM" -gt $((MEMN * 9 / 10)) ] && MEM=$((MEMN * 9 / 10))
  # a lista vai pelo ambiente (o sbatch exporta o ambiente todo); em --export a virgula de "3,7" separaria variaveis
  A=$(CASOS="$LISTA" sbatch --parsable ${P[@]+"${P[@]}"} ${MA[@]+"${MA[@]}"} ${DEP[@]+"${DEP[@]}"} --cpus-per-task="$K" \
      --mem="${MEM}M" --time="${TEMPO:-1-00:00:00}" job_casos.sbatch); A="${A%%;*}"
  echo "casos: job $A ($NCASOS casos, ate ${PARALELO:-$K} ao mesmo tempo em $K CPUs; no de ${CPN:-?} CPUs e ${MEMN:-?} MB)"
fi
echo "$(date '+%F %T') compilacao=${B:--} casos=$A lista=$LISTA" >> logs/jobs.txt
echo "acompanhe depois com: ./status.sh   (ou ./status.sh --slurm para ver tambem duracao e memoria no sacct)"
