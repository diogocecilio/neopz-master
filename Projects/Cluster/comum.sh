# comum.sh: funcoes comuns dos scripts de Projects/Cluster (compilar.sh, os submeter_*.sh, os jobs e os status).
# Uso: source "<Projects/Cluster>/comum.sh". Define:
#   VERIF_PKG    a pasta Projects/Cluster          VERIF_NEOPZ  a raiz do NeoPZ (dois niveis acima)
#   VERIF_MAQ    o nome curto desta maquina        VERIF_BOOT   o boot_id (um pid de antes de um reinicio nao conta)
#   VERIF_EU     o usuario
# Os jobs deste pacote tem nomes verif-*: verif-mc-<V2|V4> (Monte Carlo), verif-caso-<Vn> (casos completos) e
# verif-testes. compilar.sh nao recompila com algum deles na fila (os executaveis carregam a libpz.so de build/).
VERIF_PKG="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VERIF_NEOPZ="$(cd "$VERIF_PKG/../.." && pwd)"
VERIF_MAQ="$(hostname -s 2>/dev/null || uname -n)"
VERIF_MAQ="${VERIF_MAQ%%.*}"
VERIF_MAQ="${VERIF_MAQ:-sem_nome}"
VERIF_BOOT="$(cat /proc/sys/kernel/random/boot_id 2>/dev/null || echo -)"
VERIF_EU="${USER:-$(id -un)}"

# projetos_verif: os projetos de Projects/CMakeLists.txt (as linhas add_subdirectory(V<n>_...)), um por linha
projetos_verif() {
    sed -n 's/^[[:space:]]*add_subdirectory([[:space:]]*\(V[0-9][0-9]*_[A-Za-z0-9_]*\)[[:space:]]*).*/\1/p' \
        "$VERIF_NEOPZ/Projects/CMakeLists.txt" 2>/dev/null
}

# projeto_completo <nome>: o nome completo do projeto (V2, v2 ou V2_TaludesCho -> V2_TaludesCho); erro se nao existir
projeto_completo() {
    local p achado=""
    for p in $(projetos_verif); do
        if [ "$p" = "${1:-}" ] || [ "${p%%_*}" = "${1^^}" ]; then achado="$p"; fi
    done
    if [ -z "$achado" ]; then
        echo "projeto desconhecido: '${1:-}' (os de Projects/CMakeLists.txt: $(projetos_verif | tr '\n' ' '))" >&2
        return 1
    fi
    echo "$achado"
}

# alvo_exe <projeto>, alvo_teste <projeto>: os alvos do CMakeLists.txt do projeto (add_verif_exe e add_verif_test)
alvo_exe() {
    local a
    a="$(sed -n 's/^[[:space:]]*add_verif_exe([[:space:]]*\([A-Za-z0-9_]*\).*/\1/p' \
        "$VERIF_NEOPZ/Projects/$1/CMakeLists.txt" 2>/dev/null | head -n 1)"
    echo "${a:-$1}"
}
alvo_teste() {
    local a
    a="$(sed -n 's/^[[:space:]]*add_verif_test([[:space:]]*\([A-Za-z0-9_]*\).*/\1/p' \
        "$VERIF_NEOPZ/Projects/$1/CMakeLists.txt" 2>/dev/null | head -n 1)"
    echo "${a:-${1%%_*}_testes}"
}

# versao_de <executavel>: a versao gravada por compilar.sh em <pasta do executavel>/versao.txt ("<commit> <data>"),
# ou - se nao houver
versao_de() {
    local v
    v="$(awk -v e="${1##*/}" '$1 == e { $1 = ""; sub(/^ +/, ""); print; exit }' "$(dirname "$1")/versao.txt" 2>/dev/null)"
    echo "${v:--}"
}

# jobs_verif: "<id> <nome> <estado>" dos jobs deste usuario com nome verif-*; codigo 1 se o squeue nao respondeu (a
# resposta dele vai para a saida de erros). Fora do cluster (sem squeue), nenhum job.
jobs_verif() {
    command -v squeue > /dev/null 2>&1 || return 0
    local out
    if ! out="$(squeue -h -u "$VERIF_EU" -o '%i %j %T' 2>&1)"; then
        echo "$out" >&2
        return 1
    fi
    echo "$out" | awk '$2 ~ /^verif-/'
}

# processos_verif <build>: "<pid> <executavel>" dos processos deste usuario nesta maquina que rodam um executavel de
# bin/ ou de <build>/Projects (os testes), inclusive um executavel ja trocado (o "(deleted)" do /proc)
processos_verif() {
    local b="$1" p e
    if [ -d /proc/self ]; then
        for p in /proc/[0-9]*; do
            e="$(readlink "$p/exe" 2>/dev/null)" || continue
            case "$e" in "$VERIF_PKG/bin/"* | "$b/Projects/"*) echo "${p#/proc/} $e" ;; esac
        done
    else
        pgrep -fl "$VERIF_PKG/bin/|$b/Projects/" 2>/dev/null || true
    fi
}

# fonte_nova <projeto> <arquivo>: o primeiro fonte do executavel do projeto mais novo que <arquivo> (por exemplo
# depois de um git pull): a biblioteca (tudo fora de Projects/), Projects/CMakeLists.txt, Projects/Comum e o projeto
fonte_nova() {
    local p="$1" ref="$2"
    (cd "$VERIF_PKG" && find ../.. \( -name .git -o -name build -o -name 'build_*' -o -name 'cmake-build*' -o -name _deps \
        -o -name deps -o -name bin -o -name 'resultados*' -o -name .claude -o -name logs \) -prune -o \
        \( -path '../../Projects/*' ! -path '../../Projects/Comum' ! -path '../../Projects/Comum/*' \
        ! -path "../../Projects/$p" ! -path "../../Projects/$p/*" ! -path '../../Projects/CMakeLists.txt' \) -prune -o \
        -newer "$ref" \( -name '*.cpp' -o -name '*.h' -o -name '*.hpp' -o -name CMakeLists.txt -o -name '*.cmake' \) \
        -print -quit 2>/dev/null) || true
}

# precisa_compilar <projeto> [...]: o motivo para compilar antes de submeter (vazio: o executavel de bin/ esta em dia).
# RECOMPILAR=1 sempre compila; FORCAR_SEM_COMPILAR=1 nao olha os fontes (so exige o executavel).
precisa_compilar() {
    local p e f
    for p in "$@"; do
        e="$VERIF_PKG/bin/$(alvo_exe "$p")"
        if ! [ -x "$e" ]; then echo "bin/${e##*/} ainda nao existe"; return; fi
    done
    if [ -n "${RECOMPILAR:-}" ]; then echo "RECOMPILAR=1"; return; fi
    [ -n "${FORCAR_SEM_COMPILAR:-}" ] && return
    for p in "$@"; do
        e="$VERIF_PKG/bin/$(alvo_exe "$p")"
        f="$(fonte_nova "$p" "$e")"
        if [ -n "$f" ]; then echo "fontes mais novos que bin/${e##*/} (por exemplo ${f#../../})"; return; fi
    done
}

# compilado_agora <projeto>: o executavel de bin/ saiu da ultima compilacao em segundo plano (ela o deixa com a data
# do seu inicio, a de logs/compilar.pid)
compilado_agora() {
    local e="$VERIF_PKG/bin/$(alvo_exe "$1")" r="$VERIF_PKG/logs/compilar.pid"
    [ -x "$e" ] && [ -f "$r" ] && [ "$(stat -c %Y "$e" 2>/dev/null)" = "$(stat -c %Y "$r" 2>/dev/null)" ]
}

# compilando_aqui: a compilacao em segundo plano (logs/compilar.pid: "pid maquina boot inicio") esta rodando nesta
# maquina
compilando_aqui() {
    local pid no boot
    [ -f "$VERIF_PKG/logs/compilar.pid" ] || return 1
    read -r pid no boot _ < "$VERIF_PKG/logs/compilar.pid" || return 1
    [ "$no" = "$VERIF_MAQ" ] && [ "$boot" = "$VERIF_BOOT" ] && kill -0 "$pid" 2>/dev/null &&
        tr '\0' ' ' 2>/dev/null < "/proc/$pid/cmdline" | grep -q 'compilar.sh'
}

# recusa_se_compilando: sai (codigo 1) se a compilacao em segundo plano esta rodando, a menos que este seja o pedido
# que ela submete ao fim (CADEIA): um job submetido agora rodaria com a libpz.so que ela esta trocando
recusa_se_compilando() {
    [ -n "${CADEIA:-}" ] && return 0
    compilando_aqui || return 0
    echo "a compilacao em segundo plano esta rodando neste no (pid $(cut -d' ' -f1 "$VERIF_PKG/logs/compilar.pid"), log em" >&2
    echo "logs/compilar.log). Ao fim ela submete so o pedido original ($(cat "$VERIF_PKG/logs/compilar.pedido" 2>/dev/null));" >&2
    echo "repita este pedido depois dela (acompanhe com ./status.sh)." >&2
    exit 1
}

# recusa_se_jobs <motivo>: sai (codigo 1) se ha jobs verif-* na fila, ou se o squeue nao responde (sem saber, nao
# recompila); FORCAR=1 passa por cima
recusa_se_jobs() {
    local fila
    [ -n "${FORCAR:-}" ] && return 0
    if ! fila="$(jobs_verif 2>&1)"; then
        echo "$1, mas o squeue nao respondeu ($fila): sem saber se ha jobs usando a libpz.so, nao recompila." >&2
        echo "Tente de novo daqui a pouco (ou FORCAR=1)." >&2
        exit 1
    fi
    if [ -n "$fila" ]; then
        echo "$1, mas ha $(echo "$fila" | wc -l) job(s) destes projetos na fila: recompilar agora trocaria a libpz.so que" >&2
        echo "eles usam. Espere terminarem (ou scancel; as realizacoes do Monte Carlo continuam depois) e rode de novo:" >&2
        echo "$fila" | sed 's/^/  /' >&2
        echo "Para so submeter com os executaveis atuais: FORCAR_SEM_COMPILAR=1 (os fontes mais novos ficam de fora)." >&2
        exit 1
    fi
}

# compila_em_segundo_plano <pedido> <comando> [args...]: roda compilar.sh (todos os projetos) neste no de login em
# segundo plano (nohup + setsid: continua depois que voce desconectar; os nos de computacao do fermi nao tem cmake) e,
# ao fim, de novo <comando> [args...] com CADEIA=<n> (o pedido de submissao, que entao so submete). Os executaveis que
# ela copia para bin/ ficam com a data do inicio da compilacao: um fonte trocado durante ela (um git pull) continua mais
# novo, e a cadeia compila de novo (ate 3 vezes). logs/compilar.pid guarda "pid maquina boot inicio",
# logs/compilar.pedido o pedido e logs/compilar.log a saida. Escreve o pid. O "exit $?" no fim impede o bash de trocar
# o seu processo pelo do ultimo comando (exec): a cadeia continua sendo "a compilacao" ate a submissao terminar.
compila_em_segundo_plano() {
    local pedido="$1"
    shift
    mkdir -p "$VERIF_PKG/logs"
    echo "$pedido" > "$VERIF_PKG/logs/compilar.pedido"
    # uma nova rodada da cadeia (um fonte mudou durante a compilacao) continua o mesmo log
    [ -n "${CADEIA:-}" ] || : > "$VERIF_PKG/logs/compilar.log"
    nohup setsid bash -c 'P="$1"; shift
        VERIF_EM_SEGUNDO_PLANO=1 nice -n 10 "$P/compilar.sh"
        find "$P/bin" -maxdepth 1 -type f -perm -u+x -newer "$P/logs/compilar.pid" \
            -exec touch -r "$P/logs/compilar.pid" {} + 2>/dev/null
        CADEIA=$(( ${CADEIA:-0} + 1 )) RECOMPILAR= "$@"
        exit $?' _ "$VERIF_PKG" "$@" \
        >> "$VERIF_PKG/logs/compilar.log" 2>&1 < /dev/null &
    echo "$! $VERIF_MAQ $VERIF_BOOT $(date '+%F %T')" > "$VERIF_PKG/logs/compilar.pid"
    echo "$!"
}

# compilar_antes "<projeto> [...]" <pedido> <comando> [args...]: o passo de compilacao comum dos submeter_*.sh. Sem
# motivo para compilar, volta (o chamador submete). Com motivo, compila em segundo plano e sai (o pedido e repetido ao
# fim, com CADEIA). Na cadeia, depois da compilacao, ainda com motivo: com CADEIA_SEGUE=1 volta assim mesmo (o chamador
# submete; o job avisa o que ficou velho); senao compila de novo se os executaveis sairam desta compilacao (um fonte
# mudou durante ela; ate 3 vezes), ou avisa (tambem em AVISOS_CADEIA, se definido) e sai.
compilar_antes() {
    local ps="$1" pedido="$2" motivo pid p todos=1
    shift 2
    # shellcheck disable=SC2086
    motivo="$(precisa_compilar $ps)"
    [ -n "$motivo" ] || return 0
    if [ -n "${CADEIA:-}" ]; then
        [ "${CADEIA_SEGUE:-0}" = 1 ] && { echo "aviso: $motivo (a compilacao de agora nao resolveu)" >&2; return 0; }
        for p in $ps; do compilado_agora "$p" || todos=0; done
        if [ "$todos" = 0 ] || [ "$CADEIA" -ge 3 ]; then
            local m="compilacao terminada, mas $motivo: nada submetido (veja logs/compilar.log; corrija e repita o pedido)"
            echo "$(date '+%F %T') $m" >&2
            [ -n "${AVISOS_CADEIA:-}" ] && echo "$(date '+%F %T') $m" >> "$AVISOS_CADEIA"
            exit 1
        fi
        echo "$motivo (um fonte mudou durante a compilacao): compilando de novo"
    fi
    recusa_se_jobs "$motivo"
    pid="$(compila_em_segundo_plano "$pedido" "$@")"
    echo "$motivo: compilando neste no de login em segundo plano (pid $pid, log em logs/compilar.log; a biblioteca"
    echo "inteira leva uns 15 a 30 min, so os projetos alguns minutos). Ao terminar, o pedido e submetido sozinho."
    echo "Pode desconectar; acompanhe com ./status.sh. Para desistir: kill -- -$pid"
    exit 0
}

# estado_compilacao: uma linha sobre a compilacao em segundo plano, se ela esta rodando ou se terminou mal
estado_compilacao() {
    local L="$VERIF_PKG/logs" cpid cno cboot cini
    [ -f "$L/compilar.pid" ] || return 0
    read -r cpid cno cboot cini < "$L/compilar.pid"
    if compilando_aqui; then
        echo "compilacao: rodando desde $cini (pid $cpid): $(grep -o '\[ *[0-9]*%\]' "$L/compilar.log" 2>/dev/null |
            tail -n 1); ao fim submete: $(cat "$L/compilar.pedido" 2>/dev/null)"
    elif [ -f "$L/compilar.log" ] && ! grep -q '^Compilado:' "$L/compilar.log"; then
        echo "compilacao: terminou SEM gerar executaveis (comecou em $cini); veja: tail -30 $L/compilar.log"
    else
        [ -f "$L/compilar.log" ] && grep -q '^Falhou:' "$L/compilar.log" &&
            echo "compilacao: $(grep '^Falhou:' "$L/compilar.log" | tail -n 1 | sed 's/ (veja.*//') (veja $L/compilar.log)"
        [ -f "$L/compilar.log" ] && grep -q 'recusou\|nada submetido\|nenhum job submetido\|nao submetido' "$L/compilar.log" &&
            echo "compilacao: terminou, mas a submissao falhou; veja: tail -5 $L/compilar.log"
    fi
    return 0
}

# nucleos_fisicos: os nucleos fisicos desta maquina (com hyper-threading, metade das CPUs)
nucleos_fisicos() {
    local n
    n="$(lscpu -p=CORE,SOCKET 2>/dev/null | grep -v '^#' | sort -u | wc -l)"
    [ "${n:-0}" -ge 1 ] 2>/dev/null || n="$(nproc 2>/dev/null || echo 1)"
    echo "$n"
}

# paralelo_job: os processos de 1 CPU de um job: as CPUs dele (SLURM_CPUS_PER_TASK) divididas pelas threads por nucleo
# (um no inteiro do fermi: 24 CPUs, 12 nucleos fisicos -> 12); fora do SLURM, os nucleos fisicos
paralelo_job() {
    local tpc
    if [ -n "${SLURM_CPUS_PER_TASK:-}" ]; then
        tpc="$(LC_ALL=C lscpu 2>/dev/null | awk -F: '/^Thread\(s\) per core/ { gsub(/ /, "", $2); print $2 }')"
        case "$tpc" in '' | *[!0-9]* | 0) tpc=1 ;; esac
        echo $((SLURM_CPUS_PER_TASK / tpc > 0 ? SLURM_CPUS_PER_TASK / tpc : 1))
    else
        nucleos_fisicos
    fi
}

# trocas_caso <entrada>: as variaveis de uma entrada "<caso>[:VAR=valor,...]" dos casos completos, uma por linha
# (VAR=valor); um pedaco sem '=' continua o valor anterior, que pode ser uma lista (NU_LISTA=0.2,0.25,0.3 do V5)
trocas_caso() {
    local t atual="" pedacos
    [[ "$1" == *:* ]] || return 0
    IFS=, read -ra pedacos <<< "${1#*:}"
    for t in "${pedacos[@]}"; do
        if [[ "$t" == *=* ]]; then
            [ -n "$atual" ] && echo "$atual"
            atual="$t"
        elif [ -n "$t" ]; then
            atual="$atual,$t"
        fi
    done
    [ -n "$atual" ] && echo "$atual"
    return 0
}

# valida_caso <entrada>: a mensagem de erro de uma entrada dos casos completos (vazia: valida). O caso e o nome da pasta
# (e o CASO do programa, a menos que as trocas definam CASO): letras, numeros e . _ + -
valida_caso() {
    local c="${1%%:*}" t
    [[ "$c" =~ ^[A-Za-z0-9._+-]+$ ]] || { echo "\"$1\": o caso \"$c\" (o nome da pasta) so pode ter letras, numeros e . _ + -"; return 0; }
    if [[ "$1" == *:* ]]; then
        [[ "${1#*:}" == *=* ]] || { echo "\"$1\": nada como VAR=valor depois de ':'"; return 0; }
        [[ "${1#*:}" =~ ^[A-Za-z_][A-Za-z0-9_]*= ]] || { echo "\"$1\": o que vem depois de ':' deve comecar por VAR=valor"; return 0; }
        while IFS= read -r t; do
            [[ "${t%%=*}" =~ ^[A-Z_][A-Z0-9_]*$ ]] || { echo "\"$1\": \"${t%%=*}\" nao e um nome de variavel (maiusculas, numeros e _)"; return 0; }
        done < <(trocas_caso "$1")
    fi
    return 0
}

# aviso_ldd <executavel>: vazio se o executavel existe e acha as bibliotecas dele (libpz.so de build/, o MKL do modulo);
# senao, o que falta
aviso_ldd() {
    if ! [ -x "$1" ]; then echo "$1 nao existe (rode compilar.sh no no de login)"; return; fi
    local f
    f="$(ldd "$1" 2>&1 | grep 'not found' | tr -s ' \t' ' ' | tr '\n' ';')"
    [ -z "$f" ] || echo "$1: bibliotecas faltando:$f"
}
