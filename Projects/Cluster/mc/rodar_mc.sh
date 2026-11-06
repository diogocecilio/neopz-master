#!/bin/bash
# rodar_mc.sh <primeira> <ultima>: roda as realizacoes primeira..ultima do Monte Carlo de um projeto de verificacao
# (PROJETO=V2_TaludesCho ou V4_TaludeLiu, ou so V2 ou V4; ver projetos.sh), cada uma num processo de 1 CPU, PARALELO ao
# mesmo tempo (padrao: os nucleos fisicos da maquina; com hyper-threading, dois processos por nucleo deixam cada um ~2x
# mais lento). Os casos vao numa fila so, CASOS="<caso>[:VAR=valor,...] ..." (padrao: a campanha do artigo, em
# projetos.sh), em ordem intercalada (a realizacao primeira de todos, depois a seguinte, ...); DET=1 poe antes o caso
# deterministico (CAMPOS=0, as medias) de cada entrada, e CAMPOS=0 roda so os deterministicos. VAR=valor troca um
# parametro do programa so nessa entrada, que vai para a sua propria pasta (estudos de sensibilidade):
# CASOS="cho1 cho1:NUCLEO=exp cho1:COV_C=0.5,LX=40" (VAR: os parametros de projetos.sh; valor: letras, numeros, . _ + -).
# Cada realizacao fica em $SAIDA/<projeto>/<caso>[_trocas]/r<nnnnn>/ (as trocas em ordem alfabetica, cada uma com o nome
# da variavel em minusculas e sem _ seguido do valor: cho1:LX=40,COV_C=0.5 -> cho1_covc0.5_lx40), e o deterministico em
# <caso>[_trocas]_det/r00000/, sem as trocas que so mudam os campos (que ele ignora: cho1:COV_C=0.5 usa cho1_det). No
# fim, junta_mc.sh junta as concluidas em mc_resumo.csv e mc_fs.csv de cada pasta; ./resumo_mc.sh <pasta> mostra P_f, e
# posprocessa_mc.py faz o pos-processamento completo.
# Retomavel, e seguro com varios processos ou jobs (em nos diferentes) na mesma fila: cada realizacao em andamento tem um
# dono (r<nnnnn>/.ativo: maquina, processo e job do SLURM). As concluidas (resultado.txt) e as de um dono vivo sao
# puladas; a de um dono morto (job que acabou, processo que caiu, maquina reiniciada) volta para a fila (e continua do
# seu checkpoint.bin, se o programa gravar um). Uma realizacao que termina com erro (sem resultado.txt) MAX_FALHAS vezes
# (padrao 3) fica marcada (falhou.txt) e e pulada. Com o arquivo $SAIDA/PARAR (todos os projetos) ou
# $SAIDA/<projeto>/PARAR, nenhuma realizacao nova comeca (as em andamento terminam).
# rodar_mc.sh --faltam <primeira> <ultima>: so conta as realizacoes ainda nao concluidas nem marcadas como falhas;
# --pastas: "<entrada> <pasta> <pasta do deterministico ou ->" (separados por tabulacao) de cada entrada de CASOS;
# --verifica: so confere CASOS e os parametros.txt das pastas.
# Avisos (parametros recusados, disjuntor) tambem vao para $SAIDA/<projeto>/AVISOS, que status_mc.sh mostra.
# Ambiente: PROJETO, CASOS, DET, CAMPOS, SAIDA (padrao ../resultados_mc), EXE (padrao ../bin/<executavel>), PARALELO,
# MAX_FALHAS, DISJUNTOR, INATIVO_MIN, ABANDONO_H, ESPERA_S, SO_CAMPOS=1 (so os campos das realizacoes pedidas, sem
# simular, nos projetos que o tem) e qualquer parametro do programa (vale para todas as entradas; vai para parametros.txt).
set -u
export LC_ALL=C # numeros com ponto decimal em qualquer locale (awk, printf)
MC="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PKG="$(cd "$MC/.." && pwd)"
MODO=rodar
case "${1:-}" in --faltam) MODO=faltam; shift ;; --pastas) MODO=pastas; set -- 0 0 ;; --verifica) MODO=verifica; set -- 0 0 ;; esac
[ $# -eq 2 ] || { sed -n '2,/^set -u/{/^#/p}' "$0"; exit 2; }
[[ "$1" =~ ^[0-9]+$ && "$2" =~ ^[0-9]+$ ]] || { echo "<primeira> e <ultima> devem ser inteiros >= 0 (recebi '$1' e '$2')" >&2; exit 2; }
# shellcheck source=../comum.sh
source "$PKG/comum.sh"
# shellcheck source=projetos.sh
source "$MC/projetos.sh"
[ -n "${PROJETO:-}" ] || { echo "defina PROJETO (os de Monte Carlo: $(mc_projetos); V2 e V4 bastam)" >&2; exit 2; }
mc_projeto "$PROJETO" || exit 2
export EXE="${EXE:-$PKG/bin/$MC_EXE}"
# a versao do programa (gravada por compilar.sh em bin/versao.txt), que ele grava no resumo.csv de cada realizacao
[ -n "${!MC_VAR_VERSAO:-}" ] || export "$MC_VAR_VERSAO=$(versao_de "$EXE")"
if [ "$MODO" = rodar ] && ! [ -x "$EXE" ]; then
    echo "executavel $EXE nao encontrado (compile com ../compilar.sh no no de login; ver ../README.md)" >&2
    exit 2
fi
export MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 MAX_FALHAS="${MAX_FALHAS:-3}"
SAIDA_BASE="${SAIDA:-$PKG/resultados_mc}"
export SAIDA_BASE SAIDA_MC="$SAIDA_BASE/$MC_PROJETO"
CASOS="$(echo ${CASOS:-$MC_CASOS_PADRAO})"
export CAMPOS_TODOS="${CAMPOS:-1}" HOST_MC="$VERIF_MAQ" BOOT_MC="$VERIF_BOOT"
[ "$CAMPOS_TODOS" = 0 ] || [ "$CAMPOS_TODOS" = 1 ] || { echo "CAMPOS deve ser 0 ou 1" >&2; exit 2; }

# ---------------------------------------------------------------- as entradas de CASOS
# pasta <caso> <campos> <trocas>: poe em PASTA a pasta da entrada ($SAIDA/<projeto>/<caso>[_trocas][_det])
pasta() {
    local c="$SAIDA_MC/$1" e v
    if [ "${3:--}" != - ]; then
        for e in ${3//,/ }; do v="${e%%=*}"; v="${v//_/}"; c="${c}_${v,,}${e#*=}"; done
    fi
    [ "$2" = 0 ] && c="${c}_det"
    PASTA="$c"
}
# sem_campo <trocas>: as trocas sem as que so mudam os campos (as do deterministico), ou -
sem_campo() {
    local e d=""
    if [ "${1:--}" != - ]; then
        for e in ${1//,/ }; do [[ " $MC_VARS_CAMPO " == *" ${e%%=*} "* ]] || d="${d:+$d,}$e"; done
    fi
    echo "${d:--}"
}
# valida_entrada <entrada>: a mensagem de erro (vazia: a entrada e valida)
valida_entrada() {
    local e="$1" c="${1%%:*}" x="" t v vistos=" "
    [[ "$c" =~ $MC_CASOS_RE ]] || { echo "\"$e\": caso \"$c\" desconhecido em $MC_PROJETO (aceitos: $MC_CASOS_RE)"; return; }
    [[ "$e" == *:* ]] || return 0
    x="${e#*:}"
    [ -n "$x" ] || { echo "\"$e\": nada depois de ':'"; return; }
    for t in ${x//,/ }; do
        v="${t%%=*}"
        [[ "$t" == *=* ]] || { echo "\"$e\": \"$t\" nao e VAR=valor"; return; }
        [[ " $MC_VARS " == *" $v "* ]] || { echo "\"$e\": $v nao e um parametro de $MC_PROJETO (aceitos: $MC_VARS)"; return; }
        [[ "${t#*=}" =~ ^[A-Za-z0-9._+-]+$ ]] || { echo "\"$e\": valor \"${t#*=}\" de $v (so letras, numeros e . _ + -)"; return; }
        [[ "$vistos" == *" $v "* ]] && { echo "\"$e\": $v repetido"; return; }
        vistos="$vistos$v "
    done
}
# as entradas, uma vez: E_C (caso), E_X (trocas em ordem alfabetica, ou -), E_XD (as do deterministico), E_D (pasta) e
# E_DD (pasta do deterministico)
E_C=() E_X=() E_XD=() E_D=() E_DD=() E_ENT=()
for e in $CASOS; do
    m="$(valida_entrada "$e")"
    if [ -n "$m" ]; then
        echo "CASOS: $m" >&2
        echo "(formato: CASOS=\"<caso>[:VAR=valor,...] ...\", por exemplo CASOS=\"$MC_CASOS_PADRAO\")" >&2
        exit 2
    fi
    x=-
    [[ "$e" == *:* ]] && x="$(printf '%s\n' $(echo "${e#*:}" | tr ',' ' ') | sort -t= -k1,1 | paste -sd, -)"
    xd="$(sem_campo "$x")"
    pasta "${e%%:*}" 1 "$x"
    for k in "${!E_D[@]}"; do
        [ "${E_D[$k]}" = "$PASTA" ] && { echo "CASOS: \"${E_ENT[$k]}\" e \"$e\" vao para a mesma pasta ($PASTA)" >&2; exit 2; }
    done
    E_ENT+=("$e") E_C+=("${e%%:*}") E_X+=("$x") E_XD+=("$xd") E_D+=("$PASTA")
    pasta "${e%%:*}" 0 "$xd"
    E_DD+=("$PASTA")
done
[ ${#E_C[@]} -gt 0 ] || { echo "CASOS vazio" >&2; exit 2; }
if [ "${SO_CAMPOS:-0}" = 1 ] && [ "$MC_SO_CAMPOS" != 1 ]; then
    echo "SO_CAMPOS=1: $MC_PROJETO nao tem esse modo (o programa rodaria a realizacao inteira)" >&2
    exit 2
fi
if [ "$MODO" = pastas ]; then # "<entrada> <pasta> <pasta do deterministico ou ->", para status_mc.sh e submeter_mc.sh
    for k in "${!E_C[@]}"; do
        dd=-
        { [ "${DET:-0}" = 1 ] || [ "$CAMPOS_TODOS" = 0 ]; } && dd="${E_DD[$k]}"
        printf '%s\t%s\t%s\n' "${E_ENT[$k]}" "${E_D[$k]}" "$dd"
    done
    exit 0
fi

# tarefas <primeira> <ultima>: uma linha por realizacao, "<caso> <campos> <n> <trocas>": antes os deterministicos (com
# DET=1, um por pasta), depois as realizacoes intercaladas
tarefas() {
    local k n vistos=" "
    if [ "${DET:-0}" = 1 ] || [ "$CAMPOS_TODOS" = 0 ]; then
        for k in "${!E_C[@]}"; do
            case "$vistos" in *" ${E_C[$k]}:${E_XD[$k]} "*) continue ;; esac
            vistos="$vistos${E_C[$k]}:${E_XD[$k]} "
            echo "${E_C[$k]} 0 0 ${E_XD[$k]}"
        done
    fi
    [ "$CAMPOS_TODOS" = 0 ] && return
    for n in $(seq "$1" "$2"); do
        for k in "${!E_C[@]}"; do echo "${E_C[$k]} 1 $n ${E_X[$k]}"; done
    done
}
# pendentes <primeira> <ultima>: as tarefas sem resultado.txt nem falhou.txt
pendentes() {
    local c cp n x r
    tarefas "$1" "$2" | while read -r c cp n x; do
        pasta "$c" "$cp" "$x"
        printf -v r '%s/r%05d' "$PASTA" "$n"
        [ -s "$r/resultado.txt" ] || [ -e "$r/falhou.txt" ] || echo "$c $cp $n $x"
    done
}
aviso_mc() { echo "$(date '+%F %T') $HOST_MC${SLURM_JOB_ID:+ job $SLURM_JOB_ID}: $*" >> "$SAIDA_MC/AVISOS"; }
parar() {
    [ -e "$SAIDA_BASE/PARAR" ] || [ -e "$SAIDA_MC/PARAR" ] || [ -e "$INST_MC/parar" ] ||
        { [ -n "${SLURM_JOB_ID:-}" ] && [ -e "$SAIDA_MC/.parar.$SLURM_JOB_ID" ]; }
}
# job_vivo <id>: 0 se o job esta vivo, 1 se acabou, 2 se esta acabando (COMPLETING: os processos estao sendo mortos),
# 3 se o SLURM nao respondeu
job_vivo() {
    local out rc
    out="$(squeue -h -j "$1" -o %T 2>&1)"
    rc=$?
    if [ $rc -eq 0 ]; then
        case "$out" in
        '' | COMPLETED* | FAILED* | CANCELLED* | TIMEOUT* | NODE_FAIL* | PREEMPTED* | OUT_OF_MEMORY* | BOOT_FAIL* | DEADLINE*) return 1 ;;
        COMPLETING* | STOPPED* | SIGNALING*) return 2 ;;
        esac
        return 0
    fi
    case "$out" in *"Invalid job id"* | *"invalid job id"*) return 1 ;; esac
    return 3
}
# pid_vivo <pid>: o processo existe, nao e um zumbi (morto, esperando o pai o recolher) e e um trabalhador desta fila (o
# "bash -c uma ..." abaixo): um pid reaproveitado por outro programa nao conta
pid_vivo() {
    kill -0 "$1" 2>/dev/null || return 1
    [ -r "/proc/$1/stat" ] || return 0 # sem /proc: so o kill -0
    case "$(sed -n 's/^.*) \([A-Za-z]\).*/\1/p' "/proc/$1/stat" 2>/dev/null)" in Z | X | x) return 1 ;; esac
    tr '\0' ' ' 2>/dev/null < "/proc/$1/cmdline" | grep -q 'uma "$@"'
}
# ativa <r> <min>: algum arquivo da realizacao (nao as travas) mudou nos ultimos <min> minutos
ativa() { [ -n "$(find "$1" -mindepth 1 -maxdepth 1 ! -name .ativo ! -name .tomada -mmin "-$2" 2>/dev/null | head -n 1)" ]; }
# dono_vivo <r/.ativo>: 0 se o dono pode estar rodando a realizacao, 1 se morreu, 2 se o job dele esta acabando (ou
# acabou ha pouco), 3 se o SLURM nao respondeu. Um job que o SLURM da como acabado so conta como morto depois de
# INATIVO_MIN minutos (padrao 30) sem mudanca na realizacao: um no que perdeu contato com o controlador (NODE_FAIL) pode
# estar com os processos vivos, gravando na pasta
dono_vivo() {
    local a="$1" r="${1%/.ativo}" h p j b t job vez e
    if ! read -r h p j b t 2>/dev/null < "$a/dono"; then
        # sem dono legivel (o processo caiu logo depois de criar a trava, ou uma falha de leitura): viva enquanto a
        # trava for recente ou a realizacao estiver mudando
        [ -n "$(find "$a" -maxdepth 0 -mmin -10 2>/dev/null)" ] && return 0
        ativa "$r" 60
        return
    fi
    if [[ "${j:--}" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
        job="${j%%.*}" vez="${j#*.}"
        [ "$vez" = "$j" ] && vez=0
        if [ -n "${SLURM_JOB_ID:-}" ] && [ "$job" = "$SLURM_JOB_ID" ] && [ "$vez" = "${SLURM_RESTART_COUNT:-0}" ]; then
            # o proprio job, desta vez: vale o processo; em outro no do mesmo job (job de varios nos), vivo
            [ "$h" = "$HOST_MC" ] || return 0
            pid_vivo "$p"
            return
        fi
        if [ -n "${SLURM_JOB_ID:-}" ] && [ "$job" = "$SLURM_JOB_ID" ]; then
            e=1 # o proprio job numa vez anterior (o SLURM o recolocou na fila): acabou
        else
            job_vivo "$job"
            e=$?
        fi
        if [ $e -eq 1 ] && ativa "$r" "${INATIVO_MIN:-30}"; then return 2; fi
        return $e
    fi
    if [ "$h" = "$HOST_MC" ]; then [ "$b" = "$BOOT_MC" ] && pid_vivo "$p"; return; fi
    # outra maquina, fora do SLURM: viva enquanto a realizacao muda
    ativa "$r" "$((${ABANDONO_H:-6} * 60))"
}
# escreve_dono <r/.ativo>: "<maquina> <pid> <job>.<vez> <boot> <inicio>" (job "-" fora do SLURM; vez: SLURM_RESTART_COUNT)
escreve_dono() {
    local j=-
    [ -n "${SLURM_JOB_ID:-}" ] && j="$SLURM_JOB_ID.${SLURM_RESTART_COUNT:-0}"
    echo "$HOST_MC $$ $j $BOOT_MC $(date +%s)" > "$1/dono.tmp" && mv -f "$1/dono.tmp" "$1/dono"
}
# libera <r>: apaga a trava, se ainda for deste processo (um dono que a perdeu nao apaga a do novo)
libera() {
    local h p j b t
    if read -r h p j b t 2>/dev/null < "$1/.ativo/dono" && { [ "$h" != "$HOST_MC" ] || [ "$p" != "$$" ]; }; then
        echo "$(date '+%F %T') $1: a trava agora e de outro dono ($h $p $j); nao apagada" >&2
        return 0
    fi
    rm -rf "$1/.ativo"
}
# reivindica <r>: 0 se este processo passou a ser o dono da realizacao. A trava e um mkdir (atomico, tambem no NFS);
# a de um dono morto e trocada sob uma segunda trava (.tomada), para que so um processo a assuma
reivindica() {
    local r="$1" ok=1
    if mkdir "$r/.ativo" 2>/dev/null; then
        escreve_dono "$r/.ativo" && return 0
        rm -rf "$r/.ativo" # sem o dono gravado, a trava nao protege: desiste
        return 1
    fi
    dono_vivo "$r/.ativo"
    [ $? -eq 1 ] || return 1
    if ! mkdir "$r/.tomada" 2>/dev/null; then
        # outro processo esta trocando; uma .tomada velha (o processo caiu durante a troca) e apagada para a proxima vez
        [ -d "$r/.tomada" ] && [ -z "$(find "$r/.tomada" -maxdepth 0 -mmin -10 2>/dev/null)" ] && rm -rf "$r/.tomada"
        return 1
    fi
    if [ -d "$r/.ativo" ]; then
        dono_vivo "$r/.ativo"
        [ $? -eq 1 ] && rm -rf "$r/.ativo"
    fi
    if mkdir "$r/.ativo" 2>/dev/null; then
        if escreve_dono "$r/.ativo"; then ok=0; else rm -rf "$r/.ativo"; fi
    fi
    rmdir "$r/.tomada" 2>/dev/null
    return $ok
}
# uma <caso> <campos> <n> <trocas>: uma realizacao, na sua pasta; a saida do programa vai para run.log e run.err
uma() {
    local caso="$1" campos="$2" n="$3" x="${4:--}" d r rc nf nome dica
    local -a ev=() ex=()
    [ "$x" != - ] && ev=(${x//,/ })
    pasta "$caso" "$campos" "$x"
    d="$PASTA"
    printf -v r '%s/r%05d' "$d" "$n"
    nome="${d##*/}/${r##*/}"
    if [ "${SO_CAMPOS:-0}" = 1 ]; then # so os campos desta realizacao (campos.csv), sem simulacao
        [ "$campos" = 1 ] || return 0
        mkdir -p "$r" && cd "$r" || return 1
        if env ${ev[@]+"${ev[@]}"} CASO="$caso" CAMPOS=1 REAL="$n" CAMPOS_CSV=1 SO_CAMPOS=1 "$EXE" > campos.log 2>&1; then
            echo "$nome: campos em $r/campos.csv"
        else
            echo "$nome: erro ao gerar os campos; veja $r/campos.log" >&2
        fi
        return 0
    fi
    [ -s "$r/resultado.txt" ] && return 0
    parar && return 0
    mkdir -p "$r" || return 1
    [ -e "$r/falhou.txt" ] && return 0
    reivindica "$r" || return 0
    # terminou, ou foi marcada como falha, entre o teste e a trava
    if [ -s "$r/resultado.txt" ] || [ -e "$r/falhou.txt" ]; then libera "$r"; return 0; fi
    [ -s "$r/checkpoint.bin" ] && echo "$(date '+%F %T') $nome: retomando do checkpoint"
    # variaveis a mais do projeto (projetos.sh: mc_extra), com a pasta do deterministico desta entrada
    pasta "$caso" 0 "$(sem_campo "$x")"
    mapfile -t ex < <(mc_extra "$caso" "$campos" "$PASTA")
    (cd "$r" && env ${ev[@]+"${ev[@]}"} ${ex[@]+"${ex[@]}"} CASO="$caso" CAMPOS="$campos" REAL="$n" "$EXE" \
        >> run.log 2>> run.err)
    rc=$?
    if [ -s "$r/resultado.txt" ]; then
        rm -f "$r/checkpoint.bin" "$r/checkpoint.bin.tmp" # sobra de uma queda entre o resultado e a limpeza
        echo 0 > "$INST_MC/seguidas"
        echo "$(date '+%F %T') $nome: $(cat "$r/resultado.txt")"
    elif [ $rc -eq 129 ] || [ $rc -eq 130 ]; then
        # Ctrl-C ou terminal fechado (o trabalhador sobreviveu): interrupcao, nao falha
        echo "$(date '+%F %T') $nome: interrompida (sinal $((rc - 128))); volta para a fila" >&2
    else
        # falha: inclusive um SIGKILL/SIGTERM so no programa (o limite de tempo, o scancel e o Ctrl-C matam tambem o
        # trabalhador, que entao nem chega aqui): em geral o OOM killer
        nf=$(($(cat "$r/falhas" 2>/dev/null || echo 0) + 1))
        echo "$nf" > "$r/falhas"
        dica=""
        [ $rc -eq 137 ] && dica=" (SIGKILL: falta de memoria?)"
        if [ "$nf" -ge "$MAX_FALHAS" ]; then
            echo "$(date '+%F %T') $HOST_MC codigo $rc em $nf tentativas$dica" > "$r/falhou.txt"
            echo "$(date '+%F %T') $nome: falhou $nf vezes (codigo $rc$dica), marcada em falhou.txt; veja $r/run.err" >&2
        else
            echo "$(date '+%F %T') $nome: terminou com codigo $rc$dica sem resultado.txt (tentativa $nf de $MAX_FALHAS)" >&2
        fi
        # disjuntor: falhas seguidas demais nesta instancia, sem nenhum sucesso no meio, indicam um problema do ambiente
        # (executavel, bibliotecas, malha, disco) e nao das realizacoes: para de pegar realizacoes antes de marcar a fila
        local k
        k="$( (flock 9; k=$(($(cat "$INST_MC/seguidas" 2>/dev/null || echo 0) + 1)); echo $k > "$INST_MC/seguidas"; echo $k) 9> "$INST_MC/trava")"
        if [ "${k:-0}" -ge "$DISJUNTOR_MC" ] && ! [ -e "$INST_MC/parar" ]; then
            touch "$INST_MC/parar"
            {
                echo "$(date '+%F %T') $k falhas seguidas sem nenhum sucesso nesta maquina: parando de pegar realizacoes."
                echo "  Provavelmente um problema do ambiente, nao das realizacoes. O fim do run.err da ultima ($r):"
                tail -n 5 "$r/run.err" 2>/dev/null | sed 's/^/    /'
                echo "  Depois de corrigir, rode de novo; para zerar as tentativas: find $SAIDA_MC -name falhas -delete"
            } >&2
            aviso_mc "disjuntor: $k falhas seguidas sem nenhum sucesso; parou de pegar realizacoes (veja $r/run.err)"
        fi
    fi
    libera "$r"
}
export -f pasta sem_campo aviso_mc parar job_vivo pid_vivo ativa dono_vivo escreve_dono libera reivindica uma
# disponiveis <primeira> <ultima>: "<livres> <esperando>": as tarefas pendentes sem dono ou com dono morto (que uma nova
# passada pode rodar) e as que ficam livres em breve (jobs acabando ou que acabaram ha pouco, troca de dono em andamento,
# SLURM sem responder)
disponiveis() {
    local c cp n x r e
    pendentes "$1" "$2" | while read -r c cp n x; do
        pasta "$c" "$cp" "$x"
        printf -v r '%s/r%05d' "$PASTA" "$n"
        e=1
        if [ -d "$r/.tomada" ]; then e=2
        elif [ -d "$r/.ativo" ]; then dono_vivo "$r/.ativo"; e=$?
        fi
        echo "$e"
    done | awk '$1 == 1 { l++ } $1 >= 2 { a++ } END { print l + 0, a + 0 }'
}

if [ "$MODO" = faltam ]; then
    pendentes "$1" "$2" | wc -l
    exit 0
fi
P="${PARALELO:-$(nucleos_fisicos)}"
[ "$P" -ge 1 ] 2>/dev/null || P=1
# estado desta instancia (as falhas seguidas, para o disjuntor), apagado no fim
INST_MC="$(mktemp -d)"
export INST_MC DISJUNTOR_MC="${DISJUNTOR:-$((2 * P > MAX_FALHAS + 2 ? 2 * P : MAX_FALHAS + 2))}"
trap 'rm -rf "$INST_MC"' EXIT
# parametros.txt de cada pasta: os parametros do programa (MC_VARS de projetos.sh; o caso e as trocas tambem estao no nome
# da pasta) com que as suas realizacoes foram rodadas. Outra rodada na mesma pasta com outros parametros misturaria
# resultados: recusada (use outra SAIDA, ou uma entrada de CASOS com as trocas). O deterministico ignora os parametros
# dos campos, que ficam fora do parametros.txt dele.
VARS_PARAM_NOVAS=" " # variaveis acrescentadas depois: so entram quando definidas, para as pastas de antes valerem
parametros() { # parametros <trocas> <campos>: uma linha por variavel, VAR=valor (- se nao definida: o padrao do programa)
    local v val e
    for v in $MC_VARS; do
        [ "$2" = 0 ] && [[ " $MC_VARS_CAMPO " == *" $v "* ]] && continue
        val="${!v:--}"
        for e in ${1//,/ }; do [ "${e%%=*}" = "$v" ] && val="${e#*=}"; done
        [[ "$val" == */* ]] && val="${val##*/}" # um caminho: so o nome do arquivo (o mesmo no PC e no cluster)
        [ "$val" = - ] && [[ "$VARS_PARAM_NOVAS" == *" $v "* ]] && continue
        echo "$v=$val"
    done
}
if [ "${SO_CAMPOS:-0}" != 1 ]; then
    [ "$MODO" = rodar ] && mkdir -p "$SAIDA_MC"
    for k in "${!E_C[@]}"; do
        for cp in 1 0; do
            [ $cp = 0 ] && [ "${DET:-0}" != 1 ] && [ "$CAMPOS_TODOS" != 0 ] && continue
            [ $cp = 1 ] && [ "$CAMPOS_TODOS" = 0 ] && continue
            if [ $cp = 1 ]; then d="${E_D[$k]}" xx="${E_X[$k]}"; else d="${E_DD[$k]}" xx="${E_XD[$k]}"; fi
            novo="$(parametros "$xx" $cp)"
            if [ -f "$d/parametros.txt" ]; then
                if [ "$(cat "$d/parametros.txt")" != "$novo" ]; then
                    dif="$(diff <(cat "$d/parametros.txt") <(echo "$novo") | grep '^[<>]' | sed 's/^</antes:/; s/^>/agora:/' | tr '\n' ' ')"
                    echo "$d ja tem realizacoes rodadas com outros parametros: $dif" >&2
                    echo "Use outra SAIDA (ou uma entrada de CASOS com as trocas, como ${E_C[$k]}:COV_C=0.5)." >&2
                    [ "$MODO" = rodar ] && aviso_mc "recusado: ${d##*/} ja tem realizacoes com outros parametros: $dif"
                    exit 2
                fi
            elif [ "$MODO" = rodar ]; then
                mkdir -p "$d" && echo "$novo" > "$d/parametros.txt"
            fi
            # a pasta do deterministico desta entrada (para posprocessa_mc.py)
            if [ "$MODO" = rodar ] && [ $cp = 1 ] && [ "${DET:-0}" = 1 ] &&
                [ "$(cat "$d/deterministico.txt" 2>/dev/null)" != "${E_DD[$k]##*/}" ]; then
                echo "${E_DD[$k]##*/}" > "$d/deterministico.txt"
            fi
        done
    done
fi
[ "$MODO" = verifica ] && exit 0
mkdir -p "$SAIDA_MC"
echo "$(date '+%F %T') $MC_PROJETO $CASOS$([ "${DET:-0}" = 1 ] && echo ' (e deterministicos)'): realizacoes $1..$2, $P ao" \
    "mesmo tempo, $HOST_MC${SLURM_JOB_ID:+, job $SLURM_JOB_ID}" | tee -a "$SAIDA_MC/execucoes.txt"
# passadas pela fila: depois de cada uma, as realizacoes de donos que morreram nesse meio tempo (ou que falharam menos
# de MAX_FALHAS vezes) voltam para a fila. As de um job que esta acabando (por exemplo o job anterior, de que esta
# continuacao depende, cujas realizacoes so ficam livres INATIVO_MIN minutos depois de pararem de mudar) sao esperadas
# por ate INATIVO_MIN + 10 minutos. Termina quando o que falta esta com donos vivos (outros jobs)
espera=0
max_espera=$(( ((${INATIVO_MIN:-30} + 10) * 60 + ${ESPERA_S:-30} - 1) / ${ESPERA_S:-30} ))
passada=1
while :; do
    parar && break
    if [ "${SO_CAMPOS:-0}" = 1 ]; then
        tarefas "$1" "$2" | xargs -r -P "$P" -L 1 bash -c 'uma "$@"' _
        break
    fi
    pendentes "$1" "$2" | xargs -r -P "$P" -L 1 bash -c 'uma "$@"' _
    parar && break
    read -r livres acabando <<< "$(disponiveis "$1" "$2")"
    if [ "$livres" -gt 0 ] && [ "$passada" -lt 60 ]; then
        passada=$((passada + 1))
        echo "$(date '+%F %T') passada $passada: $livres realizacoes voltaram para a fila"
    elif [ "$livres" -eq 0 ] && [ "$acabando" -gt 0 ] && [ "$espera" -lt "$max_espera" ]; then
        [ "$espera" -eq 0 ] && echo "$(date '+%F %T') esperando $acabando realizacoes (de jobs que estao acabando ou trocando de dono)"
        espera=$((espera + 1))
        sleep "${ESPERA_S:-30}"
    else
        [ "$livres" -gt 0 ] && echo "$(date '+%F %T') $livres realizacoes ainda livres depois de 60 passadas: rode de novo" >&2
        [ "$acabando" -gt 0 ] && [ "$livres" -eq 0 ] &&
            echo "$(date '+%F %T') $acabando realizacoes de outros donos nao ficaram livres a tempo: rode de novo depois" >&2
        break
    fi
done
if [ "${SO_CAMPOS:-0}" != 1 ]; then
    printf '%s\n' "${E_D[@]}" "${E_DD[@]}" | sort -u | while IFS= read -r d; do
        if [ -d "$d" ]; then "$MC/junta_mc.sh" "$d"; fi
    done
fi
echo "$(date '+%F %T') $MC_PROJETO $CASOS: fim${SLURM_JOB_ID:+ do job $SLURM_JOB_ID}" | tee -a "$SAIDA_MC/execucoes.txt"
