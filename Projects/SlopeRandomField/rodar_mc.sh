#!/bin/bash
# rodar_mc.sh <primeira> <ultima>: roda as realizacoes primeira..ultima do Monte Carlo, cada uma num processo de 1 CPU,
# PARALELO ao mesmo tempo (padrao: os nucleos fisicos da maquina; com hyper-threading, dois processos por nucleo deixam
# cada um ~2x mais lento). Uma combinacao (SOLO e I_MM_H, padrao f50 e 36) ou varias numa fila so, com
# COMBS="f10:9 f50:80 ..." (solo:chuva), em ordem intercalada (a realizacao 0 de todas, depois a 1, ...); DET=1 poe antes
# o caso deterministico (CAMPOS=0, as medias) de cada uma, e CAMPOS=0 roda so os deterministicos. Uma entrada pode trocar
# parametros so para ela (estudos de sensibilidade), em pasta propria: COMBS="f50:36:COV_C=0.5 f50:36:LX=40,LY=4
# f50:36:ALFA=2" (COV_C, COV_PHI, COV_K, LX, LY, ALFA, KL_M, F_TRANS; o deterministico so muda com ALFA e F_TRANS).
# F_TRANS (variantes acopladas i e ii): a resistencia do transiente e c'/F_TRANS, tan(phi')/F_TRANS. Com 1 (o padrao) o
# FS do SRM nao cai abaixo de 1 (o transiente esta sempre em equilibrio, e a dilatancia segura o talude no limite);
# F_TRANS=0.5 deixa o esqueleto elastico no transiente e o FS livre: COMBS="f50:36:F_TRANS=0.5" (ver README.md).
# Cada realizacao fica em $SAIDA/<combinacao>/r<nnnnn>/ (resumo.csv, fs.csv, passos.csv, ...; ver README.md); no fim,
# junta_mc.sh junta as concluidas em mc_resumo.csv e mc_fs.csv de cada combinacao. ./resumo_mc.sh <pasta> calcula P_f.
# Retomavel, e seguro com varios processos ou jobs (em nos diferentes) na mesma fila: cada realizacao em andamento tem um
# dono (r<nnnnn>/.ativo: maquina, processo e job do SLURM). As concluidas (resultado.txt) e as de um dono vivo sao
# puladas; a de um dono morto (job que acabou, processo que caiu, maquina reiniciada) continua do seu checkpoint.bin.
# Uma realizacao que termina com erro MAX_FALHAS vezes (padrao 3) fica marcada (falhou.txt) e e pulada. Com o arquivo
# $SAIDA/PARAR, nenhuma realizacao nova comeca (as em andamento terminam).
# rodar_mc.sh --faltam <primeira> <ultima>: so conta as realizacoes ainda nao concluidas nem marcadas como falhas;
# --pastas: as pastas de cada entrada de COMBS; --verifica: so confere COMBS e os parametros.txt das combinacoes.
# Avisos (parametros recusados, disjuntor) tambem vao para $SAIDA/AVISOS, que status_mc.sh mostra.
# Ambiente: VARIANTE (iii, ii, i), MALHA (rahardjo_a45.msh ou rahardjo_a45_fina.msh), FS_H, COV_C, COV_PHI, COV_K, LX, LY,
# ALFA, KL_M, SEMENTE, EXE, SAIDA, PARALELO, COMB (nome da pasta, com uma combinacao so). CSV opcionais: PW_CSV,
# MEC_CSV, CAMPOS_CSV (1 = grava, o padrao). VTK: CAMPOS_VTK=1 grava campos.vtk, VTK_H=24 o mecanismo de ruptura em
# t = 24 h; SO_CAMPOS=1 so gera os campos (campos.vtk e campos.csv) das realizacoes pedidas, sem simular.
set -u
export LC_ALL=C # numeros com ponto decimal em qualquer locale (awk, printf)
PKG="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODO=rodar
case "${1:-}" in --faltam) MODO=faltam; shift ;; --pastas) MODO=pastas; set -- 0 0 ;; --verifica) MODO=verifica; set -- 0 0 ;; esac
[ $# -eq 2 ] || { sed -n '2,/^set -u/{/^#/p}' "$0"; exit 2; }
export EXE="${EXE:-$PKG/bin/SlopeRandomField}"
# a versao do programa (gravada por compilar.sh), que vai para o resumo.csv de cada realizacao
export SRF_VERSAO="$(cat "$(dirname "$EXE")/versao.txt" 2>/dev/null || echo -)"
if [ "$MODO" = rodar ] && ! [ -x "$EXE" ]; then
    echo "executavel $EXE nao encontrado (compile o alvo SlopeRandomField; ver README.md)" >&2
    exit 2
fi
export MESH="${MESH:-$PKG/../LearningMohrCoulomb2/cluster/malhas/${MALHA:-rahardjo_a45.msh}}"
export VARIANTE="${VARIANTE:-iii}" MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 MAX_FALHAS="${MAX_FALHAS:-3}"
export SAIDA_MC="${SAIDA:-$PKG/resultados_mc}"
# entradas de COMBS: <solo>:<chuva>[:VAR=valor,...], separadas por espaco (ou por virgula antes do solo); VAR troca um
# parametro do programa so nessa entrada, que vai para a sua propria pasta (f50:36:COV_C=0.5 -> f50_I36_iii_covc0.5)
COMBS="$(echo ${COMBS:-${SOLO:-f50}:${I_MM_H:-36}} | sed -E 's/,(f10|f50|f100):/ \1:/g')"
VARS_MC='COV_C|COV_PHI|COV_K|LX|LY|ALFA|KL_M|F_TRANS'
NUM_MC='[0-9]+([.][0-9]+)?'
RE_MC="^(f10|f50|f100):$NUM_MC(:($VARS_MC)=$NUM_MC(,($VARS_MC)=$NUM_MC)*)?$"
for sc in $COMBS; do
    [[ "$sc" =~ $RE_MC ]] && continue
    echo "COMBS: \"$sc\" nao e <solo>:<chuva em mm/h>[:VAR=valor,...] com solo f10, f50 ou f100 e VAR entre" \
        "${VARS_MC//|/, } (exemplos: COMBS=\"f10:9 f50:80\", COMBS=\"f50:36:COV_C=0.5 f50:36:LX=40,LY=4\")" >&2
    exit 2
done
# COMB (nome da pasta) so vale com uma combinacao
COMB_UNICA=""
[ "$(echo $COMBS | wc -w)" -eq 1 ] && COMB_UNICA="${COMB:-}"
HOST_MC="$(hostname -s 2>/dev/null || uname -n)"
HOST_MC="${HOST_MC%%.*}"
export COMB_UNICA CAMPOS_TODOS="${CAMPOS:-1}" HOST_MC="${HOST_MC:-sem_nome}"
BOOT_MC="$(cat /proc/sys/kernel/random/boot_id 2>/dev/null || echo -)"
export BOOT_MC
[ "$MODO" = rodar ] && mkdir -p "$SAIDA_MC"

# partes <entrada>: s (solo), i (chuva), x (as trocas VAR=valor, ou -) e xd (as que mudam o deterministico: ALFA e
# F_TRANS)
partes() {
    local r="${1#*:}" e d=""
    s="${1%%:*}" i="${r%%:*}" x=- xd=-
    [ "$r" != "$i" ] && x="${r#*:}"
    for e in ${x//,/ }; do case "${e%%=*}" in ALFA | F_TRANS) d="${d:+$d,}$e" ;; esac; done
    xd="${d:--}"
}
# pasta <solo> <chuva> <campos> [<trocas>]: a pasta da combinacao (o deterministico em <combinacao>_det)
pasta() {
    local c="${COMB_UNICA:-${1}_I${2}_${VARIANTE}}" e
    if [ "${4:--}" != - ]; then
        for e in ${4//,/ }; do c="${c}_$(echo "${e%%=*}" | tr -d _ | tr 'A-Z' 'a-z')${e#*=}"; done
    fi
    if [ "$3" = 0 ] && ! { [ "$CAMPOS_TODOS" = 0 ] && [ -n "$COMB_UNICA" ]; }; then c="${c}_det"; fi
    echo "$SAIDA_MC/$c"
}
# tarefas <primeira> <ultima>: uma linha por realizacao, "<solo> <chuva> <campos> <n> <trocas>"; os deterministicos,
# um por solo, chuva e ALFA (as outras trocas nao mudam o caso com as medias)
tarefas() {
    local sc n s i x xd vistos=" "
    if [ "$CAMPOS_TODOS" = 0 ]; then
        for n in $(seq "$1" "$2"); do for sc in $COMBS; do partes "$sc"; echo "$s $i 0 $n $xd"; done; done
        return
    fi
    if [ "${DET:-0}" = 1 ]; then
        for sc in $COMBS; do
            partes "$sc"
            case "$vistos" in *" $s:$i:$xd "*) continue ;; esac
            vistos="$vistos$s:$i:$xd "
            echo "$s $i 0 0 $xd"
        done
    fi
    for n in $(seq "$1" "$2"); do for sc in $COMBS; do partes "$sc"; echo "$s $i 1 $n $x"; done; done
}
# pendentes <primeira> <ultima>: as tarefas sem resultado.txt nem falhou.txt
pendentes() {
    local s i c n x r
    tarefas "$1" "$2" | while read -r s i c n x; do
        r="$(pasta "$s" "$i" "$c" "$x")/$(printf 'r%05d' "$n")"
        [ -s "$r/resultado.txt" ] || [ -e "$r/falhou.txt" ] || echo "$s $i $c $n $x"
    done
}
aviso_mc() { echo "$(date '+%F %T') $HOST_MC${SLURM_JOB_ID:+ job $SLURM_JOB_ID}: $*" >> "$SAIDA_MC/AVISOS"; }
parar() {
    [ -e "$SAIDA_MC/PARAR" ] || [ -e "$INST_MC/parar" ] || { [ -n "${SLURM_JOB_ID:-}" ] && [ -e "$SAIDA_MC/.parar.$SLURM_JOB_ID" ]; }
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
    # outra maquina, fora do SLURM: viva enquanto a realizacao muda (passos.csv, a cada passo)
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
# uma <solo> <chuva> <campos> <n> <trocas>: uma realizacao, na sua pasta; a saida do programa vai para run.log e run.err
uma() {
    local solo="$1" imm="$2" campos="$3" n="$4" x="${5:--}" d r rc nf nome
    local -a ev=()
    [ "$x" != - ] && ev=(${x//,/ })
    d="$(pasta "$solo" "$imm" "$campos" "$x")"
    r="$d/$(printf 'r%05d' "$n")"
    nome="$(basename "$d")/$(basename "$r")"
    if [ "${SO_CAMPOS:-0}" = 1 ]; then # so os campos desta realizacao, em VTK e CSV, sem simulacao
        [ "$campos" = 1 ] || return 0
        mkdir -p "$r" && cd "$r" || return 1
        if env ${ev[@]+"${ev[@]}"} SOLO="$solo" I_MM_H="$imm" CAMPOS=1 REAL="$n" CAMPOS_VTK=1 CAMPOS_CSV=1 SO_CAMPOS=1 \
            "$EXE" > campos.log 2>&1; then
            echo "$nome: campos em $r/campos.vtk"
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
    (cd "$r" && env ${ev[@]+"${ev[@]}"} SOLO="$solo" I_MM_H="$imm" CAMPOS="$campos" REAL="$n" "$EXE" >> run.log 2>> run.err)
    rc=$?
    if [ -s "$r/resultado.txt" ]; then
        rm -f "$r/checkpoint.bin" "$r/checkpoint.bin.tmp" # sobra de uma queda entre o resultado e a limpeza
        echo 0 > "$INST_MC/seguidas"
        echo "$(date '+%F %T') $nome: $(cat "$r/resultado.txt")"
    elif [ $rc -eq 129 ] || [ $rc -eq 130 ]; then
        # Ctrl-C ou terminal fechado (o trabalhador sobreviveu): interrupcao, nao falha
        echo "$(date '+%F %T') $nome: interrompida (sinal $((rc - 128))); continua do checkpoint" >&2
    else
        # falha: inclusive um SIGKILL/SIGTERM so no programa (o limite de tempo, o scancel e o Ctrl-C matam tambem o
        # trabalhador, que entao nem chega aqui): em geral o OOM killer
        nf=$(($(cat "$r/falhas" 2>/dev/null || echo 0) + 1))
        echo "$nf" > "$r/falhas"
        dica=""
        [ $rc -eq 137 ] && dica=" (SIGKILL: falta de memoria? veja memoria_MB nos resumo.csv)"
        if [ "$nf" -ge "$MAX_FALHAS" ]; then
            echo "$(date '+%F %T') $HOST_MC codigo $rc em $nf tentativas$dica" > "$r/falhou.txt"
            echo "$(date '+%F %T') $nome: falhou $nf vezes (codigo $rc$dica), marcada em falhou.txt; veja $r/run.err" >&2
        else
            echo "$(date '+%F %T') $nome: terminou com codigo $rc$dica sem resultado.txt (tentativa $nf de $MAX_FALHAS)" >&2
        fi
        # disjuntor: falhas seguidas demais nesta instancia, sem nenhum sucesso no meio, indicam um problema do ambiente
        # (executavel, malha, bibliotecas, disco) e nao das realizacoes: para de pegar realizacoes antes de marcar a fila
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
export -f aviso_mc partes pasta parar job_vivo pid_vivo ativa dono_vivo escreve_dono libera reivindica uma
# disponiveis <primeira> <ultima>: "<livres> <esperando>": as tarefas pendentes sem dono ou com dono morto (que uma nova
# passada pode rodar) e as que ficam livres em breve (jobs acabando ou que acabaram ha pouco, troca de dono em andamento,
# SLURM sem responder)
disponiveis() {
    local s i c n x r e
    pendentes "$1" "$2" | while read -r s i c n x; do
        r="$(pasta "$s" "$i" "$c" "$x")/$(printf 'r%05d' "$n")"
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
if [ "$MODO" = pastas ]; then # "<entrada> <pasta> <pasta do deterministico ou ->", para status_mc.sh
    for sc in $COMBS; do
        partes "$sc"
        dd=-
        [ "${DET:-0}" = 1 ] && dd="$(pasta "$s" "$i" 0 "$xd")"
        echo "$sc $(pasta "$s" "$i" 1 "$x") $dd"
    done
    exit 0
fi
NUC="$(lscpu -p=CORE,SOCKET 2>/dev/null | grep -v '^#' | sort -u | wc -l)"
P="${PARALELO:-${NUC:-1}}"
[ "$P" -ge 1 ] 2>/dev/null || P=1
# estado desta instancia (as falhas seguidas, para o disjuntor), apagado no fim
INST_MC="$(mktemp -d)"
export INST_MC DISJUNTOR_MC="${DISJUNTOR:-$((2 * P > MAX_FALHAS + 2 ? 2 * P : MAX_FALHAS + 2))}"
trap 'rm -rf "$INST_MC"' EXIT
# parametros.txt de cada combinacao: os parametros do programa (fora solo, chuva e variante, que estao no nome da pasta)
# com que as suas realizacoes foram rodadas. Outra rodada na mesma pasta com outros parametros misturaria resultados:
# recusada (use outra SAIDA, ou uma entrada de COMBS com as trocas)
VARS_PARAM="MESH TCHUVA_H DT_H FS_H TF_H XL XR HW BETA_DEG SUCMAX XTOE YTOE HL HR CHUVA_VERTICAL FS_TOL FS_F0 MAX_CORTES E NU
COES PHI_DEG GAMMA KS POROS SEMENTE COV_C COV_PHI COV_K LX LY ALFA KL_M KL_FRAC SOLVER TOLU TOLP MAXIT MAXIT_FROZEN F_TRANS"
# variaveis acrescentadas depois: so entram em parametros.txt quando definidas, para as pastas de antes continuarem valendo
VARS_PARAM_NOVAS=" F_TRANS "
parametros() { # parametros <trocas>: uma linha por variavel, VAR=valor (- se nao definida, o padrao do programa)
    local v val e
    for v in $VARS_PARAM; do
        val="${!v:--}"
        for e in ${1//,/ }; do [ "${e%%=*}" = "$v" ] && val="${e#*=}"; done
        [ "$v" = MESH ] && val="${val##*/}"
        [ "$val" = - ] && [[ "$VARS_PARAM_NOVAS" == *" $v "* ]] && continue
        echo "$v=$val"
    done
}
if [ "${SO_CAMPOS:-0}" != 1 ] && { [ "$MODO" = rodar ] || [ "$MODO" = verifica ]; }; then
    mkdir -p "$SAIDA_MC"
    for sc in $COMBS; do
        partes "$sc"
        for c in 1 0; do
            [ $c = 0 ] && [ "${DET:-0}" != 1 ] && [ "$CAMPOS_TODOS" != 0 ] && continue
            [ $c = 1 ] && [ "$CAMPOS_TODOS" = 0 ] && continue
            [ $c = 0 ] && xx="$xd" || xx="$x"
            d="$(pasta "$s" "$i" "$c" "$xx")"
            novo="$(parametros "$xx")"
            if [ -f "$d/parametros.txt" ]; then
                if [ "$(cat "$d/parametros.txt")" != "$novo" ]; then
                    dif="$(diff <(cat "$d/parametros.txt") <(echo "$novo") | grep '^[<>]' | sed 's/^</antes:/; s/^>/agora:/' | tr '\n' ' ')"
                    echo "$d ja tem realizacoes rodadas com outros parametros: $dif" >&2
                    echo "Use outra SAIDA (ou uma entrada de COMBS com as trocas, como f50:36:COV_C=0.5)." >&2
                    [ "$MODO" = rodar ] && aviso_mc "recusado: $(basename "$d") ja tem realizacoes com outros parametros: $dif"
                    exit 2
                fi
            else
                mkdir -p "$d" && echo "$novo" > "$d/parametros.txt"
            fi
        done
    done
fi
[ "$MODO" = verifica ] && exit 0
echo "$(date '+%F %T') $(echo $COMBS)$([ "${DET:-0}" = 1 ] && echo ' (e deterministicos)'): realizacoes $1..$2, $P ao" \
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
    for sc in $COMBS; do
        partes "$sc"
        for d in "$(pasta "$s" "$i" 1 "$x")" "$(pasta "$s" "$i" 0 "$xd")"; do
            [ -d "$d" ] && "$PKG/junta_mc.sh" "$d"
        done
    done
fi
echo "$(date '+%F %T') $(echo $COMBS): fim${SLURM_JOB_ID:+ do job $SLURM_JOB_ID}" | tee -a "$SAIDA_MC/execucoes.txt"
