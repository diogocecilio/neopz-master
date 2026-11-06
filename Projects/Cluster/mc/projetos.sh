# projetos.sh: os projetos de Monte Carlo da fila (rodar_mc.sh, job_mc.sbatch, submeter_mc.sh, status_mc.sh). Uso:
# source projetos.sh; mc_projeto <V2|V2_TaludesCho|V4|V4_TaludeLiu>.
# Cada executavel segue o contrato do Monte Carlo: uma realizacao por execucao, na pasta corrente, com REAL=<n>,
# CAMPOS=0|1 (0 = o caso deterministico, com as medias), SEMENTE, CASO=<nome do exemplo> e as trocas de parametros
# (VAR=valor); grava resumo.csv (cabecalho e uma linha: real, FS_min, ...), fs.csv (t_h, FS_SRM, ...) e, por ultimo,
# resultado.txt (um .tmp renomeado: "real=... FS_min=... situacao=ok seg=..."). Codigo de saida 0 em sucesso; com
# resultado.txt a realizacao esta concluida mesmo com outro codigo (situacao != ok: ver SITUACOES_RUPTURA).
# mc_projeto define:
#   MC_PROJETO       o nome do projeto (a pasta em Projects/ e a dos resultados em $SAIDA)
#   MC_ID            o prefixo (V2): os jobs se chamam verif-mc-<MC_ID>
#   MC_EXE           o executavel (em Projects/Cluster/bin)
#   MC_CASOS_RE      os valores de CASO aceitos (expressao regular)
#   MC_CASOS_PADRAO, MC_PRIMEIRA, MC_ULTIMA   a campanha padrao (a do artigo)
#   MC_VARS          os parametros do programa que mudam o resultado: as trocas aceitas em CASOS e o conteudo de
#                    parametros.txt (que recusa misturar realizacoes rodadas com valores diferentes na mesma pasta)
#   MC_VARS_CAMPO    os que so mudam os campos aleatorios: ficam fora do deterministico (CAMPOS=0 os ignora), entao
#                    cho1:COV_C=0.5 usa o mesmo deterministico que cho1
#   MC_VARS_CONTROLE as outras variaveis que o programa le (o contrato e as saidas opcionais): nao mudam o resultado
#   MC_VAR_VERSAO    a variavel com que o programa grava a versao (bin/versao.txt) no resumo.csv
#   MC_SO_CAMPOS     1 se o programa tem SO_CAMPOS=1 (so grava os campos de uma realizacao, sem simular)
# e a funcao mc_extra <caso> <campos> <pasta do deterministico>: variaveis a mais de uma execucao (VAR=valor, uma por
# linha; vazio: nenhuma). Para acrescentar um projeto: um bloco aqui, com as variaveis lidas pelo main.cpp/problema.h.

mc_projetos() { echo "V2_TaludesCho V4_TaludeLiu"; }

mc_extra() { :; }

mc_projeto() {
    case "${1:-}" in
    V2 | v2 | V2_TaludesCho)
        MC_PROJETO=V2_TaludesCho MC_ID=V2 MC_EXE=V2_TaludesCho
        # Cho (2010) / Liao et al. (2025): cho1 (argila nao drenada, campo de c_u) e cho2 (c-phi, campos de c e phi)
        MC_CASOS_RE='^(cho1|cho2)$'
        MC_CASOS_PADRAO="cho1 cho2"
        MC_PRIMEIRA=1 MC_ULTIMA=10000
        MC_VARS="NUCLEO SEMENTE COV_C COV_PHI LX LY RHO KL_FRAC KL_M KL_MARGEM H_MALHA ORDEM FS_TOL FS_F0 FS_MAX
            MAXIT_FROZEN COES PHI_DEG E NU GAMMA"
        MC_VARS_CAMPO="NUCLEO SEMENTE COV_C COV_PHI LX LY RHO KL_FRAC KL_M KL_MARGEM"
        MC_VARS_CONTROLE="CASO REAL CAMPOS VERBOSO CAMPOS_CSV SO_CAMPOS V2_VERSAO"
        MC_VAR_VERSAO=V2_VERSAO MC_SO_CAMPOS=1
        mc_extra() { :; }
        ;;
    V4 | v4 | V4_TaludeLiu)
        MC_PROJETO=V4_TaludeLiu MC_ID=V4 MC_EXE=V4_TaludeLiu
        # Liu et al. (2024): s1 (campo de k_s), s2 (c' e phi'), s3 (os tres); ANGULO=90 e o acamamento vertical
        MC_CASOS_RE='^(s1|s2|s3)$'
        MC_CASOS_PADRAO="s1 s2 s3 s1:ANGULO=90 s2:ANGULO=90 s3:ANGULO=90"
        MC_PRIMEIRA=1 MC_ULTIMA=2000
        MC_VARS="ANGULO SEMENTE COES PHI_DEG PHIB_DEG GAMMA E NU KS POROS A_FX N_FX M_FX P_KR Q_KS H_ESQ H_DIR
            CHUVA_VERTICAL DT0_H FATOR_DT DT_MAX_H DT_PERM_H TOL_PERM T_MAX_H MAX_CORTES SUC_INI COV_C COV_PHI COV_K
            THETA1 THETA2 RHO LX LY KL_M KL_FRAC FS_F0 FS_TOL TOLU TOLP MAXIT MAXIT_FROZEN MASSA_LIMIAR MALHA SOLVER"
        MC_VARS_CAMPO="ANGULO SEMENTE COV_C COV_PHI COV_K THETA1 THETA2 RHO LX LY KL_M KL_FRAC"
        MC_VARS_CONTROLE="CASO REAL CAMPOS VERBOSO CAMPOS_CSV PW_CSV MEC_CSV VTK PW_PERMANENTE V4_VERSAO"
        MC_VAR_VERSAO=V4_VERSAO MC_SO_CAMPOS=0
        # s2 (k_s constante): o regime permanente e o do deterministico, com os mesmos parametros hidraulicos (as
        # trocas do deterministico sao as da entrada sem as dos campos): depois que ele termina, as realizacoes leem o
        # pw.csv dele (PW_PERMANENTE) em vez de refazer o transiente. O deterministico sempre roda o transiente e grava
        # o pw.csv. Antes dele terminar, ou com outro caso, nada muda (o programa roda o transiente).
        mc_extra() {
            local r="$3/r00000"
            if [ "$2" = 0 ]; then
                printf '%s\n' "PW_PERMANENTE=" "PW_CSV=1"
            elif [ "$1" = s2 ] && [ -s "$r/pw.csv" ] && [ -s "$r/resultado.txt" ] &&
                grep -q 'situacao=ok' "$r/resultado.txt"; then
                echo "PW_PERMANENTE=$r/pw.csv"
            fi
        }
        ;;
    *)
        echo "PROJETO=${1:-} nao e um projeto de Monte Carlo (os de Monte Carlo: $(mc_projetos); V2 e V4 bastam)" >&2
        return 1
        ;;
    esac
    MC_VARS="$(echo $MC_VARS)" # uma linha so
    export MC_PROJETO MC_ID MC_EXE MC_CASOS_RE MC_CASOS_PADRAO MC_PRIMEIRA MC_ULTIMA MC_VARS MC_VARS_CAMPO \
        MC_VARS_CONTROLE MC_VAR_VERSAO MC_SO_CAMPOS
    export -f mc_extra
}

# mc_confere_vars <pasta do projeto>: as variaveis de ambiente que o main.cpp e os cabecalhos do projeto leem (os
# testes.cpp ficam de fora: as variaveis deles nao sao do executavel) e que nao estao em MC_VARS nem em
# MC_VARS_CONTROLE (uma lista desatualizada aqui deixaria um parametro novo fora de parametros.txt); vazio se confere
mc_confere_vars() {
    local v
    for v in $(cat "$1"/main.cpp "$1"/*.h 2>/dev/null | grep -o -E '(Env[SIDB]|getenv)\("[A-Z][A-Z0-9_]*"' |
        sed -E 's/.*\("//; s/"$//' | sort -u); do
        [[ " $MC_VARS $MC_VARS_CONTROLE " == *" $v "* ]] || echo "$v"
    done
}
