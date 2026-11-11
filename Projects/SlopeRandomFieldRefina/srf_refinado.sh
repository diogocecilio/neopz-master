#!/bin/bash
# srf_refinado.sh: uma realizacao do SlopeRandomField em duas etapas de malha, sem mudar o programa: o fator na malha
# grossa, uma malha conforme refinada so em volta da banda de ruptura dessa realizacao e o fator de novo nela. E um EXE
# para ../SlopeRandomField/rodar_mc.sh (EXE=<caminho deste arquivo>), que o roda na pasta da realizacao com o ambiente do
# programa (SOLO, I_MM_H, CAMPOS, REAL, MESH e as variaveis do modelo); a realizacao esta concluida quando resultado.txt
# existe nessa pasta. As etapas:
#   1. grossa/: a realizacao na malha MESH (a grossa, com o seu .geo ao lado), com mecanismo.csv so da banda;
#   2. refina_banda.py: malha_refinada.geo e malha_refinada.msh, o .geo da malha grossa com elementos de REFINA_H ate
#      REFINA_FOLGA dos pontos da banda (sqrt(J2) >= REFINA_LIMIAR do maximo, no instante do fator) encontrada em 1. O
#      contorno nao muda: o retangulo do Karhunen-Loeve e o campo aleatorio da realizacao sao os mesmos nas duas malhas;
#   3. a mesma realizacao (mesmos REAL, CAMPOS e SEMENTE) na malha refinada, nesta pasta: resultado.txt, resumo.csv,
#      fs.csv, ... como os do programa. O seg de resultado.txt e resumo.csv e o desta etapa; o total esta em refino.csv.
# Sem banda na malha grossa (fator nao positivo, ou o transiente nao chegou ao instante do fator), o resultado e o da
# malha grossa, copiado para esta pasta (refino.csv com metodo=grossa).
# refino.csv (cabecalho e uma linha): o fator, o tempo e as equacoes na malha grossa, as celulas de h x h da banda, os
# triangulos da malha refinada e o tempo do gmsh (escrito antes da etapa 3) e, ao fim, o fator, o tempo e as equacoes na
# malha refinada e o tempo total das tres etapas.
# Retomavel: cada etapa pula o que ja esta feito, e uma rodada do programa interrompida continua do seu checkpoint.bin.
# Um sinal no programa ou no refina_banda.py (Ctrl-C, scancel, kill) sai com o codigo dele (128 + sinal), sem mensagem;
# rodar_mc.sh o trata como o do programa: 129 e 130 (terminal fechado, Ctrl-C) nao contam como falha.
# Depois do sucesso, apaga malha_refinada.msh (refeita do .geo, se a realizacao for rodada de novo), grossa/mecanismo.csv e
# um checkpoint que tenha ficado em grossa/ (de uma queda entre o resultado.txt e a limpeza do programa).
# Ambiente: SRF_BIN (padrao ../SlopeRandomField/bin/SlopeRandomField), REFINA_H (0.125), REFINA_FOLGA (0.375),
# REFINA_LIMIAR (0.1), PYTHON (python3; com o modulo python gmsh ou, sem ele, um executavel gmsh no PATH que rode). Os
# REFINA_* de uma combinacao ficam em ../refino_parametros.txt; outra rodada na mesma pasta de combinacao com outros
# valores e recusada.
set -u
export LC_ALL=C # numeros com ponto decimal em qualquer locale (awk, printf)
dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRF_BIN="${SRF_BIN:-$dir/../SlopeRandomField/bin/SlopeRandomField}"
py="${PYTHON:-python3}"
refina="$dir/refina_banda.py"
if ! [ -x "$SRF_BIN" ]; then
    echo "srf_refinado: executavel $SRF_BIN nao encontrado (compile o SlopeRandomField ou defina SRF_BIN)" >&2
    exit 2
fi
# rodar_mc.sh le a versao em versao.txt ao lado do EXE, que aqui e este script: a versao e a do programa
if [ -z "${SRF_VERSAO:-}" ] || [ "$SRF_VERSAO" = - ]; then
    SRF_VERSAO="$(cat "$(dirname "$SRF_BIN")/versao.txt" 2>/dev/null || echo -)"
fi
export SRF_VERSAO
# SO_CAMPOS=1 (rodar_mc.sh): so os campos, na malha grossa, sem as etapas
[ "${SO_CAMPOS:-0}" = 1 ] && exec "$SRF_BIN"

[ -n "${MESH:-}" ] || { echo "srf_refinado: defina MESH (a malha grossa)" >&2; exit 2; }
case "$MESH" in /*) malha_g="$MESH" ;; *) malha_g="$PWD/$MESH" ;; esac
if ! [ -f "$malha_g" ] || ! [ -f "${malha_g%.*}.geo" ]; then
    echo "srf_refinado: faltam a malha grossa $malha_g ou o seu .geo, ${malha_g%.*}.geo" >&2
    exit 2
fi
# num <valor> <minimo> <maximo>: o numero normalizado (%g), ou codigo 1 se nao for um numero nesse intervalo
num() { awk -v v="$1" -v lo="$2" -v hi="$3" 'BEGIN {
    if (v !~ /^[0-9]*[.]?[0-9]+([eE][-+]?[0-9]+)?$/ || v + 0 < lo || v + 0 > hi) exit 1; printf "%g", v }'; }
if ! h="$(num "${REFINA_H:-0.125}" 1e-6 1e6)" || ! folga="$(num "${REFINA_FOLGA:-0.375}" 0 1e6)" ||
    ! limiar="$(num "${REFINA_LIMIAR:-0.1}" 1e-9 1)"; then
    echo "srf_refinado: REFINA_H (> 0), REFINA_FOLGA (>= 0) ou REFINA_LIMIAR (entre 0 e 1) invalido" >&2
    exit 2
fi
# os parametros do refinamento sao os mesmos em toda a combinacao (a pasta acima desta), como os de parametros.txt
par="REFINA_H=$h REFINA_FOLGA=$folga REFINA_LIMIAR=$limiar"
arq_par="$(dirname "$PWD")/refino_parametros.txt"
if ! [ -f "$arq_par" ]; then
    echo "$par" > "$arq_par.$$" && mv -n "$arq_par.$$" "$arq_par" 2>/dev/null
    rm -f "$arq_par.$$"
fi
if [ "$(cat "$arq_par" 2>/dev/null)" != "$par" ]; then
    echo "srf_refinado: $(dirname "$PWD") ja tem realizacoes com $(cat "$arq_par" 2>/dev/null); agora $par." \
        "Use outra SAIDA." >&2
    exit 2
fi

# campo <csv> <coluna>: o valor da coluna num CSV de cabecalho e uma linha
campo() { awk -F, -v c="$2" 'NR == 1 { for (i = 1; i <= NF; i++) if ($i == c) k = i } NR == 2 && k { print $k }' "$1" 2>/dev/null; }
# chave <arquivo> <chave>: o valor de chave=valor numa linha de texto
chave() { tr ' ' '\n' < "$1" 2>/dev/null | sed -n "s/^$2=//p" | head -n 1; }
agora() { date +%s.%N; }
soma() { awk -v a="${1:-0}" -v b="${2:-0}" -v c="${3:-0}" 'BEGIN { printf "%.6g", a + b + c }'; }
cab="real,campos,malha_grossa,fator_grossa,seg_grossa,neq_grossa,celulas_banda,h,folga,limiar,triangulos,seg_malha"
cab="$cab,metodo,fator,seg_fina,neq,seg_total"
# escreve_refino: refino.csv com os valores das variaveis (as da etapa 3 vazias ate o fim), por troca atomica
escreve_refino() {
    {
        echo "$cab"
        echo "${REAL:-0},${CAMPOS:-1},$(basename "$malha_g"),$fator_g,$seg_g,$neq_g,$celulas,$h,$folga,$limiar,$tri,$seg_m,$metodo,${fator:-},${seg_f:-},${neq:-},${seg_t:-}"
    } > refino.csv.tmp && mv -f refino.csv.tmp refino.csv
}

# etapa2_feita: refino.csv com os mesmos parametros e, com banda, o .geo da malha refinada
etapa2_feita() {
    [ -s refino.csv ] && [ "$(campo refino.csv h),$(campo refino.csv folga),$(campo refino.csv limiar)" = "$h,$folga,$limiar" ] &&
        { [ "$(campo refino.csv metodo)" = grossa ] || [ -s malha_refinada.geo ]; }
}
# a banda da malha grossa (grossa/mecanismo.csv) e apagada ao fim; se a etapa 2 tiver de ser refeita sem ela, a etapa 1
# tambem e
if ! etapa2_feita && [ -s grossa/resultado.txt ] && ! [ -f grossa/mecanismo.csv ]; then rm -f grossa/resultado.txt; fi

# ---- etapa 1: a realizacao na malha grossa, em grossa/ (so a banda em mecanismo.csv; sem pw.csv, campos.csv e VTK)
mkdir -p grossa || exit 1
if ! [ -s grossa/resultado.txt ]; then
    # sem gmsh, a etapa 2 nao teria como rodar: o erro sai ja aqui
    if ! "$py" "$refina" --gmsh > /dev/null; then
        echo "srf_refinado: sem gmsh para $py: instale o modulo ($py -m pip install --user gmsh) ou ponha no PATH" \
            "um executavel gmsh que rode" >&2
        exit 2
    fi
    echo "$(date '+%F %T') etapa 1: malha grossa $(basename "$malha_g")"
    (cd grossa && exec env MESH="$malha_g" MEC_CSV=1 MEC_LIMIAR="$limiar" PW_CSV=0 CAMPOS_CSV=0 CAMPOS_VTK=0 VTK_H= \
        "$SRF_BIN" >> run.log 2>> run.err)
    rc=$?
    if ! [ -s grossa/resultado.txt ]; then
        [ $rc -ge 128 ] ||
            echo "srf_refinado: a rodada na malha grossa terminou com codigo $rc sem resultado.txt; veja grossa/run.err" >&2
        exit $((rc ? rc : 1))
    fi
fi
fator_g="$(chave grossa/resultado.txt FS_min)"
seg_g="$(chave grossa/resultado.txt seg)"
neq_g="$(chave grossa/realizacao.txt neq)"

# ---- etapa 2: a malha refinada em volta da banda (refino.csv marca a etapa como feita)
if etapa2_feita; then
    celulas="$(campo refino.csv celulas_banda)" tri="$(campo refino.csv triangulos)" seg_m="$(campo refino.csv seg_malha)"
    metodo="$(campo refino.csv metodo)"
    if [ "$metodo" = refinada ] && ! [ -s malha_refinada.msh ]; then
        # apagada ao fim de uma rodada completa (e resultado.txt apagado para rodar de novo): refeita do .geo
        "$py" "$refina" --malha malha_refinada > /dev/null || exit 1
    fi
else
    echo "$(date '+%F %T') etapa 2: malha refinada em volta da banda (fator $fator_g na malha grossa)"
    t0="$(agora)"
    info="$("$py" "$refina" grossa "$malha_g" malha_refinada "$h" "$folga" "$limiar")"
    rc=$?
    seg_m="$(awk -v a="$t0" -v b="$(agora)" 'BEGIN { printf "%.3f", b - a }')"
    if [ $rc -eq 0 ]; then
        metodo=refinada
        celulas="$(echo "$info" | tr ' ' '\n' | sed -n 's/^celulas_banda=//p')"
        tri="$(echo "$info" | tr ' ' '\n' | sed -n 's/^triangulos=//p')"
        echo "$info"
    elif [ $rc -eq 3 ]; then
        metodo=grossa celulas=0 tri=""
        echo "$(date '+%F %T') sem banda de ruptura na malha grossa: o resultado e o da malha grossa"
    else
        [ $rc -ge 128 ] || echo "srf_refinado: refina_banda.py terminou com codigo $rc" >&2
        exit $((rc ? rc : 1))
    fi
    escreve_refino
fi

# ---- etapa 3: a mesma realizacao na malha refinada, nesta pasta (ou, sem banda, o resultado da malha grossa)
if [ "$metodo" = refinada ]; then
    echo "$(date '+%F %T') etapa 3: malha refinada ($tri triangulos)"
    env MESH="$PWD/malha_refinada.msh" "$SRF_BIN"
    rc=$?
    if ! [ -s resultado.txt ]; then
        [ $rc -ge 128 ] ||
            echo "srf_refinado: a rodada na malha refinada terminou com codigo $rc sem resultado.txt" >&2
        exit $((rc ? rc : 1))
    fi
    fator="$(chave resultado.txt FS_min)" seg_f="$(chave resultado.txt seg)" neq="$(chave realizacao.txt neq)"
    seg_t="$(soma "$seg_g" "$seg_m" "$seg_f")"
    escreve_refino
else
    rc=0
    for f in resumo.csv fs.csv fs_tentativas.csv passos.csv realizacao.txt triangulos.csv; do
        if [ -f "grossa/$f" ]; then cp -f "grossa/$f" "$f.tmp" && mv -f "$f.tmp" "$f"; fi
    done
    fator="$fator_g" seg_f=0 neq="$neq_g" seg_t="$(soma "$seg_g" "$seg_m")"
    escreve_refino
    # resultado.txt por ultimo: marca a realizacao como concluida
    cp -f grossa/resultado.txt resultado.txt.tmp && mv -f resultado.txt.tmp resultado.txt
fi
rm -f malha_refinada.msh grossa/mecanismo.csv grossa/checkpoint.bin grossa/checkpoint.bin.tmp
echo "$(date '+%F %T') fim: fator $fator ($metodo), $fator_g na malha grossa; $seg_t s nas tres etapas"
exit $rc
