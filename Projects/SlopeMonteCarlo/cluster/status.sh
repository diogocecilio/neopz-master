#!/bin/bash
# status.sh: os jobs slopeMC na fila, a ultima linha do log de cada processo e a tabela do Monte Carlo ate agora (junta os
# mc_resumo.csv e compara cada caso com o artigo). NREF e SAIDA como em submeter.sh: sem eles, a campanha de NREF = 4
# (SAIDA=<pasta> ./status.sh mostra outra, por exemplo a de 3 niveis em $DATA/ceron_mc).
cd "$(dirname "${BASH_SOURCE[0]}")"
NREF="${NREF:-4}"
SAIDA="${SAIDA:-${DATA:+$DATA/ceron_mc_h$NREF}}"
SAIDA="${SAIDA:-$PWD/resultados_h$NREF}"
command -v squeue > /dev/null && squeue -u "${USER:-$(id -un)}" -n slopeMC
if [ -d "$SAIDA/logs" ]; then
    echo "processos: $(ls "$SAIDA"/logs/p*.txt 2>/dev/null | wc -l); tarefas feitas: $(cat "$SAIDA"/logs/p*.txt 2>/dev/null | grep -c 'Gamma =')"
    for f in "$SAIDA"/logs/p*.txt; do [ -e "$f" ] && echo "  $(basename "$f" .txt): $(tail -n 1 "$f")"; done | sort -V | tail -n 5
    # tarefas que terminaram numa excecao e que nenhuma tentativa posterior concluiu: sem resultado
    erros="$(cat "$SAIDA"/logs/p*.txt 2>/dev/null | awk -F' [|] ' '
        / [|] ERRO: /   { k = $1; gsub(/ +/, " ", k); erro[k] = $2 }
        / [|] Gamma = / { k = $1; gsub(/ +/, " ", k); feita[k] = 1 }
        END { for (k in erro) if (!(k in feita)) print k " | " erro[k] }' | sort -V)"
    if [ -n "$erros" ]; then
        echo "ATENCAO: $(echo "$erros" | wc -l) tarefa(s) com erro, sem resultado (as ultimas; a proxima submissao tenta de novo):"
        echo "$erros" | tail -n 5 | sed 's/^/  /'
    fi
fi
if [ -s "$SAIDA/CONTINUACAO_PENDENTE" ]; then
    echo "ATENCAO: o SLURM recusou a continuacao de um job; quando os jobs acabarem, rode ./submeter.sh <nos>:"
    sed 's/^/  /' "$SAIDA/CONTINUACAO_PENDENTE"
fi
bin/SlopeMonteCarloCluster resumo "$SAIDA"
