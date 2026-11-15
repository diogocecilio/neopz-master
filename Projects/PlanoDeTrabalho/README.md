# Plano de trabalho

Este plano serve a dois artigos: o da amostragem por importância (AI) e o da chuva em solo não saturado. Ele valida as
ferramentas estatísticas antes de usá-las na chuva. Versão 1, de 15/11/2026; a versão final vem depois da revisão
independente.

## Princípios

- **Uma etapa por vez.** Cada etapa tem uma hipótese, um teste barato e um critério de aceite. Só se avança quando ela
  passa.
- **Tudo roda no fermi.** O PC só termina o que já estava rodando (etapa 0).
- **Modelo barato para validar métodos.** Todos os métodos usam o mesmo modelo numérico, então o viés de malha não
  atrapalha a comparação. A malha convergida fica para os números finais dos artigos.
- **Referência.** É o Monte Carlo direto no mesmo modelo, onde ele é barato (P_f de 1 a 20%). Com P_f pequena, a
  referência passa a ser a concordância entre dois métodos independentes e os valores publicados.
- **Problema de teste: o Cho2010.** É o talude seco do Cho (2010), com dois campos e o determinístico validado (1,2047
  contra 1,204), a ~80 s por realização no fermi.

## Etapas

| Etapa | O quê | Jobs no fermi | Custo | Pasta |
|---|---|---|---|---|
| 0 | Fechar o que já roda: a `chuva_ai` no fermi e o Cho2010 no PC (prévia, com o SRM antigo) | — | análise | `ETAPA0/` |
| 1 | Corrigir o SRM para começar abaixo de F = 1 e escolher a malha de teste; montar a infraestrutura do Cho2010 no fermi | 200 realizações pareadas antes e depois; determinísticos | ~15 min | `ETAPA1/` |
| 2 | FORM com M = 10…100 termos por campo e a AI com o ponto de projeto (M = 10, 20, 50, 100) | 1 job por M, 4 de cada vez; campanha na fila | ~8–10 h | `ETAPA2/` |
| 3 | AI com o deslocamento da entropia cruzada: tamanho do piloto, P_f menores (r = −0,7 e Wang 2020) e a comparação com a AI do FORM | campanhas na fila | ~5–6 h | `ETAPA3/` |
| 4 | Validar as ferramentas na chuva: f50 a 31,05 mm/h, com direto de 2000, AI e FORM com g = min_t FS − 1 | campanhas e FORM | ~1 dia | `ETAPA4/` |
| 5 | Investigar a chuva: P_f(t), curva de fragilidade P_f(I), solos, efeito do k e duração | campanhas | dias | `ETAPA5/` |

Fica parado até lá:
- a campanha grande do Ceron no fermi;
- os exemplos da literatura (van den Eijnden & Hicks 2017, Huang 2017), que são a etapa final do artigo da AI;
- a malha e o passo de tempo do rebaixamento;
- o ar ativo.

Cada pasta `ETAPAn/` guarda os resultados resumidos (`resultados/`) e o relatório da etapa. Os dados brutos ficam em
`~/projects/runs`.
