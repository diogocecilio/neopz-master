# Guia de estudo

Um resumo objetivo de tudo o que foi feito no talude não saturado sob chuva com campos aleatórios no NeoPZ:
formulação, verificação, fator de segurança, chuva, Karhunen-Loève, Monte Carlo e confiabilidade.

O guia está em `guia_estudo.pdf`, com 18 páginas. Ele resume os documentos do zip `formulacao.zip` (14/11) e os
projetos do branch `UnsaturatedSoil`, e liga cada documento ao código. Contém:
- o mapa da cadeia de trabalho;
- os números de referência;
- as equações essenciais de cada etapa;
- os resultados com números;
- as pendências por prioridade;
- um roteiro de leitura com 17 perguntas de estudo;
- a evolução da formulação;
- uma folha de fórmulas.

Para compilar:

    latexmk -pdf guia_estudo.tex

Cada seção fica em `secoes/`.
