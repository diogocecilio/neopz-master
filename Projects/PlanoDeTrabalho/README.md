# Plano de trabalho

Este plano serve a dois artigos: o da amostragem por importância (AI) e o da chuva em solo não saturado. Ele valida as
ferramentas estatísticas antes de usá-las na chuva. É a versão 2, de 15/11/2026. Ela junta a revisão independente (uma
crítica estatística e uma de custo) e a decisão de rodar tudo no fermi.

## Princípios

- **Uma etapa por vez.** Cada etapa tem uma hipótese, um teste barato e um critério de aceite. Só se avança quando ela
  passa.
- **Tudo roda no fermi.** O PC parou em 15/11 às 11:33. O que ele tinha feito fica como prévia, com o SRM antigo.
- **Mesmo modelo em cada comparação:** a mesma versão do código, o mesmo SRM e a mesma malha. Um caso novo nunca é
  comparado com um direto de outra versão.
- **Modelo barato para validar métodos.** A malha convergida fica para os números finais dos artigos.
- **Referência:**
  - onde ele é barato (P_f de 1 a 20%), o Monte Carlo direto no mesmo modelo;
  - com P_f pequena, a concordância entre métodos independentes;
  - o valor publicado serve só como plausibilidade (mesma ordem de grandeza): ele vem de outro modelo (equilíbrio
    limite).
- **Problema de teste: o Cho2010.** É o talude seco do Cho (2010), com dois campos e o determinístico validado (1,2047
  contra 1,204).

## Etapas

| Etapa | O quê | No fermi | Situação | Pasta |
|---|---|---|---|---|
| 0 | O Cho2010 que rodava no PC: direto, AI pela banda e AI pelo ponto de projeto do FORM M = 10 | campanha `etapa1` | rodando | `ETAPA0/` |
| 1 | SRM corrigido (`srm_f0` = 0,8): o antigo contra o novo nas mesmas realizações, no Cho2010 e na chuva | campanha `etapa1` | rodando | `ETAPA1/` |
| 2 | FORM com M = 10, 20, 50 e 100 termos por campo; AI pelo ponto de projeto | 4 FORM em paralelo e uma campanha | depois da 1 | `ETAPA2/` |
| 3 | AI pela entropia cruzada: tamanho do piloto e P_f menor (r = −0,7) | uma campanha | depois da 2 | `ETAPA3/` |
| 4 | Validar FORM e AI na chuva (f50, FS em 24 h) | código novo e uma campanha | depois da 3 | `ETAPA4/` |
| 5 | Investigar a chuva: P_f(t), fragilidade P_f(I), solos, k e duração | campanhas | depois da 4 | `ETAPA5/` |

As etapas 0 e 1 rodam juntas na campanha `etapa1` (`SlopeRandomFieldRefinement/cluster/casos_etapa1.txt`), com ~290 h de
núcleo, ou 6 a 7 h em 4 nós. O primeiro nó roda antes os FORM com M = 10 e 20 (`form.sh`).

### Etapa 0: o Cho2010 com o SRM novo

- **Casos:**
  - direto com r = −0,5 (4000) e r = −0,7 (1000);
  - AI pela banda (1000 cada);
  - AI pelo ponto de projeto do FORM M = 10 (1000);
  - FORM com M = 10 e 20.
- **Critérios:**
  - a AI pela banda e a AI pelo FORM concordam entre si e, juntas, com o direto (teste das diferenças, uma vez só);
  - o FORM fica perto da AI do mesmo modelo: |ln(P_f,FORM / P_f,AI)| < ln 1,3.
- **Prévia do PC (SRM antigo):** 600 realizações diretas com r = −0,5 dão P_f = 1,67% (10 falhas). O FORM M = 10 deu
  β ≈ 2,13 e o M = 20, β ≈ 2,01, ainda sem convergir.

### Etapa 1: o SRM corrigido

- **Hipótese:** com F = 1 no peso, um talude com FS logo abaixo de 1 às vezes convergia, e o FS saía ≥ 1 (o buraco
  logo abaixo de 1). Com o caminho de F começando em 0,8, o FS logo abaixo de 1 é achado.
- **O que já se sabe:**
  - **Os solos (SRM antigo, todos os instantes):**
    - o f50 tem 13 FS em [0,995; 1,000), 47 em [1,000; 1,005) e nenhum em [0,999; 1,000);
    - o f10 tem 26 e 41.

    Numa faixa tão estreita, a densidade devia ser quase plana. O SRM antigo empurrava FS logo abaixo de 1 para logo
    acima, e a P_f saía baixa.
  - **Quatro realizações da chuva no PC** (f50, 24 h), com as mesmas opções nas duas versões:
    - três passaram de 1,0001–1,0010 (antigo) para 0,9985–0,9998 (novo);
    - a outra deu 0,9989 nas duas. No antigo, o peso não convergiu com F = 1, e ele recomeçou de F = 0,5.
  - **O fator ficou mais caro perto de 1:** 54 a 71 s no novo, contra 23 a 27 s no antigo, cujo caminho parava logo
    em F = 1. No determinístico do Cho2010, o custo subiu 16% e o FS mudou de 1,2047 para 1,2044. O fermi mede o custo
    nas realizações pareadas (o `compara_pares.py` mostra o tempo do fator).
- **No fermi**, pareado, só o `srm_f0` muda:
  - Cho2010: as 1000 primeiras (`r05_antigo` contra `r05_direto`);
  - chuva: 600 com o FS em 24 h (`chuva_antigo` contra `chuva_novo`).
- **Critérios:**
  - longe de 1 (|FS − 1| > 0,03), a diferença é só ruído (< 0,002);
  - o histograma novo não tem a faixa vazia em [0,999; 1,000);
  - os pares discordantes (McNemar) medem o viés antigo na P_f.
- **Malha de teste: sem rodadas extras.** A coluna `niveis` de cada realização tri2 + 4 já traz o FS dos níveis 0 a 3.
  O tri2 + 3 são os 4 primeiros níveis do tri2 + 4. Na prévia, o FS3 − FS4 = 0,005 ± 0,002, e a P_f é a mesma. A
  malha barata só serve para o FORM: o ponto de projeto é só a densidade de sorteio, e a AI continua sem viés.
- **Conferências de graça:**
  - `r05_antigo` (fermi) contra o direto do PC (SRM antigo) nas 600 em comum: as duas máquinas dão o mesmo FS;
  - `chuva_antigo` contra a campanha dos solos em 24 h: a lista de instantes não muda o FS.

### Etapa 2: FORM e AI pelo ponto de projeto

- **Antes, um teste sem elementos finitos.** Trocar o FS por g = β − α·u e sortear 10⁶ amostras com o próprio
  `RandomField` (`SetShift`, `LogWeight`). O teste exige P_f = Φ(−β) com β = 2, 3 e 4. Ele confere a semente, a ordem
  c/φ, os K primeiros termos e o peso em segundos, também com P_f pequena.
- **FORM independentes, saindo da origem, um por job:** M = 10, 20, 50 e 100. Os outros M só se β(M) não estabilizar.
  - Parada: |g| < 0,001 e mudança de β < 0,01.
  - Teste: a KL de K termos são os K primeiros termos da de 2000, então o β exato não pode crescer com M. Se crescer,
    é falta de convergência ou ruído.
- **AI:** com o ponto de projeto de M = 10 e com o do maior M, com a mesma semente (números aleatórios comuns).
  - Métrica: o custo total para um CoV alvo, com o FORM incluído: N_FORM + N_AI (CoV/CoV_alvo)².
  - Previsão a testar: Var(M)/Var(M') ≈ exp(β_M² − β_M'²).

### Etapa 3: AI pela entropia cruzada

- **Tamanho do piloto, de graça:** sub-pilotos disjuntos do direto de 4000 (40 de 100, 16 de 250, 8 de 500, 4 de
  1000). Só o deslocamento é recalculado. Depois, uma AI só, com o menor piloto que serve.
- **P_f pequena (r = −0,7, ~0,4%):**
  - estender o direto para ~6000 (~25 falhas, CoV ≈ 20%);
  - FORM M = 20 e a AI pelo ponto de projeto, contra a AI pela banda.
- **Fica para depois:** o Wang 2020 (precisa da média de c' como opção e da entropia cruzada em dois níveis).

### Etapa 4: a chuva

- **Código antes:**
  - levar `M=`, `xi=` e `desloc_xi=` do `cho2010.cpp` para o `montecarlo.cpp`;
  - o `form_kl.py` com 3 campos e o menor FS dos instantes.
- **Validação:** o f50 com o FS em 24 h, onde o mínimo cai em quase todas as realizações. As 200 do piloto do
  `chuva_ai` são as 200 primeiras do direto. O FORM só em c' e φ': o k deu α² ≈ 0, e a AI continua sem viés com o k
  sem deslocamento.
- **O FORM por instante** sai das mesmas análises. O teste do mínimo no tempo fica para a etapa 5, num caso em que o
  instante crítico varia (no f10 ele cai em 12, 18 e 24 h).

## Critérios gerais (da revisão)

- **O teste "< 2 desvios" contra o direto só pega erros grosseiros.** Com P_f ≈ 2% e 4000 diretas (CoV ≈ 11%), um viés
  de 25 a 38% passa despercebido metade das vezes. O relatório de cada etapa diz o viés mínimo detectável.
- **Com menos de ~20 falhas no direto,** usar o intervalo exato (Clopper–Pearson).
- **Diagnósticos da AI que medem a P_f:**
  - ESS_F = (Σ w I)² / Σ (w I)², que dá o CoV;
  - a maior parcela, max(w I) / Σ(w I), com alerta acima de 10%;
  - a fração de falhas sorteadas.
  - O ESS de Kong (sobre todos os pesos) não serve para escolher o método: com |μ|² ≈ 4,5, ele dá ~11 em 1000.
- **Outro modo de ruptura,** o erro que a AI com um centro só não mostra no CoV. Três verificações:
  - comparar a banda das falhas (colunas `banda_*`) entre o direto e a AI ponderada;
  - nas falhas do direto, t = α·u: se mais de ~10% tiverem t < 1, há outro modo;
  - o cosseno entre a média dos ξ das falhas e α tem de ser > 0,9.

## Fica parado até lá

- a campanha grande do Ceron no fermi;
- os exemplos da literatura (van den Eijnden & Hicks 2017, Huang 2017, Wang 2020), a etapa final do artigo da AI;
- a malha e o passo de tempo do rebaixamento;
- o ar ativo.

Cada pasta `ETAPAn/` guarda os resultados resumidos (`resultados/`) e o relatório da etapa. Os dados brutos ficam em
`~/projects/runs` (os do fermi em `~/projects/runs/fermi_<campanha>`).
