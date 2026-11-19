# Plano de trabalho

Este plano serve a dois artigos: o da amostragem por importância (AI) e o da chuva em solo não saturado. Ele valida as
ferramentas estatísticas antes de usá-las na chuva. É a versão 2, de 15/11/2026. Ela junta a revisão independente (uma
crítica estatística e uma de custo) e a decisão de rodar tudo no fermi. O índice de todos os documentos (formulação e
resultados) está no fim, na seção [Documentos](#documentos).

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
| 2b | Base de Nyström contra a analítica: β(M), ganho da AI e direto de referência próprio | uma campanha e FORM | junto com a 2 | `ETAPA2/` |
| 2c | Gradiente do FORM pela projeção modal (sensibilidade no colapso do SRM), validado contra as diferenças centrais | FORM | junto com a 2 | `ETAPA2/` |
| 3 | AI pela entropia cruzada: tamanho do piloto e P_f menor (r = −0,7) | uma campanha | depois da 2 | `ETAPA3/` |
| 3b | ISMC até CoV de 5% no campo de Nyström com 2000 modos, com o centro de um FORM com poucos termos (M = 20 e 100 por diferenças centrais; M = 20 pelo GIM modal) | campanha `etapa3b` | pronta para rodar | `ETAPAS1a3/` |
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

### Etapa 2b: base de Nyström contra a analítica

- **Hipótese:** a KL de Nyström do próprio talude precisa de menos termos que a base do retângulo para a mesma variância.
  Com isso, ela deixa o FORM mais barato e faz o deslocamento pelos primeiros modos funcionar.
  - No Mathematica, com a falha simplificada, os 10 primeiros modos de Nyström deram ganho de 30×.
  - No Ceron, com a base do retângulo, deram 1,2×.
- **Protótipo no Mathematica, já feito:**
  - `docs/confiabilidade_campos/ai_nystrom`: a AI com Nyström e uma falha simplificada;
  - `~/Dropbox/UnsaturatedSoil/mathematica/talude_nystrom/talude_nystrom.wl`: o talude do `poro_plastic.nb` (o do Cho)
    com c′ e φ′ por ponto de integração e o SRM com bissecção. Dá FS determinístico de 1,258 e ~1 min por
    realização.
- **Pontos abertos do protótipo:**
  - com M = 60, só 64% da variância é retida: subir M ou restringir a KL à região do mecanismo;
  - a malha é grossa: FS 1,258 contra 1,205 no NeoPZ.
- **No NeoPZ:**
  - uma classe de campo de Nyström ao lado do `RandomField`, sem mexer nele, com a mesma interface (`At`, `SetShift`,
    `LogWeight`, `Direction`);
  - a interpolação de Nyström nos pontos de integração da malha adaptativa.
- **Comparação no Cho2010** (r = −0,5), cada base com o seu próprio direto de referência no fermi, porque a mesma
  semente dá outras realizações:
  - variância retida por M;
  - β(M) do FORM;
  - ganho da AI pelo ponto de projeto, pela banda e pelos primeiros modos.
- **Critério:** as duas bases dão a mesma P_f, dentro do erro dos diretos; a melhor é a de menor custo total para um CoV
  alvo, com o FORM incluído.

### Etapa 2c: o gradiente do FORM pela projeção modal

- **O problema:** o `form_kl.py` calcula o gradiente por diferenças centrais (h = 0,5), com 4M + 1 análises do SRM por
  iteração (121 com M = 30). Os termos fracos (|∂FS/∂u| de 0,001 a 0,003) ficam dominados pelo ruído do FS (~0,001), e a
  direção sai torta: a do M = 10 fica a 17,6° da do M = 20.
- **Por que importa:** o teste linear (etapa 2) mostrou que a AI com o deslocamento fora de direção sai baixa, com −3%,
  −5% e −9% para β = 2, 3 e 4. Uma direção mais certa melhora a AI pelo FORM, justo onde a P_f é pequena.
- **A ideia:** no colapso, o FS do SRM é um ponto limite, e a sua derivada usa só o modo de colapso φ, o incremento do
  último passo do arco. Com fluxo associado (o caso do Cho2010), a tangente é simétrica:
  ∂FS/∂ξₖ = − (φᵀ ∂R/∂ξₖ) / (φᵀ ∂R/∂F).
  - ∂R/∂F já é calculado no caminho de F (`ReductionLoad`).
  - ∂R/∂ξₖ é a variação do resíduo quando o termo k da KL muda c′ e φ′ nos pontos de integração, com os deslocamentos
    congelados: uma montagem por termo, sem resolver nada.
  - É a projeção da sensibilidade nos modos, que o guia de estudo descreve (seção do FORM com a KL).
- **Custo por gradiente:** 1 análise do SRM e 2M montagens do resíduo, contra 4M + 1 análises. Fica viável até com os 2000
  termos de cada campo, e sem o ruído das diferenças.
- **Implementação, sem mexer no material não saturado:**
  1. o `UnsatAnalysis` guarda o modo de colapso φ e ∂R/∂F no fim do `SRMArcLength`;
  2. o Cho2010 ganha a opção `gradiente=`, que monta ∂R/∂ξₖ por diferença no resíduo, com o campo perturbado no termo k;
  3. o `form_kl.py` ganha a opção de pedir o gradiente ao programa, em vez das diferenças centrais.
- **Validação:**
  - no ponto de projeto do M = 10 e do M = 20 (r = −0,5), comparar o gradiente modal com o das diferenças centrais, que
    já está no `form.json`. Os termos fortes (φ′₀, φ′₂, c′₀) têm de bater a menos de 5%;
  - o FORM com o gradiente modal tem de dar o mesmo β, com a mesma tolerância;
  - depois, FORM com M = 30, 100 e 2000, a curva β(M) e a AI com cada ponto.
- **Limitação:** é uma aproximação de análise limite. Ela ignora a dependência do caminho elastoplástico, que a validação
  mede.

### Etapa 3: AI pela entropia cruzada

- **Tamanho do piloto, de graça:** sub-pilotos disjuntos do direto de 4000 (40 de 100, 16 de 250, 8 de 500, 4 de
  1000). Só o deslocamento é recalculado. Depois, uma AI só, com o menor piloto que serve.
- **P_f pequena (r = −0,7, ~0,4%):**
  - estender o direto para ~6000 (~25 falhas, CoV ≈ 20%);
  - FORM M = 20 e a AI pelo ponto de projeto, contra a AI pela banda.
- **Fica para depois:** o Wang 2020 (precisa da média de c' como opção e da entropia cruzada em dois níveis).

### Etapa 3b: ISMC até CoV de 5% no campo de Nyström com 2000 modos

- **Pergunta:** um FORM com gradiente numérico e poucos termos dá uma função de amostragem que leve a ISMC, com o campo
  completo, à P_f do Monte Carlo direto, sem viés e com CoV baixo? Quantas realizações isso custa?
- **Por que pode dar certo:** o peso w = φ(ξ)/h(ξ) é exato para qualquer centro, então o centro de poucos termos não cria
  viés; ele só muda a variância. O risco prático é um centro que erre uma região de falha: os pesos ficam com cauda
  pesada, a P_f converge devagar, e o CoV estimado sai otimista. Num exemplo de viga em Mathematica
  (`~/Dropbox/UnsaturatedSoil/mathematica/viga_estudo_M_ISMC.txt`, seção 23), a ISMC no campo de 500 modos com o centro
  do FORM de 20 deu a P_f do direto (z = 1,2 nas 10 sementes juntas), com o CoV verdadeiro igual ao estimado (5,35 contra
  5,10%).
- **A base:** Nyström na grade de 0,2 × 0,1 m (15 000 nós), com 2000 modos e 99,7% da variância (`nystrom_base.py` com
  os autopares do `nystrom_autopares.wl`).
  - A grade antiga (0,5 × 0,25 m, 2400 nós) erra os autovalores acima de ~50 modos: o λ₁₀₀ em 23% e o λ₁₀₀₀ em 57%,
    contra a grade de 0,2 × 0,1 m.
  - Os FORM leem a cópia com os 100 primeiros modos (os mesmos autovetores, com os mesmos sinais).
  - Os pontos de projeto de Nyström antigos não servem na base nova, porque os sinais dos autovetores mudam.
- **Referências:** os diretos da KL analítica com 2000 termos, `r05_direto` (etapa1, 24 000, P_f = 1,70%, CoV 4,9%) e
  `r07_direto` (etapa3, 110 000, CoV ≈ 5%). É o mesmo campo, em outra base (a analítica retém 98,8% da variância).
- **Os centros (6 FORM, antes da fila, divididos entre os jobs pelas travas; `FORM_LISTA`):**
  - diferenças centrais (h = 0,5) com M = 20 e M = 100 por campo: 81 e 401 análises por iteração;
  - o GIM com o gradiente modal, M = 20: uma análise por iteração;
  - cada um para r = −0,5 e r = −0,7.
- **A ISMC (30 casos):**
  - o SRM (falha: FS < 1, como nos diretos), com os 2000 modos;
  - o deslocamento é o u* nos M primeiros termos de cada campo, e os outros termos ficam sem deslocamento;
  - 5 sementes (21 a 25) por centro e por r, as mesmas nos três centros (números aleatórios comuns).
- **A parada:** `cov_alvo=0.05` no `tarefas.sh`.
  - Antes de lançar um bloco, ele calcula a P_f e o CoV nas realizações já prontas, em ordem. Com CoV ≤ 5%, k ≥ 500 e
    ≥ 100 falhas sorteadas, para.
  - O N do caso é só o teto: 5000 com r = −0,5 e 8000 com r = −0,7.
  - O `Cho2010/ismc_parada.py` dá o N do alvo na ordem das realizações e o N a partir do qual o CoV não volta a passar.
    Ele também compara a P_f no ponto de parada com a P_f com todos os blocos que já rodavam: o viés da parada opcional.
- **A conferência de versão:** as 10 primeiras realizações do `r07_ai_form20` (etapa3) refeitas com o binário novo. O FS
  tem de ser o mesmo.
- **A previsão (sementes 3 a 5 da etapa 3, base analítica, FORM M = 20):**
  - N₅% ≈ 1600 com r = −0,5 e ≈ 2700 com r = −0,7;
  - ganho ≈ 13× e ≈ 38× sobre o direto, contando as 243 análises do FORM.
  - O FORM M = 100 custa ~1200 a 1600 análises a mais. Ele só compensa se o N₅% cair mais que isso.
- **Critérios:**
  - as sementes juntas contra o direto: |z| < 2. Com CoV ≈ 2,2% nas 5 juntas e ≈ 5% no direto, o viés mínimo
    detectável é ≈ 11%;
  - a P_f no ponto de parada e a final concordam;
  - ESS_F e a maior parcela abaixo de 10%;
  - o custo total (FORM + N₅%) por centro e por r.
- **Custo:** ~64 000 realizações (~1350 h de núcleo), cerca de 1 dia nos 4 nós. Os FORM M = 100 levam umas 2 h nos
  primeiros jobs, enquanto os outros já fazem a ISMC dos centros de M = 20.
- **Depois:** a fase 4 pausada (`CAMPANHA=etapa3`, `CASOS=fase4,...`) entra na fila de reserva. Os casos de Nyström
  dela usam a base antiga, de 2400 nós.

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

## Documentos

Cada documento mora num lugar só, conforme a natureza:

- `docs/<tema>/`: a **formulação**, que muda pouco e é a base dos artigos. Os temas são `nao_saturado`, `talude_srm` e
  `confiabilidade_campos`.
- `Projects/PlanoDeTrabalho/`: os **resultados**, gerados das rodadas:
  - `situacao/`: a situação, com os resultados principais, as hipóteses e os próximos passos;
  - `ETAPA0/`, `ETAPAS1a3/`, ...: os resultados completos de cada grupo de etapas (o `ETAPAS1a3` tem todos os do
    Cho2010 e a referência da chuva);
  - `campanhas/`: chuva, rebaixamento e solos;
  - `arquivo/`: as cópias datadas.
- `Projects/<programa>/README.md`: como usar cada programa. Fica junto do código.
- Fora do repositório:
  - `~/Dropbox/UnsaturatedSoil/referencias/`: os artigos de terceiros;
  - `~/projects/runs/<campanha>/`: os dados brutos.

Regras:

1. O relatório de resultados não repete a formulação. Ele cita a seção do documento de formulação.
2. Cada documento é uma pasta com:
   - o gerador (`gera_*.py` ou `atualiza.sh`);
   - o `.tex` e o `.pdf`;
   - a pasta `figs/`.

   O PDF vai para o git; os arquivos auxiliares do LaTeX não vão.
3. O documento vivo não tem data no nome; a capa traz a data e o commit. Antes de uma reunião ou de um envio, guarde uma
   cópia datada em `arquivo/` (`AAAA-MM-DD_nome.pdf`).
4. O índice abaixo sai da mesma lista que o "Mapa dos documentos" do relatório de situação (`documentos.py`). Para
   mudar, edite a lista e rode `python3 -I documentos.py`. Os relatórios de situação e das etapas se
   atualizam juntos com `situacao/atualiza.sh`.

<!-- indice: gerado por documentos.py; nao edite a mao -->

### Formulação: o solo não saturado e o fator de segurança (docs/)

| Documento | Páginas | O que tem |
|---|---:|---|
| `docs/nao_saturado/biot_nao_saturado_agua_ar5.pdf` | 49 | formulação poroelástica não saturada (esqueleto, água e ar), a implementação no NeoPZ, Rahardjo et al. (2007) e benchmarks; a versão 4 está ao lado |
| `docs/nao_saturado/formulacao/formulacao.pdf` | 25 | formulação acoplada u–pw–pa com Mohr–Coulomb associado: discretização, resíduos e jacobiano |
| `docs/nao_saturado/resumo/resumo_nao_saturado.pdf` | 24 | resumo da formulação u–pw–pa |
| `docs/nao_saturado/relatorio_discussao.pdf` | 7 | contraexemplos a Rahardjo et al. (2007): discussão dos resultados, argumentos e verificação |
| `docs/nao_saturado/artigo/artigo.pdf` | 19 | artigo: a intensidade crítica de chuva em taludes não saturados, e quanto as análises desacopladas e com ar passivo a superestimam |
| `docs/talude_srm/slope_mohr_coulomb/relatorio.pdf` | 18 | talude Mohr–Coulomb no NeoPZ: correções da plasticidade, FS e rebaixamento |
| `docs/talude_srm/slope_mohr_coulomb/nota_uqp.pdf` | 6 | nota: formulação de três campos u–q–p |
| `docs/talude_srm/fs_arc_length/relatorio_fs.pdf` | 34 | FS de talude não saturado sob chuva: SRM e GIM com e sem arc-length, sensibilidade e parâmetros |

### Formulação: campos aleatórios, FORM e amostragem por importância (docs/confiabilidade_campos)

| Documento | Páginas | O que tem |
|---|---:|---|
| `docs/confiabilidade_campos/formulacao/formulacao.pdf` | 20 | o modelo, as duas bases de Karhunen–Loève (analítica e Nyström) e a mudança de base, o Monte Carlo, a AI, a variável de controle e o FORM com o gradiente modal (SRM e GIM); no apêndice, um exemplo de uma variável em Mathematica (exemplo_1d.wl) |
| `docs/confiabilidade_campos/implementacao_form_ai/implementacao_form_ai.pdf` | 5 | como o FORM e a AI foram implementados no NeoPZ |
| `docs/confiabilidade_campos/ai_nystrom/ai_nystrom.pdf` | 5 | AI com a KL de Nyström, um roteiro testado em Mathematica |
| `docs/confiabilidade_campos/confiabilidade/confiabilidade.pdf` | 19 | FORM, Monte Carlo e AI em variáveis: a revisão de um código Mathematica e um talude infinito sob chuva (os números saem de Projects/Confiabilidade) |
| `docs/confiabilidade_campos/kl_form_placas/relatorio.pdf` | 26 | estudo do Ghanem & Spanos (Seção 5.3): KL analítica, Nyström e Galerkin, Monte Carlo e FORM nas placas quadrada e curva |
| `docs/confiabilidade_campos/guia_estudo/guia_estudo.pdf` | 22 | guia de estudo que liga tudo (inclui o FORM com a base da KL, Eq. 14) |

### Resultados: o plano de trabalho (Projects/PlanoDeTrabalho)

| Documento | Páginas | O que tem |
|---|---:|---|
| `Projects/PlanoDeTrabalho/README.md` |  | o plano de trabalho (etapas 0 a 5, 2b e 2c), as regras dos documentos e este índice |
| `Projects/PlanoDeTrabalho/situacao/situacao.pdf` | 6 | a situação: o plano, os resultados principais (com as seções do relatório das etapas), as hipóteses e os próximos passos |
| `Projects/PlanoDeTrabalho/ETAPA0/relatorio_etapa0.pdf` | 9 | etapa 0: a AI da chuva no fermi e a prévia do Cho2010 (SRM antigo; superado pelo das etapas 1 a 3) |
| `Projects/PlanoDeTrabalho/ETAPAS1a3/relatorio_etapas.pdf` | 16 | todos os resultados das etapas 0 a 3 e a referência da chuva: SRM corrigido, AI contra o direto, teste linear, as duas bases, o gradiente modal (SRM e GIM), o viés, a variável de controle, a fase 4 parcial, o piloto e o custo |
| `Projects/PlanoDeTrabalho/viga_mathematica/relatorio_viga.pdf` | 9 | a viga com campo aleatório em Mathematica (código em viga_mathematica/codigo): MC, FORM, SORM e ISMC com M, centros baratos (fórmulas de viga), o ponto ótimo do M do FORM e o SORM como função de amostragem |
| `Projects/PlanoDeTrabalho/arquivo/` |  | cópias datadas (AAAA-MM-DD_nome.pdf), nunca editadas: 2026-11-17_relatorio_parcial.pdf |

### Resultados: as campanhas

| Documento | Páginas | O que tem |
|---|---:|---|
| `Projects/PlanoDeTrabalho/campanhas/relatorio_chuva.pdf` | 16 | Monte Carlo da chuva (Rahardjo); a Seção 7 é a AI |
| `Projects/PlanoDeTrabalho/campanhas/relatorio_rebaixamento.pdf` | 11 | rebaixamento do Ceron: Monte Carlo direto e AI (validação: 15,70% contra 15,87%) |
| `Projects/PlanoDeTrabalho/campanhas/relatorio_solos.pdf` | 7 | os três solos do Rahardjo sob chuva |
| `Projects/SlopeRandomFieldRefina/relatorio/relatorio.pdf` | 36 | refinamento por realização e Karhunen–Loève no talude do Ceron (2025); fica com o programa, porque os scripts dele e os do KLCilindroMika gravam nesta pasta |

### Fora do repositório

| Documento | Páginas | O que tem |
|---|---:|---|
| `~/Dropbox/UnsaturatedSoil/referencias/` |  | artigos de terceiros (Cho 2010, Ceron 2025, a tese do Santos, ...) |
| `~/projects/runs/` |  | os dados brutos das rodadas, uma pasta por campanha |

Como usar cada programa: o README.md da pasta dele, em Projects/<programa>/ (fica com o código).

<!-- fim do indice -->
