# SlopeRandomField: talude sob chuva com campos aleatórios e Monte Carlo

Talude de Rahardjo et al. (2007), Série A com α = 45°, sob chuva, com a coesão c', o ângulo de atrito φ' e a
permeabilidade k variando no espaço. Os campos são gerados por Karhunen-Loève com a estatística de Vargas Ceron,
Cecílio, Linn e Maghous (2025, IJNAMG 49:2459). O fator de segurança é calculado só por redução de resistência (SRM),
nos instantes pedidos. Cada execução do programa é uma realização do Monte Carlo. Tudo o que ela calcula vai para
arquivos CSV, para o pós-processamento. Uma realização interrompida continua de onde parou.

## Modelo

- Malha, solos, condição inicial, contornos e chuva são os do driver de Rahardjo do projeto LearningMohrCoulomb2
  (`../LearningMohrCoulomb2/cluster/fonte/rahardjo.cpp`). O `main.cpp` daquele pacote é incluído para reaproveitar a
  malha, o material e as rotinas de FS.
- Os solos vêm da Tabela 1 de Rahardjo: f10, f50 e f100, com a de Fredlund-Xing = 10/50/100 kPa, m = n = 1,
  θs = 0,45, k_w = k_s Θ⁴ e k_s = 10⁻⁴/10⁻⁵/10⁻⁶ m/s. A resistência é c' = 10 kPa, φ' = φb = 26° e γ = 20 kN/m³.
- `VARIANTE=iii` (padrão) usa as hipóteses do artigo: esqueleto rígido e só a percolação, como no SEEP/W. `VARIANTE=i`
  usa o modelo completo, com u, p_w e p_a acoplados; `VARIANTE=ii`, o mesmo com o ar passivo (p_a = 0).
- **Nas variantes acopladas (i e ii), o FS do SRM não cai abaixo de 1 com a resistência real no transiente.**
  - O transiente resolve a mecânica junto com o fluxo, com c' e φ' reais, a cada passo. Todo estado aceito está em
    equilíbrio com F = 1 para as suas poropressões. Se não estivesse, o passo não convergiria e a realização terminaria
    com `situacao=transiente_nao_convergiu`.
  - O Mohr-Coulomb é associado (ψ = φ'). Pelo teorema do limite inferior, o SRM com essas mesmas pressões dá FS ≥ 1, a
    menos da resolução da busca (0,00625).
  - Quando a chuva levaria o talude abaixo de 1, ele plastifica. A dilatância (ψ = φ' = 26°) entra nos balanços de massa,
    baixa p_w na zona de ruptura e sobe a sucção; com χ = 1, a resistência sobe. O talude fica no limite, em vez de romper.
  - O FS_min fica censurado em 1, e P(FS_min < 1) não mede a ruptura. No Monte Carlo f50, 36 mm/h, variante i: 48 de 91
    realizações com FS_min entre 0,99375 e 1,00625, nenhuma abaixo; o determinístico vai de 1,1125 (12 h) a 1,01875
    (18 h) e 1 (24 h).
  - `F_TRANS` (padrão 1) é o fator de redução de resistência só do transiente (c'/F, tan φ'/F). Com `F_TRANS=0.5`:
    - c' e tan φ' dobram no transiente e o esqueleto não plastifica (`pl_pts` = 0 em `passos.csv`);
    - o acoplamento hidromecânico (elástico) e o ar continuam;
    - o SRM usa a resistência real, então o FS pode ficar abaixo de 1.
  - Determinístico f50, 36 mm/h, variante i (FS em 0, 6, 12, 18, 24, 30 e 36 h):

    | | 0 | 6 | 12 | 18 | 24 | 30 | 36 |
    |---|---|---|---|---|---|---|---|
    | `F_TRANS=1` | 2,119 | 1,475 | 1,1125 | 1,019 | **1** | 1,169 | 1,3125 |
    | `F_TRANS=0.5` | 2,119 | 1,469 | 1,094 | 1 | **0,9625** | 1,156 | 1,306 |
    | iii (rígido) | 2,119 | 1,75 | 1,425 | 1,206 | 1,0625 | 1,131 | 1,1875 |
    | ii (ar livre), `F_TRANS=0.5` | 2,119 | 1,756 | 1,4375 | 1,225 | 1,081 | 1,156 | 1,206 |

    - Quase toda a diferença entre i e iii vem do ar, não da mecânica. A ii, acoplada mas com p_a = 0, fica com a
      iii (+0,01 a +0,02).
    - Na i, a chuva fecha a superfície, a única saída do ar, e ele se comprime: p_a de 30 a 50 kPa. Base e laterais são
      fechadas ao ar.
    - Com χ = 1, isso põe p_w positivo com o solo não saturado: na zona de ruptura em 24 h, p_w = +10,5 kPa com
      S = 0,86 na i, contra +2,3 kPa com S = 1,0 na ii.
    - Esse resultado depende das condições de contorno do ar e do k_ra. O SEEP/W de Rahardjo e Ceron não tem ar.
  - Sensibilidade do ar preso (determinístico, variante i, `F_TRANS=0.5`; p_a na zona de ruptura em 24 h):

    | | 6 h | 12 h | 18 h | 24 h | p_a |
    |---|---|---|---|---|---|
    | base: `KRA_MIN=1e-4`, laterais fechadas | 1,469 | 1,094 | 1 | 0,9625 | 35,5 kPa |
    | `KRA_MIN=0.01` (o piso de Laloui et al. 2003) | 1,55 | 1,175 | 1,056 | 1,006 | 19,9 kPa |
    | `AR_LATERAIS=1` (laterais abertas ao ar) | 1,4625 | 1,094 | 1 | 0,9625 | 35,5 kPa |
    | as duas | 1,55 | 1,181 | 1,0625 | 1,0125 | 19,2 kPa |
    | ii (p_a = 0) | 1,756 | 1,4375 | 1,225 | 1,081 | 0 |

    - Abrir as laterais não muda nada: o ar fica preso entre a frente de molhamento e o nível d'água, longe delas.
    - O piso de k_ra muda o tamanho do efeito: com 1% o ar sai mais fácil, p_a cai à metade e o FS em 24 h sobe
      0,044. O efeito do ar diminui de −0,12 para −0,075, mas continua grande.
    - Na coluna de Liakopoulos (projeto V1 do branch `Verificacao`), o piso de 1e-4 é o que reproduz o benchmark.
    - Com `F_TRANS=1`, o transiente plastifica desde as primeiras horas da chuva: 126 pontos no passo em 5 h, 240 em
      6 h. No fim da chuva são 893 pontos, √J2 da deformação plástica acumulada de 0,16 e 19 cm de deslocamento,
      acelerando (10,8 cm em 18 h, 14,0 cm em 21 h, 18,8 cm em 24 h): o talude escoa no limite.
    - Com `F_TRANS=0.5`, nenhum ponto plastifica (√J2 da plástica ~1e-16), e os 6 cm de deslocamento elástico voltam
      depois da chuva.
    - `F_TRANS=0.3` dá os mesmos FS de 0,5: o transiente já é elástico, e o valor não importa abaixo disso.
  - Numa fila: `COMBS="f50:36:F_TRANS=0.5"`, com resultados na pasta `f50_I36_i_ftrans0.5`.
  - Com dados de `F_TRANS=1`, a ruptura é "chegou ao limite": FS_min ≤ 1,00625, ou seja, FS < 1,01 na grade da busca.
    - Os valores 0,99375, 1 e 1,00625 são uma classe só, a dos taludes no limite. Os 0,99375 não são FS < 1: são
      falhas do Newton do SRM em F = 1, num estado que o transiente já mostrou estar em equilíbrio.
    - Use `posprocessa_mc.py <pasta> --limite 1.01`, ou `LIMITE_FS=1.01` em `resumo_mc.sh`, `plot_mc.py`,
      `sensibilidade_mc.py` e `posprocessa_mc.py`.
    - Uma realização com `situacao` diferente de ok (o transiente parou de convergir) pode ser uma ruptura durante o
      transiente. Ela fica fora dessas estatísticas; `status_mc.sh` conta essas realizações ao lado de P(FS_min < 1,01).
- O FS vem do SRM com nova tentativa do F que falhou (a mesma `FactorSRMRetry` do driver de Rahardjo), nos instantes
  `FS_H` (padrão 0, 12, 24, 36 e 48 h).
- Um passo do transiente que não converge é dividido ao meio até `MAX_CORTES` vezes (padrão 4). Com permeabilidade
  heterogênea, isso evita perder a realização inteira.
- **Regressão:** com `CAMPOS=0`, o programa reproduz o driver de Rahardjo (f50, 36 mm/h): FS = 2,11875 em t = 0 e
  2,025 em t = 1 h na variante i, e o mesmo balanço de água por passo. Com os três CoV iguais a 0, os campos passam pelo
  mesmo caminho por ponto de integração e dão os mesmos FS que `CAMPOS=0`.

## Campos aleatórios

- **Covariância** exponencial separável, C = exp(−|Δx|/Lx − |Δy|/Ly) (Ceron, eq. 64), a mesma para os três campos.
  Os campos são independentes entre si, como no Ceron.
- **Karhunen-Loève:** como o núcleo é o produto de dois núcleos exponenciais 1D, os autopares no retângulo que contém
  a malha são produtos dos autopares 1D, que têm forma fechada (Ghanem e Spanos 1991; a mesma abordagem de Cho 2010,
  com quem o Ceron compara). `KLSeparavel.h` calcula esses autopares, guarda os M maiores produtos e avalia os modos
  diretamente nas coordenadas dos pontos de integração. O Ceron resolve o problema de autovalores por elementos finitos
  no domínio do talude; aqui não é preciso resolver autovalores na malha. Com Lx = 20 m, Ly = 2 m no retângulo
  70 m × 30 m, M = 2000 termos retêm 96,4% da variância (o Ceron relata erro abaixo de 6% com os mesmos 2000).
- **Transformação lognormal** (Ceron, eqs. 65-66): H = exp(μ_ln + σ_ln Z), com σ_ln = √ln(1 + CoV²) e
  μ_ln = ln(μ/√(1 + CoV²)), para c', φ' e k_v.
- **Aplicação aos pontos de integração** (3 por triângulo):
  - c' e φ' vão para `fmatprop = {c, φ, ψ = φ}` da memória do ponto. O critério de Mohr-Coulomb os lê antes de
    aplicar o fator de redução do SRM, então o SRM divide c' e tan φ' de cada ponto. Como χ = 1, φb = φ' em cada ponto.
  - k vai para `fKMult`, um multiplicador por ponto do tensor de condutividade do material,
    K = k_v diag(α, 1) com k_h = α k_v.
- **Estatística padrão** (Ceron, Tabela 2): CoV(c') = 0,3, CoV(φ') = 0,1, CoV(k_v) = 0,6, Lx = 20 m, Ly = 2 m, α = 1.
  As médias são as de Rahardjo: c' = 10 kPa, φ' = 26° e o k_s do solo.

### Mudanças na biblioteca

Para os campos funcionarem, a biblioteca mudou em três pontos. Sem campos, o comportamento é o mesmo de antes.
- `TPZUnsaturatedMem` ganhou `fKMult` (padrão 1).
- `TPZMatUnsaturatedUPwPa` multiplica o tensor de condutividade por esse valor em cada ponto.
- `TPZUnsatAnalysis::ResetMemory`, usado no início de cada busca de FS, passa a preservar `fmatprop` e `fKMult`.
  Antes ele apagava as propriedades de cada ponto.

## Compilar

```bash
./compilar.sh
```

O script compila a biblioteca e o alvo `SlopeRandomField` em `./build` e copia o executável para `./bin`. Opções:
`USAR_MKL` (padrão ON: PARDISO com ordenação METIS), `MKL_THREAD_MODEL` (padrão seq) e `JOBS`. No fermi, compile no nó
de login, depois de `source ../LearningMohrCoulomb2/cluster/modulos.sh`. A troca do executável é atômica: realizações
que já estão rodando continuam com o executável antigo.

## Rodar

```bash
SOLO=f50 I_MM_H=36 ./rodar_mc.sh 0 9            # 10 realizacoes, uma por nucleo fisico
CAMPOS=0 SOLO=f50 I_MM_H=36 ./rodar_mc.sh 0 0   # caso deterministico (medias), para comparar
COMBS="f10:9 f50:80" DET=1 ./rodar_mc.sh 0 99   # varias combinacoes numa fila so, com os deterministicos
./resumo_mc.sh resultados_mc/f50_I36_iii        # P_f, estatistica do FS e correlacoes
python3 plot_mc.py resultados_mc/f50_I36_iii resultados_mc/f50_I36_iii_det
python3 plot_realizacao.py resultados_mc/f50_I36_iii/r00003
```

Com `COMBS`, as realizações das combinações vão intercaladas (a 0 de todas, depois a 1...): se a rodada for cortada,
todas as combinações ficam com o mesmo número de realizações. `DET=1` põe antes o caso determinístico de cada uma.

Os casos de Rahardjo (Série A, 45°):

| `SOLO` | `I_MM_H` |
|---|---|
| f10 | 9, 80, 360 (= 1·k_s) |
| f50 | 9, 80, 36 (= 1·k_s) |
| f100 | 9, 80, 3.6 (= 1·k_s) |

Cada combinação tem o seu próprio Monte Carlo, porque cada uma é um problema diferente.

### Interrupções

- `rodar_mc.sh` pula as realizações que já têm `resultado.txt`. Para continuar um Monte Carlo interrompido (Ctrl-C,
  queda de energia, reinício da máquina), rode o mesmo comando de novo.
- Cada realização em andamento tem um dono: `r<nnnnn>/.ativo/dono` guarda a máquina, o processo e o job do SLURM.
  - Uma realização com dono vivo é pulada, então vários `rodar_mc.sh` ou vários jobs, em nós diferentes, podem dividir a
    mesma fila sem repetir nada.
  - A de um dono morto volta para a fila e continua do checkpoint: um processo que caiu (inclusive um zumbi ainda não
    recolhido), uma máquina reiniciada ou um job que acabou.
  - Um job que o SLURM dá como acabado só conta como morto depois de `INATIVO_MIN` minutos (padrão 30) sem mudança na
    pasta da realização. Um nó que perdeu contato com o controlador (`NODE_FAIL`) pode estar com os processos vivos,
    gravando.
  - Um job que o SLURM recoloca na fila (mesmo número, outra vez) retoma as próprias realizações.
  - Sem resposta do `squeue`, nenhuma trava é tomada, e a fila espera.
  - Cada processo só apaga a própria trava.
  - A trava é um `mkdir`, atômico também no NFS e no GlusterFS do fermi.
- Uma realização que termina com erro (sem `resultado.txt`) volta para a fila. Na `MAX_FALHAS`-ésima vez (padrão 3),
  fica marcada em `falhou.txt` e não é mais tentada.
  - Ctrl-C, o limite de tempo do job e o `scancel` não contam como falha.
  - Um SIGKILL só no programa conta, em geral o OOM killer, por falta de memória.
  - Para tentar de novo: `find <SAIDA> -name falhou.txt -delete; find <SAIDA> -name falhas -delete`.
- **Disjuntor:** muitas falhas seguidas numa máquina, sem nenhum sucesso (2 × `PARALELO`, ou `DISJUNTOR`), indicam um
  problema do ambiente (executável, malha, bibliotecas, disco), não das realizações. A instância para de pegar
  realizações, antes de marcar a fila, e mostra o fim do `run.err`.
- **Parâmetros:** cada combinação guarda em `parametros.txt` os parâmetros do programa com que foi rodada. Outra rodada
  na mesma pasta com outros parâmetros (malha, `FS_H`, CoV, semente...) misturaria resultados e é recusada. Use outra
  `SAIDA`, ou uma entrada de `COMBS` com as trocas (abaixo).
- **Recompilar:** o executável em `bin/` carrega a `libpz.so` de `build/`, então não recompile com realizações rodando.
  `compilar.sh` e `submeter_mc.sh` recusam enquanto houver (`FORCAR=1` passa por cima).
- Com o arquivo `resultados_mc/PARAR`, nenhuma realização nova começa: as que estão rodando terminam, e o
  `rodar_mc.sh` acaba. Para retomar, apague-o e rode o mesmo comando de novo. No fermi, o `submeter_mc.sh` apaga sozinho.
- Uma realização interrompida no meio continua do seu `checkpoint.bin`, sem refazer o que já foi calculado. O
  checkpoint é gravado depois de cada FS e a cada `CHECKPOINT_MIN` minutos (padrão 5). Ele guarda a solução, a memória
  de todos os pontos de integração (tensões, deformações plásticas, c', φ', k), a face de percolação, o tempo, os
  contadores e o tamanho de cada CSV. Na retomada, os CSV voltam ao tamanho que tinham no checkpoint, e o que foi
  gravado depois dele é calculado de novo. Por isso não sobram linhas repetidas.
- A gravação é atômica (um `.tmp` renomeado): uma queda durante a gravação deixa o checkpoint anterior intacto.
- O checkpoint só é aceito pela mesma realização com os mesmos parâmetros. Com outros parâmetros, a realização começa
  do zero. Ao terminar, ele é apagado.
- Teste: uma realização interrompida com `kill -9` no meio do transiente e de novo no meio de um FS, e depois
  retomada duas vezes, dá os mesmos CSV que a rodada sem interrupção.

### Variáveis de ambiente

| Variável | Padrão | O que é |
|---|---|---|
| `SOLO`, `I_MM_H`, `VARIANTE` | f50, 36, iii | solo, chuva (mm/h) e variante |
| `F_TRANS` | 1 | variantes i e ii: redução de resistência só no transiente; 0.5 o deixa elástico e o FS livre abaixo de 1 (ver Modelo) |
| `KRA_MIN` | 1e-4 | variante i: piso da permeabilidade relativa do ar (1e-4, a lei do benchmark ALERT; 0.01, a de Laloui et al. 2003) |
| `AR_LATERAIS` | 0 | variante i: 1 abre as laterais ao ar (p_a = 0), além da superfície |
| `MALHA` | `rahardjo_a45.msh` | ou `rahardjo_a45_fina.msh` (elementos de 0,375 m perto da face) |
| `NIVEL_FINAL` | vazio | reservatório: com ele definido, o nível d'água externo cai de `NIVEL_INICIAL` (padrão `HR`) a `NIVEL_FINAL` (m) em `T_REBAIX_H` horas (padrão 0: de uma vez, logo depois de t = 0); ver "Exemplo do Ceron" |
| `XCRISTA`, `YCRISTA`, `XTOE`, `YTOE` | 30, 30, 40, 20 | reservatório: a crista e o pé, vértices da superfície da malha; a carga da água na face é normal a ela |
| `DIREITA_RES` | 1 | reservatório: 1 põe a lateral direita, abaixo do nível, com p_w hidrostático do reservatório; 0 a deixa com a carga `HR` do início |
| `FS_H` | 0,12,24,36,48 | instantes do FS (h); o transiente vai até o maior |
| `COV_C`, `COV_PHI`, `COV_K` | 0.3, 0.1, 0.6 | coeficientes de variação de c', φ' e k_v |
| `LX`, `LY`, `ALFA` | 20, 2, 1 | distâncias de autocorrelação (m) e k_h/k_v |
| `KL_M`, `KL_FRAC` | 2000, 1 | termos de KL (para antes, ao reter essa fração da variância) |
| `SEMENTE`, `REAL` | 20250416, n | semente e número da realização; a mesma dupla gera o mesmo campo |
| `CAMPOS` | 1 | 0 roda o caso determinístico |
| `CHECKPOINT_MIN` | 5 | minutos entre checkpoints, além do checkpoint depois de cada FS; 0 desliga |
| `CAMPOS_CSV`, `PW_CSV`, `MEC_CSV` | 1, 1, 1 | 0 deixa de gravar `campos.csv`, `pw.csv` ou `mecanismo.csv` (os maiores arquivos) |
| `BANDA_LIMIAR`, `MEC_LIMIAR` | 0.1, 0.01 | a banda de ruptura e os pontos de `mecanismo.csv`: √J2(ε_p) ≥ limiar × máximo |
| `CAMPOS_VTK` | 0 | 1 grava `campos.vtk` (ParaView): c', φ', k_v, log₁₀ k_v e k_h por triângulo (média dos seus 3 pontos) |
| `VTK_H` | vazio | instantes (h) em que o mecanismo de ruptura do SRM vai para `ruptura_pontos.scal_vec.<passo>.vtk` |
| `SO_CAMPOS` | 0 | 1 só gera os campos (`campos.vtk` e `campos.csv`), sem simular; com a mesma `SEMENTE` e o mesmo `REAL`, o campo é o mesmo |
| `COMBS`, `DET` | `SOLO:I_MM_H`, 0 | combinações (`solo:chuva[:VAR=valor,...]`, separadas por espaço) numa fila só; 1 inclui os determinísticos |
| `MASSA_LIMIAR` | 0.5 | a massa que desliza: pontos com \|Δu\| ≥ limiar × máximo (Δu: o último incremento do SRM) |
| `LIMITE_FS` | 1 | pós-processamento: P_f = P(FS < limiar); 1.01 com dados de `F_TRANS=1` (ver Modelo) |
| `PARALELO` | núcleos físicos | realizações ao mesmo tempo |
| `MAX_CORTES` | 4 | divisões de um passo que não converge |
| `MAX_FALHAS` | 3 | tentativas de uma realização que termina com erro, antes de marcá-la em `falhou.txt` |
| `ABANDONO_H` | 6 | dono em outra máquina, fora do SLURM: abandonado depois de tantas horas sem mudança na pasta |
| `INATIVO_MIN` | 30 | dono de um job que o SLURM dá como acabado: morto depois de tantos minutos sem mudança na pasta |

## Rodar no fermi

Os casos de Rahardjo vão todos numa fila só, que os jobs (um nó inteiro cada, 12 realizações ao mesmo tempo) dividem.
Por padrão são as 8 combinações além de f50:36, com 1000 realizações cada, mais os determinísticos.
`submeter_mc.sh` faz tudo e volta em segundos:
- compila no nó de login, em segundo plano, se o executável não existe ou é mais velho que os fontes; os nós de
  computação não têm cmake;
- depois submete os jobs.

```bash
cd ~/neopz-master && git pull origin UnsaturatedSoil && cd Projects/SlopeRandomField
```

```bash
./submeter_mc.sh 4
```

- **Quantos nós:** o argumento é o número de jobs (nós), padrão 4. Os que a cota de CPUs não deixar rodar ficam
  esperando (`QOSMaxCpuPerUserLimit`) e entram quando outro acabar.
- **Mais nós depois:** rodar `submeter_mc.sh` de novo acrescenta jobs à mesma fila. Se um `git pull` tiver trazido
  fontes novos, ele recusa recompilar com jobs na fila (a `libpz.so` deles mudaria). Para acrescentar nós com o
  executável atual: `FORCAR_SEM_COMPILAR=1 ./submeter_mc.sh 2`.
- **Na primeira vez:** ele compila em segundo plano. Um segundo pedido durante a compilação é recusado (só o primeiro
  é submetido ao fim): repita-o depois. Um fonte trocado durante a compilação faz a cadeia compilar de novo antes de
  submeter.
- **Avisos:** o que impediu jobs de rodar fica em `$DATA/monte_carlo/AVISOS`, que o `status_mc.sh` mostra:
  - a cota recusou a submissão ou a continuação;
  - os parâmetros diferem dos de `parametros.txt`;
  - faltam bibliotecas num nó;
  - o disjuntor parou uma instância.
- **Cota de jobs:** se o SLURM recusar um job (`QOSMaxSubmitJobPerUserLimit`), os já aceitos seguem.
- **Limite de tempo:** cada job vai até 4 dias (`TEMPO`, o máximo da partição `short`). Quando faltam 15 min, deixa de
  começar realizações e submete a sua continuação, que espera ele acabar e retoma as interrompidas dos checkpoints.
  - Uma interrompida só fica livre 30 min depois de parar (`INATIVO_MIN`). A continuação roda primeiro o resto da fila
    e a pega na passada seguinte; até lá, o status a mostra como interrompida.
  - Se a cota recusar a continuação, o `status_mc.sh` avisa. Basta rodar `submeter_mc.sh` de novo quando houver vaga.
- **Estimativa:** uma realização leva cerca de 20 min no fermi. Um nó faz cerca de 36 por hora, e as 8000 realizações
  dão uns 9 dias-nó: 2,3 dias com 4 nós rodando juntos e 4,6 dias com 2.
- **Onde ficam os resultados:** em `$DATA/monte_carlo`. O manual do fermi manda rodar os jobs a partir de `$DATA`, a
  partição de dados. Ela é temporária e sem backup: copie os resultados para o PC (ou para o `$HOME`) ao terminar. Sem
  `$DATA`, ficam em `resultados_mc` nesta pasta, e `SAIDA` troca o lugar.
- **Disco:** cerca de 10 GB, 1,2 MB por realização. A cota é de 512 GB em `$DATA` e em `$HOME` (`quota -s`).

Acompanhar (por combinação: concluídas, rodando, esperando, falhas, tempo por realização, FS mínimo médio, P_f até
agora e o fim estimado):

```bash
./status_mc.sh
```

Parar sem perder nada:
- `touch $DATA/monte_carlo/PARAR`: nenhuma realização nova começa, e os jobs terminam as que estão rodando e acabam;
- ou `scancel -n slopeRF`: as interrompidas continuam do checkpoint;
- e, se uma compilação em segundo plano ainda vai submeter jobs, `kill -- -<pid>` com o pid de `logs/compilar_mc.pid`.

Para retomar, ou para trocar de campanha: `./submeter_mc.sh <nós>`. Ele apaga o `PARAR`, também quando precisa
compilar antes; só um `PARAR` criado durante a compilação faz a cadeia não submeter.

Baixar para o PC: o `status_mc.sh` mostra no fim o comando `rsync`, com o caminho de `$DATA` no fermi. Rode-o no PC,
não no fermi. Pode ser no meio da rodada; o `rsync` só traz o que mudou. O destino no PC tem o nome da pasta de
`SAIDA` (`~/projects/runs/sensibilidade/` com `SAIDA=$DATA/sensibilidade`): as pastas das combinações de uma `SAIDA`
própria podem ter o mesmo nome que as da campanha principal (`f50_I36_i`) e não podem cair na mesma pasta. Ele tem esta
forma:

```bash
rsync -av --exclude checkpoint.bin --exclude .ativo --exclude .tomada --exclude '.parar.*' --exclude /PARAR --exclude /execucoes.txt --exclude /campanha.txt diogo@fermi.cesup.ufrgs.br:<$DATA>/monte_carlo/ ~/projects/runs/monte_carlo/
```

A realização n usa os mesmos números aleatórios em todas as combinações: os mesmos campos de c' e φ' e o mesmo
multiplicador k/k_s em f10, f50 e f100, com qualquer chuva. As diferenças de P_f entre combinações vêm então do solo e
da chuva, não de amostras diferentes.

Variáveis (vão para os jobs): `COMBS`, `PRIMEIRA`, `ULTIMA` (padrão 0..999), `DET` (padrão 1), `SAIDA` (padrão
`$DATA/monte_carlo`), `TEMPO`, `PARTICAO`, `RECOMPILAR=1` e as do programa (`VARIANTE`, `MALHA`, `FS_H`, ...). A mão, um job:
`COMBS="f10:9 f50:80" PRIMEIRA=0 ULTIMA=999 DET=1 sbatch job_mc.sbatch`.

## Exemplo do Ceron: rebaixamento do reservatório

O caso de referência da Seção 6 de Vargas Ceron et al. (2025), com o material `TPZMatUnsaturatedUPwPa` e a variante iii
(esqueleto rígido, acoplamento de uma via, como no artigo):
- **Geometria** (`ceron_h5_a45.msh`, de `LearningMohrCoulomb2/cluster/malhas/ceron_h5_a45.geo`): talude de H = 5 m e
  45°, com a crista em (10, 10) e o pé em (15, 5). O domínio vai de x = 0 a 25 m, com a base em y = 0; é o dos exemplos
  da Seção 5.3 do artigo (Figs. 15 e 19), porque a Seção 6 não dá o tamanho do domínio.
- **Solo:** o f50 com φ' = 30° (`PHI_DEG=30`): c' = 10 kPa, γ = 20 kN/m³ e k_s = 1e-5 m/s, com Fredlund-Xing a = 50 kPa,
  n = m = 1, k_rw = Θ⁴, χ = 1. Os campos aleatórios têm a estatística padrão, a da Tabela 2 do artigo: CoV de 30% em c',
  10% em φ' e 60% em k; L_x = 20 m; L_y = 2 m; α = 1.
  - O k_v do Ceron, 1e-6 m⁴/(kN s), dá k = 9,81e-6 m/s; os 2% de diferença só mudam a escala de tempo do transiente.
- **Água:**
  - No início, o nível d'água está na crista (`HL=HR=10`): pressão hidrostática em todo o solo, e a água do reservatório
    sobre a face e o terreno do pé.
  - Em t = 0, o reservatório cai até o pé (`NIVEL_FINAL=5`, h_w = H). Com `T_REBAIX_H=0` a queda é de uma vez; com
    `T_REBAIX_H` > 0, linear nesse número de horas.
  - Abaixo do nível, a superfície fica com p_w = γ_w (nível − y) e com a carga da água, que entra também no SRM. Acima,
    é face de percolação, sem chuva (`I_MM_H=0`).
  - A lateral direita, abaixo do nível, acompanha o reservatório. A esquerda fica com a carga inicial (10 m), como o
    campo distante do artigo (u = 0). A base é impermeável.
- **Diferenças em relação ao artigo:** o Ceron resolve o fluxo permanente saturado, com p = 0 imposto na crista e na
  face. Aqui o fluxo é transiente e não saturado: a crista é face de percolação, e o solo acima da superfície freática
  dessatura.
- **Fator:** o Γ do Ceron é um multiplicador das cargas, Γ = P_mr / (P_γ + P_u). Por isso o caso roda com `FATOR=gim`:
  o aumento da gravidade multiplica o peso, a pressão equivalente positiva (as forças de percolação na zona saturada) e a
  carga do reservatório, e a sucção fica como está, como no `FactorGIM` do projeto. O padrão, `FATOR=srm`, é o fator de
  segurança por redução de resistência das campanhas de chuva.
  - Com o GIM, as colunas `FS_SRM` e `FS_parada` de `fs.csv` e o FS de `resumo.csv` e `resultado.txt` guardam o fator do
    GIM: os nomes ficam, para o pós-processamento. `resumo.csv` ganha a coluna `fator`.
- **Instante:** um só, como no Ceron, `FS_H=0.1` (0,1 h depois da queda, com o passo de tempo padrão de 0,1 h).
  - Logo depois da queda, o fator sobe sem patamar. No caso determinístico (valores médios), com o SRM: 2,325 antes do
    rebaixamento; 1,063 (0,001 h), 1,156 (0,003 h), 1,206 (0,01 h), 1,225 (0,03 h), 1,244 (0,1 h), 1,306 (1 h),
    1,413 (6 h) e 1,469 (24 h). Cada valor até 0,1 h foi calculado com um passo só até o instante; com passos de
    0,01 h, o de 0,1 h dá 1,25.
  - Com o GIM, o determinístico dá 1,40 em 0,01 h e 1,525 em 0,1 h; o Ceron tem Γ = 1,336.
  - Até uns 0,01 h, a pressão no solo saturado ainda está se ajustando ao nível novo, e o volume de água quase não
    muda. Depois, a face drena, e o volume cai dez vezes mais rápido de 0,01 a 0,1 h.
  - Para ter também 0,01 h: `FS_H=0.01,0.1` com `DT_H=0.01`, ao dobro do custo.
- **Resultados:** `passos.csv` ganha no fim `nivel_m` (o nível do reservatório) e `res_m3` (a água que entra pelo
  reservatório no passo, negativa se sai), que entra também no resíduo do balanço.

### Verificação: malha, fator e instante

Casos determinísticos (valores médios), rodados numa máquina sem MKL, nas malhas de 0,375 m (`ceron_h5_a45`), 0,25 m e
0,1875 m (`ceron_h5_a45_fina`) perto da superfície. `FATOR=gim_arco` é o GIM com arc-length de LearningMohrCoulomb2
(commit 2f92365, parâmetros de `ParamsGIM`); `FATOR=gim`, por lista e bissecção.

| caso | 0,375 m | 0,25 m | 0,1875 m | referência |
|---|---|---|---|---|
| seco, c/(γH) = 0,05: GIM com arc-length | 2,030 | 1,945 | 1,907 | Γ = 1,777 (Ceron); 1,81 (Bishop, fator só na coesão) |
| seco, c/(γH) = 0,05: GIM por bissecção | 2,019 | | | |
| seco, c/(γH) = 0,05: SRM | 1,244 | 1,225 | 1,219 | F_s = 1,203 (Ceron), 1,204 (Cho), 1,209 (Bishop) |
| seco, c/(γH) = 0,1: GIM por bissecção | 4,07 | 3,89 | 3,82 | 3,62 (Bishop, fator só na coesão) |
| seco, c/(γH) = 0,1: SRM | 1,644 | 1,625 | 1,619 | 1,60 (Bishop) |
| rebaixamento, 0,01 h: GIM com arc-length | 1,402 | 1,378 | 1,363 | Γ = 1,336 (Ceron) |
| rebaixamento, 0,01 h: GIM por bissecção | 1,40 | 1,381 | 1,369 | |
| rebaixamento, 0,1 h: GIM com arc-length | 1,525 | 1,497 | 1,483 | |

- O exemplo seco da Seção 5.3.2 do Ceron (Γ = 1,777 e F_s = 1,203) corresponde a c/(γH) = 0,05: com c' = 10 kPa,
  γ = 20 kN/m³ e H = 5 m (c/(γH) = 0,1), o Bishop dá 1,60 e 3,62. Os casos secos com c/(γH) = 0,05 rodaram na malha de
  H = 5 m com c' = 5 kPa, o mesmo problema por semelhança.
- Bishop simplificado: superfícies circulares, busca em grade; o "fator só na coesão" (c/F com φ' inteiro) é o
  multiplicador da gravidade.
- O GIM e o SRM caem com a malha. No seco com c/(γH) = 0,05, a extrapolação das três malhas (ordem ajustada de 1 a 1,4)
  dá Γ de 1,79 a 1,83, e o SRM fica perto de 1,21.
- O GIM com arc-length e o por bissecção diferem menos de 0,6% nesses casos.
- O Γ = 1,336 do Ceron fica perto do nosso fator em 0,01 h, extrapolado na malha. O de 0,1 h é 8% maior: a face já
  drenou.

### Malhas focadas: tempo e precisão

As malhas `ceron_h5_a45_foco*` são grossas fora da região das superfícies críticas. Dentro dela, mais uma folga de
0,75 m, usam o tamanho hFoco. A região são as bandas de ruptura (√J2 ≥ 0,1 do máximo) do caso determinístico e das
realizações 0 a 31, calculadas com o GIM com arc-length em 0,01 h na malha `ceron_h5_a45`: x de 8,15 a 16,76 m e y de
4,26 a 9,95 m. As `ceron_h5_a45_h025` e `_h0125` são as uniformes perto da superfície, como a base.

Caso determinístico do exemplo (rebaixamento até o pé, Γ em 0,01 h), rodado com `tempos_malhas.sh` numa máquina sem MKL
(solver skyline), duas rodadas ao mesmo tempo. Os tempos servem para comparar as linhas entre si; no fermi, com o PARDISO,
rode o mesmo script num job (o comando está no cabeçalho dele).

| malha | triângulos | equações | GIM | Γ (0,01 h) | erro vs Ceron (1,336) | tempo do fator (s) | iterações de Newton | soluções | tempo total (s) |
|---|---|---|---|---|---|---|---|---|---|
| focada 0,375 m (`foco375`) | 720 | 3808 | bissecção | 1,4000 | +4,8% | 19,3 | 162 | 16 | 20,5 |
| | | | arc-length | 1,3981 | +4,6% | 10,0 | 69 | 12 | 11,2 |
| uniforme 0,375 m (base) | 1103 | 5858 | bissecção | 1,4000 | +4,8% | 37,7 | 188 | 16 | 39,4 |
| | | | arc-length | 1,4022 | +5,0% | 21,0 | 80 | 13 | 22,8 |
| focada 0,25 m (`foco25`) | 1520 | 7898 | bissecção | 1,3625 | +2,0% | 66,9 | 179 | 15 | 69,9 |
| | | | arc-length | 1,3626 | +2,0% | 46,7 | 92 | 14 | 49,6 |
| uniforme 0,25 m (`h025`) | 2458 | 12798 | bissecção | 1,3813 | +3,4% | 121,8 | 161 | 15 | 128,3 |
| | | | arc-length | 1,3781 | +3,2% | 121,9 | 109 | 16 | 128,6 |
| focada 0,1875 m (`foco1875`) | 2906 | 14960 | bissecção | 1,3562 | +1,5% | 184,3 | 170 | 15 | 193,1 |
| | | | arc-length | 1,3576 | +1,6% | 160,8 | 101 | 15 | 169,6 |
| uniforme 0,1875 m (`fina`) | 4377 | 22558 | bissecção | 1,3687 | +2,5% | 329,7 | 155 | 15 | 344,6 |
| | | | arc-length | 1,3632 | +2,0% | 367,5 | 109 | 14 | 381,4 |
| focada 0,125 m (`foco125`) | 5426 | 27650 | bissecção | 1,3375 | +0,1% | 469,7 | 132 | 15 | 492,0 |
| | | | arc-length | 1,3352 | −0,1% | 525,0 | 93 | 12 | 546,5 |

- **Precisão:** quem manda é o tamanho dos elementos na região das bandas. A focada de 0,375 m dá o mesmo Γ da malha base,
  com 35% menos equações; a de 0,25 m dá o mesmo da uniforme de 0,1875 m (1,3626 e 1,3632), com 35% das equações.
- **Arc-length contra bissecção:** de 30 a 60% menos iterações de Newton. Mais rápido até a focada de 0,1875 m; nas
  maiores (uniforme de 0,1875 m e focada de 0,125 m), 10% mais lento, porque cada iteração do arco custa mais. O valor
  difere menos de 0,6%, e o da bissecção vem na resolução de 0,00625.

### As três fases no fermi (Seções 6.1 a 6.5 do Ceron)

A malha, o instante e o número de realizações ficam a escolher pela verificação acima: o Γ cai com a malha, e o de
0,01 h é o que fica perto do estado do Ceron. Os comandos abaixo usam a malha de 0,375 m e 0,1 h; para a malha fina, troque
`MALHA=ceron_h5_a45_fina.msh` (de 4 a 10 vezes mais cara por realização), e para 0,01 h, `FS_H=0.01 DT_H=0.01`.

Cada fase numa `SAIDA` própria, 1000 realizações por caso. A realização n usa os mesmos números aleatórios em todos os
casos. Ao fim de cada fase, `status_mc.sh` com a mesma `SAIDA` mostra o andamento.

| fase | `SAIDA` | casos | seções do Ceron |
|---|---|---|---|
| 1 | `$DATA/ceron_fase1` | referência (h_w/H = 1) e o determinístico | Fig. 24, Fig. 25 e 6.1 (Spearman) |
| 2 | `$DATA/ceron_fase2` | CoV(k_v) 0; 0,75; 0,9; 1. CoV(c) 0,1; 0,5; 0,7. CoV(φ) 0,05; 0,15; 0,2. s = 1,5; 2; 5; 10; 20; 400. α = 2; 3; 4; 5 (com os determinísticos) | 6.2, 6.3 e 6.4 |
| 3 | `$DATA/ceron_fase3` | h_w/H = 0,5; 0,6; 0,7; 0,8; 0,9 (nível final de 7,5 a 5,5 m, com os determinísticos) | 6.5 |

```bash
cd ~/neopz-master && git pull origin UnsaturatedSoil && cd Projects/SlopeRandomField
```

Compilar antes, no nó de login (1 a 2 min: só o programa mudou), e conferir que `bin/versao.txt` tem o commit novo:

```bash
source ../LearningMohrCoulomb2/cluster/modulos.sh && ./compilar.sh && cat bin/versao.txt
```

As variáveis comuns às três fases:

```bash
export VARIANTE=iii PRIMEIRA=0 ULTIMA=999 DET=1 MALHA=ceron_h5_a45.msh PHI_DEG=30 HL=10 HR=10 XL=0 XR=25 XTOE=15 YTOE=5 XCRISTA=10 YCRISTA=10 NIVEL_FINAL=5 T_REBAIX_H=0 FATOR=gim FS_H=0.1
```

```bash
SAIDA=$DATA/ceron_fase1 COMBS="f50:0" FORCAR_SEM_COMPILAR=1 ./submeter_mc.sh 4
```

```bash
SAIDA=$DATA/ceron_fase2 COMBS="f50:0:COV_K=0 f50:0:COV_K=0.75 f50:0:COV_K=0.9 f50:0:COV_K=1 f50:0:COV_C=0.1 f50:0:COV_C=0.5 f50:0:COV_C=0.7 f50:0:COV_PHI=0.05 f50:0:COV_PHI=0.15 f50:0:COV_PHI=0.2 f50:0:LX=30,LY=3 f50:0:LX=40,LY=4 f50:0:LX=100,LY=10 f50:0:LX=200,LY=20 f50:0:LX=400,LY=40 f50:0:LX=8000,LY=800 f50:0:ALFA=2 f50:0:ALFA=3 f50:0:ALFA=4 f50:0:ALFA=5" FORCAR_SEM_COMPILAR=1 ./submeter_mc.sh 4
```

```bash
SAIDA=$DATA/ceron_fase3 COMBS="f50:0:NIVEL_FINAL=7.5 f50:0:NIVEL_FINAL=7 f50:0:NIVEL_FINAL=6.5 f50:0:NIVEL_FINAL=6 f50:0:NIVEL_FINAL=5.5" FORCAR_SEM_COMPILAR=1 ./submeter_mc.sh 4
```

```bash
SAIDA=$DATA/ceron_fase1 ./status_mc.sh
```

- **Compilação:** com o executável já compilado, as três fases são submetidas com `FORCAR_SEM_COMPILAR=1`. Assim,
  nenhum job começa com um executável velho, que ignoraria o reservatório e o GIM. A coluna `versao` do `resumo.csv` de
  cada realização mostra o commit com que ela rodou.
- **Pastas:** `f50_I0_iii` (referência), `f50_I0_iii_<troca>` (por exemplo `f50_I0_iii_covk0.75`,
  `f50_I0_iii_lx40_ly4`, `f50_I0_iii_nivelfinal7.5`) e os determinísticos em `<pasta>_det`. As trocas de 6.2 e 6.3 não
  mudam o determinístico, que é o de `f50_I0_iii_det`.
- **Custo:** uma realização (um passo de 0,1 h e uma busca do GIM) levou de 50 a 60 s numa máquina sem MKL. Com 48 ao
  mesmo tempo (4 nós), estimam-se 20 min para a fase 1, 7 h para a fase 2 e 2 h para a fase 3. Os jobs que a cota não
  deixar rodar esperam os anteriores.
- **Se a queda instantânea não convergir em muitas realizações** (`situacao` diferente de ok): rode com
  `T_REBAIX_H=0.1` e `FS_H=0.2` numa outra `SAIDA`.

## Análise de sensibilidade

A seção 6 de Vargas Ceron et al. (2025), adaptada à chuva:

| Seção | Estudo | Como |
|---|---|---|
| 6.1 | correlação de Spearman entre o FS e c', φ', k | `sensibilidade_mc.py`, com os CSV de qualquer campanha |
| 6.2 | efeito do CoV de c' (10 a 70%), φ' (5 a 20%) e k_v (0 a 100%) | uma entrada de `COMBS` por caso, `COV_C`, `COV_PHI`, `COV_K` |
| 6.3 | distâncias de autocorrelação (L_x, L_y) = s × (20 m, 2 m), s de 1 a 400 | `LX`, `LY` |
| 6.4 | anisotropia α = k_h/k_v | `ALFA` |

**Correlação de Spearman (6.1):**
- Como em Hu et al. e no Ceron, cada campo vira um valor por realização:
  - c' e φ' médios ao longo da superfície de ruptura, que é a banda de deformação plástica (`banda_c`, `banda_phi_deg`);
  - k médio da massa acima dela, a massa que desliza (`massa_k_mult`).
- A massa são os pontos com o maior incremento de deslocamento entre os dois últimos estados convergidos do SRM.
  Reduzir a resistência não muda a parte elástica, então esse incremento é o do mecanismo.
- `sensibilidade_mc.py <pasta>` dá:
  - a correlação de postos de 1ª ordem e a parcial de 2ª ordem, controlando as outras duas (eqs. 70-71, Fig. 26), com
    o intervalo de 95%;
  - para comparar, o mesmo com as médias no solo inteiro.
- `--t 24` usa o FS em t = 24 h, em vez do FS mínimo.

**Casos de CoV, L e α (6.2 a 6.4):**
- Cada caso é uma entrada de `COMBS` com as suas trocas e vai para a sua pasta. Os casos de uma combinação dividem uma
  fila só:

```bash
SAIDA=$DATA/sensibilidade COMBS="f50:36 f50:36:COV_C=0.1 f50:36:COV_C=0.5 f50:36:COV_C=0.7 f50:36:COV_PHI=0.05 f50:36:COV_PHI=0.15 f50:36:COV_PHI=0.2 f50:36:COV_K=0 f50:36:COV_K=0.75 f50:36:COV_K=0.9 f50:36:COV_K=1" PRIMEIRA=0 ULTIMA=499 DET=0 FS_H=24 ./submeter_mc.sh 4
```

```bash
SAIDA=$DATA/sensibilidade ./status_mc.sh
```

- Use uma `SAIDA` própria, para não se misturar com a campanha principal no `status_mc.sh`.
- O caso de referência (a entrada sem trocas, `f50:36`) vai junto, com o mesmo `FS_H`.
- `sensibilidade_mc.py <pasta> <pasta> ...` dá a tabela por caso (N, média, desvio e CoV do FS, P_f e CoV(P_f), como nas
  tabelas 3 a 5) e a densidade do FS de cada caso (Fig. 27).
- A realização n usa os mesmos números aleatórios em todos os casos, então as diferenças entre casos têm muito menos
  ruído que cada P_f sozinho.
- **Custo:**
  - Uma realização leva cerca de 20 min no fermi, contra segundos na análise limite do Ceron. Os 10 000 ciclos por caso
    dele não cabem aqui.
  - Com `FS_H=24`, o transiente vai só até 24 h e há um FS só: cerca de 8 min. Vale se o FS mínimo da combinação cai
    em 24 h, o fim da chuva; confira `t_min_h` no `mc_resumo.csv` da campanha principal. Foi o caso em todos os testes
    de f50 com 36 mm/h.
  - 500 realizações por caso bastam para as correlações (IC de ±0,09) e para P_f a partir de ~5% (CoV(P_f) de ~20%).
  - Os 20 casos de 6.2 a 6.4 dão então ~4,6 dias-nó: pouco mais de 1 dia com 4 nós.

### Resultados: f50, 36 mm/h, variante i, `F_TRANS=0.5` (fase 2)

17 casos × 300 realizações (0 a 299), FS só em 24 h (`FS_H=24`; que o mínimo cai em 24 h só foi verificado nas
realizações 0 a 142 da base, na campanha principal), 9 min por realização no fermi. A base é idêntica à campanha
principal: nas realizações 0 a 142, os mesmos campos e o mesmo FS em 143 de 143.
Comparação pareada com a base (as mesmas normais em todos os casos; na família s os campos mudam, porque a base da KL é
outra, e o pareamento é mais fraco: correlação de postos do FS com a base de 0,68, 0,61 e 0,42 em s = 2, 5 e 400,
contra 0,82 a 0,99 nos outros casos): McNemar exato com correção de Holm nas 16 comparações, IC de ΔP_f pelo método 10 de
Newcombe. Os números do Ceron são das Tabelas 3 a 5 (10 000 a 90 000 ciclos; o α = 10 não está lá).

| caso | P_f (IC 95%) | ΔP_f (IC 95%) | McNemar, Holm | FS médio | dp do FS | ΔFS médio | P(FS < 0,9) | Ceron: P_f, μ |
|---|---|---|---|---|---|---|---|---|
| base (CoV 0,3 / 0,1 / 0,6; L 20 × 2 m; α = 1) | 63,7% [58,1; 68,9] | — | — | 0,972 | 0,090 | — | 21,3% | 11,5%; 1,353 |
| CoV c' 0,1 | 64,3% [58,8; 69,5] | +0,7 [-3,9; +5,3] pp | 1,00 | 0,978 | 0,063 | +0,005 | 8,7% | 0,45%; 1,367 |
| CoV c' 0,5 | 64,7% [59,1; 69,9] | +1,0 [-2,1; +4,1] pp | 1,00 | 0,959 | 0,122 | -0,014 | 31,3% | 25,78%; 1,324 |
| CoV c' 0,7 | 67,3% [61,8; 72,4] | +3,7 [+0,2; +7,2] pp | 0,55 | 0,939 | 0,148 | -0,034 | 42,3% | 36,34%; 1,296 |
| CoV φ' 0,05 | 64,0% [58,4; 69,2] | +0,3 [-2,9; +3,6] pp | 1,00 | 0,976 | 0,085 | +0,004 | 15,3% | 10,2%; 1,351 |
| CoV φ' 0,15 | 65,7% [60,1; 70,8] | +2,0 [-0,3; +4,3] pp | 0,88 | 0,967 | 0,098 | -0,005 | 25,0% | 12,71%; 1,357 |
| CoV φ' 0,2 | 66,3% [60,8; 71,4] | +2,7 [-0,7; +6,0] pp | 0,88 | 0,960 | 0,108 | -0,012 | 28,0% | 15,18%; 1,361 |
| CoV k_v 0 | 74,7% [69,5; 79,3] | +11,0 [+6,6; +15,4] pp | 2,2e-05 | 0,952 | 0,075 | -0,021 | 24,7% | 11,21%; 1,329 |
| CoV k_v 0,75 | 59,0% [53,4; 64,4] | -4,7 [-7,1; -2,3] pp | 1,2e-03 | 0,984 | 0,096 | +0,011 | 18,7% | 11,41%; 1,375 |
| CoV k_v 0,9 | 52,7% [47,0; 58,2] | -11,0 [-14,5; -7,4] pp | 2,8e-09 | 0,995 | 0,103 | +0,023 | 16,0% | 11,63%; 1,378 |
| CoV k_v 1 | 49,7% [44,0; 55,3] | -14,0 [-17,9; -10,0] pp | < 1e-9 | 1,004 | 0,107 | +0,031 | 14,0% | 11,44%; 1,391 |
| s = 2 | 60,0% [54,4; 65,4] | -3,7 [-9,2; +1,9] pp | 0,97 | 0,975 | 0,117 | +0,002 | 27,0% | 14,99%; 1,361 |
| s = 5 | 58,7% [53,0; 64,1] | -5,0 [-10,6; +0,7] pp | 0,77 | 0,978 | 0,135 | +0,006 | 28,0% | 18,93%; 1,35 |
| s = 400 | 57,3% [51,7; 62,8] | -6,3 [-12,7; +0,1] pp | 0,55 | 0,993 | 0,157 | +0,020 | 27,0% | 21,93%; 1,348 |
| α = 2 | 45,7% [40,1; 51,3] | -18,0 [-22,2; -13,6] pp | < 1e-9 | 1,012 | 0,086 | +0,039 | 8,0% | 3,9%; 1,562 |
| α = 5 | 25,3% [20,7; 30,5] | -38,3 [-43,6; -32,6] pp | < 1e-9 | 1,055 | 0,083 | +0,082 | 2,3% | 0,63%; 1,91 |
| α = 10 | 16,0% [12,3; 20,6] | -47,7 [-53,0; -41,7] pp | < 1e-9 | 1,077 | 0,082 | +0,104 | 1,3% | —; — |

Resultados (análise pareada dos `mc_resumo.csv`, mais rodadas determinísticas para o mecanismo de α):
- **Anisotropia α = k_h/k_v: o efeito mais forte, como no Ceron.** O FS sobe em 300 de 300 realizações, e nenhuma
  realização passa a romper. A causa é a anisotropia da água: o material usa um só tensor para a água e o ar, mas um
  binário de diagnóstico com α só na água reproduz todo o ganho (FS determinístico 0,9625 → 1,08125), e com α só no ar
  ganha 0,0125. Com α = 10, o fluxo horizontal entre as laterais de carga fixa leva para a lateral sob o pé a água que, com
  α = 1, sobe o lençol sob o talude: a lateral sob a crista passa de saída (1,35 m³/m em 24 h) a entrada (10,07 m³/m), e
  a saída pela lateral sob o pé vai de 5,31 a 22,83 m³/m. A carga na base sob o talude sobe de 0,8 a 0,9 m com α = 10,
  contra 3,0 a 4,8 m com α = 1, nas três variantes, e o p_w médio na banda de ruptura cai (o FS das 12 rodadas
  determinísticas segue FS = 1,109 − 0,0184 p̄_w, R² = 0,997). O ar preso descomprime por consequência. O ganho é de
  +0,12 com ar (i), +0,09 em ii e +0,10 em iii. Ele depende do domínio: sem chuva, α = 10 baixa o FS (o estado inicial
  não é permanente, e o fluxo entre as laterais aumenta a água no domínio em 3,2 m³/m, contra 0,3 com α = 1), e as
  laterais de carga fixa (x = 0 e 70 m) devem influir no tamanho do
  efeito (não testado).
- **CoV de k_v: efeito forte aqui, nulo no Ceron.** Com a média fixa, mais variabilidade deixa entrar menos água
  (38,2 → 33,2 m³/m em 24 h), e o FS sobe. A infiltração fica perto da de um talude uniforme com k igual à média geométrica do campo
  (0,2 a 0,7 m³/m abaixo da observada) e longe da média harmônica (2,6 a 5,5 m³/m abaixo). O efeito médio esconde uma
  divisão: 29% das realizações perdem FS com a heterogeneidade (22% por mais de dois passos da bissecção); em média, o
  FS sobe 0,040 onde o k médio da massa que desliza é menor que k_s e cai 0,005 onde é maior. No Ceron o fluxo é
  permanente com pressões prescritas em todo o contorno, e o nível de k não importa; ainda assim, a média de Γ sobe em
  proporção parecida (+4,7%, contra +5,4% do FS aqui), pelo padrão espacial de k. A diferença de P_f vem de onde o
  limiar fica em cada distribuição: com as variações relativas da média e do desvio do Ceron aplicadas à nossa amostra
  de k uniforme, P_f cairia de 74,7% para 61,1; 54,1; 52,9 e 50,8% com CoV de 0,6 a 1 (observados: 63,7; 59,0; 52,7
  e 49,7%); nos níveis intermediários o erro, de até 4,9 pontos, é o mesmo do cálculo com as nossas próprias variações
  (a forma da distribuição também muda).
- **CoV de c' e de φ': P_f quase não muda, mas a distribuição muda muito.** Com CoV de c' 0,7, o desvio do FS passa de
  0,090 a 0,148, a média cai 0,034 e P(FS < 0,9) dobra (21% → 42%), enquanto P(FS < 1) só vai de 63,7% a 67,3%, sem
  significância depois de Holm. A superfície de ruptura procura as zonas fracas (c' médio da banda 7,98 kPa contra 9,93 do
  campo), e a massa encolhe 21%. P(FS < 1) é insensível porque a mediana do FS fica abaixo de 1: a dispersão maior empurra
  P_f para 50%, a média menor o empurra para cima, e os dois se compensam (decomposição de Shapley: a dispersão cancela
  57% da parte da média, de 50 a 64% conforme a versão). No Ceron a referência fica do outro lado (P_f = 11,5%), e com o
  CoV de c' os dois efeitos somam (aproximação lognormal: +5,0 pontos da média e +19,8 da dispersão, contra +10,2 e
  −5,5 aqui). Com o CoV de φ' de 0,1 para 0,2, a média cai aqui (−1,2%) e sobe no Ceron (+0,6%).
- **Distâncias de autocorrelação: o desvio do FS quase dobra (0,090 → 0,157).** P_f cai de 3,7 a 6,3 pontos, sem
  resolver (p ≈ 0,05 a 0,07), e em s = 400 a média sobe 0,020 [0,005; 0,036]. A causa não é identificável: uma superfície de
  resposta ajustada a s = 400 dá +0,021 a +0,036 (convexidade em ln k, menos o ganho de infiltração da heterogeneidade
  interna, mais o fim do caminho fraco medido pela massa ou pela banda). No Ceron P_f sobe com s (11,5 → 21,9%): a
  dispersão cresce com a média
  fixa, e a referência dele está abaixo de 50%. Aqui a dispersão maior empurra P_f para 50%, e em s = 400 a média também
  sobe.
  Numa métrica de cauda comparável à do Ceron (P(FS ≤ 0,875), 12,7% na base), a classe de cada nível (sobe, cai ou
  n.s., McNemar sem ajuste) coincide com a do P_f dele em 11 dos 15 níveis, 13 com Holm do nosso lado. Nos outros, um
  lado é significativo e o outro não: CoV de φ' 0,05 (só no Ceron) e CoV de k_v 0, 0,9 e 1 (só aqui). Com limiares de
  0,85 a 0,95, nenhum nível tem os dois lados significativos em sentidos opostos.
- **Spearman parcial de 2ª ordem (banda e massa):** c' +0,80, φ' +0,50, k_v −0,54, contra +0,88, +0,51 e −0,04 no
  Ceron (Fig. 26, digitalizado). Os sinais e o primeiro lugar de c' na base coincidem (com CoV de c' 0,1, c' passa a ser a menos
  correlacionada). O k pesa de 13 a 17 vezes mais na chuva, porque o
  k absoluto da massa controla a água que entra. O que `sensibilidade_mc.py` chama de "1ª ordem" é a correlação simples,
  e não a parcial de 1ª ordem das eqs. 70-71.
- **A aproximação normal subestima P_f** em 15 de 17 casos (+1,6 pontos em média): o FS tem assimetria positiva
  (0,15 a 0,73). A lognormal erra menos.

Ressalvas:
- **Amostra:** com N = 300 o erro-padrão de P_f é ~2,8 pontos. Pareando, efeitos menores que ~3 a 4 pontos (CoV de c'
  e de φ') ou ~6 (família s) não se resolvem. "Não significativo" não quer dizer sem efeito.
- **Resolução do FS:** a bissecção do SRM (`FactorSRMRetry`, resolução 0,00625, com 1 na lista) não enviesa
  P(FS < 1), mas a média do FS é uma cota inferior (~−0,003). FS = 1 exato aparece em 1 a 4% das realizações, e
  "< 1" contra "≤ 1" muda P_f nessa medida.
- **Escala:** o talude de Rahardjo tem H = 10 m, o do Ceron 5 m. Em relação a H, a nossa base corresponde ao s = 0,5 dele,
  e o nosso s = 2 à referência dele. A KL retém 96,4% da variância na base e 99,999% em s = 400 (o campo não é reescalado).
- **Só em 24 h:** que o mínimo do FS cai em 24 h foi verificado só nas realizações 0 a 142 da base (campanha
  principal); nos outros 16 casos e nas realizações 143 a 299 da base, o FS só foi calculado em 24 h.
- **SRM e a nova tentativa:** nas 11 rodadas determinísticas do f50 com FS em vários instantes, a nova tentativa
  mudou o FS em 18 dos 31 instantes entre 6 e 18 h (até 0,375) e em nenhum de 24 h ou depois. No f100 iii ela muda o FS
  de 24 h de 1,094 para 1,475 (abaixo). Nas 5100 realizações da fase 2 o efeito ainda não foi conferido: `FS_parada`
  está no `mc_fs.csv` de cada caso, no fermi.

### Outros solos com 36 mm/h (determinístico, FS mínimo em 24 h)

| solo | i (com ar), `F_TRANS=0.5` | iii (rígido, sem ar) | i − iii |
|---|---|---|---|
| f10 (k_s = 360 mm/h) | 1,5125 | 1,51875 | ≈ 0 |
| f50 (k_s = 36 mm/h) | 0,9625 | 1,0625 | −0,10 |
| f100 (k_s = 3,6 mm/h) | 1,2625 | 1,475 | −0,21 |

O efeito do ar e do esqueleto (i − iii) parece crescer com I/k_s, mas isso não está estabelecido: os três solos mudam
k_s e a juntos, e no f100 iii o FS de 24 h depende da nova tentativa do SRM (o Newton falha em F = 1,1, 1,2 e 1,4 e
converge na nova tentativa; sem ela, o FS seria 1,094, e i − iii, +0,169). Uma explicação possível, a conferir nos campos
de p_a: no f10 a chuva é menor que k_s, a superfície não fecha e o ar sai; no f100 a superfície encharca logo e prende o
ar. A campanha do f100 com 36 mm/h (os mesmos 17 casos, mais 9, 80 e 3,6 mm/h e
`KRA_MIN=0.01`) está em `$DATA/f100` desde 7/11.

## Saídas

Volumes em m³ por metro de talude (o problema é plano), pressões em kPa, c' em kPa, φ' em graus, k_mult = k/k_s,
tempos em horas e durações (`seg`) em segundos.

### Uma pasta por realização: `resultados_mc/<solo>_I<chuva>_<variante>/r<nnnnn>/`

| Arquivo | Uma linha por | Conteúdo |
|---|---|---|
| `resumo.csv` | realização | tudo o que resume a realização (abaixo) |
| `fs.csv` | instante do FS | FS, busca do SRM, mecanismo e estado hidráulico (abaixo) |
| `fs_tentativas.csv` | solução da busca do SRM | `t_h`, `ordem`, `fase` (`gravidade`: os 4 passos de carga com c'/0,6; `srm`: cada F), `g`, `F`, `convergiu`, `newton`, `nova_tentativa`, `seg` |
| `passos.csv` | passo do transiente | balanço de água e convergência (abaixo) |
| `pw.csv` | nó e instante do FS | `t_h`, `no`, `x`, `y`, `pw_kPa`, `pa_kPa`, `S` (saturação), `superficie` (1 = nó da superfície), `face` (1 = na face de percolação, p_w = 0 imposto) |
| `mecanismo.csv` | ponto plastificado e instante do FS | o último estado convergido do SRM (F = FS): `t_h`, `FS`, `ip`, `x`, `y`, `epsp_raizJ2` (√J2 da deformação plástica), `epsp_xx`, `epsp_yy`, `epsp_xy` (distorção de engenharia), `plastico` (tipo de retorno no último passo; 0 = elástico), `du_m` (\|Δu\| do último incremento do SRM), `c`, `phi_deg`, `k_mult`; só os pontos com √J2 ≥ `MEC_LIMIAR` × máximo |
| `campos.csv` | ponto de integração | `ip`, `x`, `y`, `area` (peso × det J, m²), `c`, `phi_deg`, `k_mult` |
| `triangulos.csv` | triângulo da malha | `n0`, `n1`, `n2`: os nós (de `pw.csv`) de cada triângulo |
| `resultado.txt` | | FS mínimo e situação, em `chave=valor`; marca a realização como concluída |
| `realizacao.txt` | | versão do programa, parâmetros, termos de KL, variância retida e médias dos campos |
| `run.log`, `run.err` | | saída do programa |

Com a malha padrão e 5 instantes de FS, uma realização ocupa cerca de 1,2 MB, ou 1,2 GB a cada 1000 realizações:
`mecanismo.csv` tem 0,5 MB, `campos.csv` 0,3 MB e `pw.csv` 0,2 MB. Para gastar menos disco:
- `MEC_LIMIAR=0.1` grava só a banda de ruptura;
- `CAMPOS_CSV=0`, `PW_CSV=0` ou `MEC_CSV=0` deixam de gravar esses arquivos. Sem os três, uma realização cai para
  0,1 MB.

Os campos de qualquer realização podem ser refeitos depois com `SO_CAMPOS=1`.

**`resumo.csv`**: `real`, `solo`, `I_mm_h`, `variante`, `campos`, `semente`, `KL_M`, `KL_variancia`, `cov_c`, `cov_phi`,
`cov_k`, `Lx`, `Ly`, `alfa`, `area_solo_m2`; média, desvio, mínimo e máximo de cada campo no solo, ponderados pela
área (`c_med`, `c_dp`, `c_min`, `c_max`, e o mesmo para `phi` e `k_mult`); `FS_<t>h` em cada instante; `FS_min` e
`t_min_h`; a banda de ruptura e a massa que desliza no instante do FS mínimo (`banda_area_m2`, `banda_xc`, `banda_yc`,
`banda_ymin`, `banda_c`, `banda_phi_deg`, `banda_k_mult`, `massa_area_m2`, `massa_xc`, `massa_yc`, `massa_c`,
`massa_phi_deg`, `massa_k_mult`); `inf_total_m3`; `situacao`; `passos`, `cortes`, `newton` (transiente) e
`newton_fs`; `retomadas`; `seg` (o tempo de cálculo, somando as execuções retomadas); `memoria_MB` (pico); `versao`
(o commit compilado, de `bin/versao.txt`, com `+mod` se havia alterações locais). Pela `versao` dá para ver se todas
as realizações de uma combinação saíram do mesmo programa.

**`fs.csv`**:
- `t_h`, `FS_SRM` (−1: sem equilíbrio nem com c'/0,6 e tan φ'/0,6), `FS_parada` (onde o SRM sem nova tentativa
  pararia), `tentativas` e `newton_fs` da busca.
- O mecanismo no último estado convergido do SRM:
  - `u_max_m`: maior deslocamento. É só um indicador, porque cresce sem limite perto do colapso.
  - `pl_pts`: pontos em plastificação. `epsp_max`: maior √J2(ε_p).
  - A banda de ruptura, os pontos com √J2 ≥ `BANDA_LIMIAR` × máximo:
    - número de pontos (`banda_pts`), área (`banda_area_m2`) e centroide (`banda_xc`, `banda_yc`);
    - extensão (`banda_xmin`, `banda_xmax`, `banda_ymin`, `banda_ymax`);
    - c', φ' e k médios na banda (`banda_c`, `banda_phi_deg`, `banda_k_mult`).
  - A massa que desliza, os pontos com \|Δu\| ≥ `MASSA_LIMIAR` × máximo. Δu é o incremento de deslocamento entre os dois
    últimos estados convergidos do SRM, que fica concentrado no bloco em rotação acima da banda:
    - `du_max_m`, o maior \|Δu\|;
    - número de pontos (`massa_pts`), área (`massa_area_m2`) e centroide (`massa_xc`, `massa_yc`);
    - c', φ' e k médios na massa (`massa_c`, `massa_phi_deg`, `massa_k_mult`).
- O estado hidráulico no instante: `inf_acum_m3`, `V_agua_m3` (água no solo), `pw_max_kPa` e `face_nos`.
- `seg`.

**`passos.csv`**:
- Identificação do passo:
  - `passo` (de Δt; um passo dividido aparece em várias linhas) e `t_h`, no fim do passo;
  - `dt_h` e `nivel`, o número de divisões (0 = passo inteiro).
- Chuva e balanço de água do passo; a entrada de água é positiva:
  - `chuva_mm_h`; `chuva_m3`, a chuva que chega à superfície;
  - `inf_fluxo_m3`, que entra pelos nós com fluxo imposto; `inf_face_m3`, pela face de percolação (negativo = sai);
  - `esq_m3` e `dir_m3`, pelas laterais com p_w imposto;
  - `dV_m3`, a variação do armazenamento; `residuo_m3` = dV − entradas, que deve ser ≈ 0.
- Estado acumulado: `inf_acum_m3` e `V_agua_m3`.
- Convergência:
  - `face_nos`, nós na face de percolação;
  - `its_face`, iterações da face de percolação; `newton`, iterações de Newton;
  - `normp` e `eu`, os resíduos finais.
- `pw_min_kPa`, `pw_max_kPa` e `seg`.
- Mecânica do transiente (nas variantes i e ii):
  - `u_max_m`, o maior |u| desde o estado inicial;
  - `pl_pts`, os pontos que plastificaram no passo;
  - `epsp_max`, o maior √J2(ε_p) acumulado.
  - Com `F_TRANS=1`, um talude no limite segue plastificando passo a passo.

### A combinação inteira: `resultados_mc/<solo>_I<chuva>_<variante>/`

- `junta_mc.sh <pasta>` junta as realizações concluídas. Os arquivos são refeitos das pastas a cada chamada, então
  podem ser gerados durante o Monte Carlo. `rodar_mc.sh`, `resumo_mc.sh` e `plot_mc.py` o chamam.
  - `mc_resumo.csv`: os `resumo.csv`, uma linha por realização.
  - `mc_fs.csv`: os `fs.csv`, uma linha por realização e instante.
  - Com `TUDO=1`, também `mc_passos.csv` e `mc_tentativas.csv`.
  - A primeira coluna é a realização. Realizações de uma versão anterior do programa, sem `resumo.csv`, entram com as
    colunas que tiverem, e o resto fica vazio.
- `./resumo_mc.sh <pasta>` mostra:
  - por instante: a média, o desvio e o CoV do FS, P_f = P(FS < 1) e o CoV de P_f, √((1 − P_f)/(N P_f)) (Ceron, eq. 68);
  - o mesmo para o FS mínimo de cada realização;
  - a correlação do FS mínimo com as médias dos campos, no solo e na banda de ruptura.
- `plot_mc.py` desenha `mc.png`: FS(t) de cada realização, o histograma do FS mínimo e P_f(t). Desenha também
  `mc_correlacoes.png`: o FS mínimo contra c', φ' e k médios, no solo e na banda.
- `plot_realizacao.py <pasta r<nnnnn>> [t_h]` desenha uma realização:
  - FS(t), a entrada de água por hora e os acumulados;
  - p_w com o lençol, a face de percolação e a banda;
  - c' com a banda e o mecanismo (√J2).
  - O instante padrão é o do FS mínimo.
- `plot_campos.py <campos.csv>` desenha os três campos de uma realização.
- `sensibilidade_mc.py <pasta> [<pasta> ...]`: a análise de sensibilidade (seção acima).
- `posprocessa_mc.py <pasta> [--mapas] [--limite <FS>]`: o pós-processamento completo de uma combinação, em
  `<pasta>/posproc/`. P_f = P(FS < limite), com limite 1, `LIMITE_FS` ou `--limite`. O pós-processamento inclui:
  - convergência do Monte Carlo: P_f com IC de 95%, CoV(P_f), média e desvio contra o número de realizações;
  - estatísticas por instante (`estatisticas.csv`);
  - densidade e distribuição acumulada do FS, também em cada instante;
  - FS(t) com faixas de percentis e o determinístico, e P_f(t);
  - Spearman de 1ª e 2ª ordem, no FS mínimo e em cada instante;
  - dispersão do FS contra c', φ' e k;
  - mecanismo de ruptura, hidráulica e diagnóstico numérico;
  - `relatorio.md`.
  - Com `--mapas`, lê os CSV de cada realização e desenha a frequência de cada ponto na banda de ruptura e o p_w médio
    e o desvio no instante do FS mínimo. Nos resultados de uma versão anterior, sem as colunas da massa, estima o k da
    massa pela geometria da banda.
  - O caso determinístico (`<pasta>_det`) entra sozinho, se existir.
- `compara_mc.py <pasta de referência> <pasta> [...] [--rotulos A,B,...]`: compara Monte Carlos realização a
  realização (a mesma realização tem os mesmos campos em todas as pastas; muda só o modelo), em `comparacao/`:
  - FS_min de cada pasta contra o da referência, com os quadrantes de ruptura;
  - a diferença do FS_min e se ela depende do campo (correlação com c', φ' e k médios);
  - FS(t) e P(FS(t) < 1) médios nas realizações comuns;
  - P(FS_min < 1) e P(FS_min < 1,01) por pasta, e a tabela de rupturas pareada com o teste de McNemar.
  - Exemplo: `compara_mc.py f50_I36_iii f50_I36_i_ftrans0.5 f50_I36_i --rotulos "iii,i F_TRANS=0.5,i F_TRANS=1"`.
- **ParaView:**
  - `campos.vtk` mostra os campos aleatórios de uma realização: colorir por `coesao_kPa`, `atrito_graus` ou `log10_k_v`.
  - `ruptura_pontos.scal_vec.<passo>.vtk` mostra o último estado que convergiu na busca do SRM, com c'/F e tan φ'/F.
    `StrainPlasticJ2` mostra a superfície de ruptura.
  - Para gerar os campos de realizações já rodadas sem simular de novo: `SO_CAMPOS=1 SOLO=f50 I_MM_H=36 ./rodar_mc.sh 0 6`.
