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
  usa o modelo completo, com u, p_w e p_a acoplados.
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
./resumo_mc.sh resultados_mc/f50_I36_iii        # P_f, estatistica do FS e correlacoes
python3 plot_mc.py resultados_mc/f50_I36_iii resultados_mc/f50_I36_iii_det
python3 plot_realizacao.py resultados_mc/f50_I36_iii/r00003
```

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
| `MALHA` | `rahardjo_a45.msh` | ou `rahardjo_a45_fina.msh` (elementos de 0,375 m perto da face) |
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
| `PARALELO` | núcleos físicos | realizações ao mesmo tempo |
| `MAX_CORTES` | 4 | divisões de um passo que não converge |

No fermi, cada combinação vai num job de um nó:

```bash
SOLO=f50 I_MM_H=36 PRIMEIRA=0 ULTIMA=99 sbatch job_mc.sbatch
```

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
| `mecanismo.csv` | ponto plastificado e instante do FS | o último estado convergido do SRM (F = FS): `t_h`, `FS`, `ip`, `x`, `y`, `epsp_raizJ2` (√J2 da deformação plástica), `epsp_xx`, `epsp_yy`, `epsp_xy` (distorção de engenharia), `plastico` (tipo de retorno no último passo; 0 = elástico), `c`, `phi_deg`, `k_mult`; só os pontos com √J2 ≥ `MEC_LIMIAR` × máximo |
| `campos.csv` | ponto de integração | `ip`, `x`, `y`, `area` (peso × det J, m²), `c`, `phi_deg`, `k_mult` |
| `triangulos.csv` | triângulo da malha | `n0`, `n1`, `n2`: os nós (de `pw.csv`) de cada triângulo |
| `resultado.txt` | | FS mínimo e situação, em `chave=valor`; marca a realização como concluída |
| `realizacao.txt` | | parâmetros, termos de KL, variância retida e médias dos campos |
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
`t_min_h`; a banda de ruptura no instante do FS mínimo (`banda_area_m2`, `banda_xc`, `banda_yc`, `banda_ymin`,
`banda_c`, `banda_phi_deg`, `banda_k_mult`); `inf_total_m3`; `situacao`; `passos`, `cortes`, `newton` (transiente) e
`newton_fs`; `retomadas`; `seg` (o tempo de cálculo, somando as execuções retomadas); `memoria_MB` (pico).

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
- **ParaView:**
  - `campos.vtk` mostra os campos aleatórios de uma realização: colorir por `coesao_kPa`, `atrito_graus` ou `log10_k_v`.
  - `ruptura_pontos.scal_vec.<passo>.vtk` mostra o último estado que convergiu na busca do SRM, com c'/F e tan φ'/F.
    `StrainPlasticJ2` mostra a superfície de ruptura.
  - Para gerar os campos de realizações já rodadas sem simular de novo: `SO_CAMPOS=1 SOLO=f50 I_MM_H=36 ./rodar_mc.sh 0 6`.
