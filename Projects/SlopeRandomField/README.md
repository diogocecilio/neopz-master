# SlopeRandomField: talude sob chuva com campos aleatórios e Monte Carlo

Talude de Rahardjo et al. (2007), Série A com α = 45°, sob chuva, com a coesão c', o ângulo de atrito φ' e a
permeabilidade k variando no espaço. Os campos são gerados por Karhunen-Loève com a estatística de Vargas Ceron,
Cecílio, Linn e Maghous (2025, IJNAMG 49:2459). O fator de segurança é calculado só por redução de resistência (SRM),
nos instantes pedidos. Cada execução do programa é uma realização do Monte Carlo.

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
de login, depois de `source ../LearningMohrCoulomb2/cluster/modulos.sh`.

## Rodar

```bash
SOLO=f50 I_MM_H=36 ./rodar_mc.sh 0 9        # 10 realizacoes, uma por nucleo fisico
CAMPOS=0 SOLO=f50 I_MM_H=36 ./rodar_mc.sh 0 0   # caso deterministico (medias), para comparar
./resumo_mc.sh resultados_mc/f50_I36_iii    # P_f e estatistica do FS por instante
```

Os casos de Rahardjo (Série A, 45°):

| `SOLO` | `I_MM_H` |
|---|---|
| f10 | 9, 80, 360 (= 1·k_s) |
| f50 | 9, 80, 36 (= 1·k_s) |
| f100 | 9, 80, 3.6 (= 1·k_s) |

Cada combinação tem o seu próprio Monte Carlo, porque cada uma é um problema diferente. As realizações são retomáveis:
`rodar_mc.sh` pula as que já têm `resultado.txt`.

Variáveis de ambiente:

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
| `CAMPOS_CSV` | 0 | 1 grava os campos de cada ponto em `campos.csv` (`plot_campos.py` desenha) |
| `PARALELO` | núcleos físicos | realizações ao mesmo tempo |
| `MAX_CORTES` | 4 | divisões de um passo que não converge |

No fermi, cada combinação vai num job de um nó:

```bash
SOLO=f50 I_MM_H=36 PRIMEIRA=0 ULTIMA=99 sbatch job_mc.sbatch
```

## Saídas

Em `resultados_mc/<solo>_I<chuva>_<variante>/`:
- `r<nnnnn>/`, uma pasta por realização:
  - `fs.csv`: FS por instante;
  - `resultado.txt`: FS mínimo, instante do mínimo, situação, passos, cortes, iterações de Newton e tempo;
  - `realizacao.txt`: parâmetros, número de termos de KL, variância retida e médias do campo;
  - `run.log`.
- `mc_fs.csv`: o FS de todas as realizações, uma linha por realização e instante.
- `mc_resumo.csv`: o resumo de todas as realizações, uma linha por realização.
- `./resumo_mc.sh` calcula, por instante, a média, o desvio e o CoV do FS, P_f = P(FS < 1) e o CoV de P_f,
  √((1 − P_f)/(N P_f)) (Ceron, eq. 68).
- `plot_mc.py` desenha FS(t) de cada realização, o histograma do FS mínimo e P_f(t).
