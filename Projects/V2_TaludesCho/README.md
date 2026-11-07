# V2_TaludesCho — os dois taludes de Cho (2010), secos: SRM, campos de Karhunen-Loève e Monte Carlo

Verifica o fator de segurança por redução de resistência (SRM) do modelo `TPZMatUnsaturatedUPwPa` (Mohr-Coulomb
associado PV2, aqui com poropressões nulas), o gerador de campos aleatórios de Karhunen-Loève (KL) e a cadeia
campo → SRM → Monte Carlo → P_f, em dois exemplos clássicos sem água.

**Referência:** Liao, Zhao, Wu, Miao, Pan e Beer (2025), *Efficient reliability analysis of slopes in spatially variable
soils with active learning-assisted bootstrap polynomial chaos expansion*, Computers and Geotechnics 179:107022. Os
dois taludes são os de Cho (2010), JGGE 136(7):975-984, e de Jiang et al. (2014, 2015).

- **cho1** (Liao, seção 4.1): argila saturada em condição não drenada, talude de 5 m a 1V:2H, c_u médio 23 kPa.
  FS = 1,369 (Liao, SIGMA/W + SLOPE/W: tensões do MEF em superfícies de deslizamento) e 1,366 (Jiang et al. 2014).
  Campo lognormal de c_u (CoV 0,3; h_x = 25 m, h_z = 2,5 m): **P_f = 10,12%** (Monte Carlo simples, 10 000
  realizações do MEF, Tabela 2).
- **cho2** (Liao, seção 4.2): talude c-φ de 10 m a 1V:1H, c = 10 kPa, φ = 30°. FS = 1,205 (Bishop; 1,204 em Cho 2010)
  e 1,284 (elementos finitos, método não informado). Campos lognormais de c (CoV 0,3) e φ (CoV 0,2), h_x = 20 m,
  h_z = 2 m, correlação cruzada r = −0,5: P_f = 0,12% (MEF, 100 000 realizações, Tabela 4); 1,86% por equilíbrio
  limite.

## Geometria (reconstruída das Figs. 2 e 6)

As coordenadas só aparecem nas figuras. Lidas da grade de 0,5 m das figuras e conferidas pela contagem de nós e
elementos do texto: a grade de 0,5 m com estes contornos tem exatamente 981 nós e 910 elementos (cho1) e 1281 nós e
1210 elementos (cho2), os números das seções 4.1 e 4.2.

| caso | contorno (x, z) em m, sentido anti-horário | altura | inclinação |
|---|---|---|---|
| cho1 | (0, 0), (30, 0), (30, 5), (20, 5), (10, 10), (0, 10) | 5 m | 1V:2H (26,57°) |
| cho2 | (0, 0), (30, 0), (30, 5), (20, 5), (10, 15), (0, 15) | 10 m | 1V:1H (45°) |

Apoios (seções 4.1 e 4.2): base (z = 0) fixa nas duas direções, laterais (x = 0 e x = 30 m) fixas na horizontal.

Malha (`MalhaTalude` em `problema.h`): estruturada, de triângulos com u quadrático (T6) e pressões lineares (que aqui
ficam nulas). Linhas verticais a cada h (as quinas são nós); em cada linha, round(z_s/h) intervalos iguais da base à
superfície. Com h = 0,5 m (a grade de Liao): 1810 triângulos no cho1 (os 981 nós de Liao, mais os nós de meio de
lado do T6) e 2400 no cho2. Liao usou quadriláteros e triângulos do GeoStudio (ordem não informada).

## Dados

| grandeza | cho1 | cho2 | fonte / suposição |
|---|---|---|---|
| peso específico γ | 20 kN/m³ | 20 kN/m³ | Tabelas 1 e 3 (γ_sat), constante |
| E, ν | 100 MPa, 0,3 | 100 MPa, 0,3 | Tabelas 1 e 3 (não afetam a carga de colapso, só o caminho do Newton) |
| coesão | c_u = 23 kPa | c = 10 kPa | Tabelas 1 e 3 (médias) |
| atrito φ | 0 | 30° | cho1: **suposto** φ_u = 0 (análise não drenada em tensões totais; o artigo não escreve φ) |
| dilatância ψ | = φ | = φ | **suposto**: o modelo só tem fluxo associado; Liao não informa ψ |
| poropressões | 0 | 0 | seco / tensões totais: p_w = p_a = 0 congeladas (só o equilíbrio é resolvido) |
| distribuição | lognormal | lognormal | eqs. (5)-(7) |
| CoV | c_u: 0,3 | c: 0,3; φ: 0,2 | Tabelas 1 e 3 |
| h_x, h_z | 25 m, 2,5 m | 20 m, 2 m | Tabelas 1 e 3 (sem unidade na tabela; metros implícitos), as distâncias da eq. (1) |
| correlação cruzada | — | r = −0,5 | Tabela 3; Cholesky nos coeficientes da KL, eq. (12) |
| núcleo | ρ2 gaussiano | ρ2 gaussiano | seção 2 ("the latter is chosen"); **confirmado** pelos autovalores das Figs. 3 e 8 (abaixo) |
| truncamento da KL | ε ≥ 95% (M = 9) | ε ≥ 95% (M = 18) | eq. (4), Figs. 3 e 8 |
| domínio da KL | retângulo 30 × 10 m | 30 × 15 m | **inferido** (ver "KL" abaixo); KL_MARGEM = 0,25 dá o 30,5 × 10,5 equivalente ao de Liao |
| realizações | 10 000 | 100 000 | Tabelas 2 e 4 |

Escolhas nossas, que o artigo não fixa: SRM com c/F e tan φ/F (o FS de Liao vem de outro método), estado geostático
com F0 = 0,5, lista de F de 0,1 em 0,1 com bissecção até 0,01 e nova tentativa do F que falhou (`FatorSRM` em
`problema.h`, o `FactorSRMRetry` do `SlopeRandomField`), 25 iterações de Newton por F, semente 20241220. Como cada uma
afeta a comparação está na discussão. Com `FS_METODO=arco`, o FS vem do SRM por arc-length do `LearningMohrCoulomb2`
(branch `UnsaturatedSoil`, commit 2f92365), com os parâmetros padrão dele (`ParamsSRM`); os testes conferem também o
GIM por arc-length.

## O que o código faz

- `problema.h`: os casos (`Cho1`, `Cho2`), a malha, o modelo (`CreateSlopeCompMesh` e `UnsatMaterial` de
  `Comum/RahardjoDriver.h`, os pontos de integração com coordenadas e área), o SRM seco por bissecção (`FatorSRM`) e
  por arc-length (`FatorArco`, que chama `FactorSRMArc` ou `FactorGIM` de `Comum/RahardjoDriver.h`), o gerador de campos
  (`GeradorCampos`: KL no retângulo, lognormal, correlação cruzada), o `Simulador` (monta tudo uma vez e roda
  realizações), a execução do contrato (`RodaExecucao`) e o intervalo de Wilson.
- `KLGauss.h`: KL separável com bases 1D quaisquer. O núcleo **gaussiano** ρ2 = exp(−(dx/h_x)² − (dz/h_z)²) =
  exp(−(dx/h_x)²)·exp(−(dz/h_z)²) é separável; os autopares 1D vêm do método de **Nystrom** (Gauss-Legendre e Jacobi
  cíclico; a autofunção fora dos pontos pela própria equação integral). O **exponencial** ρ1 usa a forma fechada de
  `Comum/KLSeparavel.h`. Os autopares 2D são os produtos dos 1D, escolhidos como em `KL2D`.
- `main.cpp`: o executável do contrato de Monte Carlo (uma realização por execução).
- `testes.cpp`: os testes Catch2.
- `scripts/equilibrio_limite.py`: equilíbrio limite em superfícies circulares (Bishop simplificado; para φ = 0, o
  método sueco, exato para o círculo) na mesma geometria, como referência independente do FS.
- `cluster/`: compilar, rodar a fila de realizações (local ou SLURM) e calcular P_f.
- `dados/`: `autovalores_liao.csv` (Figs. 3 e 8, digitalizados, aproximados), `referencias.csv` (os valores
  publicados) e `mc_cho1_local.csv` (o FS censurado de cada realização do Monte Carlo local do cho1).
- `scripts/plot_campos.py` e `figuras/`: uma realização de cada caso (`SO_CAMPOS=1`), com as escalas das Figs. 4 e 9
  de Liao: `campos_cho1_real1.png` (c_u com os dois núcleos: o gaussiano tem o aspecto estratificado e suave da
  Fig. 4; o exponencial, com a mesma variância e as mesmas distâncias, é bem mais rugoso) e `campos_cho2_real1.png`
  (c e φ com r = −0,5).

## Como compilar e rodar

Na pasta de build já configurada (`-DBUILD_PLASTICITY_MATERIALS=ON -DBUILD_PROJECTS=ON`, com MKL):

```bash
cmake --build <build> --target V2_TaludesCho V2_testes -j 2
cd <build>/Projects/V2_TaludesCho
./V2_testes                      # os testes padrão (~6 min)
V2_NMC=200 ./V2_testes "[.lento]" # P_f do cho1 com os dois núcleos (22 min com 200 realizações cada)
mkdir r1 && cd r1 && CASO=cho1 REAL=1 ../V2_TaludesCho   # uma realização (~1 min com h = 0,5 m)
CASO=cho2 CAMPOS=0 ../V2_TaludesCho                      # o caso determinístico (as médias)
CASO=cho1 CAMPOS=0 FS_METODO=arco ../V2_TaludesCho       # o mesmo, com o SRM por arc-length
```

### Executável: variáveis de ambiente e saídas (o contrato de Monte Carlo)

| variável | padrão | significado |
|---|---|---|
| `CASO` | cho1 | `cho1` ou `cho2` |
| `REAL` | 0 | a realização (as normais da KL vêm de SEMENTE, REAL e do campo) |
| `CAMPOS` | 1 | 0 = caso determinístico com as médias |
| `SEMENTE` | 20241220 | semente base |
| `NUCLEO` | gauss | `gauss` (ρ2, o do artigo) ou `exp` (ρ1, o de `Comum/KLSeparavel.h`) |
| `COV_C`, `COV_PHI` | do caso | CoV de c e de φ (φ só é aleatório no cho2) |
| `LX`, `LY` | do caso | h_x e h_z da eq. (1) (m) |
| `RHO` | do caso | correlação cruzada de ln c e ln φ (cho2: −0,5) |
| `KL_FRAC`, `KL_M` | 0,95; 2000 | fração de variância da KL e máximo de termos |
| `KL_MARGEM` | 0 | o retângulo da KL vira [−m, L + m] × [−m, z_c + m] (0,25: o de Liao) |
| `H_MALHA`, `ORDEM` | 0,5; 2 | espaçamento da malha (m) e ordem de u |
| `FS_TOL`, `FS_F0`, `MAXIT_FROZEN` | 0,01; 0,5; 25 | bissecção, fator do estado geostático e iterações de Newton por F |
| `FS_MAX` | 0 | > 0: busca censurada em FS_MAX (com 1, a indicadora FS < 1 é a mesma; o FS informado é min(FS, 1)); no arc-length, o caminho termina ao passar de FS_MAX |
| `FS_ARCO_TOL` | 0 | tolerância da parada do arc-length (0: a de `ParamsSRM`, 0,001); ver a discussão: 0,01 não serve |
| `FS_METODO` | bissecao | `bissecao` (lista, bissecção e nova tentativa) ou `arco` (SRM por arc-length, `ParamsSRM`: F0 = 0,6, parada com dF < 0,001, 12 iterações por incremento; não usa `FS_TOL`, `FS_F0` nem `MAXIT_FROZEN`; se nem com F0 = 0,4 houver equilíbrio geostático, `situacao=sem_equilibrio_geostatico`, contada como ruptura) |
| `COES`, `PHI_DEG`, `E`, `NU`, `GAMMA` | do caso | o solo |
| `CAMPOS_CSV` | 0 | 1 grava `campos.csv` (c e φ em cada ponto de integração) |
| `SO_CAMPOS` | 0 | 1 só grava `campos.csv` da realização (sem SRM nem `resultado.txt`), para ver os campos |
| `VERBOSO` | 1 | 0 silencia a linha de cada tentativa do SRM |

Saídas na pasta corrente: `fs.csv` (`t_h = 0`, `FS_SRM`, `FS_parada`, `F0`, tentativas, iterações de Newton, s),
`fs_tentativas.csv` (cada solução da busca; vazio com `FS_METODO=arco`), `resumo.csv` (cabeçalho e uma linha: `real`,
`caso`, `nucleo`, `KL_M`, estatística dos campos, `FS_min`, `metodo`, os incrementos, cortes e o motivo da parada do
arc-length, `situacao`, `seg`, ...) e, por último, `resultado.txt` (gravado num `.tmp` renomeado:
`real=… FS_min=… situacao=ok seg=… caso=… nucleo=… metodo=…`). Código de saída 0 em sucesso. Se nem com
c/F0 e tan φ/F0 (F0 = 0,5, 0,25 e 0,1) houver equilíbrio geostático, `situacao=sem_equilibrio_geostatico`, `FS_min=-1`
e código 1: é uma ruptura (FS < 0,1), e `cluster/pf_v2.py` a conta como tal (não ocorreu em nenhuma realização).

### No cluster (SLURM)

O kit em `cluster/` segue o padrão de `Projects/Cluster` (módulos em `Projects/Cluster/modulos.sh`, compilação no nó
de login, um processo de 1 CPU por realização, retomável):

```bash
cd Projects/V2_TaludesCho/cluster
source ../../Cluster/modulos.sh && ./compilar_v2.sh         # build em ./build, executável em ./bin
mkdir -p logs
CASO=cho1 NUCLEO=gauss PRIMEIRA=1 ULTIMA=10000 DET=1 sbatch job_v2.sbatch
CASO=cho1 NUCLEO=exp   PRIMEIRA=1 ULTIMA=2000  sbatch job_v2.sbatch
CASO=cho2 NUCLEO=gauss PRIMEIRA=1 ULTIMA=10000 FS_MAX=1 sbatch job_v2.sbatch   # só P_f (busca censurada)
./pf_v2.py resultados/cho1_gauss resultados/cho1_exp         # P_f, Wilson 95%, estatística e FDA do FS
```

Cada realização fica em `resultados/<CASO>_<NUCLEO>/r<nnnnn>/` (COMB troca o nome da pasta, para conjuntos de
parâmetros diferentes); as concluídas (com `resultado.txt`) são puladas, então o mesmo job pode ser submetido de novo
para completar as que faltam. `./rodar_v2.sh 1 100` roda o mesmo fora do SLURM. O custo medido aqui (1 núcleo
compartilhado): ~50-65 s por realização com h = 0,5 m e a busca completa do FS; ~3-4 s com h = 1 m e FS_MAX = 1.
Para 10 000 realizações do cho1 com h = 0,5 m: ~170 h de CPU (cerca de 7 h num nó de 24 núcleos) com a busca completa
do FS, ou ~40 h (menos de 2 h num nó) com FS_MAX = 1, que basta para a P_f (com h = 0,5 m e FS_MAX = 1 medimos 7 a
28 s, média ~14 s, por realização). O kit foi conferido aqui com `rodar_v2.sh` (sem SLURM); `job_v2.sbatch` não foi
testado num cluster.

## Testes (`testes.cpp`, marca `[V2]`)

| teste | o que compara | tolerância e justificativa | resultado (tempo) |
|---|---|---|---|
| malha | nós da malha com h = 0,5 m contra os da grade de Liao (981 e 1281); área dos triângulos e dos pontos de integração contra a da seção (225 e 300 m²) | contagem exata; áreas a 1e-9 | passou (0,3 s) |
| Nystrom × forma fechada | autopares do núcleo exponencial (Nystrom com 200 pontos) contra `KLSeparavel.h`; convergência no gaussiano | 1e-3 relativo (o bico do núcleo dá erro O(n⁻²): medido 8e-5 a 9e-4); gaussiano: 80 e 160 pontos iguais a 1e-9 | passou (0,3 s) |
| autovalores × Figs. 3 e 8 | M(ε ≥ 95%) = 9 e 18; λ digitalizados | M exato; λ₁ a 3% no retângulo 30 × 10; todos a 4% + 0,5 (leitura) no 30,5 × 10,5; ε a 0,02 | passou (0,1 s) |
| campos: média, CoV, correlação | 20 000 realizações em 6 pontos (defasagens dx = 5, 10 m; dz = 0,5, 1, 2,5 m), núcleos gauss e exp, cho1 e cho2 | 5 erros-padrão contra a expansão truncada; contra o núcleo exato, a cota de Cauchy-Schwarz do resto truncado √((1−v_p)(1−v_q)) mais 5 erros-padrão; r = −0,5 a 5 erros-padrão | passou (0,5 s) |
| cho1 determinístico | FS do SRM (h = 0,5 m) contra 1,369 e 1,366 | 3%: ~1% método (círculo × SRM), ~1% malha, ~0,7% bissecção; contra o nosso LEM 1,354: 2% | passou: 1,344 (63 s) |
| cho2 determinístico | FS do SRM (h = 0,5 m) contra 1,205 (Bishop); abaixo de 1,284 (MEF) | 3% em torno de 1,205; 2% em torno do nosso Bishop 1,215 | passou: 1,219 (51 s) |
| cho1 por arc-length | SRM e GIM por arc-length (h = 0,5 m) contra 1,369 e 1,366, contra a bissecção na mesma malha e entre si | 3% (método e malha); arco − bissecção entre −0,5% e 0,01 + 0,5% (a bissecção fica até 0,01 abaixo do limite discreto, o arco até 0,46% abaixo, o pior caso do estudo de 2f92365); GIM = SRM a 1% (φ = 0: o mesmo limite, e 0,46% + 0,41% das duas paradas) | passou: SRM 1,3414 (129 iterações no caminho, 47 s), GIM 1,3422 (108 iterações, 30 s), bissecção 1,34375 (192 iterações, 12 delas no estado geostático; 61 s) |
| cho2 por arc-length | o mesmo no c-φ, contra 1,205 (Bishop) e abaixo de 1,284; GIM > SRM | 3%; as faixas do cho1; com atrito e FS > 1, F_SRM < λ_GIM (Taylor: o colapso só depende de c/(γH) e φ) | passou: SRM 1,2177 (82 iterações no caminho, 29 s), GIM 1,8728 (99 iterações, 21 s), bissecção 1,21875 (173 iterações, 8 no estado geostático; 47 s) |
| φ = 0 (Tresca) | FS com φ = 0 e φ = 0,001°; FS(2 c_u) = 2 FS(c_u) (h = 1 m) | 0,01 (bissecção); 2,5% na razão | passou: 1,35625, 1,350, 2,7125 (razão 2,000) (44 s) |
| fumaça do Monte Carlo | contrato do executável numa pasta temporária; 2 realizações do cho1 e 1 do cho2 (h = 1 m) | arquivos e chaves; FS entre 0,5 e 3; campos reprodutíveis | passou: FS = 1,331, 0,919, 1,481; cho2 1,419 (49 s) |
| Monte Carlo com `FS_METODO=arco` | o contrato com o arc-length (cho1, REAL = 1, h = 1 m) e o FS da mesma realização contra o da bissecção; com FS_MAX = 1, a mesma indicadora FS < 1 | as faixas do cho1 determinístico | passou: arco 1,3303 contra 1,3313 da bissecção na mesma realização (−0,07%); censurado em 1, o arco para em 1,088 (37 iterações) |
| tolerância do arc-length `[.lento]` | `tol` = 0,001 (padrão) × 0,01 × a bissecção: cho1 e cho2 com h = 0,5 m e as realizações 1 a 30 do cho1 com h = 1 m (`V2_NTOL`); grava `tolerancia_arco.csv` | só a identidade de prefixo: a tolerância só decide onde o caminho para, então FS(0,01) ≤ FS(0,001) e há menos iterações, exatamente; o resto é informado | passou (19 min): ver a discussão |
| P_f do cho1 `[.lento]` | V2_NMC realizações (200) por núcleo, h = 1 m, FS_MAX = 1 | 10,12% dentro do intervalo de Wilson de 95% (gauss); exp só informado | passou: gauss 9,5% (19/200), Wilson [6,2%; 14,4%]; exp 6,0% (12/200), [3,5%; 10,2%] (1319 s) |

Total do padrão: 11 testes, 345 asserções, todos passaram em 6 min 17 s (no `ctest -j 4`, com os outros projetos ao mesmo
tempo).

## Artigo × nosso resultado

| grandeza | artigo | nosso | diferença | fonte do valor do artigo |
|---|---|---|---|---|
| FS cho1 | 1,369 (Liao); 1,366 (Jiang) | 1,344 (h = 0,5 m) | −1,8%; −1,6% | Liao, seção 4.1 |
| FS cho1, malhas h = 1; 0,5; 0,25 m | — | 1,356; 1,344; 1,3375 (extrapolado ≈ 1,33) | −2,8% (h = 0,25) | — |
| FS cho1, equilíbrio limite circular (nosso script) | 1,369 | 1,354 (círculo tangente à base, centro (15, 15), R = 15 m) | −1,1% | — |
| FS cho1, SRM por arc-length | 1,369 | 1,341 (h = 0,5 m; bissecção 1,344) | −2,0% | Liao, seção 4.1 |
| FS cho1, GIM por arc-length | — (φ = 0: igual ao SRM) | 1,342 (h = 0,5 m) | +0,06% contra o nosso SRM por arc-length | — |
| FS cho2 | 1,205 (Bishop); 1,284 (MEF) | 1,219 (h = 0,5 m) | +1,2%; −5,1% | Liao, seção 4.2 |
| FS cho2, SRM por arc-length | 1,205 (Bishop) | 1,218 (h = 0,5 m; bissecção 1,219) | +1,1% | Liao, seção 4.2 |
| FS cho2, GIM por arc-length | — | 1,873 (h = 0,5 m) | 1,54 vezes o SRM (com atrito, GIM > SRM) | — |
| FS cho2, malhas h = 1; 0,5; 0,25 m | — | 1,2375; 1,219; 1,206 | +0,1% (h = 0,25, contra Bishop) | — |
| FS cho2, Bishop (nosso script) | 1,205; 1,204 (Cho) | 1,215 (centro (22,2; 20,0), R = 15,0 m) | +0,8% | — |
| autovalores da KL, cho1 | 100,2; 74,7; 46,3; 23,9; 20,2 (Fig. 3) | 100,0; 74,7; 46,3; 24,0; 20,4 (30,5 × 10,5); 98,1; 71,6; 42,6; 20,9; 19,5 (30 × 10) | < 2% (30,5 × 10,5) | Fig. 3, digitalizado |
| autovalores da KL, cho2 | 78,2; 70,8; 59,9; 47,4; 35,2 (Fig. 8) | 78,2; 70,7; 59,8; 47,4; 35,2 (30,5 × 15,5); 77,3; 69,5; 58,2; 45,5; 33,1 (30 × 15) | < 1% (30,5 × 15,5) | Fig. 8, digitalizado |
| M para ε ≥ 95% | 9 (cho1); 18 (cho2) | 9; 18 (gauss); 133; 295 (exp) | 0 (gauss) | seções 4.1 e 4.2 |
| M do núcleo exponencial (LEM do cho2) | 130 | 295 (30 × 15) | +127% | seção 4.2 |
| P_f cho1 (gauss, M = 9) | 10,12% (MCS 10 000) | 9,38% (75/800, h = 1 m), Wilson 95% [7,54%; 11,59%] | −0,74 p.p. (−7%), dentro do intervalo | Tabela 2 |
| P_f cho1 (exp, M = 133) | — | 6,0% (12/200, h = 1 m), Wilson 95% [3,5%; 10,2%] | — | — |
| P_f cho2 | 0,12% (MEF); 1,86% (LEM) | não rodado aqui (ver pendências) | — | Tabela 4 |

## Discussão

**FS do cho1.** O SRM dá 1,344 com a malha de Liao (h = 0,5 m), 1,8% abaixo de 1,369. O FS cai com o refinamento
(1,356 → 1,344 → 1,3375, convergência de primeira ordem, limite ≈ 1,33), como esperado de um método de deslocamentos
(a malha grossa é rígida demais). O nosso equilíbrio limite na mesma geometria dá 1,354 no melhor círculo, que é
tangente à base rígida (z = 0): para φ = 0 o bloco que gira num círculo é um mecanismo cinematicamente admissível,
então 1,354 é um limite superior do FS de colapso, e o SRM refinado fica abaixo dele (1,3375 com h = 0,25 m),
provavelmente por um mecanismo não circular junto à base rígida (não examinamos o mecanismo). Os valores publicados
(1,369 e 1,366) ficam 1,1% acima do nosso melhor círculo: a busca de superfícies do SLOPE/W com as tensões do SIGMA/W
não é o mesmo problema nem a mesma busca. Os −1,8% se dividem em −1,1% (do valor publicado ao nosso círculo) e −0,7%
(do círculo ao SRM com h = 0,5 m, que ainda está ~1% acima do limite da malha).

**FS do cho2.** O SRM dá 1,219 (h = 0,5 m), 1,2% acima do Bishop de Liao (1,205) e 0,3% acima do nosso Bishop na
mesma geometria (1,215); com h = 0,25 m, 1,206. Para um talude homogêneo com fluxo associado, o SRM e o equilíbrio
limite concordam (Griffiths e Lane 1999), e é o que se vê. O valor de 1,284 que Liao atribui aos "elementos finitos"
**não é reproduzido**: o nosso SRM fica 5,1% abaixo dele. O artigo não diz como esse FS foi obtido (no exemplo 1 foi
SIGMA/W + SLOPE/W, tensões do MEF em superfícies de deslizamento, que num talude íngreme costuma dar FS acima do
Bishop); não é um FS de redução de resistência. A consequência esperada para o Monte Carlo do cho2 (não rodado): com
o FS médio no nível do Bishop, a P_f do nosso SRM deve ficar muito acima de 0,12%, na ordem da do equilíbrio limite
(1,86%, que porém usou o núcleo ρ1).

**SRM e GIM por arc-length.** O `FactorSRMArc` e o `FactorGIM` do commit 2f92365 (`LearningMohrCoulomb2`, branch
`UnsaturatedSoil`; copiados sem mudanças em `Comum/RahardjoDriver.h`) seguem o caminho de equilíbrio com o
arc-length de de Souza Neto et al. (Box 4.4) e dão o maior fator do caminho, em vez do último F em que a Newton
convergiu com F prescrito. Na malha de Liao, o SRM por arc-length fica 0,17% (cho1: 1,3414 contra 1,34375) e 0,09%
(cho2: 1,2177 contra 1,21875) abaixo da bissecção. Como a bissecção dá um F convergido, o limite discreto está acima
dele, e a diferença é o erro da parada do arc-length (dF < 0,001 por incremento), dentro do pior caso de −0,46% do
estudo de 43 casos de 2f92365. Contra Liao, o arc-length dá −2,0% (cho1) e +1,1% (cho2), os mesmos desvios da
bissecção. O ganho é de custo: depois do estado geostático, 28% e 50% menos iterações de Newton (129 contra 180; 82
contra 165) e, no total, 23% e 38% menos tempo (cada iteração do arc-length monta um resíduo a mais, a derivada em F
por diferença finita). A bissecção gasta 79 e 73 dessas iterações nos F que falham (4 e 5 tentativas), e o arc-length
não tem F que falham. O GIM tem aqui uma verificação que
independe do artigo: sem atrito, se σ equilibra o peso λγ com c_u, σ/λ equilibra γ com c_u/λ, e o critério de
Tresca é homogêneo de grau 1 em σ e c_u. Logo o GIM e o SRM têm o mesmo limite, também no problema discreto. No cho1,
o GIM por arc-length dá 1,3422, 0,06% acima do SRM por arc-length (o caminho do GIM parte do estado sem tensões, o do
SRM do geostático com F0 = 0,6: o limite não depende do caminho, como deve ser num material perfeitamente plástico
associado). No cho2 o GIM dá 1,873, 54% acima do SRM: com atrito, aumentar o peso aumenta também a resistência por
atrito, e só a coesão resiste ao aumento (c/(λγH) = N(φ) no GIM e (c/F)/(γH) = N(φ_F) no SRM, com N decrescente
em φ). O GIM por bissecção antigo (agora `FactorGIMBissecao`) não foi rodado aqui: no estudo de 2f92365 ele chegou a
−26,9% do limite.

**Tolerância da parada do arc-length.** O `ParamsSRM` para quando um incremento sobe menos que `tol` em F
(`criterio = 0`). Com 0,001, o padrão do estudo de 2f92365, e com 0,01:

| | bissecção | arco, tol 0,001 | arco, tol 0,01 |
|---|---|---|---|
| cho1, h = 0,5 m | 1,34375 (180 it., 60 s) | 1,3414 (129 it., 48 s) | 1,3365 (108 it., 40 s) |
| cho2, h = 0,5 m | 1,21875 (165 it., 47 s) | 1,2177 (82 it., 29 s) | 1,2158 (72 it., 26 s) |
| cho1, h = 1 m, 30 realizações: FS − bissecção | — | −0,30% em média (−1,00% a +0,36%) | −9,5% em média (−54,6% a −0,2%) |
| idem: tempo total; iterações | 458 s; 5143 | 278 s; 2872 | 176 s; 1809 |
| idem: realizações que trocam de lado em FS = 1 | — | 0 | 3 |

(As iterações da bissecção são contadas sem as do estado geostático, como as do arc-length.) Com 0,01 o caminho para
cedo demais em 6 das 30 realizações, em F ≈ 0,61, logo depois do estado geostático (F0 = 0,6): o FS sai 34 a 55% abaixo,
e três realizações que não rompem passam a romper. O critério confunde um incremento pequeno com o patamar do limite, e
um incremento pode ser pequeno só porque o comprimento l é pequeno: no início do caminho, ou depois de cortes de l
quando a plastificação começa (no talude de Liu, V4, o mesmo erro dá FS = 0,81 em vez de 1,59). O ganho de tempo de
0,01 não compensa. Com 0,001, nenhum caso falhou, mas a margem é estreita: no V4 os incrementos depois dos cortes
tinham dF de 0,006 a 0,008, só 6 a 8 vezes acima de 0,001. O `criterio = 1` do `ParamsArco` não resolve no SRM: ele
para no primeiro incremento, porque no trecho elástico o deslocamento quase não muda com F. Um critério que não
dependa do tamanho do passo, como a rigidez corrente relativa (Bergan): (dF/l)/(dF/l) do primeiro incremento < tol,
evitaria os dois problemas.

**φ = 0.** O Mohr-Coulomb PV2 associado funciona no limite de Tresca: o retorno ao plano principal e às arestas não
usa cot φ (só o retorno ao ápice usa, e o ápice não existe no Tresca), o FS com φ = 0 e com φ = 0,001° coincide a menos
de um passo da bissecção, e o FS é exatamente proporcional a c_u (2,7125/1,35625 = 2,000), como deve ser sem atrito.

**Núcleo da KL.** Os autovalores das Figs. 3 e 8 só são reproduzidos pelo núcleo gaussiano ρ2: o exponencial ρ1 com as
mesmas distâncias dá λ₁ = 81 e 70 (contra 100 e 78) e precisa de 133 e 295 termos para 95%. Com o gaussiano e o
retângulo que contém o talude, M = 9 e 18 (os de Liao) e λ₁ a 2%; a fração ε de Liao implica um traço de ~317 e ~468
(não 300 e 450), e com o retângulo aumentado em meia célula da grade de 0,5 m (30,5 × 10,5 e 30,5 × 15,5, o traço de
uma discretização nodal de 61 × 21 e 61 × 31 nós com 0,25 m² cada) os autovalores coincidem com os digitalizados à
precisão da leitura. Usamos o retângulo 30 × 10 (30 × 15) por padrão; KL_MARGEM = 0,25 reproduz o de Liao (o mesmo M
e quase a mesma variância pontual dentro do talude: esperamos efeito pequeno na P_f, mas não o medimos). Liao diz que as
rodadas de equilíbrio limite do
cho2 usaram ρ1 com M = 130 por campo; com ρ1 no retângulo 30 × 15 obtemos 295 para 95% (o 130 não é reproduzido; ele
coincide com o M do cho1 com ρ1, 133). O truncamento em 95% reduz a variância pontual (v_p ≈ 0,95-0,97) e aumenta um
pouco a correlação a curta distância (até +0,03 contra o núcleo exato, ver o teste dos campos): é o critério do artigo.

**P_f do cho1.** Com o núcleo gaussiano (o do artigo, M = 9), a malha de 1 m e a busca censurada em F = 1, as 200
realizações do teste lento dão P_f = 9,5% (19/200, Wilson 95% [6,2%; 14,4%]); com as 600 seguintes (201 a 800, mesmas
semente e configuração, rodadas com `cluster/rodar_v2.sh`), 56/600 = 9,33%, e no total **P_f = 9,38% (75/800, Wilson
95% [7,54%; 11,59%])**. O intervalo contém os 10,12% de Liao (Monte Carlo simples com 10 000 realizações do MEF) e
também os outros valores da Tabela 2 (10,17%, 10,4%, 10,5%, 11,0%; os 11,56% de Jiang et al. 2014 ficam no limite).
Os resultados de cada realização estão em `dados/mc_cho1_local.csv`. A malha de 1 m no lugar da de 0,5 m não muda
o resultado: nas realizações 1 a 60, rodadas também com h = 0,5 m, romperam as mesmas 5 (2, 9, 22, 31 e 60), com FS
até 0,0125 menor. A concordância tem duas razões: o mesmo núcleo e o mesmo truncamento (ρ2 com 9 termos, os
autovalores da Fig. 3) e um FS determinístico próximo (1,356 com h = 1 m e 1,344 com h = 0,5 m, contra 1,369). Com o
desvio de ln FS da ordem de 0,2 (o que a Fig. 5b sugere: P_f = 10% em FS = 1 com mediana ~1,3), o nosso FS 0,9% (h = 1
m) a 1,8% (h = 0,5 m) menor faria esperar uma P_f 0,5 a 1 ponto percentual maior que a de Liao; a nossa ficou 0,74 ponto
menor, dentro da incerteza da amostra (o intervalo de 95% tem ±2 pontos).

Com o núcleo exponencial (o de `Comum/KLSeparavel.h`, M = 133) e as mesmas sementes, P_f = 6,0% (12/200, Wilson 95%
[3,5%; 10,2%]). Em 200 pares, 10 realizações romperam com os dois núcleos, 9 só com o gaussiano e 2 só com o
exponencial (teste de McNemar exato, p = 0,065): a P_f do exponencial é menor, mas a diferença não é significativa a 5%
com 200 realizações. O sentido é o esperado: para os mesmos h_x e h_z, a média do campo exponencial num comprimento da
ordem da escala de flutuação (~5 m na vertical, a altura do talude) tem variância um pouco menor que a do gaussiano
(Γ² ≈ 0,57 contra 0,64 em 5 m), e a resistência mobilizada na superfície de ruptura é uma média dessas. Só o núcleo
gaussiano é comparável aos 10,12%.

**Custo e robustez do SRM.** Cada F é resolvido por Newton com a tangente consistente e sem busca linear (a análise
não a usa com pressões congeladas). Um passo grande de F pode divergir mesmo havendo equilíbrio: no cho1 (h = 0,5 m),
F = 1,2 a partir de 1,1 não convergiu, mas convergiu depois da bissecção (incrementos menores). O `FactorSRM` de
`Comum/RahardjoDriver.h`, sem a nova tentativa do F que falhou, pararia em 1,194 (`FS_parada` em `fs.csv`) em vez de
1,344, 11% abaixo. Por isso usamos a busca com nova tentativa (`FatorSRM` em `problema.h`). As tentativas que falham
custam 19-25 iterações e dominam o tempo (~50-65 s por realização com h = 0,5 m). O SRM por arc-length
(`FS_METODO=arco`) não tem tentativas que falham com F prescrito: no determinístico com h = 0,5 m, 29 a 47 s.

## O que precisaria mudar fora do projeto (não alterado aqui)

- `Comum/RahardjoDriver.h`, `FactorSRM`: sem a nova tentativa do F que falhou, para numa divergência espúria do
  Newton num passo grande de F (cho1, h = 0,5 m: 1,194 em vez de 1,344, −11%). Incorporar a nova tentativa (como
  `FatorSRM` aqui ou `FactorSRMRetry` do `SlopeRandomField`).
- `Analysis/TPZUnsatAnalysis.cpp`, `Solve`: com pressões congeladas não há busca linear (ela exige dt > 0 e
  `!frozen`) nem parada por divergência ou estagnação; os passos grandes de F divergem mesmo havendo equilíbrio, e cada
  F que falha gasta 19-25 iterações, a maior parte do tempo do SRM.
- `FactorSRMArc`, `FactorGIM` e `SegueCaminho` estão em dois lugares: no `LearningMohrCoulomb2/main.cpp` do branch
  `UnsaturatedSoil` e na cópia de `Comum/RahardjoDriver.h` deste branch. Uma mudança num deles precisa ir para o outro;
  o melhor seria levá-los para a biblioteca (junto de `ArcLengthStepSRM` em `TPZUnsatAnalysis`) ou para um cabeçalho
  comum aos dois branches.
- `Projects/Cluster/mc`: a fila é a do `SlopeRandomField` (executável `bin/SlopeRandomField`, combinações SOLO:chuva,
  `../LearningMohrCoulomb2/cluster/modulos.sh`, que não existe neste branch); para rodar o V2 (e o V4) precisa receber o
  executável, o CASO e trocas VAR=valor livres pelo contrato. Até lá, `V2_TaludesCho/cluster` tem um kit mínimo.
- `Comum/UnsatModel.h` e `Comum/RahardjoDriver.h` definem os mesmos `UnsatParams`, `TMC`, `TUnsat` e `UnsatMaterial` e
  não podem ser incluídos juntos: as malhas estruturadas de `UnsatModel.h` não servem com o `FactorSRM`, e o V2 tem a
  sua `MalhaTalude`.

## Pendências

- P_f do cho1 com 10 000 realizações e h = 0,5 m (a malha de Liao) e a FDA do FS contra a Fig. 5b: rodar no cluster
  (`cluster/`); aqui 800 realizações com o núcleo gaussiano e 200 com o exponencial, h = 1 m (60 também com h = 0,5 m)
  e a busca censurada em F = 1 (sem a FDA acima de 1).
- cho2 estocástico: a correlação cruzada está implementada (eq. 12, nos coeficientes da KL; conferida no teste dos
  campos) e uma realização roda no teste de fumaça, mas a P_f de 0,12% exigiria ~10⁴-10⁵ realizações, fora do
  alcance local; no cluster, ~10⁴ realizações com FS_MAX = 1 dariam a nossa P_f (que deve ficar perto de 1-2%, pois o
  nosso FS determinístico é o do Bishop, não o de 1,284).
- FDA do FS (Figs. 5b e 10b): só curvas no artigo; digitalizar e comparar depois da rodada do cluster.
- O mecanismo de ruptura do SRM (VTK do último estado convergido) não foi examinado; a explicação do FS do cho1 abaixo
  do melhor círculo (mecanismo não circular junto à base) não está verificada.
- O efeito de KL_MARGEM = 0,25 (o retângulo de Liao) na P_f não foi medido.
