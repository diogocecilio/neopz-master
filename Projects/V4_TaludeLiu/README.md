# V4_TaludeLiu — talude não saturado de 5 m (1:2) sob chuva permanente: poropressões, FS e Monte Carlo

Verifica a cadeia **percolação não saturada → poropressões congeladas → FS por redução de resistência (SRM)** do modelo
`TPZMatUnsaturatedUPwPa` (variante iii: esqueleto rígido no transiente, ar à pressão atmosférica, só a equação de
Richards), o parâmetro de Bishop constante (χ = tan φ_b / tan φ'), e o Monte Carlo com campos lognormais de c', φ' e
k_s gerados por Karhunen-Loève (KL), num exemplo publicado.

**Referência:** Lei-Lei Liu, Han Han, Shi-Ya Huang, Lei Huang e Yung-Ming Cheng (2024), *Combined influence of rotated
transverse anisotropy in mechanical parameters and hydraulic conductivity of soils on reliability evaluation of an
unsaturated slope under rainfalls*, Computers and Geotechnics 165:105938. GeoStudio 2012: SEEP/W (percolação
permanente) e SLOPE/W (Bishop simplificado, busca "Entry and Exit"), 2 000 realizações de Monte Carlo com campos de
Cholesky. O exemplo é o mesmo de Jiang et al. (2014).

O que está aqui:

- **L1 (determinístico):** o regime permanente com q/k_s = 0,01 (os perfis de p_w das seções X-X e Y-Y da Fig. 3) e o
  FS (1,525 por Bishop no artigo; o nosso por SRM).
- **L4 (estocástico):** os cenários 1 (só k_s aleatório), 2 (c' e φ') e 3 (os três) da Tabela 3 com acamamento
  horizontal (α = 0°) e vertical (α = 90°, Lx e Ly trocados): média e desvio do FS (Fig. 11), a área da massa que
  desliza (Figs. 7 e 8) e o modo de ruptura (Figs. 9 e 10). O executável segue o contrato do Monte Carlo (uma
  realização por execução) e a fila do cluster está em `cluster/`.

## Geometria e contornos (Fig. 1 e seção 3.1)

Pontos (m): C(0, 0), D(30, 0), E(30, 1), F(30, 5), G(20, 5) (pé), H(10, 10) (crista), A(0, 10), B(0, 6). Altura 5 m,
inclinação 1:2, comprimento 30 m. Seção X-X: a vertical pela crista (x = 10 m); Y-Y: a vertical pelo ponto I da meia
altura da face (x = 15 m, terreno em z = 7,5 m).

Malha (`malhas/talude_liu.geo`, Gmsh): triângulos de tamanho h, u quadrático (T6) e p_w, p_a lineares.
`talude_liu_h05.msh` (h = 0,5 m, o tamanho dos quadriláteros do artigo: 1168 nós de vértice, 2181 triângulos, 11 368
equações) é a malha do executável; `talude_liu_h10.msh` (h = 1 m, 325 nós, 571 triângulos, 3 090 equações) é a dos testes padrão.

Percolação (variante iii): chuva q = 0,01 k_s = 1×10⁻⁷ m/s em toda a superfície AHIGF, normal à superfície e por
comprimento real do contorno (o "unit flux" do SEEP/W; a Fig. 2 desenha as setas normais à face; `CHUVA_VERTICAL=1`
usa a projeção horizontal, 10 % a menos de água na face); carga total de 6 m em BC (p_w = γ_w (6 − z)) e de 1 m em ED;
AB, FE e a base impermeáveis; p_a = 0 em todos os nós; um nó da superfície com p_w > 0 vira face de percolação
(p_w = 0 enquanto a entrada exigida for menor que a chuva; aqui nenhum nó fica saturado no permanente). Mecânica (SRM):
base fixa, laterais fixas na horizontal.

## Dados

| grandeza | valor | fonte / suposição |
|---|---|---|
| c', φ' | 10 kPa, 15° | Tabela 1 (médias) |
| φ_b | 10° = 2/3 φ' | Tabela 1 e seção 3.1 (duas variáveis de estado independentes, Fredlund e Rahardjo 1993) |
| χ de Bishop | tan 10°/tan 15° = 0,658, constante | nosso: c' + (σ − u_a) tan φ' + (u_a − u_w) tan φ_b = c' + (σ − u_a + χ s) tan φ'; χ = 1 onde p_w > 0 (a biblioteca já faz) |
| γ | 20 kN/m³, constante | seção 3.1 (independente da saturação) |
| E, ν | 100 MPa, 0,3 | **não dados** (equilíbrio limite); os de Griffiths e Fenton (2004) para esta família de taludes; o FS do SRM praticamente não depende deles |
| dilatância | ψ = φ' | **suposto** (o modelo só tem fluxo associado) |
| retenção | Fredlund-Xing, a = 5 kPa, n = 2, m = 1, θ_s = 0,4 (porosidade), θ_r = 0 | eqs. 7-8 e seção 3.1; a biblioteca usa C(ψ) = 1 (o artigo usa ψ_r = 3000 kPa: C ≥ 0,9989 até 20 kPa, diferença em S < 0,002, testado) |
| permeabilidade | k_r = S^p, **p = 5,57** | o artigo usa a FXH (Fredlund, Xing e Huang 1994, eq. 9), que a biblioteca **não tem**: p ajustado à FXH por mínimos quadrados de log₁₀ k_r entre ψ_aev e 20 kPa (`scripts/ajuste_kr.py`); resíduo máximo 9 % |
| ψ_aev da FXH | 2,58 kPa | **não dado**: a construção gráfica usual (tangente no ponto de inflexão da curva θ × log ψ, interceptada com θ = θ_s) |
| k_s | 1×10⁻⁵ m/s, isotrópico | Tabela 1 |
| K_w | 2,2×10⁶ kPa | nosso (só o armazenamento do transiente; não muda o permanente) |
| q | 1×10⁻⁷ m/s (q/k_s = 0,01) | seção 3.1 |
| campos | lognormais; COV 0,3 (c'), 0,2 (φ'), 0,6 (k_s) | Tabela 1 |
| correlação | gaussiana, θ₁ = 30 m, θ₂ = 3 m (eq. 1) → **exponencial separável** exp(−\|Δx\|/L_x − \|Δy\|/L_y) com L = √π θ/2: 26,59 m e 2,659 m | o nosso gerador (`Comum/KLSeparavel.h`) só tem o núcleo exponencial: mesma escala de flutuação (δ = √π θ no gaussiano, 2L no exponencial) |
| α = 90° | L_x = 2,659 m, L_y = 26,59 m | eq. 1 com α = 90°: o eixo maior na vertical |
| correlação cruzada c'-φ' | **0 (independentes)** por padrão; `RHO=-0.5` dá a do artigo | Tabela 1: ρ = −0,5 (eqs. 3-5); aqui Z_φ = ρ Z_c + √(1 − ρ²) Z₂, a mesma correlação das normais que a Cholesky do artigo |
| φ_b aleatório | φ_b = 2/3 φ' em cada ponto, via χ constante | o artigo gera φ_b como campo com a estatística de φ'; com χ constante, tan φ_b = 0,658 tan φ' e tan(2φ'/3)/tan φ' só varia de 0,644 a 0,664 entre 9° e 24° |
| KL | 2 000 termos (98-99 % da variância no retângulo 30 × 10 m) | nosso (`KL_M`, `KL_FRAC`) |
| realizações | 2 000 no artigo | o executável roda uma por vez; a fila em `cluster/` |

Escolhas nossas que o artigo não fixa: o regime permanente por um transiente longo com chuva constante (dt de 1 h
dobrando até 1 ano; para quando a maior variação de p_w num passo com dt ≥ 30 dias fica abaixo de 10⁻³ kPa; condição
inicial hidrostática com o lençol reto de 6 m a 1 m e a sucção limitada ao patamar k_r(ψ) = q/k_s do modelo; não muda
o permanente, só o caminho); SRM com c'/F e tan φ'/F (a sucção também: χ s tan φ'/F = s tan φ_b/F), estado geostático
com F₀ = 0,6, F de 0,1 em 0,1, bissecção até 0,01 e nova tentativa do F que falhou (`FactorSRMRetry`, o do
`SlopeRandomField`), 25 iterações de Newton por F; a massa que desliza é o conjunto dos pontos com |du| ≥ 0,5 max|du|
(du: o incremento de deslocamento entre os dois últimos estados convergidos do SRM; o padrão do `SlopeRandomField`, o
meio da banda de cisalhamento); D_f = 10 m − a cota mais baixa da massa; semente 20231122.

## Como compilar e rodar

Aqui (pasta de build já configurada, `$B`):

```bash
cmake --build $B --target V4_TaludeLiu V4_testes -j 2
mkdir -p det && cd det && CASO=det $B/Projects/V4_TaludeLiu/V4_TaludeLiu          # ~1,5 min (malha de 0,5 m)
CASO=s3 ANGULO=0 REAL=7 $B/Projects/V4_TaludeLiu/V4_TaludeLiu                      # uma realização do cenário 3
CASO=s2 REAL=7 PW_PERMANENTE=../det/pw.csv $B/Projects/V4_TaludeLiu/V4_TaludeLiu  # k_s constante: reusa o permanente
```

Variáveis do contrato: `CASO=det|s1|s2|s3`, `ANGULO=0|90`, `REAL=<n>`, `CAMPOS=0|1` (0: as médias), `SEMENTE`; trocas:
`COV_C`, `COV_PHI`, `COV_K`, `RHO`, `THETA1`, `THETA2`, `LX`, `LY`, `KL_M`, `KL_FRAC`, `MALHA` (nome em `malhas/` ou
caminho), `P_KR`, `Q_KS`, `KS`, `COES`, `PHI_DEG`, `PHIB_DEG`, `E`, `NU`, `FS_TOL`, `FS_F0`, `MAXIT_FROZEN`,
`MASSA_LIMIAR`, `TOL_PERM`, `DT0_H`, `DT_MAX_H`, `CHUVA_VERTICAL`, `PW_PERMANENTE`; saídas opcionais `PW_CSV` (padrão
1), `CAMPOS_CSV`, `MEC_CSV`, `VTK` (o mecanismo em `ruptura_pontos.scal_vec.0.vtk`), `VERBOSO` (0, 1, 2).

Saídas na pasta corrente: `passos.csv` (cada passo do transiente com o balanço de água), `pw.csv` (p_w, p_a e S em cada
nó no permanente), `secoes.csv` (p_w nas seções X-X e Y-Y a cada 0,5 m), `fs_tentativas.csv` (cada solução do SRM),
`fs.csv` (uma linha: `t_h`, `FS_SRM`, a massa que desliza, D_f, o modo e as áreas com os limiares 0,1, 0,2, 0,3 e
0,5), `resumo.csv` (cabeçalho e uma linha: `real`, `FS_min`, as estatísticas dos campos, o balanço, ...) e, por último,
`resultado.txt` (`real= caso= angulo= FS_min= massa_area_m2= Df_m= modo= situacao=ok t_perm_h= passos= seg=`, gravado
num `.tmp` renomeado). Código de saída 0 em sucesso, 1 se o permanente não foi atingido ou o FS ficou abaixo de F₀, 2
para parâmetros inválidos.

No cluster (SLURM, pasta `cluster/`; os módulos em `Projects/Cluster/modulos.sh`):

```bash
cd Projects/V4_TaludeLiu/cluster && source ../../Cluster/modulos.sh && ./compilar_v4.sh     # no nó de login
mkdir -p logs
CASO=s1 ANGULO=0 FS_TOL=0.002 PRIMEIRA=1 ULTIMA=2000 DET=1 sbatch job_v4.sbatch
CASO=s2 ANGULO=0 PRIMEIRA=1 ULTIMA=2000 sbatch job_v4.sbatch                 # roda o det e reusa o seu permanente
CASO=s3 ANGULO=0 PRIMEIRA=1 ULTIMA=2000 sbatch job_v4.sbatch
CASO=s3 ANGULO=0 RHO=-0.5 PRIMEIRA=1 ULTIMA=2000 sbatch job_v4.sbatch        # com a correlação cruzada do artigo
CASO=s2 ANGULO=90 RHO=-0.5 PRIMEIRA=1 ULTIMA=2000 sbatch job_v4.sbatch
./estatisticas_v4.py resultados/s3_a0 resultados/s3_a0_rho-0.5                # média, desvio, beta, P_f, área, modos
```

`rodar_v4.sh <primeira> <última>` (chamado pelo job) roda cada realização em `resultados/<COMB>/r<nnnnn>/`, retomável
(pula as que têm `resultado.txt`), com `PARALELO` processos de 1 CPU; `COMB` padrão `<CASO>_a<ANGULO>[_rho<RHO>]`. No
cenário 2 o determinístico roda antes e as realizações leem o seu `pw.csv`. Tempo por realização (1 CPU, malha de
0,5 m): ~1,5 min nos cenários 1 e 3, ~1 min no 2; 2 000 realizações de um cenário em 24 núcleos: ~2 h.

## Testes (`V4_testes`, Catch2)

| teste | marca | o que compara | tolerância e justificativa | resultado (tempo) |
|---|---|---|---|---|
| retenção e k_r | [V4] | S da biblioteca × eq. 7 (com C(ψ)); k_r = S^5,57 × FXH (ψ_aev = 2,58 kPa) de 3 a 20 kPa; o patamar k_r = q/k_s | \|ΔS\| ≤ 0,002 (C ≥ 0,9989); \|Δlog₁₀ k_r\| ≤ 0,05 (o resíduo do ajuste é ≤ 9 %); patamar a 5 % | passou (0,2 s) |
| seções X-X e Y-Y | [V4] | p_w do permanente (malha de 1 m) × Fig. 3 digitalizada a cada 0,5 m; balanço de água | 1,5 kPa abaixo do lençol (hidrostático; leitura e malha); 2 kPa acima (o ψ_aev desconhecido move o patamar da FXH de 12,5 a 15,6 kPa); lençol a 0,15 m; balanço a 1 % da chuva | passou (23 s) |
| FS determinístico | [V4] | FS do SRM × 1,525 (Bishop); modo e área da massa × Fig. 2 | ±6 % (SRM × equilíbrio limite em taludes c'-φ', ~5 %, mais ~1 % da sucção menor); ruptura de base; área a 20 % de 44,5 m² | passou (junto com o anterior) |
| FS por arc-length | [V4] | FS do SRM por arc-length (`FactorSRMArc` de `Comum/RahardjoDriver.h`, o código de 2f92365 do `LearningMohrCoulomb2`), com as mesmas poropressões congeladas, × 1,525 e × a bissecção na mesma malha | ±6 % (as parcelas do teste anterior); arco − bissecção entre −0,5 % e 0,01 + 0,5 % (a bissecção fica até FS_TOL = 0,01 abaixo do limite discreto, o arco até 0,46 % abaixo, o pior caso do estudo de 43 casos de 2f92365) | passou: 1,5932 (22 incrementos, 4 cortes, 122 iterações no caminho, 13 s) contra 1,59375 da bissecção (135 iterações, com as do estado geostático; 13 s) |
| campos de KL | [V4] | qual propriedade varia em cada cenário (Tabela 3); L = √π θ/2 e a troca com α = 90°; média e desvio de ln c', ln φ', ln k num ponto em 400 realizações; correlação com `RHO=-0.5`; reprodutibilidade | 4 erros padrão (desvio esperado = o alvo × √(variância de KL no ponto)); correlação ±0,15 | passou (4 s) |
| fumaça do Monte Carlo | [V4] | uma realização s3 (α = 0) e uma s2 (α = 90°, permanente lido) pelo contrato (variáveis de ambiente, arquivos e chaves) | FS na faixa dos FS dos cenários 2 e 3 da Fig. 11 (1,0 a 2,3) | passou (36 s) |
| malha de 0,5 m | [V4][.lento] | os três testes determinísticos (seções, FS por bissecção e por arc-length) na malha do artigo | as mesmas | passou: arc-length 1,5825, bissecção 1,5875 (141 s) |
| Monte Carlo curto | [V4][.lento] | 16 realizações de s1 e de s2 (com `RHO=-0.5`) a α = 0, malha de 1 m: média e desvio do FS × Fig. 11 | média: a do artigo vezes o viés determinístico (FS_det/1,525), a 3 erros padrão + 0,02; desvio entre 0,5 e 2 vezes o do artigo | passou (550 s) |

Total dos padrão: ~70 s (6 testes, 162 verificações, todos passaram). Os demorados: `V4_testes "[.lento]"` (~10,5 min,
2 testes, 123 verificações, todos passaram). Tempos nesta máquina, 1 CPU compartilhada.

## Resultados: artigo × nosso

### Determinístico (q/k_s = 0,01)

| grandeza | artigo | nosso (h = 0,5 m) | nosso (h = 1 m) | diferença |
|---|---|---|---|---|
| p_w na crista, X-X z = 10 m | −14,6 kPa (Fig. 3a) | −13,15 kPa | −13,15 kPa | sucção 1,45 kPa (10 %) menor |
| p_w no terreno, Y-Y z = 7,5 m | −14,15 kPa (Fig. 3b) | −12,91 kPa | −12,92 kPa | sucção 1,24 kPa (9 %) menor |
| patamar k_r(ψ) = q/k_s (gradiente unitário) | ~14,5 kPa (Fig. 3) | 13,34 kPa (k_r = S^5,57); FXH com ψ_aev = 2,58 kPa: 13,43 kPa | — | −8 % |
| lençol (p_w = 0) em X-X | z = 5,12 m | 5,13 m | 5,12 m | 0,01 m |
| lençol em Y-Y | z = 4,55 m | 4,54 m | 4,53 m | −0,02 m |
| p_w na base, X-X z = 0 | 49,1 kPa | 49,77 kPa | 49,74 kPa | +0,7 kPa |
| p_w na base, Y-Y z = 0 | 44,1 kPa | 43,91 kPa | 43,86 kPa | −0,2 kPa |
| maior diferença acima / abaixo do lençol, X-X | — | 1,47 / 0,67 kPa | 1,46 / 0,64 kPa | |
| maior diferença acima / abaixo do lençol, Y-Y | — | 1,24 / 0,24 kPa | 1,23 / 0,28 kPa | |
| FS | 1,525 (Bishop) | **1,588** (SRM) | 1,594 | **+4,1 %** (+4,5 %) |
| FS, SRM por arc-length | 1,525 (Bishop) | 1,582 | 1,593 | +3,8 % (+4,5 %); −0,32 % (−0,04 %) contra a bissecção |
| superfície crítica: entrada, saída, ponto mais baixo | x = 6,8 m, x ≈ 21 m, (15,9; 2,89) m (Fig. 2, digitalizado) | x = 7,1 m, 21,4 m, z = 3,44 m (|du| ≥ 0,5 max); x = 6,4 m, 21,9 m, z = 3,01 m (≥ 0,2) | — | a mesma ruptura de base |
| D_f e modo | 7,1 m, base | 6,56 m, base (6,99 m com o limiar 0,2) | 6,56 m, base | |
| área da massa que desliza | 44,5 m² (Fig. 2, digitalizado); 40,1 m² (média do cenário 1, Fig. 8a) | 40,7 m² (limiar 0,5); 46,3 (0,3); 49,7 (0,2); 54,5 (0,1) | 41,3 m² | −9 % e +2 % |
| balanço no permanente (m²/s por metro) | — | chuva 3,12×10⁻⁶, entra em BC 4,86×10⁻⁶, sai em ED 7,98×10⁻⁶, resíduo < 10⁻¹⁹ | idem | |

Discussão:

- **Lençol e zona saturada:** coincidem com a Fig. 3 (lençol a 0,02 m, p_w a menos de 0,7 kPa). A posição do lençol
  é a de um aquífero livre com recarga: Dupuit com R/K = 0,01 dá h(10 m) = 5,13 m, o nosso valor.
- **Sucção do patamar:** a nossa é ~1,3 kPa (9-10 %) menor que a do artigo, nas duas seções. Ela é fixada pela função de
  permeabilidade (no patamar o gradiente é unitário e k_r(ψ) = q/k_s): a FXH que integramos com os parâmetros do
  artigo dá 13,4 kPa com ψ_aev = 2,58 kPa e de 12,5 (ψ_aev → 0) a 15,6 kPa (ψ_aev = a = 5 kPa); o ajuste S^5,57 a
  reproduz a 1 %. A diferença vem portanto do ψ_aev (ou de outro detalhe da FXH do SEEP/W, que o artigo não informa),
  não do ajuste nem do transiente: o artigo é compatível com um ψ_aev ~4 kPa. Não ajustamos ψ_aev ao artigo.
- **FS:** 4 % acima de Bishop. O mecanismo do SRM é a mesma ruptura de base da Fig. 2 (entrada a 0,3-0,4 m, saída a
  ~0,4 m, ponto mais baixo 0,1-0,5 m mais alto, conforme o limiar). A diferença é a esperada entre SRM em elementos
  finitos (Mohr-Coulomb associado, critério de não convergência com 25 iterações) e o equilíbrio limite de Bishop
  simplificado (que despreza as forças cisalhantes entre fatias); a nossa sucção menor puxaria o FS ~1 % para baixo. A
  malha quase não muda o FS (1,594 com 1 m, 1,588 com 0,5 m). O SRM por arc-length do commit 2f92365 (o maior F do
  caminho de equilíbrio, em vez do último F convergido com F prescrito) dá 1,593 e 1,582: 0,04 % e 0,32 % abaixo da
  bissecção. Como o F da bissecção convergiu, o limite discreto está acima dele, e a diferença é o erro da parada do
  arc-length (dF < 0,001 por incremento), dentro do pior caso de −0,46 % do estudo de 2f92365. Com a malha de 1 m o
  custo é o mesmo: 13 s, e 122 iterações no caminho contra 135 da bissecção, contadas com as do estado geostático.
- **Área da massa:** depende do limiar em |du| (de 41 a 55 m² entre 0,5 e 0,1), porque a banda de cisalhamento do MEF
  tem 1-2 elementos de espessura; com 0,5 (o meio da banda) a área fica entre as duas referências do artigo, que já
  diferem entre si (44,5 m² da superfície da Fig. 2 e 40 m² da média do cenário 1).

### Monte Carlo (α = 0°, q/k_s = 0,01)

Monte Carlo curto do teste `[.lento]` (16 realizações por cenário, malha de 1 m, KL com 1 000 termos; o cenário 2
com `RHO=-0.5`, a correlação cruzada do artigo, e o permanente do determinístico). O artigo: 2 000 realizações, Bishop,
núcleo gaussiano.

| grandeza | artigo | nosso (N = 16) | diferença |
|---|---|---|---|
| cenário 1 (k_s): média do FS | 1,52 (calculada da Tabela 4 e da Fig. 11; mediana ~1,51) | 1,588 (mediana 1,594) | +4,5 %, o viés determinístico SRM/Bishop desta malha (+4,5 %) |
| cenário 1: desvio do FS | 0,041 (Fig. 11) | 0,027 | −34 % (erro amostral ~18 %; a resolução do SRM, 0,006 com `FS_TOL=0.01`, já pesa: use `FS_TOL=0.002`) |
| cenário 1: β (eq. 13) | 15,27 (Tabela 4) | 26,9 | β ≈ ln μ/COV: o desvio menor quase dobra o β |
| cenário 1: área da massa | 40,1 ± 1,0 m² (Fig. 8a) | 41,2 ± 0,8 m² | +3 % |
| cenário 2 (c', φ', ρ = −0,5): média do FS | 1,49 (calculada; mediana 1,47-1,49) | 1,597 (mediana ~1,62) | +7 % (+2,6 % descontado o viés determinístico) |
| cenário 2: desvio do FS | 0,149 (Fig. 11) | 0,115 | −23 % (erro amostral ~18 %) |
| cenário 2: β | 3,95 (Tabela 4) | 6,5 | o desvio menor e a média maior |
| cenário 2: área da massa | 48,0 ± 19,9 m² (Fig. 8a) | 48,9 ± 12,6 m² | média +2 %, desvio −37 % |
| modo de ruptura | base em 1 990 (s1) e 1 360 (s2) de 2 000 (Fig. 10a) | base em 16 de 16 nos dois | o artigo tem 32 % de rupturas de pé no s2 |
| P_f | 0 (s1), 4,08×10⁻⁵ (s2) | 0/16 e 0/16 | não resolvível (abaixo) |

Leitura: a média do cenário 1 é o determinístico (k_s mexe pouco no FS) e o desvio do FS com k_s aleatório é ~2/3 do do
artigo: os nossos campos de k_s mudam o patamar de sucção (de 10 a 14,5 kPa na crista em 4 realizações) e o lençol
(±0,1 m em X-X, compatível com Dupuit), mas o FS do SRM responde pouco. No cenário 2 a média fica acima da do artigo
mesmo descontado o viés do SRM, e o desvio abaixo; com N = 16 as duas diferenças estão perto do erro amostral (o erro
padrão da média é 0,029 e o do desvio ~18 %). Os cenários com campos independentes (o padrão) e o cenário 3 ficam para
o Monte Carlo completo do cluster (`cluster/`); com ρ = 0 o desvio do FS deve crescer (a Fig. 19 do artigo mostra P_f
25 a 50 vezes maior com ρ ≈ 0 que com ρ = −0,5).

Os P_f do artigo (0 no cenário 1, 4×10⁻⁵ e 1,2×10⁻⁴ nos cenários 2 e 3, Tabela 4) **não podem ser resolvidos** por
contagem com o número de realizações que rodamos (nem com 2 000: o intervalo de Wilson de 95 % para 0 falhas em 2 000
vai até 1,9×10⁻³); o artigo também os calcula pela lognormal ajustada (eq. 12), que `cluster/estatisticas_v4.py`
informa junto com o β da eq. 13. Como β ≈ ln μ/COV, um viés de 4 % no FS (o do SRM × Bishop) muda β em ~0,3 com
COV = 0,1.

Diferenças de modelo que afetam a comparação estocástica: (1) o núcleo exponencial no lugar do gaussiano (a mesma
escala de flutuação: a variância da média em comprimentos grandes é a mesma, mas o exponencial tem mais variação em
pequena escala); (2) campos independentes por padrão (o artigo usa ρ = −0,5 entre c' e φ', que reduz o desvio do FS:
na Fig. 19 do artigo, P_f vai de 1,2×10⁻⁴ com ρ = −0,5 a 3×10⁻³ e 6×10⁻³ com ρ = −0,1 e +0,1; `RHO=-0.5` reproduz a
correlação do artigo); (3) SRM × Bishop (viés de +4 % no determinístico); (4) a sucção ~10 % menor (a FXH aproximada).

## Arquivos

- `problema.h`: parâmetros (`Config`, do ambiente), malha e contornos (`Talude`), campos de KL (`AplicaCampos`), o
  transiente até o permanente (`Permanente`), p_w nas seções (`PwEm`), o SRM com nova tentativa e o mecanismo (`FS`) e a
  execução do contrato (`RodaExecucao`).
- `main.cpp`: o executável do contrato; `testes.cpp`: os testes.
- `malhas/`: `talude_liu.geo` e as malhas de 0,5 m e 1 m.
- `dados/`: `fig3_secao_xx.csv`, `fig3_secao_yy.csv` (Fig. 3), `fig2_superficie_critica.csv` (Fig. 2), todos
  digitalizados, aproximados; `referencias.csv` (os valores do artigo usados na comparação, com a fonte).
- `scripts/`: `ajuste_kr.py` (a FXH, o ψ_aev pela tangente e o ajuste de p), `digitaliza_fig3.py` e
  `digitaliza_fig2.py` (a digitalização das figuras).
- `cluster/`: `compilar_v4.sh`, `rodar_v4.sh`, `job_v4.sbatch` e `estatisticas_v4.py`.
