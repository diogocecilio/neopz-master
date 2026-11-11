# KLCilindroMika: Karhunen-Loève de Galerkin do NeoPZ no meio cilindro de Mika et al. (2021) e no talude de Ceron (2025)

Dois programas usam a expansão de Karhunen-Loève (KL) de Galerkin do NeoPZ (`TPZKarhunenLoeveMat` com
`TPZInterpolationSpace::CalcStiffB` e `CalcStiffC`) e comparam o resultado com uma referência:

- `KLCilindroMika` (`main.cpp`): o Exemplo 1-1 da Seção 6.2 de Mika, Hughes, Schillinger, Wriggers e Hiemstra (2021,
  CMAME 379:113730). Campo com núcleo exponencial G(x, x') = exp(-|x - x'| / 5), σ² = 1, no meio cilindro vazado
  8 ≤ ρ ≤ 10, 0 ≤ θ ≤ π, 0 ≤ z ≤ 15. Os 20 maiores autovalores são comparados com a coluna "Galerkin 2" da Tabela 4,
  a soma de todos os autovalores com 270π (o volume) e os tempos são medidos. O mesmo programa roda a malha grossa do
  notebook `randomfied-v2.nb` (Mathematica) e compara com os 20 autovalores guardados nele.
- `KLTaludeCeron` (`talude.cpp`): a mesma KL de Galerkin em 2D com o núcleo exponencial separável de Vargas Ceron,
  Cecílio, Linn e Maghous (2025, eq. 64), G = exp(-|x - x'| / Lx - |y - y'| / Ly), Lx = 20 m, Ly = 2 m. A referência é
  a forma fechada no retângulo que o SlopeRandomField usa (`../SlopeRandomField/KLSeparavel.h`, `klsep::KL2D`):
  - (a) no retângulo [0,25] × [0,10]: autovalores numéricos contra os da forma fechada (produtos dos autovalores 1D);
  - (b) no talude (0,0) (25,0) (25,5) (15,5) (10,10) (0,10), área 187,5 m²: a KL numérica no próprio talude (como no
    Ceron, "from discretization of the actual soil mass slope") contra a forma fechada do retângulo [0,25] × [0,10]
    avaliada no talude (como no SlopeRandomField), com o mesmo número de termos M.

## Método

- O problema é ∫_D G(x, x') φ(x') dx' = λ φ(x). Com φ = Σ v_j N_j (funções de forma H1 contínuas), fica C v = λ B v,
  com B_ij = ∫ N_i N_j e C_ij = ∫∫ N_i(x) G(x, x') N_j(x') dx dx'. O núcleo do `TPZKarhunenLoeveMat` não tem a
  variância: os autovalores são os de σ² = 1.
- B: `CalcStiffB` de cada elemento, com a regra do elemento (ordem 2P).
- C: para cada par de elementos (e1, e2), `TPZKarhunenLoeveMat::ContributeC` com os dados de cada elemento nos pontos
  da sua regra de integração (`TPZInterpolationSpace::KLPointData`: coordenadas, peso × |det J| e funções de forma),
  calculados uma vez por elemento. `CalcStiffC(jel, ce)` faz a mesma conta para um par; `VALIDAR=1` no `KLCilindroMika`
  confere uma contra a outra. Só os pares e1 ≤ e2 são calculados (C é simétrica).
- B e C são matrizes densas. A montagem de C usa `THREADS` threads, cada uma somando na sua cópia de C. Os autovalores
  saem do Eigen: Cholesky B = L Lᵀ e `SelfAdjointEigenSolver` de L⁻¹ C L⁻ᵀ. Os autovetores v = L⁻ᵀ y são
  B-ortonormais (∫ φ_k² = 1).
- Quadratura: cada elemento do par usa a sua regra de Gauss. Os pares sem nó em comum, os pares com nó em comum e o par
  de um elemento com ele mesmo podem ter regras diferentes (`NG`, `NG_VIZINHO`, `NG_PROPRIO` no `KLCilindroMika`, em
  pontos por direção; `ORDEM`, `ORDEM_VIZINHO`, `ORDEM_PROPRIO` no `KLTaludeCeron`, ordem da regra do triângulo).
  O núcleo exponencial não é derivável em x = x', e o separável também não em x = x' e em y = y' separadamente. A regra
  de Gauss não integra exatamente os pares em que isso ocorre dentro do domínio de integração.
- Com a integração exata, os autovalores de Galerkin seriam limites inferiores dos exatos (princípio do min-max). Onde
  a referência é exata (a forma fechada no retângulo) ou convergida bem abaixo do erro em questão (as referências de
  Nyström das seções 1, 2 e 4), autovalores numéricos acima dela vêm da quadratura.

## Defeitos corrigidos na biblioteca

`Material/TPZKarhunenLoeveMat.h/.cpp`:
- `Dimension()` devolvia 2 para qualquer malha; agora devolve a dimensão dada no construtor (ou em `SetDim`).
- O núcleo do tipo 2 (gaussiano isotrópico) não tinha `break` e devolvia o do tipo 3 (exponencial separável).
- Um tipo de núcleo não definido devolvia 0; agora para com erro.
- Os construtores sem os parâmetros do núcleo deixavam `fLx`, `fLy`, `fLz`, `ftype` e `fDim` sem valor, e o construtor de
  cópia não copiava `fDim`.
- `NSolutionVariables`, `SetData` e `SetLinearContext` chamavam a si mesmos (recursão infinita).
- `Solution` sobrescrevia `data.sol[0]` com o caso de carga pedido (variáveis `vec1` a `vec5`).
- `AutocorrelationFunc` recebia os pontos por valor (duas cópias de vetor a cada chamada); agora por referência constante.
- Novo `ContributeC(TPZKLPointData, TPZKLPointData, ek)`: a integral dupla de um par a partir dos dados dos pontos, com
  O(np1·np2·n2 + np1·n1·n2) operações.

`Mesh/pzinterpolationspace.h/.cpp`:
- `CalcStiffC` usava a regra de integração do primeiro elemento também no segundo (errado para elementos de tipos ou
  ordens diferentes), a matriz do par tinha nshape(e1) colunas em vez de nshape(e2), e os dados do segundo elemento eram
  calculados de novo para cada par de pontos (np1·np2 chamadas de `ComputeRequiredData`; agora np1 + np2).
- Novos `CalcStiffC(jel, regra1, regra2, ce)` e `KLPointData(regra, dados)`.
- `CalcStiffB` e `CalcStiffC` não calculam mais a solução nos pontos (`fNeedsSol = false`).

## Compilar

O projeto entra no build do NeoPZ com `BUILD_PROJECTS=ON` e precisa do Eigen3 (≥ 3.3, `find_package(Eigen3)`). Sem o
Eigen3, o CMake deixa o projeto de fora e os outros configuram normalmente. No build do SlopeRandomField:

```bash
cd ../SlopeRandomField
./compilar.sh                                                   # biblioteca e SlopeRandomField (uma vez)
cmake --build build --target KLCilindroMika KLTaludeCeron -j 2  # os dois programas
```

Os executáveis ficam em `../SlopeRandomField/build/Projects/KLCilindroMika/`. O `KLTaludeCeron` lê as malhas do talude de
`../LearningMohrCoulomb2/cluster/malhas` (caminho gravado na compilação; `MALHAS_DIR` muda).

## Rodar

```bash
B=../SlopeRandomField/build/Projects/KLCilindroMika
# meio cilindro (as linhas da tabela 1)
MALHAS=8x1x4,12x2x6,16x2x8 $B/KLCilindroMika          # geometria exata, ordem 2, regra do elemento (3 pontos)
MALHAS=24x3x12 $B/KLCilindroMika                      # 8575 equacoes: 6 min, 1,7 GB
NG=5 MALHAS=8x1x4,12x2x6,16x2x8,24x3x12 $B/KLCilindroMika   # 5 x 5 x 5 pontos de Gauss em todos os pares
VALIDAR=1 MALHAS=4x1x2 $B/KLCilindroMika              # confere CalcStiffC, a simetria e o caso 2D
# malha do notebook (tabela 2)
GEOMETRIA=reta MALHAS=15x1x3 $B/KLCilindroMika        # e com NG=5
python3 confere_notebook.py q2 5 15x1x3               # conferencia independente (numpy/scipy): q2 ou s20, pontos, malha
# retangulo (tabela 3) e talude (tabelas 4)
PARTE=a $B/KLTaludeCeron                              # 5 malhas, regra do elemento
PARTE=a ORDEM=10 $B/KLTaludeCeron                     # regra de 25 pontos em todos os pares (100x40: 3 min)
PARTE=b $B/KLTaludeCeron                              # 5 malhas do talude, regra do elemento
PARTE=b ORDEM=10 TALUDES=ceron_h5_a45_fina.msh:1,ceron_h5_a45.msh:2,ceron_h5_a45_h025.msh:2 $B/KLTaludeCeron
```

`KLCilindroMika`:

| variável | padrão | |
|---|---|---|
| `P` | 2 | ordem dos elementos (1 a 3) |
| `MALHAS` | `8x1x4,12x2x6,16x2x8` | elementos em θ × ρ × z, uma malha depois da outra (`16x2x8,24x3x12,32x4x16` com `P=1`; `15x1x3` com `GEOMETRIA=reta`) |
| `GEOMETRIA` | `exata` | `exata`: `TPZGeoBlend` com `TPZArc3D` nas arestas circunferenciais; `reta`: hexaedros `TPZGeoCube` com os nós em (ρ cos θ, ρ sin θ, z) |
| `NG`, `NG_VIZINHO`, `NG_PROPRIO` | 0 | pontos de Gauss por direção nos pares sem nó em comum, com nó em comum e no próprio elemento (`NG` 0: P + 1; `NG_VIZINHO` e `NG_PROPRIO` 0: o de `NG`) |
| `THREADS` | 2 | threads na montagem de C |
| `VALIDAR` | 0 | 1: confere `CalcStiffC` com o algoritmo anterior, a simetria e dois casos 2D |
| `CSV` | | arquivo CSV onde acrescentar uma linha por malha |

`KLTaludeCeron`:

| variável | padrão | |
|---|---|---|
| `PARTE` | `ab` | `a` (retângulo), `b` (talude) ou as duas |
| `LX`, `LY` | 20, 2 | comprimentos de correlação (m) |
| `RETANGULOS` | `25x10:1,50x20:1,100x40:1,25x10:2,50x20:2` | células nx × ny (dois triângulos cada, diagonais alternadas) : ordem |
| `TALUDES` | `ceron_h5_a45.msh:1,ceron_h5_a45_h025.msh:1,ceron_h5_a45_fina.msh:1,ceron_h5_a45.msh:2,ceron_h5_a45_h025.msh:2` | malhas do Gmsh 4.1 : ordem (a caixa dos nós deve ser [0,25] × [0,10] e a área dos elementos 187,5: a forma fechada no talude e os pontos aleatórios supõem o talude (0,0) (25,0) (25,5) (15,5) (10,10) (0,10)) |
| `M` | `50,200,1000,2000` | números de termos comparados no talude (os maiores que o número de equações ficam de fora) |
| `ORDEM`, `ORDEM_VIZINHO`, `ORDEM_PROPRIO` | 0 | ordem da regra do triângulo nos pares sem nó em comum, com nó em comum e no próprio elemento (`ORDEM` 0: 2P; `ORDEM_VIZINHO` 0: a de `ORDEM`; `ORDEM_PROPRIO` 0: a de `ORDEM_VIZINHO`; no máximo 21) |
| `PARES`, `SEMENTE` | 500, 1 | pares de pontos aleatórios no talude (cada conjunto) e semente |
| `REDE` | 4 | a rede de pontos de cada triângulo onde σ_M² é avaliada: (a/n, b/n), a + b ≤ n |
| `RET_MAX` | 20000 | termos da forma fechada usados no M equivalente |
| `THREADS` | 2 | threads na montagem de C |
| `CSV_RET`, `CSV_TALUDE` | | arquivos CSV dos resultados de (a) e (b) |

## Resultados

Máquina de 4 núcleos compartilhada; montagem de C com 2 threads; o Eigen usa 1 thread. A memória é o pico do processo
durante a malha.

### 1. Meio cilindro, geometria exata, ordem 2

Erro relativo contra Galerkin 2 (Mika, Tabela 4: NURBS p = 2, 3800 elementos, 6642 graus de liberdade):

| malha | elementos | equações | pontos de Gauss | λ1 | λ10 | λ20 | máx. 1-20 | soma/270π | C (s) | autovalores (s) | memória (MB) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 8x1x4 | 32 | 459 | 3³ | 5,8e-4 | 6,2e-3 | 1,8e-2 | 1,8e-2 | 0,9489 | 0,01 | 0,06 | 30 |
| 12x2x6 | 144 | 1625 | 3³ | 1,0e-4 | 1,1e-3 | 3,0e-3 | 3,0e-3 | 0,9607 | 0,16 | 1,57 | 88 |
| 16x2x8 | 256 | 2805 | 3³ | 3,8e-5 | 3,9e-4 | 1,1e-3 | 1,1e-3 | 0,9690 | 0,56 | 7,93 | 227 |
| 24x3x12 | 864 | 8575 | 3³ | 6,9e-6 | 7,2e-5 | 2,0e-4 | 2,0e-4 | 0,9782 | 6,49 | 336,98 | 1724 |
| 8x1x4 | 32 | 459 | 5³ | 9,5e-5 | 8,9e-4 | 2,0e-3 | 2,0e-3 | 0,8956 | 0,11 | 0,04 | 30 |
| 12x2x6 | 144 | 1625 | 5³ | 1,6e-5 | 1,5e-4 | 3,5e-4 | 3,5e-4 | 0,9301 | 2,12 | 1,54 | 91 |
| 16x2x8 | 256 | 2805 | 5³ | 5,3e-6 | 5,3e-5 | 1,3e-4 | 1,3e-4 | 0,9458 | 6,67 | 7,92 | 231 |
| 24x3x12 | 864 | 8575 | 5³ | 2,8e-7 | 2,8e-6 | 6,2e-6 | 6,2e-6 | 0,9634 | 78,24 | 338,92 | 1745 |

- Os autovalores convergem por cima para Galerkin 2 e para a referência de Nyström (tabela abaixo). Com a regra do
  elemento (3³ pontos), o erro máximo dos 20 modos contra Galerkin 2 cai de 1,8e-2 (459 equações) para 2,0e-4 (8575
  equações). Galerkin 1 (Mika, 1050 graus de liberdade) difere de Galerkin 2 de 1,6e-5 (λ1) a 4,3e-4 (λ20).
- Contra a referência de Nyström, com 5³ pontos em todos os pares o erro cai, na mesma malha, 6,0 a 6,3 vezes em λ1 e
  6,7 a 8,7 vezes no máximo dos 20 modos. 16x2x8 (2805 equações, 15,5 s) fica a 1,6e-4 da referência, abaixo do 24x3x12
  com 3³ pontos (2,2e-4). O 24x3x12 com 5³ pontos fica a 6,2e-6 de Galerkin 2 nos 20 modos (λ1: 2,8e-7) e a 3,3e-5 da
  referência (λ1: 1,2e-6). Galerkin 2 fica acima da referência em todos os 20 modos, de 9,7e-7 (λ1) a 2,7e-5 (λ20).
- O erro contra a referência cai com h⁴: de 16x2x8 para 24x3x12 (h dividido por 1,5 nas três direções), o de λ1 cai 5,0
  vezes com 3³ e com 5³ pontos (1,5⁴ = 5,06).
- A soma de todos os autovalores é o traço de B⁻¹C, a integral em D da diagonal do núcleo projetado no espaço da malha.
  Fica abaixo de 270π e sobe com a malha. Com 3³ pontos ela é maior que com 5³ na mesma malha (0,9690 contra 0,9458 em
  16x2x8).
- O tempo é o do problema de autovalores denso (todos os autovalores): 337 s com 8575 equações. A montagem de C levou
  6,5 s com 3³ pontos e 78 s com 5³. Mika (Tabela 4, um núcleo) dá 17,3 h de formação e montagem para Galerkin 2 e 6,21 s
  de solução (Lanczos).
- `VALIDAR=1` (4x1x2): `CalcStiffC` corrigido contra o algoritmo anterior em elementos iguais 1,4e-15; simetria
  C(e1,e2) = C(e2,e1)ᵀ 3,4e-16; a matriz montada pelos dados dos pontos contra `CalcStiffC` par a par 1,2e-16; 2D com
  triângulo e quadrilátero de ordem 2 e G ≈ 1: C_ij contra (∫ N_i)(∫ N_j) 1,1e-9 (a diferença é o 1 - G, de 1e-9).

**Referência independente (Nyström).** Conferência em numpy/scipy, fora deste diretório e sem o NeoPZ: Nyström com
Gauss-Legendre composto em (ρ, θ, z), pesos ρ dρ dθ dz, na forma simétrica W^½ G W^½. As simetrias θ → π − θ e
z → 15 − z separam o problema em 4 classes, cada uma no quarto 0 ≤ θ ≤ π/2, 0 ≤ z ≤ 7,5 com o núcleo somado nas
imagens; os maiores autovalores de cada classe saem do Lanczos (`eigsh`). O quarto tem m × 7m × 4m painéis em ρ × θ × z
(lados de 2/m, 2,02/m e 1,875/m) com q pontos por painel e direção: m = 1 a 4 com q = 2 (até 14336 pontos por classe),
m = 1 a 3 com q = 3 (até 20412) e m = 1 e 2 com q = 4. Os autovalores convergem com h⁴ (q = 2, m = 2, 3, 4: as
diferenças sucessivas de λ1 caem 5,80 vezes; h⁴ dá 5,94). A referência é a média de quatro extrapolações de Richardson
(h⁴ + h⁵ e h⁴ + h⁶, com q = 2 em m = 2, 3, 4 e com q = 3 em m = 1, 2, 3); a diferença entre as quatro é no máximo 4,2e-7
do valor (λ1: 5,9e-9). A soma dos traços das 4 classes é 270π.

| modo | referência (Nyström) | (Galerkin 2 − ref.)/ref. | (Galerkin 1 − ref.)/ref. |
|---|---|---|---|
| 1 | 162,796422 | 9,7e-7 | 1,7e-5 |
| 2 | 91,427879 | 1,8e-6 | 3,0e-5 |
| 3 | 57,564320 | 2,7e-6 | 4,7e-5 |
| 4 | 51,087472 | 3,0e-6 | 5,3e-5 |
| 5 | 38,794672 | 4,2e-6 | 7,1e-5 |
| 6 | 27,901128 | 5,6e-6 | 9,8e-5 |
| 7 | 25,053410 | 6,0e-6 | 1,1e-4 |
| 8 | 19,366437 | 8,1e-6 | 1,3e-4 |
| 9 | 16,154275 | 9,6e-6 | 1,7e-4 |
| 10 | 15,792576 | 9,9e-6 | 1,7e-4 |
| 11 | 15,143688 | 1,1e-5 | 1,7e-4 |
| 12 | 11,210756 | 1,4e-5 | 2,3e-4 |
| 13 | 10,176235 | 1,5e-5 | 2,7e-4 |
| 14 | 9,690825 | 1,6e-5 | 2,8e-4 |
| 15 | 8,052196 | 1,9e-5 | 3,2e-4 |
| 16 | 7,576465 | 2,1e-5 | 2,1e-4 |
| 17 | 6,719768 | 2,4e-5 | 3,9e-4 |
| 18 | 6,443754 | 2,5e-5 | 2,8e-4 |
| 19 | 6,174615 | 2,5e-5 | 4,4e-4 |
| 20 | 5,761164 | 2,7e-5 | 4,6e-4 |

Erro relativo do NeoPZ contra a referência (os 20 erros de cada conta são positivos):

| malha | equações | pontos de Gauss | λ1 | λ10 | λ20 | máx. 1-20 |
|---|---|---|---|---|---|---|
| 8x1x4 | 459 | 3³ | 5,8e-4 | 6,2e-3 | 1,8e-2 | 1,8e-2 |
| 12x2x6 | 1625 | 3³ | 1,0e-4 | 1,1e-3 | 3,0e-3 | 3,0e-3 |
| 16x2x8 | 2805 | 3³ | 3,9e-5 | 4,0e-4 | 1,1e-3 | 1,1e-3 |
| 24x3x12 | 8575 | 3³ | 7,8e-6 | 8,2e-5 | 2,2e-4 | 2,2e-4 |
| 8x1x4 | 459 | 5³ | 9,6e-5 | 9,0e-4 | 2,0e-3 | 2,0e-3 |
| 12x2x6 | 1625 | 5³ | 1,7e-5 | 1,6e-4 | 3,8e-4 | 3,8e-4 |
| 16x2x8 | 2805 | 5³ | 6,2e-6 | 6,2e-5 | 1,6e-4 | 1,6e-4 |
| 24x3x12 | 8575 | 5³ | 1,2e-6 | 1,3e-5 | 3,3e-5 | 3,3e-5 |

### 2. A malha do notebook (`randomfied-v2.nb`)

O notebook monta B e C com hexaedros retos de 20 nós (serendipity do NDSolve`FEM`, `MeshOrderAlteration` de ordem 2)
numa malha de 16 × 2 × 4 nós em θ × ρ × z (15 × 1 × 3 elementos), nós em (r cos θ, r sin θ, h). C usa
`ElementIntegrationPoints[HexahedronElement, 5]` em cada elemento do par, B usa a de 6. No NeoPZ, `GEOMETRIA=reta
MALHAS=15x1x3`: 45 hexaedros de ordem 2 (Q2, 27 funções), 651 equações. O domínio é o semi-anel poligonal, de volume
842,0423 (270π = 848,2300). `confere_notebook.py` refaz a conta com numpy/scipy, em Q2 e em serendipity de 20 nós
(408 equações), na mesma malha.

| modo | Galerkin 2 (Mika) | exato no semi-anel poligonal (Nyström) | NeoPZ Q2, 3³ | NeoPZ Q2, 5³ | serendipity 20 nós (Python), 5³ | notebook | (NeoPZ 5³ − exato)/exato | (notebook − exato)/exato |
|---|---|---|---|---|---|---|---|---|
| 1 | 162,7966 | 162,0687 | 162,1780 | 162,0871 | 162,0871 | 163,9705 | 1,1e-4 | 1,2e-2 |
| 2 | 91,4280 | 90,7761 | 90,8849 | 90,7945 | 90,7945 | 92,8766 | 2,0e-4 | 2,3e-2 |
| 3 | 57,5645 | 57,3273 | 57,4413 | 57,3425 | 57,3424 | 56,1743 | 2,7e-4 | -2,0e-2 |
| 4 | 51,0876 | 50,6386 | 50,7464 | 50,6569 | 50,6569 | 53,2844 | 3,6e-4 | 5,2e-2 |
| 5 | 38,7948 | 38,5583 | 38,6733 | 38,5742 | 38,5740 | 37,8918 | 4,1e-4 | -1,7e-2 |
| 6 | 27,9013 | 27,6118 | 27,7185 | 27,6300 | 27,6300 | 31,2591 | 6,6e-4 | 1,3e-1 |
| 7 | 25,0536 | 24,8632 | 24,9783 | 24,8797 | 24,8795 | 30,6434 | 6,6e-4 | 2,3e-1 |
| 8 | 19,3666 | 19,2931 | 19,4159 | 19,2944 | 19,2936 | 24,8308 | 7,2e-5 | 2,9e-1 |
| 9 | 16,1544 | 15,9718 | 16,0774 | 15,9898 | 15,9898 | 24,8281 | 1,1e-3 | 5,5e-1 |
| 10 | 15,7927 | 15,6507 | 15,7655 | 15,6677 | 15,6675 | 19,4531 | 1,1e-3 | 2,4e-1 |
| 11 | 15,1439 | 15,0666 | 15,1969 | 15,0711 | 15,0700 | 18,2605 | 3,0e-4 | 2,1e-1 |
| 12 | 11,2109 | 11,1397 | 11,2759 | 11,1479 | 11,1464 | 16,1117 | 7,3e-4 | 4,5e-1 |
| 13 | 10,1764 | 10,0747 | 10,1882 | 10,0918 | 10,0915 | 13,6925 | 1,7e-3 | 3,6e-1 |
| 14 | 9,6910 | 9,5735 | 9,6777 | 9,5913 | 9,5912 | 12,2071 | 1,9e-3 | 2,8e-1 |
| 15 | 8,0524 | 7,9915 | 8,1319 | 8,0031 | 8,0011 | 10,8495 | 1,5e-3 | 3,6e-1 |
| 16 | 7,5766 | 7,5500 | 7,4890 | 7,4691 | 7,4690 | 9,7552 | -1,1e-2 | 2,9e-1 |
| 17 | 6,7199 | 6,6541 | 6,7474 | 6,6586 | 6,6583 | 8,3436 | 6,9e-4 | 2,5e-1 |
| 18 | 6,4439 | 6,4088 | 6,3795 | 6,3464 | 6,3463 | 7,5551 | -9,7e-3 | 1,8e-1 |
| 19 | 6,1748 | 6,0982 | 6,2027 | 6,1155 | 6,1153 | 6,9425 | 2,8e-3 | 1,4e-1 |
| 20 | 5,7613 | 5,7109 | 5,8516 | 5,7254 | 5,7229 | 6,3018 | 2,6e-3 | 1,0e-1 |

- O Q2 do Python e o do NeoPZ dão os mesmos 20 autovalores e a mesma soma, em todos os dígitos impressos, com 3³ e com 5³
  pontos (soma 804,0538099 e 755,2080374).
- O serendipity de 20 nós dá os 20 autovalores a menos de 4,5e-4 dos do Q2 (5³ pontos). Nas cinco contas (Q2 com 3³ e 5³
  pontos em C; serendipity com 2³, 3³ e 5³), λ1 fica entre 162,09 e 162,40 e λ9 entre 15,99 e 16,31.
- Nessa malha, o NeoPZ fica a menos de 1,6e-2 de Galerkin 2 nos 20 modos (com 5³ pontos, a menos de 1,5e-2). Os valores
  guardados no notebook diferem de Galerkin 2 de -2,4e-2 a +5,4e-1 (modo 9: 24,83 contra 16,15). Os modos 8 e 9 do notebook
  (24,8308 e 24,8281) são quase iguais; em Galerkin 2 eles são 19,37 e 16,15.
- A coluna "exato no semi-anel poligonal" é a referência de Nyström da seção 1 refeita no domínio da malha reta (os 15
  hexaedros retos, Gauss-Legendre composto em cada um, só a simetria z → 15 − z; extrapolação com q = 2, m = 1, 2, 3 e
  q = 3, m = 1, 2; diferença entre as estimativas ≤ 6,2e-6 do valor). Esses autovalores ficam 0,35% a 1,2% abaixo dos do
  cilindro. Contra eles, o Q2 do NeoPZ fica de -8,1e-3 a +2,5e-2 com 3³ pontos e de -1,1e-2 a +2,8e-3 com 5³ (abaixo do
  exato só nos modos 16 e 18, com as duas regras); os valores guardados no notebook ficam de -2,0e-2 a +5,5e-1.
- O notebook guarda, logo depois dos autovalores, `Length[vec[[1]]]` = 92: os autovetores guardados têm 92 componentes. A
  malha 15 × 1 × 3 do código guardado tem 408 nós com hexaedros de 20 nós. Nenhuma das contas nessa malha reproduz os
  valores guardados. A malha de 8 × 2 × 2 nós (7 × 1 × 1 elementos de 20 nós) tem 92 nós, mas também não os reproduz
  (λ1 = 160,33 com 5³ pontos em C e 167,48 com 2³).

### 3. Retângulo [0,25] × [0,10], Lx = 20 m, Ly = 2 m

A forma fechada: λ = λ_i^x λ_j^y em ordem decrescente; λ1 = 57,1548, λ10 = 4,5114, λ20 = 1,6741. A fração de variância
retida, Σ λ / 250, é 0,7947 com 20 termos, 0,8849 com 50, 0,9556 com 200, 0,9865 com 1000 e 0,9920 com 2000. Nos 2000
maiores termos, o índice 1D vai até 166 em x e 477 em y.

Erro relativo dos autovalores numéricos contra a forma fechada (malhas estruturadas de triângulos; Σλ: soma dos M
maiores autovalores / 250):

| malha | P | equações | regra (ordem/pontos) | λ1 | λ10 | λ20 | máx. 1-20 | máx. 1-100 | Σλ, M=200 | Σλ, M=1000 | Σ todos | C (s) | threads | autovalores (s) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 25x10 | 1 | 286 | 2/3 | 2,9e-3 | 3,5e-3 | 1,7e-2 | 1,8e-1 | 4,6e-1 | 0,9120 | - | 0,9145 | 0,03 | 1 | 0,01 |
| 50x20 | 1 | 1071 | 2/3 | 7,1e-4 | 8,7e-4 | 1,1e-2 | 1,8e-2 | 3,1e-1 | 0,9414 | 0,9549 | 0,9550 | 0,50 | 1 | 0,46 |
| 100x40 | 1 | 4141 | 2/3 | 1,8e-4 | 2,2e-4 | 6,0e-3 | 6,0e-3 | 6,6e-2 | 0,9546 | 0,9744 | 0,9770 | 7,98 | 1 | 26,87 |
| 25x10 | 2 | 1071 | 4/6 | 1,6e-3 | 2,0e-3 | 1,4e-2 | 5,1e-2 | 2,5e-1 | 0,9543 | 0,9709 | 0,9710 | 0,10 | 1 | 0,46 |
| 50x20 | 2 | 4141 | 4/6 | 4,1e-4 | 5,0e-4 | 1,0e-2 | 1,1e-2 | 2,1e-1 | 0,9585 | 0,9816 | 0,9850 | 1,73 | 1 | 26,63 |
| 25x10 | 1 | 286 | 10/25 | 2,1e-4 | 2,4e-4 | -9,5e-3 | 1,0e-1 | 5,7e-1 | 0,8907 | - | 0,8923 | 0,54 | 2 | 0,01 |
| 50x20 | 1 | 1071 | 10/25 | 5,4e-5 | 6,5e-5 | -1,1e-2 | 1,1e-2 | 3,2e-1 | 0,9331 | 0,9443 | 0,9444 | 16,70 | 1 | 0,45 |
| 100x40 | 1 | 4141 | 10/25 | 1,4e-5 | 1,7e-5 | -1,1e-4 | 1,1e-4 | 6,5e-2 | 0,9513 | 0,9696 | 0,9718 | 142,80 | 2 | 27,42 |
| 25x10 | 2 | 1071 | 10/25 | 2,2e-4 | 2,7e-4 | -4,7e-3 | 4,7e-3 | 3,2e-1 | 0,9358 | 0,9479 | 0,9479 | 0,64 | 2 | 0,45 |
| 50x20 | 2 | 4141 | 10/25 | 5,5e-5 | 6,7e-5 | 1,6e-3 | 1,6e-3 | 4,9e-2 | 0,9520 | 0,9712 | 0,9736 | 19,34 | 1 | 26,60 |
| forma fechada | | | | | | | | | 0,9556 | 0,9865 | 1 | | | |

- Com a regra do elemento (ordem 2P), os autovalores ficam acima dos exatos e convergem com h²: λ1 tem erro 2,9e-3,
  7,1e-4 e 1,8e-4 em P1 com h = 1; 0,5; 0,25 m (lado das células). Em P2, 1,6e-3 e 4,1e-4 com h = 1 e 0,5 m.
- Com a mesma malha, a regra de 25 pontos (ordem 10) reduz o erro de λ1 13 vezes em P1 e 7,5 vezes em P2 (50x20). Em P1
  50x20, seis dos 20 primeiros autovalores (modos 8, 11, 12, 15, 17 e 20) ficam abaixo dos exatos. Em P1 100x40 com essa
  regra, o erro máximo dos 20 primeiros é 1,1e-4.
- O traço (Σ todos) fica abaixo de 250 e cresce com a malha; a regra do elemento o aumenta (em P2 50x20, 0,9850 contra
  0,9736 com 25 pontos). Com a regra do elemento, Σλ com M = 200 em P2 50x20 (0,9585) passa do valor exato (0,9556).
- Ordem da regra na malha 25x10, P1 (286 equações): o erro de λ1 vai de 2,9e-3 (ordem 2, 3 pontos) a 4,3e-4 (ordem 6),
  2,1e-4 (ordem 8), 8,0e-5 (ordem 16) e 5,5e-5 (ordem 21, 91 pontos); o traço / 250, de 0,9145 a 0,8915.
- Subir a ordem só nos pares com nó em comum e no próprio elemento muda pouco (25x10, P1, `ORDEM_VIZINHO=10`: λ1 com erro
  2,4e-3, contra 2,9e-3); com a ordem 10 em todos os pares, 2,1e-4. No núcleo separável, a dobra em x = x' (e em y = y')
  está em todos os pares de elementos cujas faixas de x (ou de y) se sobrepõem, tenham nó em comum ou não.

### 4. Talude: KL numérica no talude contra a forma fechada do retângulo

Malhas do Gmsh do talude do Ceron em `../LearningMohrCoulomb2/cluster/malhas` (triângulos; tamanho até 1 m da superfície
e a partir de 6 m dela: `ceron_h5_a45` 0,375 e 1,5 m, `ceron_h5_a45_h025` 0,25 e 1 m, `ceron_h5_a45_fina` 0,1875 e
0,75 m). A forma fechada é a do SlopeRandomField: o retângulo é a caixa dos nós da malha, [0,25] × [0,10], e os M termos
são os M maiores produtos λ_i^x λ_j^y.

Os 10 maiores autovalores da KL numérica no talude (P2, `ceron_h5_a45_h025`, ordem 10): 48,810; 26,730; 16,912; 10,825;
9,572; 6,824; 4,719; 4,519; 3,830; 3,349. Os da forma fechada no retângulo: 57,155; 36,231; 21,400; 13,120; 13,046;
8,570; 8,270; 5,939; 4,885; 4,511. Na KL exata, a soma de todos os autovalores é a área (187,5 no talude, 250 no
retângulo).

**KL exata do talude (Nyström).** Conferência independente em numpy/scipy, fora deste diretório e sem o NeoPZ: grade
uniforme de passo h = 5/n no talude (nós (ih, jh)), com as dobras do núcleo (x = x', y = y') nas arestas das células;
trapézio tensorial nos quadrados e regra dos vértices nos meios-quadrados ao longo da face x + y = 20; forma simétrica
W^½ G W^½ e Lanczos (`eigsh`) sem montar a matriz (o núcleo separável na grade do retângulo envolvente é Ky U Kxᵀ);
n = 10 a 320 (até 770 mil nós) e extrapolação de Richardson em h² e h⁴. No retângulo, onde a forma fechada é exata, o
mesmo esquema dá λ1 e λ10 com diferença menor que 1e-9 (n = 20, 40, 80) e 1 − Σ_{k≤200} λ_k / 250 com diferença de
3,6e-7 (n = 80, 160, 320). No talude, os 10 maiores autovalores são 48,81000; 26,73006; 16,91142; 10,82447; 9,57161;
6,82401; 4,71860; 4,51868; 3,82941; 3,34915 (as extrapolações com n = 40, 80, 160 e 80, 160, 320 diferem em até 1,5e-8 do
valor). Os da KL numérica P2 `ceron_h5_a45_h025` (6 algarismos, como o programa os escreve) diferem deles de -2,6e-5 a
+3,0e-5 com a regra de 25 pontos e de +4,3e-5 a +8,2e-4 com a regra do elemento.

**Erro médio quadrático de truncamento no talude**, 1 − média de σ_M²(x), σ_M²(x) = Σ_{k≤M} λ_k φ_k(x)². Na KL
numérica, a média é Σ_{k≤M} λ_k / 187,5 (∫ φ_k² = 1 no talude). Na forma fechada, é Σ_{k≤M} λ_k ∫_talude φ_k² / 187,5,
com a integral de cada modo em forma fechada na faixa 0 ≤ y ≤ 5 e em x, e Gauss-Legendre composto em y na faixa de
cima. Essa integral foi conferida contra Gauss 2D com 4000 pontos por direção (diferença ≤ 5e-13) e, no programa, contra
a quadratura nos pontos da malha P2 com a regra de ordem 10 (diferença ≤ 1,5e-5 na média de σ_M²). A coluna "todos" usa
todos os autovalores numéricos. Na KL exata (do talude, ou a do retângulo avaliada no talude), 1 − média de σ_M² é a
média no talude de E[(Z − Z_M)²]; para cada M, a KL exata do talude dá o menor valor entre as expansões de M termos
(entre elas, a forma fechada do retângulo restrita ao talude). Na KL numérica, o valor sai dos autovalores numéricos e
não inclui a diferença entre os modos numéricos e os exatos.

| KL | malha | P | equações | regra | todos | M = 50 | M = 200 | M = 1000 | M = 2000 |
|---|---|---|---|---|---|---|---|---|---|
| numérica no talude | ceron_h5_a45 | 1 | 609 | 2/3 | 0,0715 | 0,1149 | 0,0767 | - | - |
| numérica no talude | ceron_h5_a45_h025 | 1 | 1314 | 2/3 | 0,0500 | 0,1060 | 0,0608 | 0,0502 | - |
| numérica no talude | ceron_h5_a45_fina | 1 | 2301 | 2/3 | 0,0382 | 0,1021 | 0,0524 | 0,0390 | 0,0382 |
| numérica no talude | ceron_h5_a45 | 2 | 2320 | 4/6 | 0,0273 | 0,0995 | 0,0463 | 0,0286 | 0,0273 |
| numérica no talude | ceron_h5_a45_h025 | 2 | 5085 | 4/6 | 0,0190 | 0,0990 | 0,0421 | 0,0224 | 0,0198 |
| numérica no talude | ceron_h5_a45_fina | 1 | 2301 | 10/25 | 0,0444 | 0,1039 | 0,0560 | 0,0449 | 0,0444 |
| numérica no talude | ceron_h5_a45 | 2 | 2320 | 10/25 | 0,0419 | 0,1030 | 0,0547 | 0,0426 | 0,0419 |
| numérica no talude | ceron_h5_a45_h025 | 2 | 5085 | 10/25 | 0,0287 | 0,1002 | 0,0460 | 0,0305 | 0,0291 |
| forma fechada do retângulo | | | | | | 0,1150 | 0,0444 | 0,0135 | 0,0079 |
| exata do talude (Nyström) | | | | | 0 | 0,0993 | 0,0375 | - | - |

Termos da forma fechada do retângulo com o mesmo erro médio quadrático que a KL numérica com M termos:

| malha | P | regra | M = 50 | M = 200 | M = 1000 | M = 2000 |
|---|---|---|---|---|---|---|
| ceron_h5_a45 | 1 | 2/3 | 51 | 92 | - | - |
| ceron_h5_a45_h025 | 1 | 2/3 | 57 | 129 | 169 | - |
| ceron_h5_a45_fina | 1 | 2/3 | 60 | 159 | 240 | 247 |
| ceron_h5_a45 | 2 | 4/6 | 63 | 189 | 368 | 391 |
| ceron_h5_a45_h025 | 2 | 4/6 | 63 | 216 | 512 | 604 |
| ceron_h5_a45_fina | 1 | 10/25 | 59 | 144 | 198 | 200 |
| ceron_h5_a45 | 2 | 10/25 | 60 | 150 | 213 | 217 |
| ceron_h5_a45_h025 | 2 | 10/25 | 62 | 191 | 337 | 360 |
| KL exata do talude (Nyström) | | | 63 | 253 a 254 | - | - |

σ_M² na rede de 15 pontos de cada triângulo (vértices, pontos das arestas e interiores) e a covariância truncada
Σ_{k≤M} λ_k φ_k(x) φ_k(x') contra o núcleo em 500 pares de pontos uniformes no talude e em 500 pares próximos (o segundo
ponto a até 2 m em x e 0,5 m em y do primeiro). KL numérica: P2, `ceron_h5_a45_h025`, 5085 equações, ordem 10. Os valores
da forma fechada são calculados nos mesmos pontos. KL exata do talude (Nyström): σ_M² nos nós da grade com n = 160
(h = 3,125 cm; com n = 80 a diferença é ≤ 2e-3) e a covariância nos mesmos 1000 pares do programa (std::mt19937_64 com a
semente 1, reproduzido em Python), com os modos da grade n = 320 interpolados pela fórmula de Nyström (com n = 160, os
quatro erros são os mesmos com dois algarismos).

| M | KL | σ_M² mín. | σ_M² máx. | pares uniformes: erro máx. | erro médio | pares próximos: erro máx. | erro médio |
|---|---|---|---|---|---|---|---|
| 50 | numérica no talude | 0,769 | 0,992 | 9,2e-2 | 3,6e-3 | 1,1e-1 | 2,0e-2 |
| 50 | forma fechada do retângulo | 0,731 | 0,895 | 1,0e-1 | 4,2e-3 | 1,4e-1 | 2,5e-2 |
| 50 | exata do talude (Nyström) | 0,763 | 0,917 | 9,1e-2 | 3,4e-3 | 1,1e-1 | 1,9e-2 |
| 200 | numérica no talude | 0,896 | 1,069 | 4,8e-2 | 1,0e-3 | 5,4e-2 | 8,5e-3 |
| 200 | forma fechada do retângulo | 0,884 | 0,959 | 3,7e-2 | 8,0e-4 | 4,2e-2 | 4,9e-3 |
| 200 | exata do talude (Nyström) | 0,902 | 0,968 | 3,0e-2 | 6,5e-4 | 3,5e-2 | 3,8e-3 |
| 1000 | numérica no talude | 0,915 | 1,105 | 3,6e-2 | 7,6e-4 | 5,5e-2 | 7,7e-3 |
| 1000 | forma fechada do retângulo | 0,962 | 0,987 | 8,6e-3 | 1,0e-4 | 1,0e-2 | 6,9e-4 |
| 2000 | numérica no talude | 0,915 | 1,111 | 3,5e-2 | 7,5e-4 | 5,5e-2 | 7,7e-3 |
| 2000 | forma fechada do retângulo | 0,977 | 0,992 | 4,3e-3 | 3,7e-5 | 4,3e-3 | 2,8e-4 |

- Com M = 50, a KL numérica no talude tem erro médio quadrático de 0,099 a 0,104 em P1 `ceron_h5_a45_fina` e em P2, com
  as duas regras, contra 0,115 da forma fechada do retângulo: o retângulo precisa de 59 a 63 termos para o erro de 50
  termos da KL do talude. A KL exata do talude (Nyström) dá 0,0993 com M = 50; a forma fechada precisa de 63 termos
  para esse erro.
- Com M = 200, a KL numérica ainda muda com a malha (P2: 0,0547 com 2320 equações, 0,0460 com 5085); a forma fechada dá
  0,0444 e a KL exata do talude, 0,0375 (a forma fechada precisa de 253 a 254 termos para esse erro). Nenhuma das malhas
  usadas aqui chega a 0,0375 com M = 200 (o menor valor é 0,0421, P2 `ceron_h5_a45_h025` com a regra do elemento).
- Com M = 1000 e 2000, a KL numérica para no erro da malha: 1 − traço / área é 0,0287 na malha mais fina (P2, 5085
  equações, ordem 10); a forma fechada dá 0,0135 e 0,0079. A KL exata do talude não foi calculada com esses M; o seu
  erro é no máximo o da forma fechada.
- A regra do elemento (ordem 2P) dá erros menores que a de 25 pontos na mesma malha (M = 50: 0,0990 contra 0,1002;
  traço / área 0,9810 contra 0,9713 em P2 `ceron_h5_a45_h025`), como no retângulo, onde ela aumenta os autovalores. Com
  ela, o valor de M = 50 em P2 `ceron_h5_a45_h025` (0,0990) fica abaixo do da KL exata do talude (0,0993), que é o menor
  possível com 50 termos.
- A variância σ_M² da forma fechada fica abaixo de 1 em todos os pontos e o intervalo encolhe com M (0,977 a 0,992 com
  2000 termos). A da KL numérica passa de 1: com M = 2000, de 0,915 a 1,111 na malha mais fina, e de 0,848 a 1,165 em
  P1 `ceron_h5_a45_fina` com a regra do elemento; com todos os modos, σ² é a diagonal do núcleo projetado na malha.
- Os erros da covariância da forma fechada caem com M (pares próximos, média: 2,5e-2 com 50 termos, 2,8e-4 com 2000); os
  da KL numérica quase não mudam de M = 200 a 2000 (pares próximos, média: 8,5e-3 e 7,7e-3 na malha mais fina).
- Com M = 50 e M = 200, a KL exata do talude tem os quatro erros da covariância menores que os da forma fechada com o
  mesmo M (M = 200, pares próximos, média: 3,8e-3 contra 4,9e-3) e σ_M² entre 0,763 e 0,917 (M = 50) e entre 0,902 e
  0,968 (M = 200), abaixo de 1. Com M = 200, a KL numérica na malha mais fina tem erros maiores que os das duas
  (pares próximos, média: 8,5e-3).
- Tempo: a KL numérica na malha de 5085 equações levou 260 s (242 s no problema de autovalores com autovetores) e
  650 MB. A forma fechada com 20000 termos e as integrais de cada modo no talude levou 1,3 s.

## Conferências independentes

Feitas em numpy/scipy, fora deste diretório, sem o NeoPZ e sem o código dos dois programas:

- Meio cilindro, malha reta do notebook e talude: as referências de Nyström das seções 1, 2 e 4.
- Forma fechada do retângulo: as raízes 1D recalculadas com outra forma das equações (x sin x − ca cos x = 0 e
  x cos x + ca sin x = 0, x = ωa, por `brentq`) dão os mesmos autovalores; o Nyström 1D (trapézio, com a dobra do núcleo
  nos nós, extrapolado em h² e h⁴) concorda com eles a 6e-14 nos 12 primeiros modos de cada direção. Os produtos dão os
  valores da seção 3 (λ1 = 57,15478; λ10 = 4,51136; λ20 = 1,67408; índices 1D até 166 em x e 477 em y nos 2000 maiores).
- Integral de φ_k² no talude: pela transformada de Fourier exata do talude (dois retângulos e o triângulo da face, este
  pela diferença dividida de exp, Hermite-Genocchi), conferida com Gauss 2D de 400² pontos por parte nos modos 1, 8, 200
  e 2000 (diferença ≤ 3e-14). O erro médio quadrático da forma fechada no talude que sai dela difere do do programa em no
  máximo 2,7e-11 (M = 50, 200, 1000 e 2000), e os 29 valores de "termos da forma fechada com o mesmo erro" da seção 4 são
  os mesmos.
- Covariância truncada da forma fechada: os 1000 pares de pontos do programa, gerados de novo (std::mt19937_64 reproduzido
  em Python), dão os mesmos erros máximos e médios da tabela da seção 4 para M = 50, 200, 1000 e 2000.
- Os dois programas, compilados de novo, repetem os resultados das tabelas (malhas 8x1x4, 12x2x6, 16x2x8, 24x3x12 e a do notebook;
  as cinco malhas do retângulo com a regra do elemento; `ceron_h5_a45` P2 com ordem 10) em todos os dígitos dos arquivos
  CSV.

## Arquivos

- `main.cpp`: `KLCilindroMika`.
- `talude.cpp`: `KLTaludeCeron`.
- `confere_notebook.py`: a KL de Galerkin da malha do notebook em numpy/scipy (Q2 ou serendipity de 20 nós).
- `CMakeLists.txt`: os dois alvos, com o Eigen3 opcional.
