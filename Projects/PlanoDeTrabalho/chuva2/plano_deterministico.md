# Validação determinística da chuva (antes de qualquer estatística)

25/11/2026. Decisão do usuário: o determinístico precisa ser totalmente confiável antes da análise estatística; caso
contrário a estatística sai errada e distorcida. O plano estocástico (`plano_estocastico.md`) fica em stand-by até o fim
desta validação.

## Modelo e convenções (o que está sendo validado)

- Talude de Rahardjo et al. (2007): 45°, H_s = 10 m, H_w = 5 m abaixo do pé (NA inclinado de 7°), sucção inicial
  ≤ 75 kPa, Fredlund–Xing n = m = 1, θ_s = 0,45, k_w = k_s Θ⁴, c′ = 10 kPa, φ′ = φ^b = 26° (χ = 1), γ = 20 kN/m³.
- Solos: f10 (a = 10 kPa, k_s = 10⁻⁴ m/s = 360 mm/h), f50 (50 kPa, 10⁻⁵ = 36 mm/h), f100 (100 kPa, 10⁻⁶ = 3,6 mm/h).
- Variantes do transiente: (iii) desacoplado, esqueleto rígido, ar passivo (as hipóteses do artigo); (ii) acoplado,
  esqueleto elástico (E = 20 MPa, ν = 0,3), ar passivo; (i) acoplado, esqueleto elástico, ar ativo (p_a = 0 só na
  superfície).
- FS: SRM por elementos finitos (malha tri1 + 4 níveis adaptativos, refeitos do zero em cada instante) com p_w (e p_a)
  congeladas no instante.
- Chuva: `chuva_face=inclinada` (padrão: q = I por metro inclinado da superfície, o "unit flux" do SEEP/W) com a
  condição unilateral (face de percolação quando p_w chega a 0).

## Testes (casos_det_chuva.txt, 669 casos únicos, cada um nas três variantes)

Todos com o FS de hora em hora até 48 h (de meia em meia hora nas chuvas de até 6 h; até D + 24 h nas chuvas longas,
de 2 em 2 h na de 240 h), o diagnóstico do transiente em cada instante (`diag.csv`) e as pressões nos nós
(`pw.txt`, com p_a na variante i).

| teste | o quê | casos |
|---|---|---|
| T1 | Chuvas da Série B do Rahardjo, 24 h: f10 {9, 36, 54, 72, 90, 120, 180, 360, 900}, f50 {3,6; 9; 18; 27; 36; 45; 54; 90; 180}, f100 {0,9; 3,6; 9; 18; 36; 80} mm/h | 72 |
| T2 | O mesmo volume (864 mm) em 3, 6, 12, 24 e 72 h nos três solos (o A6 do relatório), mais as durações críticas t* = V/k_s fora dessa faixa (f10 em 1 h, f100 em 240 h) | 42 |
| T3 | Sensibilidade em torno da chuva base de cada solo (f10 90, f50 36, f100 9 mm/h; 24 h): c′ ∈ {5, 10, 15} kPa; φ′ ∈ {21, 26, 31}°; k_s ∈ {k_s/3, k_s, 3k_s}; anisotropia de k (k_h/k_v ∈ {2, 5, 10}, tensor diagonal); rotação do tensor (k_h/k_v = 5, θ ∈ {45, 90, 135}°; 135° = paralela à face); as duas juntas (k_h/k_v ∈ {2, 5, 10} com θ = 135°) | 126 |
| T4 | Verificações numéricas e de convenção nas chuvas base: malha do transiente (tri3, tri4, camadas de 0,25 m até 5 m), passo de tempo 0,05 h, malha do fator (nref 5) e a chuva na projeção horizontal (`chuva_face=vertical`) | 54 |
| T5 | Séries do relatório de outubro e do artigo (figuras e contraexemplos): Série A (α = 26,6; 33,7; 45; 63,4°) e A5 (α de 20 a 70° a cada 5°), Série C (H_w de 2,5 a 15 m), Série D (H_s = 5, 10, 20, 40 m), A2 (f100 a 63,4° e com H_s = 20 m), A3 (k₁ = 10 k_s girado de 0, 45 e 135° em todas as chuvas do f50), A4 (teto da sucção 75/150/300 kPa) | 375 |

A anisotropia determinística é só de k_s (decisão do usuário): a de c′ e φ′ existe só no campo aleatório (estrutura
de correlação girada) e fica para a etapa estocástica.

Código novo (montecarlo.cpp, sem mexer no material): `k_theta=` (direção de k₁ = alfa k_s, graus a partir do eixo x;
o material já tinha o tensor girado) e `diag=` (por instante: passos, Newton, nós na face, p_w e p_a extremos, |u| máx,
água infiltrada acumulada). Padrão inalterado.

Custo: ~12 s por FS no PC sozinho, ~60–90 s com 24 processos por nó; ~1 h por caso de 48 h. 669 casos em 4 nós:
~7–10 h, mais a cauda (H_s = 40 m, 240 h, α = 20°). Uma noite.

Scripts (Projects/SlopeRandomFieldRefinement/cluster): `gera_casos_det_chuva.py` → `casos_det_chuva.txt`,
`compilar_det_chuva.sh` (bin_det/), `submeter_det_chuva.sh`, `det_chuva.sbatch` (fila por travas; `resumo`),
`baixar_det_chuva.sh` (→ ~/projects/runs/fermi_det_chuva).

## O que extrair e investigar

1. **Convergência numérica (T4)**: FS(t) e FS mínimo contra a malha do transiente, o passo de tempo e a malha do
   fator. O FS inicial já muda com a malha do transiente (2,0933 com tri2 × 2,0978 com camadas de 0,25 m: as pressões
   iniciais são interpoladas da malha do transiente, e o perfil hidrostático com teto de 75 kPa tem uma quina). A malha
   padrão do transiente (tri2) tem só 319 nós; o relatório de outubro usava 0,75 m até 2 m de profundidade (8874
   equações). Critério: a escolha de malha/passo é aceita se o FS mínimo muda < 0,01 com o refinamento.
2. **A diferença de ~0,1 para o relatório de outubro na (iii)**: dois candidatos já vistos, testados em T4: a chuva por
   metro inclinado (41 % a mais de água na face a 45° do que na projeção horizontal do relatório) e a malha do
   transiente. Definir qual convenção de chuva é a física (a chuva cai na vertical: por metro inclinado chega I cos α)
   e qual reproduz o artigo (SEEP/W).
3. **Figuras do relatório de outubro** refeitas com os dados novos, lado a lado: FS mínimo × I com a Eq. 6 (as três
   variantes), FS(t) por variante, a tabela da Fig. 5b (6, 12, 24 e 48 h), campos de p_w e p_a, Séries A/C/D (Figs. 8,
   9, 10, Tabela 7), anisotropia e rotação (FS(t), FS mínimo × θ, infiltração), A3 (FS mínimo × I por θ), A4, A5, A6.
4. **Figuras novas**: tornado da sensibilidade por solo e variante; infiltração, escoamento, nós alagados, p_a máx e
   |u| máx no tempo; profundidade do mecanismo (banda do SRM) × profundidade da frente de umedecimento (de pw.txt);
   o instante do FS mínimo em todos os casos.
5. **Contraexemplos** (ser crítico com os resultados e com os artigos):
   - Rahardjo 2007: o patamar 0,81–0,82 e a intensidade-limite; "solos pouco permeáveis são seguros em chuvas curtas"
     (A2: agora com o esqueleto elástico a (i) pode passar abaixo de 1); "o controle é k_s" (com anisotropia k_s não é
     um número); "o NA inicial é secundário" (A4); "FS mínimo linear em α" (A5); intensidade × volume (A6); as
     categorias da Tabela 3 (elevação do NA).
   - O próprio relatório de outubro: o patamar FS ≈ 1,00 do modelo completo (era a dilatância do transiente plástico)
     e a explicação "SRM × Bishop" para os +0,05–0,15.
   - Os artigos lidos (notas em research_notes/Reprodução dos artigos novos/ e PlanoDeTrabalho): "o FS mínimo ocorre no
     fim da chuva" (Hu 2024, Cho 2014) contra o instante do mínimo em cada caso (no f100 o relatório achou o mínimo
     depois da chuva); a ruptura na frente de umedecimento do talude infinito (Cho 2014, X. Liu 2025, Jiang 2026)
     contra a profundidade do mecanismo 2D; a duração crítica.
6. **Critérios de aceitação antes do estocástico**: convergência numérica dentro de 0,01 no FS mínimo; a (iii)
   reproduzindo o artigo dentro da incerteza da digitalização (±0,03) ou com a diferença explicada; nenhuma
   inconsistência física sem explicação (FS subindo com mais água, transiente parado, Newton com muitas iterações,
   balanço de água); a escolha da convenção da chuva e da malha registrada.

## Estado

- 25/11: opções novas testadas (rotação e diag; fumaça de 8 casos incomuns e o sbatch local, sem erro); campanha pronta
  para o fermi.
