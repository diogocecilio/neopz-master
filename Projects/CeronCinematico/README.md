# CeronCinematico

Abordagem cinemática (limite superior, blocos rígidos em espiral logarítmica) de
Vargas Ceron et al. (2025), IJNAMG 49(11):2459-2491, em Python (numpy, scipy, matplotlib).

## Arquivos

| arquivo | o que faz |
|---|---|
| `cinematico.py` | o método: mecanismos clássicos (I) e (II) em forma fechada, espiral por trechos para campos heterogêneos, minimização, leitura das exportações do NeoPZ |
| `percolacao_analitica.py` | campo de velocidades semianalítico v'_opt do artigo (Seção 4.1, Eqs. 29-40), α = 1 |
| `permanente.py` | percolação permanente P1 (Eqs. 19-21) na malha da exportação |
| `figura.py` | figuras `figs/mecanismo_det.pdf` e `figs/mecanismo_real0.pdf` |
| `tests.py` | testes (a)-(e); imprime a tabela e gera as figuras |
| `dados/det`, `dados/real0` | exportações do `SlopeDrawdownMonteCarlo`: campos médios e realização 0. `linha_mc_nref5.csv`: a mesma realização com o GIM em NREF = 5 |

## Uso rápido

```bash
python3 -I tests.py          # cerca de 2 min, termina com "TODOS OS TESTES PASSARAM"
python3 -I figura.py         # só as figuras
```

```python
import cinematico as k
res = k.avaliar_realizacao("dados/real0")     # pasta <saida>_real<n>/ do NeoPZ
print(res.Gamma, res.mecanismo.classe, res.mecanismo.xA, res.mecanismo.xB)
print(res.n_normalidade, res.n_phi_nao_conv)  # avaliações descartadas; devem ser 0
```

`avaliar_realizacao` lê `vertices.csv`, `triangles.csv`, `snaps_pw.csv` (a última linha) e
`campos.csv`. A geometria padrão é a do rebaixamento: aresta da crista em (35, 35), H = 5 e
β = 45°, mas dá para mudar com `origem=`, `H=` e `beta_graus=`.

Para gerar uma realização nova (cerca de 30 s; `det=1 n=0` dá os campos médios, `nref=` o número de níveis do GIM):

```bash
PZ_THREADS=0 OMP_NUM_THREADS=1 SlopeDrawdownMonteCarlo n=1 primeira=<n> figuras=0 saida=r
```

## Convenções (as mesmas do artigo)

- Origem na aresta da crista, x para fora do talude e **y para baixo**.
- A crista é y = 0. O pé fica em T = (H/tan β, H).
- O excesso de poropressão é u = p − γw·y, com p total (compressão positiva).
  Na leitura: x' = x − 35, y' = 35 − y e u = p_w − 9,81·y'.
- O bloco gira em torno de C = (Cx, −Cy), e o ponto da curva é P(θ) = (Cx − r cos θ, −Cy + r sin θ).
- Γ = P_mr / (P_γ + P_u) (Eq. 54). Γ multiplica **todas** as forças de corpo:
  o peso γ' (= γ − γw) e a percolação −∇u. Isso equivale a multiplicar γ e p_w juntos, como faz o GIM do NeoPZ.

## O método geral (campos heterogêneos)

1. **Parâmetros:** (sB, θ2, rB).
   - B é o ponto onde a curva sai na superfície. sB é a posição de B medida a partir do pé:
     sB ≤ 0 na face, o que dá o mecanismo (I) (sB = 0 é "pelo pé", η = 1); sB > 0 no terreno à frente do pé, o mecanismo (II).
   - θ2 e rB = |CB| fixam o foco C. A única exigência sobre C é ficar acima da crista (Cy > 0); C pode ficar atrás da reta da face (Cy + tan β·Cx < 0), que é ar e é admissível.
2. **Marcha de B até A, para trás:** a curva é feita de arcos de log-espiral com o mesmo foco C, cada um com comprimento l = H/100.
3. **Escolha de φ_i em cada arco:** φ_i = máximo de φ(x) **ao longo do próprio arco**. É o que a normalidade pede ([U]·n ≥ |[U]| sin φ(x) em cada ponto do arco, Seção 5.2), e nada mais.
   - φ é P1 na malha. O arco é amostrado em 6 pontos (cerca de 1 cm), cada um localizado no triângulo exato (trifinder).
   - Entre duas amostras vizinhas em triângulos diferentes, o valor na aresta cruzada é limitado pelas extrapolações lineares dos dois triângulos. Soma-se ainda |∇φ| vezes a flecha do arco (h²/8r). O resultado é um limite superior de max φ(x) no arco.
   - Como o arco muda com φ_i, repete-se até φ_i parar de crescer. Sem convergência, o mecanismo é invalidado.
   - A marcha termina quando a curva cruza a crista. O último arco é cortado na interseção exata, e esse ponto é A.
4. **Dissipação:** usa c(x) e φ(x) locais em 2 pontos de Gauss por arco.
   - `dissipacao='exata'` (padrão): dP_mr = c·tan φ_i / tan φ(x) · r² dθ. É limite superior rigoroso porque φ(x) ≤ φ_i em todo o arco.
   - `'simples'`: dP_mr = c r² dθ. Com a regra do item 3, φ_i ≈ φ(x), e os dois modos ficam a 0,2% um do outro (realização 0).
   - Se em algum ponto de Gauss φ(x) > φ_i (mais de 1e-6 grau), o mecanismo é **invalidado** (motivo "normalidade violada"), e a minimização conta esses casos. Nada é cortado em silêncio.
5. **Bloco:** quadratura polar em torno de C, cobrindo todo o intervalo angular do contorno (espiral + crista, face e terreno até B). Em cada raio, os intervalos de ρ dentro do bloco saem dos cruzamentos com o contorno (um intervalo; dois quando C está atrás da reta da face).
   - P_γ = γ' ∫ ρ cos θ dA (exato em ρ).
   - P_u = ∫ −∇u·e_θ ρ dA (Gauss composto em ρ). O ∇u é o gradiente P1 de cada triângulo da malha.
   - Também se calcula P_u pela forma de contorno (div U = 0). Serve de verificação.
6. **Minimização:** primeiro uma grade de mecanismos clássicos (θ1, θ2, η|D) com o φ médio, avaliada com o campo real e arcos de H/25. Depois, Nelder-Mead em (sB, θ2, ln rB) a partir das 3 melhores partidas de cada classe.

## Resultados (tests.py)

| caso | Γ obtido | referência |
|---|---|---|
| (a) código geral × forma fechada, c e φ uniformes | dif. máx. 4e-10 | — |
| (a) mecanismo II com C atrás da reta da face (φ = 0,9°) | 2,4405 (geral e fechada); P_γ igual ao de um polígono independente | antes era recusado |
| (b) Cho coesivo, cu = 23 (mecanismo II, D = 5,97 m) | 1,3540 | artigo 1,354 (Cho LEM 1,356) |
| (b) Cho c-φ, H = 10 m (o texto do artigo diz 5 m; a Fig. 19 mostra 10 m) | Γ 1,7770; FS 1,2032 | artigo 1,777; 1,203 |
| (c) **transiente** NeoPZ (campos médios, t = 0,1 h) | **1,2791** (pelo pé, xA = −1,08 m) | GIM extrapolado 1,276 (NREF = 5: 1,2851). Sem referência do artigo |
| (c) percolação **permanente** na mesma malha | 1,414 (lateral u = 0); 1,429 (lateral sem fluxo) | artigo 1,336 (só informativo) |
| (d) realização 0 | 1,0806 ('exata'); 1,0788 ('simples') | GIM extrapolado 1,047 (NREF = 5: 1,0516) |
| (e) P_u hidrostático: u = −γw·y com γ' = 10,19 | 3,5540 | Γ seco com γ = 20: 3,5540 |
| (e) v'_opt, Fig. 9 (tracejada, α = 1), β = 45°, 60°, 75° | 1,905; 1,429; 1,025 | leitura da figura: 1,95; 1,50; 1,05 |
| (e) v'_opt, Fig. 8, H_crit com h_w/H = 1 | 4,83; 8,72; 8,13; 14,59 | leitura: 4,84; 8,6; 8,05; 14,7 |

Uma avaliação de Γ leva cerca de 0,6 ms com campo uniforme e cerca de 11 ms com a malha (6021 nós). A minimização completa leva cerca de 9 s por realização.

### Transiente × permanente × artigo

- **São problemas diferentes.** O 1,336 do artigo (Fig. 24; na Fig. 9, α = 1, a curva cheia −∇u'_FE passa por cerca de 1,33 em β = 45°) vem de percolação **permanente saturada** logo após um rebaixamento instantâneo (Eqs. 19-21, ∂Ω_v = ∅).
- A exportação `dados/det` é o fim (t = 0,1 h) de um **transiente não saturado** com armazenamento. Esse estado depende de escolhas que o artigo não tem: a duração do rebaixamento (TREBAIX = 0,1), a retenção e a compressibilidade, e a permeabilidade.
  - Na mesma malha, o campo permanente difere do transiente em até 23 kPa e dá Γ = 1,414.
  - Por isso o Γ do transiente (1,279) não é comparado com 1,336. Ele é comparado só com o GIM da mesma exportação.
- **Unidade de k_v.** O `rebaixamento_mc.cpp` usa kv = 1e-6 como condutividade hidráulica em m/s (`SetConductivity` espera m/s). A Tabela 2 do artigo dá k_v = 1e-6 m⁴/(kN·s), isto é k = 9,81e-6 m/s, 9,81 vezes mais permeável. No permanente isso não importa; no transiente muda o campo em t = 0,1 h. **Decidir antes do passo 2.**
- **Permanente × artigo.** O domínio e o contorno do MEF do artigo não são dados. Na nossa malha, as duas escolhas de lateral esquerda dão 1,414 e 1,429, ou seja, +6 a +7% sobre 1,336. Isso fica como informativo.

### Validação de P_u contra o artigo

- O teste hidrostático (u = −γw·y é o mesmo que um talude seco com γ = 20) trava o sinal e a escala de P_u.
- O campo semianalítico v'_opt (Seção 4.1) é totalmente especificado no artigo, sem ambiguidade de domínio. Por isso ele é a validação de P_u contra o artigo:
  - Fig. 9 (α = 1): 1,905 contra cerca de 1,95 lido em β = 45° (−2,3%); 60° e 75° a −5% e −2%.
  - Fig. 8: os quatro H_crit ficam a menos de 1,5% da leitura.
  - A Fig. 8 tem os **painéis trocados** em relação à Tabela 1. Sem percolação, os parâmetros de "London" (c = 6, φ = 32°) dão H_crit = 13,0 a 60°, e os de "Israeli" (c = 11,7, φ = 24,7°) dão 18,0; a figura mostra 18 no painel London e 13 no Israeli. Com a troca, tudo confere.
- A Eq. 31 impressa tem erros: o termo em e_θ traz (r/R_w)^(C/D−1)·C/D e o fator na frente é 1/(C − D). Com isso, div v ≠ 0 e a água entraria pela face exposta. O `percolacao_analitica.py` deduz o campo da classe (29) e do funcional (23); o primeiro termo da Eq. 40 confere com essa dedução (ver o cabeçalho do arquivo).
- A curva cheia do artigo (−∇u'_FE, 1,336) não foi reproduzida. Ver "Permanente × artigo".

### GIM não convergido

- Com NREF = 3, o GIM do caso determinístico dá 1,3253. Com mais níveis ele continua descendo: 1,2973 (NREF = 4) e 1,2851 (NREF = 5).
- A razão dos incrementos é cerca de 0,4. A extrapolação de Richardson dá 1,2758 (1,2788 com um nível a menos).
  - O cinemático (1,2791) fica 0,3% acima do GIM extrapolado.
  - O "acordo" antigo do GIM com 1,336 (1,3339 com 3 níveis) foi coincidência de malha grossa, num problema que nem é o do artigo.
- Realização 0: o GIM vai 1,0781 → 1,0588 → 1,0516, e o extrapolado é 1,047. O cinemático (1,0806) fica 3,2% acima.
  - É o esperado de um limite superior com um bloco rígido único num campo heterogêneo.
  - Os testes exigem cinemático/GIM extrapolado entre 0,99 e 1,10.

## Monte Carlo rápido (Seção 6, configuração de referência)

| arquivo | o que faz |
|---|---|
| `campos.py` | KL analítico do núcleo exp(−\|dx\|/Lx − \|dy\|/Ly) (o mesmo do `KLSeparavel.h`), M = 2000 termos por campo, lognormal; `Problema`: malha do Darcy, Darcy permanente com k heterogêneo, grade dos campos |
| `rapido.py` | `CampoGrade` (P1 numa grade uniforme, localização por índice), marcha da espiral vetorizada, `minimizar_rapido` |
| `mc.py` | driver com multiprocessing, um CSV com uma linha por realização, retomada, resumo (P_f, IC 95%, CoV) |
| `conferencia.py` | `campos` (médias, CoV, correlações), `det`, `busca` (rápido × `cinematico.minimizar` nas mesmas pressões e campos) |

```bash
python3 -I mc.py n=200 primeira=0 semente=2025 procs=4 saida=resultados/mc_piloto.csv
python3 -I mc.py saida=resultados/mc_piloto.csv resumo=1     # só o resumo
python3 -I conferencia.py busca n=30 procs=4
```

- **Domínio e contorno:** 14H × 2,8H × 8,4H (Le = 70, Ld = 14, D = 42 m) com o contorno do Ceron.
  - Malha do Darcy: grade tensorial graduada (0,0625 m no talude, até 0,25 m na janela, até 2 m fora), com 28 902 nós.
  - k é constante por célula: o valor do campo no centro da célula, o mesmo nos dois triângulos.
  - A fatoração é nova a cada realização (SuperLU, MMD_AT_PLUS_A), cerca de 0,09 s.
- **Determinístico:** nessa malha Γ = 1,3471, e ele não muda com o refino (1,3477, 1,3471, 1,3469).
  - A malha mapeada do `estudo_contorno.py` dá 1,3373 com h = 0,0625, mas ainda sobe com o refino: 1,3392 e 1,3426.
  - A energia discreta da malha tensorial é menor, então ela é a mais precisa. O 1,3373 era malha grossa.
- **Campos:** retângulo 89 × 47 m, 62 modos em x e 295 em y. A variância retida é 94,22%; o artigo diz erro < 6% com M = 2000.
  - Com a truncagem, a correlação a 2 m na horizontal é 0,937 (núcleo: 0,905) e a 0,5 m na vertical é 0,812 (núcleo: 0,779).
  - O CoV cai na mesma proporção: c 0,291 em vez de 0,30, phi 0,097 em vez de 0,10, k 0,58 em vez de 0,60.
- **Grade dos campos:** h = 0,05 m na janela x ∈ [−25, 19], y ∈ [0, 15].
  - Em 10 realizações, h = 0,05 × 0,025 difere no máximo 3e-4 em Γ; h = 0,1 × 0,025, até 1e-3.
  - Malha do Darcy 0,25 × 0,125 na janela: até 2e-3 (média +5e-4).
- **Rápido × exato** (30 realizações, mesmos campos): diferença máxima +7e-5 e média +2e-5. Nenhuma passa de 0,1%.
  - Tempo: 0,66 s por realização (cerca de 240 avaliações de 2,8 ms), contra 13,9 s (830 avaliações de 17 ms).

## Pontos em aberto

- O Monte Carlo (passo 2) deve guardar `n_normalidade` e `n_phi_nao_conv` de cada realização. Nos testes são zero, inclusive em 400 mecanismos sorteados na realização 0 (86 039 arcos conferidos com 200 pontos cada).
- A unidade de k_v no transiente (ver acima).
- O que é Ω_i (as células do campo) o artigo não define. Com a regra de φ_i acima, isso deixou de importar: só conta φ(x) no arco.
- O p_w exportado é só o valor nos vértices. O transiente é de ordem alta (59 548 equações para 6021 nós), mas aqui ele é interpolado linearmente.
- A forma fechada do mecanismo (I) continua recusando C atrás da reta da face, porque a f3 (Eq. 53) não foi conferida nesse caso. Em 120 000 sorteios, nenhum mecanismo (I) válido caiu ali. O código geral não tem esse limite.
- A grade de partida da minimização usa mecanismos clássicos com o φ médio. Com β = φ (por exemplo, β = 30° e φ = 30°, a ponta esquerda da Fig. 9), a grade fica vazia e `minimizar` falha.
- v'_opt está implementado só para α = 1. Para α > 1, a energia da zona 1 tem um termo cruzado que as Eqs. 36-37 não mostram.

## Piloto de 200 e retomada segura (2026-11-20)

* `mc.py` grava em cada linha `semente`, `h`, `hD`, `M` e, na retomada, aborta se o arquivo for de outra semente ou
  outra discretizacao; corta uma linha incompleta no fim do CSV antes de continuar; avisa quando uma realizacao fica
  sem mecanismo (Gamma = NaN, fora de S). Grava tambem `max_varr` e `n_max_regra` (marchas que passaram de 15
  varreduras e usaram a regra "so sobe").
* `dados/partida_det_*.json` guarda a assinatura do caso (h, hD, M, hw, gamma', estatisticas, dominio, janela) e e
  recalculado se ela nao bater.
* `rapido.minimizar_rapido`: o Gamma final de cada classe e de uma marcha FRIA (sem o cache de phi_i), com a
  normalidade conferida em todos os arcos e P_u na forma de volume.
* `figura_mc.py csv=... nome=... [parede=s] [conf=...]`: figura `figs/<nome>.png/.pdf` e `dados/<nome>_resumo.txt`.
* Piloto: `/home/diogo/projects/runs/ceron_mc_python/piloto200` (semente 2025, realizacoes 0 a 199), resumo em
  `dados/piloto200_resumo.txt`.
