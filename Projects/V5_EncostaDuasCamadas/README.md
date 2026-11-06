# V5_EncostaDuasCamadas: encosta em duas camadas sob chuva e fator de segurança local (LFS)

Referência: R. Abbasov, M. Fahs, A. Younes, H. Nowamooz, K. J. Måløy e R. Toussaint (2024), *Modeling
rainfall-induced landslide using the concept of local factor of safety: Uncertainty propagation and sensitivity
analysis*, Computers and Geotechnics 167:106102, doi:10.1016/j.compgeo.2024.106102 (secs. 2 a 4, Tabelas 1, 3 e 4,
Figs. 2, 4 a 6, 13 e 14). É o benchmark de encosta de Lu et al. (2012) e Shao et al. (2015).

## O caso

Seção de encosta 2D (Fig. 2): domínio de 42 m × 25 m, crista em z = 25 m, talude de 6 m de altura e pé em z = 19 m.
Uma camada superior de 2 m acompanha o terreno sobre a camada inferior. As duas camadas têm a mesma mecânica e
hidráulicas muito diferentes: K_s = 2,59 e 0,06 cm/h (25,9 e 0,6 mm/h).

A chuva tem duas etapas. Primeiro, 10 anos a 600 mm/ano a partir do estado hidrostático com o nível d'água em
z = 17 m; o fim dessa etapa é o "estado de equilíbrio", o instante 0. Depois, 20 mm/h por 15 h. A chuva de 20 mm/h é
menor que o K_s da camada superior e muito maior que o da inferior, então a água se acumula sobre a interface.

O fluxo é o de Richards com Brooks-Corey (eqs. 1 a 4). O campo de tensões é linear elástico, com o peso da mistura que
depende da umidade, ∇·σ + ((1 − θ_s) ρ_s + θ ρ_w) g = 0 (eq. 5). A tensão efetiva é σ' = σ − χ ρ_w g h I com
χ = S_e (eq. 7). O fator de segurança local de Lu et al. (2012) é LFS = cos φ' (c' + σ'_I tan φ') / σ_II (eq. 8),
com o centro σ'_I e o raio σ_II do círculo de Mohr. Na zona de investigação (22 m de largura, 8 m de profundidade sob
a crista, 113 m²) o artigo calcula o LFS médio (eq. 11) e a área com LFS < 1 (eqs. 12 e 13).

O que é comparado:

| figura | caso | grandeza |
|---|---|---|
| Fig. 4 | Tabela 1 | mapas de θ a 0, 7,5 e 15 h |
| Fig. 5 | Tabela 1 | mapas da tensão de sucção χ ρ_w g \|h\| |
| Fig. 6 | Tabela 1 | mapas do LFS |
| Fig. 14 e sec. 4.6 | Tabela 4 (determinístico) | LFS médio e A(LFS < 1) a 0, 7,5 e 15 h |
| Fig. 13 | média das faixas da Tabela 3 | efeitos marginais de ν, φ' e c' |

Os mapas e as curvas foram digitalizados em `dados/` (ver "Dados de referência").

Os exemplos de incerteza do artigo não foram feitos: a triagem de Plackett-Burman, o PCE com índices de Sobol, os
500 casos com ±10 % e o caso arenoso com o ramo de umedecimento.

## Dados

| dado | valor usado | fonte / suposição |
|---|---|---|
| geometria | crista z = 25 m até x = 19 m; pé (x = 34 m, z = 19 m); base z = 0; domínio 0 ≤ x ≤ 42 m | Fig. 2, posições medidas no desenho (25 px/m). O rótulo "23°" não fecha com os 22 m nem com os 113 m² da zona: usamos 6 m em 15 m (21,8°), que dá exatamente 4 + 15 + 3 = 22 m e 32 + 75 + 6 = 113 m² |
| camada superior | 2 m de espessura **vertical**; interface (0, 23) – (19, 23) – (34, 17) – (42, 17) | Figs. 2, 4 e 6: a interface é o terreno 2 m abaixo |
| zona de investigação | 15 ≤ x ≤ 37 m, z ≥ 17 m | Fig. 2 (medida) e os 113 m² da sec. 4.6 |
| θ_s, θ_r | 0,412 / 0,041 (sup.); 0,385 / 0,09 (inf.); porosidade = θ_s | Tabela 1 |
| K_s | 2,59 cm/h (sup.); 0,06 cm/h (inf., Tabela 1) ou 0,05 cm/h (Tabela 4) | Tabelas 1 e 4; isotrópica |
| Brooks-Corey | α_BC = 0,068 e 0,027 1/cm; n_BC = 0,322 e 0,131; l_BC = 1 nas duas | Tabela 1. O artigo imprime um só l_BC; supusemos o mesmo para a camada inferior |
| retenção usada | van Genuchten a = 0,3350 e 0,05630 1/kPa, b = 1,4762 e 1,5664, S_wr = θ_r/θ_s | ajuste a Brooks-Corey em 0 < h ≤ 1,5 m (`scripts/ajuste_vg.py`, `figuras/ajuste_vg.png`); ver "O nosso modelo" |
| condutividade | k_rw = S_e^p com p = l + 2 + 2/n_BC = 9,211 e 18,267 | eq. (4), exata (a lei de potência da biblioteca) |
| E, ν, φ' | 10 MPa; ν = 0,35 (Tabela 1) ou 0,3 (Tabela 4); φ' = 35° | Tabelas 1 e 4. E não muda as tensões (só forças de volume e apoios) |
| c' | 3 kPa (sup.) e 6 kPa (inf.) no caso `tabela1`; 4,5 kPa nas duas no `tabela4`; 5 kPa no `fig9` | A Tabela 1 dá 3 e 6 nas colunas das camadas, embora a sec. 3.1 diga que a mecânica é igual nas duas. Usamos a tabela como está |
| ρ_s, ρ_w, g | 2636 kg/m³; 1000 kg/m³; 9,81 m/s² | Tabela 1; g suposto |
| chuva | 600 mm/ano por 10 anos (365,25 d/ano); 20 mm/h por 15 h; fluxo vertical por área horizontal | sec. 2.3: "o fluxo normal no talude é o horizontal vezes o cosseno da inclinação" |
| condição inicial | hidrostática com nível em z = 17 m, sem limite de sucção | sec. 2.3 |
| lateral direita | carga hidrostática com nível 17 m abaixo de z = 17 m; face de percolação de 17 a 19 m | Fig. 2. O nível da carga e a cota da troca não são dados: supusemos 17 m, o nível inicial e a base da camada superior no desenho |
| superfície | chuva e face de percolação (p_w ≤ 0) quando satura | sec. 2.3. O artigo descreve a percolação como "saída crescendo com a saturação" e não dá a fórmula; usamos a condição de Signorini (p_w = 0 com saída, ou fluxo dado com p_w < 0) |
| base e lateral esquerda | impermeáveis; roletes; roletes também na lateral direita | Fig. 2 |
| armazenamento saturado | compressibilidade da água K_w = 2,2·10⁶ kPa | suposto (o modelo de armazenamento do COMSOL não é dado) |
| estado plano | deformação plana | suposto. O artigo não diz, e a sua eq. (A5), LFS_min = (1 + ν)/(1 − ν) sin φ', corresponde a K0 = ν (tensão plana). Os nossos resultados com deformação plana batem com as Figs. 6, 13 e 14; com tensão plana (ν* = ν/(1 + ν) = 0,23 para ν = 0,3) a área com LFS < 1 a 0 h seria ~32 m² (interpolada na nossa sensibilidade a ν) contra 14,8 m² |

### Dados de referência (`dados/`), digitalizados, aproximados

- `fig4_theta.csv`, `fig5_sigs.csv`, `fig6_lfs.csv`, `fig9_lfs.csv`: os mapas das Figs. 4 a 6 e do painel COMSOL
  da Fig. 9 (15 h) numa grade de 0,25 m da zona
  (x_zona = x − 15 m, z_zona = z − 17 m), nos três instantes (a Fig. 9 só a 15 h). `perfis_fig4_5_6.csv` traz os
  perfis verticais em
  x_zona = 2, 10, 17 e 20,5 m. Cada pixel vira o valor da cor mais próxima da barra de cores
  (`scripts/digitaliza_mapas.py <pdf>`, 250 dpi, ~0,02 m por pixel). As barras saturam nos extremos: θ ≤ 0,22 ou
  ≥ 0,40, tensão de sucção ≥ 8 kPa e LFS ≤ 0,5 ou ≥ 3 ficam no extremo.
- `fig14_lfs.csv`: as barras da Fig. 14 (face frontal das barras 3D, ±0,005 no LFS, ±0,05 m² na área).
- `fig13_efeitos_marginais.csv`: as curvas da Fig. 13 (±0,02 no LFS, ±1 m² na área). O instante não é dito; os
  valores na média das faixas batem com a Fig. 14 a 15 h.
- `resultados_malha_fina/`: metricas, perfis, sensibilidade e resumo das três execuções da malha fina (`tabela1`,
  `tabela4`, `fig9`) usadas nas tabelas abaixo. `comparacao_malha_fina.csv` é a tabela "artigo × nosso" completa.

## O nosso modelo

Material `TPZMatUnsaturatedUPwPa` com uma camada por grupo físico (`UnsatMesh` do `UnsatModel.h`, um material por id)
e a análise `TPZUnsatAnalysis` (Euler implícito, Newton com busca linear, PARDISO).

- **Fluxo de Richards num domínio fixo.** p_a = 0 é imposto em todos os nós e todos os deslocamentos ficam presos no
  fluxo (esqueleto rígido, `SetHydroMechanicalCoupling(false)`). Só a equação da água é resolvida (P1).
- **Tensões totais elásticas (eq. 5).** O coeficiente de Biot é α = 0, então as poropressões não entram no equilíbrio
  e o tensor da memória é a tensão total σ = D ε. O peso é o da mistura, (1 − n) γ_s + n S_w γ_w com n = θ_s, que é
  exatamente (1 − θ_s) ρ_s g + θ ρ_w g. A coesão do Mohr-Coulomb é posta em 10⁹ kPa, então o cálculo fica linear
  elástico.
- **Cálculo de tensões em cada instante pedido** (`Encosta::Tensoes`): parte do estado sem tensões, com as pressões
  do instante congeladas, e resolve só o equilíbrio. O estado do fluxo (solução e memória) é restaurado depois.
- **LFS** (`LFSdeTensao`): com as tensões totais em cada ponto de integração (3 por triângulo) e χ = S_e da camada
  (`EChiEffective`), σ'_I = (σ_1 + σ_3)/2 − S_e p_w (compressão positiva) e σ_II = (σ_1 − σ_3)/2 no plano. c' e φ'
  da camada entram só aqui.
- **Métricas** (eqs. 11 a 13): somas nos pontos de integração dos elementos da zona. A malha tem as arestas da zona e
  da interface, então a área da zona dá 113,000 m² exatos. Além da média literal da eq. (11), `LFS_medio_lim3` é a
  média com o LFS limitado a 3, o topo da escala de cores das Figs. 6 e 9. Ver a discussão.
- **Retenção.** Brooks-Corey (eq. 3) não existe na biblioteca. Ajustamos van Genuchten, S_e = [1 + (a p_c)^b]^(1/b − 1),
  por mínimos quadrados de S_e em 0 < h ≤ 1,5 m: a faixa de sucção do problema depois dos 10 anos, que com
  Brooks-Corey é ~1,1 m na camada superior e ~0,9 m na inferior. RMSE de S_e: 0,023 e 0,011. A condutividade é a de
  Brooks-Corey exatamente (lei de potência em S_e), então o θ dos trechos de gradiente unitário é exato: 0,2358 e
  0,3519 com 600 mm/ano, 0,4017 com 20 mm/h. A diferença está na sucção a S_e dado:

  | estado | h Brooks-Corey | h van Genuchten |
  |---|---|---|
  | 10 anos, camada superior (S_e = 0,525) | 1,088 m | 1,067 m |
  | 10 anos, camada inferior (S_e = 0,888) | 0,918 m | 0,991 m |
  | chuva, camada superior (S_e = 0,972) | 0,160 m (fim da entrada de ar) | 0,060 m |

  O van Genuchten não tem o patamar de entrada de ar, então perto da saturação a tensão de sucção na zona umedecida
  fica 0,57 kPa contra 1,53 kPa. O efeito no LFS é medido na discussão.
- **Contornos.** A superfície recebe a chuva pela projeção horizontal (`BoundaryNodes(..., true)`). Os nós da
  superfície e da lateral direita acima de z = 17 m são candidatos à face de percolação, com a lógica do `passo()` do
  `SlopeRandomField`:
  - um nó livre com p_w > 10⁻⁶ kPa passa a p_w = 0;
  - um nó com p_w = 0 que exigiria entrada maior que a chuva que recebe volta ao fluxo dado;
  - o passo só é aceito quando o conjunto não muda.

  Na lateral direita, abaixo de z = 17 m, p_w = γ_w (17 − z).
- **Passos.**
  - Etapa 1: dt de 1 h, vezes 1,5 enquanto o Newton leva até 8 iterações, até 30 dias; o passo que não converge é
    dividido ao meio. São 142 passos sem nenhum corte. O permanente é alcançado em ~3 anos: no fim, a saída pela
    lateral direita iguala a chuva, 7,985·10⁻⁷ m³/s por metro, com max |Δp_w| ~10⁻¹⁰ kPa por passo de 14 a 30 dias.
  - Etapa 2: dt = 0,1 h (150 passos). Os 15 h exatos e os instantes de saída são atingidos.
  - O balanço de água de cada passo está em `passos.csv`, com resíduo ~10⁻¹⁵ m³.
- **Malha** (`malhas/encosta.geo`, Gmsh; `scripts/gera_malhas.sh`). A camada superior é estruturada (blocos
  transfinitos): 16 linhas de 0,125 m e 0,25 m na horizontal na zona. A camada inferior não é estruturada: 0,3 m na
  zona, até 3 m na base. São 4320 nós, 8441 triângulos e 42 800 equações (u P2, p_w e p_a P1), parecido com os
  "~3 K nós" do artigo. A `encosta_grossa.msh` (969 nós) é a dos testes e dá as métricas a menos de 0,5 % da fina
  (0 h: A = 14,63 × 14,69 m², LFS médio limitado 1,732 × 1,731). A 15 h a diferença é de 2,5 % na área (22,20 ×
  22,76 m², com dt = 0,25 h na grossa) e de 0,3 % no LFS médio. A `encosta_muito_fina.msh` (10 037 nós) é para o
  estudo de malha no cluster.

## Compilar e rodar

```bash
cmake --build <build> --target V5_EncostaDuasCamadas V5_testes
cd <pasta de saída>
CASO=tabela1 <build>/Projects/V5_EncostaDuasCamadas/V5_EncostaDuasCamadas      # mapas das Figs. 4 a 6
CASO=tabela4 NU_LISTA=0.2,0.25,0.3,0.35,0.4 C_LISTA=3,3.5,4,4.5,5,5.5,6 PHI_LISTA=25,30,35,40,45 \
    <build>/Projects/V5_EncostaDuasCamadas/V5_EncostaDuasCamadas                 # Figs. 14 e 13
python3 <projeto>/scripts/compara_v5.py --tab1 <pasta tabela1> --tab4 <pasta tabela4>   # tabela e figuras
```

Cada caso leva ~20 min na malha fina, num núcleo de 2,1 GHz (etapa 1 ~7 min, etapa 2 ~12 min). Variáveis de
ambiente (padrões entre parênteses):

| variável | o que faz |
|---|---|
| `CASO` | `tabela1` (padrão), `tabela4` ou `fig9` |
| `MALHA` | `encosta.msh`; um nome sem `/` é procurado em `malhas/` |
| `NU`, `PHI_DEG`, `E` | ν, φ' (graus) e E (kPa) |
| `C_SUP`, `C_INF`, `COES` | c' de cada camada, ou `COES` para as duas (kPa) |
| `KS_SUP`, `KS_INF` | K_s de cada camada (cm/h) |
| `A_VG_SUP`, `B_VG_SUP`, `A_VG_INF`, `B_VG_INF` | o ajuste de van Genuchten (a em 1/kPa) |
| `RHO_S` | ρ_s |
| `ANOS`, `CHUVA1_MM_ANO`, `CHUVA2_MM_H`, `HORAS` | as duas etapas de chuva |
| `NIVEL`, `Z_FACE` | nível inicial e da carga na lateral direita, e a cota da troca para a face de percolação (17 e 17 m) |
| `DT0_H`, `FATOR_DT`, `DT_MAX_D` | passos da etapa 1 (1 h, 1,5 e 30 d) |
| `DT2_H` | passo da etapa 2 (0,1 h) |
| `MAX_CORTES` | divisões de um passo que não converge (12) |
| `SAIDAS_H` | instantes da etapa 2 (0,7.5,15) |
| `TOLU`, `TOLP`, `MAXIT` | tolerâncias e iterações do Newton |
| `ESTADO_INICIAL`, `GRAVA_ESTADO` | lê ou grava o estado do fim da etapa 1 (p_w nos nós), por exemplo para pular a etapa 1 |
| `NU_LISTA`, `C_LISTA`, `PHI_LISTA` | sensibilidade das métricas em cada instante (`sensibilidade.csv`): cada ν refaz as tensões; c' (o mesmo nas duas camadas) e φ' só refazem o LFS |
| `VTK`, `PONTOS_CSV`, `CAMPOS_CSV` | saídas opcionais (0, 0, 1) |
| `VERBOSO` | 0, 1 ou 2 |
| `SOLVER` | `pardiso` (padrão) ou `skyline` |

Saídas na pasta corrente:

- `passos.csv`: cada passo aceito, com o balanço de água;
- `metricas.csv`: LFS médio (literal e limitado), A(LFS < 1) total e por camada, faixas de θ, nível d'água e frente
  de umedecimento de cada instante;
- `perfis.csv`: p_w, θ, S_e, tensão de sucção e LFS nas verticais x = 17, 25, 32 e 35,5 m;
- `campos_t<h>.csv`, `nos.csv` e `elementos.csv`: os mapas (`scripts/mapas_v5.py <pasta>`);
- `pontos_t<h>.csv`, `encosta_t<h>.vtk` e `sensibilidade.csv`: as saídas opcionais;
- `estado_etapa1.csv`: o estado do fim da etapa 1;
- `resumo.csv`: uma linha;
- `resultado.txt`: uma linha chave=valor com `LFS_medio_<t>h`, `LFS_medio_lim3_<t>h`, `A_LFS_<t>h`, `situacao=ok` e
  `seg`; é gravado por último, de um `.tmp` renomeado. Código de saída 0 em sucesso.

**No cluster (fermi, SLURM)**, de dentro de `cluster/`:

```bash
source ../../Cluster/modulos.sh && ./compilar_v5.sh           # no nó de login
sbatch v5.sbatch                                              # as 10 execuções da lista
EXECUCOES="tabela1 tabela4" TESTES=1 sbatch v5.sbatch         # só os dois casos e os testes
```

O `v5.sbatch` roda as execuções nomeadas em paralelo, um núcleo cada, em `$SAIDA/<nome>` (padrão `$DATA/v5` ou
`./resultados_v5`):

- `tabela1`, `tabela4` (com as listas de sensibilidade) e `fig9`;
- `tabela4_ks001` a `tabela4_ks01`: o efeito de K_s(l) da Fig. 13;
- `tabela4_muito_fina`: o estudo de malha;
- `tabela4_dt`: o estudo de passo (dt de 0,05 h e no máximo 10 d);
- `tabela4_plano_tensao`: tensão plana via ν* = ν/(1 + ν).

No fim, o job junta os `resumo.csv`, roda `compara_v5.py` e, com `TESTES=1`, os testes padrão e os `[.lento]`.

## Testes

`V5_testes` (Catch2, marca `[V5]`; `ctest` roda os padrão; os demorados: `V5_testes "[.lento]"`). Resultados nesta
máquina:

| teste | o que compara | tolerância e justificativa | resultado |
|---|---|---|---|
| ajuste de van Genuchten a Brooks-Corey | RMSE de S_e em 0 < h ≤ 1,5 m; expoentes de K; θ do permanente de 600 mm/ano com o material da biblioteca (= Brooks-Corey); sucção desse estado | RMSE < 0,0252 e 0,0120 (o do ajuste mais 10 %); θ a 2·10⁻⁴; h a 10 % (2 % e 8 % obtidos) | passou, 0,000 s |
| coluna confinada lateralmente (analítica) | coluna 1 × 10 m, duas camadas, nível em 5 m, ν = 0,25: σ_v = ∫ γ(θ) dz, σ_h = K0 σ_v, LFS em forma fechada, A(LFS < 1) | erros obtidos 3·10⁻⁵ (σ_v), 1·10⁻⁶ (K0), 1·10⁻⁵ (LFS); limites 2·10⁻⁴, 10⁻⁵, 10⁻⁴; área igual (3 m²) | passou, 0,26 s |
| fluxo de 10 anos e chuva de 15 h (tabela1, malha grossa) | permanente (saída = chuva); θ e tensão de sucção na crista; nível d'água; frente a 7,5 e 15 h; saturação no pé; RMS dos mapas das Figs. 4 e 5 | ver os comentários do `testes.cpp`: θ ±0,01 (a cor da figura), sucção ±0,6 kPa, nível ±0,6 m (carga da lateral suposta, sem franja capilar), frente ±0,2/0,25 m (malha grossa), RMS de θ < 0,03, RMS da sucção < 2 kPa. Obtidos: RMS de θ 0,004 / 0,011 / 0,007, da sucção 0,47 / 0,67 / 0,63 kPa e do LFS limitado (Fig. 6) 0,10 / 0,11 / 0,16 a 0 / 7,5 / 15 h | passou, 94 s |
| LFS médio e A(LFS < 1) contra a Fig. 14 (tabela4, malha grossa) | área da zona; LFS médio limitado a 3 a ±6 %; A(LFS < 1) a 0 e 7,5 h a ±25 %; a 15 h só entre 0,5 e 2 vezes; queda do LFS médio com a chuva | ver a discussão: a média literal fica ~30 % acima e a área a 15 h ~50 % acima | passou, 97 s |
| caso completo na malha fina (tabela4 e tabela1) `[.lento]` | o mesmo com os passos do `main.cpp`: LFS médio limitado a ±5 %, área a 0 e 7,5 h a ±25 %, θ, nível e frentes (±0,15 e ±0,2 m), RMS dos mapas das Figs. 4 a 6 (LFS limitado < 0,5) | obtidos (tabela1): RMS de θ 0,004 / 0,008 / 0,006, da sucção 0,46 / 0,60 / 0,64 kPa, do LFS 0,083 / 0,079 / 0,156 | passou, 2241 s (as duas execuções de ~19 min) |

O `ctest` padrão (os quatro primeiros) leva 3 min 13 s; `V5_testes -d yes` mostra o tempo de cada teste e os
valores comparados. O teste da Fig. 14 só exige a ordem de grandeza da área a
15 h: essa grandeza **não** concorda com o artigo (ver a discussão), e o teste não a esconde.

## Artigo × nosso resultado (malha fina)

`python3 scripts/compara_v5.py` com as execuções `tabela1` e `tabela4` da malha fina (`dados/comparacao_malha_fina.csv`
tem a tabela inteira). As profundidades são verticais, a partir da superfície. Crista = x_zona = 2 m;
meio do talude = x_zona = 10 m.

| grandeza | artigo | nosso | diferença | fonte |
|---|---|---|---|---|
| θ camada superior, crista, 1 m, 0 h | 0,236 | 0,236 | 0,000 | Fig. 4 |
| θ topo da camada inferior, crista, 2,1 m, 0 h | 0,352 | 0,353 | +0,001 | Fig. 4 |
| θ camada inferior, crista, 3 m, 0 h | 0,357 | 0,361 | +0,004 | Fig. 4 |
| tensão de sucção, camada superior, crista, 1 m, 0 h (kPa) | 5,54 | 5,49 | −0,05 | Fig. 5 |
| tensão de sucção, topo da camada inferior, crista, 2,1 m, 0 h (kPa) | 7,95 | 8,47 | +0,52 | Fig. 5 |
| profundidade do nível d'água na crista, 0 h (m) | 3,65 a 4,0 | 4,05 | +0,05 a +0,4 | Fig. 5 (a tensão de sucção cai de 3,5 kPa a 0 entre 3,63 e 3,68 m) |
| frente de umedecimento (θ = 0,32), crista, 7,5 h / 15 h (m) | 0,85 / 1,70 | 0,91 / 1,81 | +0,06 / +0,12 | Fig. 4 |
| frente de umedecimento, meio do talude, 7,5 h / 15 h (m) | 0,91 / 1,87 | 0,92 / 1,96 | +0,01 / +0,09 | Fig. 4 |
| tensão de sucção na zona umedecida, crista, 0,5 m, 15 h (kPa) | 1,52 | 0,57 | −0,95 | Fig. 5 |
| LFS crista, 7,6 m / 4 m, 0 h | 1,122 / 1,912 | 1,124 / 1,931 | +0,002 / +0,018 | Fig. 6 |
| LFS crista, 7,6 m / 4 m, 15 h | 1,126 / 1,881 | 1,116 / 1,901 | −0,010 / +0,020 | Fig. 6 |
| LFS meio do talude, 1 m / 4 m, 0 h | 1,885 / 1,603 | 1,913 / 1,621 | +0,03 / +0,02 | Fig. 6 |
| LFS meio do talude, 1 m / 4 m, 15 h | 1,724 / 1,597 | 1,675 / 1,559 | −0,05 / −0,04 | Fig. 6 |
| A(LFS < 1) caso Tabela 1, 0 / 7,5 / 15 h (m²) | "só na superfície" / "cresce sob o talude" / "vai para o pé" | 0,01 / 3,0 / 4,7 (toda na camada superior) | qualitativo | Fig. 6 |
| mapa de θ, RMS em ~1670 pontos, 0 / 7,5 / 15 h | – | 0,004 / 0,008 / 0,006 (média +0,001 a +0,003) | – | Fig. 4 |
| mapa da tensão de sucção, RMS, 0 / 7,5 / 15 h (kPa) | – | 0,47 / 0,60 / 0,65 | – | Fig. 5 |
| mapa do LFS limitado a 3, RMS, 0 / 7,5 / 15 h | – | 0,074 / 0,079 / 0,156 (média +0,016 / +0,001 / −0,058) | – | Fig. 6 |
| caso fig9 (c' = 5 kPa): LFS crista 1 / 4 / 7,5 m, 15 h | 2,98 / 1,84 / 1,12 | 3,00 / 1,86 / 1,11 | +0,02 / +0,02 / −0,02 | Fig. 9 (COMSOL) |
| caso fig9: LFS meio do talude 1 / 4 / 5 m, 15 h | 1,87 / 1,57 / 1,41 | 1,84 / 1,52 / 1,33 | −0,04 / −0,05 / −0,08 | Fig. 9 |
| caso fig9: mapa do LFS limitado a 3, RMS, 15 h | – | 0,152 (média −0,062) | – | Fig. 9 |
| LFS médio (eq. 11), 0 / 7,5 / 15 h | 1,684 / 1,626 / 1,512 | 2,244 / 2,130 / 1,866 | +33 / +31 / +23 % | Fig. 14 |
| LFS médio limitado a 3, 0 / 7,5 / 15 h | 1,684 / 1,626 / 1,512 | 1,731 / 1,675 / 1,527 | +2,8 / +3,0 / +1,0 % | Fig. 14 |
| A(LFS < 1), 0 / 7,5 / 15 h (m²) | 14,77 / 15,19 / 14,9 | 14,69 / 15,83 / 22,76 | −0,5 / +4,2 / **+53 %** | Fig. 14, sec. 4.6 |
| área da zona (m²) | ~113 | 113,000 | 0 | sec. 4.6 |
| Fig. 13, LFS médio limitado, ν = 0,2 / 0,25 / 0,3 / 0,35 / 0,4 (15 h) | 0,99 / 1,27 / 1,55 / 1,77 / 1,95 | 1,18 / 1,33 / 1,53 / 1,64 / 1,81 | +0,19 a −0,14 | Fig. 13 |
| Fig. 13, LFS médio limitado, φ' = 25 / 30 / 35 / 40 / 45° | 1,22 / 1,41 / 1,53 / 1,67 / 1,76 | 1,28 / 1,41 / 1,53 / 1,63 / 1,73 | ≤ 5 % | Fig. 13 |
| Fig. 13, LFS médio limitado, c' = 3 / 4 / 5 / 6 kPa | 1,44 / 1,51 / 1,57 / 1,635 | 1,45 / 1,50 / 1,55 / 1,60 | ≤ 2,1 % | Fig. 13 |
| Fig. 13, A(LFS < 1), ν = 0,2 / 0,25 / 0,3 / 0,35 / 0,4 (m², 15 h) | 63,8 / 33,6 / 12,4 / 1,8 / 0,4 | 59,7 / 44,6 / 22,8 / 2,7 / 3,8 | mesma ordem; ν domina | Fig. 13 |
| Fig. 13, A(LFS < 1), φ' = 25 / 35 / 45° (m², 15 h) | 36,8 / 12,4 / 1,9 | 51,2 / 22,8 / 3,0 (a 0 h: 34,8 / 14,7 / 1,8) | 1,4 a 1,8 vezes | Fig. 13 |
| Fig. 13, A(LFS < 1), c' = 3 / 6 kPa (m², 15 h) | 17,3 / 7,3 | 32,1 / 16,3 (a 0 h: 17,8 / 11,8) | ~2 vezes | Fig. 13 |

O caso `fig9` usa `python3 scripts/compara_v5.py ... --fig9 <pasta fig9>`.

Figuras: `figuras/mapas_tabela1.png`, `figuras/perfis_v5.png` (linhas nossas, pontos do artigo),
`figuras/fig14_v5.png`, `figuras/fig13_v5.png` e `figuras/ajuste_vg.png`.

## Discussão

**Fluxo.** θ e a tensão de sucção depois dos 10 anos batem com as Figs. 4 e 5 em toda a zona. Nos trechos de
gradiente unitário, θ é exato por construção, porque a condutividade é a de Brooks-Corey. A frente de umedecimento
fica 6 a 12 cm mais funda que a do artigo (~7 %). O balanço de massa dá 0,15 m / (0,402 − 0,236) = 0,90 m a 7,5 h,
que é o nosso valor; a do artigo (0,85 m) corresponde a uma zona umedecida mais perto de θ_s. O nível d'água na
crista fica de 0,05 a 0,4 m mais fundo, conforme a leitura da Fig. 5. Três coisas o afetam: a carga da lateral
direita e a cota da face (supostas), a ausência da franja capilar de Brooks-Corey (0,37 m na camada inferior) no van
Genuchten, e K_s(l), já que com 0,05 cm/h o nível sobe 0,57 m.

O desvio mais nítido é a tensão de sucção na camada superior umedecida: 0,57 kPa contra 1,52 kPa. É o efeito do
ajuste de van Genuchten perto da saturação (h = 0,06 m contra 0,16 m em S_e = 0,972). Uma curva de van Genuchten que
passa pelos dois estados de Brooks-Corey da camada superior (a = 0,1526 1/kPa, b = 1,993: mesma sucção em S_e = 0,525
e 0,972) dá 1,62 kPa. Na malha grossa, caso `tabela4`, ela muda:

| grandeza (0 / 7,5 / 15 h) | ajuste padrão | van Genuchten de dois pontos |
|---|---|---|
| A(LFS < 1) (m²) | 14,6 / 15,9 / 22,2 | 13,2 / 13,9 / 19,1 |
| LFS médio limitado | 1,732 / 1,681 / 1,531 | 1,757 / 1,701 / 1,560 |

Essa curva não substitui a padrão: ela piora o ajuste entre os dois pontos (S_e até 0,1 acima em h de 0,2 a 0,5 m).

**Tensões e LFS.** Com o caso da Fig. 6 (Tabela 1), os perfis do LFS na crista e no meio do talude batem com o mapa
do artigo a 0,02–0,05 em toda a profundidade, a 0 e a 15 h. O LFS vale 1,12 no fundo da zona, onde fica abaixo do
nível d'água, e passa de 3 logo abaixo da crista. Isso confirma a abordagem do artigo, reproduzida aqui:

- tensão total elástica com o peso da umidade do instante, deformação plana, χ = S_e;
- a tensão de sucção que endurece a camada superior;
- o "afrouxamento" da parte saturada, onde o LFS fica perto de 1 abaixo da crista porque a poropressão positiva se
  soma a um estado K0 de tensão total.

O teste da coluna confere o pós-processamento com a solução analítica a 10⁻⁵. A dependência forte de ν (Fig. 13) é
reproduzida. Ela vem só da razão K0 = ν/(1 − ν) no cálculo elástico: o SRM do projeto, por exemplo, quase não depende
de ν.

**LFS médio: literal × limitado.** A média literal da eq. (11) fica 23 a 33 % acima da Fig. 14. O LFS ~ 1/σ_II vai a
infinito onde o círculo de Mohr se reduz a um ponto. Isso acontece numa faixa a ~1,2 m de profundidade sob a crista,
onde a tensão horizontal passa da tração da borda da crista à compressão (LFS até ~25), e em pontos isolados no pé.
Temos ~14 m² com LFS > 3 a 0 h. Limitando o LFS a 3, o topo da escala das Figs. 6 e 9, a média cai para +1 a +3 % da
Fig. 14 e bate com as curvas da Fig. 13: ≤ 5 % em φ' e c', de −8 a +20 % em ν. O artigo não diz que limitou o LFS.
Isto é uma **hipótese** que explica a concordância; as duas médias estão nas saídas.

**A(LFS < 1).** A área bate a 0 e 7,5 h (−0,5 % e +4,2 %). A 0 h ela está toda na camada inferior, no fundo da zona
sob a crista e o talude: com ν = 0,3, a poropressão positiva sobre um estado K0 leva o LFS abaixo de 1. As curvas da
Fig. 13 a 0 h também batem (φ': 34,8 / 14,7 / 1,8 m² contra 36,8 / 12,4 / 1,9).

A 15 h, a nossa área cresce para 22,8 m², enquanto a do artigo fica em 14,9 m² (+53 %). Entre 7,5 e 15 h, a água
empoleirada sobre a interface no fim do talude se liga ao nível d'água. p_w sobe 10 a 15 kPa na camada inferior
saturada perto do pé, e o LFS cai 0,3 a 0,5 ali. O efeito existe no artigo: na Fig. 6, o LFS perto do pé cai de ~1,6
para ~1,1 entre 0 e 15 h, com ν = 0,35. Mas com os dados da Tabela 4 a área do artigo quase não muda (14,8, 15,2,
14,9 m²).

O ajuste de van Genuchten explica ~1/5 da diferença (−3,1 m² com a curva de dois pontos). O resto depende de como a
água empoleirada drena pelo pé e pela lateral direita: da cota e da forma da face de percolação ("saída crescendo
com a saturação", sem fórmula no artigo, contra a nossa condição de Signorini) e da carga suposta na lateral direita.
Não ajustamos esses dados para fechar a diferença.

**Outras escolhas.**

- Talude de 21,8° (não os 23° do rótulo).
- c' = 3 e 6 kPa por camada no caso `tabela1`, como na Tabela 1. A 0 h a área com LFS < 1 é ~0, e a média depende
  pouco de c' (Fig. 13).
- E sem efeito nas tensões.
- Deformação plana confirmada pelos resultados. A eq. (A5) do artigo é a de K0 = ν; com deformação plana, o LFS a
  grande profundidade num terreno plano seco é sin φ'/(1 − 2ν).

## Pendências

- **Brooks-Corey na biblioteca** (`TPZUnsaturatedHydraulics`): eliminaria o ajuste de van Genuchten e a diferença da
  tensão de sucção perto da saturação. Ver o resultado da tarefa.
- O efeito de K_s(l) da Fig. 13, o estudo de malha (`encosta_muito_fina.msh`), o de passo e o de tensão plana estão
  no `v5.sbatch` e não foram rodados aqui. O caso `fig9` foi rodado a partir do estado de 10 anos do `tabela1`
  (`ESTADO_INICIAL`, o mesmo K_s(l)); a 15 h, A(LFS < 1) = 1,03 m² e o LFS médio limitado vale 1,740.
- A diferença de A(LFS < 1) a 15 h no caso Tabela 4: falta testar a forma da face de percolação do artigo (saída
  crescente com a saturação) e a carga da lateral direita, que não são dadas.
- A parte de incerteza do artigo não foi feita: triagem de Plackett-Burman com 14 parâmetros, PCE com índices de
  Sobol, os 500 casos com ±10 % e a probabilidade de 21 %, e o caso arenoso com o ramo de umedecimento.
- Uma variante com a tensão efetiva no cálculo elástico (α = 1, o modelo acoplado do projeto), para contrastar com a
  abordagem de tensão total do artigo.
