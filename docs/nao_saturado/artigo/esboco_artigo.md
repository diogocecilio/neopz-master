# Artigo principal: resumo, esboço e cálculos que faltam

Estado em 29/10/2026, sobre o branch `claude/relaxed-volta-1r4gog` (commit 9e33be7). Os números vêm do documento
`biot_nao_saturado_agua_ar5.pdf` (seções de contraexemplos, literatura, benchmarks e patamar) e do relatório
`relatorio_discussao.pdf` (`secao_discussao.tex`). Valores marcados com † são anteriores às correções do avanço no
tempo e precisam ser refeitos antes de entrar no artigo.

**Título de trabalho:** *Critical rainfall intensity of unsaturated slopes: how much do uncoupled, passive-air analyses
overestimate it?*

Alternativa: *Trapped air, coupling and dilatancy lower the critical rainfall intensity of unsaturated slopes.*

**Revistas-alvo:** Engineering Geology, Computers and Geotechnics, Landslides, Acta Geotechnica.

---

## 1. Resumo

### Abstract (EN, rascunho; ~280 palavras)

Rainfall thresholds and early-warning criteria for unsaturated slopes are usually derived from uncoupled seepage
analyses with atmospheric pore-air pressure, followed by a limit-equilibrium or strength-reduction stability analysis.
We quantify how these hypotheses bias the rainfall that brings a slope to failure. A fully coupled hydro-mechanical
formulation with three fields (displacement, pore-water and pore-air pressures; Bishop effective stress, Mohr–Coulomb
plasticity) is implemented in the open-source finite element library NeoPZ, and the factor of safety (FS) is
evaluated by strength reduction at successive instants with the pressures frozen. Since values below one do not
describe physical states, models are compared by the critical intensity *I_c*, the smallest intensity of a 24-h
rainfall that brings the FS to one, and by the time to failure. For the homogeneous slope of Rahardjo et al. (2007),
the uncoupled, passive-air model gives *I_c* ≈ 48 mm/h (1.3 *k_s*), whereas the coupled model with active air gives
31 mm/h (0.9 *k_s*) with an associated flow rule and 26 mm/h (0.7 *k_s*) without dilatancy: the usual hypotheses
overestimate the failure-triggering rainfall by 55–86 % and delay the predicted failure by 2–10 h in intense
rainfalls. Air trapped between the wetting front and the water table is compressed to 36–38 kPa; less water
infiltrates and the matric suction increases, but the pore-water pressure follows the air pressure and the effective
stress drops, so that a low-permeability slope classified as safe under short rainfalls reaches limit equilibrium
(minimum FS 1.009–1.015, decreasing with mesh refinement). Dilatancy does not decide whether a critical intensity
exists, but shifts it by about 20 % and controls the post-failure response: without dilatancy the transient loses
equilibrium, with ψ ≥ φ'/2 it stops at limit equilibrium. The lower FS plateau of intense rainfalls is the FS of the
slope saturated to its surface, a boundary-value property independent of *k_s*. Implications for physically based
rainfall thresholds are discussed.

*Números que dependem dos cálculos pendentes:* a decomposição entre acoplamento e ar (P1.1), a generalização para
outros solos e para χ = S_e (P1.3) e a convergência de malha de *I_c* (P1.4). Se a generalização mudar a faixa de
55–86 %, o resumo passa a dar a faixa geral e o talude de Rahardjo vira o caso de referência.

### Resumo (PT)

Limiares de chuva e critérios de alerta para taludes não saturados costumam vir de análises de percolação
desacoplada com pressão do ar atmosférica, seguidas de equilíbrio-limite ou redução de resistência. Quantificamos o
quanto essas hipóteses deslocam a chuva que leva o talude à ruptura. Uma formulação hidromecânica acoplada com três
campos (deslocamentos, pressões da água e do ar; tensão efetiva de Bishop e Mohr–Coulomb) foi implementada na
biblioteca de elementos finitos NeoPZ, e o FS é calculado por redução de resistência em instantes sucessivos com as
pressões congeladas. Como valores abaixo de 1 não descrevem estados físicos, os modelos são comparados pela
intensidade crítica *I_c* (a menor intensidade de uma chuva de 24 h que leva o FS a 1) e pelo instante da ruptura. No
talude homogêneo de Rahardjo et al. (2007), o modelo desacoplado com ar passivo dá *I_c* ≈ 48 mm/h (1,3 *k_s*); o
acoplado com ar ativo, 31 mm/h (0,9 *k_s*) com escoamento associativo e 26 mm/h (0,7 *k_s*) sem dilatância. As
hipóteses usuais superestimam em 55 a 86 % a chuva que rompe o talude e atrasam a ruptura em 2 a 10 h nas chuvas
fortes. O ar aprisionado entre a frente de umedecimento e o NA é comprimido a 36–38 kPa: infiltra menos água e a sucção
matricial aumenta, mas a poropressão acompanha a pressão do ar e a tensão efetiva cai. Assim, um talude pouco
permeável, classificado como seguro em chuvas curtas, chega ao equilíbrio limite (FS mínimo 1,009–1,015, caindo com o
refinamento da malha). A dilatância não decide se existe intensidade crítica, mas a desloca em cerca de 20 % e controla
o que acontece depois: sem dilatância o transiente perde o equilíbrio; com ψ ≥ φ'/2 ele para no equilíbrio limite. O
patamar inferior do FS nas chuvas intensas é o FS do talude saturado até a superfície, uma propriedade do problema de
contorno que não depende de *k_s*. Discutem-se as consequências para limiares de chuva de base física.

---

## 2. Mensagem e contribuições

**Mensagem única:** as hipóteses da prática (percolação desacoplada, ar passivo, escoamento associativo) superestimam
a intensidade crítica e atrasam a ruptura; o efeito maior parece vir do ar aprisionado (a confirmar com P1.1).

1. **Métrica.** Com o critério FS ≤ 1, comparar modelos pela intensidade crítica *I_c* e pelo instante da ruptura
   *t_f*, e não pelo FS mínimo. O FS mínimo abaixo de 1 é o patamar do talude saturado (0,84–0,85 nos três solos),
   uma propriedade do problema de contorno.
2. **Viés do modelo usual,** quantificado: *I_c* superestimada em 55–86 % (28–54 % contra a Eq. 6 do artigo, que dá
   FS = 1 com 40 mm/h) e ruptura 2–10 h mais tarde nas chuvas de 54 a 180 mm/h.
3. **Mecanismo do ar aprisionado.** A sucção matricial sobe (38 contra 24 kPa entre 2 e 5 m de profundidade), mas
   *p_w* sobe com *p_a* (−6 contra −24 kPa), e com χ = 1 a tensão efetiva só depende de *p_w*. Resultado: menos água
   (10,4–11,6 contra 13,5–17,9 m³/m em 24 h) e FS menor. A categoria "solo pouco permeável é seguro em chuva curta"
   de Rahardjo et al. não se sustenta: G1, f100, α = 63,4°, FS mínimo 1,009–1,015.
4. **Dilatância** em ruptura por chuva com modelo acoplado (não encontramos trabalho que a varie). Ela desloca *I_c*
   em ≈ 20 % (de 26 para 31 mm/h) e decide entre colapso (ψ = 0) e parada no equilíbrio limite (ψ ≥ φ'/2).
5. **Secundário,** aproveitável numa seção curta:
   - o patamar calculado direto, com a superfície alagada;
   - a anisotropia: com k₁ paralelo à face a ruptura vem 8–13 h antes nas chuvas de 90 e 180 mm/h;
   - a duração crítica *t\** = *V*/*k_s*.

### O que já existe (posicionamento; conferir numa revisão completa)

- **Ar.** A pressão do ar reduz o FS e o modelo monofásico superestima a segurança: Cho (2016, Engineering Geology);
  *Effects of airflow induced by rainfall infiltration on unsaturated soil slope stability* (Transport in Porous
  Media, 2015); ensaio de campo com modelo bifásico (2021); loess (CATENA, 2026). Em geral fluxo bifásico
  **desacoplado** da deformação, com equilíbrio-limite.
- **Acoplamento.** Análises acopladas dão FS menor que as desacopladas (Qi e Vanapalli 2015, Computers and
  Geotechnics; Yang et al. 2017, Engineering Geology).
- **Os nove trabalhos já confrontados** (seção de literatura do documento): todos usam *p_a* = 0; nenhum varia ψ.
- **Lacunas que o artigo ocupa:**
  - quantificar o deslocamento de *I_c* e de *t_f*;
  - o modelo de três campos acoplado à deformação, com SRM;
  - a dilatância;
  - a explicação quantitativa do "menos água, FS menor".

---

## 3. Esboço

| seção | conteúdo | material existente | falta |
|---|---|---|---|
| 1. Introdução | limiares de chuva e alerta; prática SEEP/W + SLOPE/W; literatura de ar e de acoplamento; lacuna; objetivos | seção de literatura (T0–T5) | revisão bibliográfica de fluxo bifásico e acoplamento |
| 2. Formulação | u–p_w–p_a; tensão de Bishop com χ; Mohr–Coulomb associativo e não associativo (PV); Fredlund–Xing; condição unilateral de chuva; ar só pela superfície; Euler implícito, Newton consistente, busca linear e reinício de passo; SRM com pressões congeladas e tolerância | seções de formulação do documento ar5 e "Material e análise" | versão em inglês, enxuta |
| 3. Verificação e validação | Rahardjo (série B (iii)); Zhang 2014 (0–5 %); Bishop × SRM (2,092 × 2,119); malha (G4) | tab:B, tab:zhang, Tang/Bishop, `fase_ar_G4_malha.png` | **Liakopoulos (1965)** e **Touma e Vauclin (1986)** (P1.2) |
| 4. Critério e métrica | FS ≤ 1; FS < 1 não físico; definição de *I_c* (interpolação ou bissecção) e *t_f*; resolução (tolfs 0,002; discretização ≈ 0,01); o patamar como FS do talude saturado | `secao_discussao.tex` 1.1–1.2; `patamar_alagado.png` | figura FS(t) × *I_r* com o cruzamento de 1 |
| 5. Decomposição do viés (f50, 45°) | quatro modelos: (iii) artigo, (ii) acoplado com ar passivo, (iv) desacoplado com ar ativo, (i) completo; *I_c* e *t_f* de cada um | A1c: (i) com ψ = 0 e ψ = φ', (iii) com 48 mm/h; série B (i) e (ii)† | P1.1 |
| 6. Generalização | *I_c*(i)/*I_c*(iii) para f10, f50, f100; α = 30°, 45°, 60°; χ = 1 e χ = S_e | série B (iii), A5 | P1.3, P2.3 |
| 7. Ar aprisionado | perfis de *p_a*, *p_a − p_w* e *p_w*; tempos de alagamento; volume infiltrado; G1 (f100, 63,4°); k_ra,min; malha | `fase_ar_campos.png`, `fase_ar_G1_fs_tempo.png`, tab:A2, G3, G4 | P1.4 (malha 0,25 m), P2.1 (laterais abertas ao ar) |
| 8. Dilatância | ψ = 0, φ'/2 e φ': *I_c*, *t_f*, colapso × parada; FS (SRM com ψ = φ') de 1,05–1,06 na perda de equilíbrio | tab:A1, `contra_A1_criterio.png` | P2.2 (passo de tempo e *t_f*), P3.1 (SRM não associativo) |
| 9. Consequências | limiares intensidade–duração (I–D) dos modelos (i) e (iii); categorias de Rahardjo; quando o modelo usual basta | A6 (*t\** = *V*/*k_s*) | P2.4 |
| 10. Limitações | 2D, talude homogêneo, lei de χ, sem histerese, porosidade constante, domínio fechado ao ar (limite superior do aprisionamento), SRM associativo, malha | seção "Limitações" do documento | — |
| 11. Conclusões | — | — | — |

### Figuras previstas

1. Talude, malha e condições de contorno; parâmetros (tabela).
2. Validação: Touma e Vauclin (*p_a* e frente de umedecimento) e Liakopoulos (drenagem, recalque).
3. FS(t) para várias intensidades nos modelos (i) e (iii), com o cruzamento de FS = 1.
4. *I_c* dos quatro modelos (barras) e generalização *I_c*(i)/*I_c*(iii) por solo, ângulo e lei de χ.
5. Instante da ruptura × *I_r* nos modelos (i) com ψ = 0, (i) com ψ = φ' e (iii).
6. Campos de *p_w* e *p_a* com ar ativo e passivo; perfis médios por profundidade.
7. Malha × FS mínimo (G4) e malha × *I_c*.
8. Limiares I–D dos modelos (i) e (iii).

---

## 4. Cálculos

### Rodados em 29/10/2026 (`scripts/artigo.py`, 31 casos; versão do artigo em `artigo.tex`)

| grupo | casos | resultado |
|---|---|---|
| D, decomposição (f50, 45°) | (ii) ψ = φ' com 54, 66, 78; (ii) ψ = 0 com 42, 54; (iv) com 20, 25, 30, 36, 42, 48 mm/h | *I_c* = 78 (62 com FS ≤ 1,02), 53 e 25 mm/h; com (iii) 48, (i) ψ = φ' 31 e (i) ψ = 0 26: o ar reduz *I_c* à metade, o acoplamento o aumenta 4–63 % |
| G, generalização | f10 em (i) com 100, 120, 140, 165 e em (iii) com 120, 140; f50 com χ = S_e em (i) com 24, 30, 36 e em (iii) com 36, 42, 48, 54 mm/h | f10: 135 × 139 mm/h (sem viés, *I_c* ≈ 0,4 *k_s*); χ = S_e: 48 × 30 mm/h (o viés persiste) |
| M, malha (h_sup = 0,35 m) | (i) com 30, (iii) com 48 mm/h | FS mínimo cai 0,007 e 0,024; *I_c* ≈ 30 e 46 mm/h (razão 1,5 mantida) |
| T, passo de tempo | (i) ψ = 0, 36 mm/h, passo 0,125 h | perda de equilíbrio às 15,6 h (15,5 h com 0,25 h) |
| Ar, ar aberto na base e nas laterais | G1 (f100, 63,4°, 36 mm/h); f50 (i) com 30 mm/h | FS mínimo 1,017 × 1,015 e 1,011 × 1,011: o ar é vedado pela superfície alagada |
| L, validação | coluna de Liakopoulos, ar ativo e passivo | ar passivo × ensaio (vazão e p_w, Callari e Abati 2009, Figs. 3a, b); a comparação com Khoei e Mohammadnejad (2011) foi retirada na revisão de 30/10/2026 (PDF não disponível; ver auditoria_referencias.md) |

Código acrescentado: leis de Liakopoulos em `TPZUnsaturatedHydraulics` (`swcc=liak`, `krw=liak`); opções `arbase`,
`arlateral`, `tsnap`, `parar` e a vazão pela base em `hist.csv` (`SlopeRainfall`). O caso `M_i` foi encerrado depois do
último FS pedido (26 h), sem `summary.txt` (`PARADO.txt`).

### Ainda em aberto

- Validação da fase de ar em infiltração com pressão do ar medida (Touma e Vauclin 1986): o artigo não está na pasta
  de referências.
- Malha: o FS mínimo do G1 cai 0,003 por refinamento sem convergir (h_sup = 0,25 m não foi rodado).
- Limiares intensidade–duração (I_c para 6, 12 e 48 h) e outros ângulos (P2.3, P2.4): não rodados.
- SRM não associativo (P3.1) e anisotropia no modelo completo (P3.2): não feitos.
- Conferir as referências marcadas com %VERIFICAR em `referencias.tex`; afiliação e coautores.

## 5. Referências a conferir

- Rahardjo, H., Ong, T.H., Rezaur, R.B., Leong, E.C. (2007). Factors controlling instability of homogeneous soil
  slopes under rainfall. J. Geotech. Geoenviron. Eng. 133(12).
- Cho, S.E. (2016). Stability analysis of unsaturated soil slopes considering water-air flow caused by rainfall
  infiltration. Engineering Geology.
- *Effects of airflow induced by rainfall infiltration on unsaturated soil slope stability.* Transport in Porous
  Media (2015), doi:10.1007/s11242-015-0469-x.
- Qi, S., Vanapalli, S.K. (2015). Hydro-mechanical coupling effect on surficial layer stability of unsaturated
  expansive soil slopes. Computers and Geotechnics.
- Yang, K.H. et al. (2017). Coupled hydro-mechanical analysis of two unstable unsaturated slopes subject to rainfall
  infiltration. Engineering Geology.
- Touma, J., Vauclin, M. (1986). Experimental and numerical analysis of two-phase infiltration in a partially
  saturated soil. Transport in Porous Media 1.
- Liakopoulos, A.C. (1965). Transient flow through unsaturated porous media. PhD thesis, University of California,
  Berkeley.
- Tschuchnigg, F., Schweiger, H.F., Sloan, S.W. (2015). Slope stability analysis by means of finite element limit
  analysis and finite element strength reduction techniques (Partes I e II). Computers and Geotechnics.
- Os nove trabalhos da seção de literatura (Ali 2014, Zhang 2014, Dou 2015, Tang 2018, Wang, Vardon e Hicks 2018,
  Ran 2018, Liu e Wang 2021, Wang, Idinger e Wu 2021, Pradhan 2022).
