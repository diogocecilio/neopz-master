# Plano estocástico da chuva ("Seção 6 da chuva") — EM STAND-BY

Estado em 25/11/2026: **em stand-by** até a validação determinística completa (`plano_deterministico.md`). Motivo
(decisão do usuário, 25/11): se o determinístico não for totalmente confiável, a estatística sai errada e distorcida.

## Objetivo

O análogo das Seções 6.1 a 6.5 de Vargas Ceron et al. (2025) no talude de Rahardjo et al. (2007) (45°, H = 10 m,
H_w = 5 m), com campos aleatórios de c′, φ′ e k_s (estatísticas da Tabela 2 do Ceron: CoV 0,3 / 0,1 / 0,6, distâncias
de autocorrelação 20 m × 2 m) e as três formas do modelo:

- (iii) desacoplado, esqueleto rígido, ar passivo (as hipóteses do artigo; a forma principal);
- (ii) acoplado (esqueleto elástico), ar passivo: só nas referências (fica a ≤ 0,03 da (iii) no determinístico);
- (i) acoplado (esqueleto elástico), ar ativo: a comparação principal é (iii) × (i).

## O que já foi feito (até 25/11)

| item | resultado |
|---|---|
| Referência do f50, (iii), 31,05 mm/h × 24 h, MC direto de 5000 (fermi, 532c105) | P_f = 5,58 % (CoV 5,8 %); menor FS em 24 h em 4981/5000; Spearman com o FS: c′ da banda +0,66, φ′ +0,46, k_s −0,01 |
| Varreduras determinísticas (runs/chuva2_det_solos) | chuvas de referência com o mesmo FS das médias (≈ 1,12): f10 80 mm/h × 24 h; f100 9,4 mm/h × 72 h |
| Reaproveitamento do transiente no LS (f5d1cbe) | linhas iguais bit a bit, metade do tempo |
| Parada pelo CoV no MC (e4f07f4) | `cov_alvo=` conta o menor FS de cada realização |
| Variante (i) com campo aleatório (runs/variante_i_campo, 16 realizações da ref) | FS da (i) 0,12–0,25 abaixo da (iii) (média −0,18); 13 de 16 rompem; P_f(i) estimada 50–80 % na chuva da referência |

Comando do fermi preparado e **não lançado** (MC das referências do f10/f100, complemento do f50 e LS do f50): ver
`cluster/casos_chuva2.txt` e a memória da campanha.

## Regras do método

- Falha: o menor FS entre o fim da chuva e 6 h depois < 1.
- **P_f ≥ 5 %: Monte Carlo direto** até CoV(P_f) = 5 % (parada automática `cov_alvo=0.05`;
  n ≈ (1 − P_f)/(0,0025 P_f): 7600 com 5 %, 3600 com 10 %, 1600 com 20 %, 300–400 com 50–60 %).
- **P_f < 5 %: line sampling** na direção do FORM modal em (c′, φ′) até CoV 5 %.
- Decisão por piloto: as 500 primeiras realizações do MC; com P_f do piloto ≥ 4 % o MC segue, se não o caso vai para LS.
- Custo: o MC com P_f de 5–10 % custa ~15–26 nós·hora por caso; um LS ~3 nós·hora.

## Etapas

| etapa | o quê | onde | tempo |
|---|---|---|---|
| E1 | Referências (iii): completar o f50; f10 (80 mm/h × 24 h) e f100 (9,4 mm/h × 72 h) | fermi, 3 nós | 1 noite |
| E2 | LS × MC na referência do f50 (P_f 5,6 % está na fronteira da regra: é a validação do LS) | fermi, 1 nó | mesma noite |
| E3 | (i) nas 3 referências (MC direto, P_f alta, ~300 realizações); (ii) só no f50 | fermi | ~1 noite |
| E4 | Rotação dos campos: núcleo de Nyström em coordenadas giradas (`nystrom_base.py`) e uma base de KL por campo no driver (girar um campo sem girar os outros). Não mexe no material | PC | ~1 dia |
| E5 | Famílias no f50, (iii) (lista abaixo), ~40 casos | fermi | ~4–5 noites |
| E6 | (i) na família da intensidade: curvas P_f(I) da (iii) e da (i) | fermi | 1–2 noites |
| E7 | Verificações vindas do relatório de outubro (a causa da diferença de ~0,1 na (iii); o item A2 com a (i) elástica) — **absorvidas no plano determinístico** | PC/fermi | — |
| E8 | Relatório no formato da Seção 6 do Ceron + comparação das variantes | — | 2 dias |

### Famílias (E5), valores do Ceron

- 6.2 CoV: c ∈ {0,1; 0,5; 0,7}; φ ∈ {0,05; 0,15; 0,2}; k ∈ {0; 0,75; 0,9; 1} (10 casos).
- 6.3 s: {1,5; 2; 5; 10; 20; 400} (6 casos).
- 6.4 anisotropia: (a) α = k_h/k_v ∈ {2, 3, 4, 5} (4 casos); (b) rotação da estrutura de correlação de cada campo
  isolado (c′, φ′, k_s), θ ∈ {−45° (paralela à face), −22,5°, +22,5°, +45°} (12 casos), e um caso com os três girados
  juntos no pior ângulo.
- 6.5 intensidade: I/k_s ∈ {0,5; 0,75; 1; 1,5; 2} (5 casos).

Total: ~1 semana de fermi com 4 nós, dominada pelos casos com P_f > 5 % (MC).

## Decisões em aberto

- **D1 (importante): Nyström ou KL analítica.** A referência do f50 rodou com a KL analítica (M = 2000 num retângulo);
  a regra do usuário de 22/11 é Nyström em tudo, e a rotação precisa de Nyström. Recomendação: Nyström em tudo e
  refazer a referência do f50 (+1 noite). Resolver antes de lançar o fermi.
- D2: os ângulos da rotação.
- D3: incluir também a rotação do tensor de permeabilidade nas famílias (no determinístico ela entra no plano
  determinístico).
- D4: famílias do f10 e do f100 (só a intensidade, ou outras). No f50 o k_s quase não pesa (I < k_s); a rotação de k_s
  só deve aparecer no f100 (I > k_s).

## Comparação com o relatório de outubro (docs/nao_saturado/biot_nao_saturado_agua_ar5.pdf)

FS mínimo determinístico, 24 h de chuva:

| caso | Eq. 6 (artigo) | relatório (iii) | hoje (iii) | relatório (i) | hoje (i) |
|---|---|---|---|---|---|
| f10, 90 mm/h | 1,05 | 1,15 | 1,07 | 1,15 | 1,08 |
| f10, 180 mm/h | 0,85 | 0,86 | 0,84 | 1,00 | 0,90 |
| f50, 18 mm/h | 1,49 | 1,61 | 1,47 | 1,32 | 1,18 |
| f50, 36 mm/h | 1,05 | 1,18 | 1,03 | 1,01 | 0,92 |
| f100, 9 mm/h | 1,55 | 1,69 | 1,55 | 1,36 | 1,27 |

1. A (iii) de hoje fica 0,08–0,15 abaixo do relatório e cai sobre a Eq. 6 (≤ 0,02). O relatório atribuía a diferença
   para o artigo a "SRM × Bishop e digitalização"; com o mesmo SRM ela some. O FS inicial muda só 0,026
   (2,119 → 2,093): a maior parte deve vir das pressões do transiente antigo. Causa não identificada.
2. O patamar FS ≈ 1,00 do modelo completo não existe mais ((i) desce a 0,90–0,92): era a dilatância plástica no
   transiente acoplado (o próprio item A1); hoje o esqueleto do transiente é elástico.
3. O efeito do ar ativo é robusto: (i) − (iii) = −0,29 no f50 a 18 mm/h nos dois; −0,33 × −0,28 no f100.
4. (ii) ≈ (iii) é robusto (≤ 0,03 nos dois).
5. Possível contraexemplo novo (A2): no f100 a 63,4° com ar ativo o relatório deu 1,01–1,03 (no patamar plástico); com
   o esqueleto elástico pode cair abaixo de 1.
6. Intensidade crítica do f50: (iii) 48 mm/h no relatório, hoje entre 36 e 54 (provavelmente ≈ 37–38); (i) 26–31,
   hoje ≈ 30.
7. Com campo aleatório (não estava no relatório): a (i) baixa o FS de cada realização em 0,12–0,25.

## Proposta (25/11): curva de fragilidade P_f(I_r) a partir da sigmoide do Rahardjo

Ideia do usuário: em vez da sigmoide do FS mínimo (Eq. 6 e Figs. 6–7 do Rahardjo), uma sigmoide da probabilidade de
FS ≤ 1 em função da chuva, ligando as constantes A, B, C, D às incertezas do Monte Carlo.

1. **Determinístico.** Com F(I) = A + B / (1 + e^{-(I − C)/D}) (A: patamar superior, A + B: patamar inferior, B < 0),
   a ruptura FS ≤ 1 ocorre a partir da intensidade crítica I_c = C − D ln(B/(1 − A) − 1), se A + B < 1 < A. O
   análogo da Fig. 6 para FS ≤ 1 é a curva t₁(I) (primeira passagem com chuva constante): a duração que rompe o
   talude com intensidade I, ou seja, o limiar intensidade–duração (I–D) do talude.
2. **Probabilístico.** Com o FS mínimo lognormal de mediana F(I) (a sigmoide) e dispersão β(I) = desvio de ln FS:
   P_f(I) = Φ(−ln F(I) / β(I)). Ela também é uma sigmoide: 50 % em I_c e, linearizando perto de I_c,
   P_f ≈ Φ((I − I_c)/σ_I) com σ_I = β / |d ln F/dI| em I_c. A incerteza "borra" o degrau determinístico em I_c numa
   rampa de largura σ_I. São 5 constantes: A, B, C, D (determinístico ou mediana do MC) e β (do MC).
3. **Teste com o que já existe** (referência do f50, 31,05 mm/h, MC de 5000): β = 0,0768, ln FS quase simétrico
   (assimetria 0,08). Com a mediana do MC (1,1271): Φ(−ln F/β) = 5,96 %, contra 5,58 % do MC (IC 95 % 5,0–6,3 %).
   Com o FS das médias (1,1239, o determinístico): 6,41 %, 0,8 ponto acima (a mediana do FS é um pouco maior que o FS
   das médias).
4. **Duração de graça, pelo t₁.** Com chuva constante, uma realização que rompe em t₁ rompe também em qualquer chuva de
   duração D ≥ t₁ com a mesma I. Então um MC por intensidade dá P_f(I, D) = P(t₁ ≤ D) para todas as durações até o
   horizonte: a superfície de fragilidade e os limiares I–D probabilísticos (ex.: a duração com P_f = 5 %). Ressalva:
   vale para as rupturas durante a chuva; as que ocorrem na drenagem depois da chuva dependem do fim da chuva.
5. **O que testar:**
   - β(I) muda com I: abaixo de k_s o k_s não pesa (Spearman ≈ 0 na referência), acima dele (alagamento) deve pesar;
   - a mediana do FS contra o FS das médias (efeito da variabilidade espacial);
   - a cauda lognormal;
   - a normalização I/k_s entre os solos (a do próprio Rahardjo).
6. **Como fazer:**
   - A, B, C, D vêm da varredura determinística (T1 da campanha do fermi, 6–9 intensidades por solo);
   - β vem do MC em 3–5 intensidades (a família 6.5);
   - a validação é a P_f prevista contra um MC novo numa intensidade intermediária.
   - Produto: uma fórmula fechada P_f(I_r) (e P_f(I, D)) para um talude e um material, com 5 constantes.

### Campanha frag (25/11, foco total; o resto do fermi parado)

- **Sigmoides determinísticas** (FS das médias no fim da chuva):
  - f50 (iii) em runs/sigmoide_f50: FS = 1 em ~39,8 mm/h, com uma quina em I = k_s = 36 mm/h (alagamento);
  - f50 (i), no mesmo lugar: FS = 1 em ~23 mm/h, com patamar de 0,91 a partir de ~30 mm/h, de modo que a P_f satura em ~85–90 %;
  - f100 com chuva de 72 h em runs/sigmoide_f100_72h.
- **Método: primeira passagem t₁ por realização e intensidade** (`crpc=` do SlopeRainMonteCarlo, estendido às variantes acopladas). Em cada intensidade, decrescente, a chuva é constante até crpc_dmax. O FS é avaliado em 6, 12, 24, 48 e dmax h, com bissecção nas horas guardadas. Quem não rompe numa intensidade pula as menores. Isso dá P_f(I, D) = P(t₁(I) ≤ D) para todo D ≤ dmax.
- **Desenho:** as mesmas realizações 0..1999 (os campos da referência chuva2), N = 2000.
  - f50: 24 h, intensidades (iii) 45, 40, 36, 33, 31,05, 30, 27 e (i) 36, 27, 24, 22, 20, 18.
  - f100: 72 h.
  - β de calibração: (iii) f50 na referência; (i) f50 num MC a 20 mm/h; f100 num MC por variante.
- **Conferências embutidas:** o crpc a 31,05 (iii) contra a referência, e a 20 (i) contra o MC de calibração, nas mesmas realizações.
- **Análise:** cluster/resumo_frag.py compara a P_f do MC (IC de Wilson) com Φ(−ln F/β), usando os pontos e a Eq. 6, e dá o z por intensidade e a distribuição de t₁.
