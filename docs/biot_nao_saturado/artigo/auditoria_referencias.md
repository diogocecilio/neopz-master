# Auditoria do artigo (revisão de 30/10/2026)

Pedido: reler as referências anexadas, corrigir o necessário e verificar os erros do artigo, sem afirmar nada
sem citação e sem inferir nada, indicando os nossos resultados e as contradições da literatura.

Método:
- **Literatura:** cada afirmação atribuída à literatura foi conferida no texto do PDF anexado. O texto foi extraído
  página a página e as figuras foram vistas como imagem quando necessário.
- **Página citada:** a do periódico. Nos artigos com número de artigo, é a página do PDF.
- **Referências sem PDF:** a referência bibliográfica foi conferida na lista de referências de um trabalho anexado.
  O conteúdo dessas obras não foi usado.
- **Nossos números:** foram conferidos nos arquivos de resultados (`Projects/SlopeRainfall/resultados`).

## 1. Erros do manuscrito e correções

| # | Onde | Erro | Correção / evidência |
|---|---|---|---|
| 1 | Eq. (3), massa de ar | Faltava o termo $\mathbf q_a\cdot\nabla p_a/(p_a+p_{atm})$, que o código usa (`TPZMatUnsaturatedUPwPa.cpp`, termo `-dt*cag*la*adv`) | Termo incluído; corresponde a $(\dot{\mathbf w}^g/\rho_g)\cdot\nabla\rho_g$ de Ghorbani et al. (2018, Eq. 35) e ao último termo da Eq. 36 de Li et al. (1999, p. 139), ambos vindos de $\nabla\cdot(\rho_g\ldots)$ (Li et al. 1999, Eq. 33) |
| 2 | Discretização | "forma secante conservativa" valia só para a água; o termo do ar é $nS_a^{k+1}\Delta p_a/(p_a^{k+1}+p_{atm})$, que não conserva exatamente a massa de ar | Texto corrigido; balanço global de água e ar calculado (seção 2.3) |
| 3 | Introdução | Pradhan et al. (2022) agrupado na cadeia desacoplada | Pradhan usa análise "fully coupled flow-deformation" no PLAXIS 2D, Biot, χ = S_e, p_a desprezada e SRM (p. 4) — corresponde ao nosso modelo (ii) |
| 4 | Introdução | Cho (2016), Touma e Vauclin (1986), Qi e Vanapalli (2015), Yang et al. (2017), Ng e Shi (1998), Tschuchnigg et al. (2015) citados por conteúdo sem PDF | Retirados; as afirmações foram refeitas com citações conferidas (Sheng 2008, Fredlund 2006, Ghorbani 2018, Khoei e Saeedmonir 2021, Callari e Abati 2009, Liu e Wang 2021, Pedroso 2015) |
| 5 | Introdução | "o SRM é usualmente associado"; "como nos outros estudos de infiltração com ar" | Retirados (sem fonte) |
| 6 | Seção 2.4 | Zienkiewicz et al. (1975) e Tschuchnigg et al. (2015) sem PDF | Definição do FS citada de Pradhan et al. (2022, p. 4) e Griffiths e Lane (1999); a escolha ψ = φ' no SRM ficou como escolha nossa, sem justificativa bibliográfica |
| 7 | Liakopoulos | Comparação numérica com Khoei e Mohammadnejad (2011) (`scripts/ref/liakopoulos_khoei2011.csv`, "Fig. 1, leitura"): o PDF não foi anexado e nenhum trabalho anexado contém essa solução | Comparação, tabela e marcadores retirados; os valores ficaram no repositório, mas não são usados |
| 8 | Liakopoulos | "o ensaio verifica a fase ar"; "característica da lei de permeabilidade do ar" | Retirados: o ensaio verifica só a implementação com ar passivo; a solução bifásica drena mais devagar que o ensaio |
| 9 | Liakopoulos | Propriedades atribuídas a Lewis e Schrefler (1998), sem PDF | Propriedades citadas de Callari e Abati (2009, Tabela 1 e Eq. 58). O coeficiente de retenção usado antes (0,10152 com γ_w = 9,81) era 4,9% maior que o de Callari (1,9722·10⁻¹¹ Pa⁻²·⁴²⁷⁹); os dois casos foram refeitos com o de Callari (`resultados/artigo_rev/L_liakopoulos_*_callari`): o erro RMS da vazão passou de 0,30 para 0,29·10⁻⁶ m/s |
| 10 | Critério | "FS < 1 não é um estado físico" | Redação do revisor: "a value below one indicates that the available strength is insufficient for the analysed configuration…" |
| 11 | Critério | I_c definido só com 1,01; nos modelos acoplados associados o FS não cai abaixo de 1 (o transiente usa a resistência real; nos desacoplados usa F = 0,25 — `main.cpp`) | I_c,δ para δ = 0,02; 0,01; 0,005; 0 (tabela 3); I_c(δ=0) não existe no modelo (i) associado (FS ≥ 1,0002 até 180 mm/h) |
| 12 | Critério | I_r(thr) de Rahardjo tratado como intensidade crítica | Separados: I_r(thr) é a inflexão inferior da sigmoide (p. 1539), 51 mm/h = 1,42 k_s no f50 |
| 13 | Critério | Patamar alcançado "após 7, 24 e 170 h" | Os instantes calculados só permitem dizer: entre 6 e 12 h, 12 e 24 h, 96 e 192 h (`resultados/benchmarks/patamar`) |
| 14 | Critério | Erro de discretização "cerca de 0,01" | Medido: 0,003–0,023 |
| 15 | Verificação | Malha padrão "1684 triângulos" vs "1834 elementos" (seção 8) e "1764, 3814, 7622 triângulos" (seção 6) — eram elementos incluindo linhas de contorno | Triângulos: 1684 e 7732 (45°); 1620, 3597, 7314 (63,4°) |
| 16 | Verificação | "segue a sigmoide de Rahardjo dentro de 0,01–0,14" | Fica 0,01–0,14 **acima** dela, nunca abaixo |
| 17 | Verificação | Zhang et al. (2014) "seis leis, 0–5%" | Zhang compara oito leis; reproduzimos as três lineares (φ_b = 0, 15°, φ') nos taludes íngreme e plano (Tabelas 4 e 5) com 0,05–3,9%. A lei χ = S_w foi comparada com Fredlund et al. (1996), mas o κ usado por Zhang não consta no texto: comparação retirada |
| 18 | Seção 5 | "atraso de até 10 h" | 4,5 h (ψ = φ') e 9,9 h (ψ = 0) com 54 mm/h; FS calculado a cada 3 h nesses casos |
| 19 | Seção 5 | "A dilatância antecipa a ruptura" | Sem dilatância a ruptura ocorre 5,4–6,3 h antes |
| 20 | Seção 5 | Mecanismos (inchamento elástico, dilatação abrindo poros, χ → 1 na zona de ruptura, sucção gerada pela dilatação) apresentados como fato | Retirados; o texto diz o que foi calculado e que o mecanismo não foi isolado |
| 21 | Seção 5 | ψ = φ'/2 com I_c = 36 mm/h | Só 18 e 36 mm/h foram calculados: I_c entre 18 e 36 |
| 22 | Seção 6 | "a partir de I_r = 5 k_s o talude atinge o equilíbrio limite" | Com 18 mm/h (5 k_s) FS = 1,028; dentro de 0,02 de 1 só a partir de 36 mm/h (1,015; 1,009 na malha mais fina): ruptura com δ = 0,02, não com δ = 0,01 |
| 23 | Seção 6 | "a sucção aumenta" (soava temporal) | Comparação entre modelos: 38 kPa com ar ativo contra 24 kPa com ar passivo |
| 24 | Seção 6 | Volume infiltrado com ar passivo 13,5–17,9 | 13,5–18,0 m³/m |
| 25 | Seção 6 | FS inicial "1,748 com 0,25 m" | Sem arquivo de resultado no repositório (só no relatório de discussão): retirado |
| 26 | Seção 8 | f10: "a superfície não alaga antes da ruptura e o ar escapa" | **Falso nos nossos resultados:** com 140 mm/h (modelo i) o alagamento começa às 6,75 h, 87 vértices alagados quando FS = 1,01 (13 h), p_a = 40 kPa. O motivo da ausência de viés no f10 não foi estabelecido |
| 27 | Seção 8 | Rodadas com `parar=1.01` usadas para FS mínimo (G_i f10 140/165, χ = S_e 30/36, D_iv 25–48, D_ii 78) | O FS mínimo dessas rodadas é só o primeiro valor ≤ 1,01; refeitas sem parada as necessárias (`resultados/artigo_rev`); onde não refeito, a tabela marca "≤" (limite superior) ou "n.c." |
| 28 | Limitações | "k_ra,min = 10⁻² eleva o FS em 0,035" | 0,031 (1,0153 → 1,0459) |
| 29 | Seção 9 (implicações) | "evento muito mais frequente" sem curva IDF | Retirado; só a afirmação qualitativa (menor intensidade, mesmo período ⇒ menor período de retorno) |
| 30 | Código | URL labmec/neopz | https://github.com/diogocecilio/neopz-master |

## 2. Afirmações da literatura usadas (todas conferidas no PDF)

| Referência | Página | Trecho / dado |
|---|---|---|
| Rahardjo et al. 2007 | 1532 | "the initial water table location and slope geometry only played a secondary role"; "there was a threshold rainfall intensity which would produce the global minimum factor of safety" |
| | 1533 | c' = 10 kPa, φ' = 26°, φ_b = 26°, γ = 20 kN/m³; Fredlund e Xing com C(ψ) = 1; k_w = k_s Θ^p (Leong e Rahardjo 1997b) |
| | 1534 | Tabela 1 (a, m = n = 1, θ_s = 0,45, p = 4); condições de contorno ("a boundary flux, q, equal to the desired rainfall intensity") |
| | 1535 | "limiting pore-water pressure of −75 kPa"; Bishop simplificado no SLOPE/W com pressões do SEEP/W |
| | 1539 | Ir(thr) = inflexão inferior; Fig. 7: A = 0,81 / 0,82 / 1,45; Ir(thr) = 140, 51, 11 mm/h = 0,39 / 1,42 / 3 k_s; taludes de k_s baixo: chuva escoa |
| | 1541 | "slopes with low ks … (α up to 63°) will not fail under short-duration rainfall (Tr ≤ 24 h)"; Tabela 7 e Fig. 9 |
| Tang et al. 2018 | 324 | SEEP/W + SLOPE/W; "FSmin, is extracted and adopted to reflect the most unstable state of the slope" |
| Pradhan et al. 2022 | 2 | "coupled hydro-mechanical analyses have been reported to provide a better prediction of factor of safety" |
| | 4 | análise totalmente acoplada no PLAXIS 2D (Biot); χ = S_e; "pore air pressure is assumed to be negligible"; SRM |
| | Tabela 2 | Mohr–Coulomb, c' = 3 kPa, φ' = 36°; ângulo de dilatância não listado |
| Liu e Wang 2021 | 1 | "often performed using a two-stage method … or a finite element method (FEM) with hydro-mechanical coupling" |
| | 4–5 | Drucker–Prager não associado; σ' = σ + δ S_w p_w (sem ar); ψ = 0 (Tabela 1) |
| | 14 | limiar de q/k_sat "about or slightly larger than 1 (e.g., Zhan and Ng, 2004; Rahardjo et al., 2007)" |
| Wang, Vardon e Hicks 2018 | 3, 4, 6 | χ = S; "gas pressure is omitted for simplicity"; ângulo de dilatância 0,1° (Tabela 2) |
| Sheng et al. 2008 | 811 | ar atmosférico "approximately true for most in situ conditions" |
| Fredlund 2006 | 292 | ar atmosférico "in most situations where the soil is exposed to atmospheric pressure and the air phase is continuous" |
| | 298 | a forma linear τ = c' + (σ_n − u_a) tan φ' + (u_a − u_w) tan φ_b "was published by Fredlund et al. (1978)" |
| | 312 | resistência cresce com φ' "for matric suctions up to the air entry value of the soil" |
| | 315 | "the uncoupled simultaneous solution yields satisfactory results for geotechnical engineering purposes" |
| Li et al. 1999 | 136 | ar passivo, p_g = 0, para simplificar |
| | 139 | Eqs. 33 e 36 (densidade do gás dentro do divergente; termo de gradiente de densidade) |
| Ghorbani et al. 2018 | 4 | comparações ar passivo × três fases "are rare" |
| | 11 | Eq. 35 com $(\dot{\mathbf w}^g/\rho_g)\cdot\nabla\rho_g$ |
| Khoei e Saeedmonir 2021 | 2 | ar passivo "lead to unrealistic outcomes" (cita Sheng 2008 e Callari 2009) |
| | 16 | Brooks–Corey para o gás, λ = 3 (Eq. 48) |
| Callari e Abati 2009 | 490 | p_g = 0 "acceptable in many situations … avoids problems … transition from full to partial saturation"; teste usado por vários autores |
| | 490–491 | Tabela 1 e Eq. 58 (propriedades do Liakopoulos); "in good agreement with available experimental data" |
| Pedroso 2015 | 791 | "it is common practice to import data from seepage analyses into stability codes" |
| | 792 | fase gasosa desprezada (p_g = 0) em todo o artigo |
| | Seção 3 | condições de contorno unilaterais, análogas à complementaridade da plasticidade |
| Bergamaschi e Putti 1999 | 1032 | estagnação de Picard ligada à "absence of the derivative terms in the Jacobian matrix" |
| Li et al. 2024 | 9 | FS do talude saturado "identical to that without considering suction" |
| Cheng et al. 2024 | 17 | FS mínimo constante porque o grau de saturação "keeps constant when the rainfall intensity exceeds a certain value" |
| Zhang et al. 2014 | Tabelas 4, 5 | FS (GLE) das oito leis de resistência |
| Wang, Idinger e Wu 2021 | 2913 | limiar ID proposto por Caine; a e b "empirical parameters defined by means of statistical analyses" |

## 3. Contradições encontradas na literatura

1. **Rahardjo et al. (2007), Tabela 7 × Fig. 9 (p. 1541).** Cinco das dez retas da Tabela 7 não reproduzem as
   retas da Fig. 9, digitalizada a 300 dpi com incerteza de ±0,02: f10 com 80 mm/h, f50 com 9 e 80 mm/h, f100 com 9
   e 80 mm/h. As linhas parecem trocadas. Por exemplo, a tabela dá FS_min = 1,455 − 0,013α para f100 com 80 mm/h,
   ou seja, 0,87 a 45°. A Fig. 9c mostra 1,46 e a Fig. 7c tem A = 1,45. O próprio texto (p. 1541) diz que esses
   taludes não rompem.
2. **Rahardjo et al. (2007), texto × Fig. 7c (p. 1539).** O texto dá I_r(ini) = 3,8 mm/h para o f100; a figura dá
   2,8 mm/h (0,778 k_s; V_r = 67,2 mm = 2,8 × 24).
3. **Rahardjo et al. (2007), lista de referências.** Fredlund e Xing (1994) aparece com as páginas 533–546, que são
   as de Fredlund, Xing e Huang (1994), segundo a lista de Fredlund (2006). Zhang et al. (2014), Tang et al. (2018)
   e Ghorbani et al. (2018) dão 521–532. O fascículo também diverge: 31(3) em Fredlund (2006) e Sheng et al. (2008),
   31(4) em Zhang et al. (2014) e Ghorbani et al. (2018).
4. **Liu e Wang (2021, p. 14) × Rahardjo et al. (2007, p. 1539).** Liu e Wang citam Rahardjo para um limiar de
   q/k_s "about or slightly larger than 1". Os valores de Rahardjo são 0,39; 1,42 e 3 k_s.
5. **Sheng et al. (2008, p. 811) × Fredlund (2006, p. 292).** Para Sheng, o ar atmosférico é "approximately true for
   most in situ conditions"; Fredlund exige solo exposto à atmosfera e fase ar contínua. Nos nossos resultados p_a
   chega a 36–40 kPa com a superfície alagada.
6. **Khoei e Saeedmonir (2021, p. 2) × Callari e Abati (2009, p. 491).** Para Khoei, o ar passivo leva "to unrealistic
   outcomes", citando Callari; Callari obtém "good agreement" com o ensaio. Nos nossos resultados as duas afirmações
   valem para processos diferentes: drenagem (Liakopoulos) e infiltração com a superfície alagada (talude).
7. **Pradhan et al. (2022, p. 2) × Fredlund (2006, p. 315)** sobre o acoplamento. Nos nossos resultados, o efeito do
   acoplamento sobre I_c é pequeno sem dilatância e grande com fluxo associado.
8. **Erros tipográficos dos originais (agentes de leitura).**
   - Li et al. (1999): Eq. 41 com ρ_w no lugar de ρ_g.
   - Ghorbani et al. (2018): Eq. 18 com S_wρ_w no termo do gás.
   - van Genuchten (1980): α = 0,079 no texto e 0,0079 na Tabela 1.
   - Cheng et al. (2024): valores divergentes entre as pp. 15–17 e 18.
   - Cheng et al. (2024) cita Ng e Shi como "1988".

## 4. O que não pôde ser verificado

- **Perfis da seção 6** (sucção 38 × 24 kPa; p_w −6 × −24 kPa; p_a média de 32 kPa). Vêm dos instantâneos das
  rodadas `fase_ar`, que estão só na máquina local e não no repositório. A figura `ar_aprisionado.pdf` foi mantida e
  a limitação está declarada no artigo.
- **Referências citadas só como atribuição de método**, conferidas apenas pela lista de um trabalho anexado: Bishop
  1959, Caine 1980, Guzzetti et al. 2008, Fredlund et al. 1978, Fredlund e Xing 1994, Griffiths e Lane 1999, Leong
  e Rahardjo 1997, Liakopoulos 1965. Devloo 1997 foi conferida no README da NeoPZ.
- **Crossref.** A consulta a api.crossref.org foi bloqueada pela política de rede do ambiente (403). Para liberar,
  adicione o domínio em *Allowed domains* nas configurações do ambiente:
  https://code.claude.com/docs/en/cloud-environments#network-access
