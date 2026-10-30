# UnsatCoupled: testes 0 e 4 (modelo completo)

Documento: Seções 16, 17, 18 e 22.

| teste | critério | resultado |
|---|---|---|
| 0 Jacobiano completo (u, p_w, p_a) com plasticidade, no teste de Taylor, em quatro variantes: χ = S_w; χ constante; χ = S_w com peso da mistura, grãos compressíveis (η > 0) e armazenamento sem condensação; saturado (a_VG = 10⁻⁸) | ordem mediana em (1,9; 2,1); pontos plásticos > 0 | 2,000 nas quatro; 96 pontos plásticos |
| 4a Terzaghi (limite saturado) contra a série | max \|p/p₀ − ref\| < 0,01 | 0,0028; 0,0014; 0,0007; 0,0008 em T = 0,05; 0,1; 0,2; 0,5 |
| 4b Mandel (placa rígida lisa) contra a solução analítica | erro no centro < 0,02 p₀; pico de Mandel–Cryer > 1,02 do valor inicial | 0,0008 a 0,0038 em ct/a² = 0,01 a 0,6; pico 5,375 contra 5,001 kPa |
| 4c Liakopoulos com ar passivo contra os dados experimentais digitalizados de Callari e Abati (2009, Fig. 3) | RMS até 1,5 vez o da solução numérica de Callari e Abati | vazão 0,225·10⁻⁶ m/s (Callari 0,227); p_w 0,44 kPa (Callari 0,48) |
| 4d Liakopoulos com ar ativo (p_a = 0 no topo e na base, `[.lento]`) | sem critério numérico; compara com o ar passivo | vazão na base (10⁻⁶ m/s) em 300, 1200, 3600 e 7200 s: passivo 3,38; 2,01; 0,88; 0,34; ativo 2,51; 0,76; 0,41; 0,26 |

1869 verificações nos testes rápidos, cerca de 3 min. Dados em `data/`: vazão e p_w experimentais e a solução
numérica de Callari e Abati, lidos das figuras.

Mandel na base hierárquica do NeoPZ: o deslocamento uniforme da placa é imposto com o valor nas funções de vértice e
zero nas de aresta, e a força total é a soma dos resíduos só nas funções de vértice.

O coeficiente de adensamento usa o módulo confinado λ + 2G. O coeficiente do apêndice A.1 de Haghighat et al. (2021,
arXiv:2110.03049v1, p. 25) usa K_dr, definido na p. 5 como módulo volumétrico drenado; com ele o erro no centro chega a
0,23 p₀ em ct/a² = 0,6 (o teste imprime os dois).
