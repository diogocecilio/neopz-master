# UnsatWaterPhase: teste 2 (água, ar passivo, esqueleto rígido)

Equação de Richards dentro do modelo completo: todos os deslocamentos prescritos (esqueleto rígido, Δε_v = 0) e
p_a = 0 prescrito em todo o domínio. Documento: Seções 14, 16, 17 e 22.

| teste | critério | resultado |
|---|---|---|
| 2.0 Jacobiano dos blocos da água (Taylor), com e sem condensação | ordem mínima em (1,9; 2,1) | 1,99976 e 1,9997 |
| 2.1 perfil hidrostático sem fluxo | variação < 10⁻⁸ kPa | 1,8·10⁻¹⁵ kPa |
| 2.2 infiltração permanente em solo de Gardner (a_G = 0,2 kPa⁻¹, k_rw = S_e, q = 0,3 k_s, H = 3 m) contra a solução fechada (Eq. 22.1 do documento; Li et al. 2024, Eq. 9, com Q = −q) | erro < 1% da sucção máxima; razão dos erros ao dobrar a malha > 3 | 0,0050 kPa (ny = 30) e 0,0015 kPa (ny = 60); sucção máxima 5,99 kPa |
| 2.3 balanço de massa na infiltração transiente (forma secante, condensada) | erro relativo < 10⁻⁶ | 2,2·10⁻¹⁴ |
| 2.4 infiltração de Philip (Callari e Abati 2009, Seção 4.2: coluna de 0,30 m, −100 kPa inicial, −7,5 kPa no topo, Tabela 2, 1000 passos de 1 s e depois 10 s) contra a solução semianalítica em t = 21 600 s digitalizada da Fig. 7 (`data/philip_callari_fig7.csv`): elevação em que p_w = −30, −50, −70, −90 kPa | distância < 5 mm (um elemento da malha fina de Callari e Abati) | 2,0 mm no pior ponto; 1,4 mm em −50 kPa |
| 2.4 (`[.lento]`) a mesma coluna com 120 elementos e Δt/2 | mudança < 5 mm | 1,8 mm |

3137 verificações nos testes rápidos, cerca de 3 min (a coluna de Philip tem 3060 passos); o lento leva cerca de 16 min. A lei de retenção de Gardner (`EGardner`) foi acrescentada a
`TPZUnsaturatedHydraulics` para o teste 2.2.

Observação: Li et al. (2024, p. 6) dizem que Q > 0 na infiltração, mas com Q > 0 a Eq. 9 dá sucção acima da
hidrostática; com os dados do teste 2.2 o argumento do logaritmo fica negativo. A solução reproduzida pelo modelo, com a
vazão entrando pelo topo, é a Eq. 9 com Q = −q.
