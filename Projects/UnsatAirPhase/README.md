# UnsatAirPhase: teste 3 (água e ar, esqueleto rígido)

As duas equações de massa juntas, com todos os deslocamentos prescritos (esqueleto rígido, Δε_v = 0). Documento:
Seções 15, 16, 17 e 22.

| teste | critério | resultado |
|---|---|---|
| 3.0 Jacobiano dos blocos água–ar (Taylor), com e sem condensação | ordem mínima em (1,9; 2,1) | mediana 2,00003, mínima 1,99972 |
| 3a ar confinado (lei de Boyle): injeção de água igual a 30% do volume de ar em 10 h; massa M = Σ m_j n S_a P | V_a/V_a0 < 0,75; p_a > 30 kPa; max \|ΔM\|/M₀ < 10⁻⁴ | V_a/V_a0 = 0,70; p_a até 43,5 kPa; 8,1·10⁻¹⁶ |
| 3b limite de ar passivo: mobilidade do ar ×1, ×10, ×100 contra a solução com p_a = 0 (γ_a = 0) | diferenças decrescentes, a última < 1% da faixa de p_w | 0,0247; 0,00248; 0,000248 kPa (faixa 23,8 kPa) |
| 3c ar livre no topo: balanço de ar e de água | água < 10⁻⁶, ar < 10⁻³ | 1,9·10⁻¹⁴ e 5,4·10⁻¹⁶ (35% do ar sai) |

81 verificações, cerca de 3 s. A conservação exata da massa de ar vem do termo n S_a^n Δp_a / P (decisão 2) com a
condensação nos vértices (documento, Seção 15).
