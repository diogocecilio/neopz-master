# UnsatGIM: teste 1b (aumento de gravidade, talude seco)

Mesmo talude do `UnsatSRM`. A partir do estado sem tensões e com a resistência real, o fator λ multiplica o peso e a
parte positiva de p_B; a sucção não é escalada. Bisseção até a tolerância 0,002. Documento: Seção 20.

| malha | equações | FS_GIM | diferença para Bishop (1,806) |
|---|---|---|---|
| nf = 20 | 3 618 | 2,0349 | +12,7% |
| nf = 40 | 13 210 | 1,8743 | +3,8% |
| nf = 80 (`[.lento]`) | 52 422 | 1,8146 | +0,5% |

Testes: FS diminui de nf = 20 para nf = 40 e fica a 5% de Bishop; lento: nf = 80 a 2% de Bishop; com nf = 20,
GIM ≥ SRM (2,031 contra 1,250) e a sucção uniforme de 20 kPa (p_w = −20 kPa, p_a = 0, não escalada) aumenta o fator
para 3,086. Para Mohr–Coulomb, aumentar o peso equivale a reduzir só a coesão, por isso FS_GIM ≥ FS_SRM (documento,
Eq. 20.3). Referências (README do projeto `SlopeMohrCoulomb`, branch `claude/relaxed-volta-1r4gog`): Bishop
simplificado 1,806; modelo mecânico com refino adaptativo, ciclo 5, ν = 0,3: 1,789. Cerca de 6 min (testes rápidos); o
lento leva cerca de 30 min.
