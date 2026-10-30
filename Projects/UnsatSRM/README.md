# UnsatSRM: teste 1a (redução de resistência, talude seco)

Talude de 10 m a 45°, c = 10 kPa, φ = 30°, γ = 20 kN/m³, E = 2·10⁴ kPa, ν = 0,3, p_w = p_a = 0 congeladas.
Estado geostático com F₀ = 0,6, depois c/F e tan φ / F (ψ = φ) até o colapso, com bisseção (tolerância 0,002).
Malha estruturada com nf colunas na zona da face (`SlopeMeshFace`). Documento: Seção 20.

| malha | equações | FS_SRM | diferença para Bishop (1,206) |
|---|---|---|---|
| nf = 20 | 3 618 | 1,2535 | +3,9% |
| nf = 40 | 13 210 | 1,2250 | +1,6% |
| nf = 80 (`[.lento]`) | 52 422 | 1,2109 | +0,4% |

Testes: FS diminui de nf = 20 para nf = 40 e fica a 2% de Bishop; FS independe de F₀ ∈ {0,5; 0,6; 0,8} (1,25352 nos
três); lento: nf = 80 a 1% de Bishop. Referências (README do projeto `SlopeMohrCoulomb`, branch
`claude/relaxed-volta-1r4gog`): Bishop simplificado 1,206; modelo mecânico com refino adaptativo, ciclo 5, ν = 0,3:
1,207. Cerca de 4 min (testes rápidos); o lento leva cerca de 30 min.
