# UnsatMaterialPoint: ponto material

Verifica, sem malha, as partes do material que entram em cada ponto de integração. Documento: Seções 3, 10 e 22.

| teste | critério | resultado |
|---|---|---|
| leis hidráulicas: Γ = −dS_w/dp_c e dk_r/dS_w contra diferenças centrais (todas as leis) | relativo 10⁻⁶ e 10⁻⁵ | passou |
| van Genuchten e Mualem (água e ar) contra as fórmulas das notas (Seções 3 e 5) | 10⁻¹² a 10⁻¹⁴ | passou |
| Bishop com χ = S_w (padrão): p_B, χ_w = S_w − Γ p_c, χ_a = S_a + Γ p_c, diferenças centrais | 10⁻¹⁴, 10⁻¹², 10⁻⁶ | passou |
| PV2: tensões e Jacobiano contra as Tabelas 4–6 e Eqs. 57–58 do artigo; região escolhida contra a projeção de ponto mais próximo por enumeração KKT (20 000 estados) | < 10⁻¹²; mais de 1000 estados em cada região | passou |
| PV2: tangente consistente no teste de Taylor nas cinco regiões | ordem em (1,85; 2,15); ápice e elástico exatos; \|D − Dᵀ\|/\|D\| < 10⁻¹² | passou |
| redução de resistência: c/F e atan(tan φ / F), com ψ = φ mantido | 10⁻¹⁴ | passou |

408 verificações, menos de 1 s. Rodar: `./UnsatMaterialPoint` (ou um teste: `./UnsatMaterialPoint "[mc]"`).
