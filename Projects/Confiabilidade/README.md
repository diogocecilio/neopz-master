# Confiabilidade: FORM, Monte Carlo e amostragem por importância

Este projeto revisa um código Mathematica de confiabilidade e propõe um segundo problema. O documento completo é
`doc/confiabilidade.pdf`.

O código original calcula, para g = X1 X2 − X3 com três normais:
- o FORM, pelo `NMinimize`;
- Monte Carlo com N = 1000;
- a amostragem por importância no ponto de projeto.

O problema proposto é um talude infinito não saturado sob chuva, com os parâmetros do talude de Vargas Ceron et al.
(2025). O documento também mostra como aplicar os métodos ao talude de elementos finitos com campos aleatórios
(`SlopeMonteCarlo`).

## Arquivos

| arquivo | o que é |
|---|---|
| `confiabilidade.wl` | o código Mathematica revisado, para abrir e usar *Run All Code*. Comentários só em ASCII |
| `python/confiab.py` | o mesmo algoritmo em Python: Nataf, HL-RF com Armijo, Monte Carlo e importância |
| `python/problemas.py` | os dois problemas. Gera `doc/tabelas`, `doc/figuras` e `dados/p2_fragilidade.csv` (cerca de 2 min) |
| `mathics/roda_mathics.sh` | o teste do `.wl` no Mathics, um interpretador livre da linguagem (`pip install mathics3`). Usa substitutos das funções estatísticas e leva cerca de 1 min |
| `doc/confiabilidade.tex` | o documento: `cd doc && latexmk -pdf confiabilidade.tex` |

## O que foi corrigido no código original

A lista completa está na Tabela 1 do PDF.

- **`FromXtoY`.** As matrizes L e L⁻¹ estão trocadas: o código faz y = L z, e o certo é y = L⁻¹ z. Com variáveis
  independentes isso não aparece. No problema 2, com ρ(c', φ') = −0,5, o β cai de 3,54 para 2,51, e o P_f fica 30
  vezes maior.
- **A volta z → x.** O código faz x = Jxz.(z + Mneq); o certo é x = Jxz.z + Mneq. Por isso a saída do caderno é
  {160.+X1, ...}.
- **`FromYtoZ`.** Os nomes das matrizes estão trocados. O teste de Cholesky usa um determinante, que não mede o erro.
  Com uma única variável, a função devolve escalares em vez de matrizes.
- **Correlação.** Com marginais lognormais, a correlação no espaço normal tem de ser a de Nataf, não a das variáveis
  físicas.
- **Monte Carlo com N = 1000.** Para P_f ≈ 1,2e-3, o CoV é de 92 %, e 31 % das rodadas dão P_f = 0 (β = ∞).
- **Amostragem por importância.** O ponto de projeto foi digitado à mão, e falta o CoV.
- **Reprodutibilidade.** Falta o `SeedRandom`, então as rodadas não se repetem.

## Resultados

- **Problema 1.**
  - Referência por integração: P_f = 1,177e-3.
  - FORM: β = 3,04907 e P_f = 1,148e-3, 2,5 % abaixo da referência, em 6 iterações.
  - Amostragem por importância com 1000 amostras: CoV de 6 %. Para o mesmo CoV, Monte Carlo precisaria de cerca de
    250 vezes mais amostras.
- **Problema 2, antes da chuva** (μ_h = −2 m): P_f = 1,53e-4.
  - O FORM erra +30 %.
  - A amostragem por importância com 1000 amostras chega a um CoV de 7 %.
  - Monte Carlo com 1000 amostras dá 0.
- **Problema 2, depois da chuva** (μ_h = −0,5 m): P_f = 0,121, da ordem dos 11,5 % do artigo. Monte Carlo basta.
- **Curva de fragilidade P_f(μ_h).** Vai de 1e-12 a 0,6.
  - A carga de pressão domina: α² ≈ 0,65 a 0,69.
  - O HL-RF sem a busca de Armijo não converge para μ_h ≤ −2,5 m.

## No talude com campos aleatórios

- **O espaço U.** É formado pelos ξ da expansão de Karhunen-Loève, que já são normais padrão independentes, e
  G = Γ − 1.
- **Monte Carlo.** Com P_f ≈ 11 %, ele basta.
- **Amostragem por importância.** Compensa abaixo de P_f ≈ 1e-2, com o deslocamento igual à média dos ξ das
  realizações que falham, restrito aos primeiros modos.
  - A mudança no código fica em `RandomFields::Sample` (`Projects/Common/SlopeFields.cpp`): somar o deslocamento aos ξ
    e gravar ln w.
- **FORM direto.** Não é viável: o gradiente custaria 2n análises por iteração, com n = 6000.
