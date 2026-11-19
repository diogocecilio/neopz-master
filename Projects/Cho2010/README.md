# Cho2010

O Exemplo 2 de Cho (2010) no NeoPZ, para o artigo da amostragem por importância (AI):

> S. E. Cho (2010). Probabilistic assessment of slope stability that considers the spatial variability of soil
> properties. *J. Geotech. Geoenviron. Eng.* 136(7):975–984. doi:10.1061/(ASCE)GT.1943-5606.0000309

## O problema

- **Talude** (Fig. 5 do artigo): seco, 1:1, com H = 10 m, 5 m de fundação, 10 m de terreno antes do pé e 10 m de crista.
  Aqui ele está espelhado, com a crista à esquerda.
- **Solo** (Tabela 4): c' log-normal com média 10 kPa e CoV 0,3; φ' log-normal com média 30° e CoV 0,2; γ = 20 kN/m³.
  A correlação cruzada r(c, φ) entre os campos gaussianos vai de −0,7 a 0,5, com base em −0,5.
- **Campos:** correlação exponencial separável, com l_h = 20 m e l_v = 2 m.
- **Referência:** o artigo usa equilíbrio limite (Bishop, com busca da superfície crítica) e Monte Carlo com 50 000
  realizações:
  - FS determinístico = 1,204;
  - P_f = 0,0171 com r = −0,5 (Tabela 5);
  - P_f ≈ 0,004 com r = −0,7 (Fig. 9).

## O modelo

- **Material:** o Mohr–Coulomb elastoplástico dos outros programas (`TPZPlasticStepVoigt<TPZYCMohrCoulombPV2>` dentro do
  `TPZMatUnsaturatedUPwPa`), com as pressões congeladas em zero. O solo é seco, e a tensão efetiva é a total. O fluxo é
  associado, com E = 20 MPa e ν = 0,3.
- **Fator:** pelo SRM, que reduz c' e tan φ', como o FS do equilíbrio limite. Usa arc-length em F e refinamento
  adaptativo em h, refeito do zero em cada realização.
- **Campos:** Karhunen–Loève analítica (`RandomField`), com 2000 termos. A correlação cruzada é
  `RandomField::SetCrossCorrelation`.
- **Amostragem por importância:** o deslocamento sai do piloto pela entropia cruzada, na direção da banda
  (`ImportanceSampling.h`).

O programa usa as classes do `SlopeRandomFieldRefinement`; as opções estão no cabeçalho do `cho2010.cpp`.

## O determinístico (SRM)

| malha do fator | FS | contra o Cho (1,204) |
|---|---|---|
| tri1 + 4 níveis | 1,2097 | +0,5% |
| tri2 + 3 níveis | 1,2101 | +0,5% |
| tri2 + 4 níveis | 1,2047 | +0,06% |

O GIM dá Γ = 1,79: num solo c–φ ele é outro fator, mas a P_f é a mesma. Os casos usam tri2 + 4 níveis, com ~30 s
por realização em 1 thread.

## Rodar no PC

A pasta `~/projects/runs/cho2010` tem o `casos.txt` e o `rodar.sh`. Ele usa a fila do
`SlopeRandomFieldRefinement/cluster/tarefas.sh` com 14 processos e, no fim, grava a comparação em `resultado.txt`:

| caso | N | o que é |
|---|---|---|
| `r05_direto` | 4000 | Monte Carlo direto com r = −0,5: o piloto (as 1000 primeiras) e a referência |
| `r05_ai` | 1000 | amostragem por importância na direção da banda, com outra semente |
| `r07_direto` | 1000 | piloto com r = −0,7 (P_f ~0,4%): o conjunto elite são as 10% de menor FS |
| `r07_ai` | 1000 | amostragem por importância com r = −0,7 |

Um caso sozinho:

```bash
Cho2010 n=0 det=1 figuras=-1 fator=tri2 nref=4      # o determinístico
Cho2010 n=50 rcphi=-0.5 fator=tri2 nref=4 saida=mc  # 50 realizações
```

## A base de Nyström com 2000 modos

O `nystrom_base.py` grava a base que a opção `nystrom=` lê. Para milhares de nós e de modos, os autopares saem do
Mathematica (MKL, em paralelo), porque o `eigh` do Python sem threads leva horas:

```bash
python3 -I nystrom_base.py base_nys2000.txt M=2000 dx=0.2 dy=0.1 nos=base_nys2000.nos
wolframscript -file nystrom_autopares.wl base_nys2000.nos 2000 20 2 0.02 base_nys2000
python3 -I nystrom_base.py base_nys2000.txt M=2000 dx=0.2 dy=0.1 autopares=base_nys2000 copias=100
```

A grade de 0,2 × 0,1 m tem 15 000 nós e retém 99,7% da variância com 2000 modos. A cópia `base_nys2000_M100.txt` tem os
100 primeiros modos, com os mesmos autovetores, para o FORM. A grade antiga, de 0,5 × 0,25 m, erra os autovalores
acima de ~50 modos.

## ISMC com parada pelo CoV (etapa 3b do plano)

A campanha `etapa3b` do cluster (`casos_etapa3b.txt`) roda a amostragem por importância no ponto de projeto de seis
FORM (M = 20 e 100 por diferenças centrais e M = 20 pelo GIM modal, r = −0,5 e −0,7) com os 2000 modos de Nyström.
A opção de fila `cov_alvo=0.05` do `tarefas.sh` para cada caso quando o CoV da P_f chega a 5%. O `ismc_parada.py` dá,
por semente, o N em que o CoV cruza o alvo e a P_f nesse ponto e no fim, e, por centro, a média das sementes, o z contra
os diretos (`r05_direto` da etapa1 e `r07_direto` da etapa3), o custo com o FORM e o ganho:

```bash
python3 -I ismc_parada.py ~/projects/runs/fermi_etapa3b/etapa3b ../SlopeRandomFieldRefinement/cluster/casos_etapa3b.txt figuras=figs
```
