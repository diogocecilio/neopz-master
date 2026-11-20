# Cho (2010) com a análise limite do Ceron: base de comparação análise limite × MEF

Refaz os casos das Tabelas 4 e 5 do `Projects/PlanoDeTrabalho/comparacao_metodos/relatorio_comparacao.pdf` (Cho 2010,
Exemplo 2: talude c-φ seco, H = 10 m a 45°, c' e φ' lognormais com r = −0,5 e −0,7) com a **análise limite cinemática**
do Ceron (2025) no lugar do MEF (SRM, tri2 + 4 níveis), **nos mesmos campos**. A realização n daqui é a realização n do
MEF com a mesma semente, o mesmo campo e o mesmo centro. Assim, a comparação é feita realização a realização.

## Arquivos

| arquivo | o que faz |
|---|---|
| `campos_cho.py` | c' e φ' do Cho2010 em Python, iguais aos do C++ dados os ξ. Contém: o gerador do `RandomField::Sample` (xoshiro256** + Box-Muller); a KL analítica do retângulo 30 × 15 (forma tensorial na grade); a base de Nyström de 2000 modos (Ψ = C·W·F·Λ^−½ pré-calculado em `.npy` float32); a correlação cruzada r; a lognormal. `python3 campos_cho.py cache base=… cache=…` só pré-calcula Ψ. |
| `la.py` | Avaliador da análise limite no talude do Cho. Γ pela busca rápida, que repete a busca forte quando Γ < 1,5. FS pelas Eqs. 57-58, com c/F e atan(tan φ/F) nó a nó, por secante. Modo de ruptura (pé, face ou base), como na Fig. 20 do Ceron. |
| `cinematico.py`, `rapido.py` | **Cópias congeladas** de `../cinematico.py` (md5 a147c234…) e `../rapido.py` (md5 b33db758…), de 20/11/2026 12:42. O `rapido.py` de cima está sendo mexido por outro trabalho; esta pasta não depende dele. |
| `testes.py` | `campos`: c e φ contra o `campos.csv` do Cho2010. `det`: Γ = 1,777 e FS = 1,203. `busca`: rápido × `cinematico.minimizar` exato, em 30 realizações. |
| `mc_la.py` | MC direto e ISMC (com `centro=`) pela análise limite. Usa multiprocessing, é retomável (apaga uma última linha cortada antes de continuar), para pelo CoV e grava um CSV por realização, na ordem em que terminam. A coluna `fs_conv` diz se a secante do FS convergiu. `resumo=1` imprime P_f, IC, Γ médio, os modos e o tempo. |
| `form_la.py` | FORM (HL-RF) com g = ln Γ_LA e gradiente por diferenças finitas, com M termos por campo. Grava `form.json` e `ustar.txt`. |
| `ls_la.py` | Line sampling pela análise limite numa direção de FORM (do MEF ou da análise limite). Usa os mesmos z do `Projects/Cho2010/line_sampling.py`, com a mesma semente e o mesmo K. Uma linha que não cerca a raiz em 25 análises vai para `perdidas.csv` e fica fora da P_f, como no MEF. |
| `pareado.py` | `gerar`: grava os arquivos xi (e o h em `parametros.json`). `la`: Γ e FS da análise limite. `mef`: roda o Cho2010 com `xi=` e o mesmo fator das referências; sai com código 3 se alguma realização ficar sem FS. `juntar`: compara LA × MEF (também Γ_LA/FS_MEF); junta também com os CSV do MC direto, das ISMC ou do direto pelo GIM do MEF já rodados. |
| `estima_pf.py` | Prevê a P_f da análise limite e o N para CoV de 5 %. Ajusta ln FS_LA × ln FS_MEF no piloto pareado e aplica o ajuste a todo o MC direto do MEF. |
| `resumo.py` | A tabela LA × MEF de uma pasta de saídas (o fermi ou `resultados/`). |
| `form_mef/` | Os FORM do MEF (u*, β), copiados de `~/projects/runs`: `form_r0x_{M20,M100,gim20}` (fermi_etapa3b) e `form_r0x_gim2000` (etapa3b_ls). São os centros das ISMC e as direções do LS. |
| `testes_mef/` | Saídas do Cho2010 com `figuras=`, usadas no teste dos campos (comandos abaixo). |
| `resultados/` | Saídas dos testes e dos pilotos: `piloto_h0.025/` (o atual) e `piloto/` (o antigo, com h = 0,1). |
| `../../SlopeRandomFieldRefinement/cluster/la_cho.sbatch` | O job do fermi (ver abaixo). |

## Convenções

- **Coordenadas:** as do Ceron, com origem na aresta da crista, x para fora e y para baixo. A conversão é X_cho = x + 10
  e Y_cho = 15 − y.
  - O talude é `kin.Talude(10, 45°, −10, 20, 15)`, o domínio [0, 30] × [0, 15] do Cho2010.
  - Solo seco: u = 0 e γ = 20 kN/m³.
- **Grade:** passo **h = 0,025 m** (`la.H_PADRAO`; `h=` em todos os drivers e `H=` no sbatch).
  - O φ_i (máximo de φ ao longo do arco) lê o campo P1 da grade. O campo de 2000 modos com núcleo exponencial é áspero
    na escala de 0,1 m, e com h = 0,1 o máximo sai baixo: Γ ~3,5 % e FS ~1,2 % abaixo (mediana nas 200 realizações do
    piloto r05).
  - Estudo de h (semente 777, reais 3 e 8, r05):

    | h | Γ real 3 | Γ real 8 | FS real 3 | FS real 8 | tempo (com FS) |
    |---|---|---|---|---|---|
    | 0,1 | 0,9478 | 1,2198 | 0,9805 | 1,0795 | 8-10 s |
    | 0,05 | 0,9623 | 1,2433 | 0,9860 | 1,0870 | 8-12 s |
    | 0,025 | 0,9676 | 1,2543 | 0,9881 | 1,0900 | 10-12 s |
    | 0,0125 | 0,9693 | 1,2574 | 0,9885 | 1,0914 | 13-14 s |

    Com h = 0,025 o FS fica a 0,05-0,25 % e Γ a 0,3 % do extrapolado. Γ sozinho custa 2,6 s (h = 0,025) contra 1,5 s
    (h = 0,1), nas mesmas 48 realizações.
  - O MEF também tem a sua resolução: o campo é lido nos pontos de Gauss da malha tri2 + 4 níveis.
- **Sem mecanismo:** em 1 de 200 realizações do r07 (real 144, φ até 71°) a busca não acha nenhum mecanismo admissível
  com h = 0,025 (as espirais da grade de partida não chegam à superfície). Então Γ = inf, modo `nenhum`, e o FS é
  calculado assim mesmo (deu 1,38): a falha passa a ser decidida pelo FS. Qualquer outro erro numa realização vai para
  o log e a realização fica fora do CSV, sem derrubar o job.
- **Falha:** Γ < 1, que equivale a FS < 1. A P_f da análise limite compara direto com a do MEF.
  - A análise limite é um limite superior com um bloco rígido. Por isso espera-se Γ_LA e FS_LA acima do colapso
    verdadeiro, e P_f(LA) < P_f(MEF).
- **Campos:**
  - MC direto: `campo=analitico`, semente 2010. É o mesmo do MC direto do MEF (24 000 com r05; 110 000 com r07).
  - ISMC e LS: `campo=nystrom`, com 2000 modos e 99,73 % da variância. É o mesmo das ISMC (sementes 21-25) e do LS do
    MEF (semente 7).

## Testes (no PC, 20/11/2026)

```bash
cd ~/projects/neopz-master/Projects/CeronCinematico/cho
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
nice python3 -I testes.py campos     # ~5 min na 1a vez (lê a base de Nyström de 525 MB e grava o .npz ao lado)
nice python3 -I testes.py det
nice python3 -I testes.py busca n=30 procs=3
```

**Campos.** Os campos de referência saem destes comandos, rodados em `testes_mef/` (EXE = o Cho2010 da build release):

```bash
$EXE n=2 primeira=0 semente=2010 rcphi=-0.5 fator=tri2 nref=4 srm_f0=0.8 metodo=srm figuras=0,1 saida=ana05
$EXE n=1 primeira=5 semente=2010 rcphi=-0.7 fator=tri2 nref=4 srm_f0=0.8 metodo=srm figuras=5 saida=ana07
$EXE n=1 primeira=0 semente=21 rcphi=-0.5 fator=tri2 nref=4 srm_f0=0.8 metodo=srm nystrom=<base_nys2000.txt> desloc_xi=<form_r05_M20/ustar.txt> figuras=0 saida=nys05
```

- Nas 4 realizações (1210 pontos cada), c e φ batem com o `campos.csv` do C++ até 4,9e-6, que é a precisão do CSV:
  - analítica r05, realizações 0 e 1;
  - analítica r07, realização 5;
  - Nyström com o centro da ISMC.
- A grade da análise limite bate com a fórmula ponto a ponto: 3e-15 na analítica e 2e-6 na de Nyström (Ψ em float32).
- Em 400 realizações, as estatísticas amostrais conferem: média de c ≈ 10, CoV 0,29-0,30, CoV de φ ≈ 0,20 e
  corr(ln c, ln φ) ≈ r.

**Determinístico.**

| | Γ | FS |
|---|---|---|
| exato | 1,7770 | — |
| rápido | 1,7770 (pelo pé) | 1,2032 |
| forma fechada | — | 1,2032 |
| artigo | 1,777 | 1,203 |

O Cho LEM dá 1,204 e o SRM do MEF 1,2044.

**Rápido × exato** (30 realizações, r05, semente 2010; `resultados/teste_busca.txt`):
- Onde Γ < 1,5 (7 realizações, onde se decide a falha), a diferença fica abaixo de 0,05 %, com nenhuma falha classificada
  diferente.
- Em 3 de 30, com Γ de 1,8 a 2,9, a busca rápida acha o mecanismo pelo pé e não o pela face, ficando +2,5 a +6,6 % acima.
  A busca forte acerta essas três, mas custa 2,7 vezes mais, e por isso só roda quando Γ < 1,5.
  - O `pareado.py la` usa a busca forte em todas, porque ali a comparação é valor a valor.
- Tempo: 0,9 a 1,2 s por realização, contra 11 s do exato.

## Comandos no PC

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
# MC direto LA (piloto): realizações 0..199 da semente 2010, com FS em todas
nice python3 -I mc_la.py saida=resultados/piloto_h0.025/direto_r05.csv r=-0.5 n=200 cov=0 procs=4 fs_ate=200
nice python3 -I mc_la.py saida=resultados/piloto_h0.025/direto_r05.csv resumo=1
# pareado com o MC direto do MEF já feito (as mesmas realizações)
python3 -I pareado.py juntar la=resultados/piloto_h0.025/direto_r05.csv mef=$HOME/projects/runs/fermi_etapa1/etapa1/r05_direto
# previsão da P_f LA e do N para CoV 5 %
python3 -I estima_pf.py la=resultados/piloto_h0.025/direto_r05.csv mef=$HOME/projects/runs/fermi_etapa1/etapa1/r05_direto
# pareado novo pelo xi= (LA e MEF nas mesmas realizações)
python3 -I pareado.py gerar pasta=resultados/piloto_h0.025/pareado20 r=-0.5 semente=2010 n=20
nice python3 -I pareado.py la pasta=resultados/piloto_h0.025/pareado20 procs=2
nice python3 -I pareado.py mef pasta=resultados/piloto_h0.025/pareado20 exe=$HOME/projects/neopz-master-build-release/Projects/Cho2010/Cho2010 procs=2
python3 -I pareado.py juntar pasta=resultados/piloto_h0.025/pareado20
# ISMC LA no centro do FORM M20 do MEF (campo de Nyström), pareada com a ISMC do MEF da semente 21
nice python3 -I mc_la.py saida=resultados/ismc_r05_mef_M20_s21.csv r=-0.5 campo=nystrom semente=21 n=200 cov=0 procs=4 \
    centro=form_mef/form_r05_M20/form.json base=$HOME/projects/runs/cho2010/nystrom2000/base_nys2000.txt cache=cache/psi_nys2000_h0.025.npy
python3 -I pareado.py juntar la=resultados/ismc_r05_mef_M20_s21.csv mef=$HOME/projects/runs/fermi_etapa3b/etapa3b/r05_nys_M20_s21
# FORM LA e LS LA
nice python3 -I form_la.py pasta=resultados/form/r05_M20 r=-0.5 M=20 procs=4 base=… cache=cache/psi_nys2000_h0.025.npy
nice python3 -I ls_la.py direcao=form_mef/form_r05_gim2000/form.json pasta=resultados/ls/r05_mef_gim2000 r=-0.5 procs=4 base=… cache=…
```

## Fermi (`la_cho.sbatch`)

Cada etapa é um job de um nó com 24 processos. As saídas vão para `$DATA/cho_la/h0.025/` (`H=` muda o passo e a
pasta), e todo job termina com o `resumo.py`. O sbatch faz `source modulos.sh` (o MKL do `bin/Cho2010`), para com
código ≠ 0 se um passo falhar e não esconde os logs.

| etapa | o que roda | estimativa no fermi |
|---|---|---|
| `ETAPA=checa` | Python com numpy/scipy/matplotlib; `ldd bin/Cho2010`; o determinístico; Ψ de Nyström com h = 0,025 (3,9 GB, uma vez, em 24 processos); o tempo de 48 realizações com FS. Rode primeiro. | ~0,5-1 h |
| `ETAPA=direto R=r05` | MC direto LA até CoV 5 %, no mínimo 24 000 (todo o direto do MEF). FS nas 24 000 primeiras e onde Γ < 1,3. Junta com `$DATA/etapa1/r05_direto` e, se existir, com `$DATA/etapa3b_gim/r05_direto_gim`. | ~54 000 reais, ~13-16 h |
| `ETAPA=direto R=r07` | O mesmo, no mínimo 110 000; FS nas 24 000 primeiras e onde Γ < 1,1. Junta com `$DATA/etapa3/r07_direto`. | ~370 000 reais, ~35-45 h (pode passar dos 2 dias: submeta de novo, ele retoma) |
| `ETAPA=direto R=r00` | r = 0, só Γ, até CoV 5 %: a conferência com os 5,50 % do Ceron. | ~7000 reais, ~0,5 h |
| `ETAPA=amostragem R=r05` ou `r07` | FORM LA M = 20 e 100 (`FORM_2000=1` acrescenta M = 2000); LS LA (100 linhas, `N_LS=`) nas direções mef_M20, mef_gim2000, la_M20, la_M100; ISMC LA nos centros mef_M20, mef_gim20, mef_M100 e la_M20, sementes 21-25 com 2000 cada (`N_ISMC=`), FS só na semente 21. | ~5-8 h cada |
| `ETAPA=pareado` | Gera os xi de `N_PAR` = 240 realizações (`R_PAR` = r05, `SEMENTE_PAR` = 2010). LA com a busca forte em todas; MEF com `bin/Cho2010 xi=`; junta. Com a semente 2010, confere que o FS pelo `xi=` é idêntico ao do MC direto do MEF. Sai com erro se o MEF falhar. | ~1 h |
| `ETAPA=resumo` | Só a tabela. | minutos |

- **Python:** o sbatch procura `$HOME/venv_la/bin/python3` e depois o `python3` do sistema, com numpy, scipy e
  matplotlib.
  - Se nenhum servir, para com a instrução (código 2). A preparação é feita uma vez, no nó de login:
    `python3 -m venv $HOME/venv_la && $HOME/venv_la/bin/pip install numpy scipy matplotlib`.
  - `PY=<python>` força outro python.
- **Fila:** todos os jobs têm o nome `laCho`; com `--dependency=singleton` rodam um de cada vez, deixando as outras
  vagas para as campanhas do MEF.

## Piloto no PC (20/11/2026; PC dividido com outros trabalhos, tempos inflados)

Tudo fica em `resultados/piloto_h0.025/`.

### MC direto LA, realizações 0..199 da semente 2010, com FS em todas

São as mesmas realizações do MC direto do MEF, pareadas uma a uma.

| r | falhas LA / MEF (nas 200) | P_f LA (piloto) | Γ médio | FS_LA − FS_MEF | corr(FS_LA, FS_MEF) | pelo pé |
|---|---|---|---|---|---|---|
| −0,5 | 3 / 5 | 1,50 % (CoV 57 %) | 1,989 | +0,0226 ± 0,0073 (min +0,011) | 0,9975 | 87 % |
| −0,7 | 1 / 1 | 0,50 % (CoV 100 %) | 1,988 | +0,0216 ± 0,0057 (min +0,012) | 0,9978 | 92,5 % |

- FS_LA > FS_MEF em todas as 400 realizações, como se espera de um limite superior. Nenhuma realização falhou só na LA.
- A previsão da P_f da LA (`estima_pf.py`) ajusta ln FS_LA = a + b ln FS_MEF no piloto (resíduo com dp 0,005) e aplica o
  ajuste a todo o MC direto do MEF:

  | r | P_f LA prevista | P_f MEF | razão | N para CoV 5 % |
  |---|---|---|---|---|
  | −0,5 | 1,01 % | 1,700 % | 0,60 | ~39 000 |
  | −0,7 | 0,17 % | 0,366 % | 0,45 | ~240 000 |

### Pareado novo pelo `xi=`, 20 realizações (`pareado20/`)

- O MEF local, pelos xi gravados em Python, dá FS **idêntico** (diferença 0) ao do MC direto do fermi.
- FS_LA − FS_MEF = +0,0205 ± 0,0051, com corr 0,9992.
- A LA com a busca forte em todas difere no máximo 0,47 % em Γ e 3e-4 em FS da política padrão.

### Outros testes do piloto

- **ISMC LA**, centro FORM M20 do MEF, semente 21, 40 realizações: o log_peso é igual ao do MEF (5e-10).
  - corr(FS_LA, FS_MEF) = 0,998 e FS_LA − FS_MEF = +0,0135.
  - Falhas: 15 na LA e 16 no MEF.
- **LS LA** na direção gim2000 do MEF, 8 linhas: c_i(LA) − c_i(MEF) fica entre +0,11 e +0,21 nas mesmas linhas.
  - Isso dá P_f(LA) ≈ 1,14 % nessas 8 linhas, com 4,9 análises por linha.
- **FORM LA M = 20**, r05, Nyström: β = 1,9907, contra 1,9775 do MEF. cos(u_LA, u_MEF) = 0,9988.
  - Convergiu em 4 iterações, com 324 minimizações e 140 s com 4 processos.

### Custo por realização (1 processo)

| o que | tempo no PC |
|---|---|
| campos | 0,04 s |
| Γ (busca rápida, e a forte se Γ < 1,5) | 1,5 a 2 s |
| FS | +5 a 6 s |
| Ψ de Nyström (uma vez) | 8 min com 1 thread |
| MEF (SRM tri2 + 4) | 40 a 60 s no PC; 69 s no fermi |

MC direto LA até CoV 5 %, só com Γ:

| r | CPU | PC com 4 processos | fermi |
|---|---|---|---|
| −0,5 | ~39 000 × 1,7 s ≈ 18 h | ~5 h | ~1,5 h num nó de 24 |
| −0,7 | ~240 000 × 1,7 s ≈ 115 h | ~29 h | ~10 a 16 h num nó (núcleo do fermi ~1,7 vez mais lento) |

O job também calcula o FS nas 24 000 primeiras realizações, o que acrescenta uns 40 h de CPU (~3 h num nó).
