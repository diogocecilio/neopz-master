# SlopeMohrCoulomb3D

Versão 3D de `../SlopeMohrCoulomb`: o mesmo talude (70 × 40 m, altura 10 m a 45°), extrudado
**40 m** na direção z, Mohr–Coulomb associativo perfeitamente plástico, FS por **aumento de
gravidade** e por **redução de resistência (SRM)**. Mesmo solo (E = 20000 kPa, ν = 0,49,
c = 10 kPa, φ = ψ = 30°, γ = 20 kN/m³), mesmo modelo constitutivo padrão
(`TPZPlasticStepVoigt<TPZYCMohrCoulombPV2>`, opção `pv` para o modelo antigo) e o mesmo
algoritmo (Newton com tangente consistente e busca linear, continuação com divisão do passo).

**Ainda não foi executado** (apenas checagem de sintaxe da compilação).

## Malha

* Base 2D: a malha do projeto 2D com **um nível de refinamento uniforme** (`TriGMesh(1)`,
  elementos de 5 m), gerada por subdivisão pelos pontos médios (`ref=<n>` muda o nível).
* Extrusão em z com camadas do mesmo tamanho (5 m → 8 camadas em 40 m; `nz=`, `lz=`).
* **Tetraedros (padrão)**: cada prisma triangular é dividido em 3 tetraedros, com as diagonais
  das faces escolhidas pela numeração global dos nós (malha conforme).
* **Hexaedros (`hexa`, preparado)**: versão em quadriláteros da malha 2D (quadrados de 10 m
  abaixo de y = 30 e uma faixa mapeada entre y = 30 e a crista, conforme em y = 30),
  refinada do mesmo modo e extrudada.
* Contorno: base fixa (−1), laterais x = 0 (−5) e x = 70 (−2) em rolete (uₓ = 0), faces
  z = 0 (−7) e z = 40 (−8) em rolete (u_z = 0). Com isso a solução exata é a de deformação
  plana do projeto 2D (FS esperados próximos dos da malha 2D equivalente: 3,04 e 1,40 no ciclo 0).
* Ordem polinomial 2 (`p=`), como no 2D.

## Diferenças em relação ao driver 2D

* `SlopeAnalysis3D.h`: cópia de `SlopeAnalysis.h` com o material 3D
  (`TPZMatElastoPlastic`) e pós-processamento VTK em 3D.
* `TPZMatElastoPlasticGravity` (em `main.cpp`): `TPZMatElastoPlastic::Contribute` **ignora o
  peso próprio de `SetBodyForce`** (só usa uma forcing function multiplicada pela densidade,
  que vale 0 por padrão). A subclasse soma `m_force` por unidade de volume, como faz o
  `TPZMatElastoPlastic2D`. A correção na biblioteca não foi feita.
* `main.cpp` instancia `TPZMatElastoPlastic<TPZPlasticStepVoigt<TPZYCMohrCoulombPV2>>`, que a
  biblioteca só instancia para o 2D.
* O `ContributeBC` 3D usa penalidade fixa 1e16 (o 2D usa `SetBigNumber(1e12)`).
* Refinamento adaptativo (`nref`) desligado por padrão: o padrão uniforme de tetraedro do NeoPZ
  gera 4 tetraedros + 2 pirâmides, então refinar deixa de ser "só tetraedros". Hexaedros dividem
  em 8 hexaedros.

## Compilar e executar

```
cmake -DBUILD_PLASTICITY_MATERIALS=ON -DBUILD_PROJECTS=ON <neopz>
ninja SlopeMohrCoulomb3D
./SlopeMohrCoulomb3D mesh       # só gera a malha (slope3d_*_gmesh.vtk) e informa o nº de equações
./SlopeMohrCoulomb3D            # tetraedros, FS por gravidade e por SRM
./SlopeMohrCoulomb3D hexa       # hexaedros
./SlopeMohrCoulomb3D pv         # modelo antigo (verificação)
./SlopeMohrCoulomb3D nz=1 p=1   # opções para testes rápidos
```

Estimativa de tamanho (P2, 8 camadas): da ordem de 2–3·10⁴ equações; o solver é a skyline
LDLᵀ com renumeração de banda, como no 2D, e pode ficar lento em 3D.
