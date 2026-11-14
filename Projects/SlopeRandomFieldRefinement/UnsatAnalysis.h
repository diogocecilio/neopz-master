// Analise de um UnsatModel: equacoes (espaco de cada uma, as de Dirichlet e as livres), percolacao, caminho do fator de
// seguranca (GIM ou SRM, por arc-length) e pos-processamento. Generica no material TMat, como o UnsatModel: o que muda de
// um material para outro (como o fator lambda e a reducao F entram) vem do problema, pela SetLoadFactorFunction e pela
// SetStrengthReductionFunction.
// Os materiais aceitos estao no fim do UnsatAnalysis.cpp (uma linha por material).
#pragma once

#include "UnsatModel.h"

#include "TPZLinearAnalysis.h"

#include <functional>
#include <map>
#include <string>
#include <vector>

template <class TMat>
class UnsatAnalysis : public TPZLinearAnalysis {
public:
    // Matriz esparsa nao simetrica (PARDISO) com LU e as equacoes identificadas (IdentifyEquations)
    explicit UnsatAnalysis(UnsatModel<TMat> &model);

    // ---- o que muda com o material (o problema escolhe)
    // GIM: como o fator lambda entra no material. Ex. u-p:
    //   [mat, d](REAL l) { mat->SetBodyForce(TPZManVector<REAL, 3>({0., -l * d.gamma, 0.})); }
    // scaleFrozen: a analise tambem multiplica as pressoes congeladas por lambda (u-p: sim; o nao saturado escala as
    // pressoes no proprio material, com SetGravityFactor)
    void SetLoadFactorFunction(const std::function<void(REAL)> &f, bool scaleFrozen)
    {
        fSetLoadFactor = f;
        fScaleFrozen = scaleFrozen;
    }
    // Reducao de resistencia F no material (c / F e tan(phi) / F). Obrigatoria no SRM (e o parametro do caminho).
    // No GIM e opcional: o arc-length usa so no preditor tangente do primeiro passo, com F = 1 + 1e-6: no estado
    // convergido os pontos plasticos estao sobre a superficie de escoamento, e a superficie um pouco menor faz eles
    // usarem a tangente plastica em vez da elastica (como o Common, que da 0,885747 na malha grossa de H = 10 m; sem o
    // gancho, 0,885930). Ex. Mohr-Coulomb:
    //   [mat](REAL F) { mat->PlasticModel().SetStrengthReductionFactor(F); }
    // O Cam Clay nao tem reducao de resistencia: deixe sem (e sem SRM)
    void SetStrengthReductionFunction(const std::function<void(REAL)> &f) { fSetReduction = f; }
    // false: sem as linhas de cada iteracao do Newton e de cada passo do arc-length (Monte Carlo)
    void SetVerbose(bool v) { fVerbose = v; }

    // ---- solve
    // Percolacao permanente com o esqueleto parado (u = 0): resolve as pressoes livres e guarda o campo em fX0
    void SolveSeepage();
    // Pressoes congeladas dadas por uma funcao, no lugar da SolveSeepage: fX0 = p(espaco, x) em cada no das malhas de
    // pressao (espaco 1 = p_w, 2 = p_a), deslocamento zero. Ex.: as pressoes de um instante do transiente, calculado em
    // outra malha (a transferencia da malha do transiente para a malha do fator). So os connects de vertice: pressao P1
    void SetFrozenPressures(const std::function<REAL(int espaco, const TPZVec<REAL> &x)> &p);
    // Fator Gamma do GIM: caminho de lambda a partir do estado sem tensoes, com as pressoes congeladas (fX0) e o
    // arc-length cilindrico (Box 4.4 de de Souza Neto et al. 2008); devolve o lambda do colapso
    REAL GIMArcLength();
    // Fator FS do SRM (o fator do equilibrio limite: Bishop, Rahardjo et al. 2007): peso todo (lambda = 1) com as
    // pressoes congeladas e c' / F, tan(phi') / F ate o colapso. A parcela da succao, chi s tan(phi'), cai junto: e o
    // tan(phi^b) / F do Bishop. Precisa da SetStrengthReductionFunction (e da SetLoadFactorFunction, para o peso).
    //  1. peso: lambda de 0 a 1 com F = 1, a partir do estado de partida (ApplyWeight); se o talude nao fica em pe com a
    //     resistencia real (FS < 1), de novo do estado de partida com F = 0,5 (resistencia dobrada) e, se preciso, 0,25
    //  2. caminho de F por arc-length, o mesmo do GIM, com a carga tangente -dR/dF por diferenca finita (ReductionLoad).
    //     Com nenhum ou poucos pontos plasticos essa carga e nula ou pequena (o caminho fica vertical: o preditor pediria
    //     um salto de F maior que 25%): Newton com F 10% maior, e metade disso quando o Newton nao converge
    // Devolve o maior F do caminho. Como o SRM do Common (que deu 2,116 e 1,060 no f50,-5 do Rahardjo em 0 e 24 h)
    REAL SRMArcLength();
    // Estados do ultimo caminho: o aceito no fim (o que fica carregado na malha) e, no SRM, o do fim do peso, antes de
    // reduzir a resistencia (o estado geostatico: o mecanismo e o deslocamento a partir dele)
    const TPZFMatrix<STATE> &PathState() const { return fX; }
    const TPZFMatrix<STATE> &WeightState() const { return fXPeso; }
    // Equilibrio com lambda fixo: Newton nos deslocamentos livres a partir de x (chute inicial e solucao), com as
    // pressoes congeladas e o Dirichlet de u. Converge com |R| < 1e-7 |fbar| nas livres; iter = correcoes feitas.
    // A memoria dos pontos nao muda (so com AcceptSolution); a solucao convergida fica carregada na malha
    bool Newton(TPZFMatrix<STATE> &x, REAL lambda, int &iter);
    // Memoria dos pontos de integracao <- estado da solucao carregada na malha (depois de um Newton convergido)
    void AcceptSolution();

    // ---- postprocess
    // Maior sqrt(J2(eps_p)) dos pontos de integracao de cada elemento, pelo indice geometrico (o indicador do
    // refinamento); 0 nos elementos que nao sao do dominio. Le a memoria: o estado do ultimo AcceptSolution
    std::vector<REAL> PlasticIndicator();
    // Um VTK, <nome>.vtk, com os triangulos do dominio e:
    //   nos nos      os campos da solucao carregada, pelos nomes do material (u-p {"PorePressure"} e {"Displacement"};
    //                nao saturado {"WaterPressure"} e {"Displacement"})
    //   nas celulas  sqrtJ2_eps_p (o PlasticIndicator) e ordem_u (a ordem de u de cada elemento)
    // Depois do GIM, no u-p as pressoes da solucao estao multiplicadas pelo lambda final (o GIM escala tambem o peso da
    // agua); no nao saturado nao (o material escala por dentro). Para ver a percolacao, chame logo depois da SolveSeepage
    void PostProcess(const std::string &nome, const TPZVec<std::string> &scalars, const TPZVec<std::string> &vectors);
    using TPZLinearAnalysis::PostProcess; // a de cima esconderia as PostProcess(int ...) do neopz

    // ---- equacoes
    int Field(int64_t eq) const { return fField[eq]; }                          // espaco: 0 = u, 1 = p_w, 2 = p_a
    const std::map<int64_t, REAL> &DirichletValues() const { return fDirichlet; } // equacao -> valor imposto
    const std::vector<int64_t> &FreeEquations() const { return fFree; }          // u sem apoio: as do caminho do fator
    const TPZFMatrix<STATE> &SeepageSolution() const { return fX0; }             // resultado da SolveSeepage

    // Carrega a solucao na malha multifisica e nas atomicas: o material le u, p_w e p_a das atomicas, e sem isso veria a
    // solucao antiga (o TPZPoroElastoPlasticUPAnalysis e o TPZUnsatAnalysis fazem o mesmo)
    void LoadSolution(const TPZFMatrix<STATE> &sol) override;
    using TPZLinearAnalysis::LoadSolution; // a de cima esconderia as outras LoadSolution do neopz

private:
    // Espaco de cada equacao (pelo indice do connect: os da malha atomica 0 vem primeiro, depois os da 1, ...), as de
    // Dirichlet com o valor (fModel.DirichletBCs()) e as livres
    void IdentifyEquations();
    // Filtro da matriz: so as equacoes da lista entram no sistema
    void SetFreeEquations(const std::vector<int64_t> &free);
    // Assemble() que reaproveita a estrutura da matriz esparsa: o TPZLinearAnalysis so reaproveita a matriz quando ela
    // tem NEquations() linhas, e com o filtro ela tem so as livres, entao a cada montagem ele refazia o grafo e a
    // renumeracao (16% do GIM). Aqui a primeira montagem depois de SetFreeEquations guarda uma copia da matriz ainda nao
    // fatorada (fMolde); as seguintes montam numa copia dela. Copia, e nao a mesma matriz zerada, porque o PARDISO do
    // neopz nao aceita fatorar de novo uma matriz ja fatorada
    void AssembleSystem();
    // Resolve K dx = rhs nas equacoes livres (as do filtro) com a matriz montada; dx tem todas as equacoes (0 nas fixas)
    TPZFMatrix<STATE> SolveFree(const TPZFMatrix<STATE> &rhs);
    // Norma de v so nas equacoes livres (as do filtro)
    REAL FreeNorm(const TPZFMatrix<STATE> &v);
    // Produto escalar so nas equacoes livres (as do filtro)
    REAL FreeDot(const TPZFMatrix<STATE> &a, const TPZFMatrix<STATE> &b);
    // Copia de v com zero nas equacoes fixas
    TPZFMatrix<STATE> FreeOnly(const TPZFMatrix<STATE> &v);
    // Valores das equacoes fixas do caminho do fator em x: pressoes congeladas (lambda fX0 com fScaleFrozen, senao
    // fX0; sem SolveSeepage, zero: talude seco) e o Dirichlet de u
    void SetFixed(TPZFMatrix<STATE> &x, REAL lambda) const;
    // -R(x, lambda) em todas as equacoes (o vetor do neopz e -R), com as fixas de SetFixed; so o vetor, sem a matriz
    TPZFMatrix<STATE> ResidualAt(TPZFMatrix<STATE> x, REAL lambda);
    // Carga externa por unidade de lambda, fbar = R(x, 0) - R(x, 1): com as pressoes congeladas o residuo e linear em
    // lambda, R(x, lambda) = R(x, 0) - lambda fbar (fbar: o peso e as pressoes). Deixa lambda = 0 no material
    TPZFMatrix<STATE> LoadVector(const TPZFMatrix<STATE> &x);
    // Carga tangente do SRM, fbar = -dR/dF em x, por diferenca finita: rhs = -R(x) com a reducao F e o residuo com
    // F (1 + 1e-6). So os pontos plasticos contribuem (no elastico a tensao nao depende de F). Deixa F no material
    TPZFMatrix<STATE> ReductionLoad(const TPZFMatrix<STATE> &x, const TPZFMatrix<STATE> &rhs, REAL F);
    // Parametro v do caminho no material: GIM lambda = v (e, com o gancho, F = 1 + 1e-6 no preditor tangente); SRM
    // lambda = 1 e F = v (v (1 + 1e-6) no preditor)
    void SetParameter(REAL v, bool preditor);
    // Peso de lambda = 0 a 1 com a reducao que estiver no material, a partir de x (o estado de partida): passos de
    // Newton, cada um aceito na memoria; o passo comeca em 0,5, cresce 1,5 vez quando converge em ate 8 iteracoes e cai
    // pela metade quando falha (a continuacao do Common). false: o passo ficou menor que 0,002 antes de lambda = 1.
    // x: o estado em lambda = 1
    bool ApplyWeight(TPZFMatrix<STATE> &x);
    // Um incremento do arc-length a partir do ultimo estado aceito (fX, fV); l: comprimento do arco; x e v: o novo
    // estado (carregado na malha se convergiu); iter: montagens feitas; vertical: falhou porque a carga tangente e nula
    // (SRM sem pontos plasticos)
    bool ArcLengthStep(REAL l, TPZFMatrix<STATE> &x, REAL &v, int &iter, bool &vertical);
    // Caminho do arc-length a partir do estado aceito (fX, fV), o mesmo no GIM e no SRM; devolve o maior v e o motivo
    // do fim
    REAL ArcLengthPath(std::string &fim);

    UnsatModel<TMat> &fModel;
    std::vector<int> fField;            // espaco de cada equacao
    std::map<int64_t, REAL> fDirichlet; // equacao de Dirichlet -> valor
    std::vector<int64_t> fFree;         // equacoes livres de u
    TPZFMatrix<STATE> fX0;              // pressoes da percolacao com lambda = 1 (congeladas no GIM)
    std::function<void(REAL)> fSetLoadFactor;
    bool fScaleFrozen = true;
    std::function<void(REAL)> fSetReduction; // opcional (vazio: sem reducao no preditor)
    bool fVerbose = true;

    // estado do caminho do arc-length (GIMArcLength e SRMArcLength)
    bool fSRM = false;       // parametro do caminho: lambda (GIM) ou F (SRM)
    TPZFMatrix<STATE> fX;    // ultimo estado aceito (todas as equacoes)
    TPZFMatrix<STATE> fXPeso; // SRM: estado do fim do peso (lambda = 1, antes do caminho de F)
    REAL fV = 0.;            // parametro (lambda ou F) do ultimo estado aceito
    TPZFMatrix<STATE> fDu;   // incremento do ultimo passo aceito, so nas livres (vazio: ainda nenhum)
    REAL fDv = 0.;           // incremento do parametro no ultimo passo aceito
    TPZFMatrix<STATE> fFbar; // GIM: carga por unidade de lambda (LoadVector)
    REAL fFRef = 1.;         // |carga do peso| nas livres: referencia do residuo
    TPZAutoPointer<TPZMatrix<STATE>> fMolde; // matriz das livres sem fatorar (AssembleSystem)
};
