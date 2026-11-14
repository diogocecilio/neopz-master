#include "UnsatAnalysis.h"

#include "TPZSpStructMatrix.h"
#include "TPZVTKGeoMesh.h"
#include "pzstepsolver.h"
#include "pzcompel.h"
#include "pzgeoel.h"
#include "TPZMaterial.h"
#include "pzmultiphysicselement.h"
#include "pzinterpolationspace.h"
#include "pzerror.h"

#include <cstdlib>

#include <cmath>
#include <cstdio>
#include <fstream>

// materiais das instancias do fim do arquivo
#include "Plasticity/TPZElasticResponse.h"
#include "Plasticity/TPZMatPoroElastoPlasticUP.h"
#include "Plasticity/TPZMatUnsaturatedUPwPa.h"
#include "Plasticity/TPZPlasticStepModifiedCamClay.h"
#include "Plasticity/TPZPlasticStepVoigt.h"
#include "Plasticity/TPZYCMohrCoulombPV2.h"
#include "HPMaterial.h"

// false: sem renumerar as equacoes, como o TPZPoroElastoPlasticUPAnalysis
template <class TMat>
UnsatAnalysis<TMat>::UnsatAnalysis(UnsatModel<TMat> &model) : TPZLinearAnalysis(model.CompMesh(), false), fModel(model)
{
    // matriz esparsa nao simetrica (so guarda os nao nulos): o acoplamento u-p (-Q e Q^T), a tangente plastica e a
    // permeabilidade que depende da pressao tiram a simetria. A analise guarda uma copia
    TPZSpStructMatrix<STATE> strmat(model.CompMesh());
    // threads da montagem: 8, ou o valor da variavel de ambiente PZ_THREADS (0 = serial; no cluster, com um processo
    // por nucleo, use PZ_THREADS=0)
    const char *nth = std::getenv("PZ_THREADS");
    strmat.SetNumThreads(nth ? std::atoi(nth) : 8);
    SetStructuralMatrix(strmat);

    // solver direto:  fatoracao LU; com a matriz esparsa e o MKL, quem fatora e o PARDISO
    TPZStepSolver<STATE> solver;
    solver.SetDirect(ELU);
    SetSolver(solver);

    IdentifyEquations(); // espaco de cada equacao, as de Dirichlet e as livres
}

template <class TMat>
void UnsatAnalysis<TMat>::IdentifyEquations()
{
    TPZMultiphysicsCompMesh *cmesh = fModel.CompMesh();
    TPZVec<TPZCompMesh *> &meshvec = cmesh->MeshVector();
    const int nspaces = meshvec.size();
    TPZBlock &block = cmesh->Block();
    const int64_t neq = cmesh->NEquations();

    // Os connects da malha multifisica sao os das atomicas em sequencia: primeiro os da malha 0 (u), depois os da 1
    // (p_w), ... O connect ic da malha atomica s e o connect offset[s] + ic da multifisica
    TPZManVector<int64_t, 4> offset(nspaces + 1, 0);
    for (int s = 0; s < nspaces; s++) offset[s + 1] = offset[s] + meshvec[s]->NConnects();

    // 1. espaco de cada equacao: as equacoes de um connect estao em Position(seq) ... Position(seq) + Size(seq) - 1
    fField.assign(neq, -1);
    for (int s = 0; s < nspaces; s++) {
        for (int64_t ic = offset[s]; ic < offset[s + 1]; ic++) {
            const int64_t seq = cmesh->ConnectVec()[ic].SequenceNumber();
            if (seq < 0) continue; // connect sem uso
            const int64_t pos = block.Position(seq);
            for (int k = 0; k < block.Size(seq); k++) fField[pos + k] = s;
        }
    }

    // 2. Dirichlet: os connects dos elementos de contorno de cada BC, no espaco da BC. Num connect com nstate
    // componentes, a equacao pos + k e da componente k % nstate (um connect de ordem alta tem varias funcoes)
    fDirichlet.clear();
    const int matid = fModel.Material()->Id();
    for (const auto &bc : fModel.DirichletBCs()) {
        TPZCompMesh *atomic = meshvec[bc.field];
        const int nstate = atomic->FindMaterial(matid)->NStateVariables();
        for (TPZCompEl *cel : atomic->ElementVec()) {
            if (!cel || !cel->Reference() || cel->Reference()->MaterialId() != bc.id) continue;
            TPZGeoEl *gel = cel->Reference();
            for (int i = 0; i < cel->NConnects(); i++) {
                const int64_t seq = cmesh->ConnectVec()[offset[bc.field] + cel->ConnectIndex(i)].SequenceNumber();
                const int64_t pos = block.Position(seq);
                // valor: f no no para os connects de vertice (os NCornerNodes() primeiros); nos outros (meio de aresta)
                // 0, o coeficiente hierarquico de um valor linear na aresta (f de grau maior nao e reproduzido)
                TPZManVector<STATE, 3> v2(nstate, 0.);
                if (bc.f && i < gel->NCornerNodes()) {
                    TPZManVector<REAL, 3> x(3, 0.);
                    gel->NodePtr(i)->GetCoordinates(x);
                    TPZFMatrix<STATE> v1(nstate, nstate, 0.);
                    bc.f(x, v2, v1);
                }
                for (int k = 0; k < block.Size(seq); k++)
                    for (int c : bc.comps)
                        if (k % nstate == c) fDirichlet[pos + k] = (k < nstate) ? v2[c] : 0.;
            }
        }
    }

    // 3. livres: deslocamentos sem Dirichlet (as do caminho do fator; as pressoes ficam congeladas)
    fFree.clear();
    for (int64_t eq = 0; eq < neq; eq++)
        if (fField[eq] == 0 && !fDirichlet.count(eq)) fFree.push_back(eq);
}

template <class TMat>
void UnsatAnalysis<TMat>::SetFreeEquations(const std::vector<int64_t> &free)
{
    // O filtro mora na matriz estrutural: no Assemble so as equacoes ativas entram, e a matriz e o vetor ficam com
    // free.size() linhas (na ordem crescente das equacoes). As outras ficam com o valor que estiver na solucao
    TPZVec<int64_t> ativas(free.size());
    for (size_t i = 0; i < free.size(); i++) ativas[i] = free[i];
    TPZEquationFilter &filtro = StructMatrix()->EquationFilter();
    filtro.Reset(); // SetActiveEquations para se o filtro ja tiver equacoes (de uma chamada anterior)
    filtro.SetActiveEquations(ativas);
    // a matriz esparsa guarda a estrutura (onde ha nao nulos) das equacoes do filtro antigo, e o Assemble so cria outra
    // se o numero de linhas mudar: apagando, o proximo Assemble cria a matriz com a estrutura do filtro novo
    if (fSolver) MatrixSolver<STATE>().ResetMatrix();
    fMolde = nullptr; // a estrutura das livres mudou
}

template <class TMat>
void UnsatAnalysis<TMat>::LoadSolution(const TPZFMatrix<STATE> &sol)
{
    TPZLinearAnalysis::LoadSolution(sol);              // solucao na malha multifisica
    fModel.CompMesh()->LoadSolutionFromMultiPhysics(); // e nas atomicas, de onde o material le os campos
}

template <class TMat>
void UnsatAnalysis<TMat>::AssembleSystem()
{
    TPZMatrixSolver<STATE> &sol = MatrixSolver<STATE>();
    if (!fMolde) { // primeira montagem com este filtro: cria a matriz e guarda uma copia antes de fatorar
        Assemble();
        fMolde = dynamic_cast<TPZMatrix<STATE> *>(sol.Matrix()->Clone());
        return;
    }
    sol.SetMatrix(dynamic_cast<TPZMatrix<STATE> *>(fMolde->Clone()));
    sol.Matrix()->Zero();
    fRhs.Redim(fCompMesh->NEquations(), 1);
    fStructMatrix->Assemble(*sol.Matrix(), fRhs, fGuiInterface);
}

template <class TMat>
TPZFMatrix<STATE> UnsatAnalysis<TMat>::SolveFree(const TPZFMatrix<STATE> &rhs)
{
    // a matriz so tem as livres e o vetor tem todas: junta as livres do vetor (Gather), resolve e espalha a resposta de
    // volta nas posicoes delas (Scatter). Nao usa o Solve() do TPZLinearAnalysis, que carregaria dx na malha
    TPZEquationFilter &filtro = StructMatrix()->EquationFilter();
    TPZFMatrix<STATE> br(filtro.NActiveEquations(), 1, 0.), yr, y(rhs.Rows(), 1, 0.);
    filtro.Gather(rhs, br);
    MatrixSolver<STATE>().Solve(br, yr); // na primeira vez fatora (LU do PARDISO)
    filtro.Scatter(yr, y);
    return y;
}

template <class TMat>
REAL UnsatAnalysis<TMat>::FreeNorm(const TPZFMatrix<STATE> &v)
{
    TPZEquationFilter &filtro = StructMatrix()->EquationFilter();
    TPZFMatrix<STATE> vr(filtro.NActiveEquations(), 1, 0.);
    filtro.Gather(v, vr);
    return Norm(vr);
}

template <class TMat>
REAL UnsatAnalysis<TMat>::FreeDot(const TPZFMatrix<STATE> &a, const TPZFMatrix<STATE> &b)
{
    TPZEquationFilter &filtro = StructMatrix()->EquationFilter();
    TPZFMatrix<STATE> ar(filtro.NActiveEquations(), 1, 0.), br(filtro.NActiveEquations(), 1, 0.);
    filtro.Gather(a, ar);
    filtro.Gather(b, br);
    return Dot(ar, br);
}

template <class TMat>
TPZFMatrix<STATE> UnsatAnalysis<TMat>::FreeOnly(const TPZFMatrix<STATE> &v)
{
    TPZEquationFilter &filtro = StructMatrix()->EquationFilter();
    TPZFMatrix<STATE> vr(filtro.NActiveEquations(), 1, 0.), w(v.Rows(), 1, 0.);
    filtro.Gather(v, vr);
    filtro.Scatter(vr, w);
    return w;
}

// Percolacao permanente depois do rebaixamento, com o esqueleto parado:
//   livres: p_w sem Dirichlet          fixas: u = 0 (esqueleto rigido), p_w de Dirichlet, p_a = 0 (ar passivo)
// Newton nas livres: no u-p o problema e linear (converge numa iteracao); no nao saturado a permeabilidade cai onde
// p_w < 0 e o problema e nao linear
template <class TMat>
void UnsatAnalysis<TMat>::SolveSeepage()
{
    TMat *mat = fModel.Material();
    const int64_t neq = fModel.CompMesh()->NEquations();

    // estado de partida: tudo zero e os valores de Dirichlet das pressoes (u fica em zero mesmo com Dirichlet de u: o
    // esqueleto nao se move na percolacao)
    TPZFMatrix<STATE> x(neq, 1, 0.);
    for (const auto &[eq, v] : fDirichlet)
        if (fField[eq] > 0) x(eq, 0) = v;

    // livres: so as pressoes da agua (espaco 1) sem Dirichlet; todas as outras ficam com o valor de x
    std::vector<int64_t> livres;
    for (int64_t eq = 0; eq < neq; eq++)
        if (fField[eq] == 1 && !fDirichlet.count(eq)) livres.push_back(eq);
    SetFreeEquations(livres);

    // Regime permanente: o material faz um passo implicito no tempo, R = armazenamento + dt x fluxo. Com dt enorme o
    // armazenamento (que nao cresce com dt) fica desprezivel e o passo e o permanente, com ou sem armazenamento no
    // material (a biblioteca usava dt = 1 e zerava o armazenamento no material: Kw = 1e30, S_w = 1)
    mat->SetTimeStep(1.e12);
    REAL r0 = 0.;
    bool convergiu = false;
    for (int it = 0; it < 30; it++) {
        LoadSolution(x);
        AssembleSystem(); // Rhs() = -R em todas as equacoes; a matriz tangente so nas livres
        const REAL r = FreeNorm(Rhs());
        if (it == 0) r0 = r;
        if (fVerbose) std::cout << "  percolacao: iteracao " << it << "  |R|/|R0| = " << (r0 > 0. ? r / r0 : 0.) << std::endl;
        if (r <= 1.e-10 * r0) {
            convergiu = true;
            break;
        }
        x += SolveFree(Rhs()); // Newton: K dx = -R
    }
    mat->SetTimeStep(0.); // de volta ao padrao (sem fluxo): no GIM as pressoes ficam congeladas
    if (!convergiu) DebugStop(); // a percolacao nao convergiu em 30 iteracoes
    fX0 = x;
    LoadSolution(fX0); // o campo de pressoes fica na malha (pos-processamento)
}

// Os connects das malhas de pressao vem depois dos de u na malha multifisica (offset = NConnects das anteriores). Em
// cada elemento (do dominio e de contorno), os NCornerNodes() primeiros connects sao os dos vertices; num no comum a
// varios elementos o valor e o mesmo (a funcao so depende da posicao)
template <class TMat>
void UnsatAnalysis<TMat>::SetFrozenPressures(const std::function<REAL(int espaco, const TPZVec<REAL> &x)> &p)
{
    TPZMultiphysicsCompMesh *cmesh = fModel.CompMesh();
    TPZVec<TPZCompMesh *> &meshvec = cmesh->MeshVector();
    TPZBlock &block = cmesh->Block();
    fX0.Redim(cmesh->NEquations(), 1);
    int64_t offset = meshvec[0]->NConnects();
    for (int s = 1; s < meshvec.size(); s++) {
        for (TPZCompEl *cel : meshvec[s]->ElementVec()) {
            if (!cel || !cel->Reference()) continue;
            TPZGeoEl *gel = cel->Reference();
            for (int i = 0; i < gel->NCornerNodes() && i < cel->NConnects(); i++) {
                const int64_t seq = cmesh->ConnectVec()[offset + cel->ConnectIndex(i)].SequenceNumber();
                if (seq < 0) continue;
                TPZManVector<REAL, 3> x(3, 0.);
                gel->NodePtr(i)->GetCoordinates(x);
                fX0(block.Position(seq), 0) = p(s, x);
            }
        }
        offset += meshvec[s]->NConnects();
    }
}

template <class TMat>
void UnsatAnalysis<TMat>::SetFixed(TPZFMatrix<STATE> &x, REAL lambda) const
{
    const REAL s = fScaleFrozen ? lambda : 1.; // u-p: pressoes x lambda; nao saturado: o material escala
    const bool seepage = fX0.Rows() == x.Rows();
    for (int64_t eq = 0; eq < x.Rows(); eq++)
        if (fField[eq] > 0) x(eq, 0) = seepage ? s * fX0.GetVal(eq, 0) : 0.;
    for (const auto &[eq, v] : fDirichlet)
        if (fField[eq] == 0) x(eq, 0) = v;
}

template <class TMat>
TPZFMatrix<STATE> UnsatAnalysis<TMat>::ResidualAt(TPZFMatrix<STATE> x, REAL lambda)
{
    fSetLoadFactor(lambda);
    SetFixed(x, lambda);
    LoadSolution(x);
    AssembleResidual();
    return Rhs();
}

template <class TMat>
TPZFMatrix<STATE> UnsatAnalysis<TMat>::LoadVector(const TPZFMatrix<STATE> &x)
{
    // -R(x, 1) - (-R(x, 0)) = R(x, 0) - R(x, 1): a parte interna (tensoes de x) e a mesma nas duas e se cancela
    TPZFMatrix<STATE> fbar = ResidualAt(x, 1.);
    fbar -= ResidualAt(x, 0.);
    return fbar;
}

// Newton com lambda fixo: K(x) dx = -R(x) nas livres, x <- x + dx, ate |R| < 1e-7 |fbar| (a tolerancia do Common e
// do TPZUnsatAnalysis). Nas fixas, x fica com os valores de SetFixed. As tensoes de cada iteracao saem da memoria
// (o ultimo estado aceito) e da deformacao de x; a memoria so muda no AcceptSolution
template <class TMat>
bool UnsatAnalysis<TMat>::Newton(TPZFMatrix<STATE> &x, REAL lambda, int &iter)
{
    if (!fSetLoadFactor) { // o problema tem que dizer como lambda entra no material
        std::cerr << "UnsatAnalysis::Newton: falta SetLoadFactorFunction (como lambda entra no material)" << std::endl;
        DebugStop();
    }
    const REAL tol = 1.e-7;
    const int maxit = 20;
    SetFreeEquations(fFree);
    const REAL fref = FreeNorm(LoadVector(x)); // referencia: a carga externa com lambda = 1, nas livres
    fSetLoadFactor(lambda);
    SetFixed(x, lambda);
    for (iter = 0;; iter++) {
        LoadSolution(x);
        AssembleSystem(); // Rhs() = -R(x) e a tangente K(x) nas livres
        const REAL e = FreeNorm(Rhs()) / fref;
        if (fVerbose)
            std::cout << "    Newton (lambda = " << lambda << "): iteracao " << iter << "  |R|/|fbar| = " << e << std::endl;
        if (!std::isfinite(e)) return false; // divergiu
        if (e < tol) return true;            // x convergido e carregado na malha
        if (iter == maxit) return false;
        try {
            x += SolveFree(Rhs()); // dx so nas livres
        } catch (...) {           // tangente singular (colapso)
            return false;
        }
    }
}

// SRM: fbar = -dR/dF = (R(x, F) - R(x, F + h)) / h; o vetor do neopz e -R, entao fbar = (rhs(F + h) - rhs(F)) / h
template <class TMat>
TPZFMatrix<STATE> UnsatAnalysis<TMat>::ReductionLoad(const TPZFMatrix<STATE> &x, const TPZFMatrix<STATE> &rhs, REAL F)
{
    const REAL h = 1.e-6 * F;
    fSetReduction(F + h);
    TPZFMatrix<STATE> f = ResidualAt(x, 1.); // -R(x) com F + h (so o vetor: a tangente montada fica)
    fSetReduction(F);
    f -= rhs;
    f *= 1. / h;
    return f;
}

template <class TMat>
void UnsatAnalysis<TMat>::SetParameter(REAL v, bool preditor)
{
    if (fSRM) {
        fSetLoadFactor(1.);
        fSetReduction(preditor ? v * (1. + 1.e-6) : v);
        return;
    }
    fSetLoadFactor(v);
    if (fSetReduction) fSetReduction(preditor ? 1. + 1.e-6 : 1.);
}

// Um incremento do arc-length cilindrico (de Souza Neto et al. 2008, Box 4.4) a partir do ultimo estado aceito
// (fX, fV): acha x e o parametro v (lambda no GIM, F no SRM) com R(x, v) = 0 e |Du| = l, onde Du = x - fX nas livres.
//  preditor: secante (o incremento aceito anterior fDu, esticado ate o comprimento l) ou, no primeiro passo do
//            caminho, tangente (dv = l / |dwb|, Eq. 4.122, com o sinal do criterio secante, Eq. 4.123)
//  corretor: dws = K^-1 (-R) e dwb = K^-1 fbar (os dws e dwb do seu IterativeProcessArcLength); du = dws + dv dwb,
//            com dv da restricao |Du + du| = l (Eqs. 4.116-4.117) e, das duas raizes, a de menor angulo (Eq. 4.118).
//            fbar = -dR/dv: no GIM a carga fFbar, a mesma no caminho todo; no SRM a ReductionLoad, de novo a cada
//            montagem (depende do estado: so os pontos plasticos contribuem)
template <class TMat>
bool UnsatAnalysis<TMat>::ArcLengthStep(REAL l, TPZFMatrix<STATE> &x, REAL &v, int &iter, bool &vertical)
{
    const REAL tol = 1.e-7; // |R| / |fbar| nas livres, como o Newton
    const int kmax = 12;    // passos contando o preditor
    REAL ePrev = -1.;       // |R| / |fbar| da montagem anterior do passo (regra de divergencia)
    x = fX;
    v = fV;
    vertical = false;
    TPZFMatrix<STATE> Du(x.Rows(), 1, 0.); // incremento do passo, so nas livres: x = fX + Du
    iter = 0;

    // preditor secante: o proprio incremento anterior da a direcao (e o sinal); so se l nao for desproporcional a ele
    const REAL normDun = fDu.Rows() ? FreeNorm(fDu) : 0.;
    const bool secante = normDun > 0. && l < 100. * normDun;
    int k = 1;
    if (secante) {
        const REAL s = l / normDun;
        Du = fDu;
        Du *= s;
        x += Du;
        v += s * fDv;
        k = 2; // o preditor ja ocupou o primeiro passo
    }
    const int kfirst = k;
    for (; k <= kmax; k++) {
        iter = k - kfirst + 1;       // montagens ate aqui
        if (!(v > 0.)) return false; // o caminho voltou para v <= 0
        SetParameter(v, k == 1);     // superficie um pouco menor so no preditor tangente
        SetFixed(x, fSRM ? 1. : v);
        LoadSolution(x);
        AssembleSystem();                    // Rhs() = -R(x, v) e a tangente K(x) nas livres
        const TPZFMatrix<STATE> rhs = Rhs(); // copia: a carga do SRM monta outro residuo no lugar
        const REAL e = FreeNorm(rhs) / fFRef;
        // divergiu: o residuo passou da propria carga fbar (nos passos que convergem ele fica abaixo de ~0,1); antes o
        // limite era 1e6 e os passos condenados gastavam todas as montagens ate o kmax antes de o arco ser reduzido
        if (!std::isfinite(e) || e > 1.) return false;
        if (k > 1 && e < tol) return true; // convergiu (no k = 1, x ainda e o estado aceito)
        // da 5a montagem em diante, residuo subindo: o passo nao vai convergir, desiste ja (o arco e reduzido)
        if (iter >= 5 && ePrev > 0. && e > ePrev) return false;
        ePrev = e;
        TPZFMatrix<STATE> dws, dwb;
        try {
            dws = SolveFree(rhs); // K dws = -R: a correcao com v fixo
            // K dwb = fbar: quanto x anda por unidade de v (mesma fatoracao)
            dwb = SolveFree(fSRM ? ReductionLoad(x, rhs, k == 1 ? v * (1. + 1.e-6) : v) : fFbar);
        } catch (...) {
            return false; // tangente singular
        }
        const REAL a = FreeDot(dwb, dwb);
        if (!(a > 0.)) {
            vertical = a == 0.; // carga tangente nula: nenhum ponto plastico (SRM)
            return false;
        }
        REAL dv;
        if (k == 1) {
            // preditor tangente: |dv dwb| = l, com o sinal de fDu . dwb (o seu computelamda0 usava u, zerado a cada
            // passo, e o sinal saia sempre +)
            const REAL sinal = (fDu.Rows() && FreeDot(fDu, dwb) < 0.) ? -1. : 1.;
            dv = sinal * l / std::sqrt(a);
            // SRM com poucos pontos plasticos (ex.: o comeco com a resistencia dobrada): a carga tangente e pequena e o
            // preditor pede um salto enorme de F, que o corretor nao alcanca. O caminho ainda e quase vertical
            if (fSRM && std::abs(dv) > 0.25 * v) {
                vertical = true;
                return false;
            }
        } else {
            // restricao |Du + dws + dv dwb| = l: a dv^2 + b dv + c = 0 (o aa, bb e cc do seu computelamda)
            TPZFMatrix<STATE> w(Du);
            w += dws;
            const REAL b = 2. * FreeDot(w, dwb), c = FreeDot(w, w) - l * l, disc = b * b - 4. * a * c;
            if (disc < 0.) return false; // sem raiz real: o arco nao corta o caminho perto daqui
            // raizes sem cancelamento numerico: q = -(b + sinal(b) sqrt(disc)) / 2, r1 = q / a, r2 = c / q
            const REAL q = -0.5 * (b + (b >= 0. ? 1. : -1.) * std::sqrt(disc));
            const REAL r1 = q / a, r2 = q != 0. ? c / q : r1;
            // menor angulo entre o incremento novo e o atual: maximiza (Du + dws + dv dwb) . Du, isto e,
            // dv (dwb . Du). Nao precisa das duas montagens extras da escolha pelo residuo do computelamdacris
            const REAL dir = FreeDot(dwb, Du);
            dv = r1 * dir >= r2 * dir ? r1 : r2;
        }
        v += dv;
        TPZFMatrix<STATE> du(dwb);
        du *= dv;
        du += dws; // du = dws + dv dwb, so nas livres
        x += du;
        Du += du;
    }
    return false; // nao convergiu em kmax passos
}

// O caminho, com os parametros do Common (os mesmos no GIM e no SRM):
//  - o arco l comeca em 0,25 |u| do estado inicial e se ajusta para ~10 montagens por passo (entre metade e o dobro,
//    no maximo 2 |u|); se o passo falha, l pela metade, do mesmo estado aceito;
//  - passo vertical (SRM com nenhum ou poucos pontos plasticos: F quase nao muda o residuo e o arco nao anda): Newton
//    com F 10% maior a partir do estado aceito; se o Newton nao converge, o colapso esta antes: o aumento cai pela
//    metade (ate 0,1%: fim "vertical", FS com essa precisao). O passo seguinte volta ao preditor tangente;
//  - fim: patamar (v cresce menos de 0,001 num passo com o l pedido), v >= 10, l < 1e-4 |u| ou 100 passos.
template <class TMat>
REAL UnsatAnalysis<TMat>::ArcLengthPath(std::string &fim)
{
    const REAL tolPatamar = 0.001, limite = 10.;
    const int alvo = 10, maxPassos = 100;
    const REAL uRef = FreeNorm(fX); // tamanho dos deslocamentos do estado inicial: a escala do arco
    const REAL lMax = 2. * uRef, lMin = 1.e-4 * uRef;
    REAL l = 0.25 * uRef, lPedido = l, maior = fV;
    REAL aumento = 0.1; // passo vertical: F cresce 10% (metade a cada Newton que falha)
    fim = "100 passos";
    for (int passo = 1; passo <= maxPassos; passo++) {
        TPZFMatrix<STATE> xn;
        REAL vn = 0.;
        int it = 0;
        bool vertical = false;
        if (!ArcLengthStep(l, xn, vn, it, vertical)) {
            if (vertical && fSRM) {
                xn = fX;
                vn = (1. + aumento) * fV;
                SetParameter(vn, false);
                if (!Newton(xn, 1., it)) {
                    if (fVerbose)
                        std::printf("    %5d %10.5f %10s %4d %11s  falhou: aumento pela metade\n", passo, vn, "", it,
                                    "(vertical)");
                    aumento *= 0.5;
                    if (aumento < 1.e-3) {
                        fim = "vertical";
                        break;
                    }
                    continue;
                }
                AcceptSolution();
                fX = xn;
                fV = vn;
                fDu.Resize(0, 0);
                fDv = 0.;
                maior = std::max(maior, fV);
                if (fVerbose) std::printf("    %5d %10.5f %10s %4d %11s\n", passo, fV, "", it, "(vertical)");
                continue;
            }
            if (fVerbose) std::printf("    %5d %10s %10s %4d %11.3e  falhou: arco pela metade\n", passo, "", "", it, l);
            l *= 0.5;
            if (l < lMin) {
                fim = "arco minimo";
                break;
            }
            continue;
        }
        AcceptSolution(); // memoria <- estado convergido (carregado na malha pelo ArcLengthStep)
        TPZFMatrix<STATE> du(xn);
        du -= fX;
        fDu = FreeOnly(du);
        fDv = vn - fV;
        fX = xn;
        fV = vn;
        maior = std::max(maior, fV);
        if (fVerbose) std::printf("    %5d %10.5f %+10.5f %4d %11.3e\n", passo, fV, fDv, it, l);
        if (fDv * lPedido / l < tolPatamar) { // v parou de crescer: o colapso
            fim = "patamar";
            break;
        }
        if (fV >= limite) {
            fim = fSRM ? "F limite" : "lambda limite";
            break;
        }
        l = std::min(lMax, l * std::min(2., std::max(0.5, std::sqrt(REAL(alvo) / it))));
        lPedido = l;
    }
    return maior;
}

// Fator Gamma do GIM por arc-length, com os parametros do Common (que deu 0,885747 na malha grossa de H = 10 m):
//  1. primeiro ponto: Newton com lambda = 0,5 a partir do estado sem tensoes (metade se nao convergir)
//  2. o caminho de lambda (ArcLengthPath). Gamma = o maior lambda do caminho
template <class TMat>
REAL UnsatAnalysis<TMat>::GIMArcLength()
{
    if (!fSetLoadFactor) { // o problema tem que dizer como lambda entra no material
        std::cerr << "UnsatAnalysis::GIMArcLength: falta SetLoadFactorFunction (como lambda entra no material)" << std::endl;
        DebugStop();
    }
    const int64_t neq = fModel.CompMesh()->NEquations();
    fSRM = false;

    // 1. primeiro ponto do caminho
    TPZFMatrix<STATE> x(neq, 1, 0.);
    REAL lambda = 0.5;
    int it = 0;
    while (!Newton(x, lambda, it)) {
        lambda *= 0.5;
        if (lambda < 0.01) {
            std::cout << "  GIM: sem equilibrio nem com lambda pequeno" << std::endl;
            return 0.;
        }
        x.Zero(); // de novo do estado sem tensoes
    }
    AcceptSolution();
    fX = x;
    fV = lambda;
    fDu.Resize(0, 0);
    fDv = 0.;
    fFbar = LoadVector(fX); // a mesma em todo o caminho (pressoes congeladas)
    fFRef = FreeNorm(fFbar);

    // 2. caminho
    if (fVerbose) {
        std::printf("  GIM arc-length\n    %5s %10s %10s %4s %11s\n", "passo", "lambda", "dlambda", "it", "arco l");
        std::printf("    %5d %10.5f %10s %4d %11s\n", 0, fV, "", it, "(Newton)");
    }
    std::string fim;
    const REAL gamma = ArcLengthPath(fim);
    if (fVerbose) std::printf("  GIM: Gamma = %.6f (fim: %s)\n", gamma, fim.c_str());
    if (fSetReduction) fSetReduction(1.); // um passo que falhou no preditor pode ter deixado 1 + 1e-6 no material
    LoadSolution(fX); // a malha fica com o ultimo estado aceito (o mecanismo de colapso)
    return gamma;
}

template <class TMat>
bool UnsatAnalysis<TMat>::ApplyWeight(TPZFMatrix<STATE> &x)
{
    REAL lambda = 0., dl = 0.5;
    while (lambda < 1. && dl > 2.e-3 * std::max(lambda, REAL(0.01))) {
        const REAL lt = std::min(lambda + dl, REAL(1.));
        TPZFMatrix<STATE> xt(x);
        int it = 0;
        if (Newton(xt, lt, it)) {
            AcceptSolution();
            x = xt;
            lambda = lt;
            if (it <= 8) dl *= 1.5;
        } else {
            dl *= 0.5; // de novo do ultimo estado aceito (x e a memoria nao mudaram)
        }
    }
    return lambda >= 1.;
}

template <class TMat>
REAL UnsatAnalysis<TMat>::SRMArcLength()
{
    if (!fSetLoadFactor || !fSetReduction) {
        std::cerr << "UnsatAnalysis::SRMArcLength: falta SetLoadFactorFunction (o peso) ou SetStrengthReductionFunction (F)"
                  << std::endl;
        DebugStop();
    }
    TMat *mat = fModel.Material();
    const auto partida = *mat->GetMemory(); // estado de partida: sem tensoes (e as propriedades de cada ponto)
    const int64_t neq = fModel.CompMesh()->NEquations();
    fSRM = false;

    // 1. peso todo; F = 1 e, se o talude nao fica em pe (FS < 1), a resistencia dobrada
    TPZFMatrix<STATE> x(neq, 1, 0.);
    REAL F0 = 1.;
    fSetReduction(F0);
    while (!ApplyWeight(x)) {
        if (fVerbose) std::printf("  SRM: o talude nao fica em pe com o peso todo e F = %g\n", F0);
        F0 *= 0.5;
        if (F0 < 0.2) {
            std::cout << "  SRM: sem equilibrio com o peso todo nem com F = 0,25" << std::endl;
            fSetReduction(1.);
            return 0.;
        }
        *mat->GetMemory() = partida; // de novo do estado de partida
        x.Zero();
        fSetReduction(F0);
    }
    fXPeso = x;
    fX = x;
    fV = F0;
    fDu.Resize(0, 0);
    fDv = 0.;
    fFRef = FreeNorm(LoadVector(fX)); // referencia do residuo: a carga do peso, como no GIM

    // 2. caminho de F
    fSRM = true;
    if (fVerbose) {
        std::printf("  SRM arc-length\n    %5s %10s %10s %4s %11s\n", "passo", "F", "dF", "it", "arco l");
        std::printf("    %5d %10.5f %10s %4s %11s\n", 0, fV, "", "", "(peso)");
    }
    std::string fim;
    const REAL fs = ArcLengthPath(fim);
    if (fVerbose) std::printf("  SRM: FS = %.6f (fim: %s)\n", fs, fim.c_str());
    fSRM = false;
    fSetReduction(1.);
    fSetLoadFactor(1.);
    LoadSolution(fX); // a malha fica com o ultimo estado aceito (o mecanismo de colapso)
    return fs;
}

template <class TMat>
void UnsatAnalysis<TMat>::AcceptSolution()
{
    // com o fUpdateMem ligado, a montagem do residuo grava em cada ponto de integracao o estado da solucao carregada
    // (tensao, deformacao plastica, ...); depois o material volta a so ler a memoria (como o TPZPoroElastoPlasticUPAnalysis)
    TMat *mat = fModel.Material();
    mat->SetUpdateMem(true);
    AssembleResidual();
    mat->SetUpdateMem(false);
}

// Indicador plastico de cada elemento: o maior sqrt(J2) da deformacao plastica entre os pontos de integracao dele.
// J2 e o segundo invariante da parte desviadora, J2 = (1/2) e:e; sqrt(J2) e uma medida da distorcao plastica, e e
// grande onde o solo escoou (a superficie de ruptura no colapso)
template <class TMat>
std::vector<REAL> UnsatAnalysis<TMat>::PlasticIndicator()
{
    TPZMultiphysicsCompMesh *cmesh = fModel.CompMesh();
    TMat *mat = fModel.Material();
    std::vector<REAL> ind(cmesh->Reference()->NElements(), 0.);
    for (TPZCompEl *cel : cmesh->ElementVec()) {
        if (!cel || !cel->Reference() || cel->Reference()->MaterialId() != mat->Id()) continue;
        TPZManVector<int64_t> mem;
        cel->GetMemoryIndices(mem); // um item de memoria por ponto de integracao do elemento
        REAL &v = ind[cel->Reference()->Index()];
        for (int64_t im : mem) {
            if (im < 0) continue;
            TPZTensor<REAL> ep = mat->MemItem(im).m_elastoplastic_state.m_eps_p;
            // a memoria guarda o cisalhamento de engenharia (gamma = 2 eps); o J2 e do tensor: metade (como o Common)
            ep.XY() *= 0.5;
            ep.XZ() *= 0.5;
            ep.YZ() *= 0.5;
            v = std::max(v, std::sqrt(std::max<REAL>(ep.J2(), 0.)));
        }
    }
    return ind;
}

// Um VTK so, <nome>.vtk, com os elementos do dominio como sao (triangulos: tipo VTK 5) e:
//  - nos nos: os campos pedidos, pelos nomes do material, avaliados nos vertices de cada elemento;
//  - nas celulas: o indicador plastico (sqrtJ2_eps_p) e a ordem de u (ordem_u).
// Nao usa o TPZVTKGraphMesh do neopz: ele grava cada triangulo como um quadrilatero com o ultimo vertice repetido
// (a, b, c, c), e o ParaView ("Surface With Edges") esconde a diagonal a-c do quadrilatero, que ai e um lado do
// triangulo: cada triangulo aparecia com so dois lados. So os elementos computacionais entram (as folhas): numa malha
// refinada com Divide, os pais nao aparecem por cima dos filhos
template <class TMat>
void UnsatAnalysis<TMat>::PostProcess(const std::string &nome, const TPZVec<std::string> &scalars,
                                      const TPZVec<std::string> &vectors)
{
    TPZMultiphysicsCompMesh *cmesh = fModel.CompMesh();
    TPZGeoMesh *gmesh = fModel.GeoMesh();
    TMat *mat = fModel.Material();
    const int64_t nnos = gmesh->NNodes();

    // elementos do dominio
    std::vector<TPZCompEl *> els;
    for (TPZCompEl *cel : cmesh->ElementVec())
        if (cel && cel->Reference() && cel->Reference()->MaterialId() == mat->Id()) els.push_back(cel);

    // variaveis pedidas: indice no material e numero de componentes no arquivo (vetor do VTK: sempre 3)
    struct Var {
        std::string nome;
        int indice, ncomp;
        std::vector<STATE> valor; // por no: valor[no * ncomp + c]
    };
    std::vector<Var> vars;
    auto adiciona = [&](const std::string &n, int ncomp) {
        const int iv = mat->VariableIndex(n);
        if (iv < 0) {
            std::cerr << "UnsatAnalysis::PostProcess: o material nao conhece a variavel " << n << std::endl;
            return;
        }
        vars.push_back({n, iv, ncomp, std::vector<STATE>(nnos * ncomp, 0.)});
    };
    for (const std::string &n : scalars) adiciona(n, 1);
    for (const std::string &n : vectors) adiciona(n, 3);

    // campos nos vertices: Solution(qsi) do elemento multifisico, com qsi o vertice no elemento de referencia. u e p sao
    // continuos, entao qualquer elemento que tem o no da o mesmo valor
    for (TPZCompEl *cel : els) {
        TPZGeoEl *gel = cel->Reference();
        for (int i = 0; i < gel->NCornerNodes(); i++) {
            TPZManVector<REAL, 3> qsi(gel->Dimension(), 0.);
            gel->CenterPoint(i, qsi);
            const int64_t no = gel->NodeIndex(i);
            for (Var &v : vars) {
                TPZManVector<STATE, 9> sol;
                cel->Solution(qsi, v.indice, sol);
                for (int c = 0; c < v.ncomp && c < sol.size(); c++) v.valor[no * v.ncomp + c] = sol[c];
            }
        }
    }

    // indicador plastico e ordem de u de cada elemento
    const std::vector<REAL> ind = PlasticIndicator();
    auto ordemU = [](TPZCompEl *cel) {
        auto *mf = dynamic_cast<TPZMultiphysicsElement *>(cel);
        auto *sp = mf ? dynamic_cast<TPZInterpolationSpace *>(mf->Element(0)) : nullptr;
        return sp ? sp->MaxOrder() : 0;
    };

    // arquivo
    std::ofstream out(nome + ".vtk");
    out << "# vtk DataFile Version 3.0\n" << nome << "\nASCII\nDATASET UNSTRUCTURED_GRID\n";
    out << "POINTS " << nnos << " double\n";
    for (int64_t i = 0; i < nnos; i++) {
        TPZManVector<REAL, 3> x(3, 0.);
        gmesh->NodeVec()[i].GetCoordinates(x);
        out << x[0] << " " << x[1] << " " << x[2] << "\n";
    }
    int64_t tamanho = 0;
    for (TPZCompEl *cel : els) tamanho += 1 + cel->Reference()->NCornerNodes();
    out << "CELLS " << els.size() << " " << tamanho << "\n";
    for (TPZCompEl *cel : els) {
        TPZGeoEl *gel = cel->Reference();
        out << gel->NCornerNodes();
        for (int i = 0; i < gel->NCornerNodes(); i++) out << " " << gel->NodeIndex(i);
        out << "\n";
    }
    out << "CELL_TYPES " << els.size() << "\n";
    for (TPZCompEl *cel : els) out << TPZVTKGeoMesh::GetVTK_ElType(cel->Reference()) << "\n";
    out << "CELL_DATA " << els.size() << "\nSCALARS sqrtJ2_eps_p double 1\nLOOKUP_TABLE default\n";
    for (TPZCompEl *cel : els) out << ind[cel->Reference()->Index()] << "\n";
    out << "SCALARS ordem_u int 1\nLOOKUP_TABLE default\n";
    for (TPZCompEl *cel : els) out << ordemU(cel) << "\n";
    if (!vars.empty()) out << "POINT_DATA " << nnos << "\n";
    for (const Var &v : vars) {
        if (v.ncomp == 1) out << "SCALARS " << v.nome << " double 1\nLOOKUP_TABLE default\n";
        else out << "VECTORS " << v.nome << " double\n";
        for (int64_t i = 0; i < nnos; i++) {
            for (int c = 0; c < v.ncomp; c++) out << (c ? " " : "") << v.valor[i * v.ncomp + c];
            out << "\n";
        }
    }
}

// ---- materiais aceitos: os mesmos do fim do UnsatModel.cpp
using TMC = TPZPlasticStepVoigt<TPZYCMohrCoulombPV2, TPZElasticResponse>;
template class UnsatAnalysis<TPZMatPoroElastoPlasticUP<TMC, TPZElastoPlasticMem>>;                           // u-p, Mohr-Coulomb
template class UnsatAnalysis<TPZMatPoroElastoPlasticUP<TPZPlasticStepModifiedCamClay, TPZElastoPlasticMem>>; // u-p, Cam Clay
template class UnsatAnalysis<TPZMatUnsaturatedUPwPa<TMC, TPZUnsaturatedMem>>; // u-p_w-p_a, Mohr-Coulomb
// os mesmos com a regra de integracao do hp (HPMaterial.h)
template class UnsatAnalysis<HPMaterial<TPZMatPoroElastoPlasticUP<TMC, TPZElastoPlasticMem>>>;
template class UnsatAnalysis<HPMaterial<TPZMatPoroElastoPlasticUP<TPZPlasticStepModifiedCamClay, TPZElastoPlasticMem>>>;
template class UnsatAnalysis<HPMaterial<TPZMatUnsaturatedUPwPa<TMC, TPZUnsaturatedMem>>>;
