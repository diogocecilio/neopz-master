// Chuva de Rahardjo et al. (2007), deterministico: copia do main.cpp (a TriGMesh e o material nao saturado) com o
// transiente da chuva e o fator de seguranca do artigo (SRM, reducao da resistencia, como o Bishop do SLOPE/W).
//   talude: H = 10 m a 45 graus, Fig. 2 (3H de crista, 3H de pe, 3H de altura na esquerda e 2H na direita)
//   solos: os tres da Tabela 1 (f10,-4, f50,-5 e f100,-6), com c' = 10 kPa, phi' = phi^b = 26 graus, gamma = 20 kN/m3
//   chuva: I_r por 24 h; transiente ate 48 h (esqueleto rigido, ar passivo: o SEEP/W do artigo)
//   FS: de hora em hora, com refinamento adaptativo em h refeito do zero em cada hora
// Uso: SlopeRahardjo [solo] [I_r] [horas] [fator]
//   sem argumentos: os tres solos, cada um com I_r = 1 k_s (360, 36 e 3,6 mm/h; Figs. 4c e 5 do artigo)
//   solo: f10, f50 ou f100; I_r em mm/h (sem I_r: 1 k_s do solo); horas: duracao do transiente (padrao 48);
//   fator: srm (padrao, o do artigo) ou gim (aumento do peso: lambda multiplica o peso do solo e a poropressao
//   positiva, nao a succao; a pasta ganha _gim)
//   com a variavel de ambiente SRM_VERBOSE, o SRM imprime cada passo do caminho
// Saida: tabela FS(t) na tela e a pasta chuva_<solo>_<I_r> com os arquivos da animacao:
//   python3 ../LearningMohrCoulomb2/animacao.py chuva_f50_36
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include "UnsatAnalysis.h"
#include "UnsatModel.h"
#include "Plasticity/TPZElasticResponse.h"
#include "Plasticity/TPZMatUnsaturatedUPwPa.h"
#include "Plasticity/TPZPlasticStepVoigt.h"
#include "Plasticity/TPZYCMohrCoulombPV2.h"

#include "TPZGeoLinear.h"
#include "pzcmesh.h"
#include "pzgeotriangle.h"
#include "tpzgeoelrefpattern.h"
#include "TPZVTKGeoMesh.h"
#include "HPMaterial.h"
#include "AnimationFiles.h"
#include "RainTransient.h"
#include "PointLocator.h"
#include "TPZSpStructMatrix.h"
#include "pzstepsolver.h"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <memory>
#include <set>
#include <string>
#include <vector>
// material nao saturado com Mohr-Coulomb e a regra de integracao do hp (TUnsatHP)
using TMC = TPZPlasticStepVoigt<TPZYCMohrCoulombPV2, TPZElasticResponse>;
using TUnsat = TPZMatUnsaturatedUPwPa<TMC, TPZUnsaturatedMem>;
using TUnsatHP = HPMaterial<TUnsat>; // regra de integracao 2 (p_u - 1) em cada elemento
/// BC ids da TriGMesh: -1 base, -2 right, -3 topright, -4 topleft, -5 left, -6 face
enum EIds {
    matSolo = 1,
    base = -1,
    right = -2,
    topright = -3,
    topleft = -4,
    left=-5,
    face=-6
};

TPZGeoMesh *TriGMesh(const std::vector<std::vector<REAL>> &points, int ref);

// Tabela 1 de Rahardjo et al. (2007): o que muda de um solo para outro. O nome segue o artigo: f<a>,<expoente de k_s>
struct Solo {
    std::string nome;
    REAL a;  // parametro a da curva de retencao de Fredlund-Xing (kPa): ligado a entrada de ar
    REAL ks; // condutividade hidraulica saturada (m/s)
};
const Solo f10_4 = {"f10", 10., 1.e-4}, f50_5 = {"f50", 50., 1.e-5}, f100_6 = {"f100", 100., 1.e-6};

// O que e igual nos tres solos (kN, m, kPa). E e nu nao estao no artigo (equilibrio limite) e nao mudam o FS
struct SoloComum {
    REAL c = 10., phi = 26. * M_PI / 180.; // c' e phi'
    REAL gamma = 20., gammaw = 9.81;       // peso especifico do solo (acima e abaixo do nivel) e da agua
    REAL E = 2.e4, nu = 0.3;
};

// Solo, chuva e estado inicial para o RainTransient: Fredlund-Xing com n = m = 1 e theta_s = 0,45, k_w = k_s Theta^4,
// phi^b = 26 graus, nivel d'agua inicial 5 m abaixo do pe subindo 7 graus para dentro do macico e succao ate 75 kPa
// (os padroes do RainData); aqui so o solo, a chuva (I_r em mm/h por 24 h) e o pe do talude
RainData DadosChuva(const std::vector<std::vector<REAL>> &points, const Solo &s, REAL intensidade)
{
    RainData r;
    r.aFX = s.a;
    r.ks = s.ks;
    r.intensity = intensidade;
    r.duration = 24.;
    r.xToe = points[3][0];
    r.yToe = points[3][1];
    return r;
}

// Material do solo (o MaterialChuva do montecarlo.cpp): Mohr-Coulomb associado (psi = phi), retencao e permeabilidade
// do Rahardjo, ar passivo (p_a = 0) e chi = tan(phi^b) / tan(phi'): a parcela da succao (u_a - u_w) tan(phi^b) da Eq. 5
// do artigo. Sem acoplamento: no transiente so a agua percola, com o esqueleto parado; no fator as pressoes ficam
// congeladas
TUnsatHP *MaterialRahardjo(const RainData &r)
{
    const SoloComum s;
    TPZElasticResponse ER;
    ER.SetEngineeringData(s.E, s.nu);
    auto *mat = new TUnsatHP(matSolo);
    mat->SetPlasticModel(TMC(TPZYCMohrCoulombPV2(s.phi, s.phi, s.c, ER), ER));
    TPZUnsaturatedHydraulics h;
    h.fRetention = TPZUnsaturatedHydraulics::EFredlundXing; // S_w(p_c) com C(psi) = 1
    h.fA = r.aFX, h.fN = r.nFX, h.fM = r.mFX;
    h.fKrw = TPZUnsaturatedHydraulics::EPowerLaw; // k_w = k_s Theta^p
    h.fPower = r.power;
    h.fKra = TPZUnsaturatedHydraulics::EBrooksCorey; // o ar fica passivo (p_a = 0)
    mat->SetHydraulics(h);
    mat->SetBishop(TUnsat::EChiConstant, std::tan(r.phib) / std::tan(s.phi));
    mat->SetPorosity(r.n);
    mat->SetBiot(1., 0.);
    mat->SetWaterBulkModulus(r.Kw);
    mat->SetConductivity(r.ks, s.gammaw, 1.e-3 / 1.8e-5);
    mat->SetAir(1.2 * 9.81e-3, 101.325);
    mat->SetUnitWeight(s.gamma);
    mat->SetLumpedStorage(true);
    mat->SetHydroMechanicalCoupling(false);
    return mat;
}

// u: base engastada, laterais presas so em x
void InsertBC(UnsatModel<TUnsatHP> &model)
{
    int displacement=0;
    int xrestriction=0;
    int yrestriction=1;
    model.AddDirichletBC(right, TUnsat::EDirichletU, displacement, {xrestriction});
    model.AddDirichletBC(left, TUnsat::EDirichletU, displacement, {xrestriction});
    model.AddDirichletBC(base, TUnsat::EDirichletU, displacement, {xrestriction,yrestriction});
}

// Copia so com as folhas de uma malha refinada com Divide: o TPZUnsatAnalysis acha os nos da chuva pelos elementos de
// contorno da malha geometrica, e com os pais junto a chuva entraria duas vezes nos mesmos nos
TPZGeoMesh *MalhaPlana(TPZGeoMesh *gmesh)
{
    UnsatModel<TUnsatHP> tmp;
    tmp.SetGeoMesh(gmesh); // tmp fica dono da malha antiga e a apaga
    std::vector<int64_t> pai;
    return tmp.RefinedGeoMesh({}, pai);
}

// FS (SRM) com as pressoes pw (uma por no da malha do transiente) levadas para a malha do fator, com o refinamento
// adaptativo do SolveDrawDown feito do zero: TriGMesh(points, 1) e nref niveis em h (marcados: sqrt(J2(eps_p)) >= 0,1 do
// maior). O modelo e a analise do ultimo nivel ficam em mF e aF (o estado do fim do SRM, para a animacao); niveis: o FS
// de cada nivel. gim: o fator Gamma do GIM no lugar do FS do SRM
REAL FatorAdaptativo(const std::vector<std::vector<REAL>> &points, const RainData &r, const PointLocator &loc,
                     const std::vector<REAL> &pw, int nref, std::unique_ptr<UnsatModel<TUnsatHP>> &mF,
                     std::unique_ptr<UnsatAnalysis<TUnsatHP>> &aF, std::vector<REAL> &niveis, bool gim = false)
{
    TPZGeoMesh *gmesh = TriGMesh(points, 1);
    REAL fs = 0.;
    niveis.clear();
    for (int nivel = 0;; nivel++) {
        aF.reset(); // a analise antes do modelo (ela aponta para ele)
        mF = std::make_unique<UnsatModel<TUnsatHP>>();
        mF->SetGeoMesh(gmesh);
        InsertBC(*mF);
        mF->SetMaterial(MaterialRahardjo(r));
        mF->CreateCompMesh({{2, 2}, {1, 1}, {1, 1}});
        aF = std::make_unique<UnsatAnalysis<TUnsatHP>>(*mF);
        aF->SetVerbose(std::getenv("SRM_VERBOSE") != nullptr);
        TUnsatHP *mat = mF->Material();
        aF->SetLoadFactorFunction([mat](REAL l) { mat->SetGravityFactor(l); }, false); // o peso (lambda de 0 a 1)
        aF->SetStrengthReductionFunction([mat](REAL F) { mat->PlasticModel().SetStrengthReductionFactor(F); });
        aF->SetFrozenPressures([&](int espaco, const TPZVec<REAL> &x) {
            return espaco == 1 ? loc.Interpolate(x[0], x[1], pw) : 0.; // p_w do transiente; p_a = 0
        });
        if (gim) mat->SetGravityScalesSuction(false); // a succao nao e peso: fica fixa
        fs = gim ? aF->GIMArcLength() : aF->SRMArcLength();
        niveis.push_back(fs);
        if (nivel == nref) break;

        // marca os triangulos com sqrt(J2(eps_p)) >= 0,1 do maior e refina (h)
        std::vector<REAL> ind = aF->PlasticIndicator();
        REAL vmax = *std::max_element(ind.begin(), ind.end());
        std::set<int64_t> marcados;
        for (TPZGeoEl *gel : gmesh->ElementVec())
            if (gel && gel->Dimension() == 2 && ind[gel->Index()] >= 0.1 * vmax) marcados.insert(gel->Index());
        std::vector<int64_t> pai;
        gmesh = mF->RefinedGeoMesh(marcados, pai);
    }
    return fs;
}

// Transiente da chuva de Rahardjo no talude points: solo s sob I_r = intensidade (mm/h) por 24 h, ate horas (48: 24 h
// de chuva e 24 h depois). A cada hora, o FS (SRM) com as pressoes do instante. Grava a pasta chuva_<solo>_<I_r> para o
// animacao.py. gim: o fator e o Gamma do GIM (pasta chuva_<solo>_<I_r>_gim)
void TransientRain(const std::vector<std::vector<REAL>> &points, const Solo &s, REAL intensidade, REAL horas = 48.,
                   bool gim = false)
{
    const int refTransiente = 3; // malha do transiente: TriGMesh com 3 divisoes uniformes (lado de H / 8 = 1,25 m)
    const int nref = 4;          // niveis de refinamento do fator depois da TriGMesh(points, 1)
    const REAL dtGravar = 0.25;  // intervalo entre os instantes gravados (h); o FS, de hora em hora
    const RainData r = DadosChuva(points, s, intensidade);
    const SoloComum comum;

    using Clock = std::chrono::steady_clock;
    const auto t0 = Clock::now();
    std::printf("\n==== Rahardjo: solo %s (a = %g kPa, k_s = %g m/s), I_r = %g mm/h por %g h, transiente ate %g h, "
                "fator %s\n",
                s.nome.c_str(), s.a, s.ks, intensidade, r.duration, horas, gim ? "GIM" : "SRM");

    // transiente: so a agua (sem acoplamento), na malha uniforme sem hierarquia; u preso como no fator
    UnsatModel<TUnsatHP> mT;
    mT.SetGeoMesh(MalhaPlana(TriGMesh(points, refTransiente)));
    InsertBC(mT);
    mT.SetMaterial(MaterialRahardjo(r));
    mT.CreateCompMesh({{2, 2}, {1, 1}, {1, 1}});
    TPZUnsatAnalysis aT(mT.CompMesh(), {{base, {true, true}}, {right, {true, false}}, {left, {true, false}}});
    TPZSpStructMatrix<STATE> strmat(mT.CompMesh());
    strmat.SetNumThreads(std::getenv("PZ_THREADS") ? std::atoi(std::getenv("PZ_THREADS")) : 8); // como na UnsatAnalysis
    aT.SetStructuralMatrix(strmat);
    TPZStepSolver<STATE> solver;
    solver.SetDirect(ELU);
    aT.SetSolver(solver);
    RainTransient chuva(aT, r, comum.gammaw, {left, right}, {topleft, face, topright});
    chuva.InitialState();
    const PointLocator loc(mT.GeoMesh());

    char pasta[64];
    std::snprintf(pasta, sizeof pasta, "chuva_%s_%g%s", s.nome.c_str(), intensidade, gim ? "_gim" : "");
    AnimationFiles<TUnsatHP> anim(pasta, mT);
    anim.RainSummary(intensidade, r.duration);

    // uma linha da tabela final por hora
    struct Linha {
        REAL t, fs;
        int64_t neq;
        double segundos;
    };
    std::vector<Linha> tabela;
    double sTransiente = 0.;
    const int n = int(std::lround(horas / dtGravar));
    for (int k = 0; k <= n; k++) {
        const REAL t = k * dtGravar;
        auto tt = Clock::now();
        const bool ok = chuva.AdvanceTo(t);
        sTransiente += std::chrono::duration<double>(Clock::now() - tt).count();
        if (!ok) {
            std::cout << "o transiente parou em t = " << chuva.Time() << " h" << std::endl;
            break;
        }
        const std::vector<REAL> pw = chuva.NodalWaterPressure();
        anim.Snapshot(t, pw, mT.Material()->Hydraulics());
        if (std::abs(t - std::round(t)) > 1.e-9) continue; // o FS de hora em hora

        const auto tf = Clock::now();
        std::unique_ptr<UnsatModel<TUnsatHP>> mF;
        std::unique_ptr<UnsatAnalysis<TUnsatHP>> aF;
        std::vector<REAL> niveis;
        const REAL fs = FatorAdaptativo(points, r, loc, pw, nref, mF, aF, niveis, gim);
        anim.FactorOfSafety(t, fs, *mF, *aF);
        Linha l{t, fs, mF->CompMesh()->NEquations(), std::chrono::duration<double>(Clock::now() - tf).count()};
        tabela.push_back(l);
        const auto [pmin, pmax] = chuva.PressureRange();
        std::string sn; // FS de cada nivel
        for (REAL v : niveis) sn += " " + std::to_string(v).substr(0, 6);
        std::printf("  t = %4.0f h: FS = %.4f | p_w de %.1f a %.1f kPa, %zu nos na face | fator: %lld equacoes, %.1f s |"
                    " niveis:%s\n",
                    t, fs, pmin, pmax, chuva.FaceNodes(), (long long)l.neq, l.segundos, sn.c_str());
        std::fflush(stdout);
    }

    // ---- tabela final
    std::printf("\n  solo %s, I_r = %g mm/h\n   t (h)       FS  equacoes  tempo (s)\n", s.nome.c_str(), intensidade);
    const Linha *menor = nullptr;
    for (const Linha &l : tabela) {
        std::printf("  %5.0f  %7.4f  %8lld  %9.1f\n", l.t, l.fs, (long long)l.neq, l.segundos);
        if (!menor || l.fs < menor->fs) menor = &l;
    }
    if (menor) std::printf("  FS minimo: %.4f em t = %.0f h (Rahardjo, Tabela 7, 45 graus: F_ini = 2,118)\n", menor->fs, menor->t);
    std::printf("  tempo: transiente %.1f s, total %.1f s; pasta %s\n", sTransiente,
                std::chrono::duration<double>(Clock::now() - t0).count(), pasta);
}

int main(int argc, char *argv[])
{
    // cantos do talude no sentido anti-horario: base esquerda, base direita, terreno do pe a direita, pe, crista, crista
    // a esquerda. Fig. 2 de Rahardjo et al. (2007) com H = 10 m a 45 graus
    const std::vector<std::vector<REAL>> rahardjo = {{0, 0}, {70, 0}, {70, 20}, {40, 20}, {30, 30}, {0, 30}};
    const std::vector<Solo> solos = {f10_4, f50_5, f100_6};
    if (argc == 1) { // os tres solos com I_r = 1 k_s (k_s em m/s x 3,6e6 = mm/h)
        for (const Solo &s : solos) TransientRain(rahardjo, s, s.ks * 3.6e6);
        return 0;
    }
    for (const Solo &s : solos)
        if (s.nome == argv[1]) {
            TransientRain(rahardjo, s, argc > 2 ? std::atof(argv[2]) : s.ks * 3.6e6, argc > 3 ? std::atof(argv[3]) : 48.,
                          argc > 4 && std::string(argv[4]) == "gim");
            return 0;
        }
    std::fprintf(stderr, "solo desconhecido: %s (use f10, f50 ou f100)\n", argv[1]);
    return 1;
}

/// Talude de 6 cantos no sentido anti-horario: base esquerda, base direita, terreno do pe a direita, pe, crista, crista
/// a esquerda. Ex.: {{0, 0}, {70, 0}, {70, 30}, {40, 30}, {30, 40}, {0, 40}} (H = 10 m a 45 graus; a malha antiga) ou o do
/// Ceron {{0, 0}, {25, 0}, {25, 5}, {15, 5}, {10, 10}, {0, 10}} (H = 5 m). Malha estruturada com passo perto da altura H
/// do talude: cada trecho horizontal e vertical e dividido em round(comprimento / H) partes iguais, cada retangulo em
/// dois triangulos (diagonal do canto inferior direito ao superior esquerdo) e a cunha da face (pe, crista e o ponto
/// abaixo da crista) num triangulo so. BC ids: -1 base, -2 right, -3 top right, -4 top left, -5 left, -6 face.
/// ref = refinamentos uniformes.
TPZGeoMesh *TriGMesh(const std::vector<std::vector<REAL>> &points, int ref) {
    if (points.size() != 6) DebugStop();
    const REAL x0 = points[0][0], y0 = points[0][1], L = points[1][0], yPe = points[2][1];
    const REAL xPe = points[3][0], xCr = points[4][0], yCr = points[4][1];
    // a forma: base e terrenos horizontais, laterais verticais, crista a esquerda e acima do pe
    if (points[1][1] != y0 || points[2][0] != L || points[3][1] != yPe || points[5][0] != x0 || points[5][1] != yCr ||
        !(x0 < xCr && xCr < xPe && xPe < L && y0 < yPe && yPe < yCr))
        DebugStop();
    const REAL H = yCr - yPe; // altura do talude: o passo da malha
    // pontos de a ate b em partes iguais de umas H (pelo menos uma)
    auto divide = [H](REAL a, REAL b) {
        const int n = std::max(1, int(std::lround((b - a) / H)));
        std::vector<REAL> v(n + 1);
        for (int k = 0; k <= n; k++) v[k] = a + (b - a) * k / n;
        return v;
    };
    std::vector<REAL> xs = divide(x0, xCr); // ate a crista; a cunha (crista ao pe) e uma divisao so
    const int iCr = int(xs.size()) - 1, iPe = iCr + 1;
    for (REAL x : divide(xPe, L)) xs.push_back(x); // do pe ate a lateral direita
    const std::vector<REAL> ys = divide(y0, yPe);  // ate o pe; a camada da crista (pe a crista) e uma divisao so
    const int nx = int(xs.size()) - 1, ny = int(ys.size()) - 1;

    // nos: as linhas ate o pe com todas as colunas e a linha da crista com as colunas ate a crista
    auto no = [nx](int i, int j) { return int64_t(j) * (nx + 1) + i; };
    auto topo = [nx, ny](int i) { return int64_t(ny + 1) * (nx + 1) + i; };
    auto *gmesh = new TPZGeoMesh();
    gmesh->SetDimension(2);
    gmesh->NodeVec().Resize(topo(iCr) + 1);
    auto poeNo = [gmesh](int64_t id, REAL x, REAL y) {
        TPZManVector<REAL, 3> c = {x, y, 0.};
        gmesh->NodeVec()[id] = TPZGeoNode(id, c, *gmesh);
    };
    for (int j = 0; j <= ny; j++)
        for (int i = 0; i <= nx; i++) poeNo(no(i, j), xs[i], ys[j]);
    for (int i = 0; i <= iCr; i++) poeNo(topo(i), xs[i], yCr);

    // triangulos: cada retangulo (a, b, c, e anti-horario a partir do canto inferior esquerdo) em {a, b, e} e {b, c, e}
    auto triangulo = [gmesh](int64_t a, int64_t b, int64_t c) {
        TPZManVector<int64_t, 3> nodes = {a, b, c};
        new TPZGeoElRefPattern<pzgeom::TPZGeoTriangle>(nodes, matSolo, *gmesh);
    };
    for (int j = 0; j < ny; j++)
        for (int i = 0; i < nx; i++) {
            triangulo(no(i, j), no(i + 1, j), no(i, j + 1));
            triangulo(no(i + 1, j), no(i + 1, j + 1), no(i, j + 1));
        }
    for (int i = 0; i < iCr; i++) { // camada da crista
        triangulo(no(i, ny), no(i + 1, ny), topo(i));
        triangulo(no(i + 1, ny), topo(i + 1), topo(i));
    }
    triangulo(no(iCr, ny), no(iPe, ny), topo(iCr)); // cunha da face

    // contornos no sentido anti-horario
    auto linha = [gmesh](int64_t a, int64_t b, int id) {
        TPZManVector<int64_t, 2> nodes = {a, b};
        new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodes, id, *gmesh);
    };
    for (int i = 0; i < nx; i++) linha(no(i, 0), no(i + 1, 0), base);
    for (int j = 0; j < ny; j++) linha(no(nx, j), no(nx, j + 1), right);
    for (int i = nx; i > iPe; i--) linha(no(i, ny), no(i - 1, ny), topright);
    for (int i = iCr; i > 0; i--) linha(topo(i), topo(i - 1), topleft);
    linha(topo(0), no(0, ny), left);
    for (int j = ny; j > 0; j--) linha(no(0, j), no(0, j - 1), left);
    linha(no(iPe, ny), topo(iCr), face);

    gmesh->BuildConnectivity();
    for (int d = 0; d < ref; d++) {
        const int64_t nel = gmesh->NElements();
        TPZManVector<TPZGeoEl *> sub;
        for (int64_t iel = 0; iel < nel; iel++)
            if (!gmesh->Element(iel)->HasSubElement()) gmesh->Element(iel)->Divide(sub);
    }

    std::ofstream meshFile("geo.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, meshFile, true);

    return gmesh;
}
