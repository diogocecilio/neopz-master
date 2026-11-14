
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include "UnsatAnalysis.h"
#include "UnsatModel.h"
#include "Plasticity/TPZElasticResponse.h"
#include "Plasticity/TPZMatPoroElastoPlasticUP.h"
#include "Plasticity/TPZMatUnsaturatedUPwPa.h"
#include "Plasticity/TPZPlasticStepVoigt.h"
#include "Plasticity/TPZYCMohrCoulombPV2.h"

#include "TPZGeoLinear.h"
#include "pzcmesh.h"
#include "pzgeotriangle.h"
#include "tpzgeoelrefpattern.h"
#include "TPZVTKGeoMesh.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZNullMaterial.h"
#include "HPMaterial.h"
#include "AnimationFiles.h"
#include "DrawdownTransient.h"
#include "PointLocator.h"
#include "TPZSpStructMatrix.h"
#include "pzstepsolver.h"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <set>
#include <string>
#include <vector>
// materiais: u-p saturado (TMat) e nao saturado (TUnsat) com Mohr-Coulomb; o main usa o nao saturado com a regra de
// integracao do hp (TUnsatHP)
using TMC = TPZPlasticStepVoigt<TPZYCMohrCoulombPV2, TPZElasticResponse>;
using TMat = TPZMatPoroElastoPlasticUP<TMC, TPZElastoPlasticMem>;
using TUnsat = TPZMatUnsaturatedUPwPa<TMC, TPZUnsaturatedMem>;
using TUnsatHP = HPMaterial<TUnsat>; // regra de integracao 2 (p_u - 1) em cada elemento
/// Slope 70 x 40 m, height 10 m at 45 degrees. BC ids: -1 base, -2 right, -3 topright,
/// -4 topleft, -5 left, -6 face. ref = uniform refinements.
enum EIds {
    matSolo = 1,
    base = -1,
    right = -2,
    topright = -3,
    topleft = -4,
    left=-5,
    face=-6
};

#include "LayerGMesh.h"
TPZGeoMesh *Camadas(REAL h, bool quad);

TPZGeoMesh *TriGMesh(const std::vector<std::vector<REAL>> &points, int ref);

// Poropressao hidrostatica de um NA na cota nivel (kPa): gamma_w (nivel - y) abaixo do NA e 0 acima.
// So depende de y; o contorno so escolhe o nivel. v2[0] = p_w no ponto x (nivel e gammaw ficam guardados no lambda)
ForcingFunctionBCType<STATE> HydrostaticPressure(REAL nivel, REAL gammaw)
{
    return [nivel, gammaw](const TPZVec<REAL> &x, TPZVec<STATE> &v2, TPZFMatrix<STATE> &) {
        v2[0] = gammaw * std::max(REAL(0), nivel - x[1]); // x[1] = y
    };
}

/// Parameters of a soil (units kN, m, s, kPa); defaults: chi = S_w (notes), van Genuchten - Mualem
struct UnsatParams {
    REAL E = 20000., nu = 0.3, phi = 30. * M_PI / 180., c = 10.;
    REAL gamma = 20.;            ///< constant unit weight (kN/m3)
    bool mixture = false;        ///< unit weight of the mixture (1-n) gamma_s + n S_w gamma_w
    REAL gammas =30.;           ///< unit weight of the grains (mixture)
    TUnsat::EBishop bishop = TUnsat::EChiSaturation;
    REAL chi = 1.;               ///< constant chi (bishop = EChiConstant)
    REAL n = 0.45, Kw = 2.2e6, gw = 9.81, ga = 1.2 * 9.81e-3, muRatio = 1.e-3 / 1.8e-5, patm = 101.325;
    REAL ks = 1.e-5;             ///< saturated hydraulic conductivity (m/s)
    REAL alpha = 1., Ks = 0.;    ///< Biot coefficient and bulk modulus of the grains (<= 0: incompressible)
    bool coupled = true;         ///< hydro-mechanical coupling of the mass balances
    bool lumped = true;          ///< lumped storage
    std::array<REAL, 3> g = {0., -1., 0.}; ///< unit vector of gravity
    TPZUnsaturatedHydraulics hydro;
    UnsatParams() {
        hydro.fRetention = TPZUnsaturatedHydraulics::EVanGenuchten;
        hydro.fKrw = TPZUnsaturatedHydraulics::EMualem;
        hydro.fKra = TPZUnsaturatedHydraulics::EMualemAir;
    }
};
// Material TUnsat com Mohr-Coulomb associado (psi = phi) a partir dos parametros (copia de UnsatModel.h), com a regra
// de integracao do hp (HPMaterial)
TUnsatHP *UnsatMaterial(int id, const UnsatParams &p)
{
    TPZElasticResponse ER;
    ER.SetEngineeringData(p.E, p.nu);
    auto *mat = new TUnsatHP(id);
    mat->SetPlasticModel(TMC(TPZYCMohrCoulombPV2(p.phi, p.phi, p.c, ER), ER));
    mat->SetHydraulics(p.hydro);
    mat->SetPorosity(p.n);
    mat->SetBiot(p.alpha, p.Ks);
    mat->SetWaterBulkModulus(p.Kw);
    mat->SetConductivity(p.ks, p.gw, p.muRatio);
    mat->SetAir(p.ga, p.patm);
    if (p.mixture) mat->SetUnitWeightMixture(p.gammas);
    else mat->SetUnitWeight(p.gamma);
    mat->SetGravityDirection(TPZManVector<REAL, 3>({p.g[0], p.g[1], p.g[2]}));
    mat->SetBishop(p.bishop, p.chi);
    mat->SetLumpedStorage(p.lumped);
    mat->SetHydroMechanicalCoupling(p.coupled);
    // regra de integracao: a do HPMaterial, ordem 2 (p_u - 1) em cada elemento (3 pontos em P2, 6 em P3)
    return mat;
}

// points: os 6 cantos do talude (veja TriGMesh)
void SolveDrawDown(const std::vector<std::vector<REAL>> &points)
{
    const int nref = 5; // niveis depois do nivel 0
    const int pmax = 2; // maior ordem de u (2: so h)
    TPZGeoMesh *gmesh = TriGMesh(points, 1);

    std::vector<int> ordens(gmesh->NElements(), 2); // ordem de u de cada elemento da malha geo do nivel

    // uma linha da tabela final por nivel
    struct Linha {
        int nivel;
        int64_t ntri, neq;
        std::map<int, int> ordens; // ordem de u -> triangulos
        REAL gamma;
        double segundos;
    };

    std::vector<Linha> tabela;
    using Clock = std::chrono::steady_clock;
    const auto t0 = Clock::now();

    for (int nivel = 0;; nivel++) {
        const auto tn = Clock::now();
        UnsatModel<TUnsatHP> model;
        model.SetGeoMesh(gmesh); // malha geo do nivel

        // u: base engastada, laterais presas so em x
        int displacement=0;
        int xrestriction=0;
        int yrestriction=1;
        model.AddDirichletBC(right, TUnsat::EDirichletU, displacement, {xrestriction});
        model.AddDirichletBC(left, TUnsat::EDirichletU, displacement, {xrestriction});
        model.AddDirichletBC(base, TUnsat::EDirichletU, displacement, {xrestriction,yrestriction});

        // Rebaixamento rapido do Ceron (Eq. 1 de Vargas Ceron et al. 2025): antes, NA externo na crista; com h_w = H
        // (total) ele desce ate o pe. A base fica sem BC de pressao: fluxo nulo (natural).
        // Um id pode ter BC de u e de p_w (left, right)
        int waterpressure=1;
        REAL gammaw=9.81;      // kN/m3
        REAL yCrista=points[4][1];            // NA antes do rebaixamento: na crista
        REAL hw=points[4][1]-points[3][1];    // rebaixamento h_w = H: o NA desce ate o pe
        REAL nivelNA=yCrista-hw; // NA externo depois do rebaixamento
        model.AddDirichletBC(left, TUnsat::EDirichletPw, waterpressure, {0}, HydrostaticPressure(yCrista, gammaw)); // macico
        for (int id : {right, topright, face, topleft}) // lado do pe e superficie (0 acima do NA)
            model.AddDirichletBC(id, TUnsat::EDirichletPw, waterpressure, {0}, HydrostaticPressure(nivelNA, gammaw));

            // material do solo (do problema) e a malha comp: u com a ordem de cada elemento (P2 ou mais), p_w e p_a P1
        UnsatParams p;
        model.SetMaterial(UnsatMaterial(matSolo, p));
        model.CreateCompMesh({{2, 2}, {1, 1}, {1, 1}}, ordens); // malha comp

        UnsatAnalysis<TUnsatHP> an(model); // analysis
        an.SolveSeepage();                 // solve
        TUnsatHP *mat = model.Material();
        an.SetLoadFactorFunction([mat](REAL l) { mat->SetGravityFactor(l); }, false);
        an.SetStrengthReductionFunction([mat](REAL F) { mat->PlasticModel().SetStrengthReductionFactor(F); });
        REAL gamma = an.GIMArcLength();
        an.PostProcess("gim_nivel" + std::to_string(nivel),
                           {"WaterPressure", "AirPressure", "Suction", "Saturation", "RelativePermeabilityWater",
                               "BishopPressure", "DisplacementX", "DisplacementY"},
                           {"Displacement"});

        Linha l{nivel, 0, model.CompMesh()->NEquations(), {}, gamma, 0.};

        for (TPZGeoEl *gel : gmesh->ElementVec())
        {
            if (gel && gel->Dimension() == 2) l.ntri++, l.ordens[ordens[gel->Index()]]++;
        }

        l.segundos = std::chrono::duration<double>(Clock::now() - tn).count();
        tabela.push_back(l);

        std::cout << "nivel " << nivel << ": Gamma = " << gamma << ", " << l.neq << " equacoes, " << l.segundos << " s"
        << std::endl;

        if (nivel == nref) break;

        // marca os triangulos com sqrt(J2(eps_p)) >= 0,1 do maior; p: ordem + 1 nos marcados (ate pmax)
        std::vector<REAL> ind = an.PlasticIndicator();
        REAL vmax = *std::max_element(ind.begin(), ind.end());
        std::set<int64_t> marcados;
        for (TPZGeoEl *gel : gmesh->ElementVec())
        {
            if (gel && gel->Dimension() == 2 && ind[gel->Index()] >= 0.1 * vmax) marcados.insert(gel->Index());
        }

        for (int64_t i : marcados) ordens[i] = std::min(ordens[i] + 1, pmax);

        // h: malha nova; cada elemento novo herda a ordem do elemento de onde veio (pai)
        std::vector<int64_t> pai;
        gmesh = model.RefinedGeoMesh(marcados, pai);
        std::ofstream meshFile("geoc.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, meshFile, true);
        std::vector<int> novas(pai.size());
        for (size_t j = 0; j < pai.size(); j++) novas[j] = ordens[pai[j]];
        ordens = novas;
    }

    // ---- tabela final
    std::printf("\n nivel triangulos equacoes  ordens de u         Gamma  tempo (s)\n");
    for (const Linha &l : tabela)
    {
        std::string o;
        for (const auto &[ordem, n] : l.ordens) o += "P" + std::to_string(ordem) + ":" + std::to_string(n) + " ";
        std::printf(" %5d %10lld %8lld  %-17s %8.5f %9.1f\n", l.nivel, (long long)l.ntri, (long long)l.neq, o.c_str(),
                    l.gamma, l.segundos);
    }
    std::printf(" total: %.1f s\n", std::chrono::duration<double>(Clock::now() - t0).count());
}

void InsertBC(UnsatModel<TUnsatHP> &model,const std::vector<std::vector<REAL>> &points)
{

    int displacement=0;
    int xrestriction=0;
    int yrestriction=1;
    model.AddDirichletBC(right, TUnsat::EDirichletU, displacement, {xrestriction});
    model.AddDirichletBC(left, TUnsat::EDirichletU, displacement, {xrestriction});
    model.AddDirichletBC(base, TUnsat::EDirichletU, displacement, {xrestriction,yrestriction});

    int waterpressure=1;
    REAL gammaw=9.81;      // kN/m3
    REAL yCrista=points[4][1];            // NA antes do rebaixamento: na crista
    REAL hw=points[4][1]-points[3][1];    // rebaixamento h_w = H: o NA desce ate o pe
    REAL nivelNA=yCrista-hw; // NA externo depois do rebaixamento
    model.AddDirichletBC(left, TUnsat::EDirichletPw, waterpressure, {0}, HydrostaticPressure(yCrista, gammaw)); // macico
    for (int id : {right, topright, face, topleft})
    {// lado do pe e superficie (0 acima do NA)
        model.AddDirichletBC(id, TUnsat::EDirichletPw, waterpressure, {0}, HydrostaticPressure(nivelNA, gammaw));
    }
}

// points: os 6 cantos do talude (veja TriGMesh)
void SolveDrawDownTransient(const std::vector<std::vector<REAL>> &points)
{
    const int nref = 5; // niveis depois do nivel 0
    const int pmax = 2; // maior ordem de u (2: so h)
    TPZGeoMesh *gmesh = TriGMesh(points, 1);
    std::vector<int> ordens(gmesh->NElements(), 2);

    TPZVec<std::string> scalnames ={"WaterPressure", "AirPressure", "Suction", "Saturation", "RelativePermeabilityWater","BishopPressure"};
    TPZVec<std::string> vecnames ={"Displacement"};
    using Clock = std::chrono::steady_clock;
    const auto t0 = Clock::now();

    for (int nivel = 0;; nivel++) {
        const auto tn = Clock::now();

        UnsatModel<TUnsatHP> model;

        model.SetGeoMesh(gmesh); // malha geo do nivel

        InsertBC(model,points);

        // material do solo (do problema) e a malha comp: u com a ordem de cada elemento (P2 ou mais), p_w e p_a P1
        UnsatParams p;
        p.ks = 1.e-6;
        model.SetMaterial(UnsatMaterial(matSolo, p));
        model.CreateCompMesh({{2, 2}, {1, 1}, {1, 1}}, ordens); // malha comp

        UnsatAnalysis<TUnsatHP> an(model); // analysis

        an.SolveSeepage();                 // solve
        TUnsatHP *mat = model.Material();
        an.SetLoadFactorFunction([mat](REAL l) { mat->SetGravityFactor(l); }, false);
        an.SetStrengthReductionFunction([mat](REAL F) { mat->PlasticModel().SetStrengthReductionFactor(F); });

        REAL gamma = an.GIMArcLength();

        an.PostProcess("gim_nivel" + std::to_string(nivel),scalnames,vecnames);


        std::cout << "nivel " << nivel << ": Gamma = " << gamma << ", " << " time: " << std::chrono::duration<double>(Clock::now() - tn).count()
            << std::endl;

        if (nivel == nref) break;

            // marca os triangulos com sqrt(J2(eps_p)) >= 0,1 do maior; p: ordem + 1 nos marcados (ate pmax)
        std::vector<REAL> ind = an.PlasticIndicator();
        REAL vmax = *std::max_element(ind.begin(), ind.end());

        std::set<int64_t> marcados;
        for (TPZGeoEl *gel : gmesh->ElementVec())
        {
            if (gel && gel->Dimension() == 2 && ind[gel->Index()] >= 0.1 * vmax) marcados.insert(gel->Index());
        }

        // h: malha nova; cada elemento novo herda a ordem do elemento de onde veio (pai)
        std::vector<int64_t> pai;
        gmesh = model.RefinedGeoMesh(marcados, pai);
        std::ofstream meshFile("geoc.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, meshFile, true);
    }

    std::printf(" total: %.1f s\n", std::chrono::duration<double>(Clock::now() - t0).count());
}


// Malha refinada para o mecanismo do rebaixamento total (o estado mais critico): o laco do SolveDrawDown, sem tabela.
// Devolve duas copias iguais (mesma numeracao de nos): uma para o transiente e outra para o fator de seguranca
std::pair<TPZGeoMesh *, TPZGeoMesh *> MalhaDoRebaixamento(const std::vector<std::vector<REAL>> &points, int nref)
{
    TPZGeoMesh *gmesh = TriGMesh(points, 1);
    for (int nivel = 0;; nivel++) {
        UnsatModel<TUnsatHP> model;
        model.SetGeoMesh(gmesh);
        InsertBC(model, points);
        UnsatParams p;
        model.SetMaterial(UnsatMaterial(matSolo, p));
        model.CreateCompMesh({{2, 2}, {1, 1}, {1, 1}});
        UnsatAnalysis<TUnsatHP> an(model);
        an.SetVerbose(false);
        an.SolveSeepage();
        TUnsatHP *mat = model.Material();
        an.SetLoadFactorFunction([mat](REAL l) { mat->SetGravityFactor(l); }, false);
        std::cout << "malha, nivel " << nivel << ": Gamma = " << an.GIMArcLength() << std::endl;
        std::vector<int64_t> pai;
        if (nivel == nref) return {model.RefinedGeoMesh({}, pai), model.RefinedGeoMesh({}, pai)};
        std::vector<REAL> ind = an.PlasticIndicator();
        REAL vmax = *std::max_element(ind.begin(), ind.end());
        std::set<int64_t> marcados;
        for (TPZGeoEl *gel : gmesh->ElementVec())
            if (gel && gel->Dimension() == 2 && ind[gel->Index()] >= 0.1 * vmax) marcados.insert(gel->Index());
        gmesh = model.RefinedGeoMesh(marcados, pai);
    }
}

// Transiente do rebaixamento: o reservatorio desce em linha reta da crista ate o pe em horasRebaixamento e fica parado
// por horasDepois; esqueleto rigido (so a agua percola). De hora em hora, o fator Gamma (GIM) com as pressoes do
// instante. Grava a pasta "rebaixamento" para o animacao_rebaixamento.py
void TransientDrawDown(const std::vector<std::vector<REAL>> &points)
{
    const int nref = 3;                 // niveis de refinamento da malha (a mesma no transiente e no fator)
    const REAL horasRebaixamento = 0.1; // duracao da descida do reservatorio (h)
    const REAL horasDepois = 24.;       // transiente depois do rebaixamento (h)
    const REAL dtGravar = 0.1;          // intervalo entre os instantes gravados (h); o fator, de hora em hora e no fim da descida

    auto [gT, gF] = MalhaDoRebaixamento(points, nref);
    UnsatParams p;
    p.ks = 1.e-6; // k do Ceron (m/s): antes da copia pT, para valer no transiente

    // transiente: so a agua (sem acoplamento); u preso como no fator
    UnsatModel<TUnsatHP> mT;
    mT.SetGeoMesh(gT);
    mT.AddDirichletBC(base, TUnsat::EDirichletU, 0, {0, 1});
    mT.AddDirichletBC(left, TUnsat::EDirichletU, 0, {0});
    mT.AddDirichletBC(right, TUnsat::EDirichletU, 0, {0});
    UnsatParams pT = p;
    pT.coupled = false;
    mT.SetMaterial(UnsatMaterial(matSolo, pT));
    mT.CreateCompMesh({{2, 2}, {1, 1}, {1, 1}});
    TPZUnsatAnalysis aT(mT.CompMesh(), {{base, {true, true}}, {right, {true, false}}, {left, {true, false}}});
    TPZSpStructMatrix<STATE> strmat(mT.CompMesh());
    strmat.SetNumThreads(std::getenv("PZ_THREADS") ? std::atoi(std::getenv("PZ_THREADS")) : 8); // como na UnsatAnalysis
    aT.SetStructuralMatrix(strmat);
    TPZStepSolver<STATE> solver;
    solver.SetDirect(ELU);
    aT.SetSolver(solver);
    DrawdownData d;
    d.nivelInicial = points[4][1]; // reservatorio cheio: na crista
    d.nivelFinal = points[3][1];   // vazio: no pe
    d.nivelMacico = points[4][1];  // NA do macico (lateral esquerda), fixo
    d.tRebaixamento = horasRebaixamento;
    DrawdownTransient rebaixamento(aT, d, left, {right, topright, face, topleft});
    rebaixamento.InitialState();
    const PointLocator loc(mT.GeoMesh());

    // fator de seguranca: GIM com as pressoes do transiente congeladas
    UnsatModel<TUnsatHP> mF;
    mF.SetGeoMesh(gF);
    InsertBC(mF, points);
    // agua do reservatorio sobre o talude submerso: pressao gamma_w (nivel - y) normal a face e ao terreno do pe, que
    // cresce com o fator do GIM como o peso (fator e nivel mudam ao longo da conta: a funcao le os dois)
    REAL fator = 1., nivel = d.nivelInicial;
    auto pressaoReservatorio = [&](REAL nx, REAL ny) { // normal para fora do solo
        return [&fator, &nivel, nx, ny, gw = d.gammaw](const TPZVec<REAL> &x, TPZVec<STATE> &v2, TPZFMatrix<STATE> &) {
            const REAL pr = fator * gw * std::max(REAL(0), nivel - x[1]);
            v2[0] = -pr * nx, v2[1] = -pr * ny;
        };
    };
    const REAL dx = points[3][0] - points[4][0], dy = points[4][1] - points[3][1], L = std::hypot(dx, dy);
    mF.AddNaturalBC(face, TUnsat::ENeumannU, {0., 0.}, pressaoReservatorio(dy / L, dx / L));
    mF.AddNaturalBC(topright, TUnsat::ENeumannU, {0., 0.}, pressaoReservatorio(0., 1.));
    mF.SetMaterial(UnsatMaterial(matSolo, p));
    mF.CreateCompMesh({{2, 2}, {1, 1}, {1, 1}});
    UnsatAnalysis<TUnsatHP> aF(mF);
    aF.SetVerbose(false);
    TUnsatHP *mat = mF.Material();
    aF.SetLoadFactorFunction([mat, &fator](REAL l) { mat->SetGravityFactor(l), fator = l; }, false);

    AnimationFiles<TUnsatHP> anim("rebaixamento", mF);
    anim.Summary(horasRebaixamento, d.nivelInicial, d.nivelFinal);
    const int n = int(std::lround((horasRebaixamento + horasDepois) / dtGravar));
    for (int k = 0; k <= n; k++) {
        const REAL t = k * dtGravar;
        if (!rebaixamento.AdvanceTo(t)) {
            std::cout << "o transiente parou em t = " << rebaixamento.Time() << " h" << std::endl;
            break;
        }
        const std::vector<REAL> pw = rebaixamento.NodalWaterPressure();
        anim.Snapshot(t, pw, rebaixamento.Level(), mT.Material()->Hydraulics());
        // o fator de hora em hora e tambem no fim da descida
        if (std::abs(t - std::round(t)) > 1.e-9 && std::abs(t - horasRebaixamento) > 1.e-9) continue;
        mat->ResetMemory(); // estado sem tensoes
        nivel = rebaixamento.Level();
        aF.SetFrozenPressures([&](int espaco, const TPZVec<REAL> &x) {
            return espaco == 1 ? loc.Interpolate(x[0], x[1], pw) : 0.;
        });
        const REAL gamma = aF.GIMArcLength();
        anim.FactorOfSafety(t, gamma);
        std::cout << "t = " << t << " h, nivel = " << rebaixamento.Level() << " m: Gamma = " << gamma << std::endl;
    }
}

int main()
{
    // cantos do talude no sentido anti-horario: base esquerda, base direita, terreno do pe a direita, pe, crista, crista
    // a esquerda
    const std::vector<std::vector<REAL>> talude10 = {{0, 0}, {70, 0}, {70, 30}, {40, 30}, {30, 40}, {0, 40}}; // H = 10 m
    const std::vector<std::vector<REAL>> ceron = {{0, 0}, {25, 0}, {25, 5}, {15, 5}, {10, 10}, {0, 10}}; // H = 5 m, Fig. 24]
    const std::vector<std::vector<REAL>> ceron2 = {{0, 0}, {70, 0}, {70, 30}, {40, 30}, {35, 35}, {0, 35}}; // H = 5 m, Fig. 24
    // ---- niveis de refinamento: h nos triangulos marcados e, com pmax > 2, tambem p (hp)
     SolveDrawDown(ceron2);
    //TransientDrawDown(ceron2);
    return 0;
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
