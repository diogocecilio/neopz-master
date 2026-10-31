// LearningSlopeWater: talude saturado (etapa com agua do estudo do Mohr-Coulomb PV2 com o material nao saturado
// u - p_w - p_a, ar passivo). Transiente com elevacao do nivel d'agua a montante e fator de seguranca por reducao de
// resistencia (SRM) a cada 10 passos de tempo, com o balanco de massa da agua acoplado e desacoplado da deformacao.
// Saidas: fs_tempo.csv (FS ao longo do tempo) e pressao_tempo.csv (p_w em pontos de controle ao longo do tempo).

// geometria
#include "TPZGmshReader.h"
#include "pzgeoel.h"
#include "pzgmesh.h"

// malhas computacionais
#include "TPZMultiphysicsCompMesh.h"
#include "TPZNullMaterial.h"
#include "pzcmesh.h"

// material nao saturado u - p_w - p_a com Mohr-Coulomb PV2 do artigo
#include "Plasticity/TPZElasticResponse.h"
#include "Plasticity/TPZMatUnsaturatedUPwPa.h"
#include "Plasticity/TPZPlasticStepVoigt.h"
#include "Plasticity/TPZYCMohrCoulombPV2.h"

// analise incremental (Newton, filtro de equacoes, memoria do estado convergido)
#include "TPZUnsatAnalysis.h"

#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <vector>

using TMC = TPZPlasticStepVoigt<TPZYCMohrCoulombPV2, TPZElasticResponse>;
using TUnsat = TPZMatUnsaturatedUPwPa<TMC, TPZUnsaturatedMem>;

/// Parameters of a soil (units kN, m, s, kPa); defaults: chi = S_w (notes), van Genuchten - Mualem
struct UnsatParams {
    REAL E = 2.e4, nu = 0.3, phi = 30. * M_PI / 180., c = 10.;
    REAL gamma = 20.;            ///< constant unit weight (kN/m3)
    bool mixture = false;        ///< unit weight of the mixture (1-n) gamma_s + n S_w gamma_w
    REAL gammas = 26.;           ///< unit weight of the grains (mixture)
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

// Talude saturado (geometria e resistencia de de Souza Neto et al., Fig. 8.30; fluxo associado, psi = phi):
// limite saturado do van Genuchten (a_VG = 1e-8: S_w = 1 e k_rw = 1 para qualquer pressao), ar passivo
UnsatParams SaturatedSlopeParams(bool coupled)
{
    UnsatParams p;
    p.E = 2.e4, p.nu = 0.49, p.c = 50., p.phi = 20. * M_PI / 180., p.gamma = 20.;
    p.hydro.fAVG = 1.e-8; // saturado
    p.ks = 1.e-6;         // condutividade hidraulica (m/s)
    p.coupled = coupled;  // deformacao no balanco de massa da agua
    return p;
}

// Cenario hidraulico: nivel d'agua inicial no pe do talude (y = 30 m); em t = 0 o nivel no contorno esquerdo
// sobe para hL; no contorno direito fica em 30 m e o terreno a jusante do pe (y = 30 m, x >= 45 m) fica drenado
// (p_w = 0). Base, crista e face impermeaveis.
struct WaterScenario {
    REAL y0 = 30.;        // nivel inicial (m)
    REAL hL = 38.;        // nivel a montante depois de t = 0 (m)
    REAL xToe = 45.;      // abscissa do pe do talude (m)
    REAL tEnd = 2.e6;     // tempo final (s)
    int nSteps = 50;      // passos de tempo (progressao geometrica)
    REAL dt0 = 60.;       // primeiro passo de tempo (s)
    int srmEvery = 10;    // SRM a cada srmEvery passos
};

// ids de material: grupos fisicos das malhas do Gmsh (meshes/talude.geo e meshes/sapata.geo)
enum EIds {
    matSolo = 1,       // solo (dominio)
    bcBase = -1,       // base
    bcDireita = -2,    // lateral direita
    bcSuperficie = -3, // superficie livre do terreno
    bcEsquerda = -4,   // lateral esquerda (talude) ou eixo de simetria (sapata)
    bcSapata = -5      // contato da sapata
};

TPZGeoMesh *ReadGmshMesh(const std::string &file, const std::map<std::string, int> &domains,
                         const std::map<std::string, int> &boundaries);
// Talude: 75 m x 40/30 m, h = 10 m, 45 graus
TPZGeoMesh *ReadSlopeMesh();

TUnsat *UnsatMaterial(int id, const UnsatParams &p);
TPZCompMesh *CreateAtomicMesh(TPZGeoMesh *gmesh, int nstate, int order, const std::set<int> &bcids);
TPZMultiphysicsCompMesh *CreateSlopeCompMesh(TPZGeoMesh *gmesh, int porder, const UnsatParams &p);
// Apaga a malha multifisica e as suas malhas atomicas (depois da analise que as usa)
void DeleteCompMesh(TPZMultiphysicsCompMesh *mphys);

// Fator de seguranca por reducao de resistencia com as pressoes de x congeladas
REAL FactorSRM(TPZUnsatAnalysis &an, const TPZFMatrix<STATE> &x, REAL F0 = 0.6, REAL tol = 0.01,
               const std::vector<REAL> &Fs = {0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5, 3.0, 3.5, 4.0});

// FS do estado convergido atual sem altera-lo: o SRM parte do estado sem tensoes; depois a memoria e a solucao
// convergidas do transiente sao restauradas
REAL FactorSRMKeepingState(TPZUnsatAnalysis &an);

// Transiente de um caso (acoplado ou desacoplado): grava as linhas dos CSVs com o nome do caso
void RunCase(bool coupled, const WaterScenario &sc, std::ofstream &fsCsv, std::ofstream &pCsv);

int main() {
    WaterScenario sc;
    std::ofstream fsCsv("fs_tempo.csv"), pCsv("pressao_tempo.csv");
    fsCsv << "caso,passo,t_s,t_h,FS\n";
    pCsv << "caso,passo,t_s,t_h,ponto,x,y,pw_kPa\n";
    RunCase(true, sc, fsCsv, pCsv);
    RunCase(false, sc, fsCsv, pCsv);
    return 0;
}

void RunCase(bool coupled, const WaterScenario &sc, std::ofstream &fsCsv, std::ofstream &pCsv)
{
    const std::string caso = coupled ? "acoplado" : "desacoplado";
    const UnsatParams p = SaturatedSlopeParams(coupled);
    TPZGeoMesh *gmesh = ReadSlopeMesh();
    TPZMultiphysicsCompMesh *cmesh = CreateSlopeCompMesh(gmesh, 2, p);
    auto *an = new TPZUnsatAnalysis(cmesh, {{bcBase, {true, true}}, {bcDireita, {true, false}},
                                            {bcEsquerda, {true, false}}});

    // pontos de controle da poropressao: vertice mais proximo de cada ponto
    const std::vector<std::pair<std::string, std::array<REAL, 2>>> targets = {
        {"sob_crista", {20., 30.}}, {"sob_face", {38., 30.}}, {"sob_pe", {45., 25.}},
        {"fundacao", {40., 20.}},   {"jusante", {60., 25.}}};
    std::vector<std::pair<std::string, int64_t>> probes; // nome, no
    for (auto &[nome, xy] : targets) {
        int64_t best = -1;
        REAL dmin = std::numeric_limits<REAL>::max();
        for (auto &[node, eq] : an->WaterEq()) {
            const REAL d = std::hypot(an->Coord(node, 0) - xy[0], an->Coord(node, 1) - xy[1]);
            if (d < dmin) dmin = d, best = node;
        }
        probes.push_back({nome, best});
    }

    // estado inicial: hidrostatico com o nivel y0 e ar a pressao atmosferica (p_a = 0, passivo)
    TPZFMatrix<STATE> &X = an->X();
    X.Zero();
    for (auto &[node, eq] : an->WaterEq()) X(eq, 0) = p.gw * (sc.y0 - an->Coord(node, 1));
    an->Restore();
    // estado geostatico com as pressoes iniciais congeladas (peso e pressao em 4 incrementos)
    {
        const auto fixed = an->FrozenFixed(X);
        an->SetGravityScalesSuction(true);
        for (int k = 1; k <= 4; k++) {
            an->SetGravityFactor(k / 4.);
            if (!an->Solve(0., fixed, {}, true).conv) {
                std::cout << caso << ": estado geostatico nao convergiu" << std::endl;
                DebugStop();
            }
            an->Accept();
        }
    }

    // condicoes de contorno da agua no transiente; ar passivo em todo o dominio
    std::map<int64_t, REAL> fixed = an->DisplacementConstraints();
    for (auto &[node, eq] : an->AirEq()) fixed[eq] = 0.;
    for (auto &[node, len] : an->BoundaryNodes(bcEsquerda))
        fixed[an->WaterEq().at(node)] = p.gw * (sc.hL - an->Coord(node, 1));
    for (auto &[node, len] : an->BoundaryNodes(bcDireita))
        fixed[an->WaterEq().at(node)] = p.gw * (sc.y0 - an->Coord(node, 1));
    for (auto &[node, len] : an->BoundaryNodes(bcSuperficie))
        if (an->Coord(node, 0) >= sc.xToe - 1.e-6 && std::fabs(an->Coord(node, 1) - sc.y0) < 1.e-6)
            fixed[an->WaterEq().at(node)] = 0.; // terreno a jusante do pe: drenado

    // passos de tempo em progressao geometrica: dt0 r^k, soma = tEnd
    REAL r = 1.1;
    for (int it = 0; it < 100; it++) // dt0 (r^n - 1)/(r - 1) = tEnd
        r = std::pow(sc.tEnd * (r - 1.) / sc.dt0 + 1., 1. / sc.nSteps);

    auto record = [&](int step, REAL t) {
        for (auto &[nome, node] : probes)
            pCsv << caso << ',' << step << ',' << t << ',' << t / 3600. << ',' << nome << ',' << an->Coord(node, 0)
                 << ',' << an->Coord(node, 1) << ',' << an->X()(an->WaterEq().at(node), 0) << '\n';
        pCsv.flush();
        if (step % sc.srmEvery == 0) {
            const REAL FS = FactorSRMKeepingState(*an);
            fsCsv << caso << ',' << step << ',' << t << ',' << t / 3600. << ',' << FS << '\n';
            std::cout << caso << ": passo " << step << ", t = " << t / 3600. << " h, FS = " << FS << std::endl;
        }
        fsCsv.flush();
    };

    REAL t = 0., dt = sc.dt0;
    record(0, t);
    for (int step = 1; step <= sc.nSteps; step++) {
        const auto res = an->Solve(dt, fixed, {}, false);
        if (!res.conv) {
            std::cout << caso << ": passo " << step << " nao convergiu" << std::endl;
            break;
        }
        an->Accept();
        t += dt;
        record(step, t);
        dt *= r;
    }

    delete an;
    DeleteCompMesh(cmesh);
    delete gmesh;
}

// FS do estado convergido atual sem altera-lo: o SRM parte do estado sem tensoes com as pressoes atuais congeladas;
// depois a memoria dos pontos de integracao e a solucao convergidas do transiente sao restauradas
REAL FactorSRMKeepingState(TPZUnsatAnalysis &an)
{
    const auto memory = an.SaveMemory();
    const TPZFMatrix<STATE> x = an.X();
    const REAL FS = FactorSRM(an, x);
    an.RestoreMemory(memory);
    an.X() = x;
    an.Restore();
    an.SetGravityFactor(1.);
    an.SetGravityScalesSuction(true);
    an.SetStrengthReductionFactor(1.);
    return FS;
}

// Fator de seguranca por reducao de resistencia (SRM) com as pressoes de x congeladas (so o equilibrio e
// resolvido): serve as tres etapas, com x = pressoes nulas (seco), p_w de um instante (com agua) ou p_w e p_a de um
// instante (com agua e ar). Parte do estado sem tensoes (memoria inicial, u = 0); estado geostatico com a
// resistencia aumentada (F0 < 1, peso e sucao em 4 incrementos); depois F cresce pela lista Fs, sempre a partir do
// ultimo estado convergido; no primeiro F sem convergencia, bissecao entre o ultimo F convergido e ele ate a
// diferenca ser menor que tol. Retorna o ultimo F convergido (-1 se o estado geostatico nao convergir).
REAL FactorSRM(TPZUnsatAnalysis &an, const TPZFMatrix<STATE> &x, REAL F0, REAL tol, const std::vector<REAL> &Fs)
{
    // estado sem tensoes com as pressoes de x
    an.ResetMemory();
    TPZFMatrix<STATE> &X = an.X();
    X = x;
    for (int64_t i = 0; i < an.NEq(); i++)
        if (!an.IsPressure(i)) X(i, 0) = 0.;
    an.Restore();
    const auto fixed = an.FrozenFixed(X); // apoios = 0 e todas as pressoes nos valores de x

    // estado geostatico: peso (e sucao) em incrementos, resistencia aumentada
    an.SetStrengthReductionFactor(F0);
    an.SetGravityScalesSuction(true);
    const int n = 4;
    for (int k = 1; k <= n; k++) {
        an.SetGravityFactor(REAL(k) / n);
        if (!an.Solve(0., fixed, {}, true).conv) {
            std::cout << "SRM: estado geostatico nao convergiu" << std::endl;
            return -1.;
        }
        an.Accept();
    }

    // reducao da resistencia
    REAL Fok = F0, Ffail = -1.;
    size_t next = 0;
    while (true) {
        REAL F;
        if (Ffail > 0.) {
            if (Ffail - Fok < tol) break;
            F = 0.5 * (Fok + Ffail);
        } else {
            while (next < Fs.size() && Fs[next] <= Fok) next++;
            if (next == Fs.size()) break;
            F = Fs[next];
        }
        an.SetStrengthReductionFactor(F);
        const auto r = an.Solve(0., fixed, {}, true);
        std::cout << "SRM: F = " << F << (r.conv ? " convergiu" : " nao convergiu") << " (" << r.iter << " iteracoes)"
                  << std::endl;
        if (r.conv) {
            an.Accept();
            Fok = F;
        } else {
            an.Restore();
            Ffail = F;
        }
    }
    an.SetStrengthReductionFactor(1.);
    return Fok;
}

// Malha atomica H1 com material nulo: nstate = 2 (deslocamento) ou 1 (poropressao)
TPZCompMesh *CreateAtomicMesh(TPZGeoMesh *gmesh, int nstate, int order, const std::set<int> &bcids)
{
    auto *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(2);
    cmesh->SetDefaultOrder(order);
    cmesh->SetAllCreateFunctionsContinuous();
    auto *mat = new TPZNullMaterial<STATE>(matSolo, 2, nstate);
    cmesh->InsertMaterialObject(mat);
    TPZFMatrix<STATE> val1(nstate, nstate, 0.);
    TPZManVector<STATE, 2> val2(nstate, 0.);
    for (int id : bcids) cmesh->InsertMaterialObject(mat->CreateBC(mat, id, 0, val1, val2));
    cmesh->AutoBuild();
    return cmesh;
}

// Talude com o TUnsat: base engastada, laterais com rolete (u_x = 0), superficie livre.
// Espacos: u de ordem porder (2 = triangulo de 6 nos), p_w e p_a de ordem porder - 1 (1: P1).
// Apoios impostos pela analise (TPZUnsatAnalysis,
// fixU = {{bcBase, {true, true}}, {bcDireita, {true, false}}, {bcEsquerda, {true, false}}}):
// as condicoes EDirichletU do TUnsat nao contribuem, so marcam os contornos.
TPZMultiphysicsCompMesh *CreateSlopeCompMesh(TPZGeoMesh *gmesh, int porder, const UnsatParams &p)
{
    const std::set<int> bcids = {bcBase, bcDireita, bcEsquerda};
    TPZManVector<TPZCompMesh *, 3> meshvec = {CreateAtomicMesh(gmesh, 2, porder, bcids),
                                              CreateAtomicMesh(gmesh, 1, porder - 1, bcids),
                                              CreateAtomicMesh(gmesh, 1, porder - 1, bcids)};
    TUnsat *mat = UnsatMaterial(matSolo, p);
    auto *mphys = new TPZMultiphysicsCompMesh(gmesh);
    mphys->SetDimModel(2);
    mphys->InsertMaterialObject(mat);

    TPZFMatrix<STATE> val1(2, 2, 0.);
    TPZManVector<STATE, 2> val2(2, 0.);
    for (int id : bcids) mphys->InsertMaterialObject(mat->CreateBC(mat, id, TUnsat::EDirichletU, val1, val2));

    TPZManVector<int, 3> active(3, 1);
    mphys->BuildMultiphysicsSpaceWithMemory(active, meshvec, {matSolo}, bcids);
    return mphys;
}

// Apaga a malha multifisica e as suas malhas atomicas (depois da analise que as usa)
void DeleteCompMesh(TPZMultiphysicsCompMesh *mphys)
{
    TPZVec<TPZCompMesh *> meshvec = mphys->MeshVector();
    delete mphys;
    for (TPZCompMesh *m : meshvec) delete m;
}

// Material TUnsat com Mohr-Coulomb associado (psi = phi) a partir dos parametros (copia de UnsatModel.h)
TUnsat *UnsatMaterial(int id, const UnsatParams &p)
{
    TPZElasticResponse ER;
    ER.SetEngineeringData(p.E, p.nu);
    auto *mat = new TUnsat(id);
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
    mat->SetIntegrationOrder(2); // regra de 3 pontos
    return mat;
}

// Le uma malha do Gmsh (.msh 4.1) e associa os nomes dos grupos fisicos aos ids de material.
//   domains    : grupos de dimensao 2 (dominio)
//   boundaries : grupos de dimensao 1 (contornos)
// Todo elemento deve pertencer a um grupo listado nos mapas.
TPZGeoMesh *ReadGmshMesh(const std::string &file, const std::map<std::string, int> &domains,
                         const std::map<std::string, int> &boundaries)
{
    TPZManVector<std::map<std::string, int>, 4> names(4);
    for (auto &[name, id] : domains) names[2][name] = id;
    for (auto &[name, id] : boundaries) names[1][name] = id;
    TPZGmshReader reader;
    reader.SetDimNamePhysical(names);
    TPZGeoMesh *gmesh = reader.GeometricGmshMesh(file);
    for (TPZGeoEl *gel : gmesh->ElementVec()) {
        if (!gel) continue;
        const int dim = gel->Dimension();
        bool known = false;
        for (auto &[name, id] : names[dim]) known = known || (id == gel->MaterialId());
        if (!known) {
            std::cout << "ReadGmshMesh: elemento de dimensao " << dim << " com id " << gel->MaterialId()
                      << " fora dos grupos fisicos informados" << std::endl;
            DebugStop();
        }
    }
    gmesh->BuildConnectivity();
    return gmesh;
}

// Talude: 75 m x 40/30 m, h = 10 m, 45 graus
TPZGeoMesh *ReadSlopeMesh()
{
    return ReadGmshMesh(std::string(MESH_DIR) + "/talude.msh", {{"Solo", matSolo}},
                        {{"Base", bcBase}, {"Direita", bcDireita}, {"Superficie", bcSuperficie}, {"Esquerda", bcEsquerda}});
}

