// LearningMohrCoulomb: verificacao do Mohr-Coulomb (PV2) em deformacao plana, passo a passo
// (projeto de estudo): capacidade de carga de uma sapata corrida e estabilidade de um talude,
// em solo seco com o material nao saturado u - p_w - p_a (TUnsat).

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

// analise incremental do modelo u - p_w - p_a (Newton, filtro de equacoes, memoria do estado convergido)
#include "TPZUnsatAnalysis.h"

// pos-processamento
#include "TPZVTKGeoMesh.h"

#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

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
// Sapata corrida: 5 m x 5 m, meia sapata de 0,5 m com simetria em x = 0
TPZGeoMesh *ReadFootingMesh();

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

// Solo seco com o TUnsat: u P2, p_w P1 e p_a P1. As pressoes ficam congeladas em zero na analise
// (Solve com frozen = true); os apoios tambem sao impostos pela analise (UnsatAnalysis, fixU).
// Os contornos entram na malha so como elementos sem contribuicao.

// Talude (de Souza Neto et al., Fig. 8.30): fluxo associado (psi = phi), gravidade em -y
UnsatParams SlopeUnsatParams()
{
    UnsatParams p;
    p.E = 2.e4, p.nu = 0.49, p.c = 50., p.phi = 20. * M_PI / 180., p.gamma = 20.;
    return p;
}

// Sapata corrida (Tabela 1 do artigo): sem peso proprio, fluxo associado
UnsatParams FootingUnsatParams()
{
    UnsatParams p;
    p.E = 1.e7, p.nu = 0.48, p.c = 490., p.phi = 20. * M_PI / 180., p.gamma = 0.;
    return p;
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

// Talude em solo seco com o TUnsat: base engastada, laterais com rolete (u_x = 0), superficie livre.
// Espacos: u de ordem porder (2 = triangulo de 6 nos), p_w e p_a de ordem porder - 1 (1: P1).
// Apoios impostos pela analise (UnsatAnalysis,
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

// Sapata corrida (meio modelo) em solo seco com o TUnsat: base e laterais com rolete, simetria em x = 0 e
// recalque u_y prescrito sob a sapata. Apoios e recalque impostos pela analise (UnsatAnalysis,
// fixU = {{bcBase, {false, true}}, {bcDireita, {true, false}}, {bcEsquerda, {true, false}}, {bcSapata, {false, true}}};
// o recalque entra no mapa fixed do Solve, nos vertices da sapata).
TPZMultiphysicsCompMesh *CreateFootingCompMesh(TPZGeoMesh *gmesh, int porder, const UnsatParams &p)
{
    const std::set<int> bcids = {bcBase, bcDireita, bcEsquerda, bcSapata};
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

// Fator de seguranca por reducao de resistencia (SRM) em solo seco (pressoes congeladas em zero): estado
// geostatico com a resistencia aumentada (F0 < 1, peso em 4 incrementos); depois F cresce pela lista Fs, sempre a
// partir do ultimo estado convergido; no primeiro F sem convergencia, bissecao entre o ultimo F convergido e ele ate
// a diferenca ser menor que tol. Retorna o ultimo F convergido (-1 se o estado geostatico nao convergir).
REAL FactorSRM(TPZUnsatAnalysis &an, REAL F0 = 0.6, REAL tol = 0.01,
               const std::vector<REAL> &Fs = {0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5, 3.0, 3.5, 4.0})
{
    const auto fixed = an.FrozenFixed(an.X()); // apoios = 0 e todas as pressoes = 0
    an.SetStrengthReductionFactor(F0);
    const int n = 4;
    for (int k = 1; k <= n; k++) {
        an.SetGravityFactor(REAL(k) / n);
        if (!an.Solve(0., fixed, {}, true).conv) {
            std::cout << "SRM: estado geostatico nao convergiu" << std::endl;
            return -1.;
        }
        an.Accept();
    }
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

int main() {
    // talude em solo seco: base engastada, laterais com rolete
    TPZGeoMesh *gmesh = ReadSlopeMesh();
    std::ofstream out("talude.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
    TPZMultiphysicsCompMesh *cmesh = CreateSlopeCompMesh(gmesh, 2, SlopeUnsatParams());
    TPZUnsatAnalysis an(cmesh, {{bcBase, {true, true}}, {bcDireita, {true, false}}, {bcEsquerda, {true, false}}});
    const REAL FS = FactorSRM(an);
    std::cout << "FS (SRM) = " << FS << std::endl;
    return 0;
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

// Sapata corrida: 5 m x 5 m, meia sapata de 0,5 m com simetria em x = 0
TPZGeoMesh *ReadFootingMesh()
{
    return ReadGmshMesh(std::string(MESH_DIR) + "/sapata.msh", {{"Solo", matSolo}},
                        {{"Base", bcBase}, {"Direita", bcDireita}, {"Superficie", bcSuperficie},
                         {"Simetria", bcEsquerda}, {"Sapata", bcSapata}});
}
