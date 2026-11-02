// LearningMohrCoulomb: verificacao do Mohr-Coulomb (PV2) em deformacao plana, passo a passo
// (projeto de estudo): capacidade de carga de uma sapata corrida e estabilidade de um talude com o material
// nao saturado u - p_w - p_a (TUnsat), em tres etapas: solo seco, com agua e com agua e ar.

// geometria
#include "TPZGmshReader.h"
#include "pzgeoel.h"
#include "pzgmesh.h"
#include "pzfstrmatrix.h"
#include "pzpostprocanalysis.h"
// malhas computacionais
#include "TPZBndCondT.h"
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

TUnsat *UnsatMaterial(int id, const UnsatParams &p);

TPZCompMesh *CreateAtomicMesh(TPZGeoMesh *gmesh, int nstate, int order, const std::set<int> &bcids);
TPZMultiphysicsCompMesh *CreateSlopeCompMesh(TPZGeoMesh *gmesh, int porder, const UnsatParams &p);
TPZMultiphysicsCompMesh *CreateFootingCompMesh(TPZGeoMesh *gmesh, int porder, const UnsatParams &p, REAL q);
// Apaga a malha multifisica e as suas malhas atomicas (depois da analise que as usa)
void DeleteCompMesh(TPZMultiphysicsCompMesh *mphys);

// Fator de seguranca por reducao de resistencia com as pressoes de x congeladas (seco: x com pressoes nulas)
REAL FactorSRM(TPZUnsatAnalysis &an, const TPZFMatrix<STATE> &x, REAL F0 = 0.6, REAL tol = 0.01,
               const std::vector<REAL> &Fs = {0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5, 3.0, 3.5, 4.0});
// Fator de seguranca por aumento da gravidade com as pressoes de x congeladas (seco: x com pressoes nulas)
REAL FactorGIM(TPZUnsatAnalysis &an, const TPZFMatrix<STATE> &x, REAL tol = 0.01,
               const std::vector<REAL> &lams = {0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.0});
// Fator de seguranca por aumento da carga de superficie do contorno bcid (sapata) com as pressoes de x congeladas
REAL FactorLoad(TPZUnsatAnalysis &an, const TPZFMatrix<STATE> &x, int bcid, REAL tol = 0.01,
                const std::vector<REAL> &lams = {0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.0});

// Pos-processamento VTK do estado carregado na malha: campos da solucao (deslocamento e pressoes) em
// <nome>.scal_vec.<step>.vtk e variaveis dos pontos de integracao (tensoes efetivas e plasticidade) em
// <nome>_pontos.scal_vec.<step>.vtk
void PostProcessVTK(TPZUnsatAnalysis &an, const std::string &nome, int step = 0);

int main() {
    const int porder = 2;

    // ---- talude ----
    TPZGeoMesh *gmeshslope = ReadSlopeMesh();
    std::ofstream out2("talude.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmeshslope, out2, true);

    TPZMultiphysicsCompMesh *cslope = CreateSlopeCompMesh(gmeshslope, porder, SlopeUnsatParams());
    auto *anslope = new TPZUnsatAnalysis(cslope, {{bcBase, {true, true}}, {bcDireita, {true, false}},
                                                  {bcEsquerda, {true, false}}});

    // etapa 1: solo seco (todas as pressoes nulas)
    TPZFMatrix<STATE> xslope(cslope->NEquations(), 1, 0.);
    const REAL FSslope = FactorSRM(*anslope, xslope);
    std::cout << "Talude: FS (SRM, seco) = " << FSslope << std::endl;
    PostProcessVTK(*anslope, "talude_SRM");

    const REAL FGslope = FactorGIM(*anslope, xslope);
    std::cout << "Talude: FS (GIM, seco) = " << FGslope << std::endl;
    PostProcessVTK(*anslope, "talude_GIM");

    delete anslope; // a analise antes das malhas que ela usa
    DeleteCompMesh(cslope);
    delete gmeshslope;

    // ---- sapata: pressao de trabalho q (sem peso proprio, a carga e a da sapata) ----
    TPZGeoMesh *gmeshfooting = ReadFootingMesh();
    std::ofstream out("sapata.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmeshfooting, out, true);

    const REAL q = 3000.; // kPa
    const REAL B = 0.5;   // meia largura da sapata (m)
    TPZMultiphysicsCompMesh *cfoot = CreateFootingCompMesh(gmeshfooting, porder, FootingUnsatParams(), q);
    // sapata com carga prescrita: o recalque fica livre (bcSapata fora dos apoios)
    auto *anfoot = new TPZUnsatAnalysis(cfoot, {{bcBase, {false, true}}, {bcDireita, {true, false}},
                                                {bcEsquerda, {true, false}}});
    anfoot->SetForceReference(q * B); // sem peso proprio: referencia do residuo = forca da sapata (kN/m)

    // etapa 1: solo seco (todas as pressoes nulas)
    TPZFMatrix<STATE> xfoot(cfoot->NEquations(), 1, 0.);
    const REAL FSfoot = FactorSRM(*anfoot, xfoot);
    std::cout << "Sapata: FS (SRM, seco) = " << FSfoot << std::endl;
    PostProcessVTK(*anfoot, "sapata_SRM");

    // no lugar do GIM (o solo nao tem peso): aumento da carga da sapata
    const REAL FLfoot = FactorLoad(*anfoot, xfoot, bcSapata);
    std::cout << "Sapata: FS (carga, seco) = " << FLfoot << "  -> q_u = " << FLfoot * q << " kPa" << std::endl;
    PostProcessVTK(*anfoot, "sapata_carga");

    delete anfoot;
    DeleteCompMesh(cfoot);
    delete gmeshfooting;
    return 0;
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

// Fator de seguranca por aumento da gravidade (GIM) com as pressoes de x congeladas: parte do estado sem tensoes,
// com a resistencia real (F = 1), e o fator lambda multiplica o peso (a sucao fica com o valor de x); lambda cresce
// pela lista lams, sempre a partir do ultimo estado convergido; no primeiro lambda sem convergencia, bissecao entre o
// ultimo lambda convergido e ele ate a diferenca ser menor que tol. Retorna o ultimo lambda convergido.
REAL FactorGIM(TPZUnsatAnalysis &an, const TPZFMatrix<STATE> &x, REAL tol, const std::vector<REAL> &lams)
{
    // estado sem tensoes com as pressoes de x
    an.ResetMemory();
    TPZFMatrix<STATE> &X = an.X();
    X = x;
    for (int64_t i = 0; i < an.NEq(); i++)
        if (!an.IsPressure(i)) X(i, 0) = 0.;
    an.Restore();
    const auto fixed = an.FrozenFixed(X); // apoios = 0 e todas as pressoes nos valores de x

    an.SetStrengthReductionFactor(1.);
    an.SetGravityScalesSuction(false); // lambda so no peso
    REAL lok = 0., lfail = -1.;
    size_t next = 0;
    while (true) {
        REAL lam;
        if (lfail > 0.) {
            if (lfail - lok < tol) break;
            lam = 0.5 * (lok + lfail);
        } else {
            while (next < lams.size() && lams[next] <= lok) next++;
            if (next == lams.size()) break;
            lam = lams[next];
        }
        an.SetGravityFactor(lam);
        const auto r = an.Solve(0., fixed, {}, true);
        std::cout << "GIM: lambda = " << lam << (r.conv ? " convergiu" : " nao convergiu") << " (" << r.iter
                  << " iteracoes)" << std::endl;
        if (r.conv) {
            an.Accept();
            lok = lam;
        } else {
            an.Restore();
            lfail = lam;
        }
    }
    an.SetGravityScalesSuction(true);
    an.SetGravityFactor(1.);
    return lok;
}

// Fator de seguranca por aumento da carga de superficie (sapata): como o GIM, mas o fator lambda multiplica a
// tracao do contorno bcid (ENeumannU, Val2) em vez do peso. Parte do estado sem tensoes, com a resistencia real e as
// pressoes de x congeladas; lambda cresce pela lista lams, sempre a partir do ultimo estado convergido; no primeiro
// lambda sem convergencia, bissecao ate a diferenca ser menor que tol. Retorna o ultimo lambda convergido; ao final
// a tracao volta ao valor original.
REAL FactorLoad(TPZUnsatAnalysis &an, const TPZFMatrix<STATE> &x, int bcid, REAL tol, const std::vector<REAL> &lams)
{
    auto *bc = dynamic_cast<TPZBndCondT<STATE> *>(an.Mesh()->FindMaterial(bcid));
    if (!bc) DebugStop();
    const TPZManVector<STATE, 2> t0(bc->Val2());

    // estado sem tensoes com as pressoes de x
    an.ResetMemory();
    TPZFMatrix<STATE> &X = an.X();
    X = x;
    for (int64_t i = 0; i < an.NEq(); i++)
        if (!an.IsPressure(i)) X(i, 0) = 0.;
    an.Restore();
    const auto fixed = an.FrozenFixed(X); // apoios = 0 e todas as pressoes nos valores de x

    an.SetStrengthReductionFactor(1.);
    REAL lok = 0., lfail = -1.;
    size_t next = 0;
    while (true) {
        REAL lam;
        if (lfail > 0.) {
            if (lfail - lok < tol) break;
            lam = 0.5 * (lok + lfail);
        } else {
            while (next < lams.size() && lams[next] <= lok) next++;
            if (next == lams.size()) break;
            lam = lams[next];
        }
        TPZManVector<STATE, 2> t(t0);
        for (auto &v : t) v *= lam;
        bc->SetVal2(t);
        const auto r = an.Solve(0., fixed, {}, true);
        std::cout << "Carga: lambda = " << lam << (r.conv ? " convergiu" : " nao convergiu") << " (" << r.iter
                  << " iteracoes)" << std::endl;
        if (r.conv) {
            an.Accept();
            lok = lam;
        } else {
            an.Restore();
            lfail = lam;
        }
    }
    bc->SetVal2(t0);
    return lok;
}

// Pos-processamento VTK do estado carregado na malha (ultimo estado convergido): todas as variaveis projetadas
// por elemento com o TPZPostProcAnalysis (campo linear descontinuo), em <nome>_pontos.scal_vec.<step>.vtk:
// solucao (deslocamento, poropressoes, succao e saturacao) e memoria dos pontos de integracao (tensoes efetivas,
// tracao positiva, e plasticidade)
void PostProcessVTK(TPZUnsatAnalysis &an, const std::string &nome, int step)
{
    TPZPostProcAnalysis post;
    post.SetCompMesh(an.Mesh());
    TPZManVector<int, 1> matids = {matSolo};
    TPZManVector<std::string, 14> scal = {"DisplacementX", "DisplacementY", "WaterPressure", "AirPressure",
        "Suction",       "Saturation",    "StressXX",      "StressYY",
        "StressZZ",      "StressXY",      "StrainPlasticI1", "StrainPlasticJ2",
        "FailureType",   "Hardening"};
        TPZManVector<std::string, 1> vec = {"Displacement"};
        // o pos-processamento projeta todas as variaveis: escalares e vetores
        TPZManVector<std::string, 15> vars(scal.size() + vec.size());
        for (int i = 0; i < scal.size(); i++) vars[i] = scal[i];
        for (int i = 0; i < vec.size(); i++) vars[scal.size() + i] = vec[i];
        post.SetPostProcessVariables(matids, vars);
    TPZFStructMatrix<STATE> sm(post.Mesh());
    sm.SetNumThreads(0);
    post.SetStructuralMatrix(sm);
    post.DefineGraphMesh(2, scal, vec, nome + "_pontos.vtk"); // Displacement como vetor
    post.TransferSolution();
    post.SetStep(step);
    post.PostProcess(0);
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

// Sapata corrida (meio modelo) com o TUnsat: base e laterais com rolete, simetria em x = 0 e pressao q na sapata
// (ENeumannU). Apoios impostos pela analise (TPZUnsatAnalysis,
// fixU = {{bcBase, {false, true}}, {bcDireita, {true, false}}, {bcEsquerda, {true, false}}});
// com {bcSapata, {false, true}} tambem em fixU o recalque da sapata passa a ser prescrito e a pressao e ignorada.
TPZMultiphysicsCompMesh *CreateFootingCompMesh(TPZGeoMesh *gmesh, int porder, const UnsatParams &p, REAL q)
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
    for (int id : {bcBase, bcDireita, bcEsquerda})
        mphys->InsertMaterialObject(mat->CreateBC(mat, id, TUnsat::EDirichletU, val1, val2));
    val2[1] = -q; // pressao q da sapata (para baixo); ignorada se a sapata estiver nos apoios (recalque prescrito)
    mphys->InsertMaterialObject(mat->CreateBC(mat, bcSapata, TUnsat::ENeumannU, val1, val2));

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

// Sapata corrida: 5 m x 5 m, meia sapata de 0,5 m com simetria em x = 0
TPZGeoMesh *ReadFootingMesh()
{
    return ReadGmshMesh(std::string(MESH_DIR) + "/sapata.msh", {{"Solo", matSolo}},
                        {{"Base", bcBase}, {"Direita", bcDireita}, {"Superficie", bcSuperficie},
                         {"Simetria", bcEsquerda}, {"Sapata", bcSapata}});
}
