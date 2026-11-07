// LearningMohrCoulomb2: reproducao de Rahardjo et al. (2007), passo a passo (projeto de estudo): estabilidade de um
// talude sob chuva com o material nao saturado u - p_w - p_a (TUnsat) e Mohr-Coulomb (PV2) em deformacao plana.

// geometria
#include "TPZGmshReader.h"
#include "pzgeoel.h"
#include "pzgmesh.h"
#include "pzfstrmatrix.h"
#include "pzpostprocanalysis.h"
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

#include <algorithm>
#include <iomanip>
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

// Solo f50,-5 de Rahardjo et al. (2007), Tabela 1 (E e nu nao estao no artigo: hipotese)
UnsatParams Soilf505()
{
    UnsatParams p;
    //soil properties
    p.E = 2.e4;// Young's modulus (kPa)
    p.nu = 0.3; // Poisson's ratio
    p.c = 10.; // cohesion (kPa)
    p.phi = 26. * M_PI / 180.; // friction angle (rad)
    p.gamma = 20.; // unit weight (kN/m3)
    p.n=0.45; // porosity
    p.ks=1.e-5; // saturated hydraulic conductivity (m/s)
    
    p.hydro.fRetention = TPZUnsaturatedHydraulics::EFredlundXing; // retention model
    p.hydro.fA=50.;// Fredlund-Xing parameter A
    p.hydro.fN=1.; // Fredlund-Xing parameter N
    p.hydro.fM=1.; // Fredlund-Xing parameter M

    p.hydro.fKrw=TPZUnsaturatedHydraulics::EPowerLaw; // relative permeability model for water
    p.hydro.fPower = 4.; // power law exponent
    
    p.hydro.fKra=TPZUnsaturatedHydraulics::EBrooksCorey; // relative permeability model for air

    p.bishop = TUnsat::EChiConstant; // chi constante: phi^b = atan(chi tan phi')
    p.chi = 1.; // phi^b = phi' = 26 graus (artigo)

    p.coupled=false; // esqueleto rigido na percolacao, como no SEEP/W do artigo

    return p;
}

// ids de material: grupos fisicos da malha do Gmsh (meshes/rahardjo_a45.geo)
enum EIds {
    matSolo = 1,       // solo (dominio)
    bcBase = -1,       // base
    bcDireita = -2,    // lateral direita
    bcSuperficie = -3, // superficie livre do terreno
    bcEsquerda = -4    // lateral esquerda
};

TPZGeoMesh *ReadGmshMesh(const std::string &file, const std::map<std::string, int> &domains,
                         const std::map<std::string, int> &boundaries);
// Talude de Rahardjo et al. (2007), Fig. 2: 70 m x 30/20 m, H_s = 10 m, alfa = 45 graus
TPZGeoMesh *ReadSlopeMesh();

TUnsat *UnsatMaterial(int id, const UnsatParams &p);

TPZCompMesh *CreateAtomicMesh(TPZGeoMesh *gmesh, int nstate, int order, const std::set<int> &bcids);
TPZMultiphysicsCompMesh *CreateSlopeCompMesh(TPZGeoMesh *gmesh, int porder, const UnsatParams &p);
// Apaga a malha multifisica e as suas malhas atomicas (depois da analise que as usa)
void DeleteCompMesh(TPZMultiphysicsCompMesh *mphys);

// Norma nas equacoes livres
static REAL NormaLivre(const TPZFMatrix<STATE> &v, const std::map<int64_t, REAL> &fixed, int64_t neq);

// Pos-processamento do estado carregado na malha (solucao e memoria dos pontos de integracao),
// projetado em <nome>_pontos.scal_vec.<step>.vtk. O CSV associa step ao tempo fisico.
void PostProcessVTK(TPZUnsatAnalysis &an, const std::string &nome, int step = 0);

// Arquivos lidos por animacao.py: pressoes e saturacao nos vertices em cada chamada de GravaPressoesCsv e, em cada calculo
// do FS, o ultimo estado convergido do SRM (deslocamentos nos vertices e plasticidade nos pontos de integracao)
struct SaidaAnimacao {
    std::ofstream pw{"snaps_pw.csv"}, pa{"snaps_pa.csv"}, sw{"snaps_sw.csv"}, fs{"fs.csv"};
    std::ofstream srmNos{"srm_nos.csv"}, srmPontos{"srm_pontos.csv"};
    TUnsat *mat = nullptr;
    std::vector<int64_t> pontos; // indice de memoria de cada ponto de integracao, na ordem de pontos.csv
};
void IniciaSaida(SaidaAnimacao &s, TPZUnsatAnalysis *analysis);
void GravaPressoesCsv(SaidaAnimacao &s, REAL t_h, TPZUnsatAnalysis *analysis);
void GravaSRMCsv(SaidaAnimacao &s, REAL t_h, TPZUnsatAnalysis *analysis, const TPZFMatrix<STATE> &xGeo);
REAL calculaFS(const std::string &nome, TPZUnsatAnalysis *analysis, REAL t_h, SaidaAnimacao &s);
int main()
{

    std::cout << "Starting LearningMohrCoulomb2 analysis..." << std::endl;
    UnsatParams soilA = Soilf505();
    TPZGeoMesh * gmesh = ReadSlopeMesh();
    int porder=2;
    TPZMultiphysicsCompMesh * multiphysicscmesh = CreateSlopeCompMesh(gmesh,porder, soilA);
    const std::map<int, std::array<bool, 2>> bcdisplacement = {
        {bcBase,{true,true}},
        {bcDireita,{true,false}},
        {bcEsquerda,{true,false}}
    };
    TPZUnsatAnalysis *analysis= new TPZUnsatAnalysis(multiphysicscmesh,bcdisplacement);
    analysis->SetTolerances(1.e-7, 1.e-6, 40, 25);
    std::cout << "Equacoes: " << multiphysicscmesh->NEquations() << std::endl;

    
    TPZFMatrix<STATE> x = multiphysicscmesh->Solution();//vetor de estado inicial (solucao do problema)
    x.Zero();

    //presssa hidrostatica de agua no no (gw = gamma_w), com o limite de -75 kPa do artigo
    for(const auto &eq:analysis->WaterEq())
    {
        // eq.first: no; eq.second: indice da equacao de agua correspondente ao no
        REAL coordx = analysis->Coord(eq.first,0); 
        REAL coordy = analysis->Coord(eq.first,1); 
        // nivel d'agua (artigo): H_w = 5 m abaixo do pe (40, 20), com inclinacao de 7 graus; sobe para a esquerda
        REAL waterlevel = 15+(40-coordx)*tan(7.*M_PI/180.);
        REAL pw = soilA.gw * (waterlevel - coordy);
        x(eq.second,0) = std::max(pw,-75.);
    }

 

    std::map<int64_t, REAL> contorno = analysis->AllDisplacements(analysis->X());
    for (const auto &e : analysis->AirEq()) contorno[e.second] = 0.;
    for (int id : {bcEsquerda, bcDireita}) {
        for (const auto &b : analysis->BoundaryNodes(id)) {
            const REAL coordx = analysis->Coord(b.first, 0);
            const REAL coordy = analysis->Coord(b.first, 1);
            const REAL waterlevel = 15 + (40 - coordx) * tan(7. * M_PI / 180.);
        if (coordy <= waterlevel)
            contorno[analysis->WaterEq().at(b.first)] = soilA.gw * (waterlevel - coordy);
        }
    }


    REAL hour=3600;//hour in seconds
    REAL dt=0.1*hour;//time step in seconds
    REAL totaltime=48*hour; // 48 hours in seconds
    const REAL chuva = 36. / (1000. * hour); // 36 mm/h em m/s
    int timesteps = static_cast<int>(totaltime / dt);
    

    analysis->X() = x;
    analysis->Restore();
    analysis->Accept();

    SaidaAnimacao arquivos;
    IniciaSaida(arquivos, analysis);
    std::ofstream("summary.txt") << "I 36\ntrain 24\n";
    GravaPressoesCsv(arquivos, 0., analysis);
    calculaFS("fs0", analysis, 0., arquivos);

const REAL tolP = 1.e-6, tolQ = 1.e-10;
const auto superficie = analysis->BoundaryNodes(bcSuperficie);
std::set<int64_t> saida; // nos que estao na face de percolacao

for (int it = 1; it <= timesteps; it++)
{
    const REAL chuvaAgora = (it <= 240) ? chuva : 0.; // chove ate 24 h (com dt = 0.1 h)

    // ponto de partida: a face do passo anterior e os nos com p_w positiva
    auto tentativa = saida;
    for (const auto &b : superficie)
        if (analysis->X()(analysis->WaterEq().at(b.first), 0) > tolP) tentativa.insert(b.first);

    bool fechou = false;
    int iteracoes = 0;
    for (int k = 0; k < 40 && !fechou; k++)
    {
        // na face: p_w = 0; fora dela: fluxo de chuva
        auto fixed = contorno;
        std::map<int64_t, REAL> fluxo;
        for (const auto &b : superficie)
        {
            const int64_t eq = analysis->WaterEq().at(b.first);
            if (tentativa.count(b.first)) fixed[eq] = 0.;
            else fluxo[eq] = chuvaAgora * b.second;
        }

        const auto r = analysis->Solve(dt, fixed, fluxo, false);
        iteracoes += r.iter;
        if (!r.conv) break;

        // confere cada no e monta o novo conjunto
        const auto R = analysis->FullResidual();
        const TPZFMatrix<STATE> sol = multiphysicscmesh->Solution();
        auto nova = tentativa;
        for (const auto &b : superficie)
        {
            const int64_t eq = analysis->WaterEq().at(b.first);
            if (tentativa.count(b.first))
            {
                if (R.GetVal(eq, 0) / (dt * b.second) > chuvaAgora + tolQ) nova.erase(b.first);
            }
            else if (sol.GetVal(eq, 0) > tolP) nova.insert(b.first);
        }

        if (nova == tentativa)
        {
            analysis->Accept();
            saida = tentativa;
            fechou = true;
        }
        else tentativa = nova;
    }
    if (!fechou)
    {
        std::cout << "\nnao convergiu no passo " << it << std::endl;
        break;
    }

    std::cout << "it/timesteps = " << it << "/" << timesteps << " t = " << it * dt / hour << " h, " << iteracoes << " iteracoes, face: " << saida.size() << " nos" << std::endl;
    if (it % 120 == 0)
    {
        GravaPressoesCsv(arquivos, it * dt / hour, analysis);
        calculaFS("fs" + std::to_string(it / 10) + "h", analysis, it * dt / hour, arquivos);
    }
}

    delete analysis;
    DeleteCompMesh(multiphysicscmesh);
    delete gmesh;
    return 0;

}
// FS do estado atual, sem alterar o estado do transiente. O VTK e os CSV guardam o ultimo estado convergido do SRM.
REAL calculaFS(const std::string &nome, TPZUnsatAnalysis *analysis, REAL t_h, SaidaAnimacao &s) {
    const TPZFMatrix<STATE> estado = analysis->X();
    const auto memoria = analysis->SaveMemory();
    TPZFMatrix<STATE> &X = analysis->X();

    // 1. estado sem tensoes, com u = 0 e as pressoes deste instante congeladas
    analysis->ResetMemory();
    X = estado;
    for (int64_t i = 0; i < analysis->NEq(); i++)
        if (!analysis->IsPressure(i)) X(i, 0) = 0.;
    analysis->Restore();
    const auto fixed = analysis->FrozenFixed(X);

    // 2. estado geostatico: peso em 4 incrementos, com a resistencia aumentada
    REAL F = 0.6;
    bool geostatico = true;
    analysis->SetStrengthReductionFactor(F);
    analysis->SetGravityScalesSuction(true);
    for (int k = 1; k <= 4 && geostatico; k++) {
        analysis->SetGravityFactor(REAL(k) / 4.);
        geostatico = analysis->Solve(0., fixed, {}, true).conv;
        if (geostatico) analysis->Accept();
    }
    const TPZFMatrix<STATE> xGeo = X;
    const REAL uRef = NormaLivre(xGeo, fixed, analysis->NEq());

    // 3. caminho de equilibrio, um ArcLengthStepSRM por incremento. Parametros: estudo de sensibilidade,
    // docs/fs_arc_length/relatorio_fs.tex
    TPZUnsatAnalysis::ArcLengthOptions opcoes;
    opcoes.alpha = 0.1 * uRef;
    opcoes.secantPredictor = true;
    opcoes.maxit = 12;
    const REAL lMax = 0.5 * uRef, lMin = 1.e-4 * uRef, tol = 0.01;
    REAL l = 0.01 * uRef, lPedido = l; // lPedido: o l de antes dos cortes deste incremento
    REAL dF = 0., FS = geostatico ? F : -1.;
    TPZFMatrix<STATE> dun; // incremento anterior; vazio no primeiro
    for (int n = 1; geostatico && n <= 100; n++) {
        const auto r = analysis->ArcLengthStepSRM(fixed, F, l, dun, dF, opcoes);
        if (!r.conv) { // volta ao ultimo estado aceito e tenta com l menor
            std::cout << "SRM: l = " << l << " nao convergiu (" << r.iter << " iteracoes)" << std::endl;
            analysis->Restore();
            l *= 0.5;
            if (l < lMin) break;
            continue;
        }
        analysis->Accept();
        dF = r.lambda - F;
        F = r.lambda;
        dun = r.du;
        FS = std::max(FS, F);
        std::cout << "SRM: F = " << F << " dF = " << dF << " l = " << l << " (" << r.iter << " iteracoes)" << std::endl;
        if (dF * lPedido / l < tol) break; // parada: dF medido no l pedido
        l = std::min(lMax, l * std::min(2., std::max(0.5, std::sqrt(10. / r.iter)))); // cerca de 10 iteracoes por incremento
        lPedido = l;
    }

    std::cout << nome << ": FS = " << FS << std::endl;
    s.fs << t_h << ',' << FS << std::endl;
    // FS = -1: o estado geostatico nao convergiu, e a malha e os fatores ficaram como estavam nele
    analysis->SetStrengthReductionFactor(1.);
    analysis->SetGravityFactor(1.);
    if (FS > 0.) {
        GravaSRMCsv(s, t_h, analysis, xGeo);
        PostProcessVTK(*analysis, nome);
    }
    analysis->RestoreMemory(memoria);
    analysis->X() = estado;
    analysis->Restore();
    return FS;
}

// Malha da animacao: vertices.csv (no, x, y), triangles.csv (nos de cada triangulo) e pontos.csv (cada ponto de
// integracao: a linha do seu triangulo em triangles.csv, as coordenadas e a area associada, peso * det J)
void IniciaSaida(SaidaAnimacao &s, TPZUnsatAnalysis *analysis) {
    s.mat = dynamic_cast<TUnsat *>(analysis->Mesh()->FindMaterial(matSolo));
    if (!s.mat) DebugStop();

    std::ofstream fv("vertices.csv");
    fv << std::setprecision(10) << "no,x,y\n";
    for (const auto &e : analysis->WaterEq())
        fv << e.first << ',' << analysis->Coord(e.first, 0) << ',' << analysis->Coord(e.first, 1) << '\n';

    std::ofstream ft("triangles.csv"), fp("pontos.csv");
    ft << "n0,n1,n2\n";
    fp << std::setprecision(10) << "tri,x,y,peso\n";
    int64_t tri = 0;
    for (TPZCompEl *cel : analysis->Mesh()->ElementVec()) {
        if (!cel || !cel->Reference() || cel->Reference()->MaterialId() != matSolo) continue;
        TPZGeoEl *g = cel->Reference();
        ft << g->NodeIndex(0) << ',' << g->NodeIndex(1) << ',' << g->NodeIndex(2) << '\n';
        TPZManVector<int64_t, 10> idx;
        cel->GetMemoryIndices(idx);
        const TPZIntPoints &rule = cel->GetIntegrationRule();
        if (idx.size() != rule.NPoints()) DebugStop();
        TPZManVector<REAL, 3> qsi(2, 0.), x(3, 0.);
        TPZFNMatrix<9, REAL> jac, axes, jacinv;
        REAL w, detjac;
        for (int i = 0; i < rule.NPoints(); ++i) {
            rule.Point(i, qsi, w);
            g->X(qsi, x);
            g->Jacobian(qsi, jac, axes, detjac, jacinv);
            fp << tri << ',' << x[0] << ',' << x[1] << ',' << w * std::fabs(detjac) << '\n';
            s.pontos.push_back(idx[i]);
        }
        tri++;
    }

    s.fs << "t_h,FS" << std::endl;
    s.srmNos << std::setprecision(8) << "t_h,no,ux,uy,ux_geo,uy_geo\n";
    s.srmPontos << std::setprecision(8) << "t_h,ponto,raizJ2_epsp,hardening,tipo\n";
    s.pw << std::setprecision(10);
    s.pa << std::setprecision(10);
    s.sw << std::setprecision(10);
}

// Uma linha de cada snaps_*.csv: t (h) e p_w, p_a (kPa) e S_w nos vertices, na ordem de vertices.csv
void GravaPressoesCsv(SaidaAnimacao &s, REAL t_h, TPZUnsatAnalysis *analysis) {
    s.pw << t_h;
    s.pa << t_h;
    s.sw << t_h;
    for (const auto &e : analysis->WaterEq()) {
        const REAL pw = analysis->X()(e.second, 0), pa = analysis->X()(analysis->AirEq().at(e.first), 0);
        s.pw << ',' << pw;
        s.pa << ',' << pa;
        s.sw << ',' << s.mat->Hydraulics().Evaluate(pw, pa).S;
    }
    s.pw << std::endl;
    s.pa << std::endl;
    s.sw << std::endl;
}

// Ultimo estado convergido do SRM: deslocamentos nos vertices e, nos pontos de integracao, sqrt(J2) da deformacao
// plastica, a variavel de endurecimento e o tipo de resposta (0 elastico, 1 plano principal, 2 canto direito,
// 3 canto esquerdo, -1 apice)
void GravaSRMCsv(SaidaAnimacao &s, REAL t_h, TPZUnsatAnalysis *analysis, const TPZFMatrix<STATE> &xGeo) {
    const bool temGeo = xGeo.Rows() == analysis->NEq(); // o GIM nao tem estado geostatico: u_geo = 0
    for (const auto &e : analysis->WaterEq())
        s.srmNos << t_h << ',' << e.first << ',' << analysis->NodalDisplacement(analysis->X(), e.first, 0) << ','
                 << analysis->NodalDisplacement(analysis->X(), e.first, 1) << ','
                 << (temGeo ? analysis->NodalDisplacement(xGeo, e.first, 0) : 0.) << ','
                 << (temGeo ? analysis->NodalDisplacement(xGeo, e.first, 1) : 0.) << '\n';
    auto mem = s.mat->GetMemory();
    for (size_t i = 0; i < s.pontos.size(); ++i) {
        const auto &st = (*mem)[s.pontos[i]].m_elastoplastic_state;
        TPZTensor<REAL> eps(st.m_eps_p); // distorcoes de engenharia -> componentes do tensor
        eps.XY() *= 0.5;
        eps.XZ() *= 0.5;
        eps.YZ() *= 0.5;
        s.srmPontos << t_h << ',' << i << ',' << std::sqrt(std::max(REAL(0), eps.J2())) << ',' << st.m_hardening << ','
                    << st.m_m_type << '\n';
    }
    s.srmNos.flush();
    s.srmPontos.flush();
}


// Norma nas equacoes livres
static REAL NormaLivre(const TPZFMatrix<STATE> &v, const std::map<int64_t, REAL> &fixed, int64_t neq)
{
    REAL s = 0.;
    for (int64_t i = 0; i < neq; i++)
        if (!fixed.count(i)) s += v.GetVal(i, 0) * v.GetVal(i, 0);
    return std::sqrt(s);
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

// Talude com o TUnsat. bcSuperficie e um contorno de Neumann sem carga: a chuva entra como fluxo nodal, em main.
// Espacos: u de ordem porder (2 = triangulo de 6 nos), p_w e p_a de ordem porder - 1 (1: P1).
// Apoios impostos pela analise (TPZUnsatAnalysis,
// fixU = {{bcBase, {true, true}}, {bcDireita, {true, false}}, {bcEsquerda, {true, false}}}):
// as condicoes EDirichletU do TUnsat nao contribuem, so marcam os contornos.
TPZMultiphysicsCompMesh *CreateSlopeCompMesh(TPZGeoMesh *gmesh, int porder, const UnsatParams &p)
{
    const std::set<int> bcids = {bcBase, bcDireita, bcEsquerda, bcSuperficie};
    TPZManVector<TPZCompMesh *, 3> meshvec = {CreateAtomicMesh(gmesh, 2, porder, bcids),
                                              CreateAtomicMesh(gmesh, 1, porder - 1, bcids),
                                              CreateAtomicMesh(gmesh, 1, porder - 1, bcids)};
    TUnsat *mat = UnsatMaterial(matSolo, p);
    auto *mphys = new TPZMultiphysicsCompMesh(gmesh);
    mphys->SetDimModel(2);
    mphys->InsertMaterialObject(mat);

    TPZFMatrix<STATE> val1(2, 2, 0.);
    TPZManVector<STATE, 2> val2(2, 0.);
    for (int id : bcids) mphys->InsertMaterialObject(mat->CreateBC(mat, id,
        id == bcSuperficie ? TUnsat::ENeumannU : TUnsat::EDirichletU, val1, val2));

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

// Talude de Rahardjo et al. (2007), Fig. 2: 70 m x 30/20 m, H_s = 10 m, alfa = 45 graus
TPZGeoMesh *ReadSlopeMesh()
{
    return ReadGmshMesh(std::string(MESH_DIR) + "/rahardjo_a45.msh", {{"Solo", matSolo}},
                        {{"Base", bcBase}, {"Direita", bcDireita}, {"Superficie", bcSuperficie}, {"Esquerda", bcEsquerda}});
}
