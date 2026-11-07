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

#include <algorithm>
#include <iomanip>
#include <array>
#include <cmath>
#include <fstream>
#include <functional>
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

// Solo f10,-4 da Tabela 1: so a e k_s mudam em relacao ao f50,-5
UnsatParams Soilf104()
{
    UnsatParams p = Soilf505();
    p.hydro.fA = 10.; // Fredlund-Xing parameter A (kPa)
    p.ks = 1.e-4; // saturated hydraulic conductivity (m/s)
    return p;
}

// Solo f100,-6 da Tabela 1: so a e k_s mudam em relacao ao f50,-5
UnsatParams Soilf1006()
{
    UnsatParams p = Soilf505();
    p.hydro.fA = 100.; // Fredlund-Xing parameter A (kPa)
    p.ks = 1.e-6; // saturated hydraulic conductivity (m/s)
    return p;
}

// ids de material: grupos fisicos das malhas do Gmsh (meshes/rahardjo_a45.geo e meshes/sapata.geo)
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
// Talude de Rahardjo et al. (2007), Fig. 2: 70 m x 30/20 m, H_s = 10 m, alfa = 45 graus
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
// xGeo (opcional) recebe o estado geostatico, de onde parte a reducao da resistencia
REAL FactorSRM(TPZUnsatAnalysis &an, const TPZFMatrix<STATE> &x, TPZFMatrix<STATE> *xGeo = nullptr, REAL F0 = 0.6, REAL tol = 0.01,
               const std::vector<REAL> &Fs = {0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5, 3.0, 3.5, 4.0});
// Parametros do caminho de equilibrio com arc-length (FactorSRMArc e FactorGIM). Os comprimentos sao relativos a |u_ref|:
// no SRM, o deslocamento do estado geostatico; no GIM, o do primeiro passo. Os valores vem do estudo de sensibilidade
// (relatorio_fs.tex, em ~/projects/runs/arc_length_sensibilidade); ParamsSRM traz o que muda no SRM.
struct ParamsArco {
    std::vector<REAL> lista;  // SRM com alpha = 0: F prescritos antes do arco
    REAL F0 = 0.6;            // SRM: F do estado geostatico
    REAL F0Min = 0.4;         // SRM: F do estado geostatico se F0 nao convergir (FS abaixo de F0)
    REAL lam0 = 0.5;          // GIM: lambda do primeiro passo, prescrito
    int criterio = 0;         // parada. 0: dF < tol, com dF medido no l pedido; 1: (dF / l) * s < tol, s = |u - u_ref|
    REAL tol = 0.001;
    int alvo = 10;            // iteracoes de Newton desejadas por incremento: l varia com sqrt(alvo / iteracoes)
    REAL cresceMax = 2.;      // maior aumento de l de um incremento para o seguinte
    REAL reduzMax = 0.5;      // maior reducao de l de um incremento para o seguinte
    REAL corte = 0.5;         // reducao de l quando o incremento falha
    REAL l0 = 0.25;           // l inicial / |u_ref|, se nao ha incremento anterior
    REAL l0Inc = 0.5;         // l inicial / |ultimo incremento prescrito|
    REAL lMax = 2.;           // maior l / |u_ref|
    REAL lMin = 1.e-4;        // menor l / |u_ref|: abaixo dele o caminho termina
    REAL alpha = 0.;          // alpha / |u_ref|, na restricao Du^T Du + alpha^2 dF^2 = l^2 (0: cilindrica)
    bool secante = true;      // preditor pela secante, em vez de pela tangente
    int maxIt = 12;           // limite de iteracoes de Newton por incremento
    int maxPassos = 100;      // limite de incrementos
    REAL fatorMax = 10.;      // acima deste fator o caminho termina: FS >= fatorMax (GIM em talude pouco inclinado nao tem limite)
};
inline ParamsArco ParamsSRM()
{
    ParamsArco p;
    p.alpha = 0.1; // sem lista: dF <= l / (alpha |u_ref|)
    p.l0 = 0.01;   // dF de cerca de 0,1 no trecho elastico
    p.lMax = 0.5;
    return p;
}
inline ParamsArco ParamsGIM() { return ParamsArco(); }
// Resumo do caminho: incrementos aceitos, cortes de l, iteracoes de Newton e por que terminou
struct SaidaArco {
    int incrementos = 0, cortes = 0, iteracoes = 0;
    REAL uRef = 0., lFinal = 0., inclinacao = 0., s = 0.;
    std::string motivo = "maxPassos";
};
// Fator de seguranca por reducao de resistencia com o arc-length, em vez da bissecao: o maior F do caminho de equilibrio.
// xGeo (opcional) recebe o estado geostatico
REAL FactorSRMArc(TPZUnsatAnalysis &an, const TPZFMatrix<STATE> &x, TPZFMatrix<STATE> *xGeo = nullptr, const ParamsArco &p = ParamsSRM(),
                  SaidaArco *saida = nullptr);
// Fator de seguranca por aumento da gravidade com as pressoes de x congeladas (seco: x com pressoes nulas): o maior
// lambda do caminho de equilibrio, tracado com o arc-length cilindrico
REAL FactorGIM(TPZUnsatAnalysis &an, const TPZFMatrix<STATE> &x, const ParamsArco &p = ParamsGIM(), SaidaArco *saida = nullptr);
// Fator de seguranca por aumento da carga de superficie do contorno bcid (sapata) com as pressoes de x congeladas
REAL FactorLoad(TPZUnsatAnalysis &an, const TPZFMatrix<STATE> &x, int bcid, REAL tol = 0.01,
                const std::vector<REAL> &lams = {0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.0});

// Pos-processamento do estado carregado na malha (solucao e memoria dos pontos de integracao),
// projetado em <nome>_pontos.scal_vec.<step>.vtk. O CSV associa step ao tempo fisico.
void PostProcessVTK(TPZUnsatAnalysis &an, const std::string &nome, int step = 0);

// Funcao da carga externa do reservatorio: nivel e lido por referencia a cada avaliacao.
// O objeto nivel e o material devem existir enquanto a CC estiver em uso.
auto CargaReservatorio(const REAL &nivel, TUnsat *reservatorio)
{
    return [&nivel, reservatorio](const TPZVec<REAL> &x, TPZVec<STATE> &tracao,
                                  TPZFMatrix<STATE> &) {
        const REAL pressao = reservatorio->GravityFactor() * reservatorio->GammaW()
                           * std::max(nivel-x[1], REAL(0));
        // talude.geo: face a 45 graus; patamares horizontais na crista e no pe.
        const bool face = x[0] > 35. && x[0] < 45. && x[1] > 30. && x[1] < 40.;
        tracao.Resize(2);
        tracao[0] = face ? -pressao/std::sqrt(2.) : 0.;
        tracao[1] = face ? -pressao/std::sqrt(2.) : -pressao;
    };
}

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
    TPZFMatrix<STATE> xGeo;
    //const REAL FS = FactorSRM(*analysis, estado, &xGeo);
    const REAL FS = FactorSRMArc(*analysis, estado, &xGeo);
    // const REAL FS = FactorGIM(*analysis, estado);

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


// Fator de seguranca por reducao de resistencia (SRM) com as pressoes de x congeladas (so o equilibrio e
// resolvido): serve as tres etapas, com x = pressoes nulas (seco), p_w de um instante (com agua) ou p_w e p_a de um
// instante (com agua e ar). Parte do estado sem tensoes (memoria inicial, u = 0); estado geostatico com a
// resistencia aumentada (F0 < 1, peso e sucao em 4 incrementos); depois F cresce pela lista Fs, sempre a partir do
// ultimo estado convergido; no primeiro F sem convergencia, bissecao entre o ultimo F convergido e ele ate a
// diferenca ser menor que tol. Entao, nova tentativa do F que falhou, agora com incremento < tol: se convergir, a
// falha era do tamanho do salto (o Newton pode divergir num salto grande mesmo havendo equilibrio) e a busca continua
// pela lista; se falhar de novo, e a ruptura. Retorna o ultimo F convergido (-1 se o estado geostatico nao convergir).
REAL FactorSRM(TPZUnsatAnalysis &an, const TPZFMatrix<STATE> &x, TPZFMatrix<STATE> *xGeo, REAL F0, REAL tol,
               const std::vector<REAL> &Fs)
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
            std::cout << "\nSRM: estado geostatico nao convergiu" << std::endl;
            return -1.;
        }
        an.Accept();
    }
    if (xGeo) *xGeo = an.X();

    // reducao da resistencia
    REAL Fok = F0, Ffail = -1.;
    size_t next = 0;
    while (true) {
        REAL F;
        bool retry = false;
        if (Ffail > 0.) {
            if (Ffail - Fok < tol) {
                F = Ffail;
                retry = true;
            } else {
                F = 0.5 * (Fok + Ffail);
            }
        } else {
            while (next < Fs.size() && Fs[next] <= Fok) next++;
            if (next == Fs.size()) break;
            F = Fs[next];
        }
        an.SetStrengthReductionFactor(F);
        const auto r = an.Solve(0., fixed, {}, true);
        std::cout << "\nSRM: F = " << F << (r.conv ? " convergiu" : " nao convergiu") << " (" << r.iter << " iteracoes)"
                  << (retry ? " [nova tentativa]" : "") << std::endl;
        if (r.conv) {
            an.Accept();
            Fok = F;
            if (retry) Ffail = -1.;
        } else {
            an.Restore();
            if (retry) break;
            Ffail = F;
        }
    }
    an.SetStrengthReductionFactor(1.);
    return Fok;
}

// Norma nas equacoes livres
static REAL NormaLivre(const TPZFMatrix<STATE> &v, const std::map<int64_t, REAL> &fixed, int64_t neq)
{
    REAL s = 0.;
    for (int64_t i = 0; i < neq; i++)
        if (!fixed.count(i)) s += v.GetVal(i, 0) * v.GetVal(i, 0);
    return std::sqrt(s);
}

// Traca o caminho de equilibrio com arc-length a partir do estado convergido X() com fator lam, ate a parada. O comprimento
// l se ajusta para cerca de p.alvo iteracoes por incremento, entre lMin e lMax; se o incremento falha, l e cortado.
//  - passo: um incremento de arc-length (fator, l, incremento anterior, variacao anterior do fator);
//  - vertical: sem ponto plastico e com restricao cilindrica; deve avancar o fator com valor prescrito (false: termina).
// Parada: p.criterio. dF e medido no l da primeira tentativa do incremento: um corte de l nao parece um caminho plano.
// Termina tambem em p.fatorMax, em l < lMin e em p.maxPassos, e entao avisa. Retorna o maior fator do caminho.
static REAL SegueCaminho(TPZUnsatAnalysis &an, const std::map<int64_t, REAL> &fixed, REAL lam, TPZFMatrix<STATE> dun, REAL dlamn, REAL l,
                         REAL uRef, const TPZFMatrix<STATE> &xRef, const ParamsArco &p, SaidaArco &saida, const std::string &nome,
                         const std::function<TPZUnsatAnalysis::ArcLengthResult(REAL, REAL, const TPZFMatrix<STATE> &, REAL)> &passo,
                         const std::function<bool(REAL &, TPZFMatrix<STATE> &, REAL &, REAL &)> &vertical)
{
    const REAL lMax = p.lMax * uRef, lMin = p.lMin * uRef;
    l = std::min(lMax, std::max(lMin, l));
    saida.uRef = uRef;
    REAL lamMax = lam;
    REAL lPedido = l; // l da primeira tentativa do incremento
    for (int n = 1; n <= p.maxPassos; n++) {
        const auto r = passo(lam, l, dun, dlamn);
        saida.iteracoes += r.iter;
        saida.lFinal = l;
        if (!r.conv) {
            std::cout << "\n" << nome << ": l = " << l << (r.noRoot ? " sem raiz real" : (r.vertical ? " sem ponto plastico" : " nao convergiu"))
                      << " (" << r.iter << " iteracoes)" << std::endl;
            an.Restore();
            if (r.vertical) {
                if (!vertical(lam, dun, dlamn, l)) {
                    saida.motivo = "vertical";
                    break;
                }
                lamMax = std::max(lamMax, lam);
                l = std::min(lMax, std::max(lMin, l));
                lPedido = l;
                continue;
            }
            saida.cortes++;
            l *= p.corte;
            if (l < lMin) {
                saida.motivo = "lMin";
                break;
            }
            continue;
        }
        an.Accept();
        saida.incrementos++;
        dlamn = r.lambda - lam;
        const REAL rigidez = dlamn / l;
        lam = r.lambda;
        dun = r.du;
        lamMax = std::max(lamMax, lam);
        TPZFMatrix<STATE> distancia(an.X());
        distancia -= xRef;
        saida.s = NormaLivre(distancia, fixed, an.NEq());
        saida.inclinacao = rigidez;
        std::cout << "\n" << nome << ": fator = " << lam << " d = " << dlamn << " l = " << l << " (" << r.iter << " iteracoes)" << std::endl;
        const REAL medida = p.criterio == 0 ? dlamn * lPedido / l : rigidez * saida.s;
        if (medida < p.tol) {
            saida.motivo = "parada";
            break;
        }
        if (lam >= p.fatorMax) {
            saida.motivo = "fatorMax";
            break;
        }
        const REAL fator = std::min(p.cresceMax, std::max(p.reduzMax, std::sqrt(REAL(p.alvo) / r.iter)));
        l = std::min(lMax, l * fator);
        lPedido = l;
    }
    if (saida.motivo != "parada")
        std::cout << "\n" << nome << ": terminou sem atingir o criterio de parada (" << saida.motivo << ")" << std::endl;
    return lamMax;
}

// Fator de seguranca por reducao de resistencia com o arc-length (TPZUnsatAnalysis::ArcLengthStepSRM), com as pressoes
// de x congeladas. Depois do estado geostatico, SegueCaminho. Com poucos pontos plasticos o caminho e quase vertical (u
// quase nao muda com F), e sem nenhum e vertical; com alpha > 0, o padrao, isso nao e problema: dF fica limitado por
// l / (alpha |u_ref|) e a lista fica vazia. So com alpha = 0: F sobe antes pela p.lista, com valor prescrito, e, se o
// arco achar o caminho vertical, com valor prescrito de novo: bissecao com o que falhou, ou 25% se nada falhou.
// Retorna o maior F do caminho, ou -1 se o estado geostatico nao convergir.
REAL FactorSRMArc(TPZUnsatAnalysis &an, const TPZFMatrix<STATE> &x, TPZFMatrix<STATE> *xGeo, const ParamsArco &p, SaidaArco *saida)
{
    SaidaArco resumo;
    TPZFMatrix<STATE> &X = an.X();
    std::map<int64_t, REAL> fixed; // apoios = 0 e todas as pressoes nos valores de x
    // estado geostatico a partir do estado sem tensoes: peso (e sucao) em incrementos, resistencia aumentada
    auto geostatico = [&](REAL F0) {
        an.ResetMemory();
        X = x;
        for (int64_t i = 0; i < an.NEq(); i++)
            if (!an.IsPressure(i)) X(i, 0) = 0.;
        an.Restore();
        fixed = an.FrozenFixed(X);
        an.SetStrengthReductionFactor(F0);
        an.SetGravityScalesSuction(true);
        const int n = 4;
        for (int k = 1; k <= n; k++) {
            an.SetGravityFactor(REAL(k) / n);
            if (!an.Solve(0., fixed, {}, true).conv) return false;
            an.Accept();
        }
        return true;
    };
    REAL F0 = p.F0;
    bool convergiu = geostatico(F0);
    if (!convergiu && p.F0Min < F0) {
        F0 = p.F0Min;
        convergiu = geostatico(F0);
    }
    if (!convergiu) {
        std::cout << "\nSRM: estado geostatico nao convergiu" << std::endl;
        resumo.motivo = "geostatico";
        if (saida) *saida = resumo;
        return -1.;
    }
    if (xGeo) *xGeo = an.X();
    const TPZFMatrix<STATE> xRef = an.X();
    const REAL uRef = NormaLivre(xRef, fixed, an.NEq());

    REAL F = F0, Ffalhou = -1., dFn = 0.;
    TPZFMatrix<STATE> dun; // incremento anterior; vazio se ainda nao houve incremento
    // F prescrito: aceita o estado da malha e guarda o incremento. Retorna false se a Newton nao convergir
    auto prescrito = [&](REAL Fi) {
        an.SetStrengthReductionFactor(Fi);
        const auto r = an.Solve(0., fixed, {}, true);
        std::cout << "\nSRM: F = " << Fi << (r.conv ? " convergiu" : " nao convergiu") << " (" << r.iter << " iteracoes)" << std::endl;
        resumo.iteracoes += r.iter;
        if (!r.conv) {
            an.Restore();
            return false;
        }
        dun = X;
        an.Accept();
        dun *= -1.;
        dun += X;
        dFn = Fi - F;
        F = Fi;
        return true;
    };
    const int itAntes = an.MaxItFrozen();
    an.SetMaxItFrozen(p.maxIt);
    for (const REAL Fi : p.lista) {
        if (Fi <= F) continue;
        if (!prescrito(Fi)) {
            Ffalhou = Fi;
            break;
        }
    }
    auto comprimento = [&]() {
        return dun.Rows() ? p.l0Inc * NormaLivre(dun, fixed, an.NEq()) : p.l0 * uRef;
    };

    TPZUnsatAnalysis::ArcLengthOptions opcoes;
    opcoes.alpha = p.alpha * uRef;
    opcoes.secantPredictor = p.secante;
    opcoes.maxit = p.maxIt;
    const REAL Fmax = SegueCaminho(an, fixed, F, dun, dFn, comprimento(), uRef, xRef, p, resumo, "SRM",
        [&](REAL Fn, REAL l, const TPZFMatrix<STATE> &d, REAL dF) { return an.ArcLengthStepSRM(fixed, Fn, l, d, dF, opcoes); },
        [&](REAL &Fn, TPZFMatrix<STATE> &d, REAL &dF, REAL &l) {
            // sem ponto plastico: u nao depende de F
            if (Ffalhou > 0. && Ffalhou - F < p.tol) return false;
            const REAL Fi = Ffalhou > 0. ? 0.5 * (F + Ffalhou) : 1.25 * F;
            if (prescrito(Fi)) l = comprimento();
            else Ffalhou = Fi;
            Fn = F;
            d = dun;
            dF = dFn;
            return true;
        });
    an.SetMaxItFrozen(itAntes);
    an.SetStrengthReductionFactor(1.);
    if (saida) *saida = resumo;
    return Fmax;
}

// Fator de seguranca por aumento da gravidade (GIM) com as pressoes de x congeladas: parte do estado sem tensoes,
// com a resistencia real (F = 1): lambda multiplica o peso e a parte positiva da pressao equivalente;
// a sucao fica com o valor de x. CargaReservatorio tambem usa lambda na tracao externa.
// O caminho de equilibrio e tracado com o arc-length (de Souza Neto, Peric e Owen 2008, Box 4.4;
// TPZUnsatAnalysis::ArcLengthStep), que nao depende de a Newton convergir com lambda prescrito perto do limite. O primeiro
// passo tem lambda prescrito (p.lam0, reduzido se nao convergir): a sucao entra toda nele, entao ele nao e proporcional, e
// o sinal do primeiro incremento nao vem do criterio da secante: lambda cresce. Depois, SegueCaminho. Retorna o maior
// lambda do caminho.
REAL FactorGIM(TPZUnsatAnalysis &an, const TPZFMatrix<STATE> &x, const ParamsArco &p, SaidaArco *saida)
{
    SaidaArco resumo;
    // estado sem tensoes com as pressoes de x
    an.ResetMemory();
    TPZFMatrix<STATE> &X = an.X();
    X = x;
    for (int64_t i = 0; i < an.NEq(); i++)
        if (!an.IsPressure(i)) X(i, 0) = 0.;
    an.Restore();
    const auto fixed = an.FrozenFixed(X); // apoios = 0 e todas as pressoes nos valores de x
    const TPZFMatrix<STATE> xRef = X;

    an.SetStrengthReductionFactor(1.);
    an.SetGravityScalesSuction(false); // peso e pressao equivalente positiva escalados; sucao preservada
    const int itAntes = an.MaxItFrozen();
    an.SetMaxItFrozen(p.maxIt);
    auto termina = [&](REAL resultado) {
        an.SetMaxItFrozen(itAntes);
        an.SetGravityScalesSuction(true);
        an.SetGravityFactor(1.);
        an.SetStrengthReductionFactor(1.);
        if (saida) *saida = resumo;
        return resultado;
    };

    // primeiro passo, com lambda prescrito
    REAL lam = p.lam0;
    while (true) {
        an.SetGravityFactor(lam);
        const auto r = an.Solve(0., fixed, {}, true);
        std::cout << "\nGIM: lambda = " << lam << (r.conv ? " convergiu" : " nao convergiu") << " (" << r.iter << " iteracoes)" << std::endl;
        resumo.iteracoes += r.iter;
        if (r.conv) break;
        an.Restore();
        lam *= 0.5;
        if (lam < 0.01) {
            resumo.motivo = "lam0";
            return termina(0.);
        }
    }
    an.Accept();

    // comprimento inicial e carga por unidade de lambda
    const REAL uRef = NormaLivre(X, fixed, an.NEq());
    const TPZFMatrix<STATE> fbar = an.LoadVector(X);
    an.Restore();

    TPZUnsatAnalysis::ArcLengthOptions opcoes;
    opcoes.alpha = p.alpha * uRef;
    opcoes.secantPredictor = p.secante;
    opcoes.maxit = p.maxIt;
    const REAL lamMax = SegueCaminho(an, fixed, lam, TPZFMatrix<STATE>(), 0., p.l0 * uRef, uRef, xRef, p, resumo, "GIM",
        [&](REAL ln, REAL l, const TPZFMatrix<STATE> &d, REAL dl) { return an.ArcLengthStep(fixed, fbar, ln, l, d, dl, opcoes); },
        [](REAL &, TPZFMatrix<STATE> &, REAL &, REAL &) { return false; });
    return termina(lamMax);
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

// Talude com o TUnsat: mesmos apoios; bcSuperficie recebe a pressao externa do reservatorio.
// Espacos: u de ordem porder (2 = triangulo de 6 nos), p_w e p_a de ordem porder - 1 (1: P1).
// Apoios impostos pela analise (TPZUnsatAnalysis,
// fixU = {{bcBase, {true, true}}, {bcDireita, {true, false}}, {bcEsquerda, {true, false}}}):
// as condicoes EDirichletU do TUnsat nao contribuem, so marcam os contornos.
TPZMultiphysicsCompMesh *CreateSlopeCompMesh(TPZGeoMesh *gmesh, int porder, const UnsatParams &p)
{
    const std::set<int> bcids = {bcBase, bcDireita, bcEsquerda, bcSuperficie}; // ADICAO: reservatorio
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

// Sapata corrida: 5 m x 5 m, meia sapata de 0,5 m com simetria em x = 0
TPZGeoMesh *ReadFootingMesh()
{
    return ReadGmshMesh(std::string(MESH_DIR) + "/sapata.msh", {{"Solo", matSolo}},
                        {{"Base", bcBase}, {"Direita", bcDireita}, {"Superficie", bcSuperficie},
                         {"Simetria", bcEsquerda}, {"Sapata", bcSapata}});
}
