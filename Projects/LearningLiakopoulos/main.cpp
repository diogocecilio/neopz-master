// LearningLiakopoulos: dessaturacao de uma coluna de areia, passo a passo (projeto de estudo, nos moldes do
// LearningMohrCoulomb). Ensaio de drenagem de Liakopoulos (1965) como na Secao 4.1 de Callari e Abati (2009) e na
// Secao 10 do resumo "Solos nao saturados: formulacao u-p_w-p_a", com o material nao saturado u - p_w - p_a (TUnsat)
// e ar passivo (p_a = 0). Etapas: malha geometrica (Gmsh), malha computacional (Taylor-Hood P2-P1), matriz
// estrutural e solver, analise incremental (TPZUnsatAnalysis) e pos-processamento (VTK e CSV).
// Casos (Secao 10.7): coluna parcialmente saturada com os passos da referencia 1D, com metade do passo e com os
// passos do artigo; coluna saturada (c_L = 0). Saidas: historico.csv e perfis.csv, comparados em Python
// (plot_liakopoulos.py) com a referencia 1D (coluna_liakopoulos_1d.py), com a solucao analitica da coluna saturada
// (Eq. 48 do resumo) e com as curvas digitalizadas de Callari e Abati (2009, Figs. 3 e 4).

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
#include "Plasticity/TPZUnsaturatedHydraulics.h"
#include "Plasticity/TPZYCMohrCoulombPV2.h"

// matriz estrutural e solver
#include "TPZSkylineNSymStructMatrix.h"
#include "TPZSpStructMatrix.h"
#include "pzstepsolver.h"

// analise incremental (Newton, filtro de equacoes, memoria do estado convergido)
#include "TPZUnsatAnalysis.h"

// pos-processamento
#include "TPZVTKGeoMesh.h"
#include "pzfstrmatrix.h"
#include "pzpostprocanalysis.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <utility>
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

// Coluna de Liakopoulos (Callari e Abati 2009, Tabela 1 e Eq. 58; Secao 10.3 do resumo): areia Del Monte, ar
// passivo. saturada = true e a coluna hipotetica saturada: c_L = 0 da S_w = 1, C_s = 0 e k_rw = 1
UnsatParams LiakopoulosParams(bool saturada)
{
    UnsatParams p;
    p.E = 1300., p.nu = 0.4;               // E_sk (kPa) e nu_sk
    p.c = 1.e5, p.phi = 30. * M_PI / 180.; // Mohr-Coulomb com c' muito alto: o esqueleto fica elastico
    p.gamma = 1700. * 9.81e-3;             // gamma = rho g = 16,677 kN/m3, constante
    p.n = 0.2975;                          // porosidade
    p.Kw = 2.e6;                           // modulo volumetrico da agua (kPa)
    p.gw = 9.81;                           // peso especifico da agua (kN/m3)
    p.ks = 9.81 * 4.5e-7;                  // k_s = gamma_w k_w^sat = 4,4145e-6 m/s (k_w^sat = 4,5e-7 m2/(kPa s))
    p.alpha = 1., p.Ks = 0.;               // Biot b = 1 e graos incompressiveis (1/N = 0)
    p.bishop = TUnsat::EChiSaturation;     // chi = S_w: com ar passivo, p_B = S_w p_w
    // retencao de Liakopoulos: S_w = 1 - c_L p_c^b_L com p_c em kPa (no artigo, p_c em Pa e c_L = 1,9722e-11)
    p.hydro.fRetention = TPZUnsaturatedHydraulics::ELiakopoulos;
    p.hydro.fLB = 2.4279;
    p.hydro.fLC = saturada ? 0. : 1.9722e-11 * std::pow(1000., p.hydro.fLB); // 3,7901e-4 kPa^-b_L
    p.hydro.fSwr = 0.; // o limite inferior nao atua nesta faixa (S_w > 0,91)
    // permeabilidade relativa da agua: k_rw = 1 - d_L (1 - S_w)^e_L (a do ar nao entra: ar passivo)
    p.hydro.fKrw = TPZUnsaturatedHydraulics::ELiakopoulosW;
    p.hydro.fLD = 2.207, p.hydro.fLE = 1.0121;
    return p;
}

// Ensaio: geometria da coluna, vertical dos perfis, duracao e instantes dos perfis
struct ColumnTest {
    REAL W = 0.10;        // largura (m)
    REAL H = 1.00;        // altura (m)
    REAL xProfile = 0.05; // vertical dos perfis de p_w e u_y (m)
    REAL tEnd = 7500.;    // fim do transiente (s)
    // perfis de p_w em 600, 1200 e 1800 s (Fig. 3b), de u_y em 300, 1200 e 7200 s (Fig. 4b) e instantes da tabela
    // da Secao 10.6; um perfil so e gravado se o passo cair exatamente no instante
    std::vector<REAL> tProfiles = {300., 600., 1200., 1800., 3600., 6400., 7200., 7500.};
};

// Passos de tempo por fases: dt = phases[i].second enquanto t < phases[i].first; o ultimo passo termina em tEnd
struct TimeSteps {
    std::vector<std::pair<REAL, REAL>> phases; // (t limite, dt) em s
    REAL Dt(REAL t, REAL tEnd) const
    {
        REAL dt = phases.back().second;
        for (auto &[tl, d] : phases)
            if (t < tl - 1.e-9) {
                dt = d;
                break;
            }
        return std::min(dt, tEnd - t);
    }
};

// Passos da referencia 1D (Secao 10.6): 1 s ate 60 s, 5 s ate 600 s, 10 s ate 2000 s e 20 s depois (vezes factor)
TimeSteps ReferenceSteps(REAL factor = 1.)
{
    const REAL inf = std::numeric_limits<REAL>::max();
    return {{{60., factor * 1.}, {600., factor * 5.}, {2000., factor * 10.}, {inf, factor * 20.}}};
}

// Passos sugeridos pelos marcadores das Figs. 3a e 4a do artigo: 10 s ate 300 s, 100 s ate 2400 s e 1000 s depois
TimeSteps ArticleSteps()
{
    const REAL inf = std::numeric_limits<REAL>::max();
    return {{{300., 10.}, {2400., 100.}, {inf, 1000.}}};
}

// ids de material: grupos fisicos da malha do Gmsh (meshes/coluna.geo)
enum EIds {
    matSolo = 1,    // solo (dominio)
    bcBase = -1,    // base (y = 0): engastada e drenante
    bcDireita = -2, // lateral direita (x = 0,10 m): rolete e impermeavel
    bcTopo = -3,    // topo (y = 1,00 m): livre e impermeavel
    bcEsquerda = -4 // lateral esquerda (x = 0): rolete e impermeavel
};

TPZGeoMesh *ReadGmshMesh(const std::string &file, const std::map<std::string, int> &domains,
                         const std::map<std::string, int> &boundaries);
// Coluna: 0,10 m x 1,00 m, 4 x 40 quadrados de 0,025 m divididos em dois triangulos (320 elementos)
TPZGeoMesh *ReadColumnMesh();

TUnsat *UnsatMaterial(int id, const UnsatParams &p);

TPZCompMesh *CreateAtomicMesh(TPZGeoMesh *gmesh, int nstate, int order, const std::set<int> &bcids);
TPZMultiphysicsCompMesh *CreateColumnCompMesh(TPZGeoMesh *gmesh, int porder, const UnsatParams &p);
// Apaga a malha multifisica e as suas malhas atomicas (depois da analise que as usa)
void DeleteCompMesh(TPZMultiphysicsCompMesh *mphys);

// Analise incremental do modelo u - p_w - p_a com a matriz estrutural nao simetrica e o solver direto (LU)
TPZUnsatAnalysis *CreateAnalysis(TPZMultiphysicsCompMesh *cmesh);

// Estado inicial: p_w = p_a = 0 em todo o dominio; peso com as pressoes congeladas (dt = 0) e deslocamentos zerados
void InitialState(TPZUnsatAnalysis &an);

// Nos de observacao: vertical x = xProfile (y crescente), seus nos no meio (y = H/2) e no topo, e os nos da base
struct Probes {
    std::vector<int64_t> line;
    int64_t mid = -1, top = -1;
    std::map<int64_t, REAL> base; // no -> comprimento de influencia (m)
};
Probes FindProbes(TPZUnsatAnalysis &an, const ColumnTest &test);

// Vazao que sai pela base por unidade de area (m/s), pelas reacoes da equacao da agua (antes do Accept)
REAL BaseOutflow(TPZUnsatAnalysis &an, const Probes &pr, REAL dt, REAL W);

// Linha do historico.csv: passo, t, dt, iteracoes de Newton, vazao na base, recalque do topo, p_w no meio e no topo
// da vertical dos perfis e numero de pontos plasticos
void WriteHistory(std::ofstream &csv, const std::string &caso, int step, REAL t, REAL dt, int newton, REAL qb,
                  TPZUnsatAnalysis &an, const Probes &pr);
// Linhas do perfis.csv: p_w, u_y e S_w nos nos da vertical dos perfis
void WriteProfile(std::ofstream &csv, const std::string &caso, REAL t, TPZUnsatAnalysis &an, const Probes &pr);

// Pos-processamento VTK do estado carregado na malha (ultimo estado convergido): variaveis da solucao e dos pontos
// de integracao projetadas por elemento em <nome>_pontos.scal_vec.<step>.vtk
void PostProcessVTK(TPZUnsatAnalysis &an, const std::string &nome, int step = 0);

// Um caso do ensaio na malha geometrica gmesh: estado inicial e transiente ate tEnd, com as linhas dos CSVs gravadas
// com o nome do caso. Retorna false se algum passo nao convergir
bool RunCase(TPZGeoMesh *gmesh, const std::string &caso, bool saturada, const TimeSteps &steps, bool vtk,
             const ColumnTest &test, std::ofstream &hist, std::ofstream &prof);

int main(int argc, char *argv[])
{
    const ColumnTest test;

    // casos (Secao 10.7, itens 6 e 7): coluna parcialmente saturada com os passos da referencia 1D, com metade deles
    // e com os passos do artigo; coluna saturada (c_L = 0), comparada com a solucao analitica de Terzaghi (Eq. 48)
    struct Case {
        std::string nome;
        bool saturada;
        TimeSteps steps;
        bool vtk; // grava o VTK nos instantes dos perfis
    };
    const std::vector<Case> casos = {{"parcial", false, ReferenceSteps(), true},
                                     {"parcial_meio_passo", false, ReferenceSteps(0.5), false},
                                     {"parcial_passos_artigo", false, ArticleSteps(), false},
                                     {"saturada", true, ReferenceSteps(), true}};
    // casos escolhidos na linha de comando (sem argumentos: todos), por exemplo: ./LearningLiakopoulos parcial saturada
    const std::set<std::string> escolhidos(argv + 1, argv + argc);
    for (const std::string &nome : escolhidos)
        if (std::none_of(casos.begin(), casos.end(), [&](const Case &c) { return c.nome == nome; })) {
            std::cout << "caso desconhecido: " << nome
                      << " (casos: parcial, parcial_meio_passo, parcial_passos_artigo, saturada)" << std::endl;
            return 1;
        }

    // ---- malha geometrica ----
    TPZGeoMesh *gmesh = ReadColumnMesh();
    std::ofstream out("coluna.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);

    // ---- CSVs para o pos-processamento em Python (plot_liakopoulos.py) ----
    std::ofstream hist("historico.csv"), prof("perfis.csv");
    hist.precision(10);
    prof.precision(10);
    hist << "caso,passo,t_s,dt_s,newton,qb_m_s,w_topo_m,pw_meio_kPa,pw_topo_kPa,pontos_plasticos\n";
    prof << "caso,t_s,x_m,y_m,pw_kPa,uy_m,Sw\n";

    // ---- casos (Secao 10.7, itens 6 e 7) ----
    bool ok = true;
    for (const Case &c : casos)
        if (escolhidos.empty() || escolhidos.count(c.nome))
            ok = RunCase(gmesh, c.nome, c.saturada, c.steps, c.vtk, test, hist, prof) && ok;

    delete gmesh;
    return ok ? 0 : 1;
}

bool RunCase(TPZGeoMesh *gmesh, const std::string &caso, bool saturada, const TimeSteps &steps, bool vtk,
             const ColumnTest &test, std::ofstream &hist, std::ofstream &prof)
{
    const auto start = std::chrono::steady_clock::now();
    const UnsatParams p = LiakopoulosParams(saturada);

    // ---- malha computacional: u P2, p_w e p_a P1 (Taylor-Hood) ----
    const int porder = 2;
    TPZMultiphysicsCompMesh *cmesh = CreateColumnCompMesh(gmesh, porder, p);

    // ---- matriz estrutural, solver e analise ----
    TPZUnsatAnalysis *an = CreateAnalysis(cmesh);
    const Probes pr = FindProbes(*an, test);
    std::cout << caso << ": " << an->NEq() << " equacoes, " << pr.line.size() << " nos na vertical x = "
              << test.xProfile << " m" << std::endl;

    // ---- estado inicial (Secao 10.7, item 5) ----
    InitialState(*an);

    // ---- contorno do transiente (Secao 10.7, item 4) ----
    // apoios (base engastada, laterais com rolete), ar passivo (p_a = 0 em todos os vertices) e base drenante
    // (p_w = 0); topo e laterais impermeaveis: fluxo nulo, condicao natural
    std::map<int64_t, REAL> fixed = an->DisplacementConstraints();
    for (auto &[node, eq] : an->AirEq()) fixed[eq] = 0.;
    for (auto &[node, len] : pr.base) fixed[an->WaterEq().at(node)] = 0.;

    // ---- transiente (Secao 10.7, item 6) ----
    auto isProfileTime = [&](REAL t) {
        for (REAL tp : test.tProfiles)
            if (std::fabs(t - tp) < 1.e-6) return true;
        return false;
    };
    // t = 0: fluxo permanente para baixo com gradiente unitario, vazao q_b = k_s
    REAL t = 0.;
    int step = 0, vtkStep = 0, newtonTotal = 0;
    int64_t plastic = 0;
    WriteHistory(hist, caso, step, t, 0., 0, p.ks, *an, pr);
    WriteProfile(prof, caso, t, *an, pr);
    if (vtk) PostProcessVTK(*an, "coluna_" + caso, vtkStep++);
    bool ok = true;
    while (t < test.tEnd - 1.e-9) {
        const REAL dt = steps.Dt(t, test.tEnd);
        const TPZUnsatAnalysis::StepResult res = an->Solve(dt, fixed, {}, false);
        if (!res.conv) {
            std::cout << caso << ": o passo " << step + 1 << " (t = " << t + dt << " s) nao convergiu" << std::endl;
            ok = false;
            break;
        }
        // vazao na base pelas reacoes: o residuo usa o estado do passo anterior guardado na memoria
        const REAL qb = BaseOutflow(*an, pr, dt, test.W);
        an->Accept();
        t += dt;
        step++;
        // iteracoes de Newton (solucoes do sistema linear): Solve conta tambem a verificacao final
        const int newton = res.iter - 1;
        newtonTotal += newton;
        plastic = std::max(plastic, an->Materials()[0]->NPlasticPoints());
        WriteHistory(hist, caso, step, t, dt, newton, qb, *an, pr);
        if (isProfileTime(t)) {
            WriteProfile(prof, caso, t, *an, pr);
            if (vtk) PostProcessVTK(*an, "coluna_" + caso, vtkStep++);
            const TPZFMatrix<STATE> &X = an->X();
            std::cout << caso << ": t = " << std::setw(4) << t << " s  q_b = " << std::setprecision(4) << qb * 1.e6
                      << "e-6 m/s  p_w(y = 0,5 m) = " << X.GetVal(an->WaterEq().at(pr.mid), 0)
                      << " kPa  p_w(topo) = " << X.GetVal(an->WaterEq().at(pr.top), 0)
                      << " kPa  w = " << an->NodalDisplacement(X, pr.top, 1) * 1.e3 << " mm" << std::setprecision(6)
                      << std::endl;
        }
    }
    hist.flush();
    prof.flush();
    const REAL secs = std::chrono::duration<REAL>(std::chrono::steady_clock::now() - start).count();
    std::cout << caso << ": " << step << " passos, " << REAL(newtonTotal) / std::max(step, 1)
              << " iteracoes de Newton por passo, " << plastic << " pontos plasticos"
              << (plastic == 0 ? " (esqueleto elastico)" : " (AVISO: o esqueleto plastificou)") << ", " << secs
              << " s" << std::endl;

    delete an; // a analise antes das malhas que ela usa
    DeleteCompMesh(cmesh);
    return ok;
}

// Analise incremental (TPZUnsatAnalysis): Newton com o Jacobiano completo, as equacoes prescritas eliminadas por um
// filtro e a memoria do estado convergido. Apoios: base engastada (u_x = u_y = 0) e laterais com rolete (u_x = 0).
// O construtor ja escolhe a matriz estrutural e o solver; aqui as duas escolhas ficam explicitas.
TPZUnsatAnalysis *CreateAnalysis(TPZMultiphysicsCompMesh *cmesh)
{
    auto *an = new TPZUnsatAnalysis(cmesh, {{bcBase, {true, true}}, {bcDireita, {true, false}},
                                            {bcEsquerda, {true, false}}});

    // matriz estrutural: o Jacobiano do modelo acoplado nao e simetrico (Secao 8.5 do resumo)
#ifdef PZ_USING_MKL
    TPZSpStructMatrix<STATE> matrix(cmesh); // esparsa, fatorada pelo Pardiso
#else
    TPZSkylineNSymStructMatrix<STATE> matrix(cmesh); // skyline nao simetrica
#endif
    matrix.SetNumThreads(0);
    an->SetStructuralMatrix(matrix);

    // solver direto: LU nao simetrica
    TPZStepSolver<STATE> solver;
    solver.SetDirect(ELU);
    an->SetSolver(solver);

    // Newton (Secao 8.1): converge com e_u = |R_u|/|f_b| < 1e-7 e max |dp| < 1e-6 kPa; no maximo 40 iteracoes
    // (25 com as pressoes congeladas)
    an->SetTolerances(1.e-7, 1.e-6, 40, 25);
    return an;
}

// Estado inicial (Secao 10.7, item 5): p_w = p_a = 0 em todo o dominio; peso em 4 incrementos com as pressoes
// congeladas (dt = 0, so o equilibrio); no final u = 0, mantendo a tensao efetiva (os deslocamentos do peso sao
// zerados, como no artigo)
void InitialState(TPZUnsatAnalysis &an)
{
    an.ResetMemory();
    an.X().Zero();
    an.Restore();
    const auto fixed = an.FrozenFixed(an.X()); // apoios = 0 e todas as pressoes em zero
    an.SetGravityScalesSuction(true);
    const int n = 4;
    for (int k = 1; k <= n; k++) {
        an.SetGravityFactor(REAL(k) / n);
        if (!an.Solve(0., fixed, {}, true).conv) {
            std::cout << "Estado inicial: o estado geostatico nao convergiu" << std::endl;
            DebugStop();
        }
        an.Accept();
    }
    an.ResetDisplacements(); // u = 0, mantendo a tensao efetiva
    an.Accept();             // memoria com u = 0 e p_w = 0
}

// Nos de observacao a partir das coordenadas: vertical x = xProfile, ordenada de baixo para cima; no do meio e no
// mais proximo de y = H/2; no do topo e o ultimo da vertical; base: vertices do contorno bcBase
Probes FindProbes(TPZUnsatAnalysis &an, const ColumnTest &test)
{
    Probes pr;
    for (auto &[node, eq] : an.WaterEq())
        if (std::fabs(an.Coord(node, 0) - test.xProfile) < 1.e-6) pr.line.push_back(node);
    if (pr.line.empty()) DebugStop();
    std::sort(pr.line.begin(), pr.line.end(),
              [&](int64_t a, int64_t b) { return an.Coord(a, 1) < an.Coord(b, 1); });
    REAL dmin = std::numeric_limits<REAL>::max();
    for (int64_t node : pr.line) {
        const REAL d = std::fabs(an.Coord(node, 1) - 0.5 * test.H);
        if (d < dmin) dmin = d, pr.mid = node;
    }
    pr.top = pr.line.back();
    pr.base = an.BoundaryNodes(bcBase);
    return pr;
}

// Vazao na base (Secao 10.4a): com R_w o residuo da equacao da agua (Eq. 38, sem o termo de vazao prescrita),
// Q_b = -(1/dt) sum_j R_w,j nos vertices da base e a vazao que sai por unidade de espessura; q_b = Q_b / W.
// FullResidual devolve R = -ef do estado carregado na malha (o passo convergido); deve ser chamada antes do Accept,
// porque o armazenamento do residuo usa o estado do passo anterior guardado na memoria
REAL BaseOutflow(TPZUnsatAnalysis &an, const Probes &pr, REAL dt, REAL W)
{
    const TPZFMatrix<STATE> R = an.FullResidual();
    REAL sum = 0.;
    for (auto &[node, len] : pr.base) sum += R.GetVal(an.WaterEq().at(node), 0);
    return -sum / dt / W;
}

void WriteHistory(std::ofstream &csv, const std::string &caso, int step, REAL t, REAL dt, int newton, REAL qb,
                  TPZUnsatAnalysis &an, const Probes &pr)
{
    const TPZFMatrix<STATE> &X = an.X();
    csv << caso << ',' << step << ',' << t << ',' << dt << ',' << newton << ',' << qb << ','
        << an.NodalDisplacement(X, pr.top, 1) << ',' << X.GetVal(an.WaterEq().at(pr.mid), 0) << ','
        << X.GetVal(an.WaterEq().at(pr.top), 0) << ',' << an.Materials()[0]->NPlasticPoints() << '\n';
}

void WriteProfile(std::ofstream &csv, const std::string &caso, REAL t, TPZUnsatAnalysis &an, const Probes &pr)
{
    const TPZFMatrix<STATE> &X = an.X();
    const TPZUnsaturatedHydraulics &hydro = an.Materials()[0]->Hydraulics();
    for (int64_t node : pr.line) {
        const REAL pw = X.GetVal(an.WaterEq().at(node), 0), pa = X.GetVal(an.AirEq().at(node), 0);
        REAL S, Gamma;
        hydro.Retention(pa - pw, S, Gamma);
        csv << caso << ',' << t << ',' << an.Coord(node, 0) << ',' << an.Coord(node, 1) << ',' << pw << ','
            << an.NodalDisplacement(X, node, 1) << ',' << S << '\n';
    }
}

// Pos-processamento VTK do estado carregado na malha (ultimo estado convergido): todas as variaveis projetadas
// por elemento com o TPZPostProcAnalysis (campo linear descontinuo), em <nome>_pontos.scal_vec.<step>.vtk:
// solucao (deslocamento, poropressoes, succao, saturacao, k_rw e pressao de Bishop) e memoria dos pontos de
// integracao (tensoes efetivas, tracao positiva, e tipo de resposta plastica)
void PostProcessVTK(TPZUnsatAnalysis &an, const std::string &nome, int step)
{
    TPZPostProcAnalysis post;
    post.SetCompMesh(an.Mesh());
    TPZManVector<int, 1> matids = {matSolo};
    TPZManVector<std::string, 12> scal = {"DisplacementY", "WaterPressure", "AirPressure", "Suction",
                                          "Saturation", "RelativePermeabilityWater", "BishopPressure", "StressXX",
                                          "StressYY", "StressZZ", "StressXY", "FailureType"};
    TPZManVector<std::string, 1> vec = {"Displacement"};
    // o pos-processamento projeta todas as variaveis: escalares e vetores
    TPZManVector<std::string, 13> vars(scal.size() + vec.size());
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

// Coluna com o TUnsat: base engastada, laterais com rolete (u_x = 0), topo livre.
// Espacos: u de ordem porder (2 = triangulo de 6 nos), p_w e p_a de ordem porder - 1 (1: P1); Taylor-Hood P2-P1,
// estavel na condicao inf-sup (Secao 8.6 do resumo).
// Apoios e pressoes prescritas impostos pela analise (TPZUnsatAnalysis, fixU e mapa das equacoes prescritas): as
// condicoes EDirichletU do TUnsat nao contribuem, so marcam os contornos. O topo e as laterais sao impermeaveis
// (fluxo nulo, condicao natural): o topo nao precisa de elementos de contorno.
TPZMultiphysicsCompMesh *CreateColumnCompMesh(TPZGeoMesh *gmesh, int porder, const UnsatParams &p)
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

// Coluna: 0,10 m x 1,00 m, 4 x 40 quadrados de 0,025 m divididos em dois triangulos (320 elementos)
TPZGeoMesh *ReadColumnMesh()
{
    return ReadGmshMesh(std::string(MESH_DIR) + "/coluna.msh", {{"Solo", matSolo}},
                        {{"Base", bcBase}, {"Direita", bcDireita}, {"Topo", bcTopo}, {"Esquerda", bcEsquerda}});
}
