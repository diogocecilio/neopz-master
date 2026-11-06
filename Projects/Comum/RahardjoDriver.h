#pragma once
// Funcoes compartilhadas pelos projetos de verificacao (antes o main.cpp do pacote cluster de LearningMohrCoulomb2):
// parametros do solo nao saturado (UnsatParams), leitura de malha Gmsh, material TPZMatUnsaturatedUPwPa com Mohr-Coulomb
// associado, malhas computacionais, FS por reducao de resistencia (FactorSRM), aumento de gravidade (FactorGIM), carga
// limite (FactorLoad), pos-processamento VTK e o rebaixamento de reservatorio (SolveDraDown). Inclua num unico arquivo
// .cpp por executavel (as funcoes nao sao inline).
#ifndef MESH_DIR
#define MESH_DIR "."
#endif
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

int SolveDraDown()
{
    const int porder = 2;
    const REAL dia = 86400., H0 = 40., Hf = 30.; // niveis externos do reservatorio
    const REAL xl = 0., xr = 75.; // extremos horizontais de talude.geo (m)
    const REAL hl = H0, hr = H0; // cotas iniciais do lencol em xl e xr; editar aqui (m)
    // hl/hr sao cotas absolutas: hl=10 e hr=5 colocariam o lencol abaixo do pe (y=30).
    // Mantemos hl=hr=H0 para comparar o minimo GIM do mesmo caso 40 -> 30 m.
    const REAL Td = 5.*dia, Tf = 30.*dia, Dt = 0.1*dia;
    const int cadaFS = 10; // fora do intervalo detalhado: FS, CSV e VTK a cada 10 passos
    const REAL inicioDetalhe = 3.*dia, fimDetalhe = 5.*dia;
    const REAL tolTempo = 1.e-8*dia; // inclui as extremidades apesar do arredondamento
    // Dentro do intervalo detalhado, calcula FS e grava VTK em CADA passo (Dt=0.1 dia).
    if (!(xr > xl && cadaFS > 0 && Dt > 0. && Td > 0. && Tf >= Td &&
        0. <= inicioDetalhe && inicioDetalhe <= fimDetalhe &&
        30. <= Hf && Hf < H0 && H0 <= 40.)) {
        std::cerr << "Verifique xl, xr, cadaFS, Dt, Td, Tf, intervalo detalhado e H0/Hf.\n";
    return 1;
        }

        // ---- mesmas malhas, parametros e apoios do seu codigo ----
        TPZGeoMesh *gmeshslope = ReadSlopeMesh();
        std::ofstream out2("talude.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmeshslope, out2, true);
        const UnsatParams p = SlopeUnsatParams();
        TPZMultiphysicsCompMesh *cslope = CreateSlopeCompMesh(gmeshslope, porder, p);
        auto *anslope = new TPZUnsatAnalysis(cslope, {{bcBase, {true, true}}, {bcDireita, {true, false}},
                                             {bcEsquerda, {true, false}}});
        auto *reservatorio = dynamic_cast<TUnsat *>(cslope->FindMaterial(matSolo));
        auto *bcReservatorio = dynamic_cast<TPZBndCondT<STATE> *>(cslope->FindMaterial(bcSuperficie));
        if (!reservatorio || !bcReservatorio) DebugStop();
        REAL nivel = H0; // cota externa z_w(t), atualizada a cada passo do transiente
        // Neumann mecanico: p(y,t) = gamma_w max(nivel-y,0), tracao = -p n.
        // ENeumannU deve avaliar bc.ForcingFunctionBC() em ContributeBC.
        bcReservatorio->SetForcingFunctionBC(CargaReservatorio(nivel, reservatorio));
        const auto borda = anslope->BoundaryNodes(bcSuperficie);
        if (borda.empty()) DebugStop();

        // ---- estado inicial: reta entre (xl,hl) e (xr,hr); pa = 0 ----
        // hl/hr sao apenas condicoes iniciais; nao sao mantidos como CC nas laterais.
        // Para um estado hidraulico inicialmente equilibrado, resolver antes o fluxo
        // estacionario com as CC desejadas. Na superficie submersa, pw deve ser compativel
        // com gamma_w*(H0-y); um perfil arbitrario pode exigir uma etapa de acomodacao.
        TPZFMatrix<STATE> xslope(cslope->NEquations(), 1, 0.);
    for (const auto &eq : anslope->WaterEq()) {
        const REAL x = anslope->Coord(eq.first, 0), y = anslope->Coord(eq.first, 1);
        const REAL hf = hl + (hr-hl)*(x-xl)/(xr-xl);
        xslope(eq.second, 0) = p.gw*(hf-y); // negativo acima do lencol: preserva a succao
    }
    anslope->X() = xslope; // campo inicial prescrito; ainda nao e uma solucao estacionaria
    anslope->Restore();
    const auto inicial = anslope->FrozenFixed(xslope);
    anslope->SetStrengthReductionFactor(1.);
    anslope->SetGravityScalesSuction(true);
    bool convergiu = true;
    for (int k = 1; k <= 4; k++) {
        anslope->SetGravityFactor(REAL(k)/4.); // peso, poropressao e carga externa juntos
        if (!anslope->Solve(0., inicial, {}, true).conv) { convergiu = false; break; }
        anslope->Accept();
    }
    if (convergiu) anslope->ResetDisplacements(); // mantem as tensoes do estado inicial
    else std::cerr << "Estado inicial nao convergiu com a resistencia real.\n";

    // ---- transiente: altera o nivel e as condicoes de contorno ----
    std::ofstream csv("drawdown.csv");
    csv << "passo,tempo_s,tempo_dias,nivel_m,FS_SRM,FS_GIM\n" << std::setprecision(12);
    if (!csv) { convergiu = false; std::cerr << "Nao foi possivel criar drawdown.csv.\n"; }
    std::set<int64_t> saida; // parte exposta onde ha saida de agua e p_w = 0
    const REAL tolP = 1.e-6, tolQ = 1.e-10; // kPa e m/s, tolerancias do contorno
    const int nsteps = int(std::ceil(Tf/Dt));
    REAL t = 0.;
    for (int step = 0; step <= nsteps && convergiu; step++) {
        if (step > 0) {
            const REAL tn = std::min(step*Dt, Tf), dt = tn-t;
            nivel = H0 + (Hf-H0)*std::min(tn/Td, REAL(1));
            auto tentativa = saida;
            for (const auto &b : borda) {
                const int64_t node = b.first;
                if (anslope->Coord(node, 1) <= nivel) tentativa.erase(node);
                else if (anslope->X().GetVal(anslope->WaterEq().at(node), 0) > tolP) tentativa.insert(node);
            }
            convergiu = false;
            // Atualiza a face de saida sem aceitar memoria entre as tentativas.
            for (int it = 0; it < 40; it++) {
                auto fixed = anslope->DisplacementConstraints();
                for (const auto &b : borda) {
                    const int64_t node = b.first;
                    const REAL y = anslope->Coord(node, 1);
                    if (y <= nivel)
                        fixed[anslope->WaterEq().at(node)] = p.gw*(nivel-y);
                    else {
                        fixed[anslope->AirEq().at(node)] = 0.;
                        if (tentativa.count(node)) fixed[anslope->WaterEq().at(node)] = 0.;
                    }
                }
                if (!anslope->Solve(dt, fixed, {}, false).conv) break;
                const auto R = anslope->FullResidual(); // R_w/dt e positivo para entrada de agua
                const TPZFMatrix<STATE> solucao = cslope->Solution();
                auto nova = tentativa;
                for (const auto &b : borda) {
                    const int64_t node = b.first, eq = anslope->WaterEq().at(node);
                    if (anslope->Coord(node, 1) <= nivel) continue;
                    if (tentativa.count(node)) {
                        if (R.GetVal(eq, 0)/(dt*b.second) > tolQ) nova.erase(node);
                    } else if (solucao.GetVal(eq, 0) > tolP) nova.insert(node);
                }
                if (nova == tentativa) { convergiu = true; break; }
                tentativa = nova;
            }
            if (!convergiu) {
                anslope->Restore();
                nivel = H0 + (Hf-H0)*std::min(t/Td, REAL(1));
                std::cerr << "Passo " << step << " nao convergiu. Ultimo tempo aceito: " << t/dia
                << " dias. Reduza Dt e confira Newton/contorno; nao e prova isolada de ruptura.\n";
                break;
            }
            anslope->Accept();
            saida = tentativa;
            t = tn;
            std::cout << "Passo " << step << ": t = " << t/dia << " dias; H = " << nivel << " m\n";
        }

        const bool detalhar = t >= inicioDetalhe-tolTempo && t <= fimDetalhe+tolTempo;
        if (step % cadaFS == 0 || detalhar || step == nsteps) {
            // As rotinas de FS reconstroem o equilibrio com p congelado; a busca nao foi alterada.
            // Salva/restaura o estado para que ResetMemory() nao apague a historia do transiente.
            const auto estado = anslope->X();
            const auto memoria = anslope->SaveMemory();
            const REAL dtAnterior = reservatorio->TimeStep();
            const REAL FSslope = FactorSRM(*anslope, estado);
            if (detalhar && FSslope > 0.) {
                // Ultimo estado convergido da busca SRM, antes de o GIM apagar sua memoria.
                // A rotina retorna com F=1; repoe Fok somente para exportar esse estado.
                anslope->SetStrengthReductionFactor(FSslope);
                PostProcessVTK(*anslope, "talude_srm", step);
                anslope->SetStrengthReductionFactor(1.);
            }
            const REAL FGslope = FactorGIM(*anslope, estado);
            if (detalhar && FGslope > 0.) {
                // Ultimo estado convergido da busca GIM, com lambda correspondente a FGslope.
                // Estes dois VTKs mostram os ensaios de FS; nao sao estados do transiente.
                anslope->SetGravityScalesSuction(false);
                anslope->SetGravityFactor(FGslope);
                PostProcessVTK(*anslope, "talude_gim", step);
            }
            anslope->SetStrengthReductionFactor(1.);
            anslope->SetGravityFactor(1.);
            anslope->SetGravityScalesSuction(true);
            anslope->RestoreMemory(memoria);
            anslope->X() = estado;
            anslope->Restore();
            reservatorio->SetTimeStep(dtAnterior);
            reservatorio->SetFrozenPressures(false);
            csv << step << ',' << t << ',' << t/dia << ',' << nivel << ','
            << FSslope << ',' << FGslope << '\n';
            csv.flush();
            std::cout << "  FS (SRM) = " << FSslope << "; FS (GIM) = " << FGslope << '\n';
            // Limites das listas ORIGINAIS: ao atingir 4 ou 10, o valor e um limite inferior.
            if (FSslope >= 4. || FGslope >= 10.)
                std::cout << "  Limite da lista de fatores atingido; ampliar Fs/lams para localizar a ruptura.\n";
            // Estado fisico do transiente, sempre APOS restaurar solucao, memoria e parametros.
            // Em Dt=0.1 dia: 49 registros no CSV/VTK, sendo 21 entre 3 e 5 dias (inclusive).
            // Displacement e medido desde a inicializacao. Os invariantes de plasticidade do
            // transiente incluem a mudanca de referencia de ResetDisplacements(); nao os
            // interprete isoladamente como plasticidade acumulada somente pelo rebaixamento.
            PostProcessVTK(*anslope, "talude_drawdown", step);
        }
    }
    delete anslope;
    DeleteCompMesh(cslope);
    delete gmeshslope;
    return convergiu ? 0 : 1;
}

// o programa original (talude de Rahardjo e sapata); renomeado para o arquivo poder ser incluido pelos projetos
int RahardjoMain() {
    //SolveDraDown();
    const int porder = 2;
    // Ensaio de chuva inspirado em Rahardjo et al. (2007), Fig. 2: q = Ir na superficie.
    // Mantem sua malha, material e rotinas de FS; nao reproduz todos os parametros do artigo.
    const REAL hora = 3600., dia = 24.*hora;
    const REAL I_mm_h = 36., Ir = I_mm_h/(1000.*hora); // mm/h -> m/s; use I_mm_h=0 para controle
    const REAL Tchuva = 24.*hora, Tf = 10.*dia, Dt = 0.1*hora;
    const REAL xl = 0., xr = 75.; // extremos horizontais de talude.geo (m)
    const REAL Hw = 5., beta = 7.*M_PI/180., sucMax = 75.; // profundidade no pe, inclinacao, kPa
    const REAL hl = 30.-Hw+(45.-xl)*std::tan(beta), hr = 30.-Hw+(45.-xr)*std::tan(beta);
    // hl/hr sao cotas absolutas; pe da sua malha em (45,30). Pode editar hl e hr diretamente.
    const int cadaFS = 10; // FS, CSV e VTK a cada 10 passos: 1 hora com este Dt
    const REAL inicioDetalhe = 0., fimDetalhe = Tchuva;
    const REAL tolTempo = 1.e-8*hora;
    // No intervalo detalhado, salva tambem os mecanismos SRM/GIM, na mesma frequencia do FS.
    if (!(xr > xl && cadaFS > 0 && Dt > 0. && Tchuva > 0. && Tf >= Tchuva && I_mm_h >= 0. &&
        sucMax > 0. && 0. <= hl && hl < 40. && 0. <= hr && hr < 30. &&
        0. <= inicioDetalhe && inicioDetalhe <= fimDetalhe &&
        fimDetalhe <= Tf)) {
        std::cerr << "Verifique geometria, hl/hr, sucMax, chuva, tempos e cadaFS.\n";
    return 1;
        }

        // ---- mesmas malhas, parametros e apoios do seu codigo ----
        TPZGeoMesh *gmeshslope = ReadSlopeMesh();
        std::ofstream out2("talude.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmeshslope, out2, true);
        const UnsatParams p = SlopeUnsatParams();
        TPZMultiphysicsCompMesh *cslope = CreateSlopeCompMesh(gmeshslope, porder, p);
        auto *anslope = new TPZUnsatAnalysis(cslope, {{bcBase, {true, true}}, {bcDireita, {true, false}},
                                             {bcEsquerda, {true, false}}});
        auto *reservatorio = dynamic_cast<TUnsat *>(cslope->FindMaterial(matSolo));
        if (!reservatorio) DebugStop();
        // Sem reservatorio externo: ENeumannU fica com a tracao nula criada na malha.
        // A chuva entra no balanco de agua por fq, nao pela forcing function mecanica.
        const auto borda = anslope->BoundaryNodes(bcSuperficie);
    if (borda.empty()) DebugStop();
    // b.second e o comprimento tributario REAL: q=Ir normal ao contorno, como prescrito na Fig. 2.
    // Para chuva vertical medida por area horizontal, use BoundaryNodes(bcSuperficie, true).

    // ---- estado inicial: reta entre (xl,hl) e (xr,hr); pa = 0 ----
    // Limite de pw=-75 kPa do artigo. Este perfil e imposto; nao resolve o fluxo estacionario.
    // As cargas hidraulicas nas laterais abaixo do lencol serao mantidas no transiente.
    TPZFMatrix<STATE> xslope(cslope->NEquations(), 1, 0.);
    for (const auto &eq : anslope->WaterEq()) {
        const REAL x = anslope->Coord(eq.first, 0), y = anslope->Coord(eq.first, 1);
        const REAL hf = hl + (hr-hl)*(x-xl)/(xr-xl);
        xslope(eq.second, 0) = std::max(p.gw*(hf-y), -sucMax);
    }
    anslope->X() = xslope; // campo inicial prescrito; ainda nao e uma solucao estacionaria
    anslope->Restore();
    const auto inicial = anslope->FrozenFixed(xslope);
    anslope->SetStrengthReductionFactor(1.);
    anslope->SetGravityScalesSuction(true);
    bool convergiu = true;
    for (int k = 1; k <= 4; k++) {
        anslope->SetGravityFactor(REAL(k)/4.); // peso e poropressao juntos
        if (!anslope->Solve(0., inicial, {}, true).conv) { convergiu = false; break; }
        anslope->Accept();
    }
    if (convergiu) anslope->ResetDisplacements(); // mantem as tensoes do estado inicial
    else std::cerr << "Estado inicial nao convergiu com a resistencia real.\n";

    // ---- transiente: chuva por Tchuva e depois drenagem ----
    auto contorno = anslope->DisplacementConstraints();
    for (const int id : {bcEsquerda, bcDireita}) {
        const REAL h = id == bcEsquerda ? hl : hr;
        for (const auto &b : anslope->BoundaryNodes(id)) {
            const REAL y = anslope->Coord(b.first, 1);
            if (y <= h) contorno[anslope->WaterEq().at(b.first)] = p.gw*(h-y);
        }
    }
    // Base e laterais acima do lencol: fluxo natural nulo. Ar atmosferico na superficie.
    for (const auto &b : borda) contorno[anslope->AirEq().at(b.first)] = 0.;
    std::ofstream csv("chuva.csv");
    csv << "passo,tempo_s,tempo_h,tempo_dias,chuva_media_passo_mm_h,FS_SRM,FS_GIM\n" << std::setprecision(12);
    if (!csv) { convergiu = false; std::cerr << "Nao foi possivel criar chuva.csv.\n"; }
    std::set<int64_t> saida; // superficie saturada: pw=0; excesso escoa sem acumular lamina
    const REAL tolP = 1.e-6, tolQ = 1.e-10; // kPa e m/s, tolerancias do contorno
    const int nsteps = int(std::ceil(Tf/Dt));
    REAL t = 0., chuva = 0.; // chuva media no ultimo passo aceito; t=0 ainda nao tem passo
    for (int step = 0; step <= nsteps && convergiu; step++) {
        if (step > 0) {
            const REAL tn = std::min(step*Dt, Tf), dt = tn-t;
            // Media no passo: preserva o volume mesmo se Dt nao dividir Tchuva.
            const REAL chuvaPasso = Ir*std::max(REAL(0), std::min(tn,Tchuva)-t)/dt;
            auto tentativa = saida;
            for (const auto &b : borda) {
                const int64_t node = b.first;
                if (anslope->X().GetVal(anslope->WaterEq().at(node), 0) > tolP) tentativa.insert(node);
            }
            convergiu = false;
            // Alterna fluxo de chuva / pw=0 sem aceitar memoria entre as tentativas.
            for (int it = 0; it < 40; it++) {
                auto fixed = contorno;
                std::map<int64_t, REAL> fq;
                for (const auto &b : borda) {
                    const int64_t node = b.first, eq = anslope->WaterEq().at(node);
                    if (tentativa.count(node)) fixed[eq] = 0.;
                    else fq[eq] = chuvaPasso*b.second; // entrada positiva; Solve ja multiplica por dt
                }
                if (!anslope->Solve(dt, fixed, fq, false).conv) break;
                const auto R = anslope->FullResidual(); // R_w/dt e positivo para entrada de agua
                const TPZFMatrix<STATE> solucao = cslope->Solution();
                auto nova = tentativa;
                for (const auto &b : borda) {
                    const int64_t node = b.first, eq = anslope->WaterEq().at(node);
                    if (tentativa.count(node)) {
                        // pw=0 so e admissivel se a entrada exigida nao superar a chuva disponivel.
                        if (R.GetVal(eq, 0)/(dt*b.second) > chuvaPasso+tolQ) nova.erase(node);
                    } else if (solucao.GetVal(eq, 0) > tolP) nova.insert(node);
                }
                if (nova == tentativa) { convergiu = true; break; }
                tentativa = nova;
            }
            if (!convergiu) {
                anslope->Restore();
                std::cerr << "Passo " << step << " nao convergiu. Ultimo tempo aceito: " << t/dia
                << " dias. Reduza Dt e confira Newton/contorno; nao e prova isolada de ruptura.\n";
                break;
            }
            anslope->Accept();
            saida = tentativa;
            t = tn;
            chuva = chuvaPasso;
            std::cout << "Passo " << step << ": t = " << t/hora << " h; chuva media = " << chuva*1000.*hora << " mm/h\n";
        }

        const bool detalhar = t >= inicioDetalhe-tolTempo && t <= fimDetalhe+tolTempo;
        if (step % cadaFS == 0 || step == nsteps) {
            // As rotinas de FS reconstroem o equilibrio com p congelado; a busca nao foi alterada.
            // Salva/restaura o estado para que ResetMemory() nao apague a historia do transiente.
            const auto estado = anslope->X();
            const auto memoria = anslope->SaveMemory();
            const REAL dtAnterior = reservatorio->TimeStep();
            const REAL FSslope = FactorSRM(*anslope, estado);
            if (detalhar && FSslope > 0.) {
                // Ultimo estado convergido da busca SRM, antes de o GIM apagar sua memoria.
                // A rotina retorna com F=1; repoe Fok somente para exportar esse estado.
                anslope->SetStrengthReductionFactor(FSslope);
                PostProcessVTK(*anslope, "talude_chuva_srm", step);
                anslope->SetStrengthReductionFactor(1.);
            }
            const REAL FGslope = FactorGIM(*anslope, estado);
            if (detalhar && FGslope > 0.) {
                // Ultimo estado convergido da busca GIM, com lambda correspondente a FGslope.
                // Estes dois VTKs mostram os ensaios de FS; nao sao estados do transiente.
                anslope->SetGravityScalesSuction(false);
                anslope->SetGravityFactor(FGslope);
                PostProcessVTK(*anslope, "talude_chuva_gim", step);
            }
            anslope->SetStrengthReductionFactor(1.);
            anslope->SetGravityFactor(1.);
            anslope->SetGravityScalesSuction(true);
            anslope->RestoreMemory(memoria);
            anslope->X() = estado;
            anslope->Restore();
            reservatorio->SetTimeStep(dtAnterior);
            reservatorio->SetFrozenPressures(false);
            csv << step << ',' << t << ',' << t/hora << ',' << t/dia << ',' << chuva*1000.*hora << ','
            << FSslope << ',' << FGslope << '\n';
            csv.flush();
            std::cout << "  FS (SRM) = " << FSslope << "; FS (GIM) = " << FGslope << '\n';
            // Limites das listas ORIGINAIS: ao atingir 4 ou 10, o valor e um limite inferior.
            if (FSslope >= 4. || FGslope >= 10.)
                std::cout << "  Limite da lista de fatores atingido; ampliar Fs/lams para localizar a ruptura.\n";
            // Estado fisico do transiente, sempre APOS restaurar solucao, memoria e parametros.
            // Com Dt=0.1 h e cadaFS=10: registros horarios, incluindo t=0 e t=Tf.
            // Displacement e medido desde a inicializacao. Os invariantes de plasticidade do
            // transiente incluem a mudanca de referencia de ResetDisplacements(); nao os
            // interprete isoladamente como plasticidade acumulada somente pela chuva.
            PostProcessVTK(*anslope, "talude_chuva", step);
        }
    }
    delete anslope;
    DeleteCompMesh(cslope);
    delete gmeshslope;
    return convergiu ? 0 : 1;
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
// com a resistencia real (F = 1): lambda multiplica o peso e a parte positiva da pressao equivalente;
// a sucao fica com o valor de x. CargaReservatorio tambem usa lambda na tracao externa. Lambda cresce
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
    an.SetGravityScalesSuction(false); // peso e pressao equivalente positiva escalados; sucao preservada
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
