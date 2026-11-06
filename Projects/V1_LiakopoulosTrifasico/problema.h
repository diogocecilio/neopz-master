#pragma once
// V1_LiakopoulosTrifasico: o problema de verificacao. Coluna de areia de Liakopoulos (1965) drenando por gravidade,
// com ar passivo (u-p_w, p_a = 0: variante ii, "two phase model" de Laloui et al.) e com ar ativo (u-p_w-p_a: variante
// i, "three phase model"). Referencia: Laloui, Klubertanz & Vulliet (2003), Solid-liquid-air coupling in multiphase
// porous media, IJNAMG 27:183-206, doi:10.1002/nag.269 (Secao 5, Tabela I, Figs. 5, 6, 8 e 9).
//
// O caso de ar passivo e o porte do projeto LearningLiakopoulos (branch UnsaturatedSoil), validado contra Callari e
// Abati (2009) e contra a referencia 1D: mesma malha (malhas/coluna.msh), mesmos parametros, mesmo estado inicial,
// mesmos passos e mesma vazao pelas reacoes. O caso de ar ativo usa os mesmos parametros do solo e da agua e acrescenta
// o ar com os dados de Laloui et al. (Tabela I e Secoes 3.3 e 5): mu_a, rho_g, k_ra = (1 - S_e)^3 com S_r^r = 1% e
// k_ra >= 1%, e p_a = 0 no topo e na base (as duas extremidades abertas para a atmosfera; laterais impermeaveis).
// Unidades: kN, m, s, kPa. y para cima, base em y = 0.
#include "UnsatModel.h"
#include "UnsatAnalysis.h"

#include "TPZGmshReader.h"
#include "TPZVTKGeoMesh.h"
#include "pzfstrmatrix.h"
#include "pzgeoel.h"
#include "pzgmesh.h"
#include "pzpostprocanalysis.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#ifndef VERIF_DIR
#define VERIF_DIR "."
#endif

namespace v1 {

// ---------------------------------------------------------------------------------------------------------------
// Dados
// ---------------------------------------------------------------------------------------------------------------

/// ids de material: grupos fisicos da malha do Gmsh (malhas/coluna.geo) e da malha estruturada de UnsatModel.h
enum EIds {
    matSolo = 1,    // solo (dominio)
    bcBase = -1,    // base (y = 0): engastada, drenante (p_w = 0) e aberta ao ar (p_a = 0)
    bcDireita = -2, // lateral direita (x = 0,10 m): rolete, impermeavel a agua e ao ar
    bcTopo = -3,    // topo (y = 1,00 m): livre, impermeavel a agua; com ar ativo, p_a = 0 (atmosfera)
    bcEsquerda = -4 // lateral esquerda (x = 0): rolete, impermeavel a agua e ao ar
};

/// Geometria da coluna e vertical dos perfis
struct Coluna {
    REAL W = 0.10;        // largura (m): a coluna e 1D, a largura so entra na vazao por unidade de area
    REAL H = 1.00;        // altura (m)
    REAL xPerfil = 0.05;  // vertical dos perfis (m)
};

/**
 * Coluna de Liakopoulos com ar passivo (Callari e Abati 2009, Tabela 1 e Eq. 58; o projeto LearningLiakopoulos):
 * areia Del Monte. E = 1300 kPa e nu = 0,4 (os de Schrefler e Zhan 1993, os mesmos da Tabela I de Laloui et al.),
 * n = 0,2975, k_s = 9,81 x 4,5e-7 = 4,4145e-6 m/s (= k rho_w g / mu_w com k = 4,5e-13 m2 e mu_w = 1e-3 Pa s, os da
 * Tabela I), chi = S_w (Bishop, Eq. 39 de Laloui et al.), graos e agua quase incompressiveis (K_w = 2e6 kPa).
 * Leis hidraulicas de Liakopoulos: S_w = 1 - c_L p_c^b_L e k_rw = 1 - d_L (1 - S_w)^e_L.
 */
inline UnsatParams ParamsBifasico()
{
    UnsatParams p;
    p.E = 1300., p.nu = 0.4;               // E_sk (kPa) e nu_sk
    p.c = 1.e5, p.phi = 30. * M_PI / 180.; // Mohr-Coulomb com c' muito alto: o esqueleto fica elastico (linear)
    p.gamma = 1700. * 9.81e-3;             // gamma = 16,677 kN/m3 constante (so o estado geostatico, zerado depois)
    p.n = 0.2975;                          // porosidade
    p.Kw = 2.e6;                           // modulo volumetrico da agua (kPa)
    p.gw = 9.81;                           // peso especifico da agua (kN/m3)
    p.ks = 9.81 * 4.5e-7;                  // k_s = 4,4145e-6 m/s
    p.alpha = 1., p.Ks = 0.;               // Biot = 1 e graos incompressiveis
    p.bishop = TUnsat::EChiSaturation;     // chi = S_w
    // retencao de Liakopoulos: S_w = 1 - c_L p_c^b_L com p_c em kPa (no artigo de Callari, p_c em Pa e c_L = 1,9722e-11)
    p.hydro.fRetention = TPZUnsaturatedHydraulics::ELiakopoulos;
    p.hydro.fLB = 2.4279;
    p.hydro.fLC = 1.9722e-11 * std::pow(1000., p.hydro.fLB); // 3,7901e-4 kPa^-b_L
    p.hydro.fSwr = 0.;                                       // nao atua nesta faixa (S_w > 0,9)
    // permeabilidade relativa da agua: k_rw = 1 - d_L (1 - S_w)^e_L
    p.hydro.fKrw = TPZUnsaturatedHydraulics::ELiakopoulosW;
    p.hydro.fLD = 2.207, p.hydro.fLE = 1.0121;
    return p;
}

/**
 * Coluna com ar ativo. Solo e agua como em ParamsBifasico (a comparacao "duas fases x tres fases" isola o efeito do ar,
 * como na Secao 5.3 de Laloui et al.); ar com os dados de Laloui et al.:
 *  - mu_a = 1,8e-5 Pa s e mu_w = 1e-3 Pa s (Tabela I): mu_w / mu_a = 55,6;
 *  - rho_g = 1 kg/m3 (Tabela I): gamma_a = 9,81e-3 kN/m3 (so no termo de gravidade do fluxo de ar);
 *  - k_ra = (1 - S_e)^3 com S_e = (S_r - S_r^r) / (1 - S_r^r) e S_r^r = 1% (Eq. 61 e Secao 5.1): e a lei de
 *    Brooks-Corey da biblioteca, (1 - S_e)^2 (1 - S_e^((2 + lambda)/lambda)), no limite lambda -> infinito
 *    (lambda = 1e6: erro relativo < 1e-5);
 *  - k_ra >= 1% (Secao 5.2, "minimal relative air permeability to 1%");
 *  - ar gas ideal, compressibilidade 1/(p_a + p_atm) com p_atm = 101,325 kPa (a da biblioteca; ver o README sobre o
 *    M_g impresso na Eq. 28).
 * dadosLaloui = true troca ainda o que a Tabela I e a Secao 3.3 dao e a biblioteca tem: n = 0,297 e
 * k_rw = S_e^3 (Eq. 59, lei de potencia com p = 3 e S_r^r = 1%). A retencao continua a de Liakopoulos (a de Seker,
 * Eq. 55, nao existe na biblioteca).
 */
inline UnsatParams ParamsTrifasico(bool dadosLaloui = false)
{
    UnsatParams p = ParamsBifasico();
    p.muRatio = 1.e-3 / 1.8e-5; // mu_w / mu_a
    p.ga = 1.0 * 9.81e-3;       // rho_g g (kN/m3)
    p.patm = 101.325;           // pressao atmosferica (kPa): compressibilidade do ar ideal
    p.hydro.fKra = TPZUnsaturatedHydraulics::EBrooksCorey;
    p.hydro.fLambdaBC = 1.e6;   // (1 - S_e)^3
    p.hydro.fSwr = 0.01;        // S_r^r = 1%: entra so em S_e (k_ra e, com dadosLaloui, k_rw)
    p.hydro.fKraMin = 0.01;     // k_ra >= 1%
    if (dadosLaloui) {
        p.n = 0.297;
        p.hydro.fKrw = TPZUnsaturatedHydraulics::EPowerLaw;
        p.hydro.fPower = 3.;
    }
    return p;
}

/**
 * Sensibilidade (fora de Laloui et al.): o ar com a lei usual do teste de Liakopoulos com ar ativo na literatura do
 * benchmark (Schrefler e Zhan 1993; Schrefler e Scotta 2001): Brooks-Corey com lambda = 3, S_e = (S_w - 0,2)/0,8 e
 * k_ra >= 1e-4. O resto como em ParamsTrifasico(false). Serve para ver se a faixa da Fig. 6a se explica pela
 * permeabilidade efetiva do ar.
 */
inline UnsatParams ParamsTrifasicoAlert()
{
    UnsatParams p = ParamsTrifasico(false);
    p.hydro.fKra = TPZUnsaturatedHydraulics::EBrooksCorey;
    p.hydro.fLambdaBC = 3.;
    p.hydro.fSwr = 0.2;
    p.hydro.fKraMin = 1.e-4;
    return p;
}

/// Um caso: nome, parametros e se o ar e ativo (p_a incognita) ou passivo (p_a = 0 em todos os nos)
struct Caso {
    std::string nome;
    bool arAtivo = false;
    UnsatParams p;
};

/// Casos do projeto: "bifasico" (ar passivo, o porte validado), "trifasico" (ar ativo, mesma agua) e
/// "trifasico_laloui" (ar ativo com n e k_rw da Tabela I e da Eq. 59)
inline std::vector<Caso> Casos()
{
    return {{"bifasico", false, ParamsBifasico()},
            {"trifasico", true, ParamsTrifasico(false)},
            {"trifasico_laloui", true, ParamsTrifasico(true)},
            {"trifasico_alert", true, ParamsTrifasicoAlert()}};
}

inline Caso CasoPorNome(const std::string &nome)
{
    for (const Caso &c : Casos())
        if (c.nome == nome) return c;
    std::cerr << "caso desconhecido: " << nome << " (casos: bifasico, trifasico, trifasico_laloui, trifasico_alert)" << std::endl;
    std::exit(2);
}

// ---------------------------------------------------------------------------------------------------------------
// Passos de tempo e instantes dos perfis
// ---------------------------------------------------------------------------------------------------------------

/// Passos por fases: dt = fases[i].second enquanto t < fases[i].first (s), vezes fator
struct Passos {
    std::vector<std::pair<REAL, REAL>> fases;
    REAL fator = 1.;
    REAL Dt(REAL t) const
    {
        REAL dt = fases.back().second;
        for (auto &[tl, d] : fases)
            if (t < tl - 1.e-9) {
                dt = d;
                break;
            }
        return fator * dt;
    }
};

/// Passos da referencia 1D do LearningLiakopoulos (1 s ate 60 s, 5 s ate 600 s, 10 s ate 2000 s, 20 s ate 7500 s) e
/// 60 s depois (so o caso completo de 600 min passa de 7500 s)
inline Passos PassosReferencia(REAL fator = 1.)
{
    const REAL inf = std::numeric_limits<REAL>::max();
    return {{{60., 1.}, {600., 5.}, {2000., 10.}, {7500., 20.}, {inf, 60.}}, fator};
}

/// Instantes dos perfis (s): os das Figs. 5, 8 e 9 de Laloui et al. (5, 10, 20, 30, 60, 120 e 600 min) e os de
/// Callari e Abati usados na regressao (300, 600, 1200, 1800, 3600, 6400, 7200 e 7500 s). Os passos param neles.
inline std::vector<REAL> TemposPerfis()
{
    return {300., 600., 1200., 1800., 3600., 6400., 7200., 7500., 36000.};
}

// ---------------------------------------------------------------------------------------------------------------
// Malha
// ---------------------------------------------------------------------------------------------------------------

/// Le uma malha do Gmsh (.msh 4.1) com os grupos fisicos nomeados (dominio de dimensao 2, contornos de dimensao 1)
inline TPZGeoMesh *LeMalhaGmsh(const std::string &arquivo, const std::map<std::string, int> &dominios,
                               const std::map<std::string, int> &contornos)
{
    TPZManVector<std::map<std::string, int>, 4> nomes(4);
    for (auto &[nome, id] : dominios) nomes[2][nome] = id;
    for (auto &[nome, id] : contornos) nomes[1][nome] = id;
    TPZGmshReader leitor;
    leitor.SetDimNamePhysical(nomes);
    TPZGeoMesh *gmesh = leitor.GeometricGmshMesh(arquivo);
    gmesh->BuildConnectivity();
    return gmesh;
}

/**
 * Malha da coluna 0,10 x 1,00 m. ny = 0: a do LearningLiakopoulos (malhas/coluna.msh, 4 x 40 quadrados de 0,025 m
 * divididos em dois triangulos, 320 elementos; gerada de malhas/coluna.geo com gmsh -2). ny > 0: malha estruturada
 * 4 x ny de UnsatModel.h (estudo de malha), com os mesmos ids de contorno.
 */
inline TPZGeoMesh *MalhaColuna(int ny = 0, const Coluna &col = Coluna())
{
    if (ny > 0) return RectangleMesh(col.W, col.H, 4, ny);
    return LeMalhaGmsh(std::string(VERIF_DIR) + "/malhas/coluna.msh", {{"Solo", matSolo}},
                       {{"Base", bcBase}, {"Direita", bcDireita}, {"Topo", bcTopo}, {"Esquerda", bcEsquerda}});
}

// ---------------------------------------------------------------------------------------------------------------
// Resultados
// ---------------------------------------------------------------------------------------------------------------

/// Perfil na vertical x = xPerfil (y crescente)
struct Perfil {
    REAL t = 0.;
    std::vector<REAL> y, pw, pa, uy, S;
    /// valor interpolado linearmente em y
    REAL Em(const std::vector<REAL> &v, REAL yy) const
    {
        if (yy <= y.front()) return v.front();
        for (size_t i = 1; i < y.size(); ++i)
            if (yy <= y[i]) return v[i - 1] + (v[i] - v[i - 1]) * (yy - y[i - 1]) / (y[i] - y[i - 1]);
        return v.back();
    }
    /// minimo de v e a sua altura
    std::pair<REAL, REAL> Minimo(const std::vector<REAL> &v) const
    {
        size_t k = std::min_element(v.begin(), v.end()) - v.begin();
        return {v[k], y[k]};
    }
};

/// Uma linha do historico (um passo convergido)
struct Linha {
    int passo = 0, newton = 0, cortes = 0;
    REAL t = 0., dt = 0.;
    REAL qb = 0.;          // vazao de agua que sai pela base (m/s, por unidade de area)
    REAL qaTopo = 0.;      // vazao de ar que entra pelo topo (m/s, volume de ar a pressao do passo, por area)
    REAL qaBase = 0.;      // vazao de ar que entra pela base (m/s)
    REAL wTopo = 0.;       // deslocamento vertical do topo (m, negativo para baixo)
    REAL pwMeio = 0., pwTopo = 0.; // p_w em y = 0,5 m e no topo (kPa)
    REAL paMin = 0., yPaMin = 0.;  // minimo de p_a na vertical e a sua altura
};

struct Resultado {
    std::string caso;
    bool ok = true;
    int passos = 0, newtonTotal = 0, cortes = 0;
    int64_t neq = 0, pontosPlasticos = 0;
    double seg = 0.;
    std::vector<Linha> hist;
    std::vector<Perfil> perfis;
    /// perfil no instante t (s); nullptr se nao foi gravado
    const Perfil *PerfilEm(REAL t) const
    {
        for (const Perfil &p : perfis)
            if (std::fabs(p.t - t) < 1.e-6) return &p;
        return nullptr;
    }
    /// historico interpolado linearmente no tempo (campo dado por um ponteiro para membro de Linha)
    REAL HistEm(REAL Linha::*campo, REAL t) const
    {
        if (hist.empty()) return 0.;
        if (t <= hist.front().t) return hist.front().*campo;
        for (size_t i = 1; i < hist.size(); ++i)
            if (t <= hist[i].t) {
                const REAL a = hist[i - 1].*campo, b = hist[i].*campo;
                return a + (b - a) * (t - hist[i - 1].t) / (hist[i].t - hist[i - 1].t);
            }
        return hist.back().*campo;
    }
};

/// Opcoes de uma execucao
struct Opcoes {
    REAL tFim = 36000.;                       // fim do transiente (s): 600 min, como nas Figs. 5, 6, 8 e 9
    Passos passos = PassosReferencia();
    std::vector<REAL> tPerfis = TemposPerfis();
    int ny = 0;                               // 0: malha do Gmsh (4 x 40)
    bool vtk = false;                         // VTK nos instantes dos perfis
    int maxCortes = 6;                        // passo dividido por 2 ate 6 vezes se o Newton nao convergir
    bool verboso = true;
    Coluna col;
};

// ---------------------------------------------------------------------------------------------------------------
// Analise
// ---------------------------------------------------------------------------------------------------------------

/// Nos de observacao: vertical x = xPerfil (y crescente), no do meio e do topo; nos da base e do topo (contornos)
struct Sondas {
    std::vector<int64_t> linha;
    int64_t meio = -1, topo = -1;
    std::map<int64_t, REAL> base, topoBc;
};

inline Sondas AchaSondas(TPZUnsatAnalysis &an, const Coluna &col)
{
    Sondas s;
    for (auto &[node, eq] : an.WaterEq())
        if (std::fabs(an.Coord(node, 0) - col.xPerfil) < 1.e-6) s.linha.push_back(node);
    if (s.linha.empty()) DebugStop();
    std::sort(s.linha.begin(), s.linha.end(), [&](int64_t a, int64_t b) { return an.Coord(a, 1) < an.Coord(b, 1); });
    REAL dmin = std::numeric_limits<REAL>::max();
    for (int64_t node : s.linha) {
        const REAL d = std::fabs(an.Coord(node, 1) - 0.5 * col.H);
        if (d < dmin) dmin = d, s.meio = node;
    }
    s.topo = s.linha.back();
    s.base = an.BoundaryNodes(bcBase);
    s.topoBc = an.BoundaryNodes(bcTopo);
    return s;
}

/**
 * Estado inicial (como no LearningLiakopoulos e em Laloui et al., que seguem Gawin et al.: p_w uniforme): p_w = p_a = 0
 * em toda a coluna, que e o fluxo permanente para baixo com gradiente unitario (vazao k_s). Peso em 4 incrementos com
 * as pressoes congeladas e deslocamentos zerados mantendo a tensao efetiva.
 */
inline void EstadoInicial(TPZUnsatAnalysis &an)
{
    an.ResetMemory();
    an.X().Zero();
    an.Restore();
    const auto fixed = an.FrozenFixed(an.X());
    an.SetGravityScalesSuction(true);
    const int n = 4;
    for (int k = 1; k <= n; k++) {
        an.SetGravityFactor(REAL(k) / n);
        if (!an.Solve(0., fixed, {}, true).conv) {
            std::cerr << "Estado inicial: o estado geostatico nao convergiu" << std::endl;
            DebugStop();
        }
        an.Accept();
    }
    an.ResetDisplacements();
    an.Accept();
}

/// Soma das reacoes R/dt (residuo da equacao do espaco 1 = agua ou 2 = ar) nos nos dados, por unidade de area. Para a
/// agua na base, -soma e a vazao que sai; para o ar, -soma e a vazao que sai (o sinal e trocado por quem chama).
/// FullResidual usa o estado do passo anterior guardado na memoria: chamar antes do Accept.
inline REAL Reacao(TPZUnsatAnalysis &an, const TPZFMatrix<STATE> &R, const std::map<int64_t, REAL> &nos, int espaco,
                   REAL dt, REAL W)
{
    REAL soma = 0.;
    for (auto &[node, len] : nos) soma += R.GetVal((espaco == 1 ? an.WaterEq() : an.AirEq()).at(node), 0);
    return soma / dt / W;
}

/// Pos-processamento VTK do estado convergido (variaveis da solucao e dos pontos de integracao, por elemento)
inline void GravaVTK(TPZUnsatAnalysis &an, const std::string &nome, int passo)
{
    TPZPostProcAnalysis post;
    post.SetCompMesh(an.Mesh());
    TPZManVector<int, 1> matids = {matSolo};
    TPZManVector<std::string, 10> scal = {"DisplacementY", "WaterPressure", "AirPressure", "Suction",
                                          "Saturation", "RelativePermeabilityWater", "BishopPressure",
                                          "StressXX", "StressYY", "StressXY"};
    TPZManVector<std::string, 1> vec = {"Displacement"};
    TPZManVector<std::string, 11> vars(scal.size() + vec.size());
    for (int i = 0; i < scal.size(); i++) vars[i] = scal[i];
    for (int i = 0; i < vec.size(); i++) vars[scal.size() + i] = vec[i];
    post.SetPostProcessVariables(matids, vars);
    TPZFStructMatrix<STATE> sm(post.Mesh());
    sm.SetNumThreads(0);
    post.SetStructuralMatrix(sm);
    post.DefineGraphMesh(2, scal, vec, nome + ".vtk");
    post.TransferSolution();
    post.SetStep(passo);
    post.PostProcess(0);
}

/**
 * Roda um caso: malha, espacos (u P2, p_w e p_a P1, Taylor-Hood), estado inicial e o transiente ate tFim com Euler
 * implicito e Newton (TPZUnsatAnalysis). Contorno: base engastada, laterais com rolete (u_x = 0), topo livre; agua:
 * p_w = 0 na base, topo e laterais impermeaveis (condicao natural); ar passivo: p_a = 0 em todos os nos; ar ativo:
 * p_a = 0 na base e no topo, laterais impermeaveis. Um passo que nao converge e refeito com metade do dt (ate
 * maxCortes vezes); os passos param nos instantes dos perfis.
 */
inline Resultado RodaCaso(const Caso &caso, const Opcoes &op)
{
    const auto inicio = std::chrono::steady_clock::now();
    Resultado r;
    r.caso = caso.nome;
    TPZGeoMesh *gmesh = MalhaColuna(op.ny, op.col);
    TPZMultiphysicsCompMesh *cmesh = UnsatMesh(gmesh, {{matSolo, caso.p}});
    {
        TPZUnsatAnalysis an(cmesh, {{bcBase, {true, true}}, {bcDireita, {true, false}}, {bcEsquerda, {true, false}}});
        // Newton: |R_u|/|peso| < 1e-7 e max |dp| < 1e-6 kPa (as do LearningLiakopoulos)
        an.SetTolerances(1.e-7, 1.e-6, 40, 25);
        r.neq = an.NEq();
        const Sondas s = AchaSondas(an, op.col);
        if (op.verboso)
            std::cout << caso.nome << ": " << r.neq << " equacoes, " << s.linha.size() << " nos na vertical x = "
                      << op.col.xPerfil << " m, ar " << (caso.arAtivo ? "ativo" : "passivo") << std::endl;
        EstadoInicial(an);

        // condicoes de Dirichlet do transiente
        std::map<int64_t, REAL> fixed = an.DisplacementConstraints();
        for (auto &[node, len] : s.base) fixed[an.WaterEq().at(node)] = 0.;
        if (caso.arAtivo) {
            for (auto &[node, len] : s.base) fixed[an.AirEq().at(node)] = 0.;
            for (auto &[node, len] : s.topoBc) fixed[an.AirEq().at(node)] = 0.;
        } else {
            for (auto &[node, eq] : an.AirEq()) fixed[eq] = 0.;
        }

        const TPZUnsaturatedHydraulics &hydro = an.Materials()[0]->Hydraulics();
        auto perfil = [&](REAL t) {
            const TPZFMatrix<STATE> &X = an.X();
            Perfil p;
            p.t = t;
            for (int64_t node : s.linha) {
                const REAL pw = X.GetVal(an.WaterEq().at(node), 0), pa = X.GetVal(an.AirEq().at(node), 0);
                REAL S, G;
                hydro.Retention(pa - pw, S, G);
                p.y.push_back(an.Coord(node, 1));
                p.pw.push_back(pw);
                p.pa.push_back(pa);
                p.uy.push_back(an.NodalDisplacement(X, node, 1));
                p.S.push_back(S);
            }
            return p;
        };
        auto linha = [&](int passo, REAL t, REAL dt, int newton, int cortes, REAL qb, REAL qaT, REAL qaB) {
            const TPZFMatrix<STATE> &X = an.X();
            Linha l;
            l.passo = passo, l.t = t, l.dt = dt, l.newton = newton, l.cortes = cortes;
            l.qb = qb, l.qaTopo = qaT, l.qaBase = qaB;
            l.wTopo = an.NodalDisplacement(X, s.topo, 1);
            l.pwMeio = X.GetVal(an.WaterEq().at(s.meio), 0);
            l.pwTopo = X.GetVal(an.WaterEq().at(s.topo), 0);
            l.paMin = std::numeric_limits<REAL>::max();
            for (int64_t node : s.linha) {
                const REAL pa = X.GetVal(an.AirEq().at(node), 0);
                if (pa < l.paMin) l.paMin = pa, l.yPaMin = an.Coord(node, 1);
            }
            return l;
        };

        // t = 0: fluxo permanente para baixo com gradiente unitario, vazao q_b = k_s
        REAL t = 0.;
        int passo = 0, vtkPasso = 0;
        r.hist.push_back(linha(0, 0., 0., 0, 0, caso.p.ks, 0., 0.));
        r.perfis.push_back(perfil(0.));
        if (op.vtk) GravaVTK(an, "coluna_" + caso.nome, vtkPasso++);
        size_t proxPerfil = 0;
        while (proxPerfil < op.tPerfis.size() && op.tPerfis[proxPerfil] <= 1.e-9) proxPerfil++;
        while (t < op.tFim - 1.e-9) {
            REAL alvo = op.tFim; // o passo nao passa do proximo perfil nem do fim
            if (proxPerfil < op.tPerfis.size()) alvo = std::min(alvo, op.tPerfis[proxPerfil]);
            REAL dt = std::min(op.passos.Dt(t), alvo - t);
            if (alvo - t - dt < 1.e-6 * dt) dt = alvo - t;
            TPZUnsatAnalysis::StepResult sr;
            int cortes = 0;
            for (;;) {
                sr = an.Solve(dt, fixed, {}, false);
                if (sr.conv || cortes >= op.maxCortes) break;
                an.Restore();
                dt *= 0.5;
                cortes++;
            }
            if (!sr.conv) {
                std::cerr << caso.nome << ": o passo " << passo + 1 << " (t = " << t + dt << " s) nao convergiu"
                          << std::endl;
                r.ok = false;
                break;
            }
            // reacoes (antes do Accept): vazao de agua que sai pela base e de ar que entra pelo topo e pela base
            const TPZFMatrix<STATE> R = an.FullResidual();
            const REAL qb = -Reacao(an, R, s.base, 1, dt, op.col.W);
            const REAL qaT = caso.arAtivo ? Reacao(an, R, s.topoBc, 2, dt, op.col.W) : 0.;
            const REAL qaB = caso.arAtivo ? Reacao(an, R, s.base, 2, dt, op.col.W) : 0.;
            an.Accept();
            t += dt;
            passo++;
            const int newton = sr.iter - 1; // Solve conta tambem a verificacao final
            r.newtonTotal += newton;
            r.cortes += cortes;
            r.pontosPlasticos = std::max(r.pontosPlasticos, an.Materials()[0]->NPlasticPoints());
            r.hist.push_back(linha(passo, t, dt, newton, cortes, qb, qaT, qaB));
            if (proxPerfil < op.tPerfis.size() && std::fabs(t - op.tPerfis[proxPerfil]) < 1.e-6) {
                r.perfis.push_back(perfil(t));
                proxPerfil++;
                if (op.vtk) GravaVTK(an, "coluna_" + caso.nome, vtkPasso++);
                if (op.verboso) {
                    const Linha &l = r.hist.back();
                    std::cout << caso.nome << ": t = " << std::setw(6) << t / 60. << " min  q_b = "
                              << std::setprecision(4) << qb * 1.e6 << "e-6 m/s  p_w(0,5 m) = " << l.pwMeio
                              << " kPa  p_w(topo) = " << l.pwTopo << " kPa  w = " << l.wTopo * 1.e3 << " mm";
                    if (caso.arAtivo) std::cout << "  min p_a = " << l.paMin << " kPa em y = " << l.yPaMin << " m";
                    std::cout << std::setprecision(6) << std::endl;
                }
            }
        }
        r.passos = passo;
    }
    // a analise ja foi destruida: apaga as malhas
    {
        TPZVec<TPZCompMesh *> mv = cmesh->MeshVector();
        delete cmesh;
        for (TPZCompMesh *m : mv) delete m;
        delete gmesh;
    }
    r.seg = std::chrono::duration<double>(std::chrono::steady_clock::now() - inicio).count();
    if (op.verboso)
        std::cout << caso.nome << ": " << r.passos << " passos (" << r.cortes << " cortes de passo), "
                  << REAL(r.newtonTotal) / std::max(r.passos, 1) << " iteracoes de Newton por passo, "
                  << r.pontosPlasticos << " pontos plasticos, " << r.seg << " s" << std::endl;
    return r;
}

// ---------------------------------------------------------------------------------------------------------------
// CSVs
// ---------------------------------------------------------------------------------------------------------------

inline void CabecalhoHistorico(std::ostream &o)
{
    o << "caso,passo,t_s,dt_s,newton,cortes,qb_m_s,qa_topo_m_s,qa_base_m_s,w_topo_m,pw_meio_kPa,pw_topo_kPa,"
         "pa_min_kPa,y_pa_min_m\n";
}

inline void GravaHistorico(std::ostream &o, const Resultado &r)
{
    o << std::setprecision(10);
    for (const Linha &l : r.hist)
        o << r.caso << ',' << l.passo << ',' << l.t << ',' << l.dt << ',' << l.newton << ',' << l.cortes << ','
          << l.qb << ',' << l.qaTopo << ',' << l.qaBase << ',' << l.wTopo << ',' << l.pwMeio << ',' << l.pwTopo
          << ',' << l.paMin << ',' << l.yPaMin << '\n';
}

inline void CabecalhoPerfis(std::ostream &o) { o << "caso,t_s,t_min,y_m,pw_kPa,pa_kPa,succao_kPa,Sw,uy_m\n"; }

inline void GravaPerfis(std::ostream &o, const Resultado &r)
{
    o << std::setprecision(10);
    for (const Perfil &p : r.perfis)
        for (size_t i = 0; i < p.y.size(); ++i)
            o << r.caso << ',' << p.t << ',' << p.t / 60. << ',' << p.y[i] << ',' << p.pw[i] << ',' << p.pa[i] << ','
              << p.pa[i] - p.pw[i] << ',' << p.S[i] << ',' << p.uy[i] << '\n';
}

// ---------------------------------------------------------------------------------------------------------------
// Dados de referencia (dados/*.csv): figura,curva,variavel_x,x,variavel_y,y; linhas com # sao comentarios
// ---------------------------------------------------------------------------------------------------------------

struct Curva {
    std::string figura, curva;
    std::vector<REAL> x, y;
};

inline std::vector<Curva> LeCurvas(const std::string &arquivo)
{
    std::ifstream in(arquivo);
    if (!in) {
        std::cerr << "nao abriu " << arquivo << std::endl;
        return {};
    }
    std::vector<Curva> cs;
    std::string linha;
    bool cabecalho = true;
    while (std::getline(in, linha)) {
        if (linha.empty() || linha[0] == '#') continue;
        if (cabecalho) {
            cabecalho = false;
            continue;
        }
        std::vector<std::string> campos;
        std::stringstream ss(linha);
        std::string c;
        while (std::getline(ss, c, ',')) campos.push_back(c);
        if (campos.size() < 6) continue;
        if (cs.empty() || cs.back().figura != campos[0] || cs.back().curva != campos[1])
            cs.push_back({campos[0], campos[1], {}, {}});
        cs.back().x.push_back(std::stod(campos[3]));
        cs.back().y.push_back(std::stod(campos[5]));
    }
    return cs;
}

inline const Curva *AchaCurva(const std::vector<Curva> &cs, const std::string &figura, const std::string &curva)
{
    for (const Curva &c : cs)
        if (c.figura == figura && c.curva == curva) return &c;
    return nullptr;
}

/// Caminho de um arquivo de dados do projeto
inline std::string Dados(const std::string &nome) { return std::string(VERIF_DIR) + "/dados/" + nome; }

} // namespace v1
