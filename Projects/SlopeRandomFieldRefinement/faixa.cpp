// SRM numa faixa longa de talude infinito (passo 4 do plano de reproducao): o fator de seguranca por reducao de
// resistencia (o SRMArcLength do UnsatAnalysis, o mesmo do SlopeRainMonteCarlo, com as pressoes congeladas) de uma camada
// inclinada sobre base rigida, contra o FS fechado de talude infinito de Lu e Godt (o do SlopeRainColumn):
//   FS(z) = [c' + (W cos^2 beta - sigma_s) tan phi'] / (W sin beta cos beta),   FS_LEM = min_z FS(z),
//   W = integral de gamma(z) dz desde a superficie (z vertical), sigma_s = S_e p_w se p_w < 0, senao p_w.
// Com Mohr-Coulomb associado o SRM de um talude infinito e o minimo em z do equilibrio limite (o mecanismo e o bloco que
// desliza num plano paralelo a superficie).
// A faixa, nos eixos do talude: x ao longo da superficie (para baixo), y normal a ela; [0, L] x [0, T], T = H cos beta a
// espessura normal de uma camada de profundidade vertical H. A gravidade e (sin beta, -cos beta) (SetGravityDirection).
// Um ponto na profundidade normal d = T - y esta na profundidade vertical z = d / cos beta: e o mapa da coluna vertical do
// SlopeRainColumn (o perfil p_w(z) dela, em z vertical) para a faixa. Base engastada. As laterais:
//   roletes (u_x = 0 nas pontas): as pontas presas seguram a camada, e o FS da faixa fica acima do de talude infinito,
//     com um erro que cai com H/L. O caminho do SRM tem um joelho: perto do FS de talude infinito o interior comeca a
//     deslizar, e dai F ainda sobe devagar ate as pontas romperem. Com a tolerancia padrao do patamar (0,001) o caminho
//     de uma faixa longa pode parar no joelho, com F ainda subindo: use tol_patamar=1e-5 e veja fim e dF_ult na CSV;
//   periodicas (u(L, y) = u(0, y)): sem pontas, o talude infinito de verdade; uma faixa curta (nx=4) basta, e o FS nao
//     depende de L. Os connects dos lados em x = L passam a ser os de x = 0 (Periodiza).
// O material e o do SlopeRainMonteCarlo/SlopeRainColumn: TPZMatUnsaturatedUPwPa (Mohr-Coulomb PV2 associado, chi = S_e de
// van Genuchten, p_a = 0); o peso e constante (gamma) ou o da mistura gamma_d + theta gw, theta = theta_s S_w(p_w), que
// e o W do Lu-Godt do SlopeRainColumn (SetUnitWeightMixture com gamma_s = gamma_d / (1 - theta_s), porosidade theta_s).
// O FS fechado e calculado aqui com os mesmos S_w e S_e (TPZUnsaturatedHydraulics) e W por trapezios numa grade de 4000
// intervalos, com duas p_w:
//   FS_LEM  a p_w(z) dada (no perfil, a linear entre os pontos do arquivo): o alvo do talude infinito;
//   FS_P1   a p_w que a malha recebe: SetFrozenPressures da p_w(z) nos nos, e o material interpola linearmente entre
//           eles. Os nos da malha de partida estao em z_k = k H / ny, e FS_P1 usa p_w(z_k) linear entre eles. E o alvo
//           da propria malha (o resto da diferenca e o erro do MEF); com nref > 0, e o da malha de partida.
//   Com pressao=perfil, a tela tambem da o FS_min do proprio arquivo.
// Uso: SlopeStripSRM [chave=valor ...], com (padrao entre parenteses: os numeros do Cho 2014):
//   H (2)               profundidade vertical da camada (m); a espessura normal e H cos beta
//   LH (20)             L / H (L ao longo da superficie)
//   ny (8)              celulas na espessura; nx = round(L / (T / ny)) (celulas quadradas), 2 triangulos por celula
//                       com a diagonal alternada
//   nx (0)              > 0: nx celulas ao longo da faixa, quadradas (L = nx T / ny, e LH sai dai); 0: de LH
//   laterais (roletes)  roletes ou periodicas (nx >= 2, nref 0)
//   nref (0)            niveis de refinamento adaptativo em h (os elementos com sqrt(J2(eps_p)) >= 0,1 do maior,
//                       como o FatorAdaptativo do montecarlo.cpp)
//   beta (40), c (5), phi (32)   graus, kPa, graus (dilatancia = phi: associado)
//   E (1e4), nu (0.3)   kPa
//   peso (constante)    constante: gamma; mistura: gamma_d + theta_s S_w gw
//   gamma (20), gamma_d (16), gw (9.81)   kN/m3
//   vg (0.358,3.58e-4,0.149,1.289)       theta_s,theta_r,alfa (1/kPa),n: S_w e S_e (chi = S_e)
//   pressao (seco)      seco: p_w = 0;
//                       lencol: percolacao paralela ao talude com o nivel na profundidade vertical dw:
//                         p_w = gw (z - dw) cos^2 beta abaixo; acima, 0 (acima=zero) ou o mesmo, negativo (acima=hidro);
//                       perfil: p_w(z) do <saida>_perfil_<n>.csv do SlopeRainColumn (arquivo=, t=), linear em z
//   dw (1), acima (zero), arquivo (), t (0)
//   srm_f0 (0.8)        F do peso no SRM (SetSRMStart)
//   tol_patamar (0.001) fim do caminho do SRM por patamar: um passo que sobe F menos que isso (com o arco pedido)
//                       (SetPlateauTolerance; 0,001 e o padrao dos outros drivers)
//   passos (100)        passos do caminho do SRM (SetMaxPathSteps)
//   verbose (0), vtk () prefixo do VTK do estado final (u, p_w, sqrt(J2(eps_p)))
//   csv (faixa.csv)     acrescenta uma linha (para se o arquivo tiver outro cabecalho): caso,laterais,LH,nx,ny,nref,
//                       elementos,equacoes,FS_SRM,FS_LEM,z_LEM,FS_P1,z_P1,z_banda,z_banda_min,z_banda_max,erro_rel,
//                       erro_P1,fim,passos,dF_ult,dF_ant,crit_ult,tol_patamar,segundos
//   caso ()             rotulo da linha
// z_banda: a profundidade vertical do centro do elemento de maior sqrt(J2(eps_p)) no terco central (x em [L/3, 2L/3];
// com laterais periodicas, a faixa toda); z_banda_min/max: a faixa dos centros com >= 0,5 do maior la.
// erro_rel = FS_SRM / FS_LEM - 1 e erro_P1 = FS_SRM / FS_P1 - 1. O fim do caminho do SRM (PathEnd/PathHistory): fim (o
// motivo: patamar, 100 passos, ...), passos (os aceitos), dF_ult e dF_ant (o aumento de F no ultimo passo aceito e no
// anterior) e crit_ult = dF_ult lPedido / l, o numero que o patamar compara com tol_patamar. Um patamar de verdade tem
// dF caindo para zero (dF_ult << dF_ant e << tol); um dF_ult pouco abaixo da tol com dF_ant parecido e um joelho.
#include "HPMaterial.h"
#include "UnsatAnalysis.h"
#include "UnsatModel.h"

#include "Plasticity/TPZElasticResponse.h"
#include "Plasticity/TPZMatUnsaturatedUPwPa.h"
#include "Plasticity/TPZPlasticStepVoigt.h"
#include "Plasticity/TPZYCMohrCoulombPV2.h"
#include "TPZGeoLinear.h"
#include "pzgeotriangle.h"
#include "pzintel.h" // TPZInterpolatedElement (Periodiza)
#include "pzmultiphysicselement.h"
#include "tpzgeoelrefpattern.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using TMC = TPZPlasticStepVoigt<TPZYCMohrCoulombPV2, TPZElasticResponse>;
using TUnsat = TPZMatUnsaturatedUPwPa<TMC, TPZUnsaturatedMem>;
using TUnsatHP = HPMaterial<TUnsat>;

enum EIds { matSolo = 1, base = -1, direita = -2, esquerda = -5 };

[[noreturn]] void Erro(const std::string &m)
{
    std::fprintf(stderr, "%s\n", m.c_str());
    std::exit(1);
}

REAL Numero(const std::string &k, const std::string &s)
{
    size_t pos = 0;
    REAL v = 0.;
    try {
        v = std::stod(s, &pos);
    } catch (const std::exception &) {
        Erro(k + "=" + s + ": valor invalido");
    }
    if (pos != s.size()) Erro(k + "=" + s + ": valor invalido");
    return v;
}

std::vector<REAL> Lista(const std::string &k, const std::string &s)
{
    std::vector<REAL> v;
    std::stringstream ss(s);
    std::string it;
    while (std::getline(ss, it, ',')) v.push_back(Numero(k, it));
    return v;
}

struct Opcoes {
    REAL H = 2., LH = 20., beta = 40., c = 5., phi = 32., E = 1.e4, nu = 0.3, gamma = 20., gammaD = 16., gw = 9.81;
    REAL ts = 0.358, tr = 3.58e-4, a = 0.149, nvg = 1.289;
    int ny = 8, nref = 0, nx = 0, passos = 100;
    bool mistura = false, acimaHidro = false, verbose = false, periodica = false;
    std::string pressao = "seco", arquivo, vtk, csv = "faixa.csv", caso;
    REAL dw = 1., t = 0., srmF0 = 0.8, tolPatamar = 0.001;
};

Opcoes LeOpcoes(int argc, char *argv[])
{
    Opcoes o;
    const std::set<std::string> chaves = {"H",     "LH",      "ny",    "nref",    "beta",   "c",   "phi",    "E",
                                          "nu",    "peso",    "gamma", "gamma_d", "gw",     "vg",  "pressao", "dw",
                                          "acima", "arquivo", "t",     "srm_f0",  "verbose", "vtk", "csv",    "caso",
                                          "nx",    "laterais", "tol_patamar", "passos"};
    for (int i = 1; i < argc; i++) {
        const std::string s = argv[i];
        const size_t eq = s.find('=');
        if (eq == std::string::npos) Erro("argumento sem '=': " + s);
        const std::string k = s.substr(0, eq), v = s.substr(eq + 1);
        if (!chaves.count(k)) Erro("opcao desconhecida: " + k + " (as opcoes estao no cabecalho do faixa.cpp)");
        if (k == "H") o.H = Numero(k, v);
        else if (k == "LH") o.LH = Numero(k, v);
        else if (k == "ny") o.ny = int(Numero(k, v));
        else if (k == "nref") o.nref = int(Numero(k, v));
        else if (k == "beta") o.beta = Numero(k, v);
        else if (k == "c") o.c = Numero(k, v);
        else if (k == "phi") o.phi = Numero(k, v);
        else if (k == "E") o.E = Numero(k, v);
        else if (k == "nu") o.nu = Numero(k, v);
        else if (k == "peso") {
            if (v != "constante" && v != "mistura") Erro("peso: constante ou mistura");
            o.mistura = v == "mistura";
        } else if (k == "gamma") o.gamma = Numero(k, v);
        else if (k == "gamma_d") o.gammaD = Numero(k, v);
        else if (k == "gw") o.gw = Numero(k, v);
        else if (k == "vg") {
            const std::vector<REAL> x = Lista(k, v);
            if (x.size() != 4) Erro("vg: theta_s,theta_r,alfa,n");
            o.ts = x[0], o.tr = x[1], o.a = x[2], o.nvg = x[3];
        } else if (k == "pressao") {
            if (v != "seco" && v != "lencol" && v != "perfil") Erro("pressao: seco, lencol ou perfil");
            o.pressao = v;
        } else if (k == "dw") o.dw = Numero(k, v);
        else if (k == "acima") {
            if (v != "zero" && v != "hidro") Erro("acima: zero ou hidro");
            o.acimaHidro = v == "hidro";
        } else if (k == "arquivo") o.arquivo = v;
        else if (k == "t") o.t = Numero(k, v);
        else if (k == "srm_f0") o.srmF0 = Numero(k, v);
        else if (k == "verbose") o.verbose = Numero(k, v) != 0.;
        else if (k == "vtk") o.vtk = v;
        else if (k == "csv") o.csv = v;
        else if (k == "caso") o.caso = v;
        else if (k == "nx") o.nx = int(Numero(k, v));
        else if (k == "laterais") {
            if (v != "roletes" && v != "periodicas") Erro("laterais: roletes ou periodicas");
            o.periodica = v == "periodicas";
        } else if (k == "tol_patamar") o.tolPatamar = Numero(k, v);
        else if (k == "passos") o.passos = int(Numero(k, v));
    }
    if (o.ny < 1 || o.LH <= 0. || o.H <= 0.) Erro("ny, LH e H positivos");
    if (o.nx < 0 || o.passos < 1 || !(o.tolPatamar > 0.)) Erro("nx >= 0, passos >= 1 e tol_patamar > 0");
    if (o.pressao == "perfil" && o.arquivo.empty()) Erro("pressao=perfil precisa de arquivo=");
    if (o.periodica && o.nref > 0) Erro("laterais=periodicas: so com nref=0 (refine com ny)");
    return o;
}

// ---- p_w(z), z vertical desde a superficie
struct Perfil {
    std::vector<REAL> z, pw;
    REAL fsArquivo = std::numeric_limits<REAL>::quiet_NaN(), zArquivo = 0.;
    REAL operator()(REAL zz) const
    {
        if (zz <= z.front()) return pw.front();
        if (zz >= z.back()) return pw.back();
        const size_t i = std::upper_bound(z.begin(), z.end(), zz) - z.begin();
        const REAL s = (zz - z[i - 1]) / (z[i] - z[i - 1]);
        return pw[i - 1] + s * (pw[i] - pw[i - 1]);
    }
};

// As linhas do instante t (t_h,z,y,pw,Se,theta,W,FS) do perfil do SlopeRainColumn
Perfil LePerfil(const std::string &arq, REAL t)
{
    std::ifstream in(arq);
    if (!in) Erro("nao abriu " + arq);
    std::string linha;
    std::getline(in, linha);
    if (linha.rfind("t_h,z,y,pw", 0) != 0) Erro(arq + ": cabecalho inesperado");
    Perfil p;
    REAL fsMin = std::numeric_limits<REAL>::infinity();
    while (std::getline(in, linha)) {
        std::stringstream ss(linha);
        std::string cel;
        std::vector<REAL> v;
        while (std::getline(ss, cel, ',')) v.push_back(cel == "inf" ? std::numeric_limits<REAL>::infinity() : std::stod(cel));
        if (v.size() < 8 || std::abs(v[0] - t) > 1.e-9) continue;
        if (!p.z.empty() && !(v[1] > p.z.back())) Erro(arq + ": z nao crescente");
        p.z.push_back(v[1]), p.pw.push_back(v[3]);
        if (v[1] > 0. && v[7] < fsMin) fsMin = v[7], p.zArquivo = v[1];
    }
    if (p.z.size() < 2) Erro(arq + ": sem o instante t pedido");
    p.fsArquivo = fsMin;
    return p;
}

// ---- malha da faixa [0, L] x [0, T]
TPZGeoMesh *MalhaFaixa(REAL L, REAL T, int nx, int ny)
{
    auto *g = new TPZGeoMesh();
    g->SetDimension(2);
    const int64_t nxp = nx + 1;
    g->NodeVec().Resize(nxp * (ny + 1));
    auto no = [nxp](int i, int j) { return int64_t(j) * nxp + i; };
    for (int j = 0; j <= ny; j++)
        for (int i = 0; i <= nx; i++) {
            TPZManVector<REAL, 3> c = {i == nx ? L : L * i / nx, j == ny ? T : T * j / ny, 0.};
            g->NodeVec()[no(i, j)] = TPZGeoNode(no(i, j), c, *g);
        }
    auto tri = [g](int64_t a, int64_t b, int64_t c) {
        TPZManVector<int64_t, 3> v = {a, b, c};
        new TPZGeoElRefPattern<pzgeom::TPZGeoTriangle>(v, matSolo, *g);
    };
    for (int j = 0; j < ny; j++)
        for (int i = 0; i < nx; i++) {
            const int64_t bl = no(i, j), br = no(i + 1, j), tl = no(i, j + 1), tr = no(i + 1, j + 1);
            if ((i + j) % 2 == 0) tri(bl, br, tr), tri(bl, tr, tl);
            else tri(bl, br, tl), tri(br, tr, tl);
        }
    auto linha = [g](int64_t a, int64_t b, int id) {
        TPZManVector<int64_t, 2> v = {a, b};
        new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(v, id, *g);
    };
    for (int i = 0; i < nx; i++) linha(no(i, 0), no(i + 1, 0), base);
    for (int j = 0; j < ny; j++) linha(no(nx, j), no(nx, j + 1), direita);
    for (int j = ny; j > 0; j--) linha(no(0, j), no(0, j - 1), esquerda);
    g->BuildConnectivity();
    return g;
}

TPZUnsaturatedHydraulics Hidraulica(const Opcoes &o)
{
    TPZUnsaturatedHydraulics h;
    h.fRetention = TPZUnsaturatedHydraulics::EVanGenuchten;
    h.fAVG = o.a, h.fBVG = o.nvg, h.fSwr = o.tr / o.ts;
    h.fKrw = TPZUnsaturatedHydraulics::EMualem; // nao entra (pressoes congeladas)
    h.fKra = TPZUnsaturatedHydraulics::EBrooksCorey;
    return h;
}

TUnsatHP *Material(const Opcoes &o)
{
    const REAL phi = o.phi * M_PI / 180., b = o.beta * M_PI / 180.;
    TPZElasticResponse ER;
    ER.SetEngineeringData(o.E, o.nu);
    auto *mat = new TUnsatHP(matSolo);
    mat->SetPlasticModel(TMC(TPZYCMohrCoulombPV2(phi, phi, o.c, ER), ER));
    mat->SetHydraulics(Hidraulica(o));
    mat->SetBishop(TUnsat::EChiEffective); // chi = S_e, o sigma_s = S_e p_w do Lu-Godt
    mat->SetPorosity(o.ts);
    mat->SetBiot(1., 0.);
    mat->SetConductivity(1.e-6, o.gw, 1.e-3 / 1.8e-5);
    mat->SetAir(1.2 * 9.81e-3, 101.325);
    if (o.mistura) mat->SetUnitWeightMixture(o.gammaD / (1. - o.ts)); // (1 - n) gamma_s = gamma_d; + n S_w gw
    else mat->SetUnitWeight(o.gamma);
    mat->SetGravityDirection(TPZManVector<REAL, 3>({std::sin(b), -std::cos(b), 0.}));
    mat->SetHydroMechanicalCoupling(false);
    return mat;
}

// ---- laterais periodicas: u(L, y) = u(0, y) (e as pressoes, que ja sao iguais: so dependem de y)
// Cada connect de um lado de elemento (vertice ou aresta) sobre x = L passa a ser o do lado sobre x = 0 na mesma altura:
// nos elementos das malhas atomicas (u, p_w, p_a; dominio e contorno) e nos multifisicos, que guardam os connects das
// atomicas com o offset de cada malha (refeitos aqui a partir das atomicas, como no AddConnects do neopz). Os connects de
// x = L ficam sem elemento e saem da numeracao (CleanUpUnconnectedNodes: SequenceNumber -1 e bloco vazio; a vaga no vetor
// fica). Assim o connect ic da atomica s continua sendo o offset[s] + ic da multifisica, e o UnsatAnalysis (que pula os
// connects com SequenceNumber -1) funciona sem mudar. Chamar depois da CreateCompMesh e antes de criar o UnsatAnalysis.
// So com u de ordem <= 2: a funcao de aresta de ordem 2 (4 phi_a phi_b) nao depende da orientacao da aresta, e as arestas
// de x = 0 e x = L tem nos diferentes. Devolve quantos connects de x = L foram trocados (somando as malhas)
int64_t Periodiza(TPZMultiphysicsCompMesh *mf, REAL L, REAL T)
{
    const REAL tol = 1.e-9 * std::max(L, T);
    TPZVec<TPZCompMesh *> &meshvec = mf->MeshVector();
    const int nm = meshvec.size();
    // lado de elemento sobre x = 0 (0) ou x = L (1), senao -1; chave: a dimensao do lado e as alturas dos nos
    auto lado = [&](TPZGeoEl *gel, int side, std::vector<int64_t> &chave) {
        bool x0 = true, xL = true;
        chave.assign(1, gel->SideDimension(side));
        std::vector<int64_t> ys;
        for (int i = 0; i < gel->NSideNodes(side); i++) {
            TPZManVector<REAL, 3> x(3, 0.);
            gel->SideNodePtr(side, i)->GetCoordinates(x);
            x0 = x0 && std::abs(x[0]) < tol;
            xL = xL && std::abs(x[0] - L) < tol;
            ys.push_back(std::llround(x[1] / tol));
        }
        std::sort(ys.begin(), ys.end());
        chave.insert(chave.end(), ys.begin(), ys.end());
        return x0 ? 0 : (xL ? 1 : -1);
    };
    int64_t trocados = 0;
    for (int s = 0; s < nm; s++) {
        TPZCompMesh *cm = meshvec[s];
        if (s == 0 && cm->GetDefaultOrder() > 2) Erro("laterais=periodicas: so com u de ordem <= 2");
        std::map<std::vector<int64_t>, int64_t> esquerdo; // chave -> connect do lado em x = 0
        struct Troca {
            TPZCompEl *cel;
            int loc;
            std::vector<int64_t> chave;
        };
        std::vector<Troca> trocas;
        std::set<int64_t> direitos;
        for (TPZCompEl *cel : cm->ElementVec()) {
            auto *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
            if (!intel || !cel->Reference()) continue;
            TPZGeoEl *gel = cel->Reference();
            for (int side = 0; side < gel->NSides(); side++) {
                std::vector<int64_t> chave;
                const int qual = lado(gel, side, chave);
                if (qual < 0) continue;
                const int loc = intel->SideConnectLocId(intel->NSideConnects(side) - 1, side); // o connect do proprio lado
                const int64_t ci = cel->ConnectIndex(loc);
                if (qual == 0) {
                    auto [it, novo] = esquerdo.insert({chave, ci});
                    if (!novo && it->second != ci) Erro("Periodiza: dois connects no mesmo lado de x = 0");
                } else {
                    trocas.push_back({cel, loc, chave});
                    direitos.insert(ci);
                }
            }
        }
        for (const Troca &t : trocas) {
            auto it = esquerdo.find(t.chave);
            if (it == esquerdo.end()) Erro("Periodiza: um lado de x = L sem o par em x = 0 (as malhas das laterais diferem)");
            t.cel->SetConnectIndex(t.loc, it->second);
        }
        trocados += direitos.size();
        cm->ComputeNodElCon();
        cm->CleanUpUnconnectedNodes(); // os de x = L: SequenceNumber -1
    }
    // elementos multifisicos: os connects das atomicas, em sequencia, com o offset de cada malha
    TPZManVector<int64_t, 4> offset(nm + 1, 0);
    for (int s = 0; s < nm; s++) offset[s + 1] = offset[s] + meshvec[s]->NConnects();
    for (TPZCompEl *cel : mf->ElementVec()) {
        auto *mfel = dynamic_cast<TPZMultiphysicsElement *>(cel);
        if (!mfel) continue;
        TPZManVector<int64_t, 32> indices;
        int64_t n = 0;
        for (int s = 0; s < nm; s++)
            if (TPZCompEl *ref = mfel->ReferredElement(s)) n += ref->NConnects();
        if (n != cel->NConnects()) Erro("Periodiza: elemento multifisico com outros connects que os das atomicas");
        indices.Resize(n);
        n = 0;
        for (int s = 0; s < nm; s++)
            if (TPZCompEl *ref = mfel->ReferredElement(s))
                for (int i = 0; i < ref->NConnects(); i++) indices[n++] = offset[s] + ref->ConnectIndex(i);
        mfel->SetConnectIndexes(indices);
    }
    mf->ComputeNodElCon();
    mf->CleanUpUnconnectedNodes(); // os de x = L saem da numeracao (as equacoes livres ficam em sequencia)
    return trocados;
}

// FS fechado de talude infinito (Lu-Godt), min em z, com a p_w(z) dada: W por trapezios numa grade de N intervalos
struct FSz {
    REAL fs = std::numeric_limits<REAL>::infinity(), z = 0.;
};
FSz FSFechado(const Opcoes &o, const std::function<REAL(REAL)> &pw)
{
    const TPZUnsaturatedHydraulics h = Hidraulica(o);
    const REAL b = o.beta * M_PI / 180., cb = std::cos(b), sb = std::sin(b), tphi = std::tan(o.phi * M_PI / 180.);
    auto gam = [&](REAL z) {
        if (!o.mistura) return o.gamma;
        return o.gammaD + o.ts * h.Evaluate(pw(z), 0.).S * o.gw;
    };
    const int N = 4000;
    REAL W = 0.;
    FSz r;
    for (int i = 1; i <= N; i++) {
        const REAL z0 = o.H * (i - 1) / N, z1 = o.H * i / N;
        W += 0.5 * (gam(z0) + gam(z1)) * (z1 - z0);
        const REAL p = pw(z1);
        REAL Se = 1., dSe;
        if (p < 0.) h.Effective(h.Evaluate(p, 0.).S, Se, dSe);
        const REAL ss = p < 0. ? Se * p : p;
        const REAL fs = (o.c + (W * cb * cb - ss) * tphi) / (W * sb * cb);
        if (fs < r.fs) r.fs = fs, r.z = z1;
    }
    return r;
}

struct Fator {
    std::unique_ptr<UnsatModel<TUnsatHP>> m;
    std::unique_ptr<UnsatAnalysis<TUnsatHP>> a;
};

int main(int argc, char *argv[])
{
    const auto t0 = std::chrono::steady_clock::now();
    const Opcoes o = LeOpcoes(argc, argv);
    const REAL b = o.beta * M_PI / 180., cb = std::cos(b), T = o.H * cb;
    // nx de LH (celulas quase quadradas) ou dado (quadradas, e L sai dai)
    const int nx = o.nx > 0 ? o.nx : std::max(1, int(std::lround(o.LH * o.H / (T / o.ny))));
    const REAL L = o.nx > 0 ? o.nx * T / o.ny : o.LH * o.H, LH = L / o.H;
    if (o.periodica && nx < 2) Erro("laterais=periodicas: nx >= 2 (nenhum elemento toca as duas laterais)");

    // p_w(z)
    Perfil perfil;
    if (o.pressao == "perfil") perfil = LePerfil(o.arquivo, o.t);
    const std::function<REAL(REAL)> pw = [&](REAL z) -> REAL {
        if (o.pressao == "seco") return 0.;
        if (o.pressao == "lencol") {
            const REAL p = o.gw * (z - o.dw) * cb * cb;
            return (p >= 0. || o.acimaHidro) ? p : 0.;
        }
        return perfil(z);
    };

    // FS fechado (Lu-Godt): com a p_w(z) dada (FS_LEM) e com a que a malha de partida recebe (FS_P1): p_w(z) nos nos,
    // z_j = (T - y_j) / cos beta como na SetFrozenPressures, linear entre eles
    const FSz lem = FSFechado(o, pw);
    Perfil nos;
    for (int j = o.ny; j >= 0; j--) {
        const REAL y = j == o.ny ? T : T * j / o.ny, z = (T - y) / cb;
        nos.z.push_back(z), nos.pw.push_back(pw(z));
    }
    const FSz p1 = FSFechado(o, [&nos](REAL z) { return nos(z); });

    // SRM na faixa, com refinamento adaptativo opcional
    Fator f;
    TPZGeoMesh *g = MalhaFaixa(L, T, nx, o.ny);
    REAL fs = 0.;
    std::vector<REAL> niveis;
    int64_t trocados = 0;
    for (int k = 0;; k++) {
        f.a.reset();
        f.m = std::make_unique<UnsatModel<TUnsatHP>>();
        f.m->SetGeoMesh(g);
        f.m->AddDirichletBC(base, TUnsat::EDirichletU, 0, {0, 1});
        if (!o.periodica) {
            f.m->AddDirichletBC(esquerda, TUnsat::EDirichletU, 0, {0});
            f.m->AddDirichletBC(direita, TUnsat::EDirichletU, 0, {0});
        }
        TUnsatHP *mat = Material(o);
        f.m->SetMaterial(mat);
        f.m->CreateCompMesh({{2, 2}, {1, 1}, {1, 1}});
        if (o.periodica) trocados = Periodiza(f.m->CompMesh(), L, T);
        f.a = std::make_unique<UnsatAnalysis<TUnsatHP>>(*f.m);
        f.a->SetVerbose(o.verbose);
        f.a->SetLoadFactorFunction([mat](REAL l) { mat->SetGravityFactor(l); }, false);
        f.a->SetStrengthReductionFunction([mat](REAL F) { mat->PlasticModel().SetStrengthReductionFactor(F); });
        f.a->SetSRMStart(o.srmF0);
        f.a->SetPlateauTolerance(o.tolPatamar);
        f.a->SetMaxPathSteps(o.passos);
        f.a->SetFrozenPressures([&](int espaco, const TPZVec<REAL> &x) {
            if (espaco != 1) return 0.;
            return pw((T - x[1]) / cb);
        });
        fs = f.a->SRMArcLength();
        niveis.push_back(fs);
        if (k == o.nref) break;
        const std::vector<REAL> ind = f.a->PlasticIndicator();
        const REAL vmax = *std::max_element(ind.begin(), ind.end());
        std::set<int64_t> marcados;
        for (TPZGeoEl *gel : g->ElementVec())
            if (gel && gel->Dimension() == 2 && ind[gel->Index()] >= 0.1 * vmax) marcados.insert(gel->Index());
        std::vector<int64_t> pai;
        g = f.m->RefinedGeoMesh(marcados, pai);
    }

    // fim do caminho do SRM (o da ultima malha)
    const auto &caminho = f.a->PathHistory();
    const REAL nan = std::numeric_limits<REAL>::quiet_NaN();
    REAL dFult = nan, dFant = nan, crit = nan;
    if (!caminho.empty()) {
        const auto &u = caminho.back();
        dFult = u.dv;
        crit = u.l > 0. ? u.dv * u.lPedido / u.l : u.dv; // o passo vertical nao tem arco
        if (caminho.size() > 1) dFant = caminho[caminho.size() - 2].dv;
    }

    // banda no terco central (com laterais periodicas, na faixa toda: nenhum x e especial)
    const std::vector<REAL> ind = f.a->PlasticIndicator();
    REAL vmax = 0., zb = 0.;
    std::vector<std::pair<REAL, REAL>> centrais; // (indicador, z)
    int64_t nel = 0;
    for (TPZGeoEl *gel : g->ElementVec()) {
        if (!gel || gel->Dimension() != 2 || gel->HasSubElement()) continue;
        nel++;
        TPZManVector<REAL, 3> qsi(2), x(3);
        gel->CenterPoint(gel->NSides() - 1, qsi);
        gel->X(qsi, x);
        if (!o.periodica && (x[0] < L / 3. || x[0] > 2. * L / 3.)) continue;
        const REAL z = (T - x[1]) / cb, v = ind[gel->Index()];
        centrais.push_back({v, z});
        if (v > vmax) vmax = v, zb = z;
    }
    REAL zmin = o.H, zmax = 0.;
    for (const auto &[v, z] : centrais)
        if (v >= 0.5 * vmax) zmin = std::min(zmin, z), zmax = std::max(zmax, z);
    const int64_t neq = int64_t(f.a->FreeEquations().size());
    const REAL seg = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

    std::printf("Faixa: H = %g m (T = %.4f m normal), beta = %g, L/H = %g (L = %g m), %d x %d celulas, laterais %s, "
                "nref %d: %lld triangulos, %lld equacoes livres | c' %g, phi' %g, peso %s | pressao %s%s\n",
                o.H, T, o.beta, LH, L, nx, o.ny, o.periodica ? "periodicas" : "roletes", o.nref, (long long)nel,
                (long long)neq, o.c, o.phi, o.mistura ? "mistura" : "constante", o.pressao.c_str(),
                o.pressao == "perfil" ? (" (" + o.arquivo + ", t = " + std::to_string(o.t) + " h)").c_str() : "");
    if (o.periodica) std::printf("  laterais periodicas: %lld connects de x = L trocados pelos de x = 0\n", (long long)trocados);
    std::printf("  FS_SRM = %.5f", fs);
    if (niveis.size() > 1) {
        std::printf(" (niveis:");
        for (REAL v : niveis) std::printf(" %.5f", v);
        std::printf(")");
    }
    std::printf(" | FS_LEM = min_z Lu-Godt = %.5f em z = %.4f m | erro %+.3f%% | banda %s z = %.3f m (%.3f a "
                "%.3f) | %.1f s\n",
                lem.fs, lem.z, 100. * (fs / lem.fs - 1.), o.periodica ? "" : "central", zb, zmin, zmax, seg);
    std::printf("  FS_P1 (p_w dos nos, linear) = %.5f em z = %.4f m | erro %+.3f%% | caminho: fim %s, %zu passos, dF "
                "%.2e (antes %.2e), dF lPedido / l = %.2e (tol %.0e)\n",
                p1.fs, p1.z, 100. * (fs / p1.fs - 1.), f.a->PathEnd().c_str(), caminho.size(), dFult, dFant, crit,
                o.tolPatamar);
    if (o.pressao == "perfil")
        std::printf("  FS_min do arquivo (SlopeRainColumn): %.5f em z = %.4f m\n", perfil.fsArquivo, perfil.zArquivo);

    if (!o.vtk.empty()) f.a->PostProcess(o.vtk, {"WaterPressure"}, {"Displacement"});

    const std::string cab = "caso,laterais,LH,nx,ny,nref,elementos,equacoes,FS_SRM,FS_LEM,z_LEM,FS_P1,z_P1,z_banda,"
                            "z_banda_min,z_banda_max,erro_rel,erro_P1,fim,passos,dF_ult,dF_ant,crit_ult,tol_patamar,"
                            "segundos";
    std::string primeira;
    {
        std::ifstream in(o.csv);
        if (in) std::getline(in, primeira);
    }
    if (!primeira.empty() && primeira != cab) Erro("csv=" + o.csv + ": o arquivo tem outro cabecalho (de outra versao)");
    std::ofstream csv(o.csv, std::ios::app);
    if (primeira.empty()) csv << cab << "\n";
    char lin[1024];
    std::snprintf(lin, sizeof(lin),
                  "%s,%s,%g,%d,%d,%d,%lld,%lld,%.6f,%.6f,%.4f,%.6f,%.4f,%.4f,%.4f,%.4f,%.6f,%.6f,%s,%zu,%.3e,%.3e,%.3e,%g,"
                  "%.1f\n",
                  o.caso.c_str(), o.periodica ? "periodicas" : "roletes", LH, nx, o.ny, o.nref, (long long)nel,
                  (long long)neq, fs, lem.fs, lem.z, p1.fs, p1.z, zb, zmin, zmax, fs / lem.fs - 1., fs / p1.fs - 1.,
                  f.a->PathEnd().c_str(), caminho.size(), dFult, dFant, crit, o.tolPatamar, seg);
    csv << lin;
    return 0;
}
