// Rahardjo2007: a reproducao da Serie B de H. Rahardjo, T. H. Ong, R. B. Rezaur e E. C. Leong (2007), "Factors
// controlling instability of homogeneous soil slopes under rainfall", J. Geotech. Geoenviron. Eng. 133(12):1532-1543.
//
// Serie B (Tabela 2 do artigo): talude homogeneo de H_s = 10 m a 45 graus (Fig. 2), nivel d'agua inicial H_w = 5 m
// abaixo do pe (subindo 7 graus para dentro do macico), succao inicial limitada a 75 kPa, chuva constante I_r por 24 h,
// os tres solos da Tabela 1 (f10,-4, f50,-5 e f100,-6) e sete intensidades por solo; o FS ao longo de 48 h (Fig. 5), o
// FS minimo contra I_r (Figs. 6 e 7) e a posicao do nivel d'agua no FS minimo (Tabelas 3 e 5).
//
// Os modelos do transiente (os do docs/nao_saturado/secoes_rahardjo.tex e do RainTransient):
//   iii (o padrao): esqueleto rigido e ar passivo, so a percolacao da agua (o SEEP/W do artigo);
//   ii:  acoplado com ar passivo: o esqueleto elastico (E, nu) e resolvido junto, p_a = 0 em todo o dominio;
//   i:   o modelo completo: acoplado e com o ar ativo (p_a livre, p_a = 0 so nos nos da superficie; k_ra de Brooks-Corey
//        com lambda = 3 e piso 1e-4; mu_w / mu_a = 55,6).
// Nos acoplados o transiente parte do estado geostatico (o peso com as pressoes iniciais congeladas, depois u = 0 com a
// tensao mantida) e o esqueleto e elastico (o Mohr-Coulomb com c' = 1e9: o transiente so da as pressoes e a
// deformacao). O FS e sempre o SRM (o artigo usa o Bishop simplificado do SLOPE/W), com as pressoes do transiente
// congeladas (no i, p_w e p_a; com chi = 1 a tensao efetiva e sigma - p_w e o p_a nao entra na resistencia).
//
// Uso: Rahardjo2007 [solo=f10|f50|f100] [I=<mm/h>] [modelo=iii|ii|i] [E=2e4] [nu=0.3] [tempos=artigo|<k>h|t1,...]
//                   [horas=48] [dt=0.1] [nref=4] [fs0=0.8] [dx=1] [eta0=0.003] [r=1.15] [etamax=0.025] [saida=rah]
//                   [so_transiente=0|1] [grava_malha=<prefixo>] [campos=t1,t2,...]
//   I: a intensidade (mm/h; padrao 1 k_s do solo); tempos: os instantes do FS (artigo: os da Fig. 5, 0, 1, 2, 3, 6, 12,
//   24, 28, 32, 40 e 48 h; <k>h: de k em k horas ate horas; ou a lista); dt: o passo do transiente (h); nref: os niveis
//   do refinamento adaptativo do fator; dx, eta0, r, etamax: a malha do transiente (colunas a cada dx m; camadas em
//   fracoes da altura de cada coluna, a primeira eta0 junto a superficie, crescendo na razao r ate etamax); campos: os
//   instantes (entre os do FS) em que grava o campo de p_w (<pre>_campo_<t>h.csv) e a malha final do fator com o
//   indicador plastico (<pre>_fator_<t>h.csv)
//   E, nu: o esqueleto elastico dos modelos acoplados (kPa; o artigo nao da: os do det_3solos)
//   fs0: o F com que o SRM aplica o peso antes do caminho de F (SetSRMStart; 0,8 e o padrao da UnsatAnalysis; se o
//   talude nao fica em pe com ele, a UnsatAnalysis tenta fs0/2). Testado com 0,1: o mesmo FS (0,8260 contra 0,8259 no
//   patamar saturado do f10 a 900 mm/h em 12 h; 2,0983 contra 2,0980 em 0 h), 15 a 25 % mais caro por FS
// Saidas: <pre>_fs.csv (FS, niveis do refinamento, pressoes, nivel d'agua na crista, no meio da face e no pe, categoria
//   da Tabela 3, tempos, p_a) e <pre>_nivel.csv (o nivel d'agua de metro em metro em cada instante do FS), com <pre> =
//   <saida>_<solo>_<I> no modelo iii e <saida>_<modelo>_<solo>_<I> nos outros
//
// As pecas sao as do SlopeRandomFieldRefinement (UnsatModel, UnsatAnalysis, HPMaterial, LayerGMesh, PointLocator,
// RainTransient) e o material e o fator do rahardjo.cpp de la (o driver deterministico da chuva), sem mudancas.

#include "HPMaterial.h"
#include "LayerGMesh.h"
#include "PointLocator.h"
#include "RainTransient.h"
#include "UnsatAnalysis.h"
#include "UnsatModel.h"

#include "Plasticity/TPZElasticResponse.h"
#include "Plasticity/TPZMatUnsaturatedUPwPa.h"
#include "Plasticity/TPZPlasticStepVoigt.h"
#include "Plasticity/TPZYCMohrCoulombPV2.h"
#include "TPZGeoLinear.h"
#include "TPZSpStructMatrix.h"
#include "pzgeotriangle.h"
#include "pzstepsolver.h"
#include "tpzgeoelrefpattern.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using TMC = TPZPlasticStepVoigt<TPZYCMohrCoulombPV2, TPZElasticResponse>;
using TUnsat = TPZMatUnsaturatedUPwPa<TMC, TPZUnsaturatedMem>;
using TUnsatHP = HPMaterial<TUnsat>; // regra de integracao 2 (p_u - 1) em cada elemento

/// ids de contorno (os da TriGMesh): -1 base, -2 direita, -3 terreno do pe, -4 crista, -5 esquerda, -6 face
enum EIds { matSolo = 1, base = -1, right = -2, topright = -3, topleft = -4, left = -5, face = -6 };
constexpr REAL gammaW = 9.81;

// ---- o talude da Fig. 2 (H = 10 m a 45 graus): cantos no sentido anti-horario (base esquerda, base direita, terreno
// do pe a direita, pe, crista, crista a esquerda); a crista em x = 30 m, o pe em x = 40 m
const std::vector<std::vector<REAL>> gTalude = {{0, 0}, {70, 0}, {70, 20}, {40, 20}, {30, 30}, {0, 30}};
constexpr REAL xCrista = 30., xPe = 40., yCrista = 30., yPe = 20., xFim = 70.;
static REAL Superficie(REAL x) { return x <= xCrista ? yCrista : x >= xPe ? yPe : yCrista - (x - xCrista); }

// ---- os solos da Tabela 1 (o que muda de um para outro: o a de Fredlund-Xing e o k_s) e as intensidades da Serie B
struct Solo {
    std::string nome;
    REAL a;  // kPa
    REAL ks; // m/s
};
static const std::vector<Solo> gSolos = {{"f10", 10., 1.e-4}, {"f50", 50., 1.e-5}, {"f100", 100., 1.e-6}};
// o que e igual nos tres (kN, m, kPa); E e nu nao estao no artigo (equilibrio limite) e nao mudam o FS
struct SoloComum {
    REAL c = 10., phi = 26. * M_PI / 180.; // c' e phi' (phi^b = phi' = 26 graus: chi = 1)
    REAL gamma = 20.;                       // o mesmo acima e abaixo do nivel d'agua
    REAL E = 2.e4, nu = 0.3;
};

struct Opcoes {
    Solo solo = gSolos[1];
    REAL I = -1.; // mm/h (< 0: 1 k_s)
    std::string modelo = "iii"; // iii, ii ou i
    REAL E = SoloComum().E, nu = SoloComum().nu; // o esqueleto elastico dos acoplados
    bool Acoplado() const { return modelo != "iii"; }
    bool ArAtivo() const { return modelo == "i"; }
    ERainModel Modelo() const {
        return ArAtivo() ? ERainModel::ECoupledActiveAir : Acoplado() ? ERainModel::ECoupledPassiveAir
                                                                      : ERainModel::EUncoupled;
    }
    std::vector<REAL> tempos;
    std::string tempos_txt = "artigo";
    REAL horas = 48., dt = 0.1;
    int nref = 4;
    REAL fs0 = 0.8; // o F do peso no SRM
    REAL dx = 1., eta0 = 0.003, r = 1.15, etamax = 0.025;
    std::string saida = "rah", gravaMalha;
    bool soTransiente = false;
    std::vector<REAL> campos; // instantes em que grava o campo de p_w e a malha final do fator
};

static std::vector<REAL> Lista(const std::string &s) {
    std::vector<REAL> v;
    std::stringstream ss(s);
    for (std::string t; std::getline(ss, t, ',');) v.push_back(std::stod(t));
    return v;
}

static Opcoes LeOpcoes(int argc, char **argv) {
    Opcoes op;
    for (int i = 1; i < argc; i++) {
        const std::string s = argv[i];
        const size_t eq = s.find('=');
        const std::string k = s.substr(0, eq), v = eq == std::string::npos ? "" : s.substr(eq + 1);
        bool ok = true;
        if (k == "solo") {
            ok = false;
            for (const Solo &so : gSolos)
                if (so.nome == v) op.solo = so, ok = true;
        } else if (k == "I") op.I = std::stod(v);
        else if (k == "modelo") op.modelo = v, ok = v == "iii" || v == "ii" || v == "i";
        else if (k == "E") op.E = std::stod(v);
        else if (k == "nu") op.nu = std::stod(v);
        else if (k == "tempos") op.tempos_txt = v;
        else if (k == "horas") op.horas = std::stod(v);
        else if (k == "dt") op.dt = std::stod(v);
        else if (k == "nref") op.nref = std::stoi(v);
        else if (k == "fs0") op.fs0 = std::stod(v);
        else if (k == "dx") op.dx = std::stod(v);
        else if (k == "eta0") op.eta0 = std::stod(v);
        else if (k == "r") op.r = std::stod(v);
        else if (k == "etamax") op.etamax = std::stod(v);
        else if (k == "saida") op.saida = v;
        else if (k == "grava_malha") op.gravaMalha = v;
        else if (k == "so_transiente") op.soTransiente = std::stoi(v) != 0;
        else if (k == "campos") op.campos = Lista(v);
        else ok = false;
        if (!ok) {
            std::fprintf(stderr, "opcao invalida: %s\n", argv[i]);
            std::exit(1);
        }
    }
    if (op.I < 0.) op.I = op.solo.ks * 3.6e6;
    if (op.tempos_txt == "artigo") op.tempos = {0, 1, 2, 3, 6, 12, 24, 28, 32, 40, 48};
    else if (op.tempos_txt.back() == 'h') {
        const REAL k = std::stod(op.tempos_txt.substr(0, op.tempos_txt.size() - 1));
        for (int i = 0; i * k <= op.horas + 1.e-9; i++) op.tempos.push_back(i * k);
    } else op.tempos = Lista(op.tempos_txt);
    op.tempos.erase(std::remove_if(op.tempos.begin(), op.tempos.end(), [&](REAL t) { return t > op.horas + 1.e-9; }),
                    op.tempos.end());
    return op;
}

// ---- 1. as malhas
// transiente: camadas que seguem o terreno (LayerGMesh), colunas a cada dx m passando pela crista e pelo pe, camadas em
// fracoes da altura de cada coluna, finas junto a superficie (a primeira com eta0) e crescendo na razao r ate etamax
static TPZGeoMesh *MalhaTransiente(REAL dx, REAL eta0, REAL r, REAL etamax) {
    LayerGMeshData d;
    d.superficie = {{0., yCrista}, {xCrista, yCrista}, {xPe, yPe}, {xFim, yPe}};
    d.idSuperficie = {topleft, face, topright};
    d.idMaterial = matSolo, d.idBase = base, d.idEsquerda = left, d.idDireita = right;
    // colunas: cada trecho (0-30, 30-40, 40-70) em partes de ~dx
    for (auto [a, b] : std::vector<std::pair<REAL, REAL>>{{0., xCrista}, {xCrista, xPe}, {xPe, xFim}}) {
        const int n = std::max(1, int(std::lround((b - a) / dx)));
        for (int i = d.xs.empty() ? 0 : 1; i <= n; i++) d.xs.push_back(a + (b - a) * i / n);
    }
    d.eta = {1.};
    for (REAL q : GeometricSteps(1., eta0, r, etamax)) d.eta.push_back(d.eta.back() - q);
    d.eta.back() = 0.;
    std::reverse(d.eta.begin(), d.eta.end());
    d.quadrilateros = false;
    return LayerGMesh(d);
}

// fator: a malha inicial do rahardjo.cpp (TriGMesh com 1 divisao uniforme: lados de ~5 m), refinada depois pelo
// indicador plastico
TPZGeoMesh *TriGMesh(const std::vector<std::vector<REAL>> &points, int ref);

// ---- 2. o material (o MaterialRahardjo do rahardjo.cpp): Mohr-Coulomb associado; Fredlund-Xing com C(psi) = 1, n = m
// = 1 e theta_s = 0,45; k_w = k_s Theta^4; ar passivo; chi = tan(phi^b)/tan(phi') = 1 (a parcela da succao da Eq. 5 do
// artigo); sem acoplamento: no transiente so a agua percola (esqueleto parado) e no fator as pressoes ficam congeladas
static TUnsatHP *Material(const RainData &r) {
    const SoloComum s;
    TPZElasticResponse ER;
    ER.SetEngineeringData(s.E, s.nu);
    auto *mat = new TUnsatHP(matSolo);
    mat->SetPlasticModel(TMC(TPZYCMohrCoulombPV2(s.phi, s.phi, s.c, ER), ER));
    TPZUnsaturatedHydraulics h;
    h.fRetention = TPZUnsaturatedHydraulics::EFredlundXing;
    h.fA = r.aFX, h.fN = r.nFX, h.fM = r.mFX;
    h.fKrw = TPZUnsaturatedHydraulics::EPowerLaw;
    h.fPower = r.power;
    h.fKra = TPZUnsaturatedHydraulics::EBrooksCorey; // o ar fica passivo (p_a = 0)
    mat->SetHydraulics(h);
    mat->SetBishop(TUnsat::EChiConstant, std::tan(r.phib) / std::tan(s.phi));
    mat->SetPorosity(r.n);
    mat->SetBiot(1., 0.);
    mat->SetWaterBulkModulus(r.Kw);
    mat->SetConductivity(r.ks, gammaW, 1.e-3 / 1.8e-5);
    mat->SetAir(1.2 * 9.81e-3, 101.325);
    mat->SetUnitWeight(s.gamma);
    mat->SetLumpedStorage(true);
    mat->SetHydroMechanicalCoupling(false);
    return mat;
}

// o material da malha do transiente: no iii, o Material (o esqueleto fica parado: a resistencia nao entra); nos
// acoplados, o acoplamento das equacoes de massa ligado e o esqueleto elastico (o mesmo Mohr-Coulomb com c' = 1e9: a
// projecao e sempre a elastica)
static TUnsatHP *MaterialTransiente(const RainData &r, const Opcoes &op) {
    TUnsatHP *mat = Material(r);
    if (!op.Acoplado()) return mat;
    TPZElasticResponse ER;
    ER.SetEngineeringData(op.E, op.nu);
    const REAL phi = SoloComum().phi;
    mat->SetPlasticModel(TMC(TPZYCMohrCoulombPV2(phi, phi, 1.e9, ER), ER));
    mat->SetHydroMechanicalCoupling(true); // os termos S_pi alfa d(eps_v) nas equacoes de massa
    return mat;
}

// u: base engastada, laterais presas so em x
static void Apoios(UnsatModel<TUnsatHP> &m) {
    m.AddDirichletBC(right, TUnsat::EDirichletU, 0, {0});
    m.AddDirichletBC(left, TUnsat::EDirichletU, 0, {0});
    m.AddDirichletBC(base, TUnsat::EDirichletU, 0, {0, 1});
}

// ---- 3. o fator (o FatorAdaptativo do rahardjo.cpp): SRM com as pressoes p_w do transiente (uma por no da malha
// dele), refinamento adaptativo refeito do zero a cada instante: TriGMesh(1) e nref niveis em h (marcados:
// sqrt(J2(eps_p)) >= 0,1 do maior). niveis: o FS de cada nivel; pa: p_a nos nos do transiente (ar ativo; vazio: 0)
static REAL Fator(const RainData &r, const PointLocator &loc, const std::vector<REAL> &pw, const std::vector<REAL> &pa,
                  REAL fs0, int nref, std::vector<REAL> &niveis, const std::string &arqMalha = "") {
    TPZGeoMesh *g = TriGMesh(gTalude, 1);
    std::unique_ptr<UnsatModel<TUnsatHP>> m;
    std::unique_ptr<UnsatAnalysis<TUnsatHP>> a;
    REAL fs = 0.;
    niveis.clear();
    for (int k = 0;; k++) {
        a.reset(); // a analise antes do modelo (ela aponta para ele)
        m = std::make_unique<UnsatModel<TUnsatHP>>();
        m->SetGeoMesh(g);
        Apoios(*m);
        m->SetMaterial(Material(r));
        m->CreateCompMesh({{2, 2}, {1, 1}, {1, 1}});
        a = std::make_unique<UnsatAnalysis<TUnsatHP>>(*m);
        a->SetVerbose(std::getenv("SRM_VERBOSE") != nullptr);
        TUnsatHP *mat = m->Material();
        a->SetLoadFactorFunction([mat](REAL l) { mat->SetGravityFactor(l); }, false);
        a->SetStrengthReductionFunction([mat](REAL F) { mat->PlasticModel().SetStrengthReductionFactor(F); });
        a->SetSRMStart(fs0);
        a->SetFrozenPressures([&](int espaco, const TPZVec<REAL> &x) {
            if (espaco == 2 && !pa.empty()) return loc.Interpolate(x[0], x[1], pa);
            return espaco == 1 ? loc.Interpolate(x[0], x[1], pw) : 0.;
        });
        fs = a->SRMArcLength();
        niveis.push_back(fs);
        const std::vector<REAL> ind = a->PlasticIndicator();
        if (k == nref) {
            if (!arqMalha.empty()) { // a malha final do fator com o indicador plastico de cada triangulo
                std::ofstream f(arqMalha);
                f << "x1,y1,x2,y2,x3,y3,sqrtJ2_eps_p\n";
                for (TPZGeoEl *gel : g->ElementVec()) {
                    if (!gel || gel->Dimension() != 2 || gel->HasSubElement() || gel->NCornerNodes() != 3) continue;
                    for (int i = 0; i < 3; i++)
                        f << gel->NodePtr(i)->Coord(0) << "," << gel->NodePtr(i)->Coord(1) << ",";
                    f << ind[gel->Index()] << "\n";
                }
            }
            break;
        }
        const REAL vmax = *std::max_element(ind.begin(), ind.end());
        std::set<int64_t> marcados;
        for (TPZGeoEl *gel : g->ElementVec())
            if (gel && gel->Dimension() == 2 && ind[gel->Index()] >= 0.1 * vmax) marcados.insert(gel->Index());
        std::vector<int64_t> pai;
        g = m->RefinedGeoMesh(marcados, pai);
    }
    return fs;
}

// ---- 4. o nivel d'agua na abscissa x: o topo da zona saturada que vem da base (subindo, o primeiro p_w < 0); a
// superficie se p_w >= 0 ate ela
static REAL Nivel(const PointLocator &loc, const std::vector<REAL> &pw, REAL x) {
    const REAL ys = Superficie(x) - 1.e-6, h = 0.02;
    REAL y0 = 0., p0 = loc.Interpolate(x, 0., pw);
    if (p0 < 0.) return 0.;
    for (REAL y = h;; y += h) {
        const REAL y1 = std::min(y, ys), p1 = loc.Interpolate(x, y1, pw);
        if (p1 < 0.) return y0 + (y1 - y0) * p0 / (p0 - p1);
        if (y1 >= ys) return Superficie(x);
        y0 = y1, p0 = p1;
    }
}
// a categoria da Tabela 3 do artigo: I, o nivel chega a crista; II, ao meio da face; III, ao pe; IV, nao chega
static std::string Categoria(REAL nCrista, REAL nMeio, REAL nPe) {
    const REAL tol = 0.05; // m
    if (nCrista >= yCrista - tol) return "I";
    if (nMeio >= 0.5 * (yCrista + yPe) - tol) return "II";
    if (nPe >= yPe - tol) return "III";
    return "IV";
}

// os triangulos de uma malha (x1,y1,x2,y2,x3,y3 por linha), para a figura do relatorio
static void GravaMalha(TPZGeoMesh *g, const std::string &arq) {
    std::ofstream f(arq);
    f << "x1,y1,x2,y2,x3,y3\n";
    for (TPZGeoEl *gel : g->ElementVec()) {
        if (!gel || gel->Dimension() != 2 || gel->HasSubElement()) continue;
        for (int i = 0; i < gel->NCornerNodes(); i++)
            f << (i ? "," : "") << gel->NodePtr(i)->Coord(0) << "," << gel->NodePtr(i)->Coord(1);
        f << "\n";
    }
}

int main(int argc, char **argv) {
    using Clock = std::chrono::steady_clock;
    const auto t0 = Clock::now();
    const Opcoes op = LeOpcoes(argc, argv);
    if (!op.gravaMalha.empty()) {
        std::unique_ptr<TPZGeoMesh> a(MalhaTransiente(op.dx, op.eta0, op.r, op.etamax)), b(TriGMesh(gTalude, 1));
        GravaMalha(a.get(), op.gravaMalha + "_transiente.csv");
        GravaMalha(b.get(), op.gravaMalha + "_fator.csv");
        return 0;
    }
    const Solo &s = op.solo;
    RainData r;
    r.aFX = s.a, r.ks = s.ks;
    r.intensity = op.I, r.duration = 24.;
    r.dt = op.dt, r.dtMin = op.dt / 64.;
    r.xToe = xPe, r.yToe = yPe; // nivel 5 m abaixo do pe, subindo 7 graus; succao ate 75 kPa (os padroes do RainData)
    std::printf("Rahardjo (2007), Serie B: solo %s (a = %g kPa, k_s = %g m/s), I_r = %g mm/h (%.4g k_s) por 24 h, "
                "transiente ate %g h (dt = %g h), modelo (%s)%s\n",
                s.nome.c_str(), s.a, s.ks, op.I, op.I / (s.ks * 3.6e6), op.horas, op.dt, op.modelo.c_str(),
                op.ArAtivo() ? ": acoplado, ar ativo (p_a = 0 so na superficie)"
                : op.Acoplado() ? ": acoplado, ar passivo" : ": esqueleto rigido, ar passivo");
    if (op.Acoplado()) std::printf("  esqueleto elastico: E = %g kPa, nu = %g\n", op.E, op.nu);
    std::printf("  SRM: peso com F = %g, %d niveis adaptativos\n", op.fs0, op.nref);

    // transiente: no iii so a agua, com u e p_a presos; nos acoplados u (so os apoios presos) e, no i, p_a (presa so
    // na superficie)
    UnsatModel<TUnsatHP> mT;
    mT.SetGeoMesh(MalhaTransiente(op.dx, op.eta0, op.r, op.etamax));
    Apoios(mT);
    mT.SetMaterial(MaterialTransiente(r, op));
    mT.CreateCompMesh({{2, 2}, {1, 1}, {1, 1}});
    TPZUnsatAnalysis aT(mT.CompMesh(), {{base, {true, true}}, {right, {true, false}}, {left, {true, false}}});
    TPZSpStructMatrix<STATE> strmat(mT.CompMesh());
    strmat.SetNumThreads(std::getenv("PZ_THREADS") ? std::atoi(std::getenv("PZ_THREADS")) : 8);
    aT.SetStructuralMatrix(strmat);
    TPZStepSolver<STATE> solver;
    solver.SetDirect(ELU);
    aT.SetSolver(solver);
    RainTransient chuva(aT, r, gammaW, {left, right}, {topleft, face, topright}, {}, op.Modelo());
    std::printf("  malha do transiente: %lld elementos, %lld equacoes\n", (long long)mT.GeoMesh()->NElements(),
                (long long)mT.CompMesh()->NEquations());
    if (!chuva.InitialState()) {
        std::printf("  o estado inicial falhou: %s\n", chuva.LastFailure().c_str());
        return 1;
    }
    const PointLocator loc(mT.GeoMesh());

    char sI[32];
    std::snprintf(sI, sizeof sI, "%g", op.I);
    const std::string pre = op.saida + (op.modelo == "iii" ? "" : "_" + op.modelo) + "_" + s.nome + "_" + sI;
    std::ofstream ff(pre + "_fs.csv"), fn(pre + "_nivel.csv");
    ff << "t_h,FS,niveis,pw_min,pw_max,nos_face,infiltrado_m3m,nivel_crista,nivel_meio,nivel_pe,categoria,passos,"
          "s_transiente,s_fator,s_total,pa_min,pa_max\n";
    fn << "t_h,x,y_nivel\n";
    std::printf("  %5s %8s %8s %8s %8s %6s %7s %7s %7s %4s %8s %8s\n", "t(h)", "FS", "pw min", "pw max", "pa max", "face",
                "crista", "meio", "pe", "cat", "s fator", "s total");
    double sTrans = 0.;
    for (REAL t : op.tempos) {
        if (t > 0.) {
            const auto ta = Clock::now();
            const bool ok = chuva.AdvanceTo(t);
            sTrans += std::chrono::duration<double>(Clock::now() - ta).count();
            if (!ok) {
                std::printf("  o transiente parou em %.3f h: %s\n", chuva.Time(), chuva.LastFailure().c_str());
                return 1;
            }
        }
        const std::vector<REAL> pw = chuva.NodalWaterPressure();
        const std::vector<REAL> pa = op.ArAtivo() ? chuva.NodalAirPressure() : std::vector<REAL>();
        const REAL nC = Nivel(loc, pw, xCrista), nM = Nivel(loc, pw, 0.5 * (xCrista + xPe)), nP = Nivel(loc, pw, xPe);
        for (int x = 0; x <= int(xFim); x++) fn << t << "," << x << "," << Nivel(loc, pw, x) << "\n";
        const bool grava = std::any_of(op.campos.begin(), op.campos.end(), [t](REAL c) { return std::abs(c - t) < 1.e-9; });
        char tt[32];
        std::snprintf(tt, sizeof tt, "%g", t);
        if (grava) { // o campo de p_w: cada triangulo da malha do transiente com p_w nos vertices
            std::ofstream fc(pre + "_campo_" + tt + "h.csv");
            fc << "x1,y1,p1,x2,y2,p2,x3,y3,p3" << (pa.empty() ? "" : ",a1,a2,a3") << "\n";
            for (TPZGeoEl *gel : mT.GeoMesh()->ElementVec()) {
                if (!gel || gel->Dimension() != 2 || gel->NCornerNodes() != 3) continue;
                for (int i = 0; i < 3; i++)
                    fc << (i ? "," : "") << gel->NodePtr(i)->Coord(0) << "," << gel->NodePtr(i)->Coord(1) << ","
                       << pw[gel->NodeIndex(i)];
                for (int i = 0; i < 3 && !pa.empty(); i++) fc << "," << pa[gel->NodeIndex(i)];
                fc << "\n";
            }
        }
        const auto tf = Clock::now();
        std::vector<REAL> niveis;
        const REAL fs =
            op.soTransiente ? -1.
                            : Fator(r, loc, pw, pa, op.fs0, op.nref, niveis, grava ? pre + "_fator_" + tt + "h.csv" : "");
        const double sF = std::chrono::duration<double>(Clock::now() - tf).count();
        const double sT = std::chrono::duration<double>(Clock::now() - t0).count();
        std::string nv;
        for (REAL x : niveis) nv += (nv.empty() ? "" : ";") + std::to_string(x);
        const auto pr = chuva.PressureRange(), pra = chuva.AirPressureRange();
        const std::string cat = Categoria(nC, nM, nP);
        ff << t << "," << fs << "," << nv << "," << pr.first << "," << pr.second << "," << chuva.FaceNodes() << ","
           << chuva.Infiltrated() << "," << nC << "," << nM << "," << nP << "," << cat << "," << chuva.Steps() << ","
           << sTrans << "," << sF << "," << sT << "," << pra.first << "," << pra.second << "\n";
        ff.flush(), fn.flush();
        std::printf("  %5g %8.4f %8.2f %8.2f %8.2f %6zu %7.2f %7.2f %7.2f %4s %8.1f %8.1f\n", t, fs, pr.first, pr.second,
                    pra.second, chuva.FaceNodes(), nC, nM, nP, cat.c_str(), sF, sT);
        std::fflush(stdout);
    }
    std::printf("  CSV: %s_fs.csv, %s_nivel.csv; transiente %.1f s, total %.1f s\n", pre.c_str(), pre.c_str(), sTrans,
                std::chrono::duration<double>(Clock::now() - t0).count());
    return 0;
}

/// A TriGMesh do rahardjo.cpp (sem a gravacao do geo.vtk): talude de 6 cantos no sentido anti-horario; malha
/// estruturada com passo perto da altura H do talude: cada trecho horizontal e vertical em round(comprimento / H) partes
/// iguais, cada retangulo em dois triangulos e a cunha da face num triangulo so; ref = refinamentos uniformes
TPZGeoMesh *TriGMesh(const std::vector<std::vector<REAL>> &points, int ref) {
    if (points.size() != 6) DebugStop();
    const REAL x0 = points[0][0], y0 = points[0][1], L = points[1][0], yPe_ = points[2][1];
    const REAL xPe_ = points[3][0], xCr = points[4][0], yCr = points[4][1];
    if (points[1][1] != y0 || points[2][0] != L || points[3][1] != yPe_ || points[5][0] != x0 || points[5][1] != yCr ||
        !(x0 < xCr && xCr < xPe_ && xPe_ < L && y0 < yPe_ && yPe_ < yCr))
        DebugStop();
    const REAL H = yCr - yPe_;
    auto divide = [H](REAL a, REAL b) {
        const int n = std::max(1, int(std::lround((b - a) / H)));
        std::vector<REAL> v(n + 1);
        for (int k = 0; k <= n; k++) v[k] = a + (b - a) * k / n;
        return v;
    };
    std::vector<REAL> xs = divide(x0, xCr);
    const int iCr = int(xs.size()) - 1, iPe = iCr + 1;
    for (REAL x : divide(xPe_, L)) xs.push_back(x);
    const std::vector<REAL> ys = divide(y0, yPe_);
    const int nx = int(xs.size()) - 1, ny = int(ys.size()) - 1;
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
    auto triangulo = [gmesh](int64_t a, int64_t b, int64_t c) {
        TPZManVector<int64_t, 3> nodes = {a, b, c};
        new TPZGeoElRefPattern<pzgeom::TPZGeoTriangle>(nodes, matSolo, *gmesh);
    };
    for (int j = 0; j < ny; j++)
        for (int i = 0; i < nx; i++) {
            triangulo(no(i, j), no(i + 1, j), no(i, j + 1));
            triangulo(no(i + 1, j), no(i + 1, j + 1), no(i, j + 1));
        }
    for (int i = 0; i < iCr; i++) {
        triangulo(no(i, ny), no(i + 1, ny), topo(i));
        triangulo(no(i + 1, ny), topo(i + 1), topo(i));
    }
    triangulo(no(iCr, ny), no(iPe, ny), topo(iCr));
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
    return gmesh;
}
