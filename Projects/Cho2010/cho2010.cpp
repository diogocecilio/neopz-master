// Cho2010: o Exemplo 2 de Cho (2010), talude c-phi seco com campos aleatorios, no NeoPZ (Monte Carlo direto e
// amostragem por importancia).
//   S. E. Cho (2010). Probabilistic assessment of slope stability that considers the spatial variability of soil
//   properties. J. Geotech. Geoenviron. Eng. 136(7):975-984. doi:10.1061/(ASCE)GT.1943-5606.0000309
// O talude (Fig. 5): 1:1 com H = 10 m sobre 5 m de fundacao (15 m de altura total), 10 m de terreno antes do pe e 10 m de
// crista; aqui espelhado (crista a esquerda), como a TriGMesh do SlopeRandomFieldRefinement. Solo (Tabela 4): c' lognormal
// com media 10 kPa e CoV 0,3; phi' lognormal com media 30 graus e CoV 0,2; gamma = 20 kN/m3; correlacao cruzada r(c, phi)
// de -0,7 a 0,5 (base -0,5), nos campos gaussianos (Eq. 11). Correlacao exponencial separavel com l_h = 20 m e l_v = 2 m
// (Eq. 7). O artigo usa o equilibrio limite (Bishop com busca da superficie critica) e Monte Carlo com 50 000 realizacoes:
// FS deterministico 1,204; P_f = 0,0171 com r = -0,5 (Tabela 5) e cerca de 0,004 com r = -0,7 (Fig. 9).
// Aqui:
//  - elementos finitos com Mohr-Coulomb (TPZPlasticStepVoigt<TPZYCMohrCoulombPV2> dentro do TPZMatUnsaturatedUPwPa do
//    SlopeRandomFieldRefinement, com as pressoes congeladas em zero: solo seco, tensao efetiva = total); fluxo associado;
//    E = 20 MPa e nu = 0,3 (o equilibrio limite nao tem rigidez);
//  - o fator pelo SRM (c' e tan phi' divididos por F, como o FS do equilibrio limite; padrao) ou pelo GIM, com refinamento
//    adaptativo em h feito do zero em cada realizacao (malha grossa e nref niveis);
//  - c' e phi' por Karhunen-Loeve analitico (RandomField) no retangulo do dominio, com M termos;
//  - amostragem por importancia pelo piloto (ImportanceSampling.h): o deslocamento da media dos xi pela entropia cruzada.
// Uso: Cho2010 [chave=valor ...], com (padrao entre parenteses):
//   n (10), primeira (0)   realizacoes primeira..primeira+n-1
//   det (0)                1: antes, o deterministico (as medias), como a realizacao -1
//   valores ()             c,phi: o deterministico com esses valores em todo o solo (kPa, graus)
//   metodo (srm)           srm ou gim
//   covc (0.3), covphi (0.2), rcphi (-0.5)   CoV de c' e phi' e a correlacao cruzada dos campos gaussianos
//   lh (20), lv (2), M (2000)                distancias de autocorrelacao (m) e termos da KL
//   fator (tri1), nref (4) malha grossa (TriGMesh com K divisoes: lado H / 2^K) e niveis de refinamento
//   figuras ()             realizacoes com mecanismo.csv e campos.csv (pasta <saida>_real<n>), ex.: -1,0
//   saida (mc), semente (2010)
//   banda (<arquivo>)      mecanismo.csv de um deterministico: a direcao d dos xi e as projecoes proj_* no CSV
//   piloto (<pasta>), piloto_n (1000), piloto_semente (2010), ai (banda), elite (0.1), modos (30)
//                          amostragem por importancia pelo piloto (veja ImportanceSampling.h); sai com o codigo 3 se o
//                          piloto tiver menos de 90% das realizacoes. A amostra da AI precisa de outra semente= (semente).
// Saidas (<faixa> = primeira_ultima): <saida>_<faixa>.csv (uma linha por realizacao: FS_min = o fator, niveis, banda_*,
// massa_*, medias dos campos, log_peso, proj_*), <saida>_mapa_<faixa>.csv, <saida>_parametros.txt, <saida>_real<n>/ e,
// com piloto, <saida>_desloc_<faixa>.txt.
#include "HPMaterial.h"
#include "ImportanceSampling.h"
#include "RandomField.h"
#include "UnsatAnalysis.h"
#include "UnsatModel.h"

#include "Plasticity/TPZElasticResponse.h"
#include "Plasticity/TPZMatUnsaturatedUPwPa.h"
#include "Plasticity/TPZPlasticStepVoigt.h"
#include "Plasticity/TPZYCMohrCoulombPV2.h"
#include "TPZGeoLinear.h"
#include "TPZTensor.h"
#include "pzgeotriangle.h"
#include "pzmultiphysicselement.h"
#include "tpzgeoelrefpattern.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iomanip>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <vector>

using TMC = TPZPlasticStepVoigt<TPZYCMohrCoulombPV2, TPZElasticResponse>;
using TUnsat = TPZMatUnsaturatedUPwPa<TMC, TPZUnsaturatedMem>;
using TUnsatHP = HPMaterial<TUnsat>;
using Clock = std::chrono::steady_clock;
// propriedades num ponto (x, y): {c' (kPa), phi' (graus)}
using Campo = std::function<std::vector<double>(double, double)>;
static double Seg(Clock::time_point t0) { return std::chrono::duration<double>(Clock::now() - t0).count(); }

// ids da TriGMesh
enum EIds { matSolo = 1, base = -1, right = -2, topright = -3, topleft = -4, left = -5, face = -6 };

// o talude do Exemplo 2 (Fig. 5), espelhado: base esquerda, base direita, terreno do pe, pe, crista, crista a esquerda
const std::vector<std::vector<REAL>> kTalude = {{0, 0}, {30, 0}, {30, 5}, {20, 5}, {10, 15}, {0, 15}};

// ---- dados (Tabela 4 do Cho; E e nu so para os elementos finitos)
struct Solo {
    REAL c = 10., phiDeg = 30., gamma = 20., E = 2.e4, nu = 0.3;
};
struct Campos {
    REAL covC = 0.3, covPhi = 0.2, rcphi = -0.5, lh = 20., lv = 2.;
    int M = 2000;
    uint64_t semente = 2010;
};

struct Opcoes {
    int n = 10, primeira = 0, nref = 4;
    bool det = false, gim = false;
    Campos campos;
    std::string malhaF = "tri1", saida = "mc", banda;
    std::set<int> figuras;
    std::vector<REAL> valores; // deterministico com c' e phi' dados (vazio: as medias)
    ai::Opcoes ai;
};

std::vector<REAL> Lista(const std::string &s)
{
    std::vector<REAL> v;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, ','))
        if (!item.empty()) v.push_back(std::stod(item));
    return v;
}

Opcoes LeOpcoes(int argc, char *argv[])
{
    const std::set<std::string> chaves = {"n",     "primeira", "det",   "valores", "metodo", "covc",     "covphi",
                                          "rcphi", "lh",       "lv",    "M",       "fator",  "nref",     "figuras",
                                          "saida", "semente",  "banda", "piloto",  "piloto_n", "piloto_semente",
                                          "ai",    "elite",    "modos"};
    std::map<std::string, std::string> kv;
    for (int i = 1; i < argc; i++) {
        const std::string a = argv[i];
        const size_t eq = a.find('=');
        if (eq == std::string::npos || !chaves.count(a.substr(0, eq))) {
            std::fprintf(stderr, "opcao desconhecida: %s (as opcoes estao no cabecalho do cho2010.cpp)\n", a.c_str());
            std::exit(1);
        }
        kv[a.substr(0, eq)] = a.substr(eq + 1);
    }
    auto tem = [&kv](const char *k) { return kv.count(k) > 0; };
    Opcoes op;
    if (tem("n")) op.n = std::stoi(kv["n"]);
    if (tem("primeira")) op.primeira = std::stoi(kv["primeira"]);
    if (tem("det")) op.det = std::stoi(kv["det"]) != 0;
    if (tem("metodo")) {
        if (kv["metodo"] != "srm" && kv["metodo"] != "gim") {
            std::fprintf(stderr, "metodo: srm ou gim\n");
            std::exit(1);
        }
        op.gim = kv["metodo"] == "gim";
    }
    if (tem("covc")) op.campos.covC = std::stod(kv["covc"]);
    if (tem("covphi")) op.campos.covPhi = std::stod(kv["covphi"]);
    if (tem("rcphi")) op.campos.rcphi = std::stod(kv["rcphi"]);
    if (tem("lh")) op.campos.lh = std::stod(kv["lh"]);
    if (tem("lv")) op.campos.lv = std::stod(kv["lv"]);
    if (tem("M")) op.campos.M = std::stoi(kv["M"]);
    if (tem("semente")) op.campos.semente = std::stoull(kv["semente"]);
    if (tem("fator")) op.malhaF = kv["fator"];
    if (tem("nref")) op.nref = std::stoi(kv["nref"]);
    if (tem("figuras"))
        for (REAL v : Lista(kv["figuras"])) op.figuras.insert(int(v));
    if (tem("saida")) op.saida = kv["saida"];
    if (tem("banda")) op.banda = kv["banda"];
    op.ai.semente = tem("piloto_semente") ? std::stoull(kv["piloto_semente"]) : Campos{}.semente;
    if (tem("piloto")) op.ai.piloto = kv["piloto"];
    if (tem("piloto_n")) op.ai.n = std::stoi(kv["piloto_n"]);
    if (tem("ai")) op.ai.tipo = kv["ai"];
    if (tem("elite")) op.ai.elite = std::stod(kv["elite"]);
    if (tem("modos")) op.ai.modos = std::stoi(kv["modos"]);
    if (op.ai.tipo != "banda" && op.ai.tipo != "modos") {
        std::fprintf(stderr, "ai=%s: banda ou modos\n", op.ai.tipo.c_str());
        std::exit(1);
    }
    if (!op.ai.piloto.empty() && op.ai.tipo == "banda" && op.banda.empty()) {
        std::fprintf(stderr, "piloto com ai=banda precisa de banda=<mecanismo.csv>\n");
        std::exit(1);
    }
    if (!op.ai.piloto.empty() && op.ai.semente == op.campos.semente) {
        std::fprintf(stderr, "a amostragem por importancia precisa de outra semente= (a do piloto e %llu)\n",
                     (unsigned long long)op.ai.semente);
        std::exit(1);
    }
    if (!(std::abs(op.campos.rcphi) < 1.)) {
        std::fprintf(stderr, "rcphi precisa estar entre -1 e 1\n");
        std::exit(1);
    }
    if (tem("valores")) {
        op.valores = Lista(kv["valores"]);
        if (op.valores.size() != 2) {
            std::fprintf(stderr, "valores: c,phi (kPa, graus)\n");
            std::exit(1);
        }
    }
    if (op.malhaF.rfind("tri", 0) != 0) {
        std::fprintf(stderr, "malha desconhecida: %s (use triK)\n", op.malhaF.c_str());
        std::exit(1);
    }
    return op;
}

// ---- malha: a TriGMesh do SlopeRandomFieldRefinement (montecarlo.cpp), sem hierarquia
TPZGeoMesh *TriGMesh(const std::vector<std::vector<REAL>> &points, int ref)
{
    if (points.size() != 6) DebugStop();
    const REAL x0 = points[0][0], y0 = points[0][1], L = points[1][0], yPe = points[2][1];
    const REAL xPe = points[3][0], xCr = points[4][0], yCr = points[4][1];
    if (points[1][1] != y0 || points[2][0] != L || points[3][1] != yPe || points[5][0] != x0 || points[5][1] != yCr ||
        !(x0 < xCr && xCr < xPe && xPe < L && y0 < yPe && yPe < yCr))
        DebugStop();
    const REAL H = yCr - yPe;
    auto divide = [H](REAL a, REAL b) {
        const int n = std::max(1, int(std::lround((b - a) / H)));
        std::vector<REAL> v(n + 1);
        for (int k = 0; k <= n; k++) v[k] = a + (b - a) * k / n;
        return v;
    };
    std::vector<REAL> xs = divide(x0, xCr);
    const int iCr = int(xs.size()) - 1, iPe = iCr + 1;
    for (REAL x : divide(xPe, L)) xs.push_back(x);
    const std::vector<REAL> ys = divide(y0, yPe);
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

// Copia so com as folhas (o refinamento adaptativo parte de uma malha sem hierarquia)
TPZGeoMesh *MalhaPlana(TPZGeoMesh *gmesh)
{
    UnsatModel<TUnsatHP> tmp;
    tmp.SetGeoMesh(gmesh); // tmp fica dono da malha antiga e a apaga
    std::vector<int64_t> pai;
    return tmp.RefinedGeoMesh({}, pai);
}

TPZGeoMesh *Malha(const std::string &nome) { return MalhaPlana(TriGMesh(kTalude, std::stoi(nome.substr(3)))); }

// cota da superficie do talude em x
REAL Superficie(REAL x)
{
    const auto &p = kTalude;
    if (x <= p[4][0]) return p[4][1];
    if (x >= p[3][0]) return p[3][1];
    return p[4][1] + (p[3][1] - p[4][1]) * (x - p[4][0]) / (p[3][0] - p[4][0]);
}

// ---- material: Mohr-Coulomb no TPZMatUnsaturatedUPwPa com as pressoes congeladas em zero (solo seco). A lei
// hidraulica (Fredlund-Xing padrao) nao entra: sem succao, S_w = 1 e a tensao efetiva e a total
TUnsatHP *MaterialSeco()
{
    const Solo s;
    const REAL phi = s.phiDeg * M_PI / 180.;
    TPZElasticResponse ER;
    ER.SetEngineeringData(s.E, s.nu);
    auto *mat = new TUnsatHP(matSolo);
    mat->SetPlasticModel(TMC(TPZYCMohrCoulombPV2(phi, phi, s.c, ER), ER));
    TPZUnsaturatedHydraulics h;
    mat->SetHydraulics(h);
    mat->SetBishop(TUnsat::EChiConstant, 1.);
    mat->SetPorosity(0.3);
    mat->SetBiot(1., 0.);
    mat->SetWaterBulkModulus(2.2e6);
    mat->SetConductivity(1.e-7, 9.81, 1.e-3 / 1.8e-5);
    mat->SetAir(1.2 * 9.81e-3, 101.325);
    mat->SetUnitWeight(s.gamma);
    mat->SetLumpedStorage(true);
    mat->SetHydroMechanicalCoupling(false);
    return mat;
}

// c' e phi' (fmatprop) de cada ponto de integracao, do campo
void AplicaCampos(UnsatModel<TUnsatHP> &m, const Campo &f)
{
    m.ForEachIntegrationPoint([&f](const TPZVec<REAL> &x, TPZUnsaturatedMem &mem) {
        const std::vector<double> v = f(x[0], x[1]);
        const REAL phi = v[1] * M_PI / 180.;
        mem.m_elastoplastic_state.fmatprop = {v[0], phi, phi}; // {c, phi, psi}: fluxo associado
    });
}

// u: base engastada, laterais presas so em x
void InsertBC(UnsatModel<TUnsatHP> &m)
{
    m.AddDirichletBC(base, TUnsat::EDirichletU, 0, {0, 1});
    m.AddDirichletBC(left, TUnsat::EDirichletU, 0, {0});
    m.AddDirichletBC(right, TUnsat::EDirichletU, 0, {0});
}

// ---- fator de seguranca (o do montecarlo.cpp, sem agua)
struct Fator {
    std::unique_ptr<UnsatModel<TUnsatHP>> m;
    std::unique_ptr<UnsatAnalysis<TUnsatHP>> a; // declarada depois: destruida antes do modelo
};

void NovoFator(Fator &f, TPZGeoMesh *g, const Opcoes &op)
{
    f.a.reset();
    f.m = std::make_unique<UnsatModel<TUnsatHP>>();
    f.m->SetGeoMesh(g);
    InsertBC(*f.m);
    TUnsatHP *mat = MaterialSeco();
    if (op.gim) mat->SetGravityScalesSuction(false);
    f.m->SetMaterial(mat);
    f.m->CreateCompMesh({{2, 2}, {1, 1}, {1, 1}});
    f.a = std::make_unique<UnsatAnalysis<TUnsatHP>>(*f.m);
    f.a->SetVerbose(false);
    f.a->SetLoadFactorFunction([mat](REAL l) { mat->SetGravityFactor(l); }, false);
    f.a->SetStrengthReductionFunction([mat](REAL F) { mat->PlasticModel().SetStrengthReductionFactor(F); });
}

REAL CalculaFator(Fator &f, const Campo *campos, bool gim)
{
    f.m->Material()->ResetMemory(); // estado sem tensoes (apaga tambem c e phi dos pontos)
    if (campos) AplicaCampos(*f.m, *campos);
    f.a->SetFrozenPressures([](int, const TPZVec<REAL> &) { return 0.; }); // solo seco: p_w = p_a = 0
    return gim ? f.a->GIMArcLength() : f.a->SRMArcLength();
}

// Refinamento adaptativo em h do zero: malha grossa e nref niveis (marcados: sqrt(J2(eps_p)) >= 0,1 do maior)
REAL FatorAdaptativo(Fator &f, const Opcoes &op, const Campo *campos, std::vector<REAL> &niveis)
{
    TPZGeoMesh *g = Malha(op.malhaF);
    niveis.clear();
    for (int k = 0;; k++) {
        NovoFator(f, g, op);
        const REAL fs = CalculaFator(f, campos, op.gim);
        niveis.push_back(fs);
        if (k == op.nref) return fs;
        const std::vector<REAL> ind = f.a->PlasticIndicator();
        const REAL vmax = *std::max_element(ind.begin(), ind.end());
        std::set<int64_t> marcados;
        for (TPZGeoEl *gel : g->ElementVec())
            if (gel && gel->Dimension() == 2 && ind[gel->Index()] >= 0.1 * vmax) marcados.insert(gel->Index());
        std::vector<int64_t> pai;
        g = f.m->RefinedGeoMesh(marcados, pai);
    }
}

// ---- mecanismo no fim do caminho do fator: banda (sqrt(J2(eps_p)) >= 0,1 do maior), massa (|du| >= 0,5 do maior) e
// as medias dos campos
struct Mecanismo {
    REAL bArea = 0., bProf = 0., bXmin = 0., bXmax = 0., bXc = 0., bYc = 0., bYmin = 0., bC = 0., bPhi = 0.;
    REAL mArea = 0., mXc = 0., mYc = 0., mC = 0., mPhi = 0., cMed = 0., phiMed = 0.;
    struct Ponto {
        REAL x, y, j, du, c, phi;
    };
    std::vector<Ponto> pontos; // sqrt(J2) >= 1% do maior (mecanismo.csv)
    std::vector<std::array<REAL, 2>> banda;
};

Mecanismo DescreveMecanismo(Fator &f)
{
    const Solo res;
    UnsatModel<TUnsatHP> &m = *f.m;
    TUnsatHP *mat = m.Material();
    const int iu = mat->VariableIndex("Displacement");
    struct P {
        TPZCompEl *cel;
        TPZManVector<REAL, 3> qsi;
        REAL x, y, w, j, c, phi, du;
    };
    std::vector<P> ps;
    for (TPZCompEl *cel : m.CompMesh()->ElementVec()) {
        if (!cel || !cel->Reference() || cel->Reference()->MaterialId() != mat->Id()) continue;
        TPZGeoEl *gel = cel->Reference();
        auto *mf = dynamic_cast<TPZMultiphysicsElement *>(cel);
        TPZManVector<int, 3> ordens(mf->NMeshes(), 0);
        for (int i = 0; i < mf->NMeshes(); i++)
            if (auto *sp = dynamic_cast<TPZInterpolationSpace *>(mf->Element(i))) ordens[i] = sp->MaxOrder();
        const int o = mat->IntegrationRuleOrder(ordens);
        std::unique_ptr<TPZIntPoints> regra(gel->CreateSideIntegrationRule(gel->NSides() - 1, o));
        TPZManVector<int, 3> ordem(gel->Dimension(), o);
        regra->SetOrder(ordem);
        TPZManVector<int64_t> mem;
        cel->GetMemoryIndices(mem);
        for (int ip = 0; ip < regra->NPoints(); ip++) {
            P p;
            p.cel = cel;
            p.qsi.Resize(gel->Dimension());
            REAL w;
            regra->Point(ip, p.qsi, w);
            TPZManVector<REAL, 3> x(3, 0.);
            gel->X(p.qsi, x);
            TPZFNMatrix<9, REAL> jac(2, 2), axes(2, 3), jacinv(2, 2);
            REAL detjac;
            gel->Jacobian(p.qsi, jac, axes, detjac, jacinv);
            const TPZUnsaturatedMem &me = mat->MemItem(mem[ip]);
            TPZTensor<REAL> e(me.m_elastoplastic_state.m_eps_p); // cisalhamento de engenharia -> tensor
            e.XY() *= 0.5, e.XZ() *= 0.5, e.YZ() *= 0.5;
            const auto &pr = me.m_elastoplastic_state.fmatprop;
            const bool campo = pr.size() >= 2 && (pr[0] != 0. || pr[1] != 0.);
            p.x = x[0], p.y = x[1], p.w = w * std::abs(detjac), p.j = std::sqrt(std::max<REAL>(0., e.J2()));
            p.c = campo ? pr[0] : res.c, p.phi = campo ? pr[1] * 180. / M_PI : res.phiDeg, p.du = 0.;
            ps.push_back(p);
        }
    }
    auto deslocamentos = [&](std::vector<std::array<REAL, 2>> &u) {
        u.assign(ps.size(), {0., 0.});
        for (size_t i = 0; i < ps.size(); i++) {
            TPZManVector<STATE, 9> v(2, 0.);
            ps[i].cel->Solution(ps[i].qsi, iu, v);
            u[i] = {v[0], v[1]};
        }
    };
    std::vector<std::array<REAL, 2>> uf, ug(ps.size(), {0., 0.});
    deslocamentos(uf);
    if (f.a->WeightState().Rows() == f.a->PathState().Rows()) { // SRM: o deslocamento desde o estado geostatico
        f.a->LoadSolution(f.a->WeightState());
        deslocamentos(ug);
        f.a->LoadSolution(f.a->PathState());
    }
    REAL jmax = 0., dumax = 0., area = 0.;
    Mecanismo d;
    for (size_t i = 0; i < ps.size(); i++) {
        P &p = ps[i];
        p.du = std::hypot(uf[i][0] - ug[i][0], uf[i][1] - ug[i][1]);
        jmax = std::max(jmax, p.j), dumax = std::max(dumax, p.du);
        area += p.w;
        d.cMed += p.w * p.c, d.phiMed += p.w * p.phi;
    }
    d.cMed /= area, d.phiMed /= area;
    if (!(jmax > 0.)) return d;
    d.bXmin = d.bYmin = 1.e300, d.bXmax = -1.e300;
    for (const P &p : ps) {
        if (p.j >= 0.01 * jmax) d.pontos.push_back({p.x, p.y, p.j, p.du, p.c, p.phi});
        if (p.j >= 0.1 * jmax) {
            d.bArea += p.w;
            d.bXc += p.w * p.x, d.bYc += p.w * p.y;
            d.bC += p.w * p.c, d.bPhi += p.w * p.phi;
            d.bProf = std::max(d.bProf, Superficie(p.x) - p.y);
            d.bXmin = std::min(d.bXmin, p.x), d.bXmax = std::max(d.bXmax, p.x), d.bYmin = std::min(d.bYmin, p.y);
            d.banda.push_back({p.x, p.y});
        }
        if (dumax > 0. && p.du >= 0.5 * dumax) {
            d.mArea += p.w;
            d.mXc += p.w * p.x, d.mYc += p.w * p.y;
            d.mC += p.w * p.c, d.mPhi += p.w * p.phi;
        }
    }
    for (REAL *v : {&d.bXc, &d.bYc, &d.bC, &d.bPhi}) *v /= d.bArea;
    if (d.mArea > 0.)
        for (REAL *v : {&d.mXc, &d.mYc, &d.mC, &d.mPhi}) *v /= d.mArea;
    return d;
}

// ---- arquivos
const char *kCabecalho =
    "real,t_h,FS_min,situacao,niveis,eq_fator,s_fator,banda_area_m2,banda_prof,banda_xmin,banda_xmax,banda_xc,banda_yc,"
    "banda_ymin,banda_c,banda_phi_deg,massa_area_m2,massa_xc,massa_yc,massa_c,massa_phi_deg,c_med,phi_med,metodo,cov_c,"
    "cov_phi,r_cphi,lh,lv,KL_M,KL_variancia,log_peso,proj_c,proj_phi,proj_k\n";

struct Linha {
    int real = 0;
    REAL fs = 0.;
    std::string situacao = "ok", niveis;
    int64_t eqF = 0;
    double sF = 0., logPeso = 0.;
    std::vector<double> proj = {0., 0.};
    Mecanismo mec;
};

void Grava(std::ofstream &csv, const Linha &l, const Opcoes &op, const RandomField &campos)
{
    const Mecanismo &d = l.mec;
    const Campos &cp = op.campos;
    std::ostringstream s;
    s << l.real << ",0," << std::setprecision(7) << l.fs << "," << l.situacao << "," << l.niveis << "," << l.eqF << ","
      << std::setprecision(4) << l.sF << "," << std::setprecision(6) << d.bArea << "," << d.bProf << "," << d.bXmin << ","
      << d.bXmax << "," << d.bXc << "," << d.bYc << "," << d.bYmin << "," << d.bC << "," << d.bPhi << "," << d.mArea << ","
      << d.mXc << "," << d.mYc << "," << d.mC << "," << d.mPhi << "," << d.cMed << "," << d.phiMed << ","
      << (op.gim ? "gim" : "srm") << "," << cp.covC << "," << cp.covPhi << "," << cp.rcphi << "," << cp.lh << "," << cp.lv
      << "," << campos.M() << "," << campos.VarianceFraction() << "," << std::setprecision(10) << l.logPeso << ","
      << l.proj[0] << "," << (l.proj.size() > 1 ? l.proj[1] : 0.) << ",0\n";
    csv << s.str();
    csv.flush();
}

void GravaParametros(const std::string &arq, const Opcoes &op, const RandomField &campos)
{
    const Solo s;
    const Campos &cp = op.campos;
    std::ofstream p(arq);
    p << "EXEMPLO=Cho2010_2\nC=" << s.c << "\nPHI=" << s.phiDeg << "\nGAMMA=" << s.gamma << "\nE=" << s.E << "\nNU=" << s.nu
      << "\nCOV_C=" << cp.covC << "\nCOV_PHI=" << cp.covPhi << "\nR_CPHI=" << cp.rcphi << "\nLH=" << cp.lh << "\nLV="
      << cp.lv << "\nKL_M=" << campos.M() << "\nKL_VARIANCIA=" << campos.VarianceFraction() << "\nSEMENTE=" << cp.semente
      << "\nMETODO=" << (op.gim ? "gim" : "srm") << "\nMALHA_FATOR=" << op.malhaF << "\nNREF=" << op.nref
      << "\nAI=" << (op.ai.piloto.empty() ? "nenhuma" : op.ai.tipo + "_piloto") << "\n";
    if (!op.ai.piloto.empty())
        p << "PILOTO=" << op.ai.piloto << "\nPILOTO_N=" << op.ai.n << "\nPILOTO_SEMENTE=" << op.ai.semente
          << "\nELITE=" << op.ai.elite << "\nBANDA=" << op.banda << "\n";
}

// mapa da banda: celulas de 1 m no retangulo do talude, com as falhas
struct Mapa {
    static constexpr int nx = 30, ny = 15;
    std::vector<int> n = std::vector<int>(nx * ny, 0), nFalha = std::vector<int>(nx * ny, 0);
    int nReal = 0;
    void Soma(const Mecanismo &d, bool falha)
    {
        std::set<int> celulas;
        for (const auto &p : d.banda)
            celulas.insert(std::clamp(int(p[1]), 0, ny - 1) * nx + std::clamp(int(p[0]), 0, nx - 1));
        for (int c : celulas) n[c]++, nFalha[c] += falha;
    }
    void Grava(const std::string &arq) const
    {
        std::ofstream f(arq);
        f << "t_h,x,y,n,n_falha,n_real\n";
        for (int c = 0; c < nx * ny; c++)
            if (n[c]) f << "0," << c % nx + 0.5 << "," << c / nx + 0.5 << "," << n[c] << "," << nFalha[c] << "," << nReal << "\n";
    }
};

int main(int argc, char *argv[])
{
    const Opcoes op = LeOpcoes(argc, argv);
    const Solo res;
    const Campos &cp = op.campos;
    const auto &p = kTalude;
    RandomField campos(p[0][0], p[0][1], p[1][0] - p[0][0], p[4][1] - p[0][1], cp.lh, cp.lv, cp.M,
                       {{res.c, cp.covC}, {res.phiDeg, cp.covPhi}}); // c' (kPa), phi' (graus)
    campos.SetCrossCorrelation(0, 1, cp.rcphi);
    const Campo campoAleatorio = [&campos](double x, double y) { return campos.At(x, y); };
    const Campo campoFixo = [v = op.valores](double, double) { return std::vector<double>(v.begin(), v.end()); };
    std::printf("Cho (2010), Exemplo 2: talude 1:1, H = 10 m, seco | %s com %s + %d niveis | realizacoes %d..%d%s\n"
                "  campos: c' (10 kPa, CoV %g), phi' (30 graus, CoV %g), r(c, phi) = %g; l_h = %g m, l_v = %g m; KL: %d "
                "termos, %.1f%% da variancia; semente %llu\n",
                op.gim ? "GIM" : "SRM", op.malhaF.c_str(), op.nref, op.primeira, op.primeira + op.n - 1,
                op.det ? " e o deterministico" : "", cp.covC, cp.covPhi, cp.rcphi, cp.lh, cp.lv, campos.M(),
                100. * campos.VarianceFraction(), (unsigned long long)cp.semente);

    const std::string faixa = std::to_string(op.primeira) + "_" + std::to_string(op.primeira + op.n - 1);
    std::vector<double> dir; // amostragem por importancia: a direcao da banda e o deslocamento pelo piloto
    if (!op.banda.empty()) {
        dir = campos.Direction(ai::LeBanda(op.banda));
        if (dir.empty()) {
            std::fprintf(stderr, "banda sem pontos: %s\n", op.banda.c_str());
            return 1;
        }
    }
    if (!op.ai.piloto.empty() && !ai::DeslocaPeloPiloto(campos, dir, op.ai, op.saida + "_desloc_" + faixa + ".txt"))
        return 3;
    const std::string arq = op.saida + "_" + faixa + ".csv";
    std::ofstream csv(arq);
    csv << kCabecalho;
    if (op.primeira == 0) GravaParametros(op.saida + "_parametros.txt", op, campos);
    Mapa mapa;
    std::vector<REAL> fss;
    double somaPesoFalha = 0., somaF = 0.;

    std::vector<int> realizacoes;
    if (op.det) realizacoes.push_back(-1);
    for (int n = op.primeira; n < op.primeira + op.n; n++) realizacoes.push_back(n);
    for (int n : realizacoes) {
        const auto t0 = Clock::now();
        const bool deterministico = n < 0;
        if (!deterministico) campos.Sample(cp.semente, n);
        Linha l;
        l.real = n;
        if (!deterministico) l.logPeso = campos.LogWeight(), l.proj = campos.Projection(dir);
        if (l.proj.size() < 2) l.proj.resize(2, 0.);
        const Campo *cf = deterministico ? (op.valores.empty() ? nullptr : &campoFixo) : &campoAleatorio;
        Fator f;
        std::vector<REAL> niveis;
        l.fs = FatorAdaptativo(f, op, cf, niveis);
        l.sF = Seg(t0), l.eqF = f.m->CompMesh()->NEquations();
        if (!(l.fs > 0.)) l.situacao = "fator";
        for (size_t i = 0; i < niveis.size(); i++) l.niveis += (i ? ";" : "") + std::to_string(niveis[i]).substr(0, 7);
        l.mec = DescreveMecanismo(f);
        if (op.figuras.count(n)) { // o mecanismo e os campos num reticulado de 0,5 m
            const std::string pasta = op.saida + "_real" + std::to_string(n);
            mkdir(pasta.c_str(), 0775);
            std::ofstream mf(pasta + "/mecanismo.csv");
            mf << "t_h,x,y,epsp_raizJ2,du,c,phi_deg,k_mult\n";
            for (const auto &q : l.mec.pontos)
                mf << "0," << q.x << "," << q.y << "," << q.j << "," << q.du << "," << q.c << "," << q.phi << ",1\n";
            std::ofstream cfs(pasta + "/campos.csv");
            cfs << "x,y,c,phi\n";
            for (REAL y = 0.25; y < p[4][1]; y += 0.5)
                for (REAL x = 0.25; x < p[1][0]; x += 0.5) {
                    if (y > Superficie(x)) continue;
                    const std::vector<double> v = cf ? (*cf)(x, y) : std::vector<double>{res.c, res.phiDeg};
                    cfs << x << "," << y << "," << v[0] << "," << v[1] << "\n";
                }
        }
        Grava(csv, l, op, campos);
        if (!deterministico && l.situacao == "ok") {
            fss.push_back(l.fs);
            mapa.Soma(l.mec, l.fs < 1.);
            mapa.nReal++, somaF += l.sF;
            somaPesoFalha += (l.fs < 1.) * std::exp(l.logPeso);
        }
        std::printf("  realizacao %4d: FS = %.4f (%s) | niveis %s | %.1f s\n", n, l.fs, l.situacao.c_str(),
                    l.niveis.c_str(), Seg(t0));
        std::fflush(stdout);
    }
    mapa.Grava(op.saida + "_mapa_" + faixa + ".csv");
    if (!fss.empty()) {
        REAL m = 0., v = 0.;
        for (REAL x : fss) m += x;
        m /= fss.size();
        for (REAL x : fss) v += (x - m) * (x - m);
        const int falhas = std::count_if(fss.begin(), fss.end(), [](REAL x) { return x < 1.; });
        std::printf("\n  %zu realizacoes: FS medio %.4f, desvio %.4f, P(FS < 1) contada = %.2f%%", fss.size(), m,
                    fss.size() > 1 ? std::sqrt(v / (fss.size() - 1)) : 0., 100. * falhas / fss.size());
        if (!op.ai.piloto.empty()) std::printf(", P_f ponderada = %.3f%%", 100. * somaPesoFalha / fss.size());
        std::printf("\n  tempo medio por realizacao: %.1f s\n", somaF / mapa.nReal);
    }
    std::printf("  arquivos: %s e %s\n", arq.c_str(), (op.saida + "_mapa_" + faixa + ".csv").c_str());
    return 0;
}
