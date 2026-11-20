// Caso deterministico do rebaixamento rapido do Ceron (2025), Secao 6.1, com o fator Gamma (GIM) calculado com o
// material SATURADO u-p (TPZMatPoroElastoPlasticUP) no lugar do nao saturado u-p_w-p_a (TPZMatUnsaturatedUPwPa).
// E o caso det=1 do rebaixamento_mc.cpp (SlopeDrawdownMonteCarlo n=0 det=1) com uma so troca:
//  1. o transiente e o MESMO: so a agua, u-p_w-p_a (van Genuchten-Mualem), esqueleto rigido, ar passivo, na malha
//     transiente=triK, o reservatorio descendo da crista ao pe em trebaix horas (passo dt);
//  2. o fator: as pressoes do fim do transiente congeladas (PointLocator + SetFrozenPressures), o peso e a agua do
//     reservatorio (normal a face e ao terreno do pe) crescendo com lambda, GIM por arc-length, refinamento adaptativo
//     em h do zero (malha fator=triK e nref niveis, marcados os triangulos com sqrt(J2(eps_p)) >= 0,1 do maior). Tudo
//     igual, mas o material do fator e o u-p (Mohr-Coulomb associado, c' = 10 kPa, phi' = 30 graus, gamma = 20,
//     E = 20 MPa, nu = 0,3), com sigma = sigma' - p (alfa = 1).
// O u-p nao tem chi: a pressao p que entra em sigma' = sigma + p e a pressao congelada. A opcao succao diz o que fazer
// com as pressoes negativas (succao) acima do NA:
//   zero   p = max(p_w, 0): a succao e desprezada (a analise saturada usual)
//   total  p = p_w: a succao entra inteira (chi = 1)
//   sw     p = S_w(p_w) p_w nos nos: o chi = S_w do u-p_w-p_a, aplicado nos nos (o u-p_w-p_a aplica nos pontos de
//          integracao, com o p_w interpolado); serve de conferencia do u-p contra o u-p_w-p_a
// Uso: SlopeDrawdownUP [chave=valor ...] (padrao entre parenteses):
//   nref (3), fator (tri1), transiente (tri3), trebaix (0.1), dt (0.1, h), hw (1)
//   tfim (trebaix)  instante (h) do Gamma; agua (unsat): sat = solo sempre saturado (S_w = 1, k_r = 1), como no
//                   rebaixamento_mc.cpp
//   succao (zero,total,sw)  lista: um fator por item, todos com o mesmo transiente
//   saida (up)               <saida>.csv: uma linha por (succao, nivel) com Gamma, equacoes e tempo
#include "DrawdownTransient.h"
#include "HPMaterial.h"
#include "PointLocator.h"
#include "UnsatAnalysis.h"
#include "UnsatModel.h"

#include "Plasticity/TPZElasticResponse.h"
#include "Plasticity/TPZMatPoroElastoPlasticUP.h"
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
#include <iomanip>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using TMC = TPZPlasticStepVoigt<TPZYCMohrCoulombPV2, TPZElasticResponse>;
using TUnsat = TPZMatUnsaturatedUPwPa<TMC, TPZUnsaturatedMem>;
using TUnsatHP = HPMaterial<TUnsat>;
using TUP = TPZMatPoroElastoPlasticUP<TMC, TPZElastoPlasticMem>;
using TUPHP = HPMaterial<TUP>;
using Clock = std::chrono::steady_clock;
static double Seg(Clock::time_point t0) { return std::chrono::duration<double>(Clock::now() - t0).count(); }

enum EIds { matSolo = 1, base = -1, right = -2, topright = -3, topleft = -4, left = -5, face = -6 };
const std::vector<std::vector<REAL>> kTalude = {{0, 0}, {70, 0}, {70, 30}, {40, 30}, {35, 35}, {0, 35}};

struct Solo { // o Solo do rebaixamento_mc.cpp
    REAL c = 10., phiDeg = 30.;
    REAL gamma = 20., gammaw = 9.81, E = 2.e4, nu = 0.3;
    REAL kv = 1.e-6, n = 0.45, Kw = 2.2e6;
};

struct Opcoes {
    int nref = 3;
    std::string malhaT = "tri3", malhaF = "tri1", saida = "up";
    REAL trebaix = 0.1, dt = 0.1, hw = 1., tfim = -1.;
    std::string agua = "unsat";
    std::vector<std::string> succao = {"zero", "total", "sw"};
};

Opcoes LeOpcoes(int argc, char *argv[])
{
    const std::set<std::string> chaves = {"nref", "fator", "transiente", "trebaix", "dt", "hw", "succao", "saida", "tfim", "agua"};
    std::map<std::string, std::string> kv;
    for (int i = 1; i < argc; i++) {
        const std::string a = argv[i];
        const size_t eq = a.find('=');
        if (eq == std::string::npos || !chaves.count(a.substr(0, eq))) {
            std::fprintf(stderr, "opcao desconhecida: %s (as opcoes estao no cabecalho do rebaixamento_up.cpp)\n", a.c_str());
            std::exit(1);
        }
        kv[a.substr(0, eq)] = a.substr(eq + 1);
    }
    Opcoes op;
    if (kv.count("nref")) op.nref = std::stoi(kv["nref"]);
    if (kv.count("fator")) op.malhaF = kv["fator"];
    if (kv.count("transiente")) op.malhaT = kv["transiente"];
    if (kv.count("trebaix")) op.trebaix = std::stod(kv["trebaix"]);
    if (kv.count("dt")) op.dt = std::stod(kv["dt"]);
    if (kv.count("tfim")) op.tfim = std::stod(kv["tfim"]);
    if (kv.count("agua")) op.agua = kv["agua"];
    if (kv.count("hw")) op.hw = std::stod(kv["hw"]);
    if (kv.count("saida")) op.saida = kv["saida"];
    if (kv.count("succao")) {
        op.succao.clear();
        std::stringstream ss(kv["succao"]);
        for (std::string s; std::getline(ss, s, ',');)
            if (!s.empty()) op.succao.push_back(s);
    }
    for (const std::string &s : op.succao)
        if (s != "zero" && s != "total" && s != "sw") {
            std::fprintf(stderr, "succao=%s: zero, total ou sw\n", s.c_str());
            std::exit(1);
        }
    for (const std::string *m : {&op.malhaT, &op.malhaF})
        if (m->rfind("tri", 0) != 0) {
            std::fprintf(stderr, "malha desconhecida: %s (use triK)\n", m->c_str());
            std::exit(1);
        }
    return op;
}

// ---- malha: a TriGMesh do rebaixamento_mc.cpp
TPZGeoMesh *TriGMesh(const std::vector<std::vector<REAL>> &points, int ref)
{
    const REAL x0 = points[0][0], y0 = points[0][1], L = points[1][0], yPe = points[2][1];
    const REAL xPe = points[3][0], xCr = points[4][0], yCr = points[4][1];
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

TPZGeoMesh *Malha(const std::string &nome)
{
    UnsatModel<TUnsatHP> tmp;
    tmp.SetGeoMesh(TriGMesh(kTalude, std::stoi(nome.substr(3))));
    std::vector<int64_t> pai;
    return tmp.RefinedGeoMesh({}, pai);
}

static bool gAguaSaturada = false; // agua=sat

TPZUnsaturatedHydraulics Hidraulica()
{
    TPZUnsaturatedHydraulics h;
    h.fRetention = TPZUnsaturatedHydraulics::EVanGenuchten;
    h.fKrw = TPZUnsaturatedHydraulics::EMualem;
    h.fKra = TPZUnsaturatedHydraulics::EMualemAir;
    if (gAguaSaturada) h.fAVG = 1.e-12; // S_w = 1 e k_r = 1
    return h;
}

// o material do transiente: o MaterialCeron(1, false) do rebaixamento_mc.cpp
TUnsatHP *MaterialTransiente()
{
    const Solo s;
    const REAL phi = s.phiDeg * M_PI / 180.;
    TPZElasticResponse ER;
    ER.SetEngineeringData(s.E, s.nu);
    auto *mat = new TUnsatHP(matSolo);
    mat->SetPlasticModel(TMC(TPZYCMohrCoulombPV2(phi, phi, s.c, ER), ER));
    mat->SetHydraulics(Hidraulica());
    mat->SetPorosity(s.n);
    mat->SetBiot(1., 0.);
    mat->SetWaterBulkModulus(s.Kw);
    mat->SetConductivity(s.kv, s.gammaw, 1.e-3 / 1.8e-5);
    mat->SetAir(1.2 * 9.81e-3, 101.325);
    mat->SetUnitWeight(s.gamma);
    mat->SetGravityDirection(TPZManVector<REAL, 3>({0., -1., 0.}));
    mat->SetBishop(TUnsat::EChiSaturation, 1.);
    mat->SetLumpedStorage(true);
    mat->SetHydroMechanicalCoupling(false);
    return mat;
}

// ---- fator u-p: o NovoFator do rebaixamento_mc.cpp com o material saturado
struct FatorUP {
    std::shared_ptr<std::array<REAL, 2>> estado; // {lambda, nivel do reservatorio}
    std::unique_ptr<UnsatModel<TUPHP>> m;
    std::unique_ptr<UnsatAnalysis<TUPHP>> a;
};

void NovoFatorUP(FatorUP &f, TPZGeoMesh *g, REAL nivel)
{
    f.a.reset();
    f.estado = std::make_shared<std::array<REAL, 2>>(std::array<REAL, 2>{1., nivel});
    f.m = std::make_unique<UnsatModel<TUPHP>>();
    f.m->SetGeoMesh(g);
    // u: base engastada, laterais presas so em x. No u-p a BC de Dirichlet do material e por penalidade: a direcional
    // prende so as componentes dadas (EDirichletU prenderia as duas nas laterais)
    f.m->AddDirichletBC(base, TUP::EDirichletUDirectional, 0, {0, 1});
    f.m->AddDirichletBC(left, TUP::EDirichletUDirectional, 0, {0});
    f.m->AddDirichletBC(right, TUP::EDirichletUDirectional, 0, {0});
    const Solo s;
    auto pressao = [est = f.estado, gw = s.gammaw](REAL nx, REAL ny) {
        return [est, gw, nx, ny](const TPZVec<REAL> &x, TPZVec<STATE> &v2, TPZFMatrix<STATE> &) {
            const REAL pr = (*est)[0] * gw * std::max(REAL(0), (*est)[1] - x[1]);
            v2[0] = -pr * nx, v2[1] = -pr * ny;
        };
    };
    const auto &p = kTalude;
    const REAL dx = p[3][0] - p[4][0], dy = p[4][1] - p[3][1], L = std::hypot(dx, dy);
    // ENeumannUFixed: a tracao nao e multiplicada pelo fator de carga do material; o lambda ja esta na funcao
    f.m->AddNaturalBC(face, TUP::ENeumannUFixed, {0., 0.}, pressao(dy / L, dx / L));
    f.m->AddNaturalBC(topright, TUP::ENeumannUFixed, {0., 0.}, pressao(0., 1.));
    const REAL phi = s.phiDeg * M_PI / 180.;
    TPZElasticResponse ER;
    ER.SetEngineeringData(s.E, s.nu);
    auto *mat = new TUPHP(matSolo, TPZMatPoroElastoPlasticUPBase::EPlaneStrain);
    mat->SetPlasticModel(TMC(TPZYCMohrCoulombPV2(phi, phi, s.c, ER), ER));
    mat->SetBiot(1., 0.);
    mat->SetPermeability(s.kv);
    mat->SetBodyForce(TPZManVector<REAL, 3>({0., -s.gamma, 0.}));
    f.m->SetMaterial(mat);
    f.m->CreateCompMesh({{2, 2}, {1, 1}});
    f.a = std::make_unique<UnsatAnalysis<TUPHP>>(*f.m);
    f.a->SetVerbose(false);
    auto est = f.estado;
    const REAL gamma = s.gamma;
    // peso lambda gamma; as pressoes congeladas a analise multiplica por lambda (scaleFrozen = true), como o
    // SetGravityFactor do u-p_w-p_a faz com p_B
    f.a->SetLoadFactorFunction(
        [mat, est, gamma](REAL l) { mat->SetBodyForce(TPZManVector<REAL, 3>({0., -l * gamma, 0.})), (*est)[0] = l; }, true);
    f.a->SetStrengthReductionFunction([mat](REAL F) { mat->PlasticModel().SetStrengthReductionFactor(F); });
}

int main(int argc, char *argv[])
{
    const Opcoes op = LeOpcoes(argc, argv);
    gAguaSaturada = op.agua == "sat";
    const Solo res;
    const auto &p = kTalude;
    const REAL H = p[4][1] - p[3][1], nivelFinal = p[4][1] - op.hw * H;
    std::printf("Rebaixamento do Ceron, deterministico, Gamma (GIM) com o u-p | transiente %s (dt = %g h, %g h), fator %s + "
                "%d niveis, h_w/H = %g\n",
                op.malhaT.c_str(), op.dt, op.trebaix, op.malhaF.c_str(), op.nref, op.hw);

    // 1. transiente (o mesmo do rebaixamento_mc.cpp)
    const auto t0 = Clock::now();
    UnsatModel<TUnsatHP> mT;
    mT.SetGeoMesh(Malha(op.malhaT));
    mT.AddDirichletBC(base, TUnsat::EDirichletU, 0, {0, 1});
    mT.AddDirichletBC(left, TUnsat::EDirichletU, 0, {0});
    mT.AddDirichletBC(right, TUnsat::EDirichletU, 0, {0});
    mT.SetMaterial(MaterialTransiente());
    mT.CreateCompMesh({{2, 2}, {1, 1}, {1, 1}});
    TPZUnsatAnalysis aT(mT.CompMesh(), {{base, {true, true}}, {right, {true, false}}, {left, {true, false}}});
    TPZSpStructMatrix<STATE> strmat(mT.CompMesh());
    strmat.SetNumThreads(std::getenv("PZ_THREADS") ? std::atoi(std::getenv("PZ_THREADS")) : 8);
    aT.SetStructuralMatrix(strmat);
    TPZStepSolver<STATE> solver;
    solver.SetDirect(ELU);
    aT.SetSolver(solver);
    DrawdownData d;
    d.gammaw = res.gammaw;
    d.nivelInicial = p[4][1];
    d.nivelFinal = nivelFinal;
    d.nivelMacico = p[4][1];
    d.tRebaixamento = op.trebaix;
    d.dt = op.dt;
    DrawdownTransient rebaixamento(aT, d, left, {right, topright, face, topleft});
    rebaixamento.InitialState();
    const PointLocator loc(mT.GeoMesh());
    if (!rebaixamento.AdvanceTo(std::max(op.tfim, op.trebaix))) {
        std::printf("o transiente parou em %.3f h\n", rebaixamento.Time());
        return 2;
    }
    const REAL nivel = rebaixamento.Level();
    const std::vector<REAL> pw = rebaixamento.NodalWaterPressure();
    const double sT = Seg(t0);
    std::printf("  transiente: %lld equacoes, %.1f s, p_w de %.3f a %.3f kPa, nivel %.3f m\n",
                (long long)mT.CompMesh()->NEquations(), sT, *std::min_element(pw.begin(), pw.end()),
                *std::max_element(pw.begin(), pw.end()), nivel);
    const TPZUnsaturatedHydraulics hidro = Hidraulica();

    // 2. Gamma com o u-p, um caminho adaptativo por tratamento da succao
    std::ofstream csv(op.saida + ".csv");
    csv << "succao,nivel,Gamma,equacoes,triangulos,s_nivel,s_acumulado,transiente,dt,eq_transiente,s_transiente\n";
    for (const std::string &suc : op.succao) {
        auto pFator = [&](const TPZVec<REAL> &x) {
            const REAL v = loc.Interpolate(x[0], x[1], pw);
            if (suc == "zero") return std::max(REAL(0), v);
            if (suc == "total") return v;
            return hidro.Evaluate(v, 0.).S * v; // chi = S_w, p_a = 0
        };
        TPZGeoMesh *g = Malha(op.malhaF);
        FatorUP f;
        const auto tc = Clock::now();
        std::printf("  succao = %s\n", suc.c_str());
        for (int k = 0;; k++) {
            const auto tk = Clock::now();
            NovoFatorUP(f, g, nivel);
            f.a->SetFrozenPressures([&](int espaco, const TPZVec<REAL> &x) { return espaco == 1 ? pFator(x) : 0.; });
            const REAL gam = f.a->GIMArcLength();
            int64_t ntri = 0;
            for (TPZGeoEl *gel : g->ElementVec())
                if (gel && gel->Dimension() == 2) ntri++;
            const double sk = Seg(tk), sc = Seg(tc);
            const int64_t neq = f.m->CompMesh()->NEquations();
            std::printf("    nivel %d: Gamma = %.5f, %lld equacoes, %lld triangulos, %.1f s (acumulado %.1f s)\n", k, gam,
                        (long long)neq, (long long)ntri, sk, sc);
            std::fflush(stdout);
            csv << suc << "," << k << "," << std::setprecision(7) << gam << "," << neq << "," << ntri << ","
                << std::setprecision(4) << sk << "," << sc << "," << op.malhaT << "," << op.dt << ","
                << mT.CompMesh()->NEquations() << "," << sT << "\n";
            csv.flush();
            if (k == op.nref) break;
            const std::vector<REAL> ind = f.a->PlasticIndicator();
            const REAL vmax = *std::max_element(ind.begin(), ind.end());
            std::set<int64_t> marcados;
            for (TPZGeoEl *gel : g->ElementVec())
                if (gel && gel->Dimension() == 2 && ind[gel->Index()] >= 0.1 * vmax) marcados.insert(gel->Index());
            std::vector<int64_t> pai;
            g = f.m->RefinedGeoMesh(marcados, pai);
        }
    }
    std::printf("  arquivo: %s.csv\n", op.saida.c_str());
    return 0;
}
