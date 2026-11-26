// Cho2016: a reproducao do exemplo de S. E. Cho (2016), "Stability analysis of unsaturated soil slopes considering
// water-air flow caused by rainfall infiltration", Engineering Geology 211:184-197 (manuscrito aceito).
// O exemplo (Secao 4): o talude da Fig. 1 (A(0,28), B(20,26), C(28,22), D(40,12), E(60,10), F(60,0), G(0,0)), nivel
// d'agua horizontal a 6 m, chuva constante de 3 k_s na superficie ABCDE, dois solos (areia e silte, Tabelas 1 e 2), o
// fluxo so da agua (sp, ar passivo) e o bifasico agua-ar (tp, ar ativo) no FLAC acoplado, e o FS pelo SRM ao longo do
// tempo (Fig. 12).
//
// Uso: Cho2016 [solo=areia|silte] [modelo=sp|tp] [dt=<h>] [tempos=t1,t2,...] [nref=2] [tolfs=<dFS>] [malhaF=media|grossa] [psi=<graus>] [saida=cho] [so_transiente=0|1] [grava_malha=<prefixo>]
//
// As pecas (todas ja validadas no SlopeRandomFieldRefinement, so reaproveitadas):
//   1. a malha em camadas do LayerGMesh com a superficie ABCDE (um id por trecho), finas perto da superficie;
//   2. o material u - p_w - p_a da biblioteca: van Genuchten com S_r = 0 e Mualem da agua e do ar (b = 0,5), chi = S_w,
//      esqueleto elastico (E e nu de K e G da Tabela 2), peso da mistura (gamma_d + n S_w gamma_w), acoplado;
//   3. o transiente da chuva (RainTransient): nivel fixo na lateral direita (EF), as outras fronteiras sem fluxo, na
//      superficie a chuva ou a face de percolacao (p_w = 0, a condicao do artigo, Secao 4.2) e, no tp, p_a = 0;
//   4. o fator pelo SRM com as pressoes do transiente congeladas (as rotinas NovoFator, CalculaFator e FatorAdaptativo
//      do montecarlo.cpp, copiadas), com refinamento adaptativo; apoios: AG e EF sem deslocamento normal, FG fixo.

#include "HPMaterial.h"
#include "LayerGMesh.h"
#include "MualemTabelado.h"
#include "PointLocator.h"
#include "RainTransient.h"
#include "UnsatAnalysis.h"
#include "UnsatModel.h"

#include "Plasticity/TPZElasticResponse.h"
#include "Plasticity/TPZMatUnsaturatedUPwPa.h"
#include "Plasticity/TPZPlasticStepPV.h"
#include "Plasticity/TPZPlasticStepVoigt.h"
#include "Plasticity/TPZYCMohrCoulombPV.h"
#include "Plasticity/TPZYCMohrCoulombPV2.h"

#include <algorithm>
#include <type_traits>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using TMC = TPZPlasticStepVoigt<TPZYCMohrCoulombPV2, TPZElasticResponse>;
using TUnsat = TPZMatUnsaturatedUPwPa<TMC, TPZUnsaturatedMem>;
using TUnsatHP = HPMaterial<TUnsat>;
// o fator com dilatancia psi != phi (opcao psi=): o Mohr-Coulomb PV com o TPZPlasticStepPV (o PV2 so e associado)
using TMCPV = TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>;
using TUnsatPV = HPMaterial<TPZMatUnsaturatedUPwPa<TMCPV, TPZUnsaturatedMem>>;

enum EIds { matSolo = 1, base = -1, right = -2, topright = -3, topleft = -4, left = -5, face = -6 };
constexpr REAL gammaW = 9.81;

// ---- os dados do artigo (Tabelas 1 e 2)
struct SoloCho {
    std::string nome;
    REAL ks;    // m/s
    REAL P0;    // kPa (van Genuchten: S_e = [1 + (psi/P0)^(1/(1-a))]^(-a))
    REAL a;     // o m do van Genuchten
    REAL n;     // porosidade
    REAL c;     // kPa
    REAL phi;   // graus
    REAL tfim;  // h: o fim da Fig. 12
    std::vector<REAL> tempos; // os instantes do FS (h)
};
static SoloCho Areia() {
    SoloCho s{"areia", 1.e-5, 1.5, 0.445, 0.35, 3., 35., 30., {}};
    for (int t = 0; t <= 30; t += 2) s.tempos.push_back(t);
    return s;
}
static SoloCho Silte() {
    SoloCho s{"silte", 3.e-6, 4.905, 0.336, 0.5, 10., 30., 168., {}};
    for (int t = 0; t <= 168; t += 8) s.tempos.push_back(t);
    return s;
}
// comuns: K = 33,3 MPa e G = 11,1 MPa (E = 9KG/(3K+G), nu = (3K-2G)/(2(3K+G))); rho_d = 1640 kg/m3; mu_w/mu_a = 56;
// K_w = 2 GPa; rho_a = 1,25 kg/m3; o ar como gas ideal a 101,325 kPa (K_a ~ 100 kPa); dilatancia 0 no artigo (aqui o
// SRM usa o fluxo associado, como os outros drivers do projeto)
constexpr REAL Kb = 33.3e3, Gb = 11.1e3, gammaD = 1640. * 9.81e-3;

struct Opcoes {
    SoloCho solo = Areia();
    bool tp = false;
    REAL dt = -1.;      // h; < 0: 0,05 na areia e 0,25 no silte
    int nref = 2;
    std::string saida = "cho";
    bool soTransiente = false;
    REAL eta0 = 0.002, dx = 1., rT = 1.12, etamaxT = 0.06; // rT, etamaxT: crescimento das camadas do transiente
    REAL dtMon = 1.;    // h: o intervalo do registro nos pontos H, I e J // a malha do transiente: 1a camada (fracao da altura) e passo das colunas (m)
    REAL tolfs = -1.;   // > 0: refina ate |FS_k - FS_(k-1)| < tolfs, com no maximo nref niveis
    std::string malhaF = "media"; // malha inicial do fator: media (dx 2 m) ou grossa (dx 4 m)
    std::string gravaMalha; // grava_malha=<prefixo>: os triangulos das malhas (transiente e inicial do fator) e sai
    REAL psi = -1.;     // dilatancia no fator (graus); < 0: associado (psi = phi)
};

static std::vector<REAL> Lista(const std::string &s) {
    std::vector<REAL> v;
    std::stringstream ss(s);
    for (std::string t; std::getline(ss, t, ',');) v.push_back(std::stod(t));
    return v;
}

static Opcoes LeOpcoes(int argc, char **argv) {
    Opcoes op;
    std::vector<REAL> tempos;
    for (int i = 1; i < argc; i++) {
        const std::string s = argv[i];
        const size_t eq = s.find('=');
        const std::string k = s.substr(0, eq), v = eq == std::string::npos ? "" : s.substr(eq + 1);
        if (k == "solo" && (v == "areia" || v == "silte")) op.solo = v == "areia" ? Areia() : Silte();
        else if (k == "modelo" && (v == "sp" || v == "tp")) op.tp = v == "tp";
        else if (k == "dt") op.dt = std::stod(v);
        else if (k == "tempos") tempos = Lista(v);
        else if (k == "nref") op.nref = std::stoi(v);
        else if (k == "saida") op.saida = v;
        else if (k == "eta0") op.eta0 = std::stod(v);
        else if (k == "dx") op.dx = std::stod(v);
        else if (k == "rT") op.rT = std::stod(v);
        else if (k == "etamaxT") op.etamaxT = std::stod(v);
        else if (k == "dt_mon") op.dtMon = std::stod(v);
        else if (k == "tolfs") op.tolfs = std::stod(v);
        else if (k == "malhaF" && (v == "media" || v == "grossa")) op.malhaF = v;
        else if (k == "grava_malha") op.gravaMalha = v;
        else if (k == "psi") op.psi = std::stod(v);
        else if (k == "so_transiente") op.soTransiente = std::stoi(v) != 0;
        else {
            std::fprintf(stderr, "opcao desconhecida: %s\n", argv[i]);
            std::exit(1);
        }
    }
    if (!tempos.empty()) op.solo.tempos = tempos;
    if (op.dt < 0.) op.dt = op.solo.nome == "areia" ? 0.05 : 0.25;
    return op;
}

// ---- 1. a malha: colunas a cada dx m (passando pelos vertices A..E), camadas em fracoes da altura de cada coluna,
// finas perto da superficie (a primeira com eta0 da altura) e crescendo na razao r ate etamax
static TPZGeoMesh *MalhaCho(REAL dx, REAL eta0, REAL r, REAL etamax) {
    LayerGMeshData d;
    d.superficie = {{0., 28.}, {20., 26.}, {28., 22.}, {40., 12.}, {60., 10.}};
    d.idSuperficie = {topleft, face, face, topright};
    d.idMaterial = matSolo, d.idBase = base, d.idEsquerda = left, d.idDireita = right;
    const int nx = int(std::lround(60. / dx));
    for (int i = 0; i <= nx; i++) d.xs.push_back(60. * i / nx);
    d.eta = {1.};
    for (REAL q : GeometricSteps(1., eta0, r, etamax)) d.eta.push_back(d.eta.back() - q);
    d.eta.back() = 0.;
    std::reverse(d.eta.begin(), d.eta.end());
    d.quadrilateros = false;
    return LayerGMesh(d);
}

// ---- 2. o material: srm = false, o do transiente (esqueleto elastico: coesao enorme, acoplado); srm = true, o do
// fator (Mohr-Coulomb com c' e phi' do solo, as pressoes congeladas, sem acoplamento)
static REAL gPsi = -1.; // a dilatancia do fator (graus; < 0: associado)
// T: TUnsatHP (Mohr-Coulomb PV2, associado) ou TUnsatPV (Mohr-Coulomb PV, com a dilatancia gPsi)
template <class T>
static T *MaterialCho(const SoloCho &s, bool srm) {
    const REAL phi = s.phi * M_PI / 180.;
    const REAL psi = srm && gPsi >= 0. ? gPsi * M_PI / 180. : phi;
    const REAL E = 9. * Kb * Gb / (3. * Kb + Gb), nu = (3. * Kb - 2. * Gb) / (2. * (3. * Kb + Gb));
    TPZElasticResponse ER;
    ER.SetEngineeringData(E, nu);
    auto *mat = new T(matSolo);
    if constexpr (std::is_same_v<T, TUnsatPV>) {
        TMCPV pl;
        pl.fYC.SetUp(phi, psi, srm ? s.c : 1.e9, ER);
        pl.fER = ER;
        mat->SetPlasticModel(pl);
    } else
        mat->SetPlasticModel(TMC(TPZYCMohrCoulombPV2(phi, phi, srm ? s.c : 1.e9, ER), ER));
    TPZUnsaturatedHydraulics h;
    h.fRetention = TPZUnsaturatedHydraulics::EVanGenuchten;
    h.fAVG = 1. / s.P0, h.fBVG = 1. / (1. - s.a), h.fSwr = 0.;
    h.fKrw = TPZUnsaturatedHydraulics::EMualem;
    h.fKra = TPZUnsaturatedHydraulics::EMualemAir;
    // n_VG = 1/(1 - a) < 2 nos dois solos (1,80 e 1,51): a derivada de k_r de Mualem explode na saturacao e o Newton
    // diverge; k_r da agua tabelado ate S_e = 1 - 1e-3 e fechado suavemente (MualemTabelado.h, como no Kang2020)
    if (h.fBVG < 2.) UsaTabela(h, MualemTabelado(h.fBVG, h.fSwr, h.fAVG, 1.e-3, EFechoVG::C1));
    mat->SetHydraulics(h);
    mat->SetBishop(T::EChiSaturation); // chi = S_w
    mat->SetPorosity(s.n);
    mat->SetBiot(1., 0.);
    mat->SetWaterBulkModulus(2.e6);
    mat->SetConductivity(s.ks, gammaW, 56.);
    mat->SetAir(1.25 * 9.81e-3, 101.325);
    mat->SetUnitWeightMixture(gammaD / (1. - s.n)); // gamma = (1 - n) gamma_s + n S_w gamma_w = gamma_d + n S_w gamma_w
    mat->SetLumpedStorage(true);
    mat->SetHydroMechanicalCoupling(!srm);
    return mat;
}

// os apoios de u: base engastada, laterais presas so em x
template <class T>
static void InsereApoios(UnsatModel<T> &m) {
    m.AddDirichletBC(base, T::EDirichletU, 0, {0, 1});
    m.AddDirichletBC(left, T::EDirichletU, 0, {0});
    m.AddDirichletBC(right, T::EDirichletU, 0, {0});
}

// ---- 4. o fator (copiado do montecarlo.cpp, so com o que o Cho usa)
template <class T>
struct Fator {
    std::unique_ptr<UnsatModel<T>> m;
    std::unique_ptr<UnsatAnalysis<T>> a;
};
template <class T>
static void NovoFator(Fator<T> &f, TPZGeoMesh *g, const SoloCho &s) {
    f.a.reset();
    f.m = std::make_unique<UnsatModel<T>>();
    f.m->SetGeoMesh(g);
    InsereApoios(*f.m);
    T *mat = MaterialCho<T>(s, true);
    f.m->SetMaterial(mat);
    f.m->CreateCompMesh({{2, 2}, {1, 1}, {1, 1}});
    f.a = std::make_unique<UnsatAnalysis<T>>(*f.m);
    f.a->SetVerbose(false);
    f.a->SetLoadFactorFunction([mat](REAL l) { mat->SetGravityFactor(l); }, false);
    f.a->SetStrengthReductionFunction([mat](REAL F) { mat->PlasticModel().SetStrengthReductionFactor(F); });
    f.a->SetSRMStart(0.8);
}
template <class T>
static REAL CalculaFator(Fator<T> &f, const PointLocator &loc, const std::vector<REAL> &pw, const std::vector<REAL> *pa) {
    f.m->Material()->ResetMemory();
    f.a->SetFrozenPressures([&](int espaco, const TPZVec<REAL> &x) {
        if (espaco == 2 && pa) return loc.Interpolate(x[0], x[1], *pa);
        if (espaco != 1) return 0.;
        return loc.Interpolate(x[0], x[1], pw);
    });
    return f.a->SRMArcLength();
}
// refinamento adaptativo do zero: a malha grossa e nref niveis (marcados: sqrt(J2(eps_p)) >= 0,1 do maior)
// (tolfs > 0: para antes, quando o FS muda menos que tolfs de um nivel para o outro; malha: media ou grossa)
template <class T>
static REAL FatorAdaptativo(const SoloCho &s, int nref, REAL tolfs, const std::string &malha, const PointLocator &loc,
                            const std::vector<REAL> &pw, const std::vector<REAL> *pa, std::vector<REAL> &niveis) {
    Fator<T> f;
    TPZGeoMesh *g = malha == "grossa" ? MalhaCho(4., 0.1, 1.3, 0.25) : MalhaCho(2., 0.04, 1.2, 0.12);
    niveis.clear();
    for (int k = 0;; k++) {
        NovoFator(f, g, s);
        const REAL fs = CalculaFator(f, loc, pw, pa);
        niveis.push_back(fs);
        if (k == nref) return fs;
        if (tolfs > 0. && k > 0 && std::abs(fs - niveis[k - 1]) < tolfs) return fs;
        const std::vector<REAL> ind = f.a->PlasticIndicator();
        const REAL vmax = *std::max_element(ind.begin(), ind.end());
        std::set<int64_t> marcados;
        for (TPZGeoEl *gel : g->ElementVec())
            if (gel && gel->Dimension() == 2 && ind[gel->Index()] >= 0.1 * vmax) marcados.insert(gel->Index());
        std::vector<int64_t> pai;
        g = f.m->RefinedGeoMesh(marcados, pai);
    }
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
    gPsi = op.psi;
    const SoloCho &s = op.solo;
    if (!op.gravaMalha.empty()) {
        std::unique_ptr<TPZGeoMesh> a(MalhaCho(op.dx, op.eta0, op.rT, op.etamaxT)), b(MalhaCho(4., 0.1, 1.3, 0.25)),
            c(MalhaCho(2., 0.04, 1.2, 0.12));
        GravaMalha(a.get(), op.gravaMalha + "_transiente.csv");
        GravaMalha(b.get(), op.gravaMalha + "_fator_grossa.csv");
        GravaMalha(c.get(), op.gravaMalha + "_fator_media.csv");
        return 0;
    }
    const REAL I = 3. * s.ks * 3.6e6; // 3 k_s em mm/h
    std::printf("Cho (2016): %s, modelo %s, chuva %.1f mm/h (3 k_s) ate %g h, dt = %g h\n", s.nome.c_str(),
                op.tp ? "tp (agua e ar)" : "sp (so a agua)", I, s.tfim, op.dt);

    // 1-3: a malha, o modelo e o transiente
    UnsatModel<TUnsatHP> mT;
    mT.SetGeoMesh(MalhaCho(op.dx, op.eta0, op.rT, op.etamaxT));
    InsereApoios(mT);
    mT.SetMaterial(MaterialCho<TUnsatHP>(s, false));
    mT.CreateCompMesh({{2, 2}, {1, 1}, {1, 1}});
    TPZUnsatAnalysis aT(mT.CompMesh(), {{base, {true, true}}, {right, {true, false}}, {left, {true, false}}});
    aT.SetTolerances(1.e-7, 1.e-6, 80, 25);
    RainData r;
    r.intensity = I, r.duration = 1.e9, r.dt = op.dt, r.dtMin = op.dt / 64., r.lineSearch = 1;
    r.rainVertical = false;                                  // fluxo normal a superficie, como no FLAC
    r.xToe = 0., r.yToe = 6., r.Hw = 0., r.betaW = 0.;       // o nivel horizontal a 6 m
    r.suctionMax = 30.;                                      // o teto de succao inicial (Secao 4.1)
    RainTransient chuva(aT, r, gammaW, {right}, {topleft, face, topright}, {},
                        op.tp ? ERainModel::ECoupledActiveAir : ERainModel::ECoupledPassiveAir);
    std::printf("  malha do transiente: %lld elementos, %lld equacoes\n", (long long)mT.GeoMesh()->NElements(),
                (long long)mT.CompMesh()->NEquations());
    if (!chuva.InitialState()) {
        std::printf("  o estado inicial nao convergiu: %s\n", chuva.LastFailure().c_str());
        return 1;
    }
    const PointLocator loc(mT.GeoMesh());

    // o perfil da secao C-C' (x = 28 m) em cada instante, e o FS
    const std::string pre = op.saida + "_" + s.nome + "_" + (op.tp ? "tp" : "sp");
    std::ofstream ff(pre + "_fs.csv"), fp(pre + "_cc.csv");
    ff << "t_h,FS,niveis,pw_min,pw_max,pa_max,passos,s\n";
    fp << "t_h,y_m,pw_kPa,pa_kPa\n";
    std::printf("  %6s %8s %9s %9s %9s %6s %8s\n", "t(h)", "FS", "pw min", "pw max", "pa max", "passos", "s");
    // os pontos de monitoramento da Fig. 1 (posicoes lidas da figura, +-0,2 m): H, I e J na vertical x = 34 m
    const std::vector<std::pair<std::string, std::array<REAL, 2>>> pontos = {
        {"H", {34., 16.1}}, {"I", {34., 15.1}}, {"J", {34., 12.2}}};
    std::ofstream fm(pre + "_pontos.csv");
    fm << "t_h,ponto,x,y,pw_kPa,pa_kPa,succao_kPa\n";
    auto registra = [&](REAL t) {
        const std::vector<REAL> pw = chuva.NodalWaterPressure(), pa = chuva.NodalAirPressure();
        for (const auto &[nome, x] : pontos) {
            const REAL w = loc.Interpolate(x[0], x[1], pw), a = loc.Interpolate(x[0], x[1], pa);
            fm << t << "," << nome << "," << x[0] << "," << x[1] << "," << w << "," << a << "," << a - w << "\n";
        }
        fm.flush();
    };
    registra(0.);
    REAL tMon = op.dtMon;
    for (REAL t : s.tempos) {
        for (; t > 0. && tMon <= t + 1.e-9; tMon += op.dtMon) { // avanca de hora em hora, registrando os pontos
            if (!chuva.AdvanceTo(tMon)) {
                std::printf("  o transiente parou em %.3f h: %s\n", chuva.Time(), chuva.LastFailure().c_str());
                return 1;
            }
            registra(tMon);
        }
        if (t > 0. && chuva.Time() < t - 1.e-9 && !chuva.AdvanceTo(t)) {
            std::printf("  o transiente parou em %.3f h: %s\n", chuva.Time(), chuva.LastFailure().c_str());
            return 1;
        }
        const std::vector<REAL> pw = chuva.NodalWaterPressure(), pa = chuva.NodalAirPressure();
        for (REAL y = 0.; y <= 22. + 1e-9; y += 0.25)
            fp << t << "," << y << "," << loc.Interpolate(28., std::min(y, 21.999), pw) << ","
               << loc.Interpolate(28., std::min(y, 21.999), pa) << "\n";
        std::vector<REAL> niveis;
        const REAL fs = op.soTransiente ? -1. : (op.psi >= 0. ? FatorAdaptativo<TUnsatPV>(s, op.nref, op.tolfs, op.malhaF, loc, pw, op.tp ? &pa : nullptr, niveis)
                                                         : FatorAdaptativo<TUnsatHP>(s, op.nref, op.tolfs, op.malhaF, loc, pw, op.tp ? &pa : nullptr, niveis));
        std::string nv;
        for (REAL x : niveis) nv += (nv.empty() ? "" : ";") + std::to_string(x);
        const auto pr = chuva.PressureRange(), ar = chuva.AirPressureRange();
        const double sec = std::chrono::duration<double>(Clock::now() - t0).count();
        ff << t << "," << fs << "," << nv << "," << pr.first << "," << pr.second << "," << ar.second << ","
           << chuva.Steps() << "," << sec << "\n";
        ff.flush();
        fp.flush();
        std::printf("  %6g %8.4f %9.2f %9.2f %9.2f %6lld %8.0f\n", t, fs, pr.first, pr.second, ar.second,
                    (long long)chuva.Steps(), sec);
        std::fflush(stdout);
    }
    std::printf("  CSV: %s_fs.csv, %s_cc.csv\n", pre.c_str(), pre.c_str());
    return 0;
}
