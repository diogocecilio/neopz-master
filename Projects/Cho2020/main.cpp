// Cho2020: a reproducao de S. E. Cho (2020), "Failure distribution analysis of shallow landslides under rainfall
// infiltration based on fragility curves", Landslides 17:79-91, escrita passo a passo.
// O exemplo: um talude infinito a beta = 30 graus com a rocha impermeavel a D = 2 m (na vertical; L = D cos(beta) na
// normal), o nivel d'agua inicial na rocha e o estado inicial hidrostatico, p_w = -gamma_w (L - z*) cos(beta) (z* a
// profundidade normal), uma areia de granito alterado (van Genuchten - Mualem, Tabela 1) e chuva constante p de
// 0,2 a 1,2 k_s. O artigo calcula a infiltracao com um Green-Ampt modificado (Cho 2017) e o FS do talude infinito
// (Eq. 21) com sigma_s = S_e p_w (Eq. 22, u_a = 0); aqui a infiltracao e o Richards do NeoPZ na coluna normal ao talude
// do Kang2020 (Passo 7 de la) e o FS e o mesmo (talude_infinito.h do Kang2020; com theta_r = 0, chi = S_w = S_e).
//
// Uso: Cho2020 [passo=1] [p=<p/k_s>] [ny=80] [dt=0.02] [tfim=16] [modelo=iii|ii|i] [eps=1e-3] [saida=cho20]
//
// Passo 1 (este): a analise deterministica da Secao "Deterministic analysis" (Fig. 3): com os valores medios das
//   Tabelas 1 e 2, a profundidade da frente de molhamento (Fig. 3a) e o FS (Fig. 3b) contra o tempo, para uma chuva p.
//   A coluna normal (como no Kang2020): a malha tem a altura L = D cos(beta); a gravidade ao longo dela e g cos(beta)
//   (SetGravityDirection), a chuva vertical entra como p cos(beta) por metro de superficie inclinada, e o estado
//   inicial hidrostatico com o nivel na base sai do RainTransient com gamma_w cos(beta) (o gamma_w dele so entra no
//   estado inicial: p_w = gamma_w cos(beta) (0 - y), y a altura normal acima da rocha). No topo, a chuva ou a face de
//   percolacao (p_w = 0, o excesso escorre); a base e impermeavel. O ar e passivo (p_a = 0): modelo iii (o padrao,
//   esqueleto rigido, so a percolacao, como o modelo do artigo), ii (acoplado) ou i (ar ativo).
//   A frente: o ponto mais fundo em que S_w ja subiu metade do que faltava para saturar, S_w >= (S_w0 + 1) / 2; sai na
//   profundidade vertical z = z* / cos(beta), como a Fig. 3a. O FS: o minimo da Eq. 21 ao longo do perfil.

#include "talude_infinito.h" // do Kang2020: FSInfinito (Eqs. 14 e 15 do Kang = Eq. 21 do Cho com u_a = 0)

#include "HPMaterial.h"
#include "MualemTabelado.h"
#include "RainTransient.h"
#include "UnsatModel.h"

#include "Plasticity/TPZElasticResponse.h"
#include "Plasticity/TPZMatUnsaturatedUPwPa.h"
#include "Plasticity/TPZPlasticStepVoigt.h"
#include "Plasticity/TPZYCMohrCoulombPV2.h"

#include "TPZGeoLinear.h"
#include "pzgeotriangle.h"
#include "tpzgeoelrefpattern.h"

#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

using TMC = TPZPlasticStepVoigt<TPZYCMohrCoulombPV2, TPZElasticResponse>;
using TUnsat = TPZMatUnsaturatedUPwPa<TMC, TPZUnsaturatedMem>;
using TUnsatHP = HPMaterial<TUnsat>;
using kang::FSInfinito;
using kang::PontoPerfil;

enum EIds { matSolo = 1, base = -1, direita = -2, esquerda = -5, topo = -6 };

// ---- os dados do artigo: os valores medios das Tabelas 1 e 2
namespace cho20 {
constexpr REAL thetaS = 0.355, thetaR = 0.;   // a porosidade e theta_s
constexpr REAL alfa = 0.410, nVG = 1.120;     // 1/kPa e o n de van Genuchten (m = 1 - 1/n)
constexpr REAL ks = 1.e-5;                    // m/s (3,6 cm/h)
constexpr REAL gammaD = 17.01, c = 5., phiDeg = 35., betaDeg = 30.;
constexpr REAL D = 2.;                        // m: a rocha, na vertical
constexpr REAL gammaW = 9.81;
// o esqueleto (o artigo nao precisa: hipotese nossa; so pesa nos modelos acoplados ii e i)
constexpr REAL E = 20.e3, nu = 0.3;
} // namespace cho20
using namespace cho20;

struct Opcoes {
    REAL p = 1.;     // a chuva em multiplos de k_s
    int ny = 80;     // celulas ao longo da coluna
    REAL dt = 0.02;  // h
    REAL tfim = 16.; // h
    std::string modelo = "iii";
    std::string saida = "cho20";
    REAL eps = 1.e-3; // o k_r tabelado ate S_e = 1 - eps, fechado suavemente ate a saturacao (MualemTabelado.h)
};

static Opcoes LeOpcoes(int argc, char **argv) {
    Opcoes op;
    for (int i = 1; i < argc; i++) {
        const std::string s = argv[i];
        const size_t eq = s.find('=');
        const std::string k = s.substr(0, eq), v = eq == std::string::npos ? "" : s.substr(eq + 1);
        if (k == "passo" && v == "1") continue;
        else if (k == "p") op.p = std::stod(v);
        else if (k == "ny") op.ny = std::stoi(v);
        else if (k == "dt") op.dt = std::stod(v);
        else if (k == "tfim") op.tfim = std::stod(v);
        else if (k == "modelo" && (v == "iii" || v == "ii" || v == "i")) op.modelo = v;
        else if (k == "saida") op.saida = v;
        else if (k == "eps") op.eps = std::stod(v);
        else {
            std::fprintf(stderr, "opcao desconhecida: %s\n", argv[i]);
            std::exit(1);
        }
    }
    return op;
}

// A coluna [0, w] x [0, altura] com ny celulas, cada uma em 4 triangulos em volta do centro (a MalhaColuna do Kang2020):
// cantos 2j (x = 0) e 2j + 1 (x = w) na cota y_j = altura j / ny, centros depois
static TPZGeoMesh *MalhaColuna(REAL altura, REAL w, int ny) {
    auto *g = new TPZGeoMesh();
    g->SetDimension(2);
    const int64_t nCantos = 2 * int64_t(ny + 1);
    g->NodeVec().Resize(nCantos + ny);
    auto poe = [g](int64_t id, REAL x, REAL y) {
        TPZManVector<REAL, 3> c = {x, y, 0.};
        g->NodeVec()[id] = TPZGeoNode(id, c, *g);
    };
    auto Y = [altura, ny](int j) { return j == ny ? altura : altura * j / ny; };
    for (int j = 0; j <= ny; j++) poe(2 * j, 0., Y(j)), poe(2 * j + 1, w, Y(j));
    for (int j = 0; j < ny; j++) poe(nCantos + j, 0.5 * w, 0.5 * (Y(j) + Y(j + 1)));
    for (int j = 0; j < ny; j++) {
        const int64_t bl = 2 * j, br = 2 * j + 1, tl = 2 * j + 2, tr = 2 * j + 3, cc = nCantos + j;
        for (const auto &t : {std::array<int64_t, 3>{bl, br, cc}, std::array<int64_t, 3>{br, tr, cc},
                              std::array<int64_t, 3>{tr, tl, cc}, std::array<int64_t, 3>{tl, bl, cc}}) {
            TPZManVector<int64_t, 3> v = {t[0], t[1], t[2]};
            new TPZGeoElRefPattern<pzgeom::TPZGeoTriangle>(v, matSolo, *g);
        }
    }
    auto linha = [g](int64_t a, int64_t b, int id) {
        TPZManVector<int64_t, 2> v = {a, b};
        new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(v, id, *g);
    };
    linha(0, 1, base);
    for (int j = 0; j < ny; j++) linha(2 * j + 1, 2 * j + 3, direita);
    linha(2 * int64_t(ny) + 1, 2 * int64_t(ny), topo);
    for (int j = ny; j > 0; j--) linha(2 * j, 2 * j - 2, esquerda);
    g->BuildConnectivity();
    return g;
}

// O material: van Genuchten - Mualem (n = 1,12 < 2: k_r tabelado ate S_e = 1 - 1e-3, como no Kang2020 e no Cho2016),
// chi = S_w (= S_e, theta_r = 0), peso gamma_d + n S_w gamma_w, esqueleto elastico (coesao enorme: o FS e a Eq. 21)
static TUnsatHP *MaterialCho(REAL beta, REAL eps) {
    TPZElasticResponse ER;
    ER.SetEngineeringData(E, nu);
    auto *mat = new TUnsatHP(matSolo);
    const REAL phi = phiDeg * M_PI / 180.;
    mat->SetPlasticModel(TMC(TPZYCMohrCoulombPV2(phi, phi, 1.e9, ER), ER));
    TPZUnsaturatedHydraulics h;
    h.fRetention = TPZUnsaturatedHydraulics::EVanGenuchten;
    h.fAVG = alfa, h.fBVG = nVG, h.fSwr = thetaR / thetaS;
    h.fKrw = TPZUnsaturatedHydraulics::EMualem;
    h.fKra = TPZUnsaturatedHydraulics::EMualemAir;
    UsaTabela(h, MualemTabelado(h.fBVG, h.fSwr, h.fAVG, eps, EFechoVG::C1));
    mat->SetHydraulics(h);
    mat->SetBishop(TUnsat::EChiSaturation);
    mat->SetPorosity(thetaS);
    mat->SetBiot(1., 0.);
    mat->SetWaterBulkModulus(2.2e6);
    mat->SetConductivity(ks, gammaW, 56.);
    mat->SetAir(1.2 * 9.81e-3, 101.325);
    mat->SetUnitWeightMixture(gammaD / (1. - thetaS)); // gamma_d + n S_w gamma_w
    mat->SetLumpedStorage(true);
    mat->SetHydroMechanicalCoupling(true);
    mat->SetGravityDirection(TPZManVector<REAL, 3>{0., -std::cos(beta), 0.}); // a coluna normal: g cos(beta)
    return mat;
}

int main(int argc, char **argv) {
    using Clock = std::chrono::steady_clock;
    const auto t0 = Clock::now();
    const Opcoes op = LeOpcoes(argc, argv);
    const REAL beta = betaDeg * M_PI / 180., cb = std::cos(beta), L = D * cb;
    const REAL pmmh = op.p * ks * 3.6e6; // a chuva (mm/h, por metro horizontal)
    std::printf("Cho (2020), Passo 1: coluna normal ao talude, beta = %g graus, L = %.4f m (D = %g m), p = %g k_s = %.2f "
                "mm/h, modelo %s, ny = %d, dt = %g h\n", betaDeg, L, D, op.p, pmmh, op.modelo.c_str(), op.ny, op.dt);

    UnsatModel<TUnsatHP> m;
    m.SetGeoMesh(MalhaColuna(L, L / op.ny, op.ny));
    m.AddDirichletBC(base, TUnsat::EDirichletU, 0, {0, 1});
    m.AddDirichletBC(esquerda, TUnsat::EDirichletU, 0, {0});
    m.AddDirichletBC(direita, TUnsat::EDirichletU, 0, {0});
    TUnsatHP *mat = MaterialCho(beta, op.eps);
    m.SetMaterial(mat);
    m.CreateCompMesh({{2, 2}, {1, 1}, {1, 1}});
    TPZUnsatAnalysis an(m.CompMesh(), {{base, {true, true}}, {direita, {true, false}}, {esquerda, {true, false}}});
    an.SetTolerances(1.e-7, 1.e-6, 80, 25);
    RainData r;
    r.intensity = pmmh * cb; // a chuva vertical por metro da superficie inclinada
    r.duration = 1.e9, r.dt = op.dt, r.dtMin = op.dt / 256., r.lineSearch = 1;
    r.rainVertical = false;
    r.xToe = 0., r.yToe = 0., r.Hw = 0., r.betaW = 0.; // o nivel na base (y = 0)
    r.suctionMax = 1.e6;
    const ERainModel modelo = op.modelo == "iii" ? ERainModel::EUncoupled
                              : op.modelo == "ii" ? ERainModel::ECoupledPassiveAir
                                                  : ERainModel::ECoupledActiveAir;
    RainTransient chuva(an, r, gammaW * cb, {}, {topo}, {}, modelo); // gamma_w cos(beta): o hidrostatico normal
    if (!chuva.InitialState()) {
        std::printf("  o estado inicial nao convergiu: %s\n", chuva.LastFailure().c_str());
        return 1;
    }
    std::printf("  %lld equacoes\n", (long long)m.CompMesh()->NEquations());

    // o perfil no eixo x = 0, de cima para baixo (nos 2j), na profundidade normal z* = L - y
    std::vector<int64_t> nos;
    std::vector<REAL> zn;
    for (int j = op.ny; j >= 0; j--) nos.push_back(2 * j), zn.push_back(L - L * j / op.ny);
    const TPZUnsaturatedHydraulics &h = mat->Hydraulics();
    auto Sw = [&h](REAL succao) {
        REAL S, G;
        h.Retention(succao, S, G);
        return S;
    };
    std::vector<REAL> S0(nos.size());
    {
        const std::vector<REAL> pw = chuva.NodalWaterPressure();
        for (size_t i = 0; i < nos.size(); i++) S0[i] = Sw(-pw[nos[i]]);
    }

    const std::string pre = op.saida + "_p" + std::to_string(op.p).substr(0, 4) + "_" + op.modelo;
    std::ofstream ff(pre + "_fs.csv"), fp(pre + "_perfis.csv");
    ff << "t_h,frente_m,FS,zmin_m,pw_base_kPa,face,passos\n";
    fp << "t_h,z_m,pw_kPa,Sw\n";
    REAL tChega = -1.;
    const int nSaidas = int(std::lround(op.tfim / 0.1));
    std::printf("  %6s %8s %8s %8s %9s %4s %7s\n", "t(h)", "frente", "FS", "z_min", "pw base", "face", "passos");
    for (int k = 0; k <= nSaidas; k++) {
        const REAL t = 0.1 * k;
        if (k > 0 && !chuva.AdvanceTo(t)) {
            std::printf("  o transiente parou em %.3f h: %s\n", chuva.Time(), chuva.LastFailure().c_str());
            return 1;
        }
        const std::vector<REAL> pw = chuva.NodalWaterPressure();
        std::vector<PontoPerfil> perf;
        REAL frente = 0.;
        for (size_t i = 0; i < nos.size(); i++) {
            const REAL p = pw[nos[i]], S = Sw(-p), z = zn[i] / cb; // a profundidade vertical
            perf.push_back({z, S, p, 0.});
            if (S >= 0.5 * (S0[i] + 1.) && (i == 0 || frente >= zn[i - 1] / cb - 1.e-12)) frente = z;
        }
        REAL zmin = 0.;
        const REAL fs = FSInfinito(perf, beta, phiDeg * M_PI / 180., c, gammaD, thetaS, gammaW, nullptr, &zmin);
        if (tChega < 0. && frente >= D - 1.e-9) tChega = t;
        ff << t << "," << frente << "," << fs << "," << zmin << "," << pw[nos.back()] << "," << chuva.FaceNodes() << ","
           << chuva.Steps() << "\n";
        if (k % 5 == 0)
            for (const PontoPerfil &p : perf) fp << t << "," << p.z << "," << p.pw << "," << p.Sw << "\n";
        ff.flush();
        if (k % 10 == 0)
            std::printf("  %6.1f %8.3f %8.4f %8.3f %9.2f %4zu %7lld\n", t, frente, fs, zmin, pw[nos.back()],
                        chuva.FaceNodes(), (long long)chuva.Steps());
        std::fflush(stdout);
    }
    std::printf("  a frente chega a rocha em %s; CSV: %s_fs.csv, %s_perfis.csv; tempo %.1f s\n",
                tChega < 0. ? "(nao chegou)" : (std::to_string(tChega).substr(0, 5) + " h").c_str(), pre.c_str(),
                pre.c_str(), std::chrono::duration<double>(Clock::now() - t0).count());
    return 0;
}
