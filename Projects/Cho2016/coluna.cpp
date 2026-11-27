// Cho2016Coluna: a coluna vertical 1D do silte (ou da areia) do Cho (2016) sob a chuva de 3 k_s, para separar o erro de
// malha da fisica. A comparacao com o artigo e com o Richards 1D independente (r1d.py, ~/projects/runs/cho2016/richards1d)
// mostrou que o FLAC adianta a frente na malha de 0,27 m; aqui o NeoPZ e refinado ate convergir:
//   sp: so a agua (modelo iii, esqueleto rigido; ou ii, acoplado), tp: agua e ar (modelo i, ar ativo, p_a = 0 no topo)
// A coluna tem a altura H (10 m), camadas que crescem do topo (h0, razao r, ate hmax), cada camada em 4 triangulos em
// volta do centro (a MalhaColuna do Kang2020), largura w. Estado inicial: p_w = -30 kPa em toda a coluna (o teto de
// succao do artigo, com o nivel d'agua bem abaixo), p_a = 0. Topo: a chuva 3 k_s ou a face de percolacao; base e laterais
// impermeaveis (e estanques ao ar).
//
// Uso: Cho2016Coluna [solo=silte|areia] [modelo=iii|ii|i] [h0=0.01] [r=1.1] [hmax=0.05] [H=10] [w=0.05] [dt=0.05]
//                    [tfim=160] [kramin=1e-4] [mu=56] [eps=1e-3] [saida=col]
// Saida: <saida>_frente.csv (t_h, frente_m (p_w = -15 kPa, interpolada), pa_max_kPa, pw_topo_max_kPa, infiltrado_m,
//        passos) e <saida>_perfis.csv (t_h, z_m, pw_kPa, pa_kPa) em 16, 32, 64, 96, 128 e 160 h

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

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <set>
#include <string>
#include <vector>

using TMC = TPZPlasticStepVoigt<TPZYCMohrCoulombPV2, TPZElasticResponse>;
using TUnsat = TPZMatUnsaturatedUPwPa<TMC, TPZUnsaturatedMem>;
using TUnsatHP = HPMaterial<TUnsat>;

enum EIds { matSolo = 1, base = -1, direita = -2, esquerda = -5, topo = -6 };
constexpr REAL gammaW = 9.81, gammaD = 1640. * 9.81e-3, Kb = 33.3e3, Gb = 11.1e3;

struct Opcoes {
    std::string solo = "silte", modelo = "iii", saida = "col";
    REAL h0 = 0.01, r = 1.1, hmax = 0.05, H = 10., w = 0.05, dt = 0.05, tfim = 160., kramin = 1.e-4, mu = 56.,
         eps = 1.e-3;
};

static Opcoes LeOpcoes(int argc, char **argv) {
    Opcoes op;
    for (int i = 1; i < argc; i++) {
        const std::string s = argv[i];
        const auto e = s.find('=');
        const std::string k = s.substr(0, e), v = e == std::string::npos ? "" : s.substr(e + 1);
        if (k == "solo") op.solo = v;
        else if (k == "modelo") op.modelo = v;
        else if (k == "saida") op.saida = v;
        else if (k == "h0") op.h0 = std::stod(v);
        else if (k == "r") op.r = std::stod(v);
        else if (k == "hmax") op.hmax = std::stod(v);
        else if (k == "H") op.H = std::stod(v);
        else if (k == "w") op.w = std::stod(v);
        else if (k == "dt") op.dt = std::stod(v);
        else if (k == "tfim") op.tfim = std::stod(v);
        else if (k == "kramin") op.kramin = std::stod(v);
        else if (k == "mu") op.mu = std::stod(v);
        else if (k == "eps") op.eps = std::stod(v);
        else {
            std::fprintf(stderr, "opcao desconhecida: %s\n", argv[i]);
            std::exit(1);
        }
    }
    return op;
}

// as cotas das camadas, de baixo (0) para cima (H): espessuras h0, h0 r, ... ate hmax, a partir do topo
static std::vector<REAL> Cotas(const Opcoes &op) {
    std::vector<REAL> prof = {0.};
    REAL h = op.h0;
    while (prof.back() < op.H - 1.e-9) {
        prof.push_back(std::min(prof.back() + h, op.H));
        h = std::min(h * op.r, op.hmax);
    }
    if (prof.size() > 2 && prof.back() - prof[prof.size() - 2] < 0.3 * prof[prof.size() - 2] - 0.3 * prof[prof.size() - 3])
        prof.erase(prof.end() - 2); // a ultima camada nao fica fina demais
    std::vector<REAL> y;
    for (auto it = prof.rbegin(); it != prof.rend(); ++it) y.push_back(op.H - *it);
    return y;
}

// a MalhaColuna do Kang2020 com cotas dadas: cantos 2j (x = 0) e 2j + 1 (x = w), centros depois
static TPZGeoMesh *MalhaColuna(const std::vector<REAL> &Y, REAL w) {
    const int ny = int(Y.size()) - 1;
    auto *g = new TPZGeoMesh();
    g->SetDimension(2);
    const int64_t nCantos = 2 * int64_t(ny + 1);
    g->NodeVec().Resize(nCantos + ny);
    auto poe = [g](int64_t id, REAL x, REAL y) {
        TPZManVector<REAL, 3> c = {x, y, 0.};
        g->NodeVec()[id] = TPZGeoNode(id, c, *g);
    };
    for (int j = 0; j <= ny; j++) poe(2 * j, 0., Y[j]), poe(2 * j + 1, w, Y[j]);
    for (int j = 0; j < ny; j++) poe(nCantos + j, 0.5 * w, 0.5 * (Y[j] + Y[j + 1]));
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

int main(int argc, char **argv) {
    using Clock = std::chrono::steady_clock;
    const auto t0 = Clock::now();
    const Opcoes op = LeOpcoes(argc, argv);
    const bool silte = op.solo != "areia";
    const REAL ks = silte ? 3.e-6 : 1.e-5, P0 = silte ? 4.905 : 1.5, a = silte ? 0.336 : 0.445,
               n = silte ? 0.5 : 0.35, phi = (silte ? 30. : 35.) * M_PI / 180.;
    const std::vector<REAL> Y = Cotas(op);
    const int ny = int(Y.size()) - 1;
    std::printf("Cho (2016), coluna 1D: %s, modelo %s, %d camadas (h0 %g, r %g, hmax %g), w %g, dt %g h, k_ra min %g, "
                "mu_w/mu_a %g\n", op.solo.c_str(), op.modelo.c_str(), ny, op.h0, op.r, op.hmax, op.w, op.dt, op.kramin,
                op.mu);

    UnsatModel<TUnsatHP> m;
    m.SetGeoMesh(MalhaColuna(Y, op.w));
    m.AddDirichletBC(base, TUnsat::EDirichletU, 0, {0, 1});
    m.AddDirichletBC(esquerda, TUnsat::EDirichletU, 0, {0});
    m.AddDirichletBC(direita, TUnsat::EDirichletU, 0, {0});
    TPZElasticResponse ER;
    ER.SetEngineeringData(9. * Kb * Gb / (3. * Kb + Gb), (3. * Kb - 2. * Gb) / (2. * (3. * Kb + Gb)));
    auto *mat = new TUnsatHP(matSolo);
    mat->SetPlasticModel(TMC(TPZYCMohrCoulombPV2(phi, phi, 1.e9, ER), ER));
    TPZUnsaturatedHydraulics h;
    h.fRetention = TPZUnsaturatedHydraulics::EVanGenuchten;
    h.fAVG = 1. / P0, h.fBVG = 1. / (1. - a), h.fSwr = 0.;
    h.fKrw = TPZUnsaturatedHydraulics::EMualem;
    h.fKra = TPZUnsaturatedHydraulics::EMualemAir;
    UsaTabela(h, MualemTabelado(h.fBVG, h.fSwr, h.fAVG, op.eps, EFechoVG::C1));
    h.fKraMin = op.kramin;
    mat->SetHydraulics(h);
    mat->SetBishop(TUnsat::EChiSaturation);
    mat->SetPorosity(n);
    mat->SetBiot(1., 0.);
    mat->SetWaterBulkModulus(2.e6);
    mat->SetConductivity(ks, gammaW, op.mu);
    mat->SetAir(1.25 * 9.81e-3, 101.325);
    mat->SetUnitWeightMixture(gammaD / (1. - n));
    mat->SetLumpedStorage(true);
    mat->SetHydroMechanicalCoupling(op.modelo != "iii");
    m.SetMaterial(mat);
    m.CreateCompMesh({{2, 2}, {1, 1}, {1, 1}});
    TPZUnsatAnalysis an(m.CompMesh(), {{base, {true, true}}, {direita, {true, false}}, {esquerda, {true, false}}});
    an.SetTolerances(1.e-7, 1.e-6, 80, 25);
    RainData r;
    r.intensity = 3. * ks * 3.6e6;
    r.duration = 1.e9, r.dt = op.dt, r.dtMin = op.dt / 256., r.lineSearch = 1;
    r.rainVertical = false;
    r.xToe = 0., r.yToe = -100., r.Hw = 0., r.betaW = 0.; // o nivel bem abaixo: -30 kPa em toda a coluna
    r.suctionMax = 30.;
    const ERainModel modelo = op.modelo == "iii" ? ERainModel::EUncoupled
                              : op.modelo == "ii" ? ERainModel::ECoupledPassiveAir
                                                  : ERainModel::ECoupledActiveAir;
    RainTransient chuva(an, r, gammaW, {}, {topo}, {}, modelo);
    if (!chuva.InitialState()) {
        std::printf("  o estado inicial nao convergiu: %s\n", chuva.LastFailure().c_str());
        return 1;
    }
    std::printf("  %lld equacoes\n", (long long)m.CompMesh()->NEquations());

    std::ofstream ff(op.saida + "_frente.csv"), fp(op.saida + "_perfis.csv");
    ff << "t_h,frente_m,pa_max_kPa,pw_max_kPa,infiltrado_m,passos\n";
    fp << "t_h,z_m,pw_kPa,pa_kPa\n";
    const std::set<int> perfis = {16, 32, 64, 96, 128, 160};
    for (int k = 0; k <= int(std::lround(op.tfim)); k++) {
        const REAL t = k;
        if (k > 0 && !chuva.AdvanceTo(t)) {
            std::printf("  o transiente parou em %.3f h: %s\n", chuva.Time(), chuva.LastFailure().c_str());
            return 1;
        }
        const std::vector<REAL> pw = chuva.NodalWaterPressure(), pa = chuva.NodalAirPressure();
        REAL frente = 0., pamax = 0., pwmax = -1.e9;
        for (int j = ny; j > 0; j--) { // do topo para baixo, no eixo x = 0 (nos 2j)
            const REAL p0 = pw[2 * j], p1 = pw[2 * (j - 1)], z0 = op.H - Y[j], z1 = op.H - Y[j - 1];
            pwmax = std::max(pwmax, p0);
            if (p0 >= -15. && p1 < -15.) frente = z0 + (p0 + 15.) / (p0 - p1) * (z1 - z0);
        }
        for (int j = 0; j <= ny; j++) pamax = std::max(pamax, pa[2 * j]);
        ff << t << "," << frente << "," << pamax << "," << pwmax << "," << chuva.Infiltrated() / op.w << ","
           << chuva.Steps() << "\n";
        ff.flush();
        if (perfis.count(k)) {
            for (int j = ny; j >= 0; j--) fp << t << "," << op.H - Y[j] << "," << pw[2 * j] << "," << pa[2 * j] << "\n";
            fp.flush();
        }
        if (k % 8 == 0) {
            std::printf("  %5.0f h: frente %.3f m, p_a max %.2f kPa, infiltrado %.3f m, %lld passos\n", t, frente, pamax,
                        chuva.Infiltrated() / op.w, (long long)chuva.Steps());
            std::fflush(stdout);
        }
    }
    std::printf("  tempo %.1f s\n", std::chrono::duration<double>(Clock::now() - t0).count());
    return 0;
}
