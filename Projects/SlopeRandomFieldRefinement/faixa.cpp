// SRM numa faixa longa de talude infinito (passo 4 do plano de reproducao): o fator de seguranca por reducao de
// resistencia (o SRMArcLength do UnsatAnalysis, o mesmo do SlopeRainMonteCarlo, com as pressoes congeladas) de uma camada
// inclinada sobre base rigida, contra o FS fechado de talude infinito de Lu e Godt (o do SlopeRainColumn):
//   FS(z) = [c' + (W cos^2 beta - sigma_s) tan phi'] / (W sin beta cos beta),   FS_LEM = min_z FS(z),
//   W = integral de gamma(z) dz desde a superficie (z vertical), sigma_s = S_e p_w se p_w < 0, senao p_w.
// Com Mohr-Coulomb associado o SRM de um talude infinito e o minimo em z do equilibrio limite (o mecanismo e o bloco que
// desliza num plano paralelo a superficie): a faixa so difere pelas pontas, um erro em 1/L.
// A faixa, nos eixos do talude: x ao longo da superficie (para baixo), y normal a ela; [0, L] x [0, T], T = H cos beta a
// espessura normal de uma camada de profundidade vertical H. A gravidade e (sin beta, -cos beta) (SetGravityDirection).
// Um ponto na profundidade normal d = T - y esta na profundidade vertical z = d / cos beta: e o mapa da coluna vertical do
// SlopeRainColumn (o perfil p_w(z) dela, em z vertical) para a faixa. Base engastada; pontas com roletes (u_x = 0).
// O material e o do SlopeRainMonteCarlo/SlopeRainColumn: TPZMatUnsaturatedUPwPa (Mohr-Coulomb PV2 associado, chi = S_e de
// van Genuchten, p_a = 0); o peso e constante (gamma) ou o da mistura gamma_d + theta gw, theta = theta_s S_w(p_w), que
// e o W do Lu-Godt do SlopeRainColumn (SetUnitWeightMixture com gamma_s = gamma_d / (1 - theta_s), porosidade theta_s).
// O FS fechado e calculado aqui com a mesma p_w(z) (a que a malha recebe), os mesmos S_w e S_e (TPZUnsaturatedHydraulics) e
// W por trapezios numa grade de 4000 intervalos; com pressao=perfil, a tela tambem da o FS_min do proprio arquivo.
// Uso: SlopeStripSRM [chave=valor ...], com (padrao entre parenteses: os numeros do Cho 2014):
//   H (2)               profundidade vertical da camada (m); a espessura normal e H cos beta
//   LH (20)             L / H (L ao longo da superficie)
//   ny (8)              celulas na espessura; nx = round(L / (T / ny)) (celulas quadradas), 2 triangulos por celula
//                       com a diagonal alternada
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
//   verbose (0), vtk () prefixo do VTK do estado final (u, p_w, sqrt(J2(eps_p)))
//   csv (faixa.csv)     acrescenta uma linha: caso,LH,ny,nref,elementos,equacoes,FS_SRM,FS_LEM,z_LEM,z_banda,
//                       z_banda_min,z_banda_max,erro_rel,segundos
//   caso ()             rotulo da linha
// z_banda: a profundidade vertical do centro do elemento de maior sqrt(J2(eps_p)) no terco central (x em [L/3, 2L/3]);
// z_banda_min/max: a faixa dos centros com >= 0,5 do maior la.
#include "HPMaterial.h"
#include "UnsatAnalysis.h"
#include "UnsatModel.h"

#include "Plasticity/TPZElasticResponse.h"
#include "Plasticity/TPZMatUnsaturatedUPwPa.h"
#include "Plasticity/TPZPlasticStepVoigt.h"
#include "Plasticity/TPZYCMohrCoulombPV2.h"
#include "TPZGeoLinear.h"
#include "pzgeotriangle.h"
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
    int ny = 8, nref = 0;
    bool mistura = false, acimaHidro = false, verbose = false;
    std::string pressao = "seco", arquivo, vtk, csv = "faixa.csv", caso;
    REAL dw = 1., t = 0., srmF0 = 0.8;
};

Opcoes LeOpcoes(int argc, char *argv[])
{
    Opcoes o;
    const std::set<std::string> chaves = {"H",     "LH",      "ny",    "nref",    "beta",   "c",   "phi",    "E",
                                          "nu",    "peso",    "gamma", "gamma_d", "gw",     "vg",  "pressao", "dw",
                                          "acima", "arquivo", "t",     "srm_f0",  "verbose", "vtk", "csv",    "caso"};
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
    }
    if (o.ny < 1 || o.LH <= 0. || o.H <= 0.) Erro("ny, LH e H positivos");
    if (o.pressao == "perfil" && o.arquivo.empty()) Erro("pressao=perfil precisa de arquivo=");
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

struct Fator {
    std::unique_ptr<UnsatModel<TUnsatHP>> m;
    std::unique_ptr<UnsatAnalysis<TUnsatHP>> a;
};

int main(int argc, char *argv[])
{
    const auto t0 = std::chrono::steady_clock::now();
    const Opcoes o = LeOpcoes(argc, argv);
    const REAL b = o.beta * M_PI / 180., cb = std::cos(b), T = o.H * cb, L = o.LH * o.H;
    const int nx = std::max(1, int(std::lround(L / (T / o.ny))));

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

    // FS fechado (Lu-Godt) com a mesma p_w(z)
    const TPZUnsaturatedHydraulics h = Hidraulica(o);
    const REAL tphi = std::tan(o.phi * M_PI / 180.), sb = std::sin(b);
    auto gam = [&](REAL z) {
        if (!o.mistura) return o.gamma;
        return o.gammaD + o.ts * h.Evaluate(pw(z), 0.).S * o.gw;
    };
    const int N = 4000;
    REAL W = 0., fsLEM = std::numeric_limits<REAL>::infinity(), zLEM = 0.;
    for (int i = 1; i <= N; i++) {
        const REAL z0 = o.H * (i - 1) / N, z1 = o.H * i / N;
        W += 0.5 * (gam(z0) + gam(z1)) * (z1 - z0);
        const REAL p = pw(z1);
        REAL Se = 1., dSe;
        if (p < 0.) h.Effective(h.Evaluate(p, 0.).S, Se, dSe);
        const REAL ss = p < 0. ? Se * p : p;
        const REAL fs = (o.c + (W * cb * cb - ss) * tphi) / (W * sb * cb);
        if (fs < fsLEM) fsLEM = fs, zLEM = z1;
    }

    // SRM na faixa, com refinamento adaptativo opcional
    Fator f;
    TPZGeoMesh *g = MalhaFaixa(L, T, nx, o.ny);
    REAL fs = 0.;
    std::vector<REAL> niveis;
    for (int k = 0;; k++) {
        f.a.reset();
        f.m = std::make_unique<UnsatModel<TUnsatHP>>();
        f.m->SetGeoMesh(g);
        f.m->AddDirichletBC(base, TUnsat::EDirichletU, 0, {0, 1});
        f.m->AddDirichletBC(esquerda, TUnsat::EDirichletU, 0, {0});
        f.m->AddDirichletBC(direita, TUnsat::EDirichletU, 0, {0});
        TUnsatHP *mat = Material(o);
        f.m->SetMaterial(mat);
        f.m->CreateCompMesh({{2, 2}, {1, 1}, {1, 1}});
        f.a = std::make_unique<UnsatAnalysis<TUnsatHP>>(*f.m);
        f.a->SetVerbose(o.verbose);
        f.a->SetLoadFactorFunction([mat](REAL l) { mat->SetGravityFactor(l); }, false);
        f.a->SetStrengthReductionFunction([mat](REAL F) { mat->PlasticModel().SetStrengthReductionFactor(F); });
        f.a->SetSRMStart(o.srmF0);
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

    // banda no terco central
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
        if (x[0] < L / 3. || x[0] > 2. * L / 3.) continue;
        const REAL z = (T - x[1]) / cb, v = ind[gel->Index()];
        centrais.push_back({v, z});
        if (v > vmax) vmax = v, zb = z;
    }
    REAL zmin = o.H, zmax = 0.;
    for (const auto &[v, z] : centrais)
        if (v >= 0.5 * vmax) zmin = std::min(zmin, z), zmax = std::max(zmax, z);
    const int64_t neq = int64_t(f.a->FreeEquations().size());
    const REAL seg = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

    std::printf("Faixa: H = %g m (T = %.4f m normal), beta = %g, L/H = %g (L = %g m), %d x %d celulas, nref %d: %lld "
                "triangulos, %lld equacoes livres | c' %g, phi' %g, peso %s | pressao %s%s\n",
                o.H, T, o.beta, o.LH, L, nx, o.ny, o.nref, (long long)nel, (long long)neq, o.c, o.phi,
                o.mistura ? "mistura" : "constante", o.pressao.c_str(),
                o.pressao == "perfil" ? (" (" + o.arquivo + ", t = " + std::to_string(o.t) + " h)").c_str() : "");
    std::printf("  FS_SRM = %.5f", fs);
    if (niveis.size() > 1) {
        std::printf(" (niveis:");
        for (REAL v : niveis) std::printf(" %.5f", v);
        std::printf(")");
    }
    std::printf(" | FS_LEM = min_z Lu-Godt = %.5f em z = %.4f m | erro %+.3f%% | banda central z = %.3f m (%.3f a "
                "%.3f) | %.1f s\n",
                fsLEM, zLEM, 100. * (fs / fsLEM - 1.), zb, zmin, zmax, seg);
    if (o.pressao == "perfil")
        std::printf("  FS_min do arquivo (SlopeRainColumn): %.5f em z = %.4f m\n", perfil.fsArquivo, perfil.zArquivo);

    if (!o.vtk.empty()) f.a->PostProcess(o.vtk, {"WaterPressure"}, {"Displacement"});

    const bool novo = !std::ifstream(o.csv).good();
    std::ofstream csv(o.csv, std::ios::app);
    if (novo)
        csv << "caso,LH,ny,nref,elementos,equacoes,FS_SRM,FS_LEM,z_LEM,z_banda,z_banda_min,z_banda_max,erro_rel,"
               "segundos\n";
    char lin[512];
    std::snprintf(lin, sizeof(lin), "%s,%g,%d,%d,%lld,%lld,%.6f,%.6f,%.4f,%.4f,%.4f,%.4f,%.6f,%.1f\n", o.caso.c_str(),
                  o.LH, o.ny, o.nref, (long long)nel, (long long)neq, fs, fsLEM, zLEM, zb, zmin, zmax, fs / fsLEM - 1.,
                  seg);
    csv << lin;
    return 0;
}
