//
// aterro.cpp
//
// Exemplo da Itasca (FLAC3D) "Embankment Loading on a Cam-Clay Foundation" com a estrutura nativa do NeoPZ:
// TPZMultiphysicsCompMesh (u: H1 vetorial com memória, p: H1 linear), material TPZMatPoroElastoPlastic3DMem
// com o Cam-Clay modificado (TPZModifiedCamClay), TPZLinearAnalysis com Newton, saída VTK por
// DefineGraphMesh/PostProcess e TPZPostProcAnalysis (ver PoroCamClayNativo.h). Porte de aterro_itasca.py /
// aterro_elementos.py. Unidades: kPa, m, s.
//
// Fatia de 1 m (deformação plana) com meia simetria: 20 x 1 x 10 m, malha 20 x 1 x 10 (como no FLAC3D).
// Fundação Cam-Clay modificado: ν = 0.3, M = 0.888, λ = 0.161, κ = 0.062, p_ref = 1 kPa, v_λ = 2.858,
// p'c0 = 160 kPa; densidade seca 2000 kg/m³, n = 0.3 (γ_sat = 23 kN/m³); fluido: K_f = 2·10⁵ kPa
// (α = 1, 1/M = n/K_f), mobilidade k = 10⁻⁹ m²/(kPa·s). G hipoelástico com ν constante (como o FLAC3D).
// Estado inicial: N.A. na superfície; σ_h = 0.7 σ_v em tensões totais (σ'h/σ'v ≈ 6/13), avaliado em
// cada ponto de integração. Topo drenante (p = 0), base fixa e impermeável, u_x = 0 em x = 0 e x = 20.
// Etapa 1: sobrecarga de 50 kPa em 0 <= x <= 4 m sem fluxo (não drenada), em 10 incrementos.
// Etapa 2: adensamento acoplado até t = 10⁸ s (25 passos em escala logarítmica).
//
// Uso:  AterroCamClay [hex20|hex8 ...]    (padrão: hex20 e hex8)
//   hex8   hexaedro de 8 nós, u e p lineares (Q1-Q1, o mesmo elemento do Python)
//   hex20  hexaedro de 20 nós (TPZQuadraticCube), u quadrático H1 do NeoPZ (Q2, 27 funções), p linear,
//          integração completa 3x3x3 (com 2x2x2 o Q2 tem modos espúrios e a matriz fica singular; o hex20
//          serendipity com integração reduzida está em serendipity/aterro_serendipity.cpp)
// Saída: aterro_<elem>_up.scal_vec.N.vtk (u, p), aterro_<elem>_tensoes.scal_vec.N.vtk (σ', σ, p', q, p_c, ...),
//        aterro_<elem>_malha.vtk (TPZVTKGeoMesh), aterro_<elem>_historico.csv
//
#include <cmath>
#include <cstdio>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "PoroCamClayNativo.h"
#include "TPZVTKGeoMesh.h"
#include "pzgeoel.h"
#include "malhas.h"

namespace {

const REAL LX = 20., LY = 1., LZ = 10.;
const int NX = 20, NY = 1, NZ = 10;
const REAL NU = 0.3, M_CS = 0.888, LAM = 0.161, KAP = 0.062, VLAM = 2.858, PC0 = 160.;
const REAL RHO_DRY = 2.0, POR = 0.3, RHO_W = 1.0, GRAV = 10.;  // t/m³, m/s² -> kN/m³ = kPa/m
const REAL GAM_SAT = (RHO_DRY + POR * RHO_W) * GRAV;             // 23 kPa/m
const REAL GAM_W = RHO_W * GRAV;                                 // 10 kPa/m
const REAL K0R = 0.7;                                            // σh/σv (totais)
const REAL KF = 2.0e5;                                           // kPa
const REAL MB = KF / POR;                                        // módulo de Biot
const REAL PERM = 1.0e-9;                                        // m²/(kPa·s)
const REAL LOAD = 50., XLOAD = 4.;
const int MAT_TOPO = 6, MAT_TOPO_CARREGADO = 7;

/// σ' inicial (tração +): σ_h = 0.7 σ_v em tensões totais
TPZTensor<REAL> InitialEffectiveStress(REAL z) {
    const REAL d = LZ - z;
    const REAL svTot = -GAM_SAT * d, pw = GAM_W * d;
    TPZTensor<REAL> s;
    s[_XX_] = s[_YY_] = K0R * svTot + pw;
    s[_ZZ_] = svTot + pw;
    return s;
}

/// Parâmetros do Cam-Clay no ponto x (v0 da NCL + linha κ com o p'0 local)
TPZModifiedCamClay ParAt(const TPZVec<REAL> &x) {
    const TPZTensor<REAL> s0 = InitialEffectiveStress(x[2]);
    const REAL p0 = -(s0[_XX_] + s0[_YY_] + s0[_ZZ_]) / 3.;
    const REAL v0 = VLAM - LAM * std::log(PC0) + KAP * std::log(PC0 / p0);
    TPZModifiedCamClay m;
    m.SetUp(M_CS, LAM, KAP, VLAM, v0, PC0, p0, 0., 1., TPZModifiedCamClay::EPressureDependent,
            TPZModifiedCamClay::EHypoNu, 0., NU);
    m.SetInitialStress(s0);
    return m;
}

REAL Hydrostatic(const TPZVec<REAL> &x) { return GAM_W * (LZ - x[2]); }

/// FLAC3D (Figs. 8 e 9 da Itasca, digitalizadas): recalques (m, positivos para baixo) e poropressões (kPa)
struct Row {
    REAL t, uz0, uz2, uz4, uz6, pp1, pp2;
};
const Row kFlac[] = {{0.0, 0.1400, 0.1350, 0.0550, -0.0416, 18.10, 62.40},
                     {2.5e5, 0.1414, 0.1353, 0.0565, -0.0409, 10.06, 44.75},
                     {1.0e6, 0.1470, 0.1404, 0.0605, -0.0375, 7.45, 35.96},
                     {1.0e7, 0.1700, 0.1631, 0.0825, -0.0170, 5.71, 28.44},
                     {3.0e7, 0.1862, 0.1790, 0.0979, -0.0016, 5.19, 25.83},
                     {1.0e8, 0.1930, 0.1855, 0.1041, 0.0045, 5.05, 25.07}};
/// FE Python (fe3d_up.py / aterro_elementos.py): hex8 (Q1-Q1) e hex20 serendipity (Q2-Q1)
const Row kPyHex8[] = {{0.0, 0.1362, 0.1377, 0.0426, -0.0469, 49.87, 53.76},
                       {2.5e5, 0.1636, 0.1630, 0.0675, -0.0364, 27.44, 58.43},
                       {1.0e6, 0.1883, 0.1862, 0.0838, -0.0339, 16.76, 56.42},
                       {1.0e7, 0.2437, 0.2373, 0.1265, -0.0071, 6.50, 31.69},
                       {3.0e7, 0.2644, 0.2575, 0.1463, 0.0132, 5.37, 26.72},
                       {1.0e8, 0.2730, 0.2660, 0.1546, 0.0222, 5.02, 25.11}};
const Row kPyHex20[] = {{0.0, 0.1528, 0.1523, 0.0673, -0.0386, 33.06, 55.72},
                        {2.5e5, 0.1673, 0.1675, 0.0765, -0.0372, 26.96, 58.47},
                        {1.0e6, 0.1910, 0.1897, 0.0913, -0.0349, 17.05, 55.82},
                        {1.0e7, 0.2457, 0.2404, 0.1357, -0.0069, 6.51, 31.72},
                        {3.0e7, 0.2664, 0.2608, 0.1558, 0.0136, 5.37, 26.72},
                        {1.0e8, 0.2751, 0.2693, 0.1642, 0.0227, 5.02, 25.11}};

REAL InterpLinear(const std::vector<REAL> &x, const std::vector<REAL> &y, REAL t) {
    if (t <= x.front()) return y.front();
    if (t >= x.back()) return y.back();
    for (size_t i = 1; i < x.size(); i++)
        if (t <= x[i]) return y[i - 1] + (y[i] - y[i - 1]) * (t - x[i - 1]) / (x[i] - x[i - 1]);
    return y.back();
}

void Run(const std::string &elem) {
    const int orderU = (elem == "hex8") ? 1 : 2;
    const std::string nome = "aterro_" + elem;
    std::cout << "\n" << std::string(100, '=') << "\nAterro sobre fundação Cam-Clay (FLAC3D), " << elem << ", malha " << NX
              << "x" << NY << "x" << NZ << " — TPZMultiphysicsCompMesh + TPZMatPoroElastoPlastic3DMem\n";
    TPZGeoMesh *gmesh = BoxMesh({LX, LY, LZ}, {NX, NY, NZ}, orderU == 2);
    std::unique_ptr<TPZGeoMesh> gmeshOwner(gmesh);  // destruída depois das malhas computacionais
    // faces do topo com centro em x <= 4 ('range position-x 0 4' do FLAC3D) recebem a sobrecarga
    for (TPZGeoEl *gel : gmesh->ElementVec()) {
        if (!gel || gel->Dimension() != 2 || gel->MaterialId() != MAT_TOPO) continue;
        REAL xm = 0.;
        for (int i = 0; i < gel->NNodes(); i++) xm += gel->NodePtr(i)->Coord(0) / gel->NNodes();
        if (xm <= XLOAD + 1.e-9) gel->SetMaterialId(MAT_TOPO_CARREGADO);
    }
    {
        std::ofstream f(nome + "_malha.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, f, true);
    }

    TPoroCamClayNativo::TParams par;
    par.orderU = orderU;
    TPZManVector<REAL, 3> xc = {LX / 2, LY / 2, LZ / 2};
    par.model = ParAt(xc);
    par.modelUpdate = [](const TPZVec<REAL> &x, TPZModifiedCamClay &m) { m = ParAt(x); };
    par.p0 = Hydrostatic;
    par.hydro = Hydrostatic;
    par.alpha = 1.;
    par.Se = 1. / MB;
    par.k = PERM;
    par.mu = 1.;
    par.rhof = RHO_W;
    par.g = {0., 0., -GRAV};
    par.body = {0., 0., -GAM_SAT};
    par.vtkResolution = (orderU == 2) ? 1 : 0;
    using M = TMatUP;
    par.bcs = {{1, M::EDirectionalNullU, {0, 0, 0}, {1, 0, 0, 0}},           // u_x = 0 em x = 0
               {2, M::EDirectionalNullU, {0, 0, 0}, {1, 0, 0, 0}},           // u_x = 0 em x = 20
               {3, M::EDirectionalNullU, {0, 0, 0}, {0, 1, 0, 0}},           // u_y = 0 (fatia)
               {4, M::EDirectionalNullU, {0, 0, 0}, {0, 1, 0, 0}},
               {5, M::EDirectionalNullU, {0, 0, 0}, {1, 1, 1, 0}},           // base fixa
               {MAT_TOPO, M::ETractionDirichletP, {0, 0, 0}, {0, 0, 0, 0}},  // topo drenante (p = 0)
               {MAT_TOPO_CARREGADO, M::ETractionDirichletP, {0, 0, 0}, {0, 0, 0, 0}}};
    TPoroCamClayNativo prob(gmesh, par);
    prob.DefineVTK(nome);
    std::printf("  %lld equações (u e p); max |R_u(t = 0)| nos gdl livres = %.2e kN\n",
                (long long)prob.NEquations(), prob.FreeResidualU());

    // pontos monitorados (como no FLAC3D)
    auto nodeAt = [&](REAL x, REAL y, REAL z) {
        int64_t best = 0;
        REAL dmin = 1.e300;
        for (int64_t i = 0; i < gmesh->NNodes(); i++) {
            const auto &n = gmesh->NodeVec()[i];
            const REAL d = std::hypot(std::hypot(n.Coord(0) - x, n.Coord(1) - y), n.Coord(2) - z);
            if (d < dmin) { dmin = d; best = i; }
        }
        return best;
    };
    auto zoneAt = [&](REAL x, REAL z) -> TPZGeoEl * {
        for (TPZGeoEl *gel : gmesh->ElementVec()) {
            if (!gel || gel->Dimension() != 3) continue;
            REAL x0 = 1.e300, x1 = -1.e300, z0 = 1.e300, z1 = -1.e300;
            for (int a = 0; a < 8; a++) {
                const TPZGeoNode *n = gel->NodePtr(a);
                x0 = std::min(x0, n->Coord(0)); x1 = std::max(x1, n->Coord(0));
                z0 = std::min(z0, n->Coord(2)); z1 = std::max(z1, n->Coord(2));
            }
            if (x0 <= x && x <= x1 && z0 <= z && z <= z1) return gel;
        }
        return nullptr;
    };
    const REAL xs[4] = {0., 2., 4., 6.};
    int64_t mon[4];
    for (int i = 0; i < 4; i++) mon[i] = nodeAt(xs[i], 0., LZ);
    TPZGeoEl *zpp1 = zoneAt(0.5, 9.5), *zpp2 = zoneAt(1.5, 7.5);

    struct Rec {
        REAL t, lam, uz[4], pp1, pp2;
    };
    std::vector<Rec> hist;
    int nvtk = 0;
    auto record = [&](REAL t, REAL lam) {
        Rec r;
        r.t = t;
        r.lam = lam;
        for (int i = 0; i < 4; i++) r.uz[i] = prob.NodalU(mon[i], 2);
        r.pp1 = prob.ZonePressure(zpp1);
        r.pp2 = prob.ZonePressure(zpp2);
        hist.push_back(r);
        prob.WriteVTK(nvtk++);
    };
    record(0., 0.);

    int ncut = 0;
    std::function<int(REAL, REAL, REAL, REAL, bool, int)> advance =
        [&](REAL lam0, REAL lam1, REAL ta, REAL tb, bool flow, int level) -> int {
        prob.SetBCVal2(MAT_TOPO_CARREGADO, {0., 0., -LOAD * lam1, 0.});
        try {
            return prob.Step(std::max(tb - ta, REAL(1.)), flow, false);
        } catch (std::exception &) {
            if (level >= 10) throw;
            ncut++;
            const REAL lm = 0.5 * (lam0 + lam1), tm = 0.5 * (ta + tb);
            return advance(lam0, lm, ta, tm, flow, level + 1) + advance(lm, lam1, tm, tb, flow, level + 1);
        }
    };
    // ---- etapa 1: carga não drenada
    const int nload = 10;
    for (int k = 1; k <= nload; k++) {
        const int it = advance(REAL(k - 1) / nload, REAL(k) / nload, 0., 0., false, 0);
        record(0., REAL(k) / nload);
        std::printf("  carga %5.2f: %2d it.  uz(0) = %8.4f m\n", REAL(k) / nload, it, hist.back().uz[0]);
    }
    const Rec und = hist.back();
    // ---- etapa 2: adensamento
    REAL t = 0.;
    for (int k = 0; k < 25; k++) {
        const REAL tn = std::pow(10., 2. + 6. * k / 24.);
        const int it = advance(1., 1., t, tn, true, 0);
        t = tn;
        record(t, 1.);
        std::printf("  t = %9.3e s: %2d it.  uz(0) = %8.4f m   pp1 = %7.3f  pp2 = %7.3f kPa\n", t, it,
                    hist.back().uz[0], hist.back().pp1, hist.back().pp2);
    }
    std::printf("  subdivisões de passo: %d;  %d saídas VTK (%s_up.scal_vec.N.vtk e %s_tensoes.scal_vec.N.vtk)\n",
                ncut, nvtk, nome.c_str(), nome.c_str());

    // ---- históricos e comparação
    {
        std::ofstream f(nome + "_historico.csv");
        f << std::setprecision(10) << "t,fator_carga,uz0,uz2,uz4,uz6,pp1,pp2\n";
        for (auto &r : hist)
            f << r.t << "," << r.lam << "," << r.uz[0] << "," << r.uz[1] << "," << r.uz[2] << "," << r.uz[3] << ","
              << r.pp1 << "," << r.pp2 << "\n";
    }
    std::vector<REAL> tt, v[6];
    for (size_t i = nload; i < hist.size(); i++) {
        tt.push_back(hist[i].t);
        for (int j = 0; j < 4; j++) v[j].push_back(-hist[i].uz[j]);
        v[4].push_back(hist[i].pp1);
        v[5].push_back(hist[i].pp2);
    }
    const Row *py = (elem == "hex8") ? kPyHex8 : kPyHex20;
    const char *pyname = (elem == "hex8") ? "Python hex8" : "Python hex20s";
    std::printf("\n  Recalques (-u_z, m) e poropressões de zona (kPa): NeoPZ x FE Python x FLAC3D\n");
    std::printf("  (hex20s = hex20 serendipity do Python; o hex20 do NeoPZ usa u quadrático Q2)\n");
    std::printf("  %-14s%10s%9s%9s%9s%9s%8s%8s\n", "", "t (s)", "uz0", "uz2", "uz4", "uz6", "pp1", "pp2");
    for (int k = 0; k < 6; k++) {
        const REAL tk = kFlac[k].t;
        std::printf("  %-14s%10.1e", ("NeoPZ " + elem).c_str(), tk);
        for (int j = 0; j < 4; j++) std::printf("%9.4f", InterpLinear(tt, v[j], tk));
        for (int j = 4; j < 6; j++) std::printf("%8.2f", InterpLinear(tt, v[j], tk));
        std::printf("\n");
        for (const auto &rr : {std::make_pair(pyname, py[k]), std::make_pair("FLAC3D", kFlac[k])}) {
            const Row &r = rr.second;
            std::printf("  %-14s%10.1e%9.4f%9.4f%9.4f%9.4f%8.2f%8.2f\n", rr.first, tk, r.uz0, r.uz2, r.uz4, r.uz6,
                        r.pp1, r.pp2);
        }
        std::printf("\n");
    }
    std::printf("  fim da etapa não drenada: recalque em x = 0: %.4f m (FLAC3D ~0.14); final: %.4f m (FLAC3D ~0.19)\n",
                -und.uz[0], -hist.back().uz[0]);
}

} // namespace

int main(int argc, char **argv) {
    std::vector<std::string> tipos;
    for (int i = 1; i < argc; i++) {
        const std::string a = argv[i];
        if (a != "hex8" && a != "hex20") {
            std::cerr << "uso: " << argv[0] << " [hex20|hex8 ...]\n";
            return 1;
        }
        tipos.push_back(a);
    }
    if (tipos.empty()) tipos = {"hex20", "hex8"};
    for (auto &et : tipos) Run(et);
    return 0;
}
