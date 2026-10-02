//
// triaxial_abaqus.cpp
//
// Abaqus Benchmarks Manual 1.15.2 "Consolidation of a triaxial test specimen" com o FE u-p 3D do NeoPZ
// (TPZPoroCamClayUP), hexaedros de 20 nós (u quadrático serendipity, p trilinear nos vértices: Q2-Q1, o
// mesmo par do C3D20P / C3D20RP do Abaqus) e o Cam-Clay modificado (TPZModifiedCamClay). Porte de
// triaxial_abaqus.py.
//
// Corpo de prova cilíndrico, altura/diâmetro = 3, metade superior modelada (simetria no plano médio),
// H = 60 mm, r = 20 mm; aqui em 3D: um quarto do cilindro (simetria em x = 0 e y = 0), malha "O-grid" com
// os nós de meio de aresta da superfície lateral sobre o arco. Cam-Clay modificado com elasticidade porosa:
// ν = 0.3, κ = 0.026, λ = 0.174, M = 1.0, a0 = 58.3 kPa (p_c0 = 2 a0), e0 = 1.08 (v0 = 2.08);
// k = 1.728e-4 m/dia, γ_w = 10 kN/m³; tensão efetiva inicial isotrópica de 100 kPa, pressão confinante
// P = 100 kPa constante na face lateral; a placa desce até δ/H = 0.6 em 400 dias com drenagem livre (p = 0)
// no topo. Placa lisa: só u_z prescrito no topo; placa rugosa: também u_x = u_y = 0 no topo. Pequenas
// deformações. Resultados no ponto A (r = 5 mm, z = 7.5 mm), interpolados dos pontos de Gauss.
// Unidades: m, kPa, kN, s.
//
// Versão de referência com o elemento serendipity de 20 nós (TPZPoroCamClayUP, mesmo elemento do Python e do
// C3D20P do Abaqus); a versão com a estrutura nativa do NeoPZ é ../triaxial_abaqus.cpp.
//
// Uso:  TriaxialAbaqusCamClaySerendipity                         (lisa e rugosa, hex20 e hex20r, malha 2x2x4)
//       TriaxialAbaqusCamClaySerendipity rugosa hex20 3 3 8 150   (placa, elemento, nc nr nz, passos)
// Saída: serendipity_triaxial_<placa>_<elem>_<malha>_NNN.vtk e .csv, serendipity_triaxial_homogeneo.csv
//
#include <cmath>
#include <cstdio>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "TPZPoroCamClayUP.h"
#include "malhas.h"

using TUP = TPZPoroCamClayUP;
using TMCC = TPZModifiedCamClay;

namespace {

const REAL R = 0.020, H = 0.060;                 // raio e meia altura (m)
const REAL P_CONF = 100.0;                       // pressão confinante = tensão efetiva inicial (kPa)
const REAL V0 = 2.08;                            // 1 + e0
const REAL KAP = 0.026, LAM = 0.174, M_CS = 1.0, A0 = 58.3, NU = 0.3;
const REAL K_HYD = 1.728e-4 / 86400.0;           // condutividade hidráulica (m/s)
const REAL GAM_W = 10.0;                         // kN/m³
const REAL MOBILITY = K_HYD / GAM_W;             // k/γ_w (m²/(kPa s))
const REAL T_END = 34.56e6;                      // 400 dias (s)
const REAL DH_END = 0.6;                         // δ/H no fim
const REAL E_LIN = 15.0e3;                       // E da variante com elasticidade linear (kPa)
const std::array<REAL, 3> POINT_A = {0.005, 0.0, 0.0075};

/// Cam-Clay do benchmark. porous: p = p0 exp(-v0 ε_v^e/κ), G de ν na forma incremental (elasticidade
/// porosa do Abaqus); linear: E = 15 MPa, ν = 0.3.
TMCC Parameters(bool porous) {
    TMCC m;
    if (porous) {
        m.SetUp(M_CS, LAM, KAP, 0., V0, 2. * A0, P_CONF, 0., 1., TMCC::EPressureDependent, TMCC::EHypoNu, 0., NU);
    } else {
        m.SetUp(M_CS, LAM, KAP, 0., V0, 2. * A0, P_CONF, 0., 1., TMCC::ELinear, TMCC::EConstantG,
                E_LIN / (2. * (1. + NU)), NU);
        m.SetLinearBulkModulus(E_LIN / (3. * (1. - 2. * NU)));
    }
    return m;
}

/// Pontos digitalizados do Abaqus (Figs. 1.15.2-2 e -3): q (kPa) x δ/H
const REAL kAbaqusRugosa[11][2] = {{0.0304, 59.69}, {0.08, 91.89},   {0.1294, 111.84}, {0.1688, 125.81},
                                   {0.2201, 134.73}, {0.2697, 140.85}, {0.3191, 144.98}, {0.3702, 147.78},
                                   {0.4097, 150.08}, {0.4593, 152.19}, {0.5096, 153.09}};
const REAL kAbaqusLisa[11][2] = {{0.0297, 59.68}, {0.0791, 89.96},  {0.1301, 109.94}, {0.1695, 122.81},
                                 {0.2189, 131.09}, {0.27, 137.23},   {0.3192, 141.09}, {0.3687, 144.17},
                                 {0.4081, 145.85}, {0.4593, 147.25}, {0.5086, 148.08}};
/// FE Python (triaxial_abaqus.py, malha 2x2x4, 150 passos) nos mesmos δ/H
struct PyRef {
    const char *placa, *elem;
    REAL q[11];
};
const PyRef kPython[] = {
    {"lisa", "hex20", {60.2, 93.3, 114.0, 124.7, 133.7, 139.8, 143.5, 145.8, 147.1, 148.2, 148.8}},
    {"lisa", "hex20r", {60.2, 93.3, 114.0, 124.7, 133.8, 139.8, 143.5, 145.8, 147.1, 148.2, 148.8}},
    {"rugosa", "hex20", {61.6, 95.9, 117.2, 128.7, 138.5, 144.1, 146.9, 147.9, 147.7, 146.7, 145.2}},
    {"rugosa", "hex20r", {61.5, 95.6, 116.4, 127.4, 136.9, 142.9, 146.9, 149.9, 151.6, 153.3, 154.5}}};

/// Completa com espaços até w caracteres (contando caracteres UTF-8, não bytes)
std::string Pad(const std::string &s, size_t w) {
    size_t n = 0;
    for (unsigned char c : s)
        if ((c & 0xC0) != 0x80) n++;
    return s + std::string(n < w ? w - n : 0, ' ');
}

REAL Interp(const std::vector<REAL> &x, const std::vector<REAL> &y, REAL t) {
    if (t <= x.front()) return y.front();
    if (t >= x.back()) return y.back();
    for (size_t i = 1; i < x.size(); i++)
        if (t <= x[i]) return y[i - 1] + (y[i] - y[i - 1]) * (t - x[i - 1]) / (x[i] - x[i - 1]);
    return y.back();
}

/// Placa lisa = estado homogêneo: ε_zz = -δ/H prescrita, σ_xx = σ_yy = -P, distorções nulas (Newton nas
/// deformações laterais, com σ_n e ε_n do passo anterior para a lei incremental). Colunas: δ/H, p', q, ε_v
std::vector<std::array<REAL, 4>> Homogeneous(bool porous, int nsteps = 600) {
    const TMCC P = Parameters(porous);
    TPZTensor<REAL> eps, epsp, sig = P.InitialStress();
    REAL al = 0.;
    std::vector<std::array<REAL, 4>> out = {{0., P_CONF, 0., 0.}};
    REAL ratio = NU / (1. - NU);
    TMCC::TResult r;
    for (int k = 1; k <= nsteps; k++) {
        const REAL dea = DH_END / nsteps;
        TPZTensor<REAL> e(eps);
        e[_ZZ_] -= dea;
        e[_XX_] += ratio * dea;
        e[_YY_] += ratio * dea;
        for (int it = 0; it < 50; it++) {
            P.ReturnMapping(e, epsp, al, r, &sig, &eps);
            const REAL r0 = r.stress[_XX_] + P_CONF, r1 = r.stress[_YY_] + P_CONF;
            if (std::max(std::fabs(r0), std::fabs(r1)) < 1.e-10 * P_CONF) break;
            const REAL a = r.Dep(_XX_, _XX_), b = r.Dep(_XX_, _YY_), c = r.Dep(_YY_, _XX_), d = r.Dep(_YY_, _YY_);
            const REAL det = a * d - b * c;
            e[_XX_] += (-r0 * d + r1 * b) / det;
            e[_YY_] += (-r1 * a + r0 * c) / det;
        }
        ratio = (e[_XX_] - eps[_XX_]) / dea;
        eps = e;
        epsp = r.plastic_strain;
        al = r.alpha;
        sig = r.stress;
        REAL p, q;
        TMCC::Invariants(sig, p, q);
        out.push_back({DH_END * k / nsteps, -p, q, -(e[_XX_] + e[_YY_] + e[_ZZ_])});
    }
    return out;
}

struct Resultado {
    std::vector<REAL> dh, p, q, sa, pmax;
    std::vector<int> it;
    int cuts = 0;
    REAL res0 = 0.;
};

Resultado Run(bool rugosa, TUP::EElement et, int nc, int nr, int nz, int nsteps, bool porous = true) {
    const std::string placa = rugosa ? "rugosa" : "lisa";
    char malha[32];
    std::snprintf(malha, sizeof(malha), "%d%d%d", nc, nr, nz);
    const std::string nome = "serendipity_triaxial_" + placa + "_" + TUP::ElementName(et) + "_" + malha + (porous ? "" : "_linear");
    TPZGeoMesh *gmesh = QuarterCylinderMesh(R, H, nc, nr, nz, et != TUP::EHex8);
    const TMCC par = Parameters(porous);
    TUP::TCoupling cp;
    cp.alpha = 1.;
    cp.invBiotModulus = 0.;
    cp.perm = MOBILITY;
    TUP model(gmesh, et, [&](const TPZVec<REAL> &) { return par; }, cp, nullptr, false);

    const REAL tol = 1.e-9 * R;
    std::map<int64_t, REAL> fixed0, fixedP;
    std::vector<int64_t> top;
    for (int64_t i = 0; i < model.NNodes(); i++) {
        const auto &c = model.Coord(i);
        if (std::fabs(c[0]) < tol) fixed0[3 * i] = 0.;
        if (std::fabs(c[1]) < tol) fixed0[3 * i + 1] = 0.;
        if (std::fabs(c[2]) < tol) fixed0[3 * i + 2] = 0.;
        if (std::fabs(c[2] - H) < tol) top.push_back(i);
    }
    if (rugosa)
        for (int64_t nd : top) fixed0[3 * nd] = fixed0[3 * nd + 1] = 0.;
    for (int64_t nd : top)
        if (model.PDof(nd) >= 0) fixedP[model.PDof(nd)] = 0.;
    const std::vector<REAL> Fp = model.FacePressure([](const TUP::TFace &f) { return f.marker == ELateral; }, P_CONF,
                                                    [](const TPZVec<REAL> &x) { return std::array<REAL, 3>{x[0], x[1], 0.}; });
    auto fixedAt = [&](REAL dh) {
        std::map<int64_t, REAL> fx(fixed0);
        for (int64_t nd : top) fx[3 * nd + 2] = -dh * H;
        return fx;
    };
    int64_t kA;
    std::array<REAL, 3> xiA;
    if (!model.Locate(POINT_A, kA, xiA)) throw std::runtime_error("ponto A fora da malha");
    auto platen = [&]() {
        const std::vector<REAL> Rv = model.Reactions(Fp);
        REAL s = 0.;
        for (int64_t nd : top) s += Rv[3 * nd + 2];
        return -s / (0.25 * M_PI * R * R);
    };

    Resultado res;
    {   // etapa geostática: equilíbrio do estado inicial com a pressão confinante
        const std::map<int64_t, REAL> fx = fixedAt(0.);
        const std::vector<REAL> Rv = model.Reactions(Fp);
        for (size_t i = 0; i < Rv.size(); i++)
            if (!fx.count(int64_t(i))) res.res0 = std::max(res.res0, std::fabs(Rv[i]));
    }
    int it = model.Step(Fp, fixedAt(0.), fixedP, 1., false);
    auto record = [&](REAL dh, int its) {
        REAL p, q;
        TMCC::Invariants(model.StressAt(kA, xiA), p, q);
        REAL pmax = 0.;
        for (REAL v : model.P()) pmax = std::max(pmax, std::fabs(v));
        res.dh.push_back(dh);
        res.p.push_back(-p);
        res.q.push_back(q);
        res.sa.push_back(platen());
        res.pmax.push_back(dh == 0. ? 0. : pmax);
        res.it.push_back(its);
    };
    record(0., it);
    int nvtk = 0;
    auto vtk = [&](REAL dh, REAL t) {
        char fn[256];
        std::snprintf(fn, sizeof(fn), "%s_%03d.vtk", nome.c_str(), nvtk++);
        model.WriteVTK(fn, "Abaqus 1.15.2 - triaxial com Cam-Clay, placa " + placa + ", " + TUP::ElementName(et),
                       {{"delta_H", dh}, {"tempo_s", t}});
    };
    vtk(0., 0.);

    // preditor do 1º passo: campo homogêneo elástico (u_z = -z δ/H, expansão lateral ν)
    std::vector<REAL> rate(3 * model.NNodes());
    for (int64_t i = 0; i < model.NNodes(); i++) {
        const auto &c = model.Coord(i);
        rate[3 * i] = NU * c[0];
        rate[3 * i + 1] = NU * c[1];
        rate[3 * i + 2] = -c[2];
    }
    std::function<int(REAL, REAL, int)> advance = [&](REAL dh0, REAL dh1, int level) -> int {
        const REAL dt = T_END * (dh1 - dh0) / DH_END;
        try {
            const std::vector<REAL> U0 = model.U();
            std::vector<REAL> guess(U0);
            for (size_t i = 0; i < guess.size(); i++) guess[i] += rate[i] * (dh1 - dh0);
            const int its = model.Step(Fp, fixedAt(dh1), fixedP, dt, true, 1.e-8, 30, &guess);
            for (size_t i = 0; i < rate.size(); i++) rate[i] = (model.U()[i] - U0[i]) / (dh1 - dh0);
            return its;
        } catch (std::exception &) {
            if (level >= 8) throw;
            res.cuts++;
            const REAL dm = 0.5 * (dh0 + dh1);
            return advance(dh0, dm, level + 1) + advance(dm, dh1, level + 1);
        }
    };
    for (int k = 1; k <= nsteps; k++) {
        const REAL dh0 = DH_END * (k - 1) / nsteps, dh = DH_END * k / nsteps;
        const int its = advance(dh0, dh, 0);
        record(dh, its);
        if (k % 10 == 0 || k == nsteps) vtk(dh, T_END * dh / DH_END);
    }
    REAL itmean = 0., pmax = 0.;
    int itmax = 0;
    for (size_t i = 1; i < res.it.size(); i++) {
        itmean += res.it[i] / REAL(res.it.size() - 1);
        itmax = std::max(itmax, res.it[i]);
        pmax = std::max(pmax, res.pmax[i]);
    }
    std::printf("  %-6s placa %-6s malha nc=%d nr=%d nz=%d: %lld nós, %lld elementos; resíduo inicial %.1e kN; Newton %.1f "
                "it/passo (máx %d), subdivisões %d; |p_poro| máx %.4f kPa; fim: p' = %.1f, q = %.1f kPa (A), "
                "σ_a(placa) = %.1f kPa; %d VTK\n",
                TUP::ElementName(et).c_str(), placa.c_str(), nc, nr, nz, (long long)model.NNodes(),
                (long long)model.NElements(), res.res0, itmean, itmax, res.cuts, pmax, res.p.back(), res.q.back(),
                res.sa.back(), nvtk);
    std::ofstream f(nome + ".csv");
    f << std::setprecision(10) << "dH,p_A,q_A,sigma_a_placa,pmax,iteracoes\n";
    for (size_t i = 0; i < res.dh.size(); i++)
        f << res.dh[i] << "," << res.p[i] << "," << res.q[i] << "," << res.sa[i] << "," << res.pmax[i] << ","
          << res.it[i] << "\n";
    delete gmesh;
    return res;
}

void Tabela(const std::string &titulo, const REAL abq[11][2], const std::vector<std::pair<std::string, Resultado>> &runs,
            const std::vector<std::array<REAL, 4>> *hom, const std::string &placa) {
    std::printf("\n  %s: q (kPa) nos pontos digitalizados do Abaqus\n", titulo.c_str());
    std::printf("  %s", Pad("δ/H", 34).c_str());
    for (int i = 0; i < 11; i++) std::printf("%7.3f", abq[i][0]);
    std::printf("\n  %s", Pad("Abaqus (CAX8RP), digitalizado", 34).c_str());
    for (int i = 0; i < 11; i++) std::printf("%7.1f", abq[i][1]);
    std::printf("\n");
    if (hom) {
        std::vector<REAL> x, y;
        for (auto &r : *hom) { x.push_back(r[0]); y.push_back(r[2]); }
        std::printf("  %s", Pad("solução homogênea (600 inc.)", 34).c_str());
        for (int i = 0; i < 11; i++) std::printf("%7.1f", Interp(x, y, abq[i][0]));
        std::printf("\n");
    }
    for (auto &rr : runs) {
        std::printf("  %s", Pad("NeoPZ " + rr.first, 34).c_str());
        for (int i = 0; i < 11; i++) std::printf("%7.1f", Interp(rr.second.dh, rr.second.q, abq[i][0]));
        std::printf("\n");
        for (auto &py : kPython)
            if (placa == py.placa && rr.first.find(std::string(py.elem) + " ") == 0) {
                std::printf("  %s", Pad("Python " + rr.first, 34).c_str());
                for (int i = 0; i < 11; i++) std::printf("%7.1f", py.q[i]);
                std::printf("\n");
            }
    }
}

} // namespace

int main(int argc, char **argv) {
    std::cout << "Abaqus Benchmarks 1.15.2: adensamento de um corpo de prova triaxial (Cam-Clay modificado, u-p 3D)\n";
    const auto hom = Homogeneous(true);
    {
        std::ofstream f("serendipity_triaxial_homogeneo.csv");
        f << std::setprecision(10) << "dH,p,q,ev\n";
        for (auto &r : hom) f << r[0] << "," << r[1] << "," << r[2] << "," << r[3] << "\n";
    }
    if (argc > 1) {
        const bool rugosa = std::string(argv[1]) != "lisa";
        const std::string el = argc > 2 ? argv[2] : "hex20";
        const TUP::EElement et = el == "hex20r" ? TUP::EHex20R : (el == "hex8" ? TUP::EHex8 : TUP::EHex20);
        const int nc = argc > 3 ? std::atoi(argv[3]) : 2, nr = argc > 4 ? std::atoi(argv[4]) : 2,
                  nz = argc > 5 ? std::atoi(argv[5]) : 4, ns = argc > 6 ? std::atoi(argv[6]) : 150;
        const Resultado r = Run(rugosa, et, nc, nr, nz, ns);
        Tabela(rugosa ? "Placa rugosa, ponto A" : "Placa lisa", rugosa ? kAbaqusRugosa : kAbaqusLisa,
               {{el + " " + std::to_string(nc) + "x" + std::to_string(nr) + "x" + std::to_string(nz), r}},
               rugosa ? nullptr : &hom, rugosa ? "rugosa" : "lisa");
        return 0;
    }
    std::vector<std::pair<std::string, Resultado>> lisa, rugosa;
    lisa.push_back({"hex20 2x2x4", Run(false, TUP::EHex20, 2, 2, 4, 150)});
    lisa.push_back({"hex20r 2x2x4", Run(false, TUP::EHex20R, 2, 2, 4, 150)});
    rugosa.push_back({"hex20 2x2x4", Run(true, TUP::EHex20, 2, 2, 4, 150)});
    rugosa.push_back({"hex20r 2x2x4", Run(true, TUP::EHex20R, 2, 2, 4, 150)});
    Tabela("Placa lisa (estado homogêneo)", kAbaqusLisa, lisa, &hom, "lisa");
    Tabela("Placa rugosa, ponto A", kAbaqusRugosa, rugosa, nullptr, "rugosa");
    return 0;
}
