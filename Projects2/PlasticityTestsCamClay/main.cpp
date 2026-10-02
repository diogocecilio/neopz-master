//
// PlasticityTestsCamClay/main.cpp
//
// Verificação do modelo Cam-Clay modificado do NeoPZ (TPZModifiedCamClay), porte C++ de
// modified_cam_clay.py. Reproduz os testes dos scripts Python:
//
//   1) ponto material: return mapping e tangente consistente x diferenças finitas
//      (verificacao_rs2.py, seção 2), incluindo p_t, β ≠ 1, K = -v0 p/κ e σ0 anisotrópica;
//   2) ensaios triaxiais drenados (controle misto) das Figs. 8.5-8.8 do RS2 x solução analítica
//      em forma fechada (triaxial_drained / triaxial_analytical);
//   3) benchmark da Itasca "Drained and Undrained Triaxial Compression Test on a Cam-Clay Sample"
//      com um hexaedro e TPZMatElastoPlastic<TPZModifiedCamClay> (benchmark_itasca.py). O caso
//      não drenado usa o limite de permeabilidade nula do sistema u-p de Biot: a pressão de poros
//      é condensada no ponto de integração, p = -α M ε_v e σ = σ' + α² M ε_v I;
//   4) elasticidade hipoelástica (shear = hypo_nu, benchmark 1.15.2 do Abaqus) em ponto material
//      e em elementos finitos.
//
// Os valores de referência marcados "Python" foram obtidos com os scripts originais
// (modified_cam_clay.py e poro_camclay_fem.py). As curvas são gravadas em arquivos CSV;
// compara_python.py recalcula os mesmos casos com o código Python e compara.
//
// Convenções: tração positiva, Voigt {xx, xy, xz, yy, yz, zz}, distorções de engenharia.
//
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "TPZBndCondT.h"
#include "pzstepsolver.h"
#include "pzfstrmatrix.h"
#include "TPZMaterialDataT.h"

#include "TPZModifiedCamClay.h"
#include "TPZMatElastoPlastic.h"
#include "TPZElastoPlasticMem.h"
#include "pzelastoplasticanalysis.h"

using TMCC = TPZModifiedCamClay;
using TMatCamClay = TPZMatElastoPlastic<TPZModifiedCamClay, TPZElastoPlasticMem>;

// =====================================================================================
// utilitários
// =====================================================================================
static int gFalhas = 0;

/// Compara valor x referência (|x - ref| <= max(tolabs, tolrel |ref|))
static void Check(const std::string &nome, REAL valor, REAL ref, REAL tolrel, REAL tolabs = 0.) {
    const REAL err = std::fabs(valor - ref);
    const bool ok = err <= std::max(tolabs, tolrel * std::fabs(ref));
    std::printf("   %-52s %15.8g   ref %15.8g   %s\n", nome.c_str(), valor, ref, ok ? "OK" : "FALHOU");
    if (!ok) gFalhas++;
}

/// Verifica que valor <= limite
static void CheckBelow(const std::string &nome, REAL valor, REAL limite) {
    const bool ok = valor <= limite;
    std::printf("   %-52s %15.3e   max %15.3e   %s\n", nome.c_str(), valor, limite, ok ? "OK" : "FALHOU");
    if (!ok) gFalhas++;
}

static TPZTensor<REAL> Tensor(const std::array<REAL, 6> &v) {
    TPZTensor<REAL> t;
    for (int i = 0; i < 6; i++) t[i] = v[i];
    return t;
}

/// Interpolação linear de y(x) (pontos ordenados por x internamente; fora do intervalo: extremos)
static REAL Interp(std::vector<std::array<REAL, 2>> d, REAL x) {
    std::sort(d.begin(), d.end(), [](const std::array<REAL, 2> &a, const std::array<REAL, 2> &b) { return a[0] < b[0]; });
    if (x <= d.front()[0]) return d.front()[1];
    if (x >= d.back()[0]) return d.back()[1];
    auto it = std::upper_bound(d.begin(), d.end(), x, [](REAL v, const std::array<REAL, 2> &a) { return v < a[0]; });
    const auto &b = *it, &a = *(it - 1);
    return a[1] + (b[1] - a[1]) * (x - a[0]) / (b[0] - a[0]);
}

/// Parâmetros da Tabela 8.1 do RS2 (M = 1.2, λ = 0.077, κ = 0.0066, N = 1.788, ν = 0.3, G = 20 MPa), v0 = 1.70
static TMCC RS2(REAL p0, REAL pc0, TMCC::EShear shear, TMCC::EElasticity el = TMCC::ELinear) {
    TMCC m;
    m.SetUp(1.2, 0.077, 0.0066, 1.788, 1.70, pc0, p0, 0., 1., el, shear, 20000., 0.3);
    return m;
}

/// Tangente por diferenças finitas centradas (verificação)
static void TangentFD(const TMCC &m, const TPZTensor<REAL> &eps, const TPZTensor<REAL> &epsp, REAL al,
                      TPZFMatrix<REAL> &D, const TPZTensor<REAL> *sn = nullptr, const TPZTensor<REAL> *en = nullptr,
                      REAL h = 1.e-8) {
    D.Redim(6, 6);
    TMCC::TResult r1, r2;
    for (int j = 0; j < 6; j++) {
        TPZTensor<REAL> e1(eps), e2(eps);
        e1[j] += h;
        e2[j] -= h;
        m.ReturnMapping(e1, epsp, al, r1, sn, en);
        m.ReturnMapping(e2, epsp, al, r2, sn, en);
        for (int i = 0; i < 6; i++) D(i, j) = (r1.stress[i] - r2.stress[i]) / (2. * h);
    }
}

static REAL MaxRelDiff(const TPZFMatrix<REAL> &A, const TPZFMatrix<REAL> &B) {
    REAL dmax = 0., bmax = 0.;
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            dmax = std::max(dmax, std::fabs(A.GetVal(i, j) - B.GetVal(i, j)));
            bmax = std::max(bmax, std::fabs(B.GetVal(i, j)));
        }
    }
    return dmax / bmax;
}

// =====================================================================================
// 1) ponto material
// =====================================================================================
static void TestePontoMaterial() {
    std::cout << "\n" << std::string(100, '=') << "\n1) Ponto material\n";

    // ---- OCR = 5, incremento isocórico grande -> lado seco (verificacao_rs2.py, seção 2)
    std::cout << "  OCR = 5 (p0 = 100, pc0 = 500), eps = {0.006, 0, 0, 0.006, 0, -0.012}\n";
    TMCC P = RS2(100., 500., TMCC::EConstantNu);
    TMCC::TResult r;
    const TPZTensor<REAL> zero;
    P.ReturnMapping(Tensor({0.006, 0., 0., 0.006, 0., -0.012}), zero, 0., r);
    Check("sigma_xx (Python)", r.stress[_XX_], -64.30011619733827, 1.e-10);
    Check("sigma_zz (Python)", r.stress[_ZZ_], -331.0238456825674, 1.e-10);
    Check("alpha (Python)", r.alpha, -2.065723363352569e-03, 1.e-9);
    Check("dgamma (Python)", r.dgamma, 1.2204833056699915e-05, 1.e-9);
    Check("pc (Python)", r.pc, 475.67058865645004, 1.e-10);
    TPZFNMatrix<36, REAL> Dfd;
    TangentFD(P, Tensor({0.006, 0., 0., 0.006, 0., -0.012}), zero, 0., Dfd);
    CheckBelow("tangente consistente x DF (erro relativo)", MaxRelDiff(r.Dep, Dfd), 1.e-7);
    {
        // interface TPZPlasticBase: ApplyLoad (inversa) recupera a deformação de um passo plástico com
        // endurecimento (argila NC; no lado seco, com amolecimento, a inversa não é única)
        TMCC Pnc = RS2(200., 200., TMCC::EConstantNu);
        const TPZTensor<REAL> epsnc = Tensor({-0.002, 0.0005, 0., -0.003, -0.001, -0.01});
        TMCC::TResult rnc;
        Pnc.ReturnMapping(epsnc, zero, 0., rnc);
        TPZTensor<REAL> epsInv;
        Pnc.ApplyLoad(rnc.stress, epsInv);
        REAL err = 0.;
        for (int i = 0; i < 6; i++) err = std::max(err, std::fabs(epsInv[i] - epsnc[i]));
        CheckBelow("ApplyLoad: max |eps - eps_recuperada| (NC, plástico)", err, 1.e-12);
        Check("ApplyLoad: alpha do estado", Pnc.GetState().m_hardening, rnc.alpha, 1.e-8);
    }

    // ---- 200 estados aleatórios: pt = 20, β = 0.6, K = -v0 p/κ, ν constante, σ0 anisotrópica
    TMCC PX;
    PX.SetUp(1.2, 0.077, 0.0066, 1.788, 1.70, 200., 150., 20., 0.6, TMCC::EPressureDependent, TMCC::EConstantNu,
             20000., 0.3);
    PX.SetInitialStress(Tensor({-150., 10., -5., -140., 3., -160.}));
    std::mt19937 gen(1);
    std::uniform_real_distribution<REAL> ue(-0.005, 0.005), up(-0.001, 0.001), ua(0., 0.005);
    int nplast = 0, nbeta = 0;
    REAL worst = 0.;
    for (int k = 0; k < 200; k++) {
        TPZTensor<REAL> e, ep;
        for (int i = 0; i < 6; i++) e[i] = ue(gen);
        for (int i = 0; i < 6; i++) ep[i] = up(gen);
        const REAL al = ua(gen);
        PX.ReturnMapping(e, ep, al, r);
        nplast += r.plastic;
        nbeta += (r.plastic && r.b != 1.);
        TangentFD(PX, e, ep, al, Dfd);
        worst = std::max(worst, MaxRelDiff(r.Dep, Dfd));
    }
    std::cout << "  200 estados aleatórios (pt = 20, beta = 0.6, K = -v0 p/kappa, sigma0 anisotrópica): " << nplast
              << " plásticos (" << nbeta << " com b = beta)\n";
    CheckBelow("tangente consistente x DF (pior erro relativo)", worst, 1.e-6);

    // ---- estados determinísticos (gravados para compara_python.py)
    std::ofstream csv("estados_pontuais.csv");
    csv << std::setprecision(17);
    csv << "caso,k,plastic,b,alpha";
    for (int i = 0; i < 6; i++) csv << ",s" << i;
    for (int i = 0; i < 6; i++) csv << ",ep" << i;
    for (int i = 0; i < 36; i++) csv << ",D" << i;
    csv << "\n";
    auto grava = [&csv](const std::string &caso, int k, const TMCC::TResult &rr) {
        csv << caso << "," << k << "," << int(rr.plastic) << "," << rr.b << "," << rr.alpha;
        for (int i = 0; i < 6; i++) csv << "," << rr.stress[i];
        for (int i = 0; i < 6; i++) csv << "," << rr.plastic_strain[i];
        for (int i = 0; i < 6; i++)
            for (int j = 0; j < 6; j++) csv << "," << rr.Dep(i, j);
        csv << "\n";
    };
    for (int k = 0; k < 12; k++) {
        TPZTensor<REAL> e, ep;
        for (int i = 0; i < 6; i++) e[i] = std::sin(1.3 * k + i) * 0.005;
        for (int i = 0; i < 6; i++) ep[i] = std::cos(0.7 * k + 2 * i) * 0.001;
        const REAL al = 0.0025 * (1. + std::sin(REAL(k)));
        PX.ReturnMapping(e, ep, al, r);
        grava("px", k, r);
    }

    // ---- modo hipoelástico (parâmetros do benchmark 1.15.2 do Abaqus), σ_n e ε_n dados
    TMCC PH;
    PH.SetUp(1.0, 0.174, 0.026, 0., 2.08, 2. * 58.3, 100., 0., 1., TMCC::EPressureDependent, TMCC::EHypoNu, 0., 0.3);
    const TPZTensor<REAL> sig_n = Tensor({-95., 4., -2., -105., 1., -120.});
    const TPZTensor<REAL> eps_n = Tensor({0.001, 0.0005, 0., 0.0008, -0.0003, -0.002});
    const TPZTensor<REAL> epsp_n = Tensor({0.0001, 0., 0.0002, -0.0001, 0., 0.0003});
    const std::array<std::array<REAL, 6>, 3> deps = {{{0.0002, 0.0001, 0., 0.0002, 0., -0.0005},
                                                       {0.003, 0.002, -0.001, 0.003, 0.0005, -0.009},
                                                       {-0.002, 0., 0., -0.002, 0., -0.002}}};
    REAL worsth = 0.;
    for (int k = 0; k < 3; k++) {
        TPZTensor<REAL> e(eps_n);
        for (int i = 0; i < 6; i++) e[i] += deps[k][i];
        PH.ReturnMapping(e, epsp_n, 0.001, r, &sig_n, &eps_n);
        grava("hypo", k, r);
        TangentFD(PH, e, epsp_n, 0.001, Dfd, &sig_n, &eps_n);
        worsth = std::max(worsth, MaxRelDiff(r.Dep, Dfd));
        if (k == 1) {
            Check("hypo: sigma_zz caso plástico (Python)", r.stress[_ZZ_], -131.11032495153646, 1.e-10);
            Check("hypo: alpha caso plástico (Python)", r.alpha, 4.805889498522143e-03, 1.e-9);
        }
    }
    CheckBelow("hypo: tangente consistente x DF (pior erro relativo)", worsth, 1.e-6);
}

// =====================================================================================
// 2) ensaio triaxial drenado (ponto material, controle misto) e solução analítica
// =====================================================================================
struct TriaxState {
    TPZTensor<REAL> eps, epsp, sig;
    REAL alpha = 0.;
    REAL ratio = 0.;
};

/// Um passo: ε_zz prescrita, σ_xx = σ_yy = σ0_xx por Newton nas deformações laterais
static bool TriaxialStep(const TMCC &m, const TriaxState &st, REAL dea, REAL tol, int maxit, TriaxState &out,
                         TMCC::TResult &r, int &its) {
    const REAL sc = m.InitialStress()[_XX_];
    const bool hypo = (m.Shear() == TMCC::EHypoNu);
    TPZTensor<REAL> e(st.eps);
    e[_ZZ_] -= dea;
    e[_XX_] -= st.ratio * dea;  // preditor das deformações laterais
    e[_YY_] -= st.ratio * dea;
    for (int it = 0; it <= maxit; it++) {
        try {
            m.ReturnMapping(e, st.epsp, st.alpha, r, hypo ? &st.sig : nullptr, hypo ? &st.eps : nullptr);
        } catch (TMCC::ReturnMappingError &) {
            return false;
        }
        const REAL r0 = r.stress[_XX_] - sc, r1 = r.stress[_YY_] - sc;
        if (std::max(std::fabs(r0), std::fabs(r1)) <= tol * std::fabs(sc)) {
            out.eps = e;
            out.epsp = r.plastic_strain;
            out.sig = r.stress;
            out.alpha = r.alpha;
            out.ratio = (e[_XX_] - st.eps[_XX_]) / (-dea);
            its = it;
            return true;
        }
        const REAL a = r.Dep(_XX_, _XX_), b = r.Dep(_XX_, _YY_), c = r.Dep(_YY_, _XX_), d = r.Dep(_YY_, _YY_);
        const REAL det = a * d - b * c;
        e[_XX_] += (-r0 * d + r1 * b) / det;
        e[_YY_] += (-r1 * a + r0 * c) / det;
    }
    return false;
}

/// Passo com subdivisão recursiva (bisseção) em caso de falha
static bool TriaxialAdvance(const TMCC &m, const TriaxState &st, REAL dea, int level, REAL tol, int maxit, int maxcut,
                            int &cuts, TriaxState &out, TMCC::TResult &r, int &its) {
    if (TriaxialStep(m, st, dea, tol, maxit, out, r, its)) return true;
    if (level >= maxcut) return false;
    cuts++;
    TriaxState mid;
    int it1, it2;
    if (!TriaxialAdvance(m, st, dea / 2., level + 1, tol, maxit, maxcut, cuts, mid, r, it1)) return false;
    if (!TriaxialAdvance(m, mid, dea / 2., level + 1, tol, maxit, maxcut, cuts, out, r, it2)) return false;
    its = it1 + it2;
    return true;
}

/// Linhas (compressão positiva): ε_a, p', q, ε_v, ε_q, σ_a, p_c
using TRow7 = std::array<REAL, 7>;

static std::vector<TRow7> TriaxialDrained(const TMCC &m, REAL ea_max, int nsteps, REAL tol = 1.e-10, int maxit = 30,
                                          int maxcut = 8) {
    const REAL K0 = m.K0(), G0 = m.G0();
    TriaxState st;
    st.sig = m.InitialStress();
    st.ratio = -(3. * K0 - 2. * G0) / (2. * (3. * K0 + G0));
    REAL p, q;
    TMCC::Invariants(m.InitialStress(), p, q);
    std::vector<TRow7> out;
    out.push_back({0., -p, q, 0., 0., -m.InitialStress()[_ZZ_], m.Pc0()});
    const REAL dea = ea_max / nsteps;
    int cuts = 0;
    for (int k = 0; k < nsteps; k++) {
        TriaxState nxt;
        TMCC::TResult r;
        int its;
        if (!TriaxialAdvance(m, st, dea, 0, tol, maxit, maxcut, cuts, nxt, r, its)) {
            std::cout << "TriaxialDrained: falha no passo " << k + 1 << "\n";
            break;
        }
        st = nxt;
        TMCC::Invariants(st.sig, p, q);
        const TPZTensor<REAL> &e = st.eps;
        out.push_back({-e[_ZZ_], -p, q, -(e[_XX_] + e[_YY_] + e[_ZZ_]), 2. / 3. * std::fabs(e[_ZZ_] - e[_XX_]),
                       -st.sig[_ZZ_], r.pc});
    }
    return out;
}

/// Solução em forma fechada do ensaio triaxial drenado convencional (triaxial_analytical do Python).
/// Válida para p_t = 0 e β = 1. Linhas: ε_a, p', q, ε_v, ε_q, σ_a (compressão positiva).
static std::vector<std::array<REAL, 6>> TriaxialAnalytical(const TMCC &m, int npts, REAL ea_max) {
    const REAL M = m.M(), lam = m.Lambda(), kap = m.Kappa(), v0 = m.V0();
    const REAL p0 = -m.PIni(), pc0 = m.Pc0();
    const REAL gfac = 3. * (1. - 2. * m.Nu()) / (2. * (1. + m.Nu()));
    auto F = [M](REAL x) {
        return (1. / M) * std::log(std::fabs((M + x) / (M - x))) - (2. / M) * std::atan(x / M) -
               std::log(std::fabs(M - x)) / (3. - M) - std::log(M + x) / (3. + M) + 6. * std::log(3. - x) / (9. - M * M);
    };
    auto elastic = [&](REAL pp, REAL q, REAL &ev, REAL &eq) {
        REAL K;
        if (m.Elasticity() == TMCC::ELinear) {
            ev = (pp - p0) / m.K0();
            K = m.K0();
        } else {
            ev = (kap / v0) * std::log(pp / p0);
            K = v0 * pp / kap;
        }
        const REAL G = (m.Shear() == TMCC::EConstantG) ? m.ShearModulusParameter() : gfac * K;
        eq = q / (3. * G);
    };
    const REAL A = 9. + M * M, B = -(18. * p0 + M * M * pc0), C = 9. * p0 * p0;
    const REAL disc = B * B - 4. * A * C;
    REAL py = 1.e300;
    for (REAL root : {(-B + std::sqrt(disc)) / (2. * A), (-B - std::sqrt(disc)) / (2. * A)}) {
        if (root >= p0 - 1.e-9) py = std::min(py, root);
    }
    const REAL qy = 3. * (py - p0), etay = qy / py;
    std::vector<std::array<REAL, 4>> rows;  // p', q, ε_v, ε_q
    REAL ev, eq;
    if (py - p0 > 1.e-9) {
        for (int i = 0; i <= 40; i++) {
            const REAL pp = p0 + (py - p0) * i / 40.;
            const REAL q = 3. * (pp - p0);
            elastic(pp, q, ev, eq);
            rows.push_back({pp, q, ev, eq});
        }
    } else {
        rows.push_back({p0, 0., 0., 0.});
    }
    std::vector<REAL> etas;
    for (int i = 1; i < npts; i++) etas.push_back(etay + (M - etay) * i / npts);
    for (int i = 1; i <= npts; i++) etas.push_back(M - (M - etay) * std::exp(-12. * i / npts));
    std::sort(etas.begin(), etas.end());
    etas.erase(std::unique(etas.begin(), etas.end()), etas.end());
    if (etay > M) std::reverse(etas.begin(), etas.end());
    for (REAL eta : etas) {
        const REAL pp = 3. * p0 / (3. - eta), q = eta * pp, pc = pp * (1. + eta * eta / (M * M));
        elastic(pp, q, ev, eq);
        ev += (lam - kap) / v0 * std::log(pc / pc0);
        eq += (lam - kap) / v0 * (F(eta) - F(etay));
        rows.push_back({pp, q, ev, eq});
    }
    std::vector<std::array<REAL, 6>> out;
    for (auto &R : rows) {
        const REAL ea = R[3] + R[2] / 3.;
        if (ea <= ea_max) out.push_back({ea, R[0], R[1], R[2], R[3], R[0] + 2. * R[1] / 3.});
    }
    return out;
}

static void TesteTriaxialDrenado() {
    std::cout << "\n" << std::string(100, '=') << "\n2) Ensaios triaxiais drenados (400 passos de 0.05% até ea = 20%)"
              << " x solução analítica\n";
    struct Caso {
        std::string nome;
        TMCC P;
        REAL q20, ev20, errmax;  // valores do Python (modified_cam_clay.py)
    };
    std::vector<Caso> casos = {
        {"Fig8.5", RS2(200., 200., TMCC::EConstantNu), 387.8188418845451, 0.05109456243161903, 1.311994743837758},
        {"Fig8.6", RS2(200., 200., TMCC::EConstantG), 387.6047314402047, 0.05107056939734153, 1.2518896064013632},
        {"Fig8.7", RS2(100., 200., TMCC::EConstantNu), 196.73602832340944, 0.023020808022259892, 0.154664569118097},
        {"Fig8.8", RS2(100., 500., TMCC::EConstantNu), 202.99953424604533, -0.013545862154933225, 0.14525841124796557},
        {"OCR5_K_p", RS2(100., 500., TMCC::EConstantNu, TMCC::EPressureDependent), 202.87145141510956,
         -0.014193750345572631, 0.13973481850126745},
    };
    for (auto &c : casos) {
        std::cout << "  " << c.nome << ": p0 = " << c.P.P0() << ", pc0 = " << c.P.Pc0() << "\n";
        const auto num = TriaxialDrained(c.P, 0.2, 400);
        const auto ana = TriaxialAnalytical(c.P, 3000, 0.3);
        std::vector<std::array<REAL, 2>> aq, aev;
        for (auto &a : ana) {
            aq.push_back({a[0], a[2]});
            aev.push_back({a[0], a[3]});
        }
        REAL err = 0.;
        for (auto &n : num) err = std::max(err, std::fabs(n[2] - Interp(aq, n[0])));
        Check("q(ea = 20%) (Python)", num.back()[2], c.q20, 1.e-8);
        Check("ev(ea = 20%) (Python)", num.back()[3], c.ev20, 1.e-8);
        Check("q(ea = 20%) analítico (diferença < 0.3%)", num.back()[2], Interp(aq, 0.2), 3.e-3);
        Check("max |q_num - q_anal| (Python)", err, c.errmax, 1.e-6);
        std::ofstream f("triaxial_" + c.nome + ".csv");
        f << std::setprecision(12) << "ea,p,q,ev,eq,sa,pc\n";
        for (auto &n : num) f << n[0] << "," << n[1] << "," << n[2] << "," << n[3] << "," << n[4] << "," << n[5] << "," << n[6] << "\n";
        std::ofstream fa("triaxial_" + c.nome + "_analitico.csv");
        fa << std::setprecision(12) << "ea,p,q,ev,eq,sa\n";
        for (auto &a : ana) fa << a[0] << "," << a[1] << "," << a[2] << "," << a[3] << "," << a[4] << "," << a[5] << "\n";
    }
    // convergência com o passo (seção 4 de verificacao_rs2.py): erro ~ 1ª ordem
    std::cout << "  Convergência com o passo: max |q_num - q_anal| (kPa)\n     passos      NC (8.5)   OCR=5 (8.8)\n";
    std::vector<REAL> e85, e88;
    for (int n : {50, 100, 200, 400, 800}) {
        REAL errs[2];
        int ic = 0;
        for (int ic2 : {0, 3}) {
            const auto num = TriaxialDrained(casos[ic2].P, 0.2, n);
            const auto ana = TriaxialAnalytical(casos[ic2].P, 3000, 0.3);
            std::vector<std::array<REAL, 2>> aq;
            for (auto &a : ana) aq.push_back({a[0], a[2]});
            REAL err = 0.;
            for (auto &r : num) err = std::max(err, std::fabs(r[2] - Interp(aq, r[0])));
            errs[ic++] = err;
        }
        e85.push_back(errs[0]);
        e88.push_back(errs[1]);
        std::printf("     %6d  %12.3f  %12.3f\n", n, errs[0], errs[1]);
    }
    Check("ordem de convergência NC (50 -> 800 passos)", std::log(e85.front() / e85.back()) / std::log(16.), 1., 0.1);
    Check("ordem de convergência OCR=5 (50 -> 800 passos)", std::log(e88.front() / e88.back()) / std::log(16.), 1., 0.1);
}

// =====================================================================================
// 3) elementos finitos: amostra cúbica com TPZMatElastoPlastic<TPZModifiedCamClay>
// =====================================================================================

/// Limite não drenado (permeabilidade nula) do sistema u-p de Biot, com a pressão de poros
/// condensada no ponto de integração: p = -α M ε_v, σ_total = σ' + α² M ε_v I.
class TMatCamClayUndrained : public TMatCamClay {
public:
    TMatCamClayUndrained(int id, REAL biotModulus, REAL biotAlpha = 1.)
        : TMatCamClay(id), fBiotM(biotModulus), fBiotAlpha(biotAlpha) {}

    TPZMaterial *NewMaterial() const override { return new TMatCamClayUndrained(*this); }

    std::string Name() const override { return "TMatCamClayUndrained"; }

    void Contribute(const TPZMaterialDataT<STATE> &data, REAL weight, TPZFMatrix<STATE> &ek,
                    TPZFMatrix<STATE> &ef) override {
        const REAL ev = TotalVolStrain(data);  // antes de a memória ser atualizada
        TMatCamClay::Contribute(data, weight, ek, ef);
        AddBiot(data, weight, ev, &ek, ef);
    }

    void Contribute(const TPZMaterialDataT<STATE> &data, REAL weight, TPZFMatrix<STATE> &ef) override {
        const REAL ev = TotalVolStrain(data);
        TMatCamClay::Contribute(data, weight, ef);
        AddBiot(data, weight, ev, nullptr, ef);
    }

    REAL PorePressure(REAL ev) const { return -fBiotAlpha * fBiotM * ev; }

private:
    REAL TotalVolStrain(const TPZMaterialDataT<STATE> &data) {
        TPZFNMatrix<6, REAL> ds(6, 1, 0.);
        this->ComputeDeltaStrainVector(data, ds);
        const TPZTensor<REAL> &eps = this->MemItem(data.intGlobPtIndex).m_elastoplastic_state.m_eps_t;
        return ds(_XX_, 0) + ds(_YY_, 0) + ds(_ZZ_, 0) + eps[_XX_] + eps[_YY_] + eps[_ZZ_];
    }

    void AddBiot(const TPZMaterialDataT<STATE> &data, REAL weight, REAL ev, TPZFMatrix<STATE> *ek,
                 TPZFMatrix<STATE> &ef) {
        TPZFNMatrix<9, REAL> axesT;
        TPZFNMatrix<60, REAL> dphiXYZ;
        data.axes.Transpose(&axesT);
        axesT.Multiply(data.dphix, dphiXYZ);
        const int phr = data.phi.Rows();
        std::vector<REAL> Bm(3 * phr);  // Bᵀ m
        for (int a = 0; a < phr; a++)
            for (int d = 0; d < 3; d++) Bm[3 * a + d] = dphiXYZ(d, a);
        const REAL kw = fBiotAlpha * fBiotAlpha * fBiotM * weight;
        for (int i = 0; i < 3 * phr; i++) {
            ef(i, 0) -= kw * ev * Bm[i];
            if (ek)
                for (int j = 0; j < 3 * phr; j++) (*ek)(i, j) += kw * Bm[i] * Bm[j];
        }
    }

    REAL fBiotM;
    REAL fBiotAlpha;
};

/// Caixa [0,L]³ com n x n x n hexaedros. Contornos: -1 x=0, -2 x=L, -3 y=0, -4 y=L, -5 z=0, -6 z=L
static TPZGeoMesh *BoxMesh(int n, REAL L = 1.) {
    auto *gmesh = new TPZGeoMesh;
    gmesh->SetDimension(3);
    const int np = n + 1;
    auto id = [np](int i, int j, int k) { return int64_t(k) * np * np + int64_t(j) * np + i; };
    gmesh->NodeVec().Resize(int64_t(np) * np * np);
    for (int k = 0; k < np; k++)
        for (int j = 0; j < np; j++)
            for (int i = 0; i < np; i++) {
                TPZManVector<REAL, 3> x = {L * i / n, L * j / n, L * k / n};
                gmesh->NodeVec()[id(i, j, k)].Initialize(x, *gmesh);
            }
    int64_t index;
    for (int k = 0; k < n; k++)
        for (int j = 0; j < n; j++)
            for (int i = 0; i < n; i++) {
                TPZManVector<int64_t, 8> t = {id(i, j, k),         id(i + 1, j, k),         id(i + 1, j + 1, k),
                                              id(i, j + 1, k),     id(i, j, k + 1),         id(i + 1, j, k + 1),
                                              id(i + 1, j + 1, k + 1), id(i, j + 1, k + 1)};
                gmesh->CreateGeoElement(ECube, t, 1, index);
            }
    auto quad = [&](int64_t a, int64_t b, int64_t c, int64_t d, int mat) {
        TPZManVector<int64_t, 4> t = {a, b, c, d};
        gmesh->CreateGeoElement(EQuadrilateral, t, mat, index);
    };
    for (int a = 0; a < n; a++)
        for (int b = 0; b < n; b++) {
            quad(id(0, a, b), id(0, a + 1, b), id(0, a + 1, b + 1), id(0, a, b + 1), -1);
            quad(id(n, a, b), id(n, a + 1, b), id(n, a + 1, b + 1), id(n, a, b + 1), -2);
            quad(id(a, 0, b), id(a + 1, 0, b), id(a + 1, 0, b + 1), id(a, 0, b + 1), -3);
            quad(id(a, n, b), id(a + 1, n, b), id(a + 1, n, b + 1), id(a, n, b + 1), -4);
            quad(id(a, b, 0), id(a + 1, b, 0), id(a + 1, b + 1, 0), id(a, b + 1, 0), -5);
            quad(id(a, b, n), id(a + 1, b, n), id(a + 1, b + 1, n), id(a, b + 1, n), -6);
        }
    gmesh->BuildConnectivity();
    return gmesh;
}

/// Linhas: ε_a, p', q, v, u, p_c (médias sobre os pontos de integração)
using TRow6 = std::array<REAL, 6>;

/// Ensaio triaxial em elementos finitos: σ_xx = σ_yy = σ0_xx nas faces x = L e y = L (tração),
/// simetria em x = 0, y = 0, z = 0 e deslocamento vertical prescrito no topo.
/// biotModulus <= 0 -> drenado; > 0 -> não drenado (α = 1).
static std::vector<TRow6> TriaxialFE(const TMCC &model, REAL biotModulus, REAL ea_max, int nsteps, int n,
                                     int &maxNewton) {
    const REAL L = 1.;
    TPZGeoMesh *gmesh = BoxMesh(n, L);
    auto *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(3);
    cmesh->SetDefaultOrder(1);
    cmesh->SetAllCreateFunctionsContinuousWithMem();

    const bool drained = biotModulus <= 0.;
    TMatCamClay *mat = drained ? new TMatCamClay(1) : new TMatCamClayUndrained(1, biotModulus);
    TMCC plastic(model);
    mat->SetPlasticityModel(plastic);  // memória padrão: σ = σ0
    cmesh->InsertMaterialObject(mat);

    const REAL sc = model.InitialStress()[_XX_];
    TPZFNMatrix<9, STATE> v1(3, 3, 0.);
    TPZManVector<STATE, 3> v2(3, 0.);
    v2 = {1., 0., 0.};
    cmesh->InsertMaterialObject(mat->CreateBC(mat, -1, 3, v1, v2));
    v2 = {0., 1., 0.};
    cmesh->InsertMaterialObject(mat->CreateBC(mat, -3, 3, v1, v2));
    v2 = {0., 0., 1.};
    cmesh->InsertMaterialObject(mat->CreateBC(mat, -5, 3, v1, v2));
    v2 = {sc, 0., 0.};
    cmesh->InsertMaterialObject(mat->CreateBC(mat, -2, 1, v1, v2));
    v2 = {0., sc, 0.};
    cmesh->InsertMaterialObject(mat->CreateBC(mat, -4, 1, v1, v2));
    TPZFNMatrix<9, STATE> mask(3, 3, 0.);
    mask(2, 2) = 1.;
    v2 = {0., 0., 0.};
    auto *top = mat->CreateBC(mat, -6, 6, mask, v2);  // u_z prescrito (incremento do passo)
    cmesh->InsertMaterialObject(top);
    cmesh->AutoBuild();

    TPZElastoPlasticAnalysis an(cmesh, std::cout);
    TPZFStructMatrix<STATE> str(cmesh);
    an.SetStructuralMatrix(str);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELU);
    an.SetSolver(step);
    const int64_t neq = cmesh->NEquations();

    auto record = [&](REAL ea, std::vector<TRow6> &hist) {
        TPZTensor<REAL> sig, eps;
        REAL pc = 0.;
        int npts = 0;
        for (int64_t el = 0; el < cmesh->NElements(); el++) {
            TPZCompEl *cel = cmesh->Element(el);
            if (!cel || cel->Dimension() != 3) continue;
            TPZManVector<int64_t> idx;
            cel->GetMemoryIndices(idx);
            for (int64_t ip : idx) {
                const TPZElastoPlasticMem &mem = mat->MemItem(ip);
                sig += mem.m_sigma;
                eps += mem.m_elastoplastic_state.m_eps_t;
                pc += model.Pc(mem.m_elastoplastic_state.m_hardening);
                npts++;
            }
        }
        sig *= 1. / npts;
        eps *= 1. / npts;
        REAL p, q;
        TMCC::Invariants(sig, p, q);
        const REAL ev = eps[_XX_] + eps[_YY_] + eps[_ZZ_];
        const REAL u = drained ? 0. : dynamic_cast<TMatCamClayUndrained *>(mat)->PorePressure(ev);
        hist.push_back({ea, -p, q, model.V0() * (1. + ev), u, pc / npts});
    };

    std::vector<TRow6> hist;
    record(0., hist);
    std::vector<bool> freeEq;
    maxNewton = 0;
    const REAL dea = ea_max / nsteps;
    for (int k = 1; k <= nsteps; k++) {
        TPZManVector<STATE, 3> val = {0., 0., -dea * L};
        top->SetVal2(val);
        TPZFMatrix<STATE> x(neq, 1, 0.);
        an.LoadSolution(x);
        bool conv = false;
        int it;
        for (it = 1; it <= 25; it++) {
            an.Assemble();
            if (freeEq.empty()) {  // equações sem penalidade de Dirichlet
                auto K = an.MatrixSolver<STATE>().Matrix();
                freeEq.resize(neq);
                for (int64_t i = 0; i < neq; i++) freeEq[i] = std::fabs(K->GetVal(i, i)) < 1.e10;
            }
            const TPZFMatrix<STATE> rhs = an.Rhs();
            REAL rfree = 0.;
            for (int64_t i = 0; i < neq; i++)
                if (freeEq[i]) rfree = std::max(rfree, std::fabs(rhs.GetVal(i, 0)));
            if (it > 1 && rfree <= 1.e-10 * std::fabs(sc) * L * L) {
                conv = true;
                break;
            }
            an.Solve();
            const TPZFMatrix<STATE> dx = an.Solution();
            x += dx;
            an.LoadSolution(x);
        }
        if (!conv) {
            std::cout << "TriaxialFE: Newton global não convergiu no passo " << k << "\n";
            break;
        }
        maxNewton = std::max(maxNewton, it - 1);
        an.AcceptSolution();
        record(dea * k * 1., hist);
    }
    delete cmesh;
    delete gmesh;
    return hist;
}

static void GravaCSV(const std::string &nome, const std::vector<TRow6> &h) {
    std::ofstream f(nome);
    f << std::setprecision(12) << "ea,p,q,v,u,pc\n";
    for (auto &r : h) f << r[0] << "," << r[1] << "," << r[2] << "," << r[3] << "," << r[4] << "," << r[5] << "\n";
}

static void TesteItasca() {
    std::cout << "\n" << std::string(100, '=') << "\n3) Benchmark Itasca: ensaio triaxial drenado e não drenado em"
              << " amostra de Cam-Clay (1 hexaedro, TPZMatElastoPlastic)\n";
    // parâmetros (Itasca): M = 1.02, λ = 0.2, κ = 0.05, v_λ = 3.32 (p1 = 1 kPa), G = 250 kPa, p'0 = 5 kPa
    auto params = [](REAL R) {
        TMCC m;
        m.SetUp(1.02, 0.2, 0.05, 3.32, -1., R * 5., 5., 0., 1., TMCC::EPressureDependent, TMCC::EConstantG, 250.);
        return m;
    };
    struct Ref {
        REAL p, q, v, u;
    };
    // FE Python (poro_camclay_fem.py, 1 hexaedro u-p) e FLAC3D (tabelas da Itasca)
    const Ref py_d16 = {7.573292329731367, 7.719876989194128, 2.8111971280170205, 0.}, flac_d16 = {7.573, 7.718, 2.811, 0.};
    const Ref py_d8 = {7.583881848464304, 7.751645544348984, 2.8105109617028727, 0.}, flac_d8 = {7.583, 7.747, 2.811, 0.};
    const Ref py_u16 = {4.23437830785905, 4.318526792569606, 2.9273993413277837, 2.2051306230111782},
              flac_u16 = {4.234, 4.312, 2.927, 2.203};
    const Ref py_u8 = {14.047913734946093, 14.422435155179878, 2.6865536965569357, -4.240435349622619},
              flac_u8 = {14.05, 14.42, 2.687, -4.241};

    struct Caso {
        std::string nome;
        REAL R;
        bool drained;
        Ref py, flac;
    };
    const std::vector<Caso> casos = {{"drenado_R1.6", 1.6, true, py_d16, flac_d16},
                                     {"drenado_R8", 8.0, true, py_d8, flac_d8},
                                     {"nao_drenado_R1.6", 1.6, false, py_u16, flac_u16},
                                     {"nao_drenado_R8", 8.0, false, py_u8, flac_u8}};
    for (auto &c : casos) {
        TMCC P = params(c.R);
        const REAL n0 = (P.V0() - 1.) / P.V0();
        int maxit;
        const auto h = c.drained ? TriaxialFE(P, 0., 0.5, 500, 1, maxit) : TriaxialFE(P, 2.e4 / n0, 0.1, 400, 1, maxit);
        GravaCSV("itasca_" + c.nome + ".csv", h);
        const auto &f = h.back();
        std::printf("  %s (v0 = %.4f): %zu passos, máx. %d iterações de Newton por passo\n", c.nome.c_str(), P.V0(),
                    h.size() - 1, maxit);
        std::printf("      fim do ensaio: p' = %.3f (FLAC %.3f)  q = %.3f (FLAC %.3f)  v = %.3f (FLAC %.3f)", f[1],
                    c.flac.p, f[2], c.flac.q, f[3], c.flac.v);
        if (!c.drained) std::printf("  u = %.3f (FLAC %.3f)", f[4], c.flac.u);
        std::printf("\n");
        const REAL tol = c.drained ? 1.e-6 : 1.e-4;
        Check("p' final (FE Python)", f[1], c.py.p, tol);
        Check("q final (FE Python)", f[2], c.py.q, tol);
        Check("v final (FE Python)", f[3], c.py.v, 1.e-6);
        if (!c.drained) Check("u final (FE Python)", f[4], c.py.u, tol);
        Check("q final (FLAC3D, 0.2%)", f[2], c.flac.q, 2.e-3);
        if (c.drained) {
            REAL err = 0.;
            for (auto &r : h) err = std::max(err, std::fabs(r[2] - 3. * (r[1] - 5.)));
            CheckBelow("trajetória drenada max |q - 3(p' - p'0)|", err, 1.e-7);
        }
    }
    // malha 2x2x2 (deformação homogênea: mesmo resultado de 1 elemento)
    {
        TMCC P = params(8.);
        int maxit;
        const auto h1 = TriaxialFE(P, 0., 0.1, 100, 1, maxit);
        const auto h2 = TriaxialFE(P, 0., 0.1, 100, 2, maxit);
        Check("drenado R=8, ea = 10%: q malha 2x2x2 x 1 elemento", h2.back()[2], h1.back()[2], 1.e-8);
    }

    // ---- hipoelástico (benchmark 1.15.2 do Abaqus) em FE x ponto material
    std::cout << "  Hipoelástico (Abaqus 1.15.2: M = 1, lambda = 0.174, kappa = 0.026, v0 = 2.08, a0 = 58.3, p0 = 100)\n";
    TMCC PH;
    PH.SetUp(1.0, 0.174, 0.026, 0., 2.08, 2. * 58.3, 100., 0., 1., TMCC::EPressureDependent, TMCC::EHypoNu, 0., 0.3);
    int maxit;
    const auto hfe = TriaxialFE(PH, 0., 0.1, 200, 1, maxit);
    GravaCSV("hypo_fe.csv", hfe);
    const auto hpt = TriaxialDrained(PH, 0.1, 200);
    Check("hypo FE: p' em ea = 10% (ponto material Python)", hfe.back()[1], 134.68786122479625, 1.e-6);
    Check("hypo FE: q em ea = 10% (ponto material Python)", hfe.back()[2], 104.06358367438861, 1.e-6);
    Check("hypo FE: ev em ea = 10% (ponto material Python)", 1. - hfe.back()[3] / PH.V0(), 0.04729048128779287, 1.e-6);
    Check("hypo ponto material: q em ea = 10% (Python)", hpt.back()[2], 104.06358367438861, 1.e-8);
}

int main() {
    std::cout << std::setprecision(10);
    TestePontoMaterial();
    TesteTriaxialDrenado();
    TesteItasca();
    std::cout << "\n" << std::string(100, '=') << "\n";
    if (gFalhas) {
        std::cout << gFalhas << " verificação(ões) FALHARAM\n";
        return 1;
    }
    std::cout << "Todas as verificações passaram\n";
    return 0;
}
