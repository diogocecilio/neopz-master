// Slope stability (plane strain, associative perfectly plastic Mohr-Coulomb):
// factor of safety by gravity increase and by shear strength reduction, with h-refinement
// of the plastic zone.
//
// Constitutive model (default): closed-form closest-point projection in rotated
// Haigh-Westergaard space with consistent tangent (Lira Cecilio, J Eng Math 157:10, 2026),
// i.e. TPZPlasticStepVoigt<TPZYCMohrCoulombPV2>. "pv" selects the legacy iterative
// TPZPlasticStepPV<TPZYCMohrCoulombPV> for cross-checking.
//
// Usage: SlopeMohrCoulomb [pv] [check] [nref=<n>] [nu=<poisson>]
#include "SlopeAnalysis.h"
#include "Plasticity/TPZPlasticStepPV.h"
#include "Plasticity/TPZPlasticStepVoigt.h"
#include "Plasticity/TPZYCMohrCoulombPV.h"
#include "Plasticity/TPZYCMohrCoulombPV2.h"
#include "TPZGeoLinear.h"
#include "TPZVTKGeoMesh.h"
#include "pzgeotriangle.h"
#include "tpzgeoelrefpattern.h"

#include <cstring>
#include <random>
#include <string>
#include <vector>

using TMCVoigt = TPZPlasticStepVoigt<TPZYCMohrCoulombPV2, TPZElasticResponse>;
using TMCPV = TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>;

struct Soil {
    REAL E = 20000., nu = 0.49, c = 10., phi = 30. * M_PI / 180., gamma = 20.;
};

TMCVoigt ModelVoigt(const Soil &s) {
    TPZElasticResponse ER;
    ER.SetEngineeringData(s.E, s.nu);
    return TMCVoigt(TPZYCMohrCoulombPV2(s.phi, s.phi, s.c, ER), ER);
}

TMCPV ModelPV(const Soil &s) {
    TPZElasticResponse ER;
    ER.SetEngineeringData(s.E, s.nu);
    TMCPV pv;
    pv.fYC.SetUp(s.phi, s.phi, s.c, ER);
    pv.fER = ER;
    return pv;
}

/// Slope 70 x 40 m, height 10 m at 45 degrees. BC ids: -1 base, -2 right, -3 top right,
/// -4 top left, -5 left, -6 face. ref = uniform refinements.
TPZGeoMesh *TriGMesh(int ref) {
    auto *gmesh = new TPZGeoMesh();
    gmesh->SetDimension(2);
    const std::vector<std::vector<REAL>> co = {
        {0, 0}, {10, 0}, {20, 0}, {30, 0}, {40, 0}, {50, 0}, {60, 0}, {70, 0},
        {0, 10}, {10, 10}, {20, 10}, {30, 10}, {40, 10}, {50, 10}, {60, 10}, {70, 10},
        {0, 20}, {10, 20}, {20, 20}, {30, 20}, {40, 20}, {50, 20}, {60, 20}, {70, 20},
        {0, 30}, {10, 30}, {20, 30}, {30, 30}, {40, 30}, {50, 30}, {60, 30}, {70, 30},
        {0, 40}, {10, 40}, {20, 40}, {30, 40}};
    const std::vector<std::vector<int64_t>> tri = {
        {0, 1, 8}, {1, 9, 8}, {1, 2, 9}, {2, 10, 9}, {2, 3, 10}, {3, 11, 10}, {3, 4, 11},
        {4, 12, 11}, {4, 5, 12}, {5, 13, 12}, {5, 6, 13}, {6, 14, 13}, {6, 7, 14}, {7, 15, 14},
        {8, 9, 16}, {9, 17, 16}, {9, 10, 17}, {10, 18, 17}, {10, 11, 18}, {11, 19, 18}, {11, 12, 19},
        {12, 20, 19}, {12, 13, 20}, {13, 21, 20}, {13, 14, 21}, {14, 22, 21}, {14, 15, 22}, {15, 23, 22},
        {16, 17, 24}, {17, 25, 24}, {17, 18, 25}, {18, 26, 25}, {18, 19, 26}, {19, 27, 26}, {19, 20, 27},
        {20, 28, 27}, {20, 21, 28}, {21, 29, 28}, {21, 22, 29}, {22, 30, 29}, {22, 23, 30}, {23, 31, 30},
        {24, 25, 32}, {25, 33, 32}, {25, 26, 33}, {26, 34, 33}, {26, 27, 34}, {27, 35, 34}, {27, 28, 35}};
    const std::vector<std::pair<std::vector<int64_t>, int>> lines = {
        {{0, 1}, -1}, {{1, 2}, -1}, {{2, 3}, -1}, {{3, 4}, -1}, {{4, 5}, -1}, {{5, 6}, -1}, {{6, 7}, -1},
        {{7, 15}, -2}, {{15, 23}, -2}, {{23, 31}, -2},
        {{31, 30}, -3}, {{30, 29}, -3}, {{29, 28}, -3},
        {{35, 34}, -4}, {{34, 33}, -4}, {{33, 32}, -4},
        {{32, 24}, -5}, {{24, 16}, -5}, {{16, 8}, -5}, {{8, 0}, -5},
        {{28, 35}, -6}};
    gmesh->NodeVec().Resize(co.size());
    for (size_t i = 0; i < co.size(); i++) {
        TPZManVector<REAL, 3> x = {co[i][0], co[i][1], 0.};
        gmesh->NodeVec()[i] = TPZGeoNode(i, x, *gmesh);
    }
    for (auto &t : tri) {
        TPZManVector<int64_t, 3> nodes = {t[0], t[1], t[2]};
        new TPZGeoElRefPattern<pzgeom::TPZGeoTriangle>(nodes, 1, *gmesh);
    }
    for (auto &l : lines) {
        TPZManVector<int64_t, 2> nodes = {l.first[0], l.first[1]};
        new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodes, l.second, *gmesh);
    }
    gmesh->BuildConnectivity();
    for (int d = 0; d < ref; d++) {
        const int64_t nel = gmesh->NElements();
        TPZManVector<TPZGeoEl *> sub;
        for (int64_t iel = 0; iel < nel; iel++)
            if (!gmesh->Element(iel)->HasSubElement()) gmesh->Element(iel)->Divide(sub);
    }
    return gmesh;
}

/// H1 mesh with memory; base fixed, lateral sides on rollers (BC type 3: val2 = constrained directions)
template <class TPlastic>
TPZCompMesh *CreateCMesh(TPZGeoMesh *gmesh, int porder, const TPlastic &model, const Soil &s) {
    auto *mat = new TPZMatElastoPlastic2D<TPlastic, TPZElastoPlasticMem>(1, 1);
    TPlastic m(model);
    mat->SetPlasticityModel(m);
    mat->SetBodyForce({0., -s.gamma, 0.});
    mat->SetBigNumber(1.e12);
    auto *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(porder);
    cmesh->SetDimModel(2);
    cmesh->InsertMaterialObject(mat);
    TPZFMatrix<STATE> val1(2, 2, 0.);
    const int bcid[3] = {-1, -2, -5};
    const REAL fixy[3] = {1., 0., 0.};
    for (int i = 0; i < 3; i++) {
        TPZManVector<STATE, 2> val2 = {1., fixy[i]};
        cmesh->InsertMaterialObject(mat->CreateBC(mat, bcid[i], 3, val1, val2));
    }
    cmesh->SetAllCreateFunctionsContinuousWithMem();
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    return cmesh;
}

/// Taylor test of the material tangent (paper Eq. 63-65): order ~2 expected in every regime.
/// Strains are built from principal values in a random frame; Voigt strain = engineering shear.
template <class TPlastic>
void CheckTangent(const TPlastic &model, const char *name) {
    std::mt19937 gen(7);
    std::uniform_real_distribution<REAL> U(-1., 1.);
    const REAL cases[][3] = {{4.e-3, 1.e-3, -5.e-3},  {5.e-3, -2.5e-3, -2.5e-3}, {2.5e-3, 2.5e-3, -5.e-3},
                             {1.e-3, 1.e-3, 1.e-3},   {3.e-3, -1.e-3, -2.e-3},   {2.e-4, 1.e-4, -3.e-4}};
    const char *label[] = {"general", "e2 = e3", "e1 = e2", "tension", "general", "small"};
    std::cout << "\nTaylor test of the consistent tangent: " << name << "\n";
    for (int test = 0; test < 6; test++) {
        TPZFNMatrix<9, REAL> Q(3, 3), A(3, 3);                 // random rotation (Gram-Schmidt)
        for (int i = 0; i < 9; i++) A(i / 3, i % 3) = U(gen);
        for (int i = 0; i < 3; i++) {
            for (int k = 0; k < i; k++) {
                REAL d = 0.; for (int j = 0; j < 3; j++) d += A(j, i) * Q(j, k);
                for (int j = 0; j < 3; j++) A(j, i) -= d * Q(j, k);
            }
            REAL n = 0.; for (int j = 0; j < 3; j++) n += A(j, i) * A(j, i);
            for (int j = 0; j < 3; j++) Q(j, i) = A(j, i) / std::sqrt(n);
        }
        const REAL shift = test == 3 ? 0. : -2.e-5;            // mild volumetric compression
        TPZTensor<REAL> eps, deps;
        const int ij[6][2] = {{0, 0}, {0, 1}, {0, 2}, {1, 1}, {1, 2}, {2, 2}};
        for (int v = 0; v < 6; v++) {
            REAL e = 0.;
            for (int k = 0; k < 3; k++) e += Q(ij[v][0], k) * (cases[test][k] + shift) * Q(ij[v][1], k);
            eps[v] = (ij[v][0] == ij[v][1]) ? e : 2. * e;
            deps[v] = 1.e-3 * U(gen);
        }
        auto stress = [&](const TPZTensor<REAL> &e, TPZFMatrix<REAL> *D) {
            TPlastic m(model);
            m.SetState(TPZPlasticState<REAL>());
            TPZTensor<REAL> sig;
            m.ApplyStrainComputeSigma(e, sig, D);
            return std::make_pair(sig, m.GetState().m_m_type);
        };
        TPZFNMatrix<36, REAL> D(6, 6, 0.);
        auto [s0, type] = stress(eps, &D);
        REAL err[2], alpha[2] = {1.e-3, 2.e-3};
        for (int k = 0; k < 2; k++) {
            TPZTensor<REAL> e(eps);
            e.Add(deps, alpha[k]);
            TPZTensor<REAL> s1 = stress(e, nullptr).first;
            REAL n2 = 0.;
            for (int i = 0; i < 6; i++) {
                REAL lin = 0.;
                for (int j = 0; j < 6; j++) lin += D(i, j) * deps[j] * alpha[k];
                n2 += std::pow(s1[i] - s0[i] - lin, 2);
            }
            err[k] = std::sqrt(n2);
        }
        REAL asym = 0., dmax = 0.; // associative plasticity: symmetric tangent (skyline LDLt)
        for (int i = 0; i < 6; i++)
            for (int j = 0; j < 6; j++) { asym = std::max(asym, std::fabs(D(i, j) - D(j, i))); dmax = std::max(dmax, std::fabs(D(i, j))); }
        std::cout << "  " << label[test] << ": m_type = " << type << "  |D-D^T|/|D| = " << (dmax > 0. ? asym / dmax : 0.) << "  |E| = " << err[0];
        if (err[0] > 1.e-10 * Norm(s0)) std::cout << "  order p = " << std::log(err[1] / err[0]) / std::log(alpha[1] / alpha[0]);
        std::cout << "\n";
    }
}

template <class TPlastic>
void Run(const TPlastic &model, const Soil &s, int nref, const std::string &tag) {
    TPZGeoMesh *gmesh = TriGMesh(1);
    TPZCompMesh *cmesh = CreateCMesh(gmesh, 2, model, s);
    SlopeAnalysis<TPlastic> slope(cmesh);
    std::vector<std::vector<REAL>> table;
    for (int k = 0;; k++) {
        const int64_t neq = cmesh->NEquations();
        const REAL fsGI = slope.GravityIncrease();
        slope.PostProcess(tag + "_GI_ref" + std::to_string(k) + ".vtk");
        slope.MarkPlasticZone(0.1); // failure mechanism at collapse
        const REAL fsSRM = slope.StrengthReduction();
        slope.PostProcess(tag + "_SRM_ref" + std::to_string(k) + ".vtk");
        slope.MarkPlasticZone(0.1);
        table.push_back({REAL(k), REAL(neq), fsGI, fsSRM});
        if (k == nref) break;
        slope.Refine();
    }
    std::ofstream vtk(tag + "_mesh.vtk");
    TPZVTKGeoMesh::PrintCMeshVTK(cmesh, vtk, true);
    std::cout << "\n" << tag << ": refinement  equations  FS(gravity increase)  FS(strength reduction)\n";
    for (auto &r : table) std::cout << "  " << r[0] << "  " << r[1] << "  " << r[2] << "  " << r[3] << "\n";
    delete cmesh;
    delete gmesh;
}

int main(int argc, char *argv[]) {
    bool pv = false, check = false;
    int nref = 5;
    Soil s;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "pv")) pv = true;
        else if (!strcmp(argv[i], "check")) check = true;
        else if (!strncmp(argv[i], "nref=", 5)) nref = atoi(argv[i] + 5);
        else if (!strncmp(argv[i], "nu=", 3)) s.nu = atof(argv[i] + 3);
    }
    if (check) {
        CheckTangent(ModelVoigt(s), "TPZPlasticStepVoigt<TPZYCMohrCoulombPV2>");
        CheckTangent(ModelPV(s), "TPZPlasticStepPV<TPZYCMohrCoulombPV>");
        return 0;
    }
    if (pv) Run(ModelPV(s), s, nref, "slope_pv");
    else Run(ModelVoigt(s), s, nref, "slope_rhw");
    return 0;
}
