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
#include "SlopeModel.h"

#include <cstring>
#include <random>
#include <string>
#include <vector>

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
        slope.PostPlasticity(tag + "_GI_ref" + std::to_string(k) + ".vtk");
        slope.MarkPlasticZone(0.1); // failure mechanism at collapse
        const REAL fsSRM = slope.StrengthReduction();
        slope.PostPlasticity(tag + "_SRM_ref" + std::to_string(k) + ".vtk");
        slope.MarkPlasticZone(0.1);
        table.push_back({REAL(k), REAL(neq), fsGI, fsSRM});
        if (k == nref) break;
        slope.Refine();
    }
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
