// Diagnostic of the err2 floor of the Voigt model in the edge (and apex) cases: the library stress
// sigma = sigtr + sum_i (sproj_i - sigtr_i) n_i n_i uses TPZTensor::EigenSystem (closed form) for n_i; here the
// same return (TPZYCMohrCoulombPV2::ProjectSigma of the library) is rebuilt with a long-double cyclic
// Jacobi eigen-decomposition (orthonormal eigenvectors to round-off), with the same library tangent D.
// Same eps/deps construction as CheckTangent (mt19937 gen(7)).
// Usage: diag_floor <out.csv>
#include "SlopeModel.h" // compile with -I<neopz>/Projects/SlopeMohrCoulomb and the flags of SlopeMohrCoulomb

#include <algorithm>
#include <cstdio>
#include <random>

static const REAL cases[][3] = {{4.e-3, 1.e-3, -5.e-3},  {5.e-3, -2.5e-3, -2.5e-3}, {2.5e-3, 2.5e-3, -5.e-3},
                                {1.e-3, 1.e-3, 1.e-3},   {3.e-3, -1.e-3, -2.e-3},   {2.e-4, 1.e-4, -3.e-4}};
static const char *label[] = {"general", "e2 = e3", "e1 = e2", "tension", "general", "small"};

using LD = long double;
/// cyclic Jacobi in long double for a symmetric 3x3 (Voigt XX,XY,XZ,YY,YZ,ZZ): eigenvalues descending,
/// eigenvectors (columns of V) orthonormal to long-double round-off
static void JacobiLD(const TPZTensor<REAL> &t, LD lam[3], LD V[3][3]) {
    LD a[3][3] = {{t[0], t[1], t[2]}, {t[1], t[3], t[4]}, {t[2], t[4], t[5]}};
    for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) V[i][j] = (i == j);
    for (int sweep = 0; sweep < 50; sweep++) {
        LD off = std::fabs(a[0][1]) + std::fabs(a[0][2]) + std::fabs(a[1][2]);
        if (off == 0.L) break;
        for (int p = 0; p < 2; p++)
            for (int q = p + 1; q < 3; q++) {
                if (a[p][q] == 0.L) continue;
                LD th = (a[q][q] - a[p][p]) / (2.L * a[p][q]);
                LD tt = (th >= 0 ? 1.L : -1.L) / (std::fabs(th) + std::sqrt(th * th + 1.L));
                LD c = 1.L / std::sqrt(tt * tt + 1.L), s = tt * c;
                for (int k = 0; k < 3; k++) { LD akp = a[k][p], akq = a[k][q]; a[k][p] = c * akp - s * akq; a[k][q] = s * akp + c * akq; }
                for (int k = 0; k < 3; k++) { LD apk = a[p][k], aqk = a[q][k]; a[p][k] = c * apk - s * aqk; a[q][k] = s * apk + c * aqk; }
                for (int k = 0; k < 3; k++) { LD vkp = V[k][p], vkq = V[k][q]; V[k][p] = c * vkp - s * vkq; V[k][q] = s * vkp + c * vkq; }
            }
    }
    int idx[3] = {0, 1, 2};
    std::sort(idx, idx + 3, [&](int x, int y) { return a[x][x] > a[y][y]; });
    LD W[3][3];
    for (int i = 0; i < 3; i++) { lam[i] = a[idx[i]][idx[i]]; for (int k = 0; k < 3; k++) W[k][i] = V[k][idx[i]]; }
    for (int i = 0; i < 3; i++) for (int k = 0; k < 3; k++) V[k][i] = W[k][i];
}

int main(int argc, char *argv[]) {
    FILE *f = std::fopen(argc > 1 ? argv[1] : "diag_floor.csv", "w");
    std::fprintf(f, "case,label,alpha,m_type_pert,err2_lib,err2_jacobi,diff_lib_jacobi,orth_defect_lib,plastic_corr_norm,trial_gap_min,eigval_err_lib,eigproj_err_lib\n");
    Soil s;
    TMCVoigt model = ModelVoigt(s);
    TPZElasticResponse ER;
    ER.SetEngineeringData(s.E, s.nu);
    std::mt19937 gen(7);
    std::uniform_real_distribution<REAL> U(-1., 1.);
    for (int test = 0; test < 6; test++) {
        TPZFNMatrix<9, REAL> Q(3, 3), A(3, 3);
        for (int i = 0; i < 9; i++) A(i / 3, i % 3) = U(gen);
        for (int i = 0; i < 3; i++) {
            for (int k = 0; k < i; k++) {
                REAL d = 0.; for (int j = 0; j < 3; j++) d += A(j, i) * Q(j, k);
                for (int j = 0; j < 3; j++) A(j, i) -= d * Q(j, k);
            }
            REAL n = 0.; for (int j = 0; j < 3; j++) n += A(j, i) * A(j, i);
            for (int j = 0; j < 3; j++) Q(j, i) = A(j, i) / std::sqrt(n);
        }
        const REAL shift = test == 3 ? 0. : -2.e-5;
        TPZTensor<REAL> eps, deps;
        const int ij[6][2] = {{0, 0}, {0, 1}, {0, 2}, {1, 1}, {1, 2}, {2, 2}};
        for (int v = 0; v < 6; v++) {
            REAL e = 0.;
            for (int k = 0; k < 3; k++) e += Q(ij[v][0], k) * (cases[test][k] + shift) * Q(ij[v][1], k);
            eps[v] = (ij[v][0] == ij[v][1]) ? e : 2. * e;
            deps[v] = 1.e-3 * U(gen);
        }
        auto lib = [&](const TPZTensor<REAL> &e, TPZFMatrix<REAL> *D, int &type) {
            TMCVoigt m(model);
            m.SetState(TPZPlasticState<REAL>());
            TPZTensor<REAL> sig;
            m.ApplyStrainComputeSigma(e, sig, D);
            type = m.GetState().m_m_type;
            return sig;
        };
        // same return, eigen-decomposition by long-double Jacobi; also library eigenvector diagnostics
        auto ref = [&](const TPZTensor<REAL> &e, REAL &defect, REAL &corr, REAL &gap, REAL &everr, REAL &perr) {
            TPZTensor<REAL> sigtr;
            ER.ComputeStress(e, sigtr);
            LD lam[3], V[3][3];
            JacobiLD(sigtr, lam, V);
            TPZManVector<STATE, 3> str(3), spr(3), epstr(3);
            for (int i = 0; i < 3; i++) str[i] = (REAL)lam[i];
            TPZManVector<STATE, 2> dlambda(2, 0.);
            TPZFNMatrix<9> G3(3, 3, 0.);
            STATE h0 = 0., h1 = 0.;
            int type = 0;
            auto yc = model.LocalCriterion();
            yc.ProjectSigma(str, h0, dlambda, spr, epstr, G3, h1, type);
            TPZTensor<REAL> sig(sigtr);
            const int ij[6][2] = {{0, 0}, {0, 1}, {0, 2}, {1, 1}, {1, 2}, {2, 2}};
            if (type != 0)
                for (int v = 0; v < 6; v++) {
                    LD add = 0.L;
                    for (int i = 0; i < 3; i++) add += ((LD)spr[i] - lam[i]) * V[ij[v][0]][i] * V[ij[v][1]][i];
                    sig[v] = (REAL)((LD)sigtr[v] + add);
                }
            // library eigen-system of the same trial stress
            TPZTensor<REAL>::TPZDecomposed ed;
            sigtr.EigenSystem(ed);
            defect = 0.;
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++) {
                    REAL d = 0.; for (int k = 0; k < 3; k++) d += ed.fEigenvectors[i][k] * ed.fEigenvectors[j][k];
                    defect = std::max(defect, std::fabs(d - (i == j ? 1. : 0.)));
                }
            corr = 0.; for (int i = 0; i < 3; i++) corr = std::max(corr, std::fabs(spr[i] - str[i]));
            if (type == 0) corr = 0.;
            gap = std::min(std::fabs(str[0] - str[1]), std::fabs(str[1] - str[2]));
            everr = 0.; perr = 0.;
            for (int i = 0; i < 3; i++) {
                everr = std::max(everr, std::fabs(ed.fEigenvalues[i] - str[i]));
                for (int a = 0; a < 3; a++)
                    for (int b = 0; b < 3; b++)
                        perr = std::max(perr, (REAL)std::fabs(ed.fEigenvectors[i][a] * ed.fEigenvectors[i][b] - V[a][i] * V[b][i]));
            }
            return sig;
        };
        TPZFNMatrix<36, REAL> D(6, 6, 0.);
        int t0;
        TPZTensor<REAL> s0 = lib(eps, &D, t0);
        REAL dd, cc, gg, ee, pp;
        TPZTensor<REAL> r0 = ref(eps, dd, cc, gg, ee, pp);
        for (int k = 2; k <= 14; k++) {
            const REAL alpha = std::pow(10., -0.5 * k);
            TPZTensor<REAL> e(eps);
            e.Add(deps, alpha);
            int tp;
            TPZTensor<REAL> s1 = lib(e, nullptr, tp);
            REAL defect, corr, gap, everr, perr;
            TPZTensor<REAL> r1 = ref(e, defect, corr, gap, everr, perr);
            REAL n2l = 0., n2r = 0., nd = 0.;
            for (int i = 0; i < 6; i++) {
                REAL lin = 0.;
                for (int j = 0; j < 6; j++) lin += D(i, j) * deps[j] * alpha;
                n2l += std::pow(s1[i] - s0[i] - lin, 2);
                n2r += std::pow(r1[i] - r0[i] - lin, 2);
                nd += std::pow(s1[i] - r1[i], 2);
            }
            std::fprintf(f, "%d,%s,%.17g,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n", test + 1, label[test], alpha, tp, std::sqrt(n2l),
                         std::sqrt(n2r), std::sqrt(nd), defect, corr, gap, everr, perr);
        }
    }
    std::fclose(f);
    return 0;
}
