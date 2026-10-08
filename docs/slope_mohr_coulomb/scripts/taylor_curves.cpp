// Taylor test data for a log-log figure of the consistent tangent (Lira Cecilio, J Eng Math 157:10, 2026,
// Eq. 63-65). Construction of eps, deps copied verbatim from CheckTangent in
// Projects/SlopeMohrCoulomb/main.cpp (std::mt19937 gen(7), same cases, labels, shift, deps).
// For alpha = 10^-1, 10^-1.5, ..., 10^-7 it records
//   err2 = |sigma(eps + alpha deps) - sigma(eps) - alpha D deps|   (second-order remainder)
//   err1 = |sigma(eps + alpha deps) - sigma(eps)|                  (first-order remainder)
// Usage: taylor <taylor.csv> <taylor_full.csv> <case_info.csv>
#include "SlopeModel.h" // compile with -I<neopz>/Projects/SlopeMohrCoulomb and the flags of SlopeMohrCoulomb

#include <cstdio>
#include <random>
#include <string>
#include <vector>

static const REAL cases[][3] = {{4.e-3, 1.e-3, -5.e-3},  {5.e-3, -2.5e-3, -2.5e-3}, {2.5e-3, 2.5e-3, -5.e-3},
                                {1.e-3, 1.e-3, 1.e-3},   {3.e-3, -1.e-3, -2.e-3},   {2.e-4, 1.e-4, -3.e-4}};
static const char *label[] = {"general", "e2 = e3", "e1 = e2", "tension", "general", "small"};

template <class TPlastic>
void Taylor(const TPlastic &model, const char *mname, FILE *fcsv, FILE *ffull, FILE *finfo) {
    std::mt19937 gen(7);
    std::uniform_real_distribution<REAL> U(-1., 1.);
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
        auto remainders = [&](REAL alpha, REAL &e2, REAL &e1, int &tp) {
            TPZTensor<REAL> e(eps);
            e.Add(deps, alpha);
            auto r = stress(e, nullptr);
            const TPZTensor<REAL> &s1 = r.first;
            tp = r.second;
            REAL n2 = 0., n1 = 0.;
            for (int i = 0; i < 6; i++) {
                REAL lin = 0.;
                for (int j = 0; j < 6; j++) lin += D(i, j) * deps[j] * alpha;
                n2 += std::pow(s1[i] - s0[i] - lin, 2);
                n1 += std::pow(s1[i] - s0[i], 2);
            }
            e2 = std::sqrt(n2);
            e1 = std::sqrt(n1);
        };
        // reproduction of the original two-point check (alpha = 1e-3, 2e-3)
        REAL eo[2], dummy;
        int tdummy;
        remainders(1.e-3, eo[0], dummy, tdummy);
        remainders(2.e-3, eo[1], dummy, tdummy);
        const REAL p = std::log(eo[1] / eo[0]) / std::log(2.);
        REAL asym = 0., dmax = 0.;
        for (int i = 0; i < 6; i++)
            for (int j = 0; j < 6; j++) { asym = std::max(asym, std::fabs(D(i, j) - D(j, i))); dmax = std::max(dmax, std::fabs(D(i, j))); }
        const REAL s0n = Norm(s0);
        std::fprintf(finfo, "%s,%d,%s,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n", mname, test + 1, label[test], type,
                     cases[test][0] + shift, cases[test][1] + shift, cases[test][2] + shift, s0n, dmax > 0. ? asym / dmax : 0.,
                     eo[0], eo[1], p);
        std::printf("%-6s case %d (%-8s) m_type=%d |sig0|=%.6g |D-D^T|/|D|=%.3g |E(1e-3)|=%.6g p(1e-3,2e-3)=%.6f\n", mname, test + 1,
                    label[test], type, s0n, dmax > 0. ? asym / dmax : 0., eo[0], p);
        for (int k = 2; k <= 14; k++) {                         // alpha = 10^(-k/2), k = 2..14
            const REAL alpha = std::pow(10., -0.5 * k);
            REAL e2, e1;
            int tp;
            remainders(alpha, e2, e1, tp);
            std::fprintf(fcsv, "%s,%d,%d,%.17g,%.17g,%.17g\n", mname, test + 1, type, alpha, e2, e1);
            std::fprintf(ffull, "%s,%d,%s,%d,%d,%.17g,%.17g,%.17g,%.17g\n", mname, test + 1, label[test], type, tp, alpha, e2, e1, s0n);
        }
    }
}

int main(int argc, char *argv[]) {
    if (argc < 4) { std::fprintf(stderr, "usage: taylor taylor.csv taylor_full.csv case_info.csv\n"); return 1; }
    FILE *fcsv = std::fopen(argv[1], "w"), *ffull = std::fopen(argv[2], "w"), *finfo = std::fopen(argv[3], "w");
    std::fprintf(fcsv, "model,case,m_type,alpha,err2,err1\n");
    std::fprintf(ffull, "model,case,label,m_type,m_type_pert,alpha,err2,err1,sig0_norm\n");
    std::fprintf(finfo, "model,case,label,m_type,e1,e2,e3,sig0_norm,asym_rel,err_1e-3,err_2e-3,p_orig\n");
    Soil s;
    Taylor(ModelVoigt(s), "Voigt", fcsv, ffull, finfo);
    Taylor(ModelPV(s), "PV", fcsv, ffull, finfo);
    std::fclose(fcsv);
    std::fclose(ffull);
    std::fclose(finfo);
    return 0;
}
