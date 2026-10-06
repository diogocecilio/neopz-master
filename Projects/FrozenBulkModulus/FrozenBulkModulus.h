/**
 * @file FrozenBulkModulus.h
 * @brief Sect. 6.1 and Table 3 of the article: exact integration of the porous elastic law during the
 * plastic correction (this work) versus the bulk modulus frozen at its trial value, in drained and
 * undrained triaxial tests at a material point of the RS2 clay.
 */
#pragma once

#include "MCCPaperTools.h"

#include <array>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

/**
 * @brief Exact and frozen integration of the porous law in the local return mapping (Sect. 4.3 and 6.1,
 * Table 3).
 *
 * During the plastic correction the volumetric equation of the local system (18) is
 * \f[
 *   R_1 = \frac{\xi-\xi_{tr}\exp(-v_0\Delta\alpha/\kappa)}{a_n}\quad\text{(exact, this work)},\qquad
 *   R_1 = \frac{\xi-\xi_{tr}(1-v_0\Delta\alpha/\kappa)}{a_n}\quad\text{(bulk modulus frozen at } K_{tr}\text{)}.
 * \f]
 * The frozen form reverses the sign of \f$\xi_{tr}\f$ when \f$v_0\Delta\alpha/\kappa>1\f$; for \f$p_t=0\f$
 * this bounds the growth of the preconsolidation pressure in one increment to
 * \f$p'_c/p'_{c,n}<\exp[\kappa/(\lambda-\kappa)]=1.098\f$. The two forms are selected with
 * TPZPlasticStepModifiedCamClay::SetPorousIntegration (TPZYCModifiedCamClayRHW::EExact or EFrozen); the
 * rest of the algorithm is unchanged.
 *
 * Material: RS2 clay, \f$M=1.2\f$, \f$\lambda=0.077\f$, \f$\kappa=0.0066\f$, \f$v_0=1.70\f$, porous
 * elasticity with constant \f$G=20\f$ MPa. States: normally consolidated (NC, \f$p'_0=p'_{c0}=200\f$ kPa)
 * and OCR = 5 (\f$p'_0=100\f$, \f$p'_{c0}=500\f$ kPa). Increments
 * \f$n\in\{10,15,20,25,50,100,200,400,800,1600\}\f$ up to \f$\varepsilon_a=20\%\f$.
 *  - Drained test (mcc::TriaxialDrained, constant cell pressure): error of q against the closed form of
 *    Appendix B.1 (mcc::TriaxialDrainedClosed with 20000 points, interpolated at the computed
 *    \f$\varepsilon_a\f$), at the end of the test and the largest along the path, and the mean and largest
 *    numbers of local Newton iterations per plastic projection (mcc::TLocalStats). A failed local projection
 *    means no solution (dash in Table 3); then the program checks whether the first increment alone converges.
 *  - Undrained test (mcc::TriaxialUndrained, \f$\Delta\varepsilon_{xx}=\Delta\varepsilon_{yy}=
 *    -\Delta\varepsilon_{zz}/2\f$): at the plastic states (\f$|p'-p'_0|>10^{-9}p'_0\f$) the largest
 *    \f$|p'-p'(\eta)|\f$ with the closed-form path (B.7) at the same \f$\eta=q/p'\f$ (mcc::UndrainedClosedP),
 *    and \f$p'_{end}-p'_0(R/2)^\Lambda\f$, \f$\Lambda=(\lambda-\kappa)/\lambda\f$, \f$R=p'_{c0}/p'_0\f$.
 *  - First increment of the drained NC test with 20 increments: \f$q\f$, \f$\varepsilon_v\f$,
 *    \f$p'_c/p'_{c,n}\f$ from the consistency condition (\f$p'_c=p'+q^2/(M^2p')\f$) and
 *    \f$v_0\Delta\alpha/\kappa=\frac{\lambda-\kappa}{\kappa}\ln(p'_c/p'_{c,n})\f$.
 *  - First increment with 15 and 10 increments and the frozen form: sweep of the lateral strain
 *    \f$\varepsilon_r\in[-0.03,0.01]\f$ (4001 states) to bracket the solution \f$\sigma_r=-p'_0\f$.
 *
 * The program mirrors the functions frozen() and undrained_point() of gen_data.py (Python transcription).
 * There is no finite element mesh in this example: the "structure" is a single integration point, and the
 * methods follow the sequence of the finite element examples: set-up of the material (CreateMaterial), of
 * the reference solutions (ClosedForm), solution (RunDrained, RunUndrained, FirstIncrementSweep) and
 * post-processing (PostProcess, Print).
 */
class FrozenBulkModulus {
public:
    /** @brief Integration of the porous law during the plastic correction */
    enum EIntegration {
        EExact = 0, ///< exact integral of the porous law (this work)
        EFrozen = 1 ///< bulk modulus frozen at its trial value (Sanei et al. 2020)
    };

    /** @brief Initial isotropic state of a test */
    struct TState {
        std::string fName; ///< "NC" or "OCR5"
        REAL fP0;          ///< mean effective stress p'_0 (kPa)
        REAL fPc0;         ///< preconsolidation pressure p'_c0 (kPa)
    };

    /** @brief Result of a drained test */
    struct TDrained {
        bool fConverged = false;               ///< false if a local projection failed (no solution)
        bool fFirstIncrementConverges = false; ///< for a failed test: the first increment alone converges
        std::vector<std::array<REAL, 6>> fPath; ///< rows (eps_a, p', q, eps_v, eps_q, sigma_a), compression positive
        std::vector<REAL> fQClosed;            ///< closed-form q at the eps_a of each row
        REAL fEnd = std::numeric_limits<REAL>::quiet_NaN(); ///< q_end - q_closed(eps_a,end) (kPa)
        REAL fMax = std::numeric_limits<REAL>::quiet_NaN(); ///< max |q - q_closed(eps_a)| along the path (kPa)
        mcc::TLocalStats fStats;               ///< local Newton iterations of the plastic projections
    };

    /** @brief Result of an undrained test */
    struct TUndrained {
        bool fConverged = false;               ///< false if a local projection failed
        std::vector<std::array<REAL, 3>> fPath; ///< rows (eps_a, p', q)
        std::vector<REAL> fPClosed;            ///< p' of the path (B.7) at the eta of each row (NaN at elastic states)
        int fNPlastic = 0;                     ///< number of plastic states of the path
        REAL fEnd = std::numeric_limits<REAL>::quiet_NaN(); ///< p'_end - p'_0 (R/2)^Lambda (kPa)
        REAL fMax = std::numeric_limits<REAL>::quiet_NaN(); ///< max |p' - p'(eta)| at the plastic states (kPa)
        mcc::TLocalStats fStats;               ///< local Newton iterations of the plastic projections
    };

    /** @brief Drained and undrained tests of a state, an integration and a number of increments */
    struct TRun {
        int fState = 0;                  ///< index in fStates
        EIntegration fIntegration = EExact;
        int fN = 0;                      ///< number of increments
        TDrained fDrained;
        TUndrained fUndrained;
    };

    /** @brief Analysis of the first increment of a drained test */
    struct TFirstIncrement {
        REAL fQ = 0.;          ///< deviatoric stress (kPa)
        REAL fEpsV = 0.;       ///< volumetric strain
        REAL fPcRatio = 0.;    ///< p'_c / p'_c,n from the consistency condition
        REAL fV0DalKappa = 0.; ///< v0 Delta alpha / kappa = (lambda - kappa)/kappa ln(p'_c/p'_c,n)
    };

    /** @brief Sweep of the lateral strain in the first increment of the frozen form */
    struct TSweep {
        int fN = 0;            ///< number of increments of the test (Delta eps_a = 0.2 / n)
        REAL fDea = 0.;        ///< axial strain increment
        std::vector<std::array<REAL, 5>> fRows; ///< (eps_r, converged, sigma_r + p'_0, q, p'_c/p'_c,n)
        int fFailures = 0;     ///< states where the local projection failed
        bool fBracket = false; ///< a change of sign of sigma_r + p'_0 was found
        std::array<REAL, 4> fLeft{}, fRight{}; ///< converged states (eps_r, sigma_r + p'_0, q, p'_c/p'_c,n) around it
    };

    /** @brief Values of the Python transcription (gen_data.py frozen()) of one test */
    struct TPythonReference {
        bool fAvailable = false; ///< the test is in the reference list
        bool fNone = false;      ///< the Python test failed (no solution)
        REAL fEnd = 0., fMax = 0., fIts = 0.;
        int fItMax = 0;
    };

    /** @brief One row of Table 3 of the article (a negative value is a dash: no convergence) */
    struct TArticleRow {
        int fN;
        REAL fDrainedExact, fDrainedFrozen, fItsExact, fItsFrozen, fUndrainedNC, fUndrainedOCR5;
    };

    /** @name Parameters (Sect. 6.1) */
    /** @{ */
    REAL fM = 1.2, fLambda = 0.077, fKappa = 0.0066, fV0 = 1.70;
    REAL fG = 20000.;      ///< constant shear modulus (kPa)
    REAL fNu = 0.3;        ///< Poisson ratio (not used with constant G; argument of the closed form)
    REAL fEaMax = 0.2;     ///< final axial strain
    int fClosedPoints = 20000; ///< points of the closed form of the drained test (npts of triaxial_closed)
    std::vector<int> fIncrements = {10, 15, 20, 25, 50, 100, 200, 400, 800, 1600};
    std::vector<TState> fStates = {{"NC", 200., 200.}, {"OCR5", 100., 500.}};
    /** @} */

    /** @brief Constitutive model: MCC with porous elasticity, constant G and the given integration */
    mcc::TPlastic CreateMaterial(EIntegration integ) const;

    /** @brief Closed form of the drained test, eqs. (B.1)-(B.6), rows (eps_a, p', q, eps_v, eps_q, sigma_a) */
    std::vector<std::array<REAL, 6>> ClosedForm(const TState &st, int npts) const;

    /**
     * @brief Drained test with n increments and its errors against the closed form
     * @param model constitutive model
     * @param st initial state
     * @param n number of increments
     * @param closed closed form of the drained test of this state
     */
    TDrained RunDrained(const mcc::TPlastic &model, const TState &st, int n,
                        const std::vector<std::array<REAL, 6>> &closed) const;

    /** @brief Undrained test with n increments and its errors against the closed-form path (B.7) */
    TUndrained RunUndrained(const mcc::TPlastic &model, const TState &st, int n) const;

    /** @brief First increment of a drained test: q, eps_v, p'_c/p'_c,n and v0 Delta alpha / kappa */
    TFirstIncrement FirstIncrement(const TDrained &d, const TState &st) const;

    /**
     * @brief First increment of the drained NC test with the frozen form and n increments: the lateral strain
     * eps_r is swept in [-0.03, 0.01] (4001 states, numpy.linspace) and the stress update is evaluated from
     * the initial state; the change of sign of sigma_r + p'_0 brackets the solution of the increment
     */
    TSweep FirstIncrementSweep(int n) const;

    /** @brief Runs all the tests (both integrations, both states, all the increments) */
    std::vector<TRun> RunTests() const;

    /** @brief Writes the CSV files of Table 3, of the paths, of the closed forms and of the first increment */
    void PostProcess(const std::vector<TRun> &runs, const TFirstIncrement first[2],
                     const std::vector<TSweep> &sweeps) const;

    /** @brief Prints Table 3 and the other quantities of Sect. 6.1 next to the article and Python values */
    void Print(const std::vector<TRun> &runs, const TFirstIncrement first[2], const std::vector<TSweep> &sweeps) const;

    /** @brief Runs the complete example */
    void RunAll();

    /** @brief Value of the Python transcription of a test ("drained"/"undrained", state, integration, n) */
    static TPythonReference PythonReference(const std::string &test, const std::string &state, EIntegration integ,
                                            int n);

    /** @brief Rows of Table 3 of the article */
    static const std::vector<TArticleRow> &ArticleTable3();

    /** @brief "exact" or "frozen" */
    static std::string IntegrationName(EIntegration integ) { return integ == EExact ? "exact" : "frozen"; }

private:
    /** @brief The run of a state, an integration and a number of increments (nullptr if absent) */
    static const TRun *Find(const std::vector<TRun> &runs, int state, EIntegration integ, int n);
};

// ------------------------------------------------------------------------------------------------ set-up

inline mcc::TPlastic FrozenBulkModulus::CreateMaterial(EIntegration integ) const {
    mcc::TPlastic model;
    model.SetModifiedCamClay(fM, fLambda, fKappa); // pt = 0, omega = 1
    model.SetPorousElasticity();
    model.SetConstantShearModulus(fG);
    model.SetPorousIntegration(integ == EExact ? TPZYCModifiedCamClayRHW::EExact : TPZYCModifiedCamClayRHW::EFrozen);
    model.SetDefaultSpecificVolume(fV0);
    return model;
}

inline std::vector<std::array<REAL, 6>> FrozenBulkModulus::ClosedForm(const TState &st, int npts) const {
    return mcc::TriaxialDrainedClosed(st.fP0, st.fPc0, fV0, fM, fLambda, fKappa, fG, fNu, npts);
}

// ------------------------------------------------------------------------------------------------ solution

inline FrozenBulkModulus::TDrained FrozenBulkModulus::RunDrained(const mcc::TPlastic &model, const TState &st, int n,
                                                                 const std::vector<std::array<REAL, 6>> &closed) const {
    TDrained d;
    d.fPath = mcc::TriaxialDrained(model, st.fP0, st.fPc0, fV0, fEaMax, n, &d.fStats);
    d.fConverged = !d.fPath.empty();
    if (!d.fConverged) {
        // does the projection fail already in the first increment?
        d.fFirstIncrementConverges = !mcc::TriaxialDrained(model, st.fP0, st.fPc0, fV0, fEaMax / n, 1).empty();
        return d;
    }
    d.fMax = 0.;
    for (auto &row : d.fPath) {
        const REAL qa = mcc::Interpolate(closed, row[0], 2);
        d.fQClosed.push_back(qa);
        d.fMax = std::max(d.fMax, std::fabs(row[2] - qa));
    }
    d.fEnd = d.fPath.back()[2] - d.fQClosed.back();
    return d;
}

inline FrozenBulkModulus::TUndrained FrozenBulkModulus::RunUndrained(const mcc::TPlastic &model, const TState &st,
                                                                     int n) const {
    TUndrained u;
    u.fPath = mcc::TriaxialUndrained(model, st.fP0, st.fPc0, fV0, fEaMax, n, &u.fStats);
    u.fConverged = !u.fPath.empty();
    if (!u.fConverged) return u;
    const REAL R = st.fPc0 / st.fP0, Lam = (fLambda - fKappa) / fLambda;
    u.fMax = 0.;
    for (auto &row : u.fPath) {
        if (std::fabs(row[1] - st.fP0) > 1e-9 * st.fP0) { // plastic state (p' = p'_0 on the elastic part)
            const REAL pclosed = mcc::UndrainedClosedP(st.fP0, R, fM, fLambda, fKappa, row[2] / row[1]);
            u.fPClosed.push_back(pclosed);
            u.fMax = std::max(u.fMax, std::fabs(row[1] - pclosed));
            u.fNPlastic++;
        } else {
            u.fPClosed.push_back(std::numeric_limits<REAL>::quiet_NaN());
        }
    }
    u.fEnd = u.fPath.back()[1] - st.fP0 * std::pow(R / 2., Lam);
    return u;
}

inline FrozenBulkModulus::TFirstIncrement FrozenBulkModulus::FirstIncrement(const TDrained &d, const TState &st) const {
    TFirstIncrement f;
    if (d.fPath.size() < 2) return f;
    const auto &e = d.fPath[1];
    const REAL pc = e[1] + e[2] * e[2] / (fM * fM * e[1]); // consistency condition (pt = 0, omega = 1)
    f.fQ = e[2];
    f.fEpsV = e[3];
    f.fPcRatio = pc / st.fPc0;
    f.fV0DalKappa = (fLambda - fKappa) / fKappa * std::log(pc / st.fPc0);
    return f;
}

inline FrozenBulkModulus::TSweep FrozenBulkModulus::FirstIncrementSweep(int n) const {
    const TState &st = fStates[0];
    const mcc::TPlastic model = CreateMaterial(EFrozen);
    TSweep s;
    s.fN = n;
    s.fDea = fEaMax / n;
    const mcc::TPointState initial(mcc::IsotropicTensor(-st.fP0), st.fPc0, fV0);
    const TPZTensor<REAL> epsn; // the increment starts from zero strain
    const int npts = 4001;
    const REAL start = -0.03, stop = 0.01, step = (stop - start) / (npts - 1);
    TPZFNMatrix<36, REAL> D(6, 6, 0.);
    for (int i = 0; i < npts; ++i) {
        const REAL er = (i == npts - 1) ? stop : REAL(i) * step + start; // numpy.linspace
        TPZTensor<REAL> eps, sig;
        eps.ZZ() = -s.fDea;
        eps.XX() = er;
        eps.YY() = er;
        REAL pc = st.fPc0;
        int type = 0;
        if (mcc::ApplyStrain(model, epsn, initial, eps, sig, D, pc, type)) {
            s.fRows.push_back({er, 1., sig.XX() + st.fP0, mcc::DeviatoricStress(sig), pc / st.fPc0});
        } else {
            const REAL nan = std::numeric_limits<REAL>::quiet_NaN();
            s.fRows.push_back({er, 0., nan, nan, nan});
            s.fFailures++;
        }
    }
    // first change of sign of sigma_r + p'_0 between consecutive converged states
    auto sign = [](REAL x) { return (x > 0.) - (x < 0.); };
    const std::array<REAL, 5> *prev = nullptr;
    for (auto &r : s.fRows) {
        if (r[1] == 0.) continue;
        if (prev && sign((*prev)[2]) != sign(r[2])) {
            s.fBracket = true;
            s.fLeft = {(*prev)[0], (*prev)[2], (*prev)[3], (*prev)[4]};
            s.fRight = {r[0], r[2], r[3], r[4]};
            break;
        }
        prev = &r;
    }
    return s;
}

inline std::vector<FrozenBulkModulus::TRun> FrozenBulkModulus::RunTests() const {
    std::vector<TRun> runs;
    for (EIntegration integ : {EExact, EFrozen}) {
        const mcc::TPlastic model = CreateMaterial(integ);
        for (int is = 0; is < int(fStates.size()); ++is) {
            const TState &st = fStates[is];
            const auto closed = ClosedForm(st, fClosedPoints);
            for (int n : fIncrements) {
                TRun r;
                r.fState = is;
                r.fIntegration = integ;
                r.fN = n;
                r.fDrained = RunDrained(model, st, n, closed);
                r.fUndrained = RunUndrained(model, st, n);
                runs.push_back(std::move(r));
            }
        }
    }
    return runs;
}

inline const FrozenBulkModulus::TRun *FrozenBulkModulus::Find(const std::vector<TRun> &runs, int state,
                                                              EIntegration integ, int n) {
    for (auto &r : runs)
        if (r.fState == state && r.fIntegration == integ && r.fN == n) return &r;
    return nullptr;
}

// ------------------------------------------------------------------------------------------------ post-processing

inline void FrozenBulkModulus::PostProcess(const std::vector<TRun> &runs, const TFirstIncrement first[2],
                                           const std::vector<TSweep> &sweeps) const {
    const REAL nan = std::numeric_limits<REAL>::quiet_NaN();
    // Table 3 (all the increments of the sweep)
    {
        std::vector<std::vector<REAL>> rows;
        for (int n : fIncrements) {
            const TRun *ne = Find(runs, 0, EExact, n), *nf = Find(runs, 0, EFrozen, n);
            const TRun *oe = Find(runs, 1, EExact, n), *of = Find(runs, 1, EFrozen, n);
            rows.push_back({REAL(n), ne->fDrained.fMax, nf->fDrained.fMax,
                            ne->fDrained.fConverged ? ne->fDrained.fStats.Mean() : nan,
                            nf->fDrained.fConverged ? nf->fDrained.fStats.Mean() : nan, nf->fUndrained.fMax,
                            of->fUndrained.fMax, ne->fUndrained.fMax, oe->fUndrained.fMax});
        }
        mcc::WriteCSV("frozen_table3.csv",
                      {"n", "drained_NC_exact_max_err_q", "drained_NC_frozen_max_err_q", "drained_NC_exact_mean_its",
                       "drained_NC_frozen_mean_its", "undrained_NC_frozen_max_err_p", "undrained_OCR5_frozen_max_err_p",
                       "undrained_NC_exact_max_err_p", "undrained_OCR5_exact_max_err_p"},
                      rows);
    }
    // every test with the values of the Python transcription
    {
        std::vector<std::vector<REAL>> rows;
        for (auto &r : runs) {
            for (int test = 0; test < 2; ++test) {
                const bool drained = (test == 0);
                const TPythonReference py =
                    PythonReference(drained ? "drained" : "undrained", fStates[r.fState].fName, r.fIntegration, r.fN);
                const bool conv = drained ? r.fDrained.fConverged : r.fUndrained.fConverged;
                const mcc::TLocalStats &stats = drained ? r.fDrained.fStats : r.fUndrained.fStats;
                const bool pyok = py.fAvailable && !py.fNone;
                rows.push_back({REAL(test), REAL(r.fState), REAL(r.fIntegration), REAL(r.fN), REAL(conv),
                                drained ? r.fDrained.fEnd : r.fUndrained.fEnd,
                                drained ? r.fDrained.fMax : r.fUndrained.fMax, conv ? stats.Mean() : nan,
                                conv ? REAL(stats.fMax) : nan, pyok ? py.fEnd : nan, pyok ? py.fMax : nan,
                                pyok ? py.fIts : nan, pyok ? REAL(py.fItMax) : nan});
            }
        }
        mcc::WriteCSV("frozen_runs.csv",
                      {"test(0=drained;1=undrained)", "state(0=NC;1=OCR5)", "integration(0=exact;1=frozen)", "n",
                       "converged", "end_error", "max_error", "mean_local_its", "max_local_its", "python_end_error",
                       "python_max_error", "python_mean_local_its", "python_max_local_its"},
                      rows);
    }
    // paths, one file per test, state and integration (column n identifies the number of increments)
    for (int is = 0; is < int(fStates.size()); ++is) {
        for (EIntegration integ : {EExact, EFrozen}) {
            const std::string tag = fStates[is].fName + "_" + IntegrationName(integ);
            std::vector<std::vector<REAL>> drows, urows;
            for (int n : fIncrements) {
                const TRun *r = Find(runs, is, integ, n);
                const TDrained &d = r->fDrained;
                for (size_t i = 0; i < d.fPath.size(); ++i) {
                    const auto &w = d.fPath[i];
                    drows.push_back({REAL(n), w[0], w[1], w[2], w[3], w[4], w[5], d.fQClosed[i], w[2] - d.fQClosed[i]});
                }
                const TUndrained &u = r->fUndrained;
                for (size_t i = 0; i < u.fPath.size(); ++i) {
                    const auto &w = u.fPath[i];
                    const REAL eta = w[2] / w[1];
                    const bool plastic = !std::isnan(u.fPClosed[i]);
                    urows.push_back({REAL(n), w[0], w[1], w[2], eta, REAL(plastic), u.fPClosed[i],
                                     plastic ? w[1] - u.fPClosed[i] : nan});
                }
            }
            mcc::WriteCSV("frozen_drained_" + tag + ".csv",
                          {"n", "eps_a", "p_eff", "q", "eps_v", "eps_q", "sigma_a", "q_closed", "q_error"}, drows);
            mcc::WriteCSV("frozen_undrained_" + tag + ".csv",
                          {"n", "eps_a", "p_eff", "q", "eta", "plastic", "p_closed_B7", "p_error"}, urows);
        }
        // closed forms for plotting: drained (B.1)-(B.6) with 600 points and undrained path (B.7)
        std::vector<std::vector<REAL>> crows;
        for (auto &w : ClosedForm(fStates[is], 600)) crows.push_back({w[0], w[1], w[2], w[3], w[4], w[5]});
        mcc::WriteCSV("frozen_closed_drained_" + fStates[is].fName + ".csv",
                      {"eps_a", "p_eff", "q", "eps_v", "eps_q", "sigma_a"}, crows);
        std::vector<std::vector<REAL>> brows;
        const REAL R = fStates[is].fPc0 / fStates[is].fP0;
        for (int i = 0; i < 400; ++i) {
            const REAL eta = fM * i / 400.;
            const REAL p = mcc::UndrainedClosedP(fStates[is].fP0, R, fM, fLambda, fKappa, eta);
            brows.push_back({eta, p, eta * p});
        }
        mcc::WriteCSV("frozen_closed_undrained_" + fStates[is].fName + ".csv", {"eta", "p_eff", "q"}, brows);
    }
    // first increment of the drained NC test with 20 increments
    {
        const auto closed = ClosedForm(fStates[0], fClosedPoints);
        const REAL de = fEaMax / 20.;
        std::vector<std::vector<REAL>> rows;
        rows.push_back({0., first[0].fQ, first[0].fEpsV, first[0].fPcRatio, first[0].fV0DalKappa});
        rows.push_back({1., first[1].fQ, first[1].fEpsV, first[1].fPcRatio, first[1].fV0DalKappa});
        rows.push_back({2., mcc::Interpolate(closed, de, 2), mcc::Interpolate(closed, de, 3), nan, nan});
        mcc::WriteCSV("frozen_first_increment.csv",
                      {"integration(0=exact;1=frozen;2=closed_form)", "q", "eps_v", "pc_over_pcn", "v0_dal_over_kappa"},
                      rows);
    }
    for (auto &s : sweeps) {
        std::vector<std::vector<REAL>> rows;
        for (auto &r : s.fRows) rows.push_back({r[0], r[1], r[2], r[3], r[4]});
        mcc::WriteCSV("frozen_first_sweep_n" + std::to_string(s.fN) + ".csv",
                      {"eps_r", "converged", "sigma_r_plus_p0", "q", "pc_over_pcn"}, rows);
    }
}

inline void FrozenBulkModulus::Print(const std::vector<TRun> &runs, const TFirstIncrement first[2],
                                     const std::vector<TSweep> &sweeps) const {
    std::ostringstream os;
    // a value with a fixed number of decimals, or a dash for no solution
    auto fmt = [](REAL v, int prec) {
        std::ostringstream s;
        if (std::isnan(v)) s << "-";
        else s << std::fixed << std::setprecision(prec) << v;
        return s.str();
    };
    // "value [reference]" in a column of the given width
    auto cell = [&fmt](REAL v, REAL ref, int prec, int width) {
        std::ostringstream s;
        s << std::setw(width) << (fmt(v, prec) + " [" + fmt(ref, prec) + "]");
        return s.str();
    };
    auto ref = [](const TPythonReference &py, bool maxerr) {
        if (!py.fAvailable || py.fNone) return std::numeric_limits<REAL>::quiet_NaN();
        return maxerr ? py.fMax : py.fEnd;
    };
    auto refits = [](const TPythonReference &py) {
        if (!py.fAvailable || py.fNone) return std::numeric_limits<REAL>::quiet_NaN();
        return py.fIts;
    };
    const REAL nan = std::numeric_limits<REAL>::quiet_NaN();

    std::cout << "\nTable 3 - exact integration of the porous law and bulk modulus frozen at its trial value\n"
              << "RS2 clay (M = 1.2, lambda = 0.077, kappa = 0.0066, v0 = 1.70, G = 20 MPa), triaxial tests to eps_a = 20%\n"
              << "this work [reference]: article Table 3; rows 15 and 25 (not in the article): Python gen_data.py frozen()\n"
              << "  '-' : no convergence in the first increment\n\n";
    std::cout << "        |            Drained, NC                                       |    Undrained, frozen\n"
              << "        |  largest error in q (kPa)      |  local iterations           |  largest error in p' (kPa)\n"
              << "      n |      exact          frozen     |     exact        frozen     |       NC             OCR = 5\n";
    for (int n : fIncrements) {
        const TRun *ne = Find(runs, 0, EExact, n), *nf = Find(runs, 0, EFrozen, n), *of = Find(runs, 1, EFrozen, n);
        TArticleRow art{n, nan, nan, nan, nan, nan, nan};
        bool inArticle = false;
        for (auto &a : ArticleTable3()) {
            if (a.fN == n) {
                art = a;
                inArticle = true;
            }
        }
        if (inArticle) {
            for (REAL *v : {&art.fDrainedExact, &art.fDrainedFrozen, &art.fItsExact, &art.fItsFrozen, &art.fUndrainedNC,
                            &art.fUndrainedOCR5})
                if (*v < 0.) *v = nan;
        } else {
            art.fDrainedExact = ref(PythonReference("drained", "NC", EExact, n), true);
            art.fDrainedFrozen = ref(PythonReference("drained", "NC", EFrozen, n), true);
            art.fItsExact = refits(PythonReference("drained", "NC", EExact, n));
            art.fItsFrozen = refits(PythonReference("drained", "NC", EFrozen, n));
            art.fUndrainedNC = ref(PythonReference("undrained", "NC", EFrozen, n), true);
            art.fUndrainedOCR5 = ref(PythonReference("undrained", "OCR5", EFrozen, n), true);
        }
        const REAL itse = ne->fDrained.fConverged ? ne->fDrained.fStats.Mean() : nan;
        const REAL itsf = nf->fDrained.fConverged ? nf->fDrained.fStats.Mean() : nan;
        std::cout << std::setw(7) << n << (inArticle ? " |" : "*|") << cell(ne->fDrained.fMax, art.fDrainedExact, 2, 15)
                  << cell(nf->fDrained.fMax, art.fDrainedFrozen, 2, 16) << " |" << cell(itse, art.fItsExact, 2, 13)
                  << cell(itsf, art.fItsFrozen, 2, 14) << " |" << cell(nf->fUndrained.fMax, art.fUndrainedNC, 3, 15)
                  << cell(of->fUndrained.fMax, art.fUndrainedOCR5, 3, 16) << "\n";
    }
    std::cout << "  (*) reference from the Python transcription (row not printed in the article)\n";

    // failures of the frozen form
    std::cout << "\nDrained tests without solution (article: frozen form, 15 increments or fewer, no convergence in the "
                 "first increment):\n";
    bool anyfail = false;
    for (auto &r : runs) {
        if (r.fDrained.fConverged) continue;
        anyfail = true;
        const TPythonReference py = PythonReference("drained", fStates[r.fState].fName, r.fIntegration, r.fN);
        std::cout << "  " << fStates[r.fState].fName << " " << IntegrationName(r.fIntegration) << " n = " << r.fN
                  << ": local projection failed; the first increment alone "
                  << (r.fDrained.fFirstIncrementConverges ? "converges" : "does not converge")
                  << "   [Python: " << (py.fNone ? "no solution, first increment does not converge" : "solution") << "]\n";
    }
    if (!anyfail) std::cout << "  none\n";

    // other quantities of the drained tests quoted in the text of Sect. 6.1
    auto py = [&ref](const char *test, const char *state, EIntegration integ, int n, bool maxerr) {
        return ref(PythonReference(test, state, integ, n), maxerr);
    };
    std::cout << "\nDrained tests, other errors in q (kPa), this work [Python]:\n"
              << "        |         NC, error at eps_a = 20%        |  OCR = 5, largest error along the path  |"
                 "     OCR = 5, error at eps_a = 20%\n"
              << "      n |         exact               frozen      |         exact               frozen      |"
                 "         exact               frozen\n";
    for (int n : fIncrements) {
        const TRun *ne = Find(runs, 0, EExact, n), *nf = Find(runs, 0, EFrozen, n);
        const TRun *oe = Find(runs, 1, EExact, n), *of = Find(runs, 1, EFrozen, n);
        std::cout << std::setw(7) << n << " |" << cell(ne->fDrained.fEnd, py("drained", "NC", EExact, n, false), 4, 20)
                  << cell(nf->fDrained.fEnd, py("drained", "NC", EFrozen, n, false), 4, 20) << " |"
                  << cell(oe->fDrained.fMax, py("drained", "OCR5", EExact, n, true), 4, 20)
                  << cell(of->fDrained.fMax, py("drained", "OCR5", EFrozen, n, true), 4, 20) << " |"
                  << cell(oe->fDrained.fEnd, py("drained", "OCR5", EExact, n, false), 4, 20)
                  << cell(of->fDrained.fEnd, py("drained", "OCR5", EFrozen, n, false), 4, 20) << "\n";
    }

    // undrained tests: exact integration on the closed-form path (B.7); end values of the frozen form
    auto sci = [](REAL v, REAL r, int width) {
        std::ostringstream s, t;
        s << std::scientific << std::setprecision(1) << v << " [" << r << "]";
        t << std::setw(width) << s.str();
        return t.str();
    };
    std::cout << "\nUndrained tests (kPa), this work [Python]: largest |p' - p'(B.7)| with the exact integration "
                 "(article: < 1e-10)\nand p'_end - p'_0 (R/2)^Lambda with the frozen modulus\n"
              << "        |    exact, largest error in p'       |   frozen, p' error at eps_a = 20%\n"
              << "      n |               NC            OCR = 5 |                NC             OCR = 5\n";
    REAL maxexact[2] = {0., 0.};
    for (int n : fIncrements) {
        const TRun *ne = Find(runs, 0, EExact, n), *oe = Find(runs, 1, EExact, n);
        const TRun *nf = Find(runs, 0, EFrozen, n), *of = Find(runs, 1, EFrozen, n);
        maxexact[0] = std::max(maxexact[0], ne->fUndrained.fMax);
        maxexact[1] = std::max(maxexact[1], oe->fUndrained.fMax);
        std::cout << std::setw(7) << n << " |" << sci(ne->fUndrained.fMax, py("undrained", "NC", EExact, n, true), 18)
                  << sci(oe->fUndrained.fMax, py("undrained", "OCR5", EExact, n, true), 19) << " |"
                  << cell(nf->fUndrained.fEnd, py("undrained", "NC", EFrozen, n, false), 4, 18)
                  << cell(of->fUndrained.fEnd, py("undrained", "OCR5", EFrozen, n, false), 4, 20) << "\n";
    }
    std::cout << std::scientific << std::setprecision(1) << "  largest over all n: NC " << maxexact[0] << ", OCR5 "
              << maxexact[1] << " kPa   [article: < 1e-10 kPa; Python: 2.0e-11, 6.5e-11]\n"
              << std::defaultfloat;

    // ratios frozen / exact
    REAL rnc[2] = {1e300, -1e300}, rocm[2] = {1e300, -1e300}, roce[2] = {1e300, -1e300}, dits = 0.;
    for (int n : fIncrements) {
        const TRun *ne = Find(runs, 0, EExact, n), *nf = Find(runs, 0, EFrozen, n);
        const TRun *oe = Find(runs, 1, EExact, n), *of = Find(runs, 1, EFrozen, n);
        if (ne->fDrained.fConverged && nf->fDrained.fConverged) {
            const REAL a = nf->fDrained.fMax / ne->fDrained.fMax;
            rnc[0] = std::min(rnc[0], a);
            rnc[1] = std::max(rnc[1], a);
            if (n >= 50)
                dits = std::max(dits, std::fabs(nf->fDrained.fStats.Mean() / ne->fDrained.fStats.Mean() - 1.));
        }
        if (oe->fDrained.fConverged && of->fDrained.fConverged) {
            const REAL a = of->fDrained.fMax / oe->fDrained.fMax, b = of->fDrained.fEnd / oe->fDrained.fEnd;
            rocm[0] = std::min(rocm[0], a);
            rocm[1] = std::max(rocm[1], a);
            roce[0] = std::min(roce[0], b);
            roce[1] = std::max(roce[1], b);
        }
    }
    std::cout << std::fixed << std::setprecision(2) << "\nFrozen versus exact integration, drained tests (Sect. 6.1):\n"
              << "  drained NC, largest error frozen / exact: " << rnc[0] << " to " << rnc[1]
              << "   [article: 2.4 to 5; Python: 2.41 to 5.18]\n"
              << "  drained NC, n >= 50: mean local iterations differ by at most " << 100. * dits
              << "%   [article: less than 3%; Python: 2.43%]\n"
              << "  drained OCR5, largest error of the frozen form " << 100. * (rocm[0] - 1.) << "% to "
              << 100. * (rocm[1] - 1.) << "% larger   [article: 10 to 20%; Python: 9.52% to 20.00%]\n"
              << "  drained OCR5, error at 20% of the frozen form " << 100. * (roce[0] - 1.) << "% to "
              << 100. * (roce[1] - 1.) << "% larger   [article: 4 to 6%; Python: 3.95% to 6.19%]\n"
              << std::defaultfloat;

    // first increment of the drained NC test with 20 increments
    const auto closed = ClosedForm(fStates[0], fClosedPoints);
    const REAL de = fEaMax / 20.;
    const REAL bound = std::exp(fKappa / (fLambda - fKappa));
    std::cout << "\nFirst increment of the drained NC test with 20 increments (Delta eps_a = 1%), this work "
                 "[article; Python]:\n"
              << std::fixed << std::setprecision(4) << "  closed form at eps_a = 1%: q = " << std::setprecision(2)
              << mcc::Interpolate(closed, de, 2) << " kPa [100.9; 100.94], eps_v = " << std::setprecision(5)
              << mcc::Interpolate(closed, de, 3) << " [0.0121; 0.01209]\n";
    const char *aq[2] = {"85.0", "41.8"}, *pyq[2] = {"85.00", "41.82"};
    const char *aev[2] = {"-", "0.0247"}, *pyev[2] = {"0.00981", "0.02467"};
    const char *apc[2] = {"1.25", "1.098"}, *pypc[2] = {"1.2515", "1.0981"};
    const char *adal[2] = {"2.39", "0.998"}, *pydal[2] = {"2.3931", "0.9981"};
    for (int i = 0; i < 2; ++i) {
        std::cout << "  " << std::setw(6) << IntegrationName(EIntegration(i)) << ": q = " << std::setprecision(2)
                  << first[i].fQ << " kPa [" << aq[i] << "; " << pyq[i] << "], eps_v = " << std::setprecision(5)
                  << first[i].fEpsV << " [" << aev[i] << "; " << pyev[i] << "], p'c/p'c,n = " << std::setprecision(4)
                  << first[i].fPcRatio << " [" << apc[i] << "; " << pypc[i] << "], v0 dal/kappa = " << first[i].fV0DalKappa
                  << " [" << adal[i] << "; " << pydal[i] << "]\n";
    }
    std::cout << "  bound of the frozen form p'c/p'c,n < exp(kappa/(lambda-kappa)) = " << bound << " [1.098; 1.0983]\n";

    // sweeps of the first increment with 15 and 10 increments
    std::cout << "\nFrozen form, first increment with 15 and 10 increments: sweep of eps_r in [-0.03, 0.01] (4001 states), "
                 "this work [Python]\n";
    struct TSweepRef {
        int fN, fFail;
        REAL fE0, fE1, fQ0, fQ1, fPc0, fPc1;
    };
    const TSweepRef sref[2] = {{15, 1847, -0.01067, -0.01065, 41.7352, 42.0285, 1.098271, 1.098270},
                               {10, 2212, -0.01936, -0.01523, 10.3333, 69.7412, 1.098285, 1.098285}};
    for (auto &s : sweeps) {
        const TSweepRef *r = nullptr;
        for (auto &x : sref)
            if (x.fN == s.fN) r = &x;
        std::cout << "  n = " << s.fN << ": projection failed in " << s.fFailures << " of " << s.fRows.size()
                  << " states [" << (r ? std::to_string(r->fFail) : "?") << "; round-off sensitive, see README]";
        if (!s.fBracket) {
            std::cout << "; no change of sign\n";
            continue;
        }
        std::cout << "; sigma_r + p'0 changes sign between eps_r = " << std::setprecision(5) << s.fLeft[0] << " and "
                  << s.fRight[0];
        if (r) std::cout << " [" << r->fE0 << ", " << r->fE1 << "]";
        std::cout << "\n      q = " << std::setprecision(4) << s.fLeft[2] << " and " << s.fRight[2];
        if (r) std::cout << " [" << r->fQ0 << ", " << r->fQ1 << "]";
        std::cout << ", p'c/p'c,n = " << std::setprecision(6) << s.fLeft[3] << " and " << s.fRight[3];
        if (r) std::cout << " [" << r->fPc0 << ", " << r->fPc1 << "]";
        std::cout << ", eps_v = " << std::setprecision(4) << s.fDea - 2. * s.fLeft[0] << " and "
                  << s.fDea - 2. * s.fRight[0] << "\n";
    }
    std::cout << "  (article: with 15 increments the solution lies at the bound, q = 41.9 kPa and eps_v = 3.5%)\n"
              << std::defaultfloat;

    // agreement with the Python transcription
    REAL dend = 0., dmax = 0., dits2 = 0., dund = 0.;
    int ditmax = 0, nmissing = 0, nsame = 0, ncmp = 0;
    for (auto &r : runs) {
        for (int test = 0; test < 2; ++test) {
            const bool drained = (test == 0);
            const TPythonReference py =
                PythonReference(drained ? "drained" : "undrained", fStates[r.fState].fName, r.fIntegration, r.fN);
            const bool conv = drained ? r.fDrained.fConverged : r.fUndrained.fConverged;
            if (!py.fAvailable) {
                nmissing++;
                continue;
            }
            ncmp++;
            if (py.fNone || !conv) {
                if (py.fNone == !conv) nsame++;
                continue;
            }
            nsame++;
            const REAL end = drained ? r.fDrained.fEnd : r.fUndrained.fEnd;
            const REAL mx = drained ? r.fDrained.fMax : r.fUndrained.fMax;
            const mcc::TLocalStats &st = drained ? r.fDrained.fStats : r.fUndrained.fStats;
            if (!drained && r.fIntegration == EExact) {
                // round-off level values: absolute difference
                dund = std::max({dund, std::fabs(end - py.fEnd), std::fabs(mx - py.fMax)});
            } else {
                dend = std::max(dend, std::fabs(end - py.fEnd) / std::max(std::fabs(py.fEnd), REAL(1e-300)));
                dmax = std::max(dmax, std::fabs(mx - py.fMax) / std::max(std::fabs(py.fMax), REAL(1e-300)));
            }
            dits2 = std::max(dits2, std::fabs(st.Mean() - py.fIts));
            ditmax = std::max(ditmax, std::abs(st.fMax - py.fItMax));
        }
    }
    std::cout << std::scientific << std::setprecision(2) << "\nAgreement with the Python transcription (" << ncmp
              << " tests, file frozen_runs.csv): same convergence in " << nsame << " tests\n"
              << "  largest relative difference: end error " << dend << ", largest error " << dmax
              << " (undrained exact, absolute: " << dund << " kPa)\n"
              << "  largest difference of the mean local iterations " << dits2 << ", of the maximum " << ditmax << "\n"
              << std::defaultfloat;
    if (nmissing) std::cout << "  " << nmissing << " tests without Python reference\n";
}

inline void FrozenBulkModulus::RunAll() {
    const auto t0 = std::chrono::steady_clock::now();
    std::cout << "FrozenBulkModulus: Sect. 6.1, Table 3 of the article (exact versus frozen integration of the porous law)"
              << std::endl;
    const std::vector<TRun> runs = RunTests();

    TFirstIncrement first[2];
    for (EIntegration integ : {EExact, EFrozen}) first[integ] = FirstIncrement(Find(runs, 0, integ, 20)->fDrained, fStates[0]);
    std::vector<TSweep> sweeps = {FirstIncrementSweep(15), FirstIncrementSweep(10)};

    Print(runs, first, sweeps);
    PostProcess(runs, first, sweeps);
    const REAL secs = std::chrono::duration<REAL>(std::chrono::steady_clock::now() - t0).count();
    std::cout << "\nFiles: frozen_table3.csv, frozen_runs.csv, frozen_{drained,undrained}_{NC,OCR5}_{exact,frozen}.csv,\n"
              << "       frozen_closed_{drained,undrained}_{NC,OCR5}.csv, frozen_first_increment.csv, "
                 "frozen_first_sweep_n{15,10}.csv\n"
              << "Run time: " << std::fixed << std::setprecision(2) << secs << " s" << std::defaultfloat << std::endl;
}

// ------------------------------------------------------------------------------------------------ reference values

inline const std::vector<FrozenBulkModulus::TArticleRow> &FrozenBulkModulus::ArticleTable3() {
    // Table 3 of the article (-1: dash, no convergence in the first increment)
    static const std::vector<TArticleRow> rows = {
        {10, 31.07, -1., 9.05, -1., 1.340, 5.161},     {20, 18.15, 94.06, 7.20, 8.98, 1.135, 3.301},
        {50, 8.47, 28.11, 5.72, 5.86, 0.674, 2.383},   {100, 4.59, 12.73, 5.16, 5.15, 0.401, 1.400},
        {200, 2.43, 6.23, 4.35, 4.33, 0.234, 0.769},   {400, 1.26, 3.11, 4.05, 4.04, 0.127, 0.405},
        {800, 0.64, 1.56, 3.69, 3.65, 0.066, 0.208},   {1600, 0.32, 0.78, 3.19, 3.12, 0.034, 0.106}};
    return rows;
}

inline FrozenBulkModulus::TPythonReference FrozenBulkModulus::PythonReference(const std::string &test,
                                                                              const std::string &state,
                                                                              EIntegration integ, int n) {
    // gen_data.py frozen(): test state integration n end max mean_its max_its ("none": no solution)
    static const char *table = R"(
drained NC exact 10 -9.39763053418443 31.07422486195506 9.045454545454545 11
undrained NC exact 10 8.87611690814083e-08 5.684341886080802e-14 8.0 8
drained NC exact 15 -6.222098948169219 22.823612784084702 7.9 10
undrained NC exact 15 5.50627987649932e-10 2.2737367544323206e-13 7.066666666666666 8
drained NC exact 20 -4.65244082037708 18.15235854071122 7.197530864197531 9
undrained NC exact 20 1.0700773600547109e-11 8.526512829121202e-14 7.0 7
drained NC exact 25 -3.7170252167454123 15.180634000987425 6.777777777777778 9
undrained NC exact 25 4.973799150320701e-13 4.263256414560601e-14 7.0 7
drained NC exact 50 -1.8588713088274744 8.470467061278612 5.718232044198895 7
undrained NC exact 50 -1.4210854715202004e-14 5.684341886080802e-14 6.0 6
drained NC exact 100 -0.9326305953787823 4.590535170709046 5.158730158730159 6
undrained NC exact 100 -1.4210854715202004e-14 1.4210854715202004e-12 5.0 5
drained NC exact 200 -0.46800966543452205 2.4261716632349817 4.34789644012945 5
undrained NC exact 200 -4.263256414560601e-14 9.947598300641403e-14 4.995 5
drained NC exact 400 -0.23462053258384685 1.2560265319373372 4.051070840197694 5
undrained NC exact 400 -2.8421709430404007e-13 2.8421709430404007e-13 3.9975 4
drained NC exact 800 -0.11749682263018713 0.6406176463491988 3.6898464163822524 4
undrained NC exact 800 -4.405364961712621e-13 2.0747847884194925e-12 3.99625 4
drained NC exact 1600 -0.058799686176712385 0.3237433884126233 3.1887742260256733 4
undrained NC exact 1600 -6.0396132539608516e-12 2.020783540501725e-11 3.98875 4
drained OCR5 exact 10 1.6923763855269272 4.656380744887969 7.046511627906977 9
undrained OCR5 exact 10 -3.263662051722349e-07 2.842170943040401e-14 7.0 7
drained OCR5 exact 15 1.1197704678779985 3.277065111892796 6.885245901639344 7
undrained OCR5 exact 15 -2.2866686322231544e-09 6.508571459562518e-11 6.066666666666666 7
drained OCR5 exact 20 0.8377618771264679 2.54238259923369 6.3076923076923075 7
undrained OCR5 exact 20 -4.988010005035903e-11 2.8421709430404007e-13 6.0 6
drained OCR5 exact 25 0.670222035497062 2.090834926440664 5.947368421052632 6
undrained OCR5 exact 25 -3.637978807091713e-12 4.035882739117369e-12 5.96 6
drained OCR5 exact 50 0.33291459672162205 1.0873413438485215 5.446428571428571 6
undrained OCR5 exact 50 -2.0179413695586845e-12 4.035882739117369e-12 5.0 5
drained OCR5 exact 100 0.16621143673575034 0.5593626763038912 4.98019801980198 5
undrained OCR5 exact 100 -1.7053025658242404e-13 1.9895196601282805e-13 5.0 5
drained OCR5 exact 200 0.08287107080070655 0.2818065156957914 3.998299319727891 4
undrained OCR5 exact 200 -1.2505552149377763e-12 1.2221335055073723e-12 4.0 4
drained OCR5 exact 400 0.041383802950434756 0.14154761053515585 3.9991445680068436 4
undrained OCR5 exact 400 -4.547473508864641e-13 4.547473508864641e-13 4.0 4
drained OCR5 exact 800 0.02068705970370388 0.07104230388188171 3.999524262607041 4
undrained OCR5 exact 800 -9.947598300641403e-13 8.242295734817162e-13 4.0 4
drained OCR5 exact 1600 0.010339572289922216 0.03555424732155643 3.0 3
undrained OCR5 exact 1600 -3.979039320256561e-12 3.154809746774845e-12 3.0 3
drained NC frozen 10 none
undrained NC frozen 10 -1.2665790200311022 1.3398045830472967 8.0 8
drained NC frozen 15 none
undrained NC frozen 15 -1.1486117927259585 1.236792178072335 7.066666666666666 8
drained NC frozen 20 -15.91851989062451 94.06276884088584 8.975903614457831 14
undrained NC frozen 20 -1.04325043631259 1.1354031817685666 7.0 7
drained NC frozen 25 -10.280326755010265 70.88540328824445 7.829787234042553 12
undrained NC frozen 25 -0.949466711708638 1.0366360163370416 7.0 7
drained NC frozen 50 -3.5351563575640057 28.11410014737561 5.857142857142857 8
undrained NC frozen 50 -0.623975639633386 0.6737827062813722 6.0 6
drained NC frozen 100 -1.5592861676312282 12.73251851830446 5.153125 6
undrained NC frozen 100 -0.3714885014456968 0.40063426294365456 5.0 5
drained NC frozen 200 -0.7458741089332079 6.229339170639776 4.332797427652733 5
undrained NC frozen 200 -0.21676661094960537 0.2339847725857851 4.995 5
drained NC frozen 400 -0.36643545073235373 3.111304108971396 4.040329218106996 5
undrained NC frozen 400 -0.11778114047773158 0.12722145663202866 3.9975 4
drained NC frozen 800 -0.18183407142976193 1.5598843532506805 3.654452492543673 4
undrained NC frozen 800 -0.06143347678886357 0.06635620933825237 3.995 4
drained NC frozen 1600 -0.0906025734037712 0.781752711517214 3.1191853155644957 4
undrained NC frozen 1600 -0.03138124287363553 0.033895098620732256 3.989375 4
drained OCR5 frozen 10 1.759235888745735 5.099566646055848 7.0476190476190474 9
undrained OCR5 frozen 10 -5.161262417861025 5.161261966253932 7.0 7
drained OCR5 frozen 15 1.1683937954905161 3.6202171740002314 6.901960784313726 7
undrained OCR5 frozen 15 -4.090199505156619 4.090199502100063 6.133333333333334 7
drained OCR5 frozen 20 0.8764165674194544 2.854516965263457 5.985294117647059 7
undrained OCR5 frozen 20 -3.301402244173204 3.301402244108374 6.0 6
drained OCR5 frozen 25 0.7027056981512771 2.368588305772761 5.941860465116279 6
undrained OCR5 frozen 25 -2.7504373860056717 2.7504373860027442 5.96 6
drained OCR5 frozen 50 0.35091132714811124 1.2605896991550765 5.466257668711656 6
undrained OCR5 frozen 50 -2.3828570009573866 2.382857000957415 5.0 5
drained OCR5 frozen 100 0.17586608553691008 0.659684473426779 4.97986577181208 5
undrained OCR5 frozen 100 -1.3995049461584586 1.3995049461585722 5.0 5
drained OCR5 frozen 200 0.08783776468825977 0.33516893477059284 3.998299319727891 4
undrained OCR5 frozen 200 -0.7691721478092575 0.7691721478092575 4.0 4
drained OCR5 frozen 400 0.04390785339967351 0.16915473928887081 3.9991445680068436 4
undrained OCR5 frozen 400 -0.4051714263261772 0.40517142632609193 4.0 4
drained OCR5 frozen 800 0.021962981258752734 0.08514803842280116 3.999524262607041 4
undrained OCR5 frozen 800 -0.20824252109366626 0.20824252109352415 4.0 4
drained OCR5 frozen 1600 0.010979989719402283 0.042665071783829944 3.0 3
undrained OCR5 frozen 1600 -0.10560815342657293 0.10560815342580554 3.0 3
)";
    TPythonReference ref;
    std::istringstream in(table);
    std::string line;
    while (std::getline(in, line)) {
        std::istringstream ls(line);
        std::string t, s, i, rest;
        int nn = 0;
        if (!(ls >> t >> s >> i >> nn)) continue;
        if (t != test || s != state || i != IntegrationName(integ) || nn != n) continue;
        ref.fAvailable = true;
        ls >> rest;
        if (rest == "none") {
            ref.fNone = true;
            return ref;
        }
        ref.fEnd = std::stod(rest);
        ls >> ref.fMax >> ref.fIts >> ref.fItMax;
        return ref;
    }
    return ref;
}
