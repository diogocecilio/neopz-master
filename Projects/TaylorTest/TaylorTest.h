/**
 * @file TaylorTest.h
 * @brief Sect. 4.5 and Fig. 3 of the article: Taylor test of the consistent tangent operator of the
 * Modified Cam-Clay return mapping in rotated Haigh-Westergaard space (TPZPlasticStepModifiedCamClay).
 */
#pragma once

#include "MCCPaperTools.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <random>
#include <string>
#include <vector>
#include <array>
#include <iostream>
#include <iomanip>
#include <sstream>

/**
 * @brief Taylor test (23) of the consistent tangent at elastic, subcritical and supercritical states.
 *
 * The test compares the stress update with its first-order Taylor expansion along a random strain
 * direction \f$\Delta\varepsilon\f$:
 * \f[
 *   E(\alpha) = \|\sigma(\varepsilon_0+\alpha\Delta\varepsilon)-\sigma(\varepsilon_0)-\alpha\,\mathbb{D}\,\Delta\varepsilon\|
 *   = C\alpha^2 + O(\alpha^3), \qquad
 *   p = \frac{\log E(\alpha_2)-\log E(\alpha_1)}{\log\alpha_2-\log\alpha_1}\approx 2 .
 * \f]
 * The material is the clay of the Abaqus benchmark (Table 1): \f$M=1\f$, \f$\lambda=0.174\f$,
 * \f$\kappa=0.026\f$, \f$v_0=2.08\f$, porous elasticity with \f$\nu=0.3\f$. Every evaluation starts from an
 * isotropic converged state \f$\sigma_n=-p'_n I\f$, \f$p'_{c,n}=116.6\f$ kPa, \f$\varepsilon_n=0\f$ and applies
 * the total strain \f$\varepsilon\f$ in a single increment (mcc::ApplyStrain, i.e. the
 * ApplyStrainComputeSigma of TPZPlasticStepModifiedCamClay with the 6x6 tangent):
 *  - kind 0 (elastic): \f$p'_n=100\f$, strain components drawn in \f$[-0.002, 0.002]\f$;
 *  - kind 1 (subcritical, compaction): \f$p'_n=100\f$, components in \f$[-0.005, 0.005]\f$;
 *  - kind 2 (supercritical, dilation): \f$p'_n=50\f$, components in \f$[-0.02, 0.02]\f$.
 *
 * For each kind the state \f$\varepsilon_0\f$ (six engineering Voigt components XX, XY, XZ, YY, YZ, ZZ) is
 * drawn until the response is of that kind; \f$\sigma_0=\sigma(\varepsilon_0)\f$ and \f$\mathbb{D}_0\f$ are
 * computed together with the asymmetry \f$\|\mathbb{D}_0-\mathbb{D}_0^T\|/\|\mathbb{D}_0\|\f$. Then, for
 * \f$\mathbb{D}=\mathbb{D}_0\f$ and \f$\mathbb{D}=\mathbb{D}_0^T\f$, 300 pairs of amplitudes
 * \f$\alpha_1,\alpha_2\sim U(10^{-4},10^{-2})\f$ and directions (uniform in \f$[-1,1]^6\f$, scaled to
 * \f$\|\Delta\varepsilon\|=10^{-3}\f$) are drawn; pairs whose perturbed states change kind are discarded.
 * The slope of the least-squares line of \f$\log E\f$ against \f$\log\alpha\f$ over the 300 points
 * \f$(\alpha_1,E_1)\f$ and the median of the pairwise slopes \f$p\f$ are reported (Fig. 3).
 *
 * The random numbers are consumed in exactly the order of the function taylor() of gen_data.py. Two
 * generators are used:
 *  - TNumpyRandom: a transcription of numpy's default_rng(2026) (SeedSequence + PCG64 XSL-RR and
 *    Generator.uniform). It reproduces the draws of the Python script, so the states of Fig. 3, the
 *    asymmetries and the fitted slopes are those of the article;
 *  - TMersenneRandom: std::mt19937_64 with seed 2026, an independent sample that must be statistically
 *    equivalent (second order with \f$\mathbb{D}\f$, first order with \f$\mathbb{D}^T\f$ at plastic states).
 *
 * There is no finite element mesh in this example: the "material" is a single integration point, and
 * the methods follow the same sequence of the finite element examples (material set up, solution,
 * post-processing).
 */
class TaylorTest {
public:
    /** @brief Kind of response of the stress update (m_m_type of TPZPlasticStepModifiedCamClay) */
    enum EKind { EElastic = 0, ESubcritical = 1, ESupercritical = 2 };

    /** @name Random number generators */
    /** @{ */

    /** @brief Source of uniformly distributed random numbers */
    class TRandom {
    public:
        virtual ~TRandom() = default;
        /** @brief Uniform number in [lo, hi) */
        virtual REAL Uniform(REAL lo, REAL hi) = 0;
        /** @brief Short name, used in the names of the output files */
        virtual std::string Name() const = 0;
        /** @brief Description printed in the report */
        virtual std::string Description() const = 0;
    };

    /**
     * @brief Transcription of numpy.random.default_rng(seed).uniform(lo, hi): the seed is expanded by
     * numpy's SeedSequence (pool of four 32-bit words, hashmix/mix of O'Neill's seed_seq) into the
     * 128-bit state and increment of PCG64 (permuted congruential generator with XSL-RR output);
     * a double is \f$(x\gg 11)\,2^{-53}\f$ and Generator.uniform returns \f$lo+(hi-lo)\,u\f$.
     *
     * Implemented with portable 64-bit arithmetic (no compiler 128-bit integer). The first three raw
     * outputs of default_rng(2026) are 3300764713747675562, 11804314397344746687 and 8619580609625321962.
     */
    class TNumpyRandom : public TRandom {
    public:
        explicit TNumpyRandom(uint64_t seed);
        /** @brief Next 64-bit output of PCG64 (random_raw) */
        uint64_t Next64();
        /** @brief Compares the first raw outputs of seed 2026 with numpy's default_rng(2026).bit_generator.random_raw(3) */
        static bool SelfTest() {
            TNumpyRandom r(2026);
            const uint64_t a = r.Next64(), b = r.Next64(), c = r.Next64();
            return a == 3300764713747675562ULL && b == 11804314397344746687ULL && c == 8619580609625321962ULL;
        }
        /** @brief Next double in [0, 1) (next_double of numpy) */
        REAL NextDouble() { return REAL(Next64() >> 11) * (1.0 / 9007199254740992.0); }
        REAL Uniform(REAL lo, REAL hi) override { return lo + (hi - lo) * NextDouble(); }
        std::string Name() const override { return "pcg64"; }
        std::string Description() const override {
            return "numpy-compatible PCG64 stream, default_rng(2026): same draws as gen_data.py taylor()";
        }

    private:
        /** @brief 128-bit unsigned integer as two 64-bit words */
        struct TU128 {
            uint64_t fHi = 0, fLo = 0;
        };
        static TU128 Add(const TU128 &a, const TU128 &b);
        static TU128 Mul(const TU128 &a, const TU128 &b);
        static uint64_t MulHi(uint64_t a, uint64_t b);
        void Step() { fState = Add(Mul(fState, fMultiplier), fInc); }
        TU128 fState, fInc;
        /** @brief PCG_DEFAULT_MULTIPLIER_128 */
        const TU128 fMultiplier{0x2360ED051FC65DA4ULL, 0x4385DF649FCCF645ULL};
    };

    /** @brief std::mt19937_64 with std::uniform_real_distribution */
    class TMersenneRandom : public TRandom {
    public:
        explicit TMersenneRandom(uint64_t seed) : fGen(seed) {}
        REAL Uniform(REAL lo, REAL hi) override { return std::uniform_real_distribution<REAL>(lo, hi)(fGen); }
        std::string Name() const override { return "mt19937"; }
        std::string Description() const override {
            return "std::mt19937_64, seed 2026: independent sample (statistical comparison)";
        }

    private:
        std::mt19937_64 fGen;
    };
    /** @} */

    /** @brief Converged state from which every strain is applied (sigma_n = -p'_n I, pc_n, eps_n = 0) */
    struct TState {
        REAL fPn;     ///< mean effective stress p'_n of the isotropic state (kPa)
        REAL fPcn;    ///< preconsolidation pressure p'_c,n (kPa)
        REAL fRange;  ///< the components of eps_0 are drawn in [-range, range]
        std::string fName;
    };

    /** @brief Result of one panel of Fig. 3: a state and the operator D or D^T */
    struct TPanel {
        int fKind = 0;              ///< EKind of the state
        bool fTransposed = false;   ///< D^T instead of D
        int fDraws = 0;             ///< number of draws of eps_0 until the response was of the required kind
        TPZManVector<REAL, 6> fX0;  ///< state eps_0 (engineering Voigt components)
        TPZTensor<REAL> fSigma0;    ///< sigma(eps_0)
        TPZFNMatrix<36, REAL> fD0;  ///< consistent tangent at eps_0 (not transposed)
        REAL fPc = 0.;              ///< preconsolidation pressure after the update
        REAL fP = 0., fQ = 0.;      ///< p' and q of sigma_0 (kPa)
        REAL fAsym = 0.;            ///< ||D0 - D0^T|| / ||D0|| (Frobenius norms)
        int fRejected = 0;          ///< pairs discarded because a perturbed state changed kind
        std::vector<std::array<REAL, 4>> fPairs; ///< (alpha1, E1, alpha2, E2)
        std::vector<REAL> fSlopes;  ///< pairwise slopes p
        REAL fFitSlope = 0., fFitIntercept = 0.; ///< least-squares line log E1 = b0 + b1 log alpha1
        REAL fMedian = 0., fMinSlope = 0., fMaxSlope = 0.;
    };

    /**
     * @brief Values reported by the article (Fig. 3 and Sect. 4.5) and by the Python script
     * (gen_data.py taylor(), numpy default_rng(2026)); a negative value means "not reported"
     */
    struct TReference {
        REAL fArtP, fArtQ, fArtAsym, fArtFit, fArtMedian;
        REAL fPyP, fPyQ, fPyAsym, fPyFit, fPyMedian;
    };

    /** @name Parameters (Table 1, Taylor test) */
    /** @{ */
    REAL fM = 1.0, fLambda = 0.174, fKappa = 0.026, fV0 = 2.08, fNu = 0.3;
    int fNPairs = 300;                          ///< pairs of amplitudes per panel
    REAL fAlphaMin = 1.e-4, fAlphaMax = 1.e-2;  ///< range of the amplitudes
    REAL fDirectionNorm = 1.e-3;                ///< norm of the strain direction
    uint64_t fSeed = 2026;
    /** @} */

    /** @brief Constitutive model: MCC, porous elasticity, shear modulus from the Poisson ratio */
    mcc::TPlastic CreateMaterial() const;

    /** @brief Converged state of a kind of test (Sect. 4.5) */
    TState State(int kind) const;

    /**
     * @brief Stress update from the converged state of a kind for the total strain x (one increment)
     * @param model constitutive model (not changed)
     * @param kind which converged state (EKind)
     * @param x total strain, engineering Voigt components XX, XY, XZ, YY, YZ, ZZ
     * @param[out] sigma updated stress
     * @param[out] D consistent tangent d sigma / d eps (6x6)
     * @param[out] pc updated preconsolidation pressure
     * @param[out] type kind of the response (0 elastic, 1 subcritical, 2 supercritical)
     * @return false if the local projection failed
     */
    bool Response(const mcc::TPlastic &model, int kind, const TPZVec<REAL> &x, TPZTensor<REAL> &sigma,
                  TPZFMatrix<REAL> &D, REAL &pc, int &type) const;

    /** @brief Draws eps_0 until the response is of the given kind; computes sigma_0, D_0, p', q and the asymmetry */
    TPanel DrawState(const mcc::TPlastic &model, int kind, TRandom &rng) const;

    /** @brief Draws the pairs of amplitudes and directions of one panel and computes the slopes */
    void Perturb(const mcc::TPlastic &model, TPanel &panel, TRandom &rng) const;

    /** @brief The complete test with one generator: kinds 0, 1, 2, each with D and D^T (order of taylor()) */
    std::vector<TPanel> Run(TRandom &rng) const;

    /** @brief Writes one CSV file per panel with the points of Fig. 3 and a summary file with the fits */
    void PostProcess(const std::vector<TPanel> &panels, const std::string &prefix) const;

    /** @brief Prints the results of a run next to the values of the article and of the Python script */
    void Print(const std::vector<TPanel> &panels, const TRandom &rng, bool comparePython) const;

    /** @brief Runs the test with the numpy-compatible generator and with std::mt19937_64 */
    void RunAll();

    /** @brief Reference values of a panel */
    static TReference Reference(int kind, bool transposed);

    /** @brief Name of a kind of state */
    static std::string KindName(int kind);

private:
    /** @brief ||sigma(x0 + alpha dx) - sigma0 - alpha D dx||; returns false if the update failed */
    bool TaylorError(const mcc::TPlastic &model, const TPanel &panel, const TPZFMatrix<REAL> &D, REAL alpha,
                     const TPZVec<REAL> &dx, REAL &error, int &type) const;
};

// ------------------------------------------------------------------------------------------------ TNumpyRandom

inline TaylorTest::TNumpyRandom::TNumpyRandom(uint64_t seed) {
    // SeedSequence(seed): entropy as 32-bit words (least significant first)
    std::vector<uint32_t> entropy;
    do {
        entropy.push_back(uint32_t(seed & 0xFFFFFFFFULL));
        seed >>= 32;
    } while (seed);
    const uint32_t INIT_A = 0x43b0d7e5U, MULT_A = 0x931e8875U, INIT_B = 0x8b51f9ddU, MULT_B = 0x58f38dedU;
    const uint32_t MIX_MULT_L = 0xca01f9ddU, MIX_MULT_R = 0x4973f715U;
    uint32_t hashConst = INIT_A;
    auto hashmix = [&hashConst, MULT_A](uint32_t value) {
        value ^= hashConst;
        hashConst *= MULT_A;
        value *= hashConst;
        value ^= value >> 16;
        return value;
    };
    auto mix = [MIX_MULT_L, MIX_MULT_R](uint32_t x, uint32_t y) {
        uint32_t r = MIX_MULT_L * x - MIX_MULT_R * y;
        r ^= r >> 16;
        return r;
    };
    // mix_entropy with a pool of four words
    const size_t npool = 4;
    uint32_t pool[npool];
    for (size_t i = 0; i < npool; ++i) pool[i] = hashmix(i < entropy.size() ? entropy[i] : 0U);
    for (size_t isrc = 0; isrc < npool; ++isrc)
        for (size_t idst = 0; idst < npool; ++idst)
            if (isrc != idst) pool[idst] = mix(pool[idst], hashmix(pool[isrc]));
    for (size_t isrc = npool; isrc < entropy.size(); ++isrc)
        for (size_t idst = 0; idst < npool; ++idst) pool[idst] = mix(pool[idst], hashmix(entropy[isrc]));
    // generate_state(4, uint64): eight 32-bit words viewed as four little-endian 64-bit words
    uint32_t words[8];
    uint32_t hb = INIT_B;
    for (int i = 0; i < 8; ++i) {
        uint32_t v = pool[i % npool];
        v ^= hb;
        hb *= MULT_B;
        v *= hb;
        v ^= v >> 16;
        words[i] = v;
    }
    uint64_t s[4];
    for (int k = 0; k < 4; ++k) s[k] = uint64_t(words[2 * k]) | (uint64_t(words[2 * k + 1]) << 32);
    // pcg64_set_seed: initstate = (s0, s1), initseq = (s2, s3) as (high, low) words
    const TU128 initstate{s[0], s[1]};
    TU128 initseq{s[2], s[3]};
    // pcg_setseq_128_srandom_r
    fState = TU128{0, 0};
    fInc.fHi = (initseq.fHi << 1) | (initseq.fLo >> 63);
    fInc.fLo = (initseq.fLo << 1) | 1ULL;
    Step();
    fState = Add(fState, initstate);
    Step();
}

inline uint64_t TaylorTest::TNumpyRandom::MulHi(uint64_t a, uint64_t b) {
    const uint64_t a0 = a & 0xFFFFFFFFULL, a1 = a >> 32, b0 = b & 0xFFFFFFFFULL, b1 = b >> 32;
    const uint64_t p00 = a0 * b0, p01 = a0 * b1, p10 = a1 * b0, p11 = a1 * b1;
    const uint64_t mid = (p00 >> 32) + (p01 & 0xFFFFFFFFULL) + (p10 & 0xFFFFFFFFULL);
    return p11 + (p01 >> 32) + (p10 >> 32) + (mid >> 32);
}

inline TaylorTest::TNumpyRandom::TU128 TaylorTest::TNumpyRandom::Add(const TU128 &a, const TU128 &b) {
    TU128 r;
    r.fLo = a.fLo + b.fLo;
    r.fHi = a.fHi + b.fHi + (r.fLo < a.fLo ? 1ULL : 0ULL);
    return r;
}

inline TaylorTest::TNumpyRandom::TU128 TaylorTest::TNumpyRandom::Mul(const TU128 &a, const TU128 &b) {
    TU128 r;
    r.fLo = a.fLo * b.fLo;
    r.fHi = MulHi(a.fLo, b.fLo) + a.fHi * b.fLo + a.fLo * b.fHi;
    return r;
}

inline uint64_t TaylorTest::TNumpyRandom::Next64() {
    // pcg_setseq_128_xsl_rr_64_random_r: advance, then XSL-RR output of the new state
    Step();
    const uint64_t x = fState.fHi ^ fState.fLo;
    const unsigned rot = unsigned(fState.fHi >> 58);
    return (x >> rot) | (x << ((64U - rot) & 63U));
}

// ------------------------------------------------------------------------------------------------ TaylorTest

inline mcc::TPlastic TaylorTest::CreateMaterial() const {
    mcc::TPlastic model;
    model.SetModifiedCamClay(fM, fLambda, fKappa);
    model.SetPorousElasticity();
    model.SetPoissonRatio(fNu);
    model.SetDefaultSpecificVolume(fV0);
    return model;
}

inline TaylorTest::TState TaylorTest::State(int kind) const {
    switch (kind) {
        case EElastic: return {100., 116.6, 0.002, "elastic"};
        case ESubcritical: return {100., 116.6, 0.005, "subcritical"};
        default: return {50., 116.6, 0.02, "supercritical"};
    }
}

inline std::string TaylorTest::KindName(int kind) {
    const char *names[3] = {"elastic", "subcritical", "supercritical"};
    return names[kind];
}

inline bool TaylorTest::Response(const mcc::TPlastic &model, int kind, const TPZVec<REAL> &x, TPZTensor<REAL> &sigma,
                                 TPZFMatrix<REAL> &D, REAL &pc, int &type) const {
    const TState st = State(kind);
    const mcc::TPointState sn(mcc::IsotropicTensor(-st.fPn), st.fPcn, fV0);
    TPZTensor<REAL> eps, epsn;
    for (int i = 0; i < 6; ++i) eps[i] = x[i];
    return mcc::ApplyStrain(model, epsn, sn, eps, sigma, D, pc, type);
}

inline TaylorTest::TPanel TaylorTest::DrawState(const mcc::TPlastic &model, int kind, TRandom &rng) const {
    const TState st = State(kind);
    TPanel panel;
    panel.fKind = kind;
    panel.fX0.Resize(6, 0.);
    REAL pc = 0.;
    int type = -1;
    while (true) {
        for (int i = 0; i < 6; ++i) panel.fX0[i] = rng.Uniform(-st.fRange, st.fRange);
        panel.fDraws++;
        if (Response(model, kind, panel.fX0, panel.fSigma0, panel.fD0, pc, type) && type == kind) break;
    }
    panel.fPc = pc;
    panel.fP = mcc::MeanEffectiveStress(panel.fSigma0);
    panel.fQ = mcc::DeviatoricStress(panel.fSigma0);
    REAL nd = 0., na = 0.;
    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j) {
            const REAL d = panel.fD0(i, j), a = panel.fD0(i, j) - panel.fD0(j, i);
            nd += d * d;
            na += a * a;
        }
    panel.fAsym = std::sqrt(na) / std::sqrt(nd);
    return panel;
}

inline bool TaylorTest::TaylorError(const mcc::TPlastic &model, const TPanel &panel, const TPZFMatrix<REAL> &D,
                                    REAL alpha, const TPZVec<REAL> &dx, REAL &error, int &type) const {
    TPZManVector<REAL, 6> x(6);
    for (int i = 0; i < 6; ++i) x[i] = panel.fX0[i] + alpha * dx[i];
    TPZTensor<REAL> sigma;
    TPZFNMatrix<36, REAL> Dx(6, 6, 0.);
    REAL pc;
    if (!Response(model, panel.fKind, x, sigma, Dx, pc, type)) return false;
    REAL e2 = 0.;
    for (int i = 0; i < 6; ++i) {
        REAL lin = 0.;
        for (int j = 0; j < 6; ++j) lin += (alpha * D.GetVal(i, j)) * dx[j];
        const REAL r = sigma[i] - panel.fSigma0[i] - lin;
        e2 += r * r;
    }
    error = std::sqrt(e2);
    return true;
}

inline void TaylorTest::Perturb(const mcc::TPlastic &model, TPanel &panel, TRandom &rng) const {
    TPZFNMatrix<36, REAL> D(6, 6, 0.);
    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j) D(i, j) = panel.fTransposed ? panel.fD0(j, i) : panel.fD0(i, j);
    panel.fPairs.clear();
    panel.fSlopes.clear();
    panel.fRejected = 0;
    TPZManVector<REAL, 6> dx(6);
    while ((int)panel.fPairs.size() < fNPairs) {
        const REAL a1 = rng.Uniform(fAlphaMin, fAlphaMax);
        const REAL a2 = rng.Uniform(fAlphaMin, fAlphaMax);
        REAL norm = 0.;
        for (int i = 0; i < 6; ++i) {
            dx[i] = rng.Uniform(-1., 1.);
            norm += dx[i] * dx[i];
        }
        norm = std::sqrt(norm);
        for (int i = 0; i < 6; ++i) dx[i] = fDirectionNorm * dx[i] / norm;
        REAL e1 = 0., e2 = 0.;
        int t1 = -1, t2 = -1;
        const bool ok1 = TaylorError(model, panel, D, a1, dx, e1, t1);
        const bool ok2 = TaylorError(model, panel, D, a2, dx, e2, t2);
        if (!ok1 || !ok2 || t1 != panel.fKind || t2 != panel.fKind) {
            panel.fRejected++;
            continue;
        }
        panel.fPairs.push_back({a1, e1, a2, e2});
        panel.fSlopes.push_back(std::log(e2 / e1) / std::log(a2 / a1));
    }
    // least-squares line of log E1 against log alpha1 (np.polyfit of degree 1)
    const REAL n = panel.fPairs.size();
    REAL xm = 0., ym = 0.;
    for (auto &p : panel.fPairs) {
        xm += std::log(p[0]);
        ym += std::log(p[1]);
    }
    xm /= n;
    ym /= n;
    REAL sxy = 0., sxx = 0.;
    for (auto &p : panel.fPairs) {
        const REAL dxl = std::log(p[0]) - xm, dyl = std::log(p[1]) - ym;
        sxy += dxl * dyl;
        sxx += dxl * dxl;
    }
    panel.fFitSlope = sxy / sxx;
    panel.fFitIntercept = ym - panel.fFitSlope * xm;
    // median (mean of the two central values for an even number), extreme values
    std::vector<REAL> sorted(panel.fSlopes);
    std::sort(sorted.begin(), sorted.end());
    const size_t m = sorted.size();
    panel.fMedian = (m % 2) ? sorted[m / 2] : 0.5 * (sorted[m / 2 - 1] + sorted[m / 2]);
    panel.fMinSlope = sorted.front();
    panel.fMaxSlope = sorted.back();
}

inline std::vector<TaylorTest::TPanel> TaylorTest::Run(TRandom &rng) const {
    const mcc::TPlastic model = CreateMaterial();
    std::vector<TPanel> panels;
    for (int kind : {EElastic, ESubcritical, ESupercritical}) {
        const TPanel state = DrawState(model, kind, rng);
        for (bool transposed : {false, true}) {
            TPanel panel(state);
            panel.fTransposed = transposed;
            Perturb(model, panel, rng);
            panels.push_back(panel);
        }
    }
    return panels;
}

inline void TaylorTest::PostProcess(const std::vector<TPanel> &panels, const std::string &prefix) const {
    std::vector<std::vector<REAL>> summary;
    for (auto &p : panels) {
        std::vector<std::vector<REAL>> rows;
        for (size_t k = 0; k < p.fPairs.size(); ++k) {
            const auto &q = p.fPairs[k];
            const REAL la = std::log(q[0]);
            rows.push_back({q[0], q[1], la, std::log(q[1]), p.fFitIntercept + p.fFitSlope * la, q[2], q[3],
                            p.fSlopes[k]});
        }
        const std::string file = prefix + "_" + KindName(p.fKind) + (p.fTransposed ? "_DT" : "_D") + ".csv";
        mcc::WriteCSV(file, {"alpha1", "E1", "log_alpha1", "log_E1", "fit_log_E1", "alpha2", "E2", "pair_slope"}, rows);
        std::vector<REAL> row = {REAL(p.fKind), REAL(p.fTransposed), p.fP, p.fQ, p.fPc, p.fAsym, p.fFitSlope,
                                 p.fFitIntercept, p.fMedian, p.fMinSlope, p.fMaxSlope, REAL(p.fRejected),
                                 REAL(p.fDraws)};
        for (int i = 0; i < 6; ++i) row.push_back(p.fX0[i]);
        for (int i = 0; i < 6; ++i) row.push_back(p.fSigma0[i]);
        summary.push_back(row);
    }
    mcc::WriteCSV(prefix + "_summary.csv",
                  {"kind", "transposed", "p_eff", "q", "pc", "asym", "fit_slope", "fit_intercept", "median_slope",
                   "min_slope", "max_slope", "rejected", "draws", "x0_xx", "x0_xy", "x0_xz", "x0_yy", "x0_yz",
                   "x0_zz", "sig0_xx", "sig0_xy", "sig0_xz", "sig0_yy", "sig0_yz", "sig0_zz"},
                  summary);
}

inline TaylorTest::TReference TaylorTest::Reference(int kind, bool transposed) {
    // article: Fig. 3 and Sect. 4.5; Python: gen_data.py taylor() with numpy default_rng(2026)
    switch (kind) {
        case EElastic:
            return transposed ? TReference{-1, -1, -1, -1, 2.000, 105.254901154231, 16.8355063140829, 0.0,
                                           2.284155357681005, 1.9999917921427404}
                              : TReference{-1, -1, -1, -1, 2.000, 105.254901154231, 16.8355063140829, 0.0,
                                           2.043055635448701, 1.9999943584092064};
        case ESubcritical:
            return transposed ? TReference{103.4, 40.7, 0.016, 1.039, 1.000, 103.411561037478, 40.7038106469357,
                                           0.015835622133207647, 1.0393804933210131, 1.0002706710798948}
                              : TReference{103.4, 40.7, 0.016, 1.996, 2.000, 103.411561037478, 40.7038106469357,
                                           0.015835622133207647, 1.996015341520114, 1.9999874160215851};
        default:
            return transposed ? TReference{23.7, 45.6, 0.033, 1.013, 1.000, 23.6707002558066, 45.6311078924772,
                                           0.032875208106802306, 1.0129987302952381, 0.9998092644988404}
                              : TReference{23.7, 45.6, 0.033, 1.976, 2.000, 23.6707002558066, 45.6311078924772,
                                           0.032875208106802306, 1.976462982150101, 1.9999975521482396};
    }
}

inline void TaylorTest::Print(const std::vector<TPanel> &panels, const TRandom &rng, bool comparePython) const {
    auto ref = [](REAL v, int prec) {
        std::ostringstream s;
        if (v < 0.) s << "-";
        else s << std::fixed << std::setprecision(prec) << v;
        return s.str();
    };
    std::cout << "\n" << rng.Description() << "\n";
    std::cout << "  values: this work [article]" << (comparePython ? " {Python}" : "")
              << "; asym = ||D - D^T||/||D||; slopes of log E against log alpha\n";
    std::cout << std::left << std::setw(17) << "  state" << std::setw(5) << "op" << std::setw(17) << "p' (kPa)"
              << std::setw(17) << "q (kPa)" << std::setw(28) << "asym" << std::setw(32) << "fitted slope"
              << "median slope" << "\n";
    for (auto &p : panels) {
        const TReference r = Reference(p.fKind, p.fTransposed);
        std::ostringstream cp, cq, ca, cf, cm;
        cp << std::fixed << std::setprecision(3) << p.fP << " [" << ref(r.fArtP, 1) << "]";
        cq << std::fixed << std::setprecision(3) << p.fQ << " [" << ref(r.fArtQ, 1) << "]";
        ca << std::fixed << std::setprecision(5) << p.fAsym << " [" << ref(r.fArtAsym, 3) << "]";
        cf << std::fixed << std::setprecision(6) << p.fFitSlope << " [" << ref(r.fArtFit, 3) << "]";
        cm << std::fixed << std::setprecision(6) << p.fMedian << " [" << ref(r.fArtMedian, 3) << "]";
        if (comparePython) {
            ca << " {" << ref(r.fPyAsym, 5) << "}";
            cf << " {" << ref(r.fPyFit, 6) << "}";
            cm << " {" << ref(r.fPyMedian, 6) << "}";
        }
        std::cout << "  " << std::setw(15) << KindName(p.fKind) << std::setw(5) << (p.fTransposed ? "D^T" : "D")
                  << std::setw(17) << cp.str() << std::setw(17) << cq.str() << std::setw(28) << ca.str()
                  << std::setw(32) << cf.str() << cm.str() << "\n";
    }
    std::cout << std::right;
    if (comparePython) {
        std::cout << "  Python states (gen_data.py): p' = " << std::setprecision(15) << Reference(0, false).fPyP << ", "
                  << Reference(1, false).fPyP << ", " << Reference(2, false).fPyP
                  << "; q = " << Reference(0, false).fPyQ << ", " << Reference(1, false).fPyQ << ", "
                  << Reference(2, false).fPyQ << "\n";
        std::cout << "  this work:                   p' = ";
        for (int k = 0; k < 3; ++k) std::cout << panels[2 * k].fP << (k < 2 ? ", " : "; q = ");
        for (int k = 0; k < 3; ++k) std::cout << panels[2 * k].fQ << (k < 2 ? ", " : "\n");
    }
    std::cout << "  pairs per panel = " << fNPairs << ", rejected (state changed kind):";
    for (auto &p : panels) std::cout << " " << p.fRejected;
    std::cout << "; pairwise slopes, min/max:";
    std::cout << std::fixed << std::setprecision(3);
    for (auto &p : panels) std::cout << " " << p.fMinSlope << "/" << p.fMaxSlope;
    std::cout << std::defaultfloat << std::setprecision(6) << "\n";
}

inline void TaylorTest::RunAll() {
    std::cout << "Taylor test of the consistent tangent (Sect. 4.5, Fig. 3): M = " << fM << ", lambda = " << fLambda
              << ", kappa = " << fKappa << ", v0 = " << fV0 << ", porous elasticity, nu = " << fNu << "\n";
    std::cout << "states: sigma_n = -p'_n I, p'_c,n = 116.6 kPa; elastic p'_n = 100 (range 0.002), subcritical "
                 "p'_n = 100 (0.005), supercritical p'_n = 50 (0.02)\n";
    // 1. same random numbers as the Python script: reproduces the states and slopes of Fig. 3
    std::cout << "check of the PCG64 transcription against numpy default_rng(2026).bit_generator.random_raw(3): "
              << (TNumpyRandom::SelfTest() ? "ok" : "FAILED") << "\n";
    TNumpyRandom pcg(fSeed);
    const std::vector<TPanel> fig3 = Run(pcg);
    Print(fig3, pcg, true);
    PostProcess(fig3, "taylor_" + pcg.Name());
    // 2. independent sample with the C++ standard generator: statistical comparison
    TMersenneRandom mt(fSeed);
    const std::vector<TPanel> sample = Run(mt);
    Print(sample, mt, false);
    PostProcess(sample, "taylor_" + mt.Name());
}
