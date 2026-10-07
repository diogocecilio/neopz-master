/**
 * @file TaylorTest.h
 * @brief Sect. 4.5 and Fig. 3 of the article: Taylor test of the consistent tangent operator of the
 * Modified Cam-Clay return mapping in rotated Haigh-Westergaard space (TPZPlasticStepModifiedCamClay), and the
 * Taylor slopes of the alternative tangent operators of Sect. 6.7 (last column of Table 10).
 */
#pragma once

#include "MCCPaperTools.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <sstream>
#include <string>
#include <vector>

/**
 * @ingroup mccpaper
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
 * The pairwise slope \f$p\f$ uses one direction for both amplitudes, so the constant \f$C(\Delta\varepsilon)\f$
 * cancels; the least-squares line mixes 300 directions and its slope carries the scatter of
 * \f$\log C(\Delta\varepsilon)\f$. In elastic steps the shear modulus is that of the converged state and only
 * the porous volumetric law is nonlinear, so \f$C\propto(\mathrm{tr}\,\Delta\varepsilon)^2\f$: the scatter of
 * \f$\log C\f$ is large and the fitted slope is not a measure of the order there (the article reports only
 * the median, 2.000, for elastic steps).
 *
 * The random numbers are consumed in exactly the order of the function taylor() of gen_data.py. Two
 * generators are used:
 *  - TNumpyRandom: a transcription of numpy's default_rng(2026) (SeedSequence + PCG64 XSL-RR and
 *    Generator.uniform). It reproduces the draws of the Python script, so the states of Fig. 3, the
 *    asymmetries and the fitted slopes are those of the article;
 *  - TMersenneRandom: std::mt19937_64 with seed 2026, an independent sample that must be statistically
 *    equivalent (second order with \f$\mathbb{D}\f$, first order with \f$\mathbb{D}^T\f$ at plastic states).
 *
 * Taylor slopes of the operators of Table 10 (Sect. 6.7; tangentes() of gen_data.py, OperatorSlopes): at the
 * subcritical and supercritical states of Fig. 3 (the same \f$\varepsilon_0\f$), the test is repeated with the
 * consistent tangent \f$\mathbb{D}\f$, its transpose \f$\mathbb{D}^T\f$, its symmetric part
 * \f$(\mathbb{D}+\mathbb{D}^T)/2\f$ and the continuum operator (TPZPlasticStepModifiedCamClay::ContinuumTangent at
 * \f$\sigma_0\f$ and the updated \f$p_c\f$), with 300 pairs per operator drawn from TNumpyRandom(7), i.e.
 * numpy's default_rng(7), in the order of the Python script (states 1, 2; operators D, DT, sym, cont). Each
 * operator is also obtained from the library with TPZPlasticStepModifiedCamClay::SetTangentMode (EConsistentTangent,
 * ETransposedTangent, ESymmetricTangent, EContinuumTangent), the operators of the global iterations of Table 10,
 * and must be identical to the one tested (column diff_library_mode of taylor_operators_summary.csv).
 *
 * There is no finite element mesh in this example: the "material" is a single integration point, and
 * the methods follow the same sequence as the finite element examples (material set-up, solution,
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
        /** @brief Virtual destructor of the interface */
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
        /** @brief Seeds the generator as numpy.random.default_rng(seed) (SeedSequence, then pcg64_set_seed) */
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
        /** @brief Generator.uniform of numpy: \f$lo+(hi-lo)\,u\f$ */
        REAL Uniform(REAL lo, REAL hi) override { return lo + (hi - lo) * NextDouble(); }
        std::string Name() const override { return "pcg64"; }
        std::string Description() const override {
            return "numpy-compatible PCG64 stream, default_rng(2026): same draws as gen_data.py taylor()";
        }

    private:
        /** @brief 128-bit unsigned integer as two 64-bit words */
        struct TU128 {
            uint64_t fHi = 0; ///< most significant word
            uint64_t fLo = 0; ///< least significant word
        };
        /** @brief Sum modulo \f$2^{128}\f$ */
        static TU128 Add(const TU128 &a, const TU128 &b);
        /** @brief Product modulo \f$2^{128}\f$ */
        static TU128 Mul(const TU128 &a, const TU128 &b);
        /** @brief Most significant 64 bits of the 128-bit product of two 64-bit words */
        static uint64_t MulHi(uint64_t a, uint64_t b);
        /** @brief Linear congruential step: state = state * multiplier + increment */
        void Step() { fState = Add(Mul(fState, fMultiplier), fInc); }
        /** @brief State of the congruential generator */
        TU128 fState;
        /** @brief Increment (odd) of the congruential generator */
        TU128 fInc;
        /** @brief PCG_DEFAULT_MULTIPLIER_128 */
        const TU128 fMultiplier{0x2360ED051FC65DA4ULL, 0x4385DF649FCCF645ULL};
    };

    /**
     * @brief std::mt19937_64 with the 53-bit conversion \f$u=(x\gg 11)\,2^{-53}\f$ and \f$lo+(hi-lo)\,u\f$:
     * the engine and the conversion are fully specified, so the sample is the same with any standard
     * library (std::uniform_real_distribution is implementation defined)
     */
    class TMersenneRandom : public TRandom {
    public:
        /** @brief Seeds std::mt19937_64 */
        explicit TMersenneRandom(uint64_t seed) : fGen(seed) {}
        REAL Uniform(REAL lo, REAL hi) override {
            return lo + (hi - lo) * (REAL(fGen() >> 11) * (1.0 / 9007199254740992.0));
        }
        std::string Name() const override { return "mt19937"; }
        std::string Description() const override {
            return "std::mt19937_64, seed 2026: independent sample (statistical comparison, orders only)";
        }

    private:
        /** @brief The 64-bit Mersenne Twister engine */
        std::mt19937_64 fGen;
    };
    /** @} */

    /** @brief Converged state from which every strain is applied (sigma_n = -p'_n I, pc_n, eps_n = 0) */
    struct TState {
        REAL fPn;     ///< mean effective stress p'_n of the isotropic state (kPa)
        REAL fPcn;    ///< preconsolidation pressure p'_c,n (kPa)
        REAL fRange;  ///< the components of eps_0 are drawn in [-range, range]
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
        REAL fP = 0.;               ///< p' of sigma_0 (kPa)
        REAL fQ = 0.;               ///< q of sigma_0 (kPa)
        REAL fAsym = 0.;            ///< ||D0 - D0^T|| / ||D0|| (Frobenius norms)
        int fRejected = 0;          ///< pairs discarded because a perturbed state changed kind
        std::vector<std::array<REAL, 4>> fPairs; ///< (alpha1, E1, alpha2, E2)
        std::vector<REAL> fSlopes;  ///< pairwise slopes p
        REAL fFitSlope = 0.;        ///< slope b1 of the least-squares line log E1 = b0 + b1 log alpha1
        REAL fFitIntercept = 0.;    ///< intercept b0 of the least-squares line
        REAL fMedian = 0.;          ///< median of the pairwise slopes
        REAL fMinSlope = 0.;        ///< smallest pairwise slope
        REAL fMaxSlope = 0.;        ///< largest pairwise slope
    };

    /**
     * @brief Values reported by the article (Fig. 3 and Sect. 4.5) and by the Python script
     * (gen_data.py taylor(), numpy default_rng(2026)); a negative value means "not reported"
     */
    struct TReference {
        REAL fArtP;      ///< article: p' of the state (kPa)
        REAL fArtQ;      ///< article: q of the state (kPa)
        REAL fArtAsym;   ///< article: asymmetry ||D - D^T||/||D||
        REAL fArtFit;    ///< article: slope of the least-squares line
        REAL fArtMedian; ///< article: median of the pairwise slopes
        REAL fPyP;       ///< Python: p' of the state (kPa)
        REAL fPyQ;       ///< Python: q of the state (kPa)
        REAL fPyAsym;    ///< Python: asymmetry
        REAL fPyFit;     ///< Python: slope of the least-squares line (np.polyfit)
        REAL fPyMedian;  ///< Python: median of the pairwise slopes
    };

    /** @name Parameters (Table 1, Taylor test) */
    /** @{ */
    REAL fM = 1.0;          ///< slope of the critical state line
    REAL fLambda = 0.174;   ///< slope of the normal compression line
    REAL fKappa = 0.026;    ///< slope of the swelling line
    REAL fV0 = 2.08;        ///< specific volume
    REAL fNu = 0.3;         ///< Poisson ratio (shear modulus from the porous bulk modulus)
    int fNPairs = 300;      ///< pairs of amplitudes per panel
    REAL fAlphaMin = 1.e-4; ///< smallest amplitude
    REAL fAlphaMax = 1.e-2; ///< largest amplitude
    REAL fDirectionNorm = 1.e-3; ///< norm of the strain direction
    uint64_t fSeed = 2026;  ///< seed of both generators (as in gen_data.py)
    uint64_t fOperatorSeed = 7; ///< seed of the draws of the comparison of the operators (tangentes() of gen_data.py)
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

    /**
     * @brief Prints the results of a run next to the reference values
     * @param panels results of Run
     * @param rng generator used (for the description)
     * @param sameDraws true for the numpy-compatible stream: the states are those of the article, and the
     * article and Python values are printed; false for an independent sample: only the slopes are compared
     * with the article
     */
    void Print(const std::vector<TPanel> &panels, const TRandom &rng, bool sameDraws) const;

    /** @brief Draws the pairs of one panel with the operator D (6x6) and computes the slopes and the fit */
    void PerturbWith(const mcc::TPlastic &model, TPanel &panel, const TPZFMatrix<REAL> &D, TRandom &rng) const;

    /** @brief Taylor test of one tangent operator at a state of Fig. 3 (last column of Table 10) */
    struct TOperatorPanel {
        std::string fOperator; ///< D, DT, sym or cont (names of gen_data.py)
        TPanel fPanel;         ///< state, pairs, slopes and fit
        TPZFNMatrix<36, REAL> fMatrix; ///< the operator
        REAL fLibraryDiff = -1.; ///< largest |difference| from the tangent of the library mode (SetTangentMode)
    };

    /**
     * @brief Taylor slopes of the consistent tangent D, its transpose, its symmetric part and the continuum
     * operator at the subcritical and supercritical states of Fig. 3 (tangentes() of gen_data.py)
     * @param fig3 panels of Run with the numpy-compatible stream of seed 2026 (the states eps_0 of Fig. 3)
     * @return eight panels, in the order (kind 1: D, DT, sym, cont), (kind 2: D, DT, sym, cont); the draws come
     * from TNumpyRandom(fOperatorSeed)
     */
    std::vector<TOperatorPanel> OperatorSlopes(const std::vector<TPanel> &fig3) const;

    /** @brief Writes taylor_operators_summary.csv and one file of points per operator and state */
    void PostProcessOperators(const std::vector<TOperatorPanel> &ops) const;

    /** @brief Prints the slopes of the operators with the values of the Python code and of Table 10 */
    void PrintOperators(const std::vector<TOperatorPanel> &ops) const;

    /** @brief Python values (data_tangentes.pkl) of the fitted and median slopes of an operator at a state */
    static void OperatorReference(int kind, const std::string &op, REAL &fit, REAL &median);

    /** @brief Runs the test with the numpy-compatible generator and with std::mt19937_64, then the operators */
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
    const TU128 initseq{s[2], s[3]};
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
        case EElastic: return {100., 116.6, 0.002};
        case ESubcritical: return {100., 116.6, 0.005};
        default: return {50., 116.6, 0.02};
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
    TPZFNMatrix<36, REAL> Dpert(6, 6, 0.); // tangent at the perturbed state (not used)
    REAL pc;
    if (!Response(model, panel.fKind, x, sigma, Dpert, pc, type)) return false;
    // (sigma - sigma0) - (alpha D) dx, as r[0] - f0 - a * D @ dx in gen_data.py
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
    PerturbWith(model, panel, D, rng);
}

inline void TaylorTest::PerturbWith(const mcc::TPlastic &model, TPanel &panel, const TPZFMatrix<REAL> &D,
                                    TRandom &rng) const {
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
    // median (mean of the two central values for an even number, as np.median), extreme values
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

inline void TaylorTest::Print(const std::vector<TPanel> &panels, const TRandom &rng, bool sameDraws) const {
    // reference value with a given number of decimals ("-" if not reported); scale 100 for percentages
    auto ref = [](REAL v, int prec, REAL scale = 1., const char *unit = "") {
        std::ostringstream s;
        if (v < 0.) s << "-";
        else s << std::fixed << std::setprecision(prec) << scale * v << unit;
        return s.str();
    };
    std::cout << "\n" << rng.Description() << "\n";
    std::cout << "  values: this work [article]" << (sameDraws ? " {Python}" : "")
              << "; asym = ||D - D^T||/||D||; slopes of log E against log alpha\n";
    if (!sameDraws)
        std::cout << "  (the states differ from those of Fig. 3: only the orders 2 (D) and 1 (D^T) are comparable)\n";
    std::cout << std::left << std::setw(17) << "  state" << std::setw(5) << "op" << std::setw(17) << "p' (kPa)"
              << std::setw(17) << "q (kPa)" << std::setw(26) << "asym" << std::setw(32) << "fitted slope"
              << "median slope" << "\n";
    for (auto &p : panels) {
        const TReference r = Reference(p.fKind, p.fTransposed);
        std::ostringstream cp, cq, ca, cf, cm;
        cp << std::fixed << std::setprecision(3) << p.fP;
        cq << std::fixed << std::setprecision(3) << p.fQ;
        ca << std::fixed << std::setprecision(3) << 100. * p.fAsym << "%";
        cf << std::fixed << std::setprecision(6) << p.fFitSlope << " [" << ref(r.fArtFit, 3) << "]";
        cm << std::fixed << std::setprecision(6) << p.fMedian << " [" << ref(r.fArtMedian, 3) << "]";
        if (sameDraws) {
            cp << " [" << ref(r.fArtP, 1) << "]";
            cq << " [" << ref(r.fArtQ, 1) << "]";
            ca << " [" << ref(r.fArtAsym, 1, 100., "%") << "] {" << ref(r.fPyAsym, 3, 100., "%") << "}";
            cf << " {" << ref(r.fPyFit, 6) << "}";
            cm << " {" << ref(r.fPyMedian, 6) << "}";
        }
        std::cout << "  " << std::setw(15) << KindName(p.fKind) << std::setw(5) << (p.fTransposed ? "D^T" : "D")
                  << std::setw(17) << cp.str() << std::setw(17) << cq.str() << std::setw(26) << ca.str()
                  << std::setw(32) << cf.str() << cm.str() << "\n";
    }
    std::cout << std::right;
    if (sameDraws) {
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

inline void TaylorTest::OperatorReference(int kind, const std::string &op, REAL &fit, REAL &median) {
    // gen_data.py tangentes(), numpy default_rng(7): out[('taylor', kind, op)] = (fitted slope, median slope)
    struct TRef {
        int fKind;
        const char *fOp;
        REAL fFit, fMedian;
    };
    static const TRef refs[] = {{1, "D", 2.0139100258796616, 2.0000096402788152},
                                {1, "DT", 1.0158055249438278, 0.9998283289949064},
                                {1, "sym", 0.9167202240286861, 0.9991788058849537},
                                {1, "cont", 1.0053936737427394, 1.0000050727073335},
                                {2, "D", 2.041269916294635, 2.000002687021147},
                                {2, "DT", 1.0496453929964105, 1.0002973098453514},
                                {2, "sym", 1.0131830871801109, 1.0005421178195406},
                                {2, "cont", 0.9988922920878781, 1.0000170654756}};
    fit = median = -1.;
    for (auto &r : refs)
        if (r.fKind == kind && op == r.fOp) {
            fit = r.fFit;
            median = r.fMedian;
        }
}

inline std::vector<TaylorTest::TOperatorPanel> TaylorTest::OperatorSlopes(const std::vector<TPanel> &fig3) const {
    const mcc::TPlastic model = CreateMaterial();
    TNumpyRandom rg(fOperatorSeed);
    std::vector<TOperatorPanel> out;
    for (int kind : {ESubcritical, ESupercritical}) {
        // the state eps_0 of Fig. 3 (panel with D of this kind); sigma_0, D_0 and p_c recomputed as in tangentes()
        const TPanel *fig = nullptr;
        for (auto &p : fig3)
            if (p.fKind == kind && !p.fTransposed) fig = &p;
        if (!fig) DebugStop();
        TPanel base(*fig);
        TPZFNMatrix<36, REAL> D0(6, 6, 0.);
        REAL pc = 0.;
        int type = -1;
        if (!Response(model, kind, base.fX0, base.fSigma0, D0, pc, type) || type != kind) DebugStop();
        base.fD0 = D0;
        base.fPc = pc;
        TPZFNMatrix<36, REAL> DT(6, 6, 0.), Dsym(6, 6, 0.), Dc(6, 6, 0.);
        for (int i = 0; i < 6; ++i)
            for (int j = 0; j < 6; ++j) {
                DT(i, j) = D0(j, i);
                Dsym(i, j) = 0.5 * (D0(i, j) + D0(j, i));
            }
        model.ContinuumTangent(base.fSigma0, pc, fV0, Dc);
        const std::vector<std::pair<std::string, TPZFNMatrix<36, REAL>>> ops = {
            {"D", D0}, {"DT", DT}, {"sym", Dsym}, {"cont", Dc}};
        // the same operators returned by the library to the global iterations (SetTangentMode, Table 10)
        const std::map<std::string, TPZPlasticStepModifiedCamClay::ETangentMode> modes = {
            {"D", TPZPlasticStepModifiedCamClay::EConsistentTangent},
            {"DT", TPZPlasticStepModifiedCamClay::ETransposedTangent},
            {"sym", TPZPlasticStepModifiedCamClay::ESymmetricTangent},
            {"cont", TPZPlasticStepModifiedCamClay::EContinuumTangent}};
        for (auto &op : ops) {
            TOperatorPanel res;
            res.fOperator = op.first;
            res.fMatrix = op.second;
            mcc::TPlastic lib(model);
            lib.SetTangentMode(modes.at(op.first));
            TPZTensor<REAL> sl;
            TPZFNMatrix<36, REAL> Dl(6, 6, 0.);
            REAL pcl = 0.;
            int typel = -1;
            if (!Response(lib, kind, base.fX0, sl, Dl, pcl, typel) || typel != kind) DebugStop();
            res.fLibraryDiff = 0.;
            for (int i = 0; i < 6; ++i)
                for (int j = 0; j < 6; ++j)
                    res.fLibraryDiff = std::max(res.fLibraryDiff, std::fabs(Dl(i, j) - op.second.GetVal(i, j)));
            res.fPanel = base;
            res.fPanel.fTransposed = op.first == "DT";
            PerturbWith(model, res.fPanel, op.second, rg);
            out.push_back(res);
        }
    }
    return out;
}

inline void TaylorTest::PostProcessOperators(const std::vector<TOperatorPanel> &ops) const {
    std::vector<std::vector<REAL>> summary;
    for (auto &o : ops) {
        const TPanel &p = o.fPanel;
        std::vector<std::vector<REAL>> rows;
        for (size_t k = 0; k < p.fPairs.size(); ++k) {
            const auto &q = p.fPairs[k];
            const REAL la = std::log(q[0]);
            rows.push_back({q[0], q[1], la, std::log(q[1]), p.fFitIntercept + p.fFitSlope * la, q[2], q[3],
                            p.fSlopes[k]});
        }
        mcc::WriteCSV("taylor_operators_" + KindName(p.fKind) + "_" + o.fOperator + ".csv",
                      {"alpha1", "E1", "log_alpha1", "log_E1", "fit_log_E1", "alpha2", "E2", "pair_slope"}, rows);
        REAL pyfit, pymed;
        OperatorReference(p.fKind, o.fOperator, pyfit, pymed);
        // asymmetry of the operator itself
        REAL nd = 0., na = 0.;
        for (int i = 0; i < 6; ++i)
            for (int j = 0; j < 6; ++j) {
                const REAL d = o.fMatrix.GetVal(i, j), a = o.fMatrix.GetVal(i, j) - o.fMatrix.GetVal(j, i);
                nd += d * d;
                na += a * a;
            }
        const REAL opcode = o.fOperator == "D" ? 0. : (o.fOperator == "DT" ? 1. : (o.fOperator == "sym" ? 2. : 3.));
        summary.push_back({REAL(p.fKind), opcode, p.fP, p.fQ, p.fPc, std::sqrt(na) / std::sqrt(nd), p.fFitSlope,
                           p.fFitIntercept, p.fMedian, p.fMinSlope, p.fMaxSlope, REAL(p.fRejected), pyfit, pymed,
                           o.fLibraryDiff});
    }
    mcc::WriteCSV("taylor_operators_summary.csv",
                  {"kind", "operator", "p_eff", "q", "pc", "asym_operator", "fit_slope", "fit_intercept",
                   "median_slope", "min_slope", "max_slope", "rejected", "python_fit_slope", "python_median_slope",
                   "diff_library_mode"},
                  summary);
}

inline void TaylorTest::PrintOperators(const std::vector<TOperatorPanel> &ops) const {
    std::cout << "\nTaylor slopes of the tangent operators of Sect. 6.7 (last column of Table 10), states of Fig. 3, "
                 "numpy-compatible stream default_rng("
              << fOperatorSeed << ")\n";
    std::cout << "  values: this work {Python, data_tangentes.pkl} [Table 10, median]; operator codes in "
                 "taylor_operators_summary.csv: 0 D, 1 D^T, 2 (D+D^T)/2, 3 continuum\n";
    std::cout << std::left << std::setw(17) << "  state" << std::setw(7) << "op" << std::setw(36) << "fitted slope"
              << "median slope" << std::right << "\n";
    REAL dfit = 0., dmed = 0.;
    for (auto &o : ops) {
        REAL pyfit, pymed;
        OperatorReference(o.fPanel.fKind, o.fOperator, pyfit, pymed);
        dfit = std::max(dfit, std::fabs(o.fPanel.fFitSlope - pyfit));
        dmed = std::max(dmed, std::fabs(o.fPanel.fMedian - pymed));
        std::ostringstream f, m;
        f << std::fixed << std::setprecision(9) << o.fPanel.fFitSlope << " {" << pyfit << "}";
        m << std::fixed << std::setprecision(9) << o.fPanel.fMedian << " {" << pymed << "} ["
          << (o.fOperator == "D" ? "2.00" : "1.00") << "]";
        std::cout << "  " << std::left << std::setw(15) << KindName(o.fPanel.fKind) << std::setw(7) << o.fOperator
                  << std::setw(36) << f.str() << m.str() << std::right << "\n";
    }
    std::cout << "  largest difference with the Python values: fitted slope " << std::scientific << std::setprecision(2)
              << dfit << ", median " << dmed << std::defaultfloat << std::setprecision(6) << "; rejected pairs:";
    for (auto &o : ops) std::cout << " " << o.fPanel.fRejected;
    REAL dlib = 0.;
    for (auto &o : ops) dlib = std::max(dlib, o.fLibraryDiff);
    std::cout << "\n  the operators tested are those returned by TPZPlasticStepModifiedCamClay::SetTangentMode "
                 "(D, DT, sym, cont): largest |difference| "
              << dlib << "\n";
}

inline void TaylorTest::RunAll() {
    std::cout << "Taylor test of the consistent tangent (Sect. 4.5, Fig. 3): M = " << fM << ", lambda = " << fLambda
              << ", kappa = " << fKappa << ", v0 = " << fV0 << ", porous elasticity, nu = " << fNu << "\n";
    std::cout << "converged states sigma_n = -p'_n I, eps_n = 0:";
    for (int kind : {EElastic, ESubcritical, ESupercritical}) {
        const TState st = State(kind);
        std::cout << (kind ? "; " : " ") << KindName(kind) << " p'_n = " << st.fPn << ", p'_c,n = " << st.fPcn
                  << ", eps_0 in [-" << st.fRange << ", " << st.fRange << "]";
    }
    std::cout << "\n";
    // 1. same random numbers as the Python script: reproduces the states and slopes of Fig. 3
    std::cout << "check of the PCG64 transcription against numpy default_rng(2026).bit_generator.random_raw(3): "
              << (TNumpyRandom::SelfTest() ? "ok" : "FAILED") << "\n";
    TNumpyRandom pcg(fSeed);
    const std::vector<TPanel> fig3 = Run(pcg);
    Print(fig3, pcg, true);
    PostProcess(fig3, "taylor_" + pcg.Name());
    // 2. alternative tangent operators at the states of Fig. 3 (Table 10, last column; tangentes() of gen_data.py)
    const std::vector<TOperatorPanel> ops = OperatorSlopes(fig3);
    PrintOperators(ops);
    PostProcessOperators(ops);
    // 3. independent sample with the C++ standard engine: statistical comparison
    TMersenneRandom mt(fSeed);
    const std::vector<TPanel> sample = Run(mt);
    Print(sample, mt, false);
    PostProcess(sample, "taylor_" + mt.Name());
}
