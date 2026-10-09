/**
 * @file TPZPlasticStepModifiedCamClay.h
 * @brief Stress update and consistent tangent operator of the Modified Cam-Clay model with
 * pressure-dependent (porous) elasticity in rotated Haigh-Westergaard space.
 */

#ifndef TPZPLASTICSTEPMODIFIEDCAMCLAY_H
#define TPZPLASTICSTEPMODIFIEDCAMCLAY_H

#include "TPZPlasticBase.h"
#include "TPZPlasticState.h"
#include "TPZElasticResponse.h"
#include "TPZYCModifiedCamClayRHW.h"
#include "TPZTensor.h"
#include "pzfmatrix.h"

/**
 * @ingroup material
 * @brief Stress update of the Modified Cam-Clay model (Algorithm 1 of the article; routine
 * ProjectStressCC of camclay-perf.m with the four-unknown local system) and its consistent tangent operator.
 *
 * The class follows the interface of the NeoPZ plastic steps (TPZPlasticBase) so that it can be
 * used as the template argument of the elastoplastic materials. The steps are:
 *  -# elastic trial state: \f$\sigma_{tr}=s_n+2G\Delta e+p_{tr}I\f$ with
 *     \f$p_{tr}=p_n\exp(-v_0\Delta\varepsilon_v/\kappa)\f$ (porous law) or
 *     \f$p_{tr}=p_n+K_0\Delta\varepsilon_v\f$ (linear law), and \f$G\f$ constant or
 *     \f$G=\frac{3(1-2\nu)}{2(1+\nu)}K(p_n)\f$ frozen in the step (hypoelastic law of FLAC3D/Abaqus);
 *  -# spectral decomposition of \f$\sigma_{tr}\f$ (TPZTensor::EigenSystem);
 *  -# elastic check \f$\Phi(\xi_{tr},\rho_{tr},a_n)\le10^{-11}a_n^2\f$;
 *  -# local Newton iterations of the reduced problem in the angle of the meridian ellipse and the hardening
 *     increment (TPZYCModifiedCamClayRHW::ProjectReduced; the four-unknown system TPZYCModifiedCamClayRHW::ProjectHW
 *     is available as a cross-check, SetLocalSolver), with the region check for \f$\omega\neq1\f$;
 *  -# projected principal stresses \f$\sigma_{proj}=\xi\mathbf1/\sqrt3+\rho\,n\f$;
 *  -# Jacobian of the projection in principal stresses, by implicit differentiation of the converged local system;
 *  -# consistent tangent operator in spectral form with the rotational correction (ComputedDep of HWTools.m),
 *     with \f$\kappa_{ij}=2G\rho/\rho_{tr}\f$ (see ComputedDep),
 *     assembled column by column, i.e. already transposed with respect to the row-wise storage of
 *     the Wolfram Language routine (Sect. 3 of the article).
 *
 * The class also implements a linear elastic model (SetLinearElastic) used in the Terzaghi
 * consolidation example, so that the same u-p material serves all the examples.
 *
 * State conventions (TPZPlasticState):
 *  - m_eps_t: total strain of the last converged step (engineering shear components
 *    \f$\gamma_{xy}=2\varepsilon_{xy}\f$, as in the strain tensors of TPZElasticResponse);
 *  - m_hardening: preconsolidation pressure \f$p_c\f$ (positive);
 *  - m_m_type: 0 elastic, 1 plastic in the subcritical region, 2 plastic in the supercritical region;
 *  - m_eps_p: not used (zero), the plastic strain is not defined by the hypoelastic predictor;
 *  - fmatprop[0]: specific volume \f$v_0\f$ of the point (if fmatprop is empty, the default
 *    specific volume of the class is used).
 *
 * The stress update is incremental (hypoelastic predictor): ApplyStrainComputeSigma receives in
 * @c sigma the stress of the last converged step and returns the updated stress. Strains are stored in
 * TPZTensor objects with ENGINEERING shear components (\f$\gamma_{xy}=2\varepsilon_{xy}\f$), the
 * convention of TPZElasticResponse in this version of NeoPZ and of the Voigt map
 * \f$V_\varepsilon\f$ of eq. (6) of the article. The tangent is
 * \f$\partial\sigma/\partial\varepsilon\f$ in the Voigt order of TPZTensor (XX, XY, XZ, YY, YZ, ZZ)
 * with respect to these engineering strains (same convention of TPZElasticResponse::De): the entry
 * (l, k) is \f$\partial\sigma_l/\partial\varepsilon_k\f$, i.e. column k is the response to the unit
 * engineering strain \f$e_k\f$.
 *
 * The operator returned to the global Newton iterations is selected with SetTangentMode (comparison of
 * the tangent operators of Sect. 6.7 and Table 10 of the article, option TANGENT['mode'] of gen_data.py):
 * the consistent tangent (default), its transpose, its symmetric part, the continuum elastoplastic
 * operator or central differences of the stress update. The stress and the state do not depend on the
 * mode, only the tangent does; the one exception is the central-difference mode, in which a failure of one of
 * the extra stress updates makes the whole call fail (see ApplyStrainComputeSigma).
 */
class TPZPlasticStepModifiedCamClay : public TPZPlasticBase {
public:

    /** @brief Constitutive model */
    enum EModel {
        EModifiedCamClay = 0, ///< Modified Cam-Clay plasticity (default)
        ELinearElastic = 1    ///< linear isotropic elasticity (E, nu), no plasticity
    };

    /** @brief Shear modulus of the elastic predictor */
    enum EShear {
        EConstantG = 0, ///< constant shear modulus G
        EPoisson = 1    ///< G from a constant Poisson ratio and the bulk modulus of the start of the step
    };

    /**
     * @brief Tangent operator returned by ApplyStrainComputeSigma and ApplyStrainComputeDep
     * (Sect. 6.7 and Table 10 of the article; option TANGENT['mode'] of gen_data.py, whose names are given
     * in parentheses). In every mode the matrix is \f$\partial\sigma/\partial\varepsilon\f$ with respect
     * to engineering strains in the Voigt order of TPZTensor, entry (l, k) = row l, column k.
     *
     * The modes apply to the Modified Cam-Clay model only: the linear elastic model (SetLinearElastic)
     * always returns its elastic operator, as the elastic material of fe_user.py, which does not go through
     * the switch of gen_data.py. When the stress update fails (LastProjectionFailed), the tangent is the
     * elastic trial operator whatever the mode.
     */
    enum ETangentMode {
        /** ('D') consistent tangent (8) with the rotational correction (9) (ComputedDep); default */
        EConsistentTangent = 0,
        /** ('DT') transpose of the consistent tangent: what the column-wise Voigt assembly of the Wolfram
         * Language routine returns without the correction of Sect. 3 (also selected by SetTransposedTangent) */
        ETransposedTangent = 1,
        /** ('sym') symmetric part \f$(D+D^T)/2\f$ of the consistent tangent at every point (what a code
         * with a symmetric solver uses) */
        ESymmetricTangent = 2,
        /** ('cont') continuum elastoplastic operator at the end-of-step state at the plastic points
         * (ContinuumTangent); the elastic points keep the elastic trial operator */
        EContinuumTangent = 3,
        /** ('fd') central differences of the stress update at every point (elastic points included), with
         * step FiniteDifferenceStep() (1e-7) on the engineering strains (FiniteDifferenceTangent) */
        EFiniteDifferenceTangent = 4
    };

    /** @brief Short name of a tangent mode, as in gen_data.py: "D", "DT", "sym", "cont" or "fd" */
    static const char *TangentModeName(ETangentMode mode);

    TPZPlasticStepModifiedCamClay();

    TPZPlasticStepModifiedCamClay(const TPZPlasticStepModifiedCamClay &cp) = default;

    TPZPlasticStepModifiedCamClay &operator=(const TPZPlasticStepModifiedCamClay &cp) = default;

    virtual ~TPZPlasticStepModifiedCamClay() = default;

    /** @name Set up of the model */
    /** @{ */

    /**
     * @brief Modified Cam-Clay parameters (Table 1 of the article)
     * @param M slope of the critical state line
     * @param lambda slope of the normal compression line
     * @param kappa slope of the swelling line
     * @param pt tensile strength (default 0)
     * @param omega shape parameter of the subcritical region (default 1)
     */
    void SetModifiedCamClay(REAL M, REAL lambda, REAL kappa, REAL pt = 0., REAL omega = 1.);

    /** @brief Porous (logarithmic) volumetric law (14) (default) */
    void SetPorousElasticity() { fYC.SetVolumetricLaw(TPZYCModifiedCamClayRHW::EPorous); }

    /** @brief Linear volumetric law with bulk modulus K0 */
    void SetLinearVolumetric(REAL K0) { fYC.SetVolumetricLaw(TPZYCModifiedCamClayRHW::ELinear, K0); }

    /** @brief Constant shear modulus G (15) */
    void SetConstantShearModulus(REAL G) { fShear = EConstantG; fG = G; }

    /** @brief Shear modulus from a constant Poisson ratio and the bulk modulus at the start of the step (15) */
    void SetPoissonRatio(REAL nu) { fShear = EPoisson; fNu = nu; }

    /** @brief Exact (default) or frozen integration of the porous law (comparison with the frozen bulk modulus) */
    void SetPorousIntegration(TPZYCModifiedCamClayRHW::EPorousIntegration integ) { fYC.SetPorousIntegration(integ); }

    /** @brief Local system solved at the plastic points: reduced (two unknowns, default) or full (four unknowns) */
    void SetLocalSolver(TPZYCModifiedCamClayRHW::ELocalSolver solver) { fYC.SetLocalSolver(solver); }

    /** @brief Linear isotropic elastic model (no plasticity): \f$\sigma=\sigma_n+D_e\Delta\varepsilon\f$ */
    void SetLinearElastic(REAL E, REAL nu);

    /**
     * @brief Returns the transpose of the consistent tangent (what the column-wise Voigt assembly of
     * the Wolfram Language routine returns without the correction of Sect. 3). Used only to
     * reproduce the comparisons of Sects. 4.5 and 6.7 (Taylor test and global iterations).
     *
     * Kept for compatibility: true selects ETransposedTangent and false EConsistentTangent (SetTangentMode).
     */
    void SetTransposedTangent(bool transposed) { fTangentMode = transposed ? ETransposedTangent : EConsistentTangent; }

    /** @brief Selects the tangent operator returned by the stress update (see ETangentMode; default EConsistentTangent) */
    void SetTangentMode(ETangentMode mode) { fTangentMode = mode; }

    /**
     * @brief Step h of the central differences of EFiniteDifferenceTangent (default 1e-7, as in gen_data.py);
     * column k is \f$(\sigma(\varepsilon+he_k)-\sigma(\varepsilon-he_k))/(2h)\f$. h must be positive and
     * finite (DebugStop otherwise).
     */
    void SetFiniteDifferenceStep(REAL h);

    /** @brief Specific volume used when the plastic state does not hold one (fmatprop empty) */
    void SetDefaultSpecificVolume(REAL v0) { fV0Default = v0; }

    /** @} */

    /** @name Access */
    /** @{ */
    EModel Model() const { return fModel; }
    EShear Shear() const { return fShear; }
    REAL ShearModulus() const { return fG; }
    REAL PoissonRatio() const { return fNu; }
    /** @brief True if the tangent mode is ETransposedTangent */
    bool TransposedTangent() const { return fTangentMode == ETransposedTangent; }
    /** @brief Tangent operator returned by the stress update */
    ETangentMode TangentMode() const { return fTangentMode; }
    /** @brief Step of the central differences of EFiniteDifferenceTangent */
    REAL FiniteDifferenceStep() const { return fFDStep; }
    const TPZYCModifiedCamClayRHW &YC() const { return fYC; }
    TPZYCModifiedCamClayRHW &YC() { return fYC; }

    /** @brief Specific volume of the current state */
    REAL SpecificVolume() const;

    /** @brief Sets the specific volume of the current state (stored in fmatprop[0]) */
    void SetSpecificVolume(REAL v0);

    /** @brief Number of local Newton iterations of the last projection (0 for elastic steps) */
    int LastNewtonIterations() const { return fLastNewtonIterations; }

    /** @brief True if the last call of ApplyStrainComputeSigma failed (local projection did not converge) */
    bool LastProjectionFailed() const { return fFailed; }
    /** @} */

    /**
     * @brief Stress update for a given total strain.
     * @param epsTotal total strain at the end of the step (Voigt components with engineering shear strains)
     * @param[in,out] sigma on input the stress of the last converged step; on output the updated stress
     * @param tangent if not null, the 6x6 tangent with respect to engineering strains selected by
     * SetTangentMode (default: the consistent tangent)
     *
     * The state (fN) is updated: m_eps_t, m_hardening (pc), m_m_type, and the elastic response
     * (GetElasticResponse) receives the pair (E, nu) equivalent to the trial moduli. On failure of the
     * local projection (no convergence, wrong region or singular Jacobian of the projection) the state
     * and the elastic response are not changed, sigma is returned unchanged, LastProjectionFailed()
     * becomes true and the tangent is set to the elastic trial operator.
     *
     * Tangent modes (computed only when @c tangent is not null):
     *  - EConsistentTangent, ETransposedTangent, ESymmetricTangent: \f$D\f$, \f$D^T\f$ or
     *    \f$(D+D^T)/2\f$, with \f$D\f$ the consistent tangent (elastic trial operator at elastic points);
     *  - EContinuumTangent: ContinuumTangent at the updated stress and preconsolidation pressure at the
     *    plastic points (m_m_type 1 or 2), the elastic trial operator at the elastic points;
     *  - EFiniteDifferenceTangent: FiniteDifferenceTangent, i.e. 12 extra stress updates from the same
     *    converged state (state and stress on input), whose results are discarded: the state, the
     *    elastic response, the stress returned and LastNewtonIterations() are those of the update with
     *    epsTotal. If one of the extra updates fails, the call fails as a whole (as the exception of the
     *    Python code, caught by fe_user.py, which marks the point as failed).
     */
    void ApplyStrainComputeSigma(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma,
                                 TPZFMatrix<REAL> *tangent = NULL) override;

    /** @brief Same as ApplyStrainComputeSigma with the tangent (of the mode set by SetTangentMode) */
    void ApplyStrainComputeDep(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma, TPZFMatrix<REAL> &Dep) override;

    /** @brief Not available for the incremental (hypoelastic) model: DebugStop */
    void ApplyStrain(const TPZTensor<REAL> &epsTotal) override;

    /** @brief Not available for this model: DebugStop */
    void ApplyLoad(const TPZTensor<REAL> &sigma, TPZTensor<REAL> &epsTotal) override;

    void SetState(const TPZPlasticState<REAL> &state) override { fN = state; }

    TPZPlasticState<REAL> GetState() const override { return fN; }

    TPZPlasticCriterion &GetYC() override { return fYC; }

    /**
     * @brief Value of the yield function for the stress of the last stress update
     * (the argument is not used, the model is incremental)
     */
    void Phi(const TPZTensor<REAL> &epsTotal, TPZVec<REAL> &phi) const override;

    /** @brief Elastic response of the last step (equivalent pair E, nu of K_tr and G) */
    void SetElasticResponse(TPZElasticResponse &ER) override { fER = ER; }

    TPZElasticResponse GetElasticResponse() const override { return fER; }

    const char *Name() const override { return "TPZPlasticStepModifiedCamClay"; }

    void Print(std::ostream &out) const override;

    int ClassId() const override;

    void Write(TPZStream &buf, int withclassid) const override;

    void Read(TPZStream &buf, void *context) override;

    /**
     * @brief Elastic predictor (16) (TrialStressCC)
     * @param deps strain increment (engineering shear strains)
     * @param sigman stress of the last converged step
     * @param v0 specific volume
     * @param[out] sigtr trial stress
     * @param[out] Ktr trial bulk modulus \f$dp_{tr}/d\varepsilon_v\f$
     * @param[out] G shear modulus of the step
     */
    void TrialStress(const TPZTensor<REAL> &deps, const TPZTensor<REAL> &sigman, REAL v0, TPZTensor<REAL> &sigtr,
                     REAL &Ktr, REAL &G) const;

    /**
     * @brief Consistent tangent (8) in spectral form with the rotational correction (9) (ComputedDep)
     * @param Dproj Jacobian of the projection in principal stresses (22)
     * @param eigvec eigenvectors of the trial stress (eigvec[i] is the i-th eigenvector)
     * @param Ktr trial bulk modulus
     * @param G shear modulus
     * @param ratio \f$\rho/\rho_{tr}=1/(1+6G\Delta\gamma/M^2)\f$ (1 for an elastic state)
     * @param[out] Dep 6x6 consistent tangent (columns are the responses to unit engineering strains)
     *
     * The coefficients of the rotational correction are
     * \f$\kappa_{ij}=(\sigma_i-\sigma_j)/(\varepsilon^{tr}_i-\varepsilon^{tr}_j)\f$. Since the return of
     * the Modified Cam-Clay model is radial in the deviatoric plane and
     * \f$\varepsilon^{tr}_i-\varepsilon^{tr}_j=(\sigma^{tr}_i-\sigma^{tr}_j)/(2G)\f$, all of them are
     * equal to \f$2G\rho/\rho_{tr}\f$. The closed form is used instead of the quotient of the Wolfram
     * Language routine, which is corrupted by cancellation when two trial eigenvalues are close (the
     * routine switches to the limit only when the difference is below 1e-15 in absolute value). Both
     * forms agree to round-off for distinct eigenvalues.
     */
    static void ComputedDep(const TPZFMatrix<REAL> &Dproj, const TPZVec<TPZManVector<REAL, 3>> &eigvec, REAL Ktr,
                            REAL G, REAL ratio, TPZFMatrix<REAL> &Dep);

    /**
     * @brief Continuum (non-algorithmic) elastoplastic operator at the state (sigma, pc)
     * (continuum_tangent of camclay_hw.py, Sect. 6.7 of the article)
     * @param sigma stress (tension positive)
     * @param pc preconsolidation pressure
     * @param v0 specific volume
     * @param[out] Dc 6x6 operator, engineering strains, Voigt order of TPZTensor
     *
     * \f$D_c=\mathbb{C}-\dfrac{(\mathbb{C}n)(n^T\mathbb{C})}{n^T\mathbb{C}n+h}\f$ with:
     *  - \f$\mathbb{C}=\mathbb{C}(K,G)\f$ (ElasticOperator) with the moduli at the state: \f$K=-v_0p/\kappa\f$
     *    (porous law) or \f$K_0\f$ (linear law), \f$G\f$ constant or \f$\frac{3(1-2\nu)}{2(1+\nu)}K\f$;
     *  - \f$n=\partial\Phi/\partial\sigma=\frac{2\bar p}{3b^2}\mathbf{m}+\frac{3}{M^2}s\f$ in engineering form
     *    (shear components doubled), \f$\bar p=p-p_t+a\f$, \f$a=(p_c+p_t)/(1+\omega)\f$, \f$b=1\f$ if
     *    \f$\bar p\ge0\f$ and \f$\omega\f$ otherwise, \f$s\f$ the deviatoric stress;
     *  - \f$h=-\frac{\partial\Phi}{\partial a}\frac{da}{d\gamma}\f$ with
     *    \f$\frac{\partial\Phi}{\partial a}=\frac{2\bar p}{b^2}-2a\f$ and
     *    \f$\frac{da}{d\gamma}=-H\frac{2\bar p}{b^2}\f$, \f$H=\frac{v_0p_c}{(\lambda-\kappa)(1+\omega)}\f$.
     */
    void ContinuumTangent(const TPZTensor<REAL> &sigma, REAL pc, REAL v0, TPZFMatrix<REAL> &Dc) const;

    /**
     * @brief Tangent by central differences of the stress update (option 'fd' of gen_data.py)
     * @param epsTotal total strain at the end of the step
     * @param sigman stress of the last converged step
     * @param[out] Dfd 6x6 matrix, column k = \f$(\sigma(\varepsilon+he_k)-\sigma(\varepsilon-he_k))/(2h)\f$
     * with \f$h\f$ = FiniteDifferenceStep() and \f$e_k\f$ the unit engineering strain in the Voigt order
     * of TPZTensor
     * @return false if one of the 12 stress updates failed
     *
     * Each update starts from the current state (GetState), which is restored afterwards together with
     * the elastic response, the stress of Phi and the counters: the object is not changed.
     */
    bool FiniteDifferenceTangent(const TPZTensor<REAL> &epsTotal, const TPZTensor<REAL> &sigman,
                                 TPZFMatrix<REAL> &Dfd);

    /**
     * @brief Isotropic elastic operator \f$\mathbb{C}(K,G)\f$ of (10), engineering strains
     * @param K bulk modulus
     * @param G shear modulus
     * @param[out] C 6x6 matrix
     */
    static void ElasticOperator(REAL K, REAL G, TPZFMatrix<REAL> &C);

    /**
     * @brief Spectral decomposition of a symmetric tensor: eigenvalues in descending order and an
     * orthonormal set of eigenvectors. Uses TPZTensor::EigenSystem and falls back to cyclic Jacobi
     * rotations when the native decomposition does not return an orthonormal basis (repeated
     * eigenvalues).
     */
    static void EigenSystem(const TPZTensor<REAL> &sigma, TPZManVector<REAL, 3> &eigval,
                            TPZManVector<TPZManVector<REAL, 3>, 3> &eigvec);

protected:
    /**
     * @brief Stress update with the consistent tangent (Algorithm 1): the body of ApplyStrainComputeSigma
     * in the mode EConsistentTangent (see ApplyStrainComputeSigma for the arguments and the failure case)
     */
    void StressUpdate(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma, TPZFMatrix<REAL> *tangent);

    /** @brief Yield criterion and local projection */
    TPZYCModifiedCamClayRHW fYC;
    /** @brief Equivalent elastic response of the last step (E, nu from K_tr and G) or of the linear elastic model */
    TPZElasticResponse fER;
    /** @brief Constitutive model */
    EModel fModel;
    /** @brief Option of the shear modulus */
    EShear fShear;
    /** @brief Constant shear modulus */
    REAL fG;
    /** @brief Poisson ratio of the hypoelastic shear modulus */
    REAL fNu;
    /** @brief Specific volume used when the state does not hold one */
    REAL fV0Default;
    /** @brief Tangent operator returned by the stress update (comparison of Sect. 6.7) */
    ETangentMode fTangentMode;
    /** @brief Step of the central differences of EFiniteDifferenceTangent */
    REAL fFDStep;
    /** @brief Plastic state of the integration point (see the class description) */
    TPZPlasticState<REAL> fN;
    /** @brief Stress of the last stress update (used by Phi) */
    TPZTensor<REAL> fSigma;
    /** @brief Number of local iterations of the last projection */
    int fLastNewtonIterations;
    /** @brief Failure flag of the last projection */
    bool fFailed;
};

#endif
