/**
 * @file TPZMatPoroElastoPlasticUP.h
 * @brief Monolithic mixed u-p material for Biot consolidation with elastoplastic effective stress
 * (plane strain, axisymmetry and three dimensions with a single six-component Voigt operator).
 */

#ifndef TPZMATPOROELASTOPLASTICUP_H
#define TPZMATPOROELASTOPLASTICUP_H

#include "TPZMatBase.h"
#include "TPZMatCombinedSpaces.h"
#include "TPZMatWithMem.h"
#include "TPZMaterialDataT.h"
#include "TPZBndCondT.h"
#include "TPZElastoPlasticMem.h"
#include "TPZTensor.h"
#include "pzfmatrix.h"
#include <functional>
#include <atomic>

class TPZCompMesh;

/**
 * @ingroup matcombinedspaces
 * @brief Non-template interface of the u-p poro-elastoplastic materials, used by the incremental
 * driver TPZPoroElastoPlasticUPAnalysis (time step, load factor, boundary condition types and
 * failure of the local stress update).
 */
class TPZMatPoroElastoPlasticUPBase {
public:

    /** @brief Kinematic model of the strain-displacement operator, eq. (28) of the article */
    enum EKinematics {
        EPlaneStrain = 0,   ///< 2D plane strain (rows XZ, YZ and ZZ of B vanish)
        EAxisymmetric = 1,  ///< 2D axisymmetric, x = r and y = z; the ZZ row carries the hoop strain u_r/r
        EThreeDimensional = 2 ///< 3D
    };

    /**
     * @brief Boundary condition types of the u-p material
     *
     * The values of the boundary condition are given in Val2 (and the components in Val1 for the
     * directional condition). A forcing function of the boundary condition, if present, replaces Val2.
     *
     * The Dirichlet conditions are imposed by penalty (BigNumber) in Contribute. The driver
     * TPZPoroElastoPlasticUPAnalysis eliminates their equations by default (SetEliminateDirichlet) and
     * imposes the values directly on the solution; in that case the forcing functions of the Dirichlet
     * conditions are evaluated at the nodes.
     */
    enum EBCType {
        EDirichletU = 0,        ///< all displacement components prescribed: u = Val2[0..dim-1]
        ENeumannU = 1,          ///< traction t = lambda Val2[0..dim-1], scaled by the load factor lambda
        EDirichletP = 2,        ///< pore pressure prescribed: p = Val2[0] (drained boundary)
        EDirichletUDirectional = 3, ///< components d with Val1(d,d) != 0 prescribed: u_d = Val2[d]
        ENeumannUFixed = 4      ///< traction t = Val2[0..dim-1], not scaled by the load factor
    };

    TPZMatPoroElastoPlasticUPBase() = default;

    /** @brief Copy constructor (the failure counter is not copied) */
    TPZMatPoroElastoPlasticUPBase(const TPZMatPoroElastoPlasticUPBase &cp)
        : fKinematics(cp.fKinematics), fTimeStep(cp.fTimeStep), fLoadFactor(cp.fLoadFactor),
          fExternalOnly(cp.fExternalOnly), fBig(cp.fBig), fNFailed(0) {}

    virtual ~TPZMatPoroElastoPlasticUPBase() = default;

    /** @brief Time step \f$\Delta t\f$ of the current step (0: undrained step, no flow) */
    void SetTimeStep(REAL dt) { fTimeStep = dt; }
    REAL TimeStep() const { return fTimeStep; }

    /** @brief Load factor \f$\lambda\f$ of the boundary tractions of type ENeumannU */
    void SetLoadFactor(REAL lambda) { fLoadFactor = lambda; }
    REAL LoadFactor() const { return fLoadFactor; }

    /**
     * @brief When true, Contribute assembles only the external forces \f$f_b+\lambda f_t\f$ (body forces of
     * the material and tractions of the boundary conditions), used for the normalization of the residual
     */
    void SetAssembleExternalForcesOnly(bool flag) { fExternalOnly = flag; }
    bool AssembleExternalForcesOnly() const { return fExternalOnly; }

    /** @brief Number of integration points whose local projection failed since the last reset */
    int NFailedProjections() const { return fNFailed.load(); }
    void ResetFailedProjections() { fNFailed = 0; }

    /** @brief Penalty number of the Dirichlet conditions (used when their equations are not eliminated) */
    void SetBigNumber(REAL big) { fBig = big; }
    REAL BigNumber() const { return fBig; }

    /** @brief Kinematic model */
    EKinematics Kinematics() const { return fKinematics; }

    /** @brief Spatial dimension of the displacement field (2 or 3) */
    int DimensionU() const { return fKinematics == EThreeDimensional ? 3 : 2; }

    /** @brief True if the boundary condition type prescribes displacement components */
    static bool IsDirichletU(int bctype) { return bctype == EDirichletU || bctype == EDirichletUDirectional; }

    /** @brief True if the boundary condition type prescribes the pore pressure */
    static bool IsDirichletP(int bctype) { return bctype == EDirichletP; }

    /** @brief Material id of the domain material (the TPZMaterial id) */
    virtual int MaterialId() const = 0;

protected:
    EKinematics fKinematics = EPlaneStrain;
    REAL fTimeStep = 0.;
    REAL fLoadFactor = 1.;
    bool fExternalOnly = false;
    REAL fBig = 1.e12;
    std::atomic<int> fNFailed{0};
};

/**
 * @ingroup matcombinedspaces
 * @brief Coupled u-p material of Sect. 5 of the article (Biot consolidation, monolithic Newton).
 *
 * Unknowns: displacement \f$u\f$ (datavec[0], H1, quadratic) and pore pressure \f$p_w\f$ (datavec[1], H1,
 * linear), positive in compression. The total stress is \f$\sigma=\sigma'-\alpha_B p_w I\f$ (tension
 * positive). The material assembles the residuals of the backward-Euler step (26)
 * \f[ R_u = F_{int}(U) - QP - f_b - \lambda f_t, \qquad
 *     R_p = Q^T(U-U_n) + S(P-P_n) + \Delta t\,(HP - f_g), \f]
 * and the monolithic tangent (27)
 * \f[ \begin{bmatrix} K_T & -Q\\ Q^T & S+\Delta t H\end{bmatrix}, \qquad
 *     K_T=\int B^T\mathbb{D}B\,d\Omega, \f]
 * with \f$Q=\int\alpha_B B^T m N_p\f$, \f$S=\int M_B^{-1}N_p^TN_p\f$, \f$H=\int k\nabla N_p^T\nabla N_p\f$,
 * \f$f_g=\int k\nabla N_p^T\rho_w g\f$ and \f$f_b=\int N_u^Tb\f$. As usual in NeoPZ the element vector is
 * ef = -R and the element matrix is the tangent, so that the analysis solves ek du = ef.
 *
 * The strain-displacement operator has six rows in the Voigt order of TPZTensor (XX, XY, XZ, YY, YZ, ZZ)
 * with engineering shear strains for plane strain, axisymmetry and 3D, eq. (28); in axisymmetry
 * \f$d\Omega\f$ is multiplied by \f$2\pi r\f$.
 *
 * The mesh solution holds the TOTAL displacement and pore pressure. The memory of each integration point
 * (TMEM, by default TPZElastoPlasticMem) holds the state of the last converged step:
 *  - m_sigma: effective stress \f$\sigma'_n\f$;
 *  - m_elastoplastic_state.m_eps_t: total strain \f$\varepsilon_n\f$ (engineering shear components);
 *  - m_elastoplastic_state.m_hardening: preconsolidation pressure \f$p_{c,n}\f$;
 *  - m_elastoplastic_state.fmatprop[0]: specific volume \f$v_0\f$;
 *  - m_elastoplastic_state.m_m_type: 0 elastic, 1 subcritical plastic, 2 supercritical plastic;
 *  - m_elastoplastic_state.fpressure: pore pressure \f$p_{w,n}\f$ at the point.
 * The memory is updated only when TPZMatWithMem::fUpdateMem is set (after convergence of the step).
 *
 * @tparam T plastic step with the interface of TPZPlasticStepModifiedCamClay: ApplyStrainComputeSigma
 * receives the converged stress in @c sigma and returns the tangent with respect to engineering strains;
 * LastProjectionFailed() reports a failure of the local projection.
 * @tparam TMEM memory type (TPZElastoPlasticMem)
 */
template <class T, class TMEM = TPZElastoPlasticMem>
class TPZMatPoroElastoPlasticUP : public TPZMatBase<STATE, TPZMatCombinedSpacesT<STATE>, TPZMatWithMem<TMEM>>,
                                  public TPZMatPoroElastoPlasticUPBase {
    using TBase = TPZMatBase<STATE, TPZMatCombinedSpacesT<STATE>, TPZMatWithMem<TMEM>>;

public:

    /** @brief Post-processing variables (names accepted by VariableIndex in parentheses) */
    enum ESolutionVar {
        EDisplacement = 1,  ///< displacement vector ("Displacement")
        EPorePressure = 2,  ///< pore pressure ("PorePressure" or "Pressure")
        EDisplacementX = 3, ///< ("DisplacementX")
        EDisplacementY = 4, ///< ("DisplacementY")
        EDisplacementZ = 5  ///< ("DisplacementZ")
    };

    /** @brief Default constructor (plane strain, id 0) */
    TPZMatPoroElastoPlasticUP();

    /**
     * @brief Constructor
     * @param id material id
     * @param kinematics plane strain, axisymmetric or 3D
     */
    TPZMatPoroElastoPlasticUP(int id, EKinematics kinematics);

    /** @brief Copy constructor (the memory is copied, as in TPZMatWithMem) */
    TPZMatPoroElastoPlasticUP(const TPZMatPoroElastoPlasticUP &cp);

    virtual ~TPZMatPoroElastoPlasticUP() = default;

    /** @name Set up */
    /** @{ */

    /** @brief Plastic model; it is also used to build the default memory item */
    void SetPlasticModel(const T &model);

    /** @brief Plastic model */
    T &PlasticModel() { return fPlasticity; }
    const T &PlasticModel() const { return fPlasticity; }

    /**
     * @brief Biot coefficient and inverse of the Biot modulus
     * @param alpha Biot coefficient \f$\alpha_B\f$
     * @param invBiotModulus \f$1/M_B\f$ (0 for incompressible constituents; \f$n/K_f\f$ for incompressible grains)
     */
    void SetBiot(REAL alpha, REAL invBiotModulus) { fAlpha = alpha; fInvBiotModulus = invBiotModulus; }

    /** @brief Mobility k (hydraulic conductivity divided by the unit weight of water) */
    void SetPermeability(REAL k) { fPermeability = k; }

    /** @brief Body force of the saturated soil b (force per unit volume) */
    void SetBodyForce(const TPZVec<REAL> &b);

    /** @brief Weight of the fluid per unit volume, \f$\rho_w g\f$ (vector) */
    void SetFluidWeight(const TPZVec<REAL> &rhowg);

    /**
     * @brief Order of the integration rule passed to the Gauss rules (3: 2x2 or 2x2x2 points (reduced),
     * 4: 3x3 or 3x3x3 points (full)); 0 uses the NeoPZ default
     *
     * The order defines the number of memory items of each element, so it must be set before the
     * multiphysics space is built (TPZMultiphysicsCompMesh::BuildMultiphysicsSpaceWithMemory).
     */
    void SetIntegrationOrder(int order) {
        if (order != fIntegrationOrder && this->GetMemory() && this->GetMemory()->NElements() > 0) {
            std::cout << __PRETTY_FUNCTION__ << ": the integration order cannot change after the memory "
                      << "of the integration points was created" << std::endl;
            DebugStop();
        }
        fIntegrationOrder = order;
    }
    int IntegrationOrder() const { return fIntegrationOrder; }

    /** @} */

    REAL Alpha() const { return fAlpha; }
    REAL InvBiotModulus() const { return fInvBiotModulus; }
    REAL Permeability() const { return fPermeability; }

    int MaterialId() const override { return this->Id(); }

    /** @name Access to the integration points */
    /** @{ */

    /**
     * @brief Calls f for every integration point of the elements of this material, with the element,
     * the index of the point, its physical coordinates, its parametric coordinates, the weight times
     * the determinant of the Jacobian (without the 2 pi r factor) and the memory item.
     * The integration rule is the same used by the multiphysics element in the assembly. Only the
     * multiphysics elements found directly in the mesh are visited (not the elements inside condensed
     * elements or submeshes).
     */
    void ForEachIntegrationPoint(TPZCompMesh *mesh,
                                 const std::function<void(TPZCompEl *cel, int ip, const TPZVec<REAL> &x,
                                                          const TPZVec<REAL> &qsi, REAL weight, TMEM &mem)> &f);

    /**
     * @brief Initializes the memory of every integration point with init(x, mem), x physical coordinates.
     * The default memory item (built by SetPlasticModel) is copied first.
     */
    void InitializeMemory(TPZCompMesh *mesh, const std::function<void(const TPZVec<REAL> &x, TMEM &mem)> &init);
    /** @} */

    /** @name TPZMaterial interface */
    /** @{ */
    int Dimension() const override { return DimensionU(); }

    int NStateVariables() const override { return DimensionU(); }

    int IntegrationRuleOrder(const TPZVec<int> &elPMaxOrder) const override;

    void FillDataRequirements(TPZVec<TPZMaterialDataT<STATE>> &datavec) const override;

    void FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE>> &datavec) const override;

    void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                    TPZFMatrix<STATE> &ef) override;

    void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ef) override;

    void ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                      TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) override;

    void ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ef,
                      TPZBndCondT<STATE> &bc) override;

    int VariableIndex(const std::string &name) const override;

    int NSolutionVariables(int var) const override;

    void Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<STATE> &sol) override;

    TPZMaterial *NewMaterial() const override;

    std::string Name() const override { return "TPZMatPoroElastoPlasticUP"; }

    void Print(std::ostream &out = std::cout) const override;

    int ClassId() const override;

    void Write(TPZStream &buf, int withclassid) const override;

    void Read(TPZStream &buf, void *context) override;
    /** @} */

    /**
     * @brief Six-row strain-displacement operator (28) at a point, engineering shear strains
     * @param phi shape functions of the displacement (n x 1)
     * @param gradphi global gradients of the shape functions (dim x n)
     * @param x physical coordinates (radius x[0] in axisymmetry)
     * @param[out] B 6 x (dim n) matrix
     */
    void ComputeB(const TPZFMatrix<REAL> &phi, const TPZFMatrix<REAL> &gradphi, const TPZVec<REAL> &x,
                  TPZFMatrix<REAL> &B) const;

protected:
    /** @brief Contribution with or without the tangent */
    void ContributeInternal(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> *ek,
                            TPZFMatrix<STATE> &ef);

    /** @brief Boundary contribution with or without the tangent */
    void ContributeBCInternal(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> *ek,
                              TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc);

    /** @brief Global gradients (dim x n) of shape functions given in the axes of the element */
    void GlobalGradients(const TPZMaterialDataT<STATE> &data, TPZFMatrix<REAL> &grad) const;

    /** @brief Plastic model used to build the stress update at each point */
    T fPlasticity;
    /** @brief Biot coefficient */
    REAL fAlpha = 1.;
    /** @brief Inverse of the Biot modulus */
    REAL fInvBiotModulus = 0.;
    /** @brief Mobility k */
    REAL fPermeability = 0.;
    /** @brief Body force of the saturated soil */
    TPZManVector<REAL, 3> fBodyForce;
    /** @brief Weight of the fluid per unit volume */
    TPZManVector<REAL, 3> fFluidWeight;
    /** @brief Order of the integration rule (0: NeoPZ default) */
    int fIntegrationOrder = 0;
};

#endif
