/**
 * @file TPZUnsatAnalysis.h
 * @brief Incremental analysis of the unsaturated u - p_w - p_a model (TPZMatUnsaturatedUPwPa with the Mohr-Coulomb PV2
 * step): backward Euler step solved by Newton with the full Jacobian, Dirichlet equations eliminated by an equation
 * filter, nodal water inflows and memory of the converged state.
 *
 * Moved from Projects/UnsatCommon/UnsatAnalysis.h (port of Projects/SlopeRainfall/UnsatSlope.h), with the constraints
 * and the boundary nodes taken from the boundary ids of the mesh. Requires BUILD_PLASTICITY_MATERIALS.
 */
#ifndef TPZUNSATANALYSIS_H
#define TPZUNSATANALYSIS_H

#include "TPZLinearAnalysis.h"
#include "TPZMultiphysicsCompMesh.h"
#include "Plasticity/TPZElasticResponse.h"
#include "Plasticity/TPZMatUnsaturatedUPwPa.h"
#include "Plasticity/TPZPlasticStepVoigt.h"
#include "Plasticity/TPZYCMohrCoulombPV2.h"

#include <array>
#include <functional>
#include <map>
#include <set>
#include <vector>

/**
 * @ingroup analysis
 * @brief Newton solution of the backward-Euler steps of TPZMatUnsaturatedUPwPa (2D, triangles, u P2 and p_w, p_a P1).
 *
 * The mesh solution holds the TOTAL values (u, p_w, p_a). The constrained equations of a step (map equation -> value)
 * are eliminated with the equation filter of the structural matrix and their values are written in the solution. With
 * frozen pressures (factor of safety, dry soil) only the momentum balance is assembled and all the pressure equations
 * must be in the map of constrained equations.
 *
 * Convergence: relative momentum residual |R_u| / |weight| < tolu (|weight| replaced by SetForceReference when the
 * weight is null) and, with free pressures, largest pressure correction < tolp. The memory of the integration points
 * and the converged solution X() are updated only by Accept.
 */
class TPZUnsatAnalysis : public TPZLinearAnalysis {
public:
    /// Material of the domains
    using TMat = TPZMatUnsaturatedUPwPa<TPZPlasticStepVoigt<TPZYCMohrCoulombPV2, TPZElasticResponse>, TPZUnsaturatedMem>;

    struct StepResult {
        bool conv = false;
        int iter = 0;
        REAL eu = 0., normp = 0.; ///< last relative momentum residual and norm of the pressure residual
    };

    /// fixU: boundary id -> (fix u_x, fix u_y) for all the displacement connects of the boundary
    TPZUnsatAnalysis(TPZMultiphysicsCompMesh *mesh, const std::map<int, std::array<bool, 2>> &fixU);

    /// Linear solver of the analyses constructed afterwards: the sparse LU of MKL PARDISO when NeoPZ is built with MKL
    /// (default), or the skyline LU (SetUseSkyline(true), and always without MKL)
    static void SetUseSkyline(bool skyline);
    /// Whether the analyses constructed now use PARDISO
    static bool UsesPardiso();

    const std::vector<TMat *> &Materials() const { return fMats; }
    TPZMultiphysicsCompMesh *Mesh() { return fMPhys; }
    TPZGeoMesh *GMesh() { return fMPhys->Reference(); }

    void SetStrengthReductionFactor(REAL F);
    void SetGravityFactor(REAL lam);
    void SetGravityScalesSuction(bool b);

    std::vector<TPZAdmChunkVector<TPZUnsaturatedMem>> SaveMemory();
    void RestoreMemory(const std::vector<TPZAdmChunkVector<TPZUnsaturatedMem>> &s);
    void ResetMemory();

    /// node -> equation of p_w / p_a, node -> first displacement equation (vertices)
    const std::map<int64_t, int64_t> &WaterEq() const { return fEqW; }
    const std::map<int64_t, int64_t> &AirEq() const { return fEqA; }
    const std::map<int64_t, int64_t> &DispEq() const { return fEqU; }
    int64_t NEq() const { return fNEq; }
    bool IsPressure(int64_t eq) const { return fIsP[eq]; }
    /// displacement constraints of the boundaries (equation -> 0)
    const std::map<int64_t, REAL> &DisplacementConstraints() const { return fFixU; }
    /// constrained displacement equations of a boundary id, component comp (0: x, 1: y)
    std::set<int64_t> ConstrainedEqs(int id, int comp) const;
    /// all the displacement equations at their values in x (rigid skeleton)
    std::map<int64_t, REAL> AllDisplacements(const TPZFMatrix<STATE> &x) const;
    /// displacement constraints and all the pressures at their values in x (frozen pressures)
    std::map<int64_t, REAL> FrozenFixed(const TPZFMatrix<STATE> &x) const;

    /// vertices of the boundary elements with the id, with their tributary length (half of each edge; horizontal
    /// projection if horizontal = true)
    std::map<int64_t, REAL> BoundaryNodes(int id, bool horizontal = false) const;
    REAL Coord(int64_t node, int d) const { return fMPhys->Reference()->NodeVec()[node].Coord(d); }

    void SetTolerances(REAL tolu, REAL tolp, int maxit, int maxitfrozen) {
        fTolU = tolu;
        fTolP = tolp;
        fMaxIt = maxit;
        fMaxItFrozen = maxitfrozen;
    }
    void SetVerbose(int v) { fVerbose = v; }
    /// limit of Newton iterations with frozen pressures
    int MaxItFrozen() const { return fMaxItFrozen; }
    void SetMaxItFrozen(int maxit) { fMaxItFrozen = maxit; }

    /// converged solution (total u, p_w, p_a)
    TPZFMatrix<STATE> &X() { return fX; }

    /**
     * One backward-Euler step from the converged state: Newton with the equations of fixed eliminated and the
     * water inflow dt * fq (m3/s per unit thickness) added to the water equations. On convergence the solution is
     * loaded in the mesh; the memory and fX are updated only by Accept. A failure of the linear solver (singular
     * tangent) ends the step as not converged.
     */
    StepResult Solve(REAL dt, const std::map<int64_t, REAL> &fixed, const std::map<int64_t, REAL> &fq, bool frozen,
                     int linesearch = 4);

    /// Result of one increment of ArcLengthStep
    struct ArcLengthResult {
        bool conv = false;
        bool noRoot = false; ///< the constraint equation had no real root
        bool vertical = false; ///< null tangential solution: no point is plastic and u does not depend on the factor
        int iter = 0;
        REAL eu = 0.;        ///< last |r| / |fbar| on the free equations
        REAL lambda = 0.;    ///< load factor at the end of the increment
        TPZFMatrix<STATE> du; ///< incremental displacement of the increment (all the equations)
    };

    /// Options of an arc-length increment
    struct ArcLengthOptions {
        /// scale of the factor in the constraint Delta u^T Delta u + alpha^2 Delta lambda^2 = l^2, in units of displacement;
        /// 0 is the cylindrical constraint of Box 4.4
        REAL alpha;
        /// start from the previous increment scaled to the length l, instead of the tangent predictor
        bool secantPredictor;
        /// limit of Newton iterations; 0: the frozen limit of SetTolerances
        int maxit;
        /// evaluate the tangent predictor with the yield surface slightly shrunk: at the converged state the plastic
        /// points are exactly on the surface, and round-off alone decides between the elastic and the plastic tangent there
        bool plasticPredictor;
        ArcLengthOptions() : alpha(0.), secantPredictor(false), maxit(0), plasticPredictor(true) {}
    };

    /**
     * External load per unit of load factor with frozen pressures, fbar = R(x, lambda = 0) - R(x, lambda = 1) on all the
     * equations, so that R(x, lambda) = R(x, 0) - lambda fbar. The load factor is the gravity factor of the material:
     * fbar is the weight and the part of the equivalent pressure that the factor scales (SetGravityScalesSuction).
     * The gravity factor must be set again by the caller.
     */
    TPZFMatrix<STATE> LoadVector(const TPZFMatrix<STATE> &x);

    /**
     * One increment of the combined Newton-Raphson/cylindrical arc-length scheme, de Souza Neto, Peric and Owen,
     * Computational Methods for Plasticity (2008), Box 4.4, from the converged state X() with load factor lambdan:
     * finds u and lambda with R(u, lambda) = 0 and Delta u^T Delta u = l^2 on the free equations. All the pressures must
     * be in fixed (frozen pressures). The tangential solution is K_T dubar = +fbar, as in Eq. (4.115) (Box 4.4 prints it
     * with the opposite sign, which contradicts Eqs. 4.112 and 4.113). The sign of the predictor comes from the secant
     * criterion (4.123) with the previous incremental displacement dun and factor dlambdan; with no dun (0 rows), lambda
     * increases.
     * On convergence the solution is loaded in the mesh; the memory and X() are updated only by Accept.
     */
    ArcLengthResult ArcLengthStep(const std::map<int64_t, REAL> &fixed, const TPZFMatrix<STATE> &fbar, REAL lambdan, REAL l,
                                  const TPZFMatrix<STATE> &dun, REAL dlambdan = 0.,
                                  const ArcLengthOptions &options = ArcLengthOptions());

    /**
     * The same increment with the strength reduction factor F as the path parameter, from the converged state X() with
     * factor Fn: finds u and F with R(u, F) = 0 and Delta u^T Delta u = l^2. R is not linear in F, so the vector fbar of
     * Box 4.4 becomes fbar = -dR/dF at the current iterate, by a forward difference (one more residual assembly per
     * iteration). It only has contributions of the plastic points: where no point is plastic, u does not depend on F, the
     * path is vertical and the increment ends as not converged. Convergence as in Solve: |r| / |weight| < tolu.
     * The predictor is evaluated just above Fn: at Fn itself the points that were plastic in the converged state are
     * exactly on the yield surface and the tangent is elastic there, although they load plastically as F grows.
     */
    ArcLengthResult ArcLengthStepSRM(const std::map<int64_t, REAL> &fixed, REAL Fn, REAL l, const TPZFMatrix<STATE> &dun,
                                     REAL dFn = 0., const ArcLengthOptions &options = ArcLengthOptions());

    /// Updates the memory with the solution loaded in the mesh and stores it as the converged state
    void Accept();
    /// Loads the converged state in the mesh (after a step that did not converge)
    void Restore() { LoadSolution(fX); }

    /// Residual R (all equations, no nodal inflow) of x with the step dt (the vector of NeoPZ is -R)
    TPZFMatrix<STATE> Residual(const TPZFMatrix<STATE> &x, REAL dt, bool frozen);
    /// Product of the tangent dR/dX of x with d (all equations; no filter)
    TPZFMatrix<STATE> TangentTimes(const TPZFMatrix<STATE> &x, REAL dt, bool frozen, const TPZFMatrix<STATE> &d);
    /// Water (or air) residual of the solution loaded in the mesh: R_w / dt at a vertex with p_w prescribed is the
    /// inflow there
    TPZFMatrix<STATE> FullResidual();

    /// Loads the solution in the multiphysics mesh and in its atomic meshes
    void LoadSolution(const TPZFMatrix<STATE> &sol) override;
    using TPZLinearAnalysis::LoadSolution;

    /// u = 0 keeping the effective stress
    void ResetDisplacements();

    REAL Nodal(const TPZFMatrix<STATE> &x, int64_t node, int space) const {
        return x.GetVal((space == 1 ? fEqW : fEqA).at(node), 0);
    }
    REAL NodalDisplacement(const TPZFMatrix<STATE> &x, int64_t node, int comp) const {
        return x.GetVal(fEqU.at(node) + comp, 0);
    }
    /// reference of the momentum residual when the weight is null (kN)
    void SetForceReference(REAL f) { fForceRef = f; }

private:
    /**
     * Newton iterations of the arc-length increment (Box 4.4). setFactor sets the path parameter in the material.
     * tangentLoad(lambda, rhs, fbar) returns the vector -dR/dlambda of the current iterate in fbar, given rhs = -R at
     * lambda; it may assemble residuals, and must leave the factor at lambda. fref is the reference of the residual.
     * predictorShift: the predictor (first iteration) is evaluated at lambda_n (1 + predictorShift); setFactor is told
     * whether it is the predictor.
     */
    ArcLengthResult ArcLengthIterate(const std::map<int64_t, REAL> &fixed, REAL lambdan, REAL l, const TPZFMatrix<STATE> &dun,
                                     REAL dlambdan, const ArcLengthOptions &options, REAL fref, REAL predictorShift,
                                     const std::function<void(REAL, bool)> &setFactor,
                                     const std::function<void(REAL, const TPZFMatrix<STATE> &, TPZFMatrix<STATE> &)> &tangentLoad);
    void AssembleSystem();
    void AssembleFullResidual();
    void SetFilter(const std::map<int64_t, REAL> &fixed);
    int64_t Position(int64_t connect) const;
    void IdentifyEquations(const std::map<int, std::array<bool, 2>> &fixU);

    TPZMultiphysicsCompMesh *fMPhys;
    std::vector<TMat *> fMats;
    std::map<int64_t, int64_t> fEqW, fEqA, fEqU;
    std::map<int64_t, REAL> fFixU;
    std::map<int, std::array<std::set<int64_t>, 2>> fFixById;
    std::vector<bool> fIsP;
    std::set<int64_t> fFilterKey;
    TPZAutoPointer<TPZMatrix<STATE>> fTemplate;
    TPZFMatrix<STATE> fX;
    REAL fTolU = 1.e-7, fTolP = 1.e-6, fForceRef = 1.e-12;
    int fMaxIt = 40, fMaxItFrozen = 25, fVerbose = 0;
    int64_t fNEq = 0;
};

#endif
