/**
 * @file TPZPoroElastoPlasticUPAnalysis.h
 * @brief Incremental Newton driver for the coupled u-p consolidation problem (Sect. 5.4 of the article).
 */

#ifndef TPZPOROELASTOPLASTICUPANALYSIS_H
#define TPZPOROELASTOPLASTICUPANALYSIS_H

#include "TPZLinearAnalysis.h"
#include "TPZMatPoroElastoPlasticUP.h"
#include "pzfmatrix.h"
#include <vector>
#include <map>
#include <set>
#include <functional>
#include <array>

class TPZMultiphysicsCompMesh;

/**
 * @ingroup analysis
 * @brief Incremental solution of the backward-Euler u-p equations (26)-(27) with a monolithic Newton
 * method, step bisection and displacement control (routines SolveStepUP, AdvanceUP and
 * IterativeProcessUP of poro-camclay-fem.m).
 *
 * The loading is a sequence of states (t, lambda, u_c): time, load factor of the tractions of type
 * TPZMatPoroElastoPlasticUPBase::ENeumannU and value of the controlled displacement, applied to the
 * components registered with SetControlledDisplacement (boundary conditions of type EDirichletU or
 * EDirichletUDirectional whose Val2 is updated at each step).
 *
 * In each increment the Newton iterations start from the displacement of the previous step plus,
 * optionally, the increment of the previous step (predictor), and stop when the normalized residual
 * \f$ r = \|R\|/\max(\|f_{ext}\|,1) \f$ restricted to the free equations is below the tolerance
 * (1e-8); at least one correction is always performed. If the iterations do not converge (25) or the
 * local projection fails at an integration point, the increment is bisected (up to 8 levels).
 * The memory of the integration points is updated only after convergence (AcceptSolution).
 *
 * The equations of the Dirichlet conditions are eliminated from the linear systems with the equation
 * filter of the structural matrix (TPZEquationFilter) and their values are written in the solution before
 * the iterations, as in the reference implementation (SetEliminateDirichlet). The norms of the residual
 * and of the external forces are computed in the nodal basis of the serendipity element
 * (SetNodalResidualNorm), so that the normalized residual and the number of iterations are those of the
 * article (Table 8).
 *
 * The analysis is built without bandwidth optimization (the constructor passes
 * mustOptimizeBandwidth = false), so that the pressure equations are numbered after the displacement
 * equations: the LU decomposition of the skyline non-symmetric matrix (TPZSkylineNSymStructMatrix with
 * ELU, no pivoting) is then stable also in undrained steps with incompressible constituents (zero
 * pressure block), whose pivots are the Schur complement \f$Q^TK^{-1}Q\f$.
 *
 * The mesh solution holds the total displacement and pore pressure.
 */
class TPZPoroElastoPlasticUPAnalysis : public TPZLinearAnalysis {
public:

    /** @brief State of the loading: time, load factor and controlled displacement */
    struct TLoadState {
        REAL fTime;   ///< time t
        REAL fLambda; ///< load factor of the tractions
        REAL fUc;     ///< value of the controlled displacement
        TLoadState() : fTime(0.), fLambda(0.), fUc(0.) {}
        TLoadState(REAL t, REAL lambda, REAL uc) : fTime(t), fLambda(lambda), fUc(uc) {}
    };

    /** @brief Convergence record of a converged (sub)increment */
    struct TStepLog {
        TLoadState fState;             ///< state at the end of the increment
        std::vector<REAL> fResiduals;  ///< normalized residual of each evaluation (predictor included)
        int fLevel = 0;                ///< bisection level
    };

    /**
     * @brief Constructor
     * @param mesh multiphysics mesh (displacement space first, pressure space second)
     * @param mat the u-p material of the domain
     * @param out output stream for the log
     */
    TPZPoroElastoPlasticUPAnalysis(TPZMultiphysicsCompMesh *mesh, TPZMatPoroElastoPlasticUPBase *mat,
                                   std::ostream &out = std::cout);

    virtual ~TPZPoroElastoPlasticUPAnalysis() = default;

    /** @name Set up */
    /** @{ */
    /** @brief Registers the component (0, 1 or 2) of the boundary condition bcid as controlled displacement */
    void SetControlledDisplacement(int bcid, int component);

    /** @brief Newton tolerance on the normalized residual, maximum iterations and maximum bisection level */
    void SetNewtonParameters(REAL tol, int maxit, int maxbisections) {
        fTol = tol;
        fMaxIt = maxit;
        fMaxBisections = maxbisections;
    }

    /** @brief Enables the predictor (initial guess with the displacement increment of the previous step) */
    void SetPredictor(bool predictor) { fPredictor = predictor; }

    /** @brief Prints the residual of each iteration (verbose = 2) or a line per step (verbose = 1) */
    void SetVerbose(int verbose) { fVerbose = verbose; }

    /**
     * @brief If true (default), the equations with Dirichlet conditions are eliminated from the linear
     * systems with the equation filter of the structural matrix (TPZEquationFilter) and their values are
     * imposed exactly, as in the Python code; if false, the conditions are only imposed by the penalty
     * terms of the material (accuracy limited by the penalty number)
     */
    void SetEliminateDirichlet(bool eliminate) { fEliminateDirichlet = eliminate; }

    /**
     * @brief If true (default), the norms of the residual and of the external forces are computed with
     * the components in the nodal basis of the serendipity Q8/Hex20 element (vertex and mid-edge nodes)
     * instead of the hierarchical basis of NeoPZ, so that the normalized residual is the same of the
     * reference implementation. With N_c = phi_c - 1/2 sum_e phi_e and N_e = phi_e, the nodal components
     * are R_c - 1/2 sum_e R_e at the vertices and R_e at the edges.
     */
    void SetNodalResidualNorm(bool nodal) { fNodalNorm = nodal; }

    /** @brief Converts the displacement components of a vector from the hierarchical to the nodal basis */
    void ToNodalBasis(TPZFMatrix<STATE> &v) const;
    /** @} */

    /** @name Solution */
    /** @{ */
    /**
     * @brief Solves the sequence of states from start (IterativeProcessUP)
     * @param steps states at the end of each increment
     * @param monitor called after each converged increment with its index (1, 2, ...) and state; it is
     * also called with index 0 and the start state before the first increment
     * @param start initial state (time, load factor and controlled displacement of the current solution)
     * @return false if an increment failed after the maximum number of bisections
     */
    bool Run(const std::vector<TLoadState> &steps, const std::function<void(int, const TLoadState &)> &monitor,
             const TLoadState &start = TLoadState());

    /**
     * @brief Advances from s0 to s1 with bisection (AdvanceUP); the converged solution of the end state is
     * loaded in the mesh and the memory is updated
     * @param s0 initial state (current converged solution)
     * @param s1 final state
     * @param dupred predicted increment of the displacement equations (pressure entries are ignored)
     * @param level bisection level
     * @param[out] nbisect number of bisections performed
     */
    bool AdvanceStep(const TLoadState &s0, const TLoadState &s1, const TPZFMatrix<STATE> &dupred, int level,
                     int &nbisect);

    /**
     * @brief Newton iterations of one backward-Euler step (SolveStepUP)
     * @param s0 state of the converged solution
     * @param s1 state at the end of the step (Dt = s1.fTime - s0.fTime)
     * @param guess initial guess of the displacement equations
     * @param[out] residuals normalized residual of each evaluation
     * @return true if converged; the converged solution is then loaded in the mesh (memory not updated)
     */
    bool SolveStep(const TLoadState &s0, const TLoadState &s1, const TPZFMatrix<STATE> &guess,
                   std::vector<REAL> &residuals);

    /** @brief Updates the memory of the integration points with the solution loaded in the mesh */
    void AcceptSolution();

    /** @brief Loads a solution in the multiphysics mesh and in its atomic meshes */
    void LoadSolution(const TPZFMatrix<STATE> &sol) override;
    using TPZLinearAnalysis::LoadSolution;

    /** @brief Converged solution of the last increment (total displacement and pore pressure) */
    const TPZFMatrix<STATE> &ConvergedSolution() const { return fConverged; }

    /** @brief State of the last converged increment */
    const TLoadState &CurrentState() const { return fCurrent; }

    /** @brief Convergence records of all converged increments */
    const std::vector<TStepLog> &StepLog() const { return fStepLog; }
    void ClearStepLog() { fStepLog.clear(); }
    /** @} */

    /** @name Post-processing helpers */
    /** @{ */
    /**
     * @brief Residual of the momentum balance without the Dirichlet conditions,
     * \f$R_u = F_{int} - QP - f_b - \lambda f_t\f$ (all equations); its entries at the constrained
     * equations are the reactions
     */
    void ComputeReactionVector(TPZFMatrix<STATE> &Ru);

    /**
     * @brief Sum of the reactions in the component of the nodes of the boundary conditions bcids,
     * excluding the nodes that belong also to the boundary conditions exclude (ReactionByMarker)
     */
    REAL Reaction(const std::set<int> &bcids, int component, const std::set<int> &exclude = {});

    /** @brief Geometric node indices of the elements with the given material ids */
    std::set<int64_t> NodesOfMaterials(const std::set<int> &matids) const;

    /**
     * @brief Equation of a nodal degree of freedom (vertex connect of the node)
     * @param node geometric node index
     * @param space 0 displacement, 1 pore pressure
     * @param component component of the displacement (0 for the pressure)
     * @return equation index, or -1 if the node has no connect in the space
     */
    int64_t NodeEquation(int64_t node, int space, int component) const;

    /** @brief Value of a nodal degree of freedom in the solution loaded in the mesh */
    REAL NodalValue(int64_t node, int space, int component) const;

    /** @brief Flags of the equations with Dirichlet conditions (true = constrained) */
    const std::vector<bool> &ConstrainedEquations() const { return fConstrained; }

    /** @brief Flags of the pressure equations (true = pore pressure) */
    const std::vector<bool> &PressureEquations() const { return fIsPressure; }
    /** @} */

protected:
    /** @brief Builds the maps of equations (pressure flags, nodes) and the constrained equations */
    void IdentifyEquations();

    /** @brief Sets the values of the controlled displacement in the boundary conditions */
    void ApplyControlledDisplacement(REAL uc);

    /** @brief Writes the prescribed values in the constrained equations of a solution vector */
    void ImposeDirichletValues(TPZFMatrix<STATE> &sol);

    /** @brief External force vector \f$f_b+\lambda f_t\f$ */
    void ExternalForces(TPZFMatrix<STATE> &fext);

    /** @brief Activates the equation filter that removes the constrained equations (SetEliminateDirichlet) */
    void ApplyEquationFilter();

    /** @brief Assembles the residual of all the equations (the equation filter is suspended) */
    void AssembleFullResidual();

    TPZMultiphysicsCompMesh *fMPhys;
    TPZMatPoroElastoPlasticUPBase *fMat;
    std::vector<std::pair<int, int>> fControlled;
    REAL fTol = 1.e-8;
    int fMaxIt = 25;
    int fMaxBisections = 8;
    bool fPredictor = true;
    int fVerbose = 0;
    bool fEliminateDirichlet = true;
    bool fNodalNorm = true;
    /** @brief For each edge of the displacement mesh: first equation of the edge connect and of its two vertex connects */
    std::vector<std::array<int64_t, 3>> fEdgeEquations;
    std::vector<bool> fConstrained;
    std::vector<bool> fIsPressure;
    std::map<int64_t, std::pair<int64_t, int64_t>> fNodeConnects; ///< node -> (u connect, p connect) of the multiphysics mesh
    bool fEquationsIdentified = false;
    TPZFMatrix<STATE> fConverged;
    TLoadState fCurrent;
    std::vector<TStepLog> fStepLog;
};

#endif
